//! The live control port of a running case (L19, D198 (f)).
//!
//! A run started with `control` listens on a TCP port of the loopback
//! interface (`127.0.0.1`, never another address) so that an agent host or a
//! person can ask a running solver for its state, stop it or change its
//! run-time settings, without a window and without the control files of
//! D40. The port has the same meaning as the files, applied at the same
//! place (after an outer iteration or time step, D76):
//!
//! | operation | effect |
//! |---|---|
//! | `status` | the last iteration or step row, the stop request and the changes made through the port |
//! | `stop` | as `.stop`: restart file, final output, exit code 3 (`RUN022`) |
//! | `restart` | writes `restart.h5` at the end of the iteration and goes on |
//! | `set` | changes run-time settings ([`crate::control::RUNTIME_KEYS`], without `run.limits.*`: a limit is the person's, an agent cannot raise its own); a setting that is invalid when it is applied is reported and left out |
//!
//! Access: the run writes two files into its output directory, `control.json`
//! (the port, the process, the token file) and `control.token` (a random token;
//! on Unix readable by its user only, mode 0600). A request carries the
//! token; without it, or with another one, the request is refused. Both
//! files are removed when the run ends. The protocol is one JSON object per
//! line, one request and one answer per connection:
//!
//! ```text
//! {"token": "...", "op": "set", "set": {"solver.momentum.relaxation": 0.3}}
//! {"ok": true, "queued": ["solver.momentum.relaxation"]}
//! ```
//!
//! Every change that the run applies or cannot apply is an `event` row of
//! the run record with the equivalent `[[schedule]]` entry (D128), so a
//! person can reproduce a steered run from its case file.

use anyhow::{Context, Result};
use serde_json::{Value, json};
use std::io::{BufRead as _, Read as _, Write as _};
use std::net::{Ipv4Addr, SocketAddr, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Format of `control.json` and of the answer of `status`.
pub const CONTROL_FORMAT: &str = "foehn.control";
/// Format of the status.
pub const STATUS_FORMAT: &str = "foehn.run-status";
/// The file with the port and the process.
pub const CONTROL_FILE: &str = "control.json";
/// The file with the access token.
pub const TOKEN_FILE: &str = "control.token";

/// A setting asked for through the port, applied at the next boundary.
#[derive(Debug, Clone)]
pub struct Pending {
    /// `section.key`.
    pub key: String,
    /// The new value.
    pub value: toml::Value,
}

#[derive(Debug, Default)]
struct Inner {
    /// The last iteration, step or sweep row of the run record.
    last: Option<Value>,
    /// Settings waiting for the next boundary.
    pending: Vec<Pending>,
    /// What the run did with the requests, newest last.
    log: Vec<Value>,
    /// A stop was asked for.
    stop_requested: bool,
    /// The state word of the run.
    state: &'static str,
    started: Option<Instant>,
}

/// What the port and the loops of a run share. A clone is another handle.
#[derive(Debug, Clone)]
pub struct Live(Arc<Mutex<Inner>>);

impl Default for Live {
    fn default() -> Self {
        Self(Arc::new(Mutex::new(Inner {
            state: "not-started",
            ..Inner::default()
        })))
    }
}

impl Live {
    fn with<R>(&self, f: impl FnOnce(&mut Inner) -> R) -> R {
        match self.0.lock() {
            Ok(mut g) => f(&mut g),
            Err(p) => f(&mut p.into_inner()),
        }
    }

    /// The run is going: the clock of the status starts.
    pub(crate) fn begin(&self) {
        self.with(|s| {
            s.state = "running";
            s.started = Some(Instant::now());
        });
    }

    /// The run record row `row` went by; keeps the last iteration or step.
    pub(crate) fn observe(&self, row: &Value) {
        if matches!(row["type"].as_str(), Some("iteration" | "step")) {
            self.with(|s| s.last = Some(row.clone()));
        }
    }

    /// The settings asked for since the last call.
    pub(crate) fn take_settings(&self) -> Vec<Pending> {
        self.with(|s| std::mem::take(&mut s.pending))
    }

    /// A stop was asked for through the port.
    pub(crate) fn stop_requested(&self) -> bool {
        self.with(|s| s.stop_requested)
    }

    /// What the run did with a request (`applied`, `rejected`), for the
    /// status.
    pub(crate) fn note(&self, entry: Value) {
        self.with(|s| {
            s.log.push(entry);
            let n = s.log.len();
            if n > 200 {
                s.log.drain(..n - 200);
            }
        });
    }

    /// The run ended with `state`.
    pub(crate) fn end(&self, state: &'static str) {
        self.with(|s| s.state = state);
    }

    fn status(&self) -> Value {
        self.with(|s| {
            json!({
                "format": STATUS_FORMAT,
                "state": s.state,
                "pid": std::process::id(),
                "wall_time": s.started.map(|t| t.elapsed().as_secs_f64()),
                "stop_requested": s.stop_requested,
                "last": s.last,
                "control": s.log,
                "pending": s.pending.iter().map(|p| p.key.clone()).collect::<Vec<_>>(),
            })
        })
    }
}

/// Keys the port may change: those of a re-read (D40) except the limits.
pub fn settable(key: &str) -> bool {
    crate::control::RUNTIME_KEYS.contains(&key) && !key.starts_with("run.limits.")
}

fn to_toml(v: &Value) -> Option<toml::Value> {
    Some(match v {
        Value::Bool(b) => toml::Value::Boolean(*b),
        Value::Number(n) => n
            .as_i64()
            .map(toml::Value::Integer)
            .or_else(|| n.as_f64().map(toml::Value::Float))?,
        Value::String(s) => toml::Value::String(s.clone()),
        Value::Array(a) => toml::Value::Array(a.iter().map(to_toml).collect::<Option<_>>()?),
        Value::Null | Value::Object(_) => return None,
    })
}

/// The answer to one request line.
fn handle(line: &str, token: &str, live: &Live, interrupts: &crate::engine::Interrupts) -> Value {
    let refuse = |msg: &str| json!({"ok": false, "code": "TOOL007", "error": msg});
    let req: Value = match serde_json::from_str(line) {
        Ok(v) => v,
        Err(e) => return refuse(&format!("the request is not JSON: {e}")),
    };
    if !constant_time_eq(
        req["token"].as_str().unwrap_or("").as_bytes(),
        token.as_bytes(),
    ) {
        return refuse("the token is missing or wrong");
    }
    match req["op"].as_str().unwrap_or("") {
        "status" => {
            let mut s = live.status();
            s["ok"] = json!(true);
            s
        }
        "stop" => {
            live.with(|s| s.stop_requested = true);
            json!({"ok": true, "stop_requested": true})
        }
        "restart" => {
            interrupts.request_restart();
            json!({"ok": true, "restart_requested": true})
        }
        "set" => {
            let Some(map) = req["set"].as_object().filter(|m| !m.is_empty()) else {
                return refuse("`set` needs an object of settings, {\"key\": value}");
            };
            let mut pending = Vec::new();
            for (k, v) in map {
                if !settable(k) {
                    return refuse(&format!(
                        "`{k}` cannot be changed through the control port: only run-time settings can (the keys of a re-read, D40), and never `run.limits.*`"
                    ));
                }
                let Some(value) = to_toml(v) else {
                    return refuse(&format!(
                        "the value of `{k}` is not a number, string, boolean or list"
                    ));
                };
                pending.push(Pending {
                    key: k.clone(),
                    value,
                });
            }
            let queued: Vec<String> = pending.iter().map(|p| p.key.clone()).collect();
            live.with(|s| s.pending.extend(pending));
            json!({"ok": true, "queued": queued, "applied_at": "the end of the current iteration or step"})
        }
        other => refuse(&format!(
            "unknown operation `{other}` (status, stop, restart, set)"
        )),
    }
}

fn constant_time_eq(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).fold(0u8, |d, (x, y)| d | (x ^ y)) == 0
}

/// A random token of 128 bits as hexadecimal text.
fn new_token() -> String {
    let mut b = [0u8; 16];
    let from_os = std::fs::File::open("/dev/urandom")
        .and_then(|mut f| f.read_exact(&mut b))
        .is_ok();
    if !from_os {
        use std::hash::{BuildHasher as _, Hasher as _};
        for chunk in b.chunks_mut(8) {
            let mut h = std::collections::hash_map::RandomState::new().build_hasher();
            h.write_u128(
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_nanos()),
            );
            h.write_u32(std::process::id());
            chunk.copy_from_slice(&h.finish().to_le_bytes()[..chunk.len()]);
        }
    }
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn write_private(path: &Path, text: &str) -> Result<()> {
    let mut o = std::fs::OpenOptions::new();
    o.write(true).create(true).truncate(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt as _;
        o.mode(0o600);
    }
    let mut f = o
        .open(path)
        .with_context(|| format!("writing {}", path.display()))?;
    f.write_all(text.as_bytes())?;
    Ok(())
}

/// The listening port of a run. Dropping it closes the port and removes
/// `control.json` and `control.token`.
pub struct Port {
    dir: PathBuf,
    addr: SocketAddr,
    done: Arc<AtomicBool>,
    thread: Option<std::thread::JoinHandle<()>>,
    live: Live,
}

impl Port {
    /// Opens the port for the run whose output directory is `dir`.
    pub fn open(dir: &Path, live: &Live, interrupts: &crate::engine::Interrupts) -> Result<Self> {
        let listener = TcpListener::bind((Ipv4Addr::LOCALHOST, 0))
            .context("opening the control port on 127.0.0.1")?;
        let addr = listener.local_addr()?;
        let token = new_token();
        write_private(&dir.join(TOKEN_FILE), &token)?;
        write_private(
            &dir.join(CONTROL_FILE),
            &serde_json::to_string_pretty(&json!({
                "format": CONTROL_FORMAT,
                "format_version": 1,
                "host": "127.0.0.1",
                "port": addr.port(),
                "pid": std::process::id(),
                "token_file": TOKEN_FILE,
            }))?,
        )?;
        let done = Arc::new(AtomicBool::new(false));
        live.begin();
        let thread = {
            let (done, live, interrupts) = (done.clone(), live.clone(), interrupts.clone());
            std::thread::spawn(move || {
                for conn in listener.incoming() {
                    if done.load(Ordering::SeqCst) {
                        break;
                    }
                    let Ok(mut s) = conn else { continue };
                    let _ = s.set_read_timeout(Some(Duration::from_secs(5)));
                    let mut line = String::new();
                    let read = std::io::BufReader::new(&s)
                        .take(1 << 16)
                        .read_line(&mut line);
                    let answer = match read {
                        Ok(n) if n > 0 => handle(&line, &token, &live, &interrupts),
                        _ => json!({"ok": false, "code": "TOOL007", "error": "no request"}),
                    };
                    let _ = writeln!(s, "{answer}");
                }
            })
        };
        Ok(Self {
            dir: dir.to_path_buf(),
            addr,
            done,
            thread: Some(thread),
            live: live.clone(),
        })
    }
}

impl Drop for Port {
    fn drop(&mut self) {
        self.done.store(true, Ordering::SeqCst);
        // unblock `accept`
        let _ = TcpStream::connect_timeout(&self.addr, Duration::from_millis(500));
        if let Some(t) = self.thread.take() {
            let _ = t.join();
        }
        let _ = std::fs::remove_file(self.dir.join(CONTROL_FILE));
        let _ = std::fs::remove_file(self.dir.join(TOKEN_FILE));
        self.live.end("finished");
    }
}

/// Sends one request to the run in `dir` (the client side: used by the tools
/// `run.status`, `run.stop` and `run.set`). `req` is the request without the
/// token.
pub fn request(dir: &Path, mut req: Value) -> Result<Value> {
    let info: Value = serde_json::from_str(
        &std::fs::read_to_string(dir.join(CONTROL_FILE))
            .with_context(|| format!("no run listens in {} (no {CONTROL_FILE}): the run has ended, was started without the control port, or the directory is wrong", dir.display()))?,
    )?;
    let port = info["port"]
        .as_u64()
        .and_then(|p| u16::try_from(p).ok())
        .context("control.json has no port")?;
    let token =
        std::fs::read_to_string(dir.join(info["token_file"].as_str().unwrap_or(TOKEN_FILE)))
            .with_context(|| format!("reading the token file in {}", dir.display()))?;
    req["token"] = json!(token.trim());
    let mut s = TcpStream::connect_timeout(
        &SocketAddr::from((Ipv4Addr::LOCALHOST, port)),
        Duration::from_secs(3),
    )
    .with_context(|| {
        format!("connecting to the run on 127.0.0.1:{port} (is the run still going?)")
    })?;
    s.set_read_timeout(Some(Duration::from_secs(10)))?;
    writeln!(s, "{req}")?;
    let mut line = String::new();
    std::io::BufReader::new(&s).read_line(&mut line)?;
    serde_json::from_str(&line).context("the run answered with something that is not JSON")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn port_answers_and_cleans_up() {
        let dir = std::env::temp_dir().join(format!("foehn-live-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let live = Live::default();
        let ints = crate::engine::Interrupts::new();
        let port = Port::open(&dir, &live, &ints).unwrap();
        live.observe(&json!({"type": "iteration", "iteration": 7, "continuity": 0.5}));
        let s = request(&dir, json!({"op": "status"})).unwrap();
        assert_eq!(s["ok"], true);
        assert_eq!(s["state"], "running");
        assert_eq!(s["last"]["iteration"], 7);
        // limits are not settable, run-time keys are
        let r = request(
            &dir,
            json!({"op": "set", "set": {"run.limits.max_wall_time": 1e9}}),
        )
        .unwrap();
        assert_eq!(r["ok"], false);
        let r = request(
            &dir,
            json!({"op": "set", "set": {"solver.momentum.relaxation": 0.3}}),
        )
        .unwrap();
        assert_eq!(r["queued"][0], "solver.momentum.relaxation");
        let p = live.take_settings();
        assert_eq!(p.len(), 1);
        assert_eq!(p[0].value.as_float(), Some(0.3));
        // a wrong token is refused
        let info = std::fs::read_to_string(dir.join(CONTROL_FILE)).unwrap();
        let port_no = serde_json::from_str::<Value>(&info).unwrap()["port"]
            .as_u64()
            .unwrap() as u16;
        let mut c = TcpStream::connect((Ipv4Addr::LOCALHOST, port_no)).unwrap();
        writeln!(c, "{}", json!({"token": "nope", "op": "stop"})).unwrap();
        let mut line = String::new();
        std::io::BufReader::new(&c).read_line(&mut line).unwrap();
        assert!(line.contains("\"ok\":false"));
        assert!(!live.stop_requested());
        request(&dir, json!({"op": "stop"})).unwrap();
        assert!(live.stop_requested());
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt as _;
            let m = std::fs::metadata(dir.join(TOKEN_FILE)).unwrap();
            assert_eq!(m.permissions().mode() & 0o777, 0o600);
        }
        drop(port);
        assert!(!dir.join(CONTROL_FILE).exists());
        assert!(request(&dir, json!({"op": "status"})).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
