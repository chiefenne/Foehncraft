//! Broken-input corpus of the case readers (D202): the flow, mesh, scalar,
//! Poisson and benchmark case files of `examples/cases`, damaged in many
//! ways, must end in an error with a diagnostic code other than `GEN001`
//! (the code of "no more specific code"), never in a panic or a hang. A
//! broken mesh file named by a flow case ends in `MESH003`.
//!
//! The mesh readers have their own corpus (foehncraft-mesh,
//! `tests/broken_inputs.rs`); here the case layer is checked to give their
//! errors a code.

use foehncraft_flow::case::{FlowCase, check_case};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::Duration;

const LIMIT: Duration = Duration::from_secs(30);

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// Runs `f` on its own thread; a panic or a run past [`LIMIT`] fails, an
/// error must carry a specific code. Returns the code (`None`: no error).
fn guarded(label: &str, f: impl FnOnce() -> anyhow::Result<()> + Send + 'static) -> Option<String> {
    let (tx, rx) = mpsc::channel();
    std::thread::Builder::new()
        .stack_size(64 << 20)
        .spawn(move || {
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                f().map_err(|e| Box::new(diag::of_error(&e)))
            }));
            let _ = tx.send(r);
        })
        .expect("thread");
    match rx.recv_timeout(LIMIT) {
        Ok(Ok(Ok(()))) => None,
        Ok(Ok(Err(d))) => {
            assert_ne!(
                d.code,
                diag::GENERAL,
                "{label}: an error without a specific diagnostic code: {}",
                d.message
            );
            assert!(!d.message.trim().is_empty(), "{label}: an empty message");
            Some(d.code.to_string())
        }
        Ok(Err(p)) => {
            let m = p
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default();
            panic!("{label}: the reader panicked: {m}")
        }
        Err(_) => panic!("{label}: the reader did not finish within {LIMIT:?}"),
    }
}

fn spread(n: usize, k: usize) -> Vec<usize> {
    if n <= k {
        return (0..n).collect();
    }
    (0..k).map(|i| i * n / k).collect()
}

/// The damaged versions of a case file.
fn mutations(text: &str) -> Vec<(String, Vec<u8>)> {
    let lines: Vec<&str> = text.lines().collect();
    let join = |v: &[&str]| v.join("\n").into_bytes();
    let mut out: Vec<(String, Vec<u8>)> = vec![
        ("empty".into(), Vec::new()),
        ("not TOML".into(), b"this is not a case file\n".to_vec()),
        ("invalid UTF-8".into(), {
            let mut v = text.as_bytes().to_vec();
            v.splice(0..0, [0xff, 0xfe, 0x00, 0x80]);
            v
        }),
        ("an unclosed table".into(), {
            let mut v = text.as_bytes().to_vec();
            v.extend_from_slice(b"\n[mesh\n");
            v
        }),
        ("a key twice".into(), {
            let mut v = text.as_bytes().to_vec();
            v.extend_from_slice(b"\ncase_format = \"0.6\"\ncase_format = \"0.6\"\n");
            v
        }),
        ("an unknown table".into(), {
            let mut v = text.as_bytes().to_vec();
            v.extend_from_slice(b"\n[no_such_table]\nx = 1\n");
            v
        }),
        ("deep nesting".into(), {
            let mut v = b"x = ".to_vec();
            v.extend(std::iter::repeat_n(b'[', 20_000));
            v
        }),
    ];
    for k in spread(lines.len(), 8) {
        out.push((format!("truncated after line {k}"), join(&lines[..k])));
    }
    // the key lines: values replaced, lines dropped and repeated
    let keyed: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, l)| {
            let t = l.trim_start();
            !t.starts_with('#') && t.contains('=') && !t.starts_with('[')
        })
        .map(|(i, _)| i)
        .collect();
    const BAD: [&str; 12] = [
        "-1",
        "0",
        "1e999",
        "-1e308",
        "nan",
        "inf",
        "\"x\"",
        "[]",
        "[1, \"a\"]",
        "{}",
        "true",
        "18446744073709551616",
    ];
    for &i in &spread(keyed.len(), 16)
        .iter()
        .map(|&k| keyed[k])
        .collect::<Vec<_>>()
    {
        let (key, _) = lines[i].split_once('=').unwrap();
        for bad in BAD {
            let mut v = lines.clone();
            let line = format!("{key}= {bad}");
            v[i] = &line;
            out.push((format!("line {i}: {} = {bad}", key.trim()), join(&v)));
        }
        let mut v = lines.clone();
        v.remove(i);
        out.push((format!("line {i} dropped"), join(&v)));
        let mut v = lines.clone();
        v.insert(i, lines[i]);
        out.push((format!("line {i} twice"), join(&v)));
    }
    // the table headers dropped (keys land in the table above)
    for (i, l) in lines.iter().enumerate() {
        if l.trim_start().starts_with('[') {
            let mut v = lines.clone();
            v.remove(i);
            out.push((format!("header {} dropped", l.trim()), join(&v)));
        }
    }
    out
}

fn cases(dir: &str) -> Vec<PathBuf> {
    fn walk(d: &Path, out: &mut Vec<PathBuf>) {
        for e in std::fs::read_dir(d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                walk(&p, out);
            } else if p.extension().is_some_and(|x| x == "toml") {
                out.push(p);
            }
        }
    }
    let mut v = Vec::new();
    walk(&root().join("examples/cases").join(dir), &mut v);
    v.sort();
    assert!(!v.is_empty(), "no cases in {dir}");
    v
}

fn scratch(tag: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("broken-case-{tag}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Every damaged version of every case of `dir` through `read`; the
/// damaged file is written next to the original's name in a scratch
/// directory (so relative paths of the original no longer resolve, which
/// must also end in a code).
fn corpus(dir: &str, sample: usize, read: fn(&Path) -> anyhow::Result<()>) {
    let all = cases(dir);
    let picks: Vec<PathBuf> = spread(all.len(), sample)
        .into_iter()
        .map(|i| all[i].clone())
        .collect();
    let tmp = scratch(dir);
    for case in picks {
        let text = std::fs::read_to_string(&case).unwrap();
        let name = case.file_name().unwrap().to_owned();
        for (what, m) in mutations(&text) {
            let p = tmp.join(&name);
            std::fs::write(&p, &m).unwrap();
            let label = format!("{dir} case {}, {what}", case.display());
            guarded(&label, move || read(&p));
        }
    }
}

/// The flow cases are damaged in memory and read against their own
/// directory, so that the mesh and table paths resolve and the validation
/// is reached; the invalid UTF-8 goes through a file.
#[test]
fn flow_case_reader_survives_the_corpus() {
    let all = cases("flow");
    let tmp = scratch("flow-utf8");
    let mut errors = 0;
    for case in &all {
        let text = std::fs::read_to_string(case).unwrap();
        let dir = case.parent().unwrap().to_path_buf();
        for (what, m) in mutations(&text) {
            let label = format!("flow case {}, {what}", case.display());
            let got = match String::from_utf8(m) {
                Ok(t) => {
                    let d = dir.clone();
                    guarded(&label, move || FlowCase::from_toml(&t, &d).map(|_| ()))
                }
                Err(e) => {
                    let p = tmp.join(case.file_name().unwrap());
                    std::fs::write(&p, e.as_bytes()).unwrap();
                    guarded(&label, move || FlowCase::from_file(&p).map(|_| ()))
                }
            };
            errors += got.is_some() as usize;
        }
    }
    assert!(errors > 1000, "only {errors} errors: is the damage seen?");
}

#[test]
fn mesh_case_reader_survives_the_corpus() {
    corpus("mesh", 12, |p| {
        foehncraft_mesh::case::MeshCase::from_file(p)
            .map(|_| ())
            .map_err(diag::case_error)
    });
}

#[test]
fn scalar_case_reader_survives_the_corpus() {
    corpus("scalar", 8, |p| {
        foehncraft_fv::case::ScalarCase::from_file(p)
            .map(|_| ())
            .map_err(diag::case_error)
    });
}

#[test]
fn poisson_case_reader_survives_the_corpus() {
    corpus("poisson", 5, |p| {
        foehncraft_fv::poisson::PoissonCase::from_file(p)
            .map(|_| ())
            .map_err(diag::case_error)
    });
}

#[test]
fn bench_case_reader_survives_the_corpus() {
    corpus("bench", 10, |p| {
        foehncraft_bench::case::G1Case::from_file(p)
            .map(|_| ())
            .map_err(diag::case_error)
    });
}

/// A flow case whose mesh file is broken (the hand-written inputs of the
/// mesh corpus, damaged copies of a valid mesh, a missing file) stops with
/// `MESH003`.
#[test]
fn broken_mesh_files_end_in_mesh003() {
    let meshes = root().join("crates/foehncraft-mesh/tests/data");
    let tmp = scratch("mesh003");
    let mut files: Vec<PathBuf> = std::fs::read_dir(meshes.join("broken"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x != "md"))
        .collect();
    files.sort();
    // damaged copies of valid meshes, written with their extension
    for (src, cut) in [
        ("box_tet_ascii.msh", 3),
        ("box_tet_binary.msh", 2),
        ("channel_2d.su2", 2),
        ("vtk/channel_mixed_zlib.vtu", 2),
        ("fluent_msh/channel_hex.msh", 2),
    ] {
        let b = std::fs::read(meshes.join(src)).unwrap();
        let name = Path::new(src)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let p = tmp.join(format!("half_{name}"));
        std::fs::write(&p, &b[..b.len() / cut]).unwrap();
        files.push(p);
    }
    files.push(tmp.join("missing.msh"));
    files.push(tmp.join("no_extension"));
    assert!(files.len() >= 12);
    for f in files {
        // a literal TOML string: no escapes in Windows paths
        let text = format!(
            "case_format = \"0.6\"\n[mesh]\nfile = '{}'\n[fluid]\ndensity = 1.0\n\
             viscosity = 1e-3\n[boundary.wall]\ntype = \"wall\"\n",
            f.display()
        );
        let label = format!("flow case with the mesh {}", f.display());
        let path = tmp.join("case.toml");
        let dir = tmp.clone();
        let code = guarded(&label, move || {
            let case = FlowCase::from_toml(&text, &dir)?;
            check_case(&case, &path, |_| {})
        });
        assert_eq!(code.as_deref(), Some("MESH003"), "{label}");
    }
}
