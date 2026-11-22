//! The action table of agent calls (D198 (a), (b), (c)).
//!
//! A host that lets an assistant call the tools (the MCP server `foehn mcp`,
//! L18; later the desktop app, L7) asks this module, before each call, what
//! the call would do and whether the project allows it. The **action
//! kinds** are: start a run, exceed the estimate caps of the project (D81),
//! write outside the write roots, overwrite or delete earlier results, and
//! publish. Each kind is set to `default`, `allow`, `ask` or `refuse` in the
//! project's action table (`foehn-mcp.toml`); an answer "always" to a
//! question of the host is written back into the table, an answer "once"
//! holds for the one call.
//!
//! ```toml
//! # foehn-mcp.toml: the action table of the project (D198 (b))
//! [actions]
//! start-run = "ask"
//! exceed-estimate = "ask"
//! write-outside-root = "refuse"
//! overwrite-results = "ask"
//! publish = "refuse"
//!
//! # the estimate caps: a run whose estimate (D81) exceeds one is the kind
//! # `exceed-estimate`
//! [estimate]
//! max_cells = 2000000
//! max_host_memory = 8e9   # bytes
//! max_output = 5e9        # bytes
//! max_wall_time = 3600    # s, when the estimate has a wall time
//! ```
//!
//! A run started through a host must also bound itself: `[run.limits]
//! max_wall_time` and `max_output` are mandatory (D198 (c), `MCP002`).
//! Nothing here decides alone what a tool does: what a call writes comes
//! from the tool itself ([`crate::effects`]).

use crate::{Action, Effects, ToolSpec};
use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic, Fix};
use foehncraft_flow::engine::Roots;
use schemars::JsonSchema;
use serde::Serialize;
use std::path::{Path, PathBuf};

/// The default file name of the action table.
pub const POLICY_FILE: &str = "foehn-mcp.toml";

/// What a call may do that needs a rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Kind {
    /// The call starts a run (tools of action kind `start-run`).
    StartRun,
    /// The estimate of the run (D81) exceeds a cap of `[estimate]`.
    ExceedEstimate,
    /// The call writes outside the write roots.
    WriteOutsideRoot,
    /// The call writes into a directory with earlier results, or over a file.
    OverwriteResults,
    /// The call publishes something outside the computer (no tool does yet).
    Publish,
}

impl Kind {
    /// Every kind, in table order.
    pub const ALL: [Kind; 5] = [
        Kind::StartRun,
        Kind::ExceedEstimate,
        Kind::WriteOutsideRoot,
        Kind::OverwriteResults,
        Kind::Publish,
    ];

    /// The key in `[actions]`.
    pub fn name(self) -> &'static str {
        match self {
            Kind::StartRun => "start-run",
            Kind::ExceedEstimate => "exceed-estimate",
            Kind::WriteOutsideRoot => "write-outside-root",
            Kind::OverwriteResults => "overwrite-results",
            Kind::Publish => "publish",
        }
    }

    /// What `default` means for the kind.
    pub fn default_setting(self) -> Setting {
        match self {
            Kind::StartRun | Kind::ExceedEstimate | Kind::OverwriteResults => Setting::Ask,
            Kind::WriteOutsideRoot | Kind::Publish => Setting::Refuse,
        }
    }

    fn of(name: &str) -> Option<Kind> {
        Kind::ALL.into_iter().find(|k| k.name() == name)
    }
}

/// The rule of a kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Setting {
    /// The built-in rule of the kind ([`Kind::default_setting`]).
    Default,
    /// Allowed without a question.
    Allow,
    /// The host asks the person before the call.
    Ask,
    /// Refused (`MCP003`).
    Refuse,
}

impl Setting {
    /// Lower-case name.
    pub fn name(self) -> &'static str {
        match self {
            Setting::Default => "default",
            Setting::Allow => "allow",
            Setting::Ask => "ask",
            Setting::Refuse => "refuse",
        }
    }

    fn of(name: &str) -> Option<Setting> {
        [
            Setting::Default,
            Setting::Allow,
            Setting::Ask,
            Setting::Refuse,
        ]
        .into_iter()
        .find(|s| s.name() == name)
    }
}

/// The estimate caps of `[estimate]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, JsonSchema)]
pub struct Caps {
    /// Cells of the largest level.
    pub max_cells: usize,
    /// Host memory (bytes).
    pub max_host_memory: f64,
    /// Output volume (bytes).
    pub max_output: f64,
    /// Wall time (s), checked when the estimate has one.
    pub max_wall_time: f64,
}

impl Default for Caps {
    fn default() -> Self {
        Self {
            max_cells: 2_000_000,
            max_host_memory: 8e9,
            max_output: 5e9,
            max_wall_time: 3600.0,
        }
    }
}

/// The action table of a project.
#[derive(Debug, Clone, PartialEq, Serialize, JsonSchema)]
pub struct Policy {
    /// The rule of each kind, in table order.
    pub actions: Vec<(Kind, Setting)>,
    /// The estimate caps.
    pub estimate: Caps,
}

impl Default for Policy {
    fn default() -> Self {
        Self {
            actions: Kind::ALL
                .into_iter()
                .map(|k| (k, Setting::Default))
                .collect(),
            estimate: Caps::default(),
        }
    }
}

fn policy_error(msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::MCP_POLICY, msg))
}

impl Policy {
    /// The rule of `kind` as set (`default` included).
    pub fn setting(&self, kind: Kind) -> Setting {
        self.actions
            .iter()
            .find(|(k, _)| *k == kind)
            .map_or(Setting::Default, |(_, s)| *s)
    }

    /// The rule of `kind` with `default` resolved.
    pub fn rule(&self, kind: Kind) -> Setting {
        match self.setting(kind) {
            Setting::Default => kind.default_setting(),
            s => s,
        }
    }

    /// Sets the rule of `kind`.
    pub fn set(&mut self, kind: Kind, s: Setting) {
        match self.actions.iter_mut().find(|(k, _)| *k == kind) {
            Some(e) => e.1 = s,
            None => self.actions.push((kind, s)),
        }
    }

    /// Reads a table from TOML text.
    pub fn parse(text: &str) -> Result<Policy> {
        let v: toml::Table =
            toml::from_str(text).map_err(|e| policy_error(format!("not TOML: {e}")))?;
        let mut p = Policy::default();
        for (k, val) in &v {
            match (k.as_str(), val) {
                ("actions", toml::Value::Table(t)) => {
                    for (name, s) in t {
                        let kind = Kind::of(name).ok_or_else(|| {
                            policy_error(format!(
                                "[actions] {name}: not an action kind ({})",
                                Kind::ALL.map(Kind::name).join(", ")
                            ))
                        })?;
                        let setting = s.as_str().and_then(Setting::of).ok_or_else(|| {
                            policy_error(format!(
                                "[actions] {name}: \"default\", \"allow\", \"ask\" or \"refuse\""
                            ))
                        })?;
                        p.set(kind, setting);
                    }
                }
                ("estimate", toml::Value::Table(t)) => {
                    for (name, x) in t {
                        let num = x
                            .as_float()
                            .or_else(|| x.as_integer().map(|i| i as f64))
                            .filter(|f| *f > 0.0)
                            .ok_or_else(|| {
                                policy_error(format!("[estimate] {name}: a number > 0"))
                            })?;
                        match name.as_str() {
                            "max_cells" => p.estimate.max_cells = num as usize,
                            "max_host_memory" => p.estimate.max_host_memory = num,
                            "max_output" => p.estimate.max_output = num,
                            "max_wall_time" => p.estimate.max_wall_time = num,
                            _ => {
                                return Err(policy_error(format!(
                                    "[estimate] {name}: not a cap (max_cells, max_host_memory, max_output, max_wall_time)"
                                )));
                            }
                        }
                    }
                }
                _ => {
                    return Err(policy_error(format!(
                        "`{k}`: the table has only [actions] and [estimate]"
                    )));
                }
            }
        }
        Ok(p)
    }

    /// Reads the table at `path`; a missing file is the default table.
    pub fn load(path: &Path) -> Result<Policy> {
        match std::fs::read_to_string(path) {
            Ok(t) => Policy::parse(&t).map_err(|e| {
                let mut d = diag::of_error(&e);
                d.message = format!("{}: {}", path.display(), d.message);
                diag::err(d)
            }),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(Policy::default()),
            Err(e) => Err(policy_error(format!("reading {}: {e}", path.display()))),
        }
    }

    /// The table as TOML, with comments.
    pub fn to_toml(&self) -> String {
        let mut s = String::from(
            "# The action table of this project for agent calls (foehn mcp, D198 (b)).\n\
             # Each kind: \"default\", \"allow\", \"ask\" or \"refuse\". Answers \"always\"\n\
             # of the approval questions are written here; edit by hand at any time.\n\
             [actions]\n",
        );
        for k in Kind::ALL {
            s.push_str(&format!(
                "{} = \"{}\"   # default: {}\n",
                k.name(),
                self.setting(k).name(),
                k.default_setting().name()
            ));
        }
        let e = &self.estimate;
        s.push_str(&format!(
            "\n# A run whose estimate (D81) exceeds a cap is the kind exceed-estimate.\n\
             [estimate]\n\
             max_cells = {}\n\
             max_host_memory = {:e}   # bytes\n\
             max_output = {:e}        # bytes\n\
             max_wall_time = {}       # s, when the estimate has a wall time\n",
            e.max_cells, e.max_host_memory, e.max_output, e.max_wall_time
        ));
        s
    }

    /// Writes the table to `path`.
    pub fn save(&self, path: &Path) -> Result<()> {
        if let Some(d) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
            std::fs::create_dir_all(d)
                .map_err(|e| policy_error(format!("creating {}: {e}", d.display())))?;
        }
        std::fs::write(path, self.to_toml())
            .map_err(|e| policy_error(format!("writing {}: {e}", path.display())))
    }
}

/// One rule a call meets.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct Need {
    /// The action kind.
    pub kind: Kind,
    /// The rule of the table (`default` resolved).
    pub rule: Setting,
    /// Why the call is of this kind, for the person.
    pub reason: String,
    /// For `write-outside-root`: the directory the call would write.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub path: Option<PathBuf>,
}

/// What the table says about one call.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct Verdict {
    /// The rules the call meets, in table order.
    pub needs: Vec<Need>,
    /// A refusal that no answer changes (`MCP002`, mandatory run limits).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub refused: Option<Diagnostic>,
}

impl Verdict {
    /// The kinds that the table refuses.
    pub fn refusals(&self) -> impl Iterator<Item = &Need> {
        self.needs.iter().filter(|n| n.rule == Setting::Refuse)
    }

    /// The kinds that need the person's answer.
    pub fn questions(&self) -> impl Iterator<Item = &Need> {
        self.needs.iter().filter(|n| n.rule == Setting::Ask)
    }

    /// The directories outside the write roots that the call may write once
    /// allowed.
    pub fn extra_write_roots(&self) -> Vec<PathBuf> {
        self.needs
            .iter()
            .filter(|n| n.kind == Kind::WriteOutsideRoot)
            .filter_map(|n| n.path.clone())
            .collect()
    }
}

fn gb(b: f64) -> String {
    if b >= 1e9 {
        format!("{:.1} GB", b / 1e9)
    } else {
        format!("{:.0} MB", b / 1e6)
    }
}

/// A directory with earlier results, or an existing file.
fn holds_results(p: &Path) -> bool {
    if p.is_file() {
        return true;
    }
    std::fs::read_dir(p).is_ok_and(|mut d| d.next().is_some())
}

/// What the table says about a call of the tool `spec` with the effects
/// `fx` under the roots `roots`.
pub fn judge(policy: &Policy, spec: &ToolSpec, fx: &Effects, roots: &Roots) -> Verdict {
    let mut needs = Vec::new();
    let mut need = |kind: Kind, reason: String, path: Option<PathBuf>| {
        needs.push(Need {
            kind,
            rule: policy.rule(kind),
            reason,
            path,
        });
    };
    let mut refused = None;
    if spec.action == Action::StartRun {
        need(Kind::StartRun, format!("`{}` starts a run", spec.id), None);
    }
    if let Some(l) = &fx.limits {
        let mut missing = Vec::new();
        let mut fix = Vec::new();
        if l.max_wall_time.is_none() {
            missing.push("max_wall_time");
            let t = fx
                .estimate
                .as_ref()
                .and_then(|e| e.wall_time)
                .map_or(3600.0, |t| (2.0 * t).max(600.0).ceil());
            fix.push(format!("run.limits.max_wall_time={t}"));
        }
        if l.max_output.is_none() {
            missing.push("max_output");
            let b = fx
                .estimate
                .as_ref()
                .map_or(5e9, |e| (2.0 * e.output_bytes).max(1e8).ceil());
            fix.push(format!("run.limits.max_output={b:e}"));
        }
        if !missing.is_empty() {
            let mut d = Diagnostic::new(
                diag::MCP_LIMITS,
                format!(
                    "a run started through MCP needs [run.limits] {} (D198 (c)); add them to the case or pass them with `set`",
                    missing.join(" and ")
                ),
            );
            d.key = Some("run.limits".into());
            d.fix = Some(Fix { set: fix });
            refused = Some(d);
        }
    }
    if let Some(e) = &fx.estimate {
        let c = &policy.estimate;
        let mut over = Vec::new();
        if e.cells > c.max_cells {
            over.push(format!("{} cells > {}", e.cells, c.max_cells));
        }
        if e.host_memory > c.max_host_memory {
            over.push(format!(
                "host memory {} > {}",
                gb(e.host_memory),
                gb(c.max_host_memory)
            ));
        }
        if e.output_bytes > c.max_output {
            over.push(format!(
                "output {} > {}",
                gb(e.output_bytes),
                gb(c.max_output)
            ));
        }
        if let Some(t) = e.wall_time.filter(|t| *t > c.max_wall_time) {
            over.push(format!("wall time {t:.0} s > {:.0} s", c.max_wall_time));
        }
        if !over.is_empty() {
            need(
                Kind::ExceedEstimate,
                format!("the estimate exceeds the caps: {}", over.join(", ")),
                None,
            );
        }
    }
    for w in &fx.writes {
        if roots.check_write(w).is_err() {
            let dir = if w.is_file() {
                w.parent().map_or_else(|| w.clone(), Path::to_path_buf)
            } else {
                w.clone()
            };
            need(
                Kind::WriteOutsideRoot,
                format!("writes {}, outside the write roots", w.display()),
                Some(dir),
            );
        }
    }
    for w in &fx.writes {
        if holds_results(w) {
            need(
                Kind::OverwriteResults,
                format!(
                    "writes into {}, which holds earlier results or files",
                    w.display()
                ),
                None,
            );
        }
    }
    if fx.publishes {
        need(
            Kind::Publish,
            format!("`{}` publishes outside this computer", spec.id),
            None,
        );
    }
    needs.sort_by_key(|n| n.kind);
    needs.dedup_by_key(|n| n.kind);
    Verdict { needs, refused }
}

/// The refusal (`MCP003`) of the kinds a table refuses.
pub fn refusal(spec: &ToolSpec, needs: &[&Need], policy_file: &Path) -> Diagnostic {
    let what: Vec<String> = needs
        .iter()
        .map(|n| format!("{} ({})", n.kind.name(), n.reason))
        .collect();
    Diagnostic::new(
        diag::MCP_REFUSED,
        format!(
            "`{}` refused by the action table {}: {}; the person can change the rule there",
            spec.id,
            policy_file.display(),
            what.join("; ")
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_flow::case::Limits;

    fn spec(action: Action) -> ToolSpec {
        ToolSpec {
            id: "flow.run",
            version: 1,
            summary: "",
            category: crate::Category::Flow,
            action,
            streams: false,
        }
    }

    #[test]
    fn the_table_reads_and_writes_itself() {
        let mut p = Policy::default();
        assert_eq!(p.rule(Kind::StartRun), Setting::Ask);
        assert_eq!(p.rule(Kind::WriteOutsideRoot), Setting::Refuse);
        p.set(Kind::StartRun, Setting::Allow);
        p.estimate.max_cells = 10;
        let back = Policy::parse(&p.to_toml()).unwrap();
        assert_eq!(back, p);
        let e = Policy::parse("[actions]\nstart-runs = \"allow\"\n").unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MCP_POLICY);
        let e = Policy::parse("[actions]\nstart-run = \"yes\"\n").unwrap_err();
        assert_eq!(diag::of_error(&e).code, diag::MCP_POLICY);
    }

    #[test]
    fn a_run_needs_limits_and_meets_its_kinds() {
        let d = std::env::temp_dir().join(format!("foehn-approval-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        let inside = d.join("runs");
        std::fs::create_dir_all(inside.join("old")).unwrap();
        std::fs::write(inside.join("old/result.json"), "{}").unwrap();
        let roots = Roots {
            read: vec![d.clone()],
            write: vec![inside.clone()],
        };
        let p = Policy::default();
        let mut fx = Effects {
            writes: vec![inside.join("new")],
            limits: Some(Limits::default()),
            ..Default::default()
        };
        let v = judge(&p, &spec(Action::StartRun), &fx, &roots);
        let r = v.refused.as_ref().unwrap();
        assert_eq!(r.code, diag::MCP_LIMITS);
        assert_eq!(r.fix.as_ref().unwrap().set.len(), 2);
        assert_eq!(v.questions().count(), 1, "start-run asks");
        fx.limits = Some(Limits {
            max_wall_time: Some(60.0),
            max_output: Some(1e9),
            ..Default::default()
        });
        fx.writes = vec![inside.join("old"), d.join("elsewhere")];
        let v = judge(&p, &spec(Action::StartRun), &fx, &roots);
        assert!(v.refused.is_none());
        let kinds: Vec<Kind> = v.needs.iter().map(|n| n.kind).collect();
        assert_eq!(
            kinds,
            [
                Kind::StartRun,
                Kind::WriteOutsideRoot,
                Kind::OverwriteResults
            ]
        );
        assert_eq!(v.refusals().count(), 1);
        assert_eq!(v.extra_write_roots(), [d.join("elsewhere")]);
        let read = judge(&p, &spec(Action::Read), &Effects::default(), &roots);
        assert!(read.needs.is_empty());
    }
}
