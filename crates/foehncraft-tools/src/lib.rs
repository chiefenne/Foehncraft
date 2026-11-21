//! The tool registry of Foehncraft (D196).
//!
//! Every function of the engine is a **tool**: a Rust type implementing
//! [`Tool`], with a stable dotted id (`case.check`), a contract version,
//! typed input and output (their JSON Schemas come from the same structs,
//! so they cannot drift), a category and an action kind (D198). The
//! [`registry`] is an explicit list in catalogue order; [`call`] runs a
//! tool by id from JSON and returns the result envelope
//! (`foehn.tool-result`). The command line, the Python package and later
//! the MCP server (L18) and the desktop app call only this crate, so no
//! function exists in one front end only and an agent learns every tool
//! from one catalogue (`foehn commands`).
//!
//! A tool writes its human text (D79) through [`Ctx::line`]; the front end
//! decides where it goes (standard output, standard error, a Python
//! callback, nowhere). Errors are diagnostics with the codes of D78.

use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic};
use foehncraft_flow::engine::{BuildInfo, Interrupts, Observer, Roots};
use schemars::JsonSchema;
use serde::Serialize;
use serde::de::DeserializeOwned;
use serde_json::{Value, json};

pub mod case;
pub mod expression;
pub mod flow;
pub mod markdown;
pub mod materials;
pub mod mesh;
pub mod post;
pub mod system;

/// Format of the result envelope.
pub const RESULT_FORMAT: &str = "foehn.tool-result";
/// Format of the catalogue.
pub const CATALOGUE_FORMAT: &str = "foehn.tool-catalogue";
/// Version of the envelope and the catalogue formats.
pub const FORMAT_VERSION: u32 = 1;

/// The area a tool belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Category {
    /// Build, platform, the catalogue itself.
    System,
    /// Case files: schema, check, resolve, estimate, upgrade.
    Case,
    /// Diagnostic codes (D78).
    Diag,
    /// The material database (D94).
    Materials,
    /// Expressions (D53).
    Expression,
    /// Meshes.
    Mesh,
    /// Flow runs.
    Flow,
    /// Results of earlier runs.
    Result,
    /// Post-processing.
    Post,
    /// Verification and solver studies.
    Study,
}

impl Category {
    /// Lower-case name.
    pub fn name(self) -> &'static str {
        match self {
            Self::System => "system",
            Self::Case => "case",
            Self::Diag => "diag",
            Self::Materials => "materials",
            Self::Expression => "expression",
            Self::Mesh => "mesh",
            Self::Flow => "flow",
            Self::Result => "result",
            Self::Post => "post",
            Self::Study => "study",
        }
    }
}

/// What a tool may do (D198 (b)): the approvals of an agent host follow it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Action {
    /// Reads files and computes; writes nothing.
    Read,
    /// Writes or rewrites files.
    WriteFiles,
    /// Starts a run (writes an output directory, may take long).
    StartRun,
}

impl Action {
    /// Kebab-case name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Read => "read",
            Self::WriteFiles => "write-files",
            Self::StartRun => "start-run",
        }
    }
}

/// What the catalogue says about a tool.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ToolSpec {
    /// Stable id, dotted lowercase `area.verb`.
    pub id: &'static str,
    /// Version of the input and output contract; a breaking change raises it.
    pub version: u32,
    /// One line.
    pub summary: &'static str,
    /// Area.
    pub category: Category,
    /// Action kind (D198 (b)).
    pub action: Action,
    /// The tool streams events (rows of the run record) while it runs.
    pub streams: bool,
}

/// A tool: one function of the engine with typed input and output.
pub trait Tool {
    /// Stable id (`case.check`).
    const ID: &'static str;
    /// Contract version.
    const VERSION: u32;
    /// One line for the catalogue.
    const SUMMARY: &'static str;
    /// Area.
    const CATEGORY: Category;
    /// Action kind.
    const ACTION: Action;
    /// Streams events.
    const STREAMS: bool = false;
    /// The input (read from JSON; unknown fields are refused).
    type Input: DeserializeOwned + JsonSchema;
    /// The output.
    type Output: Serialize + JsonSchema;

    /// The catalogue entry.
    fn spec() -> ToolSpec {
        ToolSpec {
            id: Self::ID,
            version: Self::VERSION,
            summary: Self::SUMMARY,
            category: Self::CATEGORY,
            action: Self::ACTION,
            streams: Self::STREAMS,
        }
    }

    /// Runs the tool.
    fn call(input: Self::Input, ctx: &mut Ctx<'_>) -> Result<Self::Output>;

    /// The exit code of a finished call (D64): 0 finished, 2 failed checks,
    /// 3 stopped; 1 when the output reports an error (a case that cannot
    /// run).
    fn exit_code(_output: &Self::Output) -> i32 {
        0
    }
}

/// What every call gets from its front end: roots (D82), interrupts and the
/// observer of a run, the build information (D85) and the sink of the
/// human text (D79). The journal of D197 hooks in here later.
pub struct Ctx<'a> {
    /// File confinement.
    pub roots: Roots,
    /// Stop and restart requests.
    pub interrupts: Interrupts,
    /// Sees every row of the run record of a run.
    pub observer: Option<Observer>,
    /// Build and host.
    pub build: BuildInfo,
    text: Box<dyn FnMut(&str) + Send + 'a>,
}

impl<'a> Ctx<'a> {
    /// A context with the build information, no roots and no text.
    pub fn new(build: BuildInfo) -> Self {
        Self {
            roots: Roots::default(),
            interrupts: Interrupts::new(),
            observer: None,
            build,
            text: Box::new(|_| {}),
        }
    }

    /// With the read and write roots.
    pub fn roots(mut self, roots: Roots) -> Self {
        self.roots = roots;
        self
    }

    /// With the sink of the human text.
    pub fn text(mut self, f: impl FnMut(&str) + Send + 'a) -> Self {
        self.text = Box::new(f);
        self
    }

    /// One line of human text.
    pub fn line(&mut self, l: &str) {
        (self.text)(l);
    }
}

/// An input without fields.
#[derive(Debug, Clone, Default, serde::Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NoInput {}

/// A tool of the registry with its types erased.
pub struct Entry {
    /// Catalogue entry.
    pub spec: ToolSpec,
    input_schema: fn() -> Value,
    output_schema: fn() -> Value,
    run: fn(Value, &mut Ctx<'_>) -> Result<(Value, i32)>,
}

impl Entry {
    /// JSON Schema of the input.
    pub fn input_schema(&self) -> Value {
        (self.input_schema)()
    }

    /// JSON Schema of the output.
    pub fn output_schema(&self) -> Value {
        (self.output_schema)()
    }
}

fn schema_of<T: JsonSchema>() -> Value {
    schemars::SchemaGenerator::default()
        .into_root_schema_for::<T>()
        .to_value()
}

/// Reads the input of tool `T` from JSON (`null` is no fields).
pub fn input_of<T: Tool>(input: Value) -> Result<T::Input> {
    let input = if input.is_null() { json!({}) } else { input };
    serde_json::from_value(input).map_err(|e| {
        diag::err(Diagnostic::new(
            diag::TOOL_INPUT,
            format!("input of `{}`: {e}", T::ID),
        ))
    })
}

fn run_erased<T: Tool>(input: Value, ctx: &mut Ctx<'_>) -> Result<(Value, i32)> {
    let out = T::call(input_of::<T>(input)?, ctx)?;
    let code = T::exit_code(&out);
    Ok((serde_json::to_value(&out)?, code))
}

/// The registry entry of tool `T`.
pub fn entry<T: Tool>() -> Entry {
    Entry {
        spec: T::spec(),
        input_schema: schema_of::<T::Input>,
        output_schema: schema_of::<T::Output>,
        run: run_erased::<T>,
    }
}

/// Every tool, in catalogue order. An explicit list: no macros, no
/// link-time collection, so the order is fixed and `docs/tools.md` can be
/// compared with it.
pub fn registry() -> Vec<Entry> {
    vec![
        entry::<system::Info>(),
        entry::<system::Backends>(),
        entry::<system::Commands>(),
        entry::<case::Schema>(),
        entry::<case::Check>(),
        entry::<case::Resolve>(),
        entry::<case::EstimateTool>(),
        entry::<case::Upgrade>(),
        entry::<diag_codes::Codes>(),
        entry::<materials::List>(),
        entry::<materials::Show>(),
        entry::<expression::Eval>(),
        entry::<mesh::Process>(),
        entry::<flow::Run>(),
        entry::<flow::Read>(),
        entry::<post::Derive>(),
    ]
}

/// `diag.codes`.
pub mod diag_codes {
    use super::*;
    use foehncraft_flow::diag::CodeInfo;

    /// The diagnostic codes.
    #[derive(Debug, Serialize, JsonSchema)]
    pub struct CodesOutput {
        /// `foehncraft.codes`.
        pub format: &'static str,
        /// Format version of the engine.
        pub format_version: &'static str,
        /// Every code, in registry order.
        pub codes: Vec<CodeInfo>,
    }

    /// `diag.codes`: the diagnostic codes (D78).
    pub struct Codes;

    impl Tool for Codes {
        const ID: &'static str = "diag.codes";
        const VERSION: u32 = 1;
        const SUMMARY: &'static str =
            "The diagnostic codes (D78) with their usual severity and meaning.";
        const CATEGORY: Category = Category::Diag;
        const ACTION: Action = Action::Read;
        type Input = NoInput;
        type Output = CodesOutput;

        fn call(_: NoInput, ctx: &mut Ctx<'_>) -> Result<CodesOutput> {
            for c in diag::CODES {
                ctx.line(&format!(
                    "{:<8} {:<7} {}",
                    c.code,
                    c.severity.name(),
                    c.title
                ));
            }
            Ok(CodesOutput {
                format: "foehncraft.codes",
                format_version: foehncraft_flow::engine::FORMAT_VERSION,
                codes: diag::CODES.to_vec(),
            })
        }
    }
}

/// One tool in the catalogue.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ToolDescription {
    /// What the catalogue says about it.
    #[serde(flatten)]
    pub spec: ToolSpec,
    /// JSON Schema of the input.
    pub input_schema: Value,
    /// JSON Schema of the output (the `result` of the envelope).
    pub output_schema: Value,
}

/// The catalogue (`foehn commands`, `system.commands`).
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct Catalogue {
    /// [`CATALOGUE_FORMAT`].
    pub format: &'static str,
    /// [`FORMAT_VERSION`].
    pub version: u32,
    /// Version of Foehncraft.
    pub foehncraft_version: &'static str,
    /// Every tool, in catalogue order.
    pub tools: Vec<ToolDescription>,
}

/// The catalogue of every tool.
pub fn catalogue() -> Catalogue {
    Catalogue {
        format: CATALOGUE_FORMAT,
        version: FORMAT_VERSION,
        foehncraft_version: foehncraft_core::VERSION,
        tools: registry()
            .into_iter()
            .map(|e| ToolDescription {
                input_schema: e.input_schema(),
                output_schema: e.output_schema(),
                spec: e.spec,
            })
            .collect(),
    }
}

/// The result envelope of [`call`].
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ToolResult {
    /// [`RESULT_FORMAT`].
    pub format: &'static str,
    /// [`FORMAT_VERSION`].
    pub version: u32,
    /// Tool id.
    pub tool: String,
    /// Contract version of the tool (0 for an unknown tool).
    pub tool_version: u32,
    /// `finished`, `error`, `failed-checks` or `stopped` (D64).
    pub status: &'static str,
    /// Exit code (D64).
    pub exit_code: i32,
    /// The output of the tool (`null` on an error).
    pub result: Value,
    /// The diagnostic of a call that failed.
    pub diagnostics: Vec<Diagnostic>,
    /// Known limitations of the result (D199; empty until it fills them).
    pub limitations: Vec<Value>,
}

/// The status of an exit code (D64).
pub fn status_of(code: i32) -> &'static str {
    match code {
        0 => "finished",
        2 => "failed-checks",
        3 => "stopped",
        _ => "error",
    }
}

fn distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut prev = row[0];
        row[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let cur = row[j + 1];
            row[j + 1] = (prev + usize::from(ca != *cb)).min(row[j] + 1).min(cur + 1);
            prev = cur;
        }
    }
    row[b.len()]
}

/// Runs the tool `id` with the JSON `input` and returns the envelope.
/// Never fails: an unknown id, a bad input or an error of the tool is a
/// diagnostic in the envelope.
pub fn call(id: &str, input: Value, ctx: &mut Ctx<'_>) -> ToolResult {
    let reg = registry();
    let Some(e) = reg.iter().find(|e| e.spec.id == id) else {
        let near = reg
            .iter()
            .map(|e| e.spec.id)
            .min_by_key(|n| distance(id, n));
        let hint = near.map_or(String::new(), |n| format!("; did you mean `{n}`?"));
        return failed(
            id,
            0,
            Diagnostic::new(
                diag::TOOL_UNKNOWN,
                format!("unknown tool `{id}`{hint} (`foehn commands` lists the tools)"),
            ),
        );
    };
    match (e.run)(input, ctx) {
        Ok((result, code)) => ToolResult {
            format: RESULT_FORMAT,
            version: FORMAT_VERSION,
            tool: id.into(),
            tool_version: e.spec.version,
            status: status_of(code),
            exit_code: code,
            result,
            diagnostics: Vec::new(),
            limitations: Vec::new(),
        },
        Err(err) => {
            let mut d = diag::of_error(&err);
            d.message = format!("{err:#}");
            failed(id, e.spec.version, d)
        }
    }
}

fn failed(id: &str, version: u32, d: Diagnostic) -> ToolResult {
    ToolResult {
        format: RESULT_FORMAT,
        version: FORMAT_VERSION,
        tool: id.into(),
        tool_version: version,
        status: "error",
        exit_code: 1,
        result: Value::Null,
        diagnostics: vec![d],
        limitations: Vec::new(),
    }
}

/// Runs tool `T` with a typed input (what the clap commands do).
pub fn run<T: Tool>(input: T::Input, ctx: &mut Ctx<'_>) -> Result<(T::Output, i32)> {
    let out = T::call(input, ctx)?;
    let code = T::exit_code(&out);
    Ok((out, code))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_are_unique_dotted_lowercase() {
        let mut seen = std::collections::BTreeSet::new();
        for e in registry() {
            let id = e.spec.id;
            assert!(seen.insert(id), "{id} twice");
            let parts: Vec<&str> = id.split('.').collect();
            assert_eq!(parts.len(), 2, "{id}");
            assert!(
                parts
                    .iter()
                    .all(|p| !p.is_empty() && p.chars().all(|c| c.is_ascii_lowercase())),
                "{id}"
            );
            assert_eq!(
                parts[0],
                e.spec.category.name(),
                "{id}: area is the category"
            );
            assert!(e.spec.version >= 1, "{id}");
        }
    }

    #[test]
    fn every_schema_is_an_object_schema() {
        for e in registry() {
            let i = e.input_schema();
            assert_eq!(i["type"], "object", "{} input", e.spec.id);
            assert!(e.output_schema().is_object(), "{} output", e.spec.id);
        }
    }

    #[test]
    fn unknown_tool_and_bad_input_are_diagnostics() {
        let mut ctx = Ctx::new(BuildInfo::default());
        let r = call("case.chek", Value::Null, &mut ctx);
        assert_eq!(r.exit_code, 1);
        assert_eq!(r.diagnostics[0].code, diag::TOOL_UNKNOWN);
        assert!(r.diagnostics[0].message.contains("`case.check`"));
        let r = call("diag.codes", json!({"x": 1}), &mut ctx);
        assert_eq!(r.diagnostics[0].code, diag::TOOL_INPUT);
        assert_eq!(r.status, "error");
    }

    #[test]
    fn codes_through_call() {
        let mut ctx = Ctx::new(BuildInfo::default());
        let r = call("diag.codes", Value::Null, &mut ctx);
        assert_eq!(r.status, "finished");
        assert_eq!(
            r.result["codes"].as_array().map(Vec::len),
            Some(diag::CODES.len())
        );
    }

    #[test]
    fn tools_md_is_current() {
        let doc = include_str!("../../../docs/tools.md").replace("\r\n", "\n");
        assert!(
            doc == markdown::catalogue(),
            "docs/tools.md is out of date: regenerate it with\n  \
             cargo run --release -p foehncraft-cli -- commands --markdown > docs/tools.md"
        );
    }
}
