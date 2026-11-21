//! `expression.eval`: an expression (D53) compiled and evaluated with the
//! solver's own compiler and functions.

use crate::{Action, Category, Ctx, Tool};
use anyhow::Result;
use foehncraft_flow::diag::{self, Diagnostic};
use foehncraft_flow::expr::{self, Env, Scalar, TableSection, Vars};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

fn all_variables() -> String {
    "xyzt".into()
}

/// Input of `expression.eval`.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EvalInput {
    /// The expression, as in a case file.
    pub source: String,
    /// `[constants]` and `[tables]` the expression may refer to, as TOML
    /// text in the form of a case file.
    #[serde(default)]
    pub env: String,
    /// The variables it may use: a string of `x`, `y`, `z`, `t`.
    #[serde(default = "all_variables")]
    pub allowed: String,
    /// Points `[x, y, z, t]` at which it is evaluated.
    #[serde(default)]
    pub points: Vec<[f64; 4]>,
    /// Directory of relative table files (default: the working directory).
    #[serde(default)]
    pub base: Option<PathBuf>,
}

impl EvalInput {
    /// An expression without constants, tables or points.
    pub fn new(source: impl Into<String>) -> Self {
        Self {
            source: source.into(),
            env: String::new(),
            allowed: all_variables(),
            points: Vec::new(),
            base: None,
        }
    }
}

/// Values of an expression.
#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct EvalOutput {
    /// The value at each point.
    pub values: Vec<f64>,
    /// The variables the expression uses (a string of `x`, `y`, `z`, `t`).
    pub uses: String,
}

/// The constants and tables an expression may refer to, as in a case file.
#[derive(Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct EnvText {
    #[serde(default)]
    constants: BTreeMap<String, Scalar>,
    #[serde(default)]
    tables: BTreeMap<String, TableSection>,
}

fn variables(v: &str) -> Result<Vars> {
    let mut vars = Vars::NONE;
    for c in v.chars() {
        vars = vars.with(match c {
            'x' => Vars::X,
            'y' => Vars::Y,
            'z' => Vars::Z,
            't' => Vars::T,
            _ => {
                return Err(diag::err(Diagnostic::new(
                    diag::EXPR_VARIABLE,
                    format!("unknown variable `{c}` (x, y, z, t)"),
                )));
            }
        });
    }
    Ok(vars)
}

/// `expression.eval`.
pub struct Eval;

impl Tool for Eval {
    const ID: &'static str = "expression.eval";
    const VERSION: u32 = 1;
    const SUMMARY: &'static str =
        "Compile an expression (D53) with constants and tables and evaluate it at points.";
    const CATEGORY: Category = Category::Expression;
    const ACTION: Action = Action::Read;
    type Input = EvalInput;
    type Output = EvalOutput;

    fn call(a: EvalInput, ctx: &mut Ctx<'_>) -> Result<EvalOutput> {
        let mut e: EnvText = toml::from_str(&a.env).map_err(|e| {
            diag::err(Diagnostic::new(
                diag::CASE_UNKNOWN_KEY,
                format!("constants and tables: {e}"),
            ))
        })?;
        let base = a.base.clone().unwrap_or_else(|| PathBuf::from("."));
        for (name, t) in &mut e.tables {
            if let Some(f) = &t.file {
                let f = if f.is_relative() {
                    base.join(f)
                } else {
                    f.clone()
                };
                ctx.roots.check_read(&f, "the table file")?;
                t.file = Some(f);
            }
            t.load(name)?;
        }
        let env = Env::new(&e.constants, &e.tables)?;
        let p = expr::compile(&a.source, &env, variables(&a.allowed)?, "expression")?;
        let prog = expr::Value::Program(Arc::new(p.clone()));
        let values: Vec<f64> = a
            .points
            .iter()
            .map(|q| prog.at([q[0], q[1], q[2]], q[3]))
            .collect();
        let uses: String = [
            ('x', Vars::X),
            ('y', Vars::Y),
            ('z', Vars::Z),
            ('t', Vars::T),
        ]
        .iter()
        .filter(|(_, v)| p.uses.contains(*v))
        .map(|(c, _)| *c)
        .collect();
        for (q, v) in a.points.iter().zip(&values) {
            ctx.line(&format!("({}, {}, {}, {}) -> {v}", q[0], q[1], q[2], q[3]));
        }
        Ok(EvalOutput { values, uses })
    }
}
