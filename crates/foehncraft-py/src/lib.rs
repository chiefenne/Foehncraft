//! The native module `foehncraft._native` of the Python package (D61).
//!
//! A thin layer over the tool registry of `foehncraft-tools` (D196):
//! `call(id, input)` runs any tool by id, and the functions below are
//! one-line wrappers of their tools, so Python reaches exactly what the
//! command line reaches. Check and run a case file in process, see every row of the run record while the
//! run goes on (monitors and stop criteria at iteration or step level),
//! evaluate expressions (D53) with the solver's own compiler and
//! functions, list the material database (D102) and return the case-file
//! schema (D77). Everything crosses the boundary as JSON text, the same
//! versioned formats the command line writes (D80), so the Python side
//! and `foehn --json` see the same data.
//!
//! The run releases the GIL; the callbacks take it back between outer
//! iterations or time steps, never per cell, so GPU speed and
//! reproducibility are kept (D53, D31).

use foehncraft_flow::diag;
use foehncraft_flow::engine::{BuildInfo, Observer, StopSignal};
use foehncraft_tools::{self as tools, Ctx, NoInput};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

create_exception!(
    _native,
    FoehncraftError,
    PyException,
    "An error of the solver with its diagnostic code (D78): args = (code, message)."
);

/// The error of the solver as a Python exception with its code.
fn py_err(e: anyhow::Error) -> PyErr {
    let d = diag::of_error(&e);
    FoehncraftError::new_err((d.code.to_string(), format!("{e:#}")))
}

/// Build and host for the provenance (D85).
fn build_info() -> BuildInfo {
    let mut features = vec!["python".to_string()];
    if cfg!(feature = "hdf5") {
        features.push("hdf5".into());
    }
    BuildInfo {
        commit: option_env!("FOEHNCRAFT_GIT_COMMIT")
            .unwrap_or("unknown")
            .into(),
        features,
        platform: format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH),
        cpu: String::new(),
        gpu: None,
    }
}

/// The version of the solver the module was built from.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// The context of a call: the build of this module, no roots, no text.
fn ctx() -> Ctx<'static> {
    Ctx::new(build_info())
}

/// The output of tool `T` as JSON text.
fn typed<T: tools::Tool>(input: T::Input, ctx: &mut Ctx<'_>) -> PyResult<String> {
    let (o, _) = tools::run::<T>(input, ctx).map_err(py_err)?;
    serde_json::to_string(&o).map_err(|e| py_err(e.into()))
}

/// The callbacks of a call: `observer(row)` gets every row of the run
/// record as JSON text and stops the run with a true value; `log(line)`
/// gets the human text. An exception in either stops the run and is kept
/// in the returned slot; Ctrl-C stops the run the same way.
fn callbacks(
    ctx: &mut Ctx<'static>,
    observer: Option<PyObject>,
    log: Option<PyObject>,
) -> Arc<Mutex<Option<PyErr>>> {
    let failed: Arc<Mutex<Option<PyErr>>> = Arc::new(Mutex::new(None));
    if let Some(cb) = observer {
        let stop = ctx.interrupts.clone();
        let failed = failed.clone();
        ctx.observer = Some(Observer::new(move |row: &serde_json::Value| {
            Python::with_gil(|py| {
                let text = row.to_string();
                match cb.call1(py, (text,)).and_then(|r| r.is_truthy(py)) {
                    Ok(true) => stop.request_stop(StopSignal::Api),
                    Ok(false) => {}
                    Err(e) => {
                        failed.lock().unwrap().get_or_insert(e);
                        stop.request_stop(StopSignal::Api);
                    }
                }
                if py.check_signals().is_err() {
                    stop.request_stop(StopSignal::Int);
                }
            })
        }));
    }
    let lg_failed = failed.clone();
    let lg_stop = ctx.interrupts.clone();
    let c = std::mem::replace(ctx, Ctx::new(BuildInfo::default()));
    *ctx = c.text(move |line: &str| {
        Python::with_gil(|py| {
            if let Some(cb) = &log
                && let Err(e) = cb.call1(py, (line,))
            {
                lg_failed.lock().unwrap().get_or_insert(e);
                lg_stop.request_stop(StopSignal::Api);
            }
            if py.check_signals().is_err() {
                lg_stop.request_stop(StopSignal::Int);
            }
        })
    });
    failed
}

/// Runs the tool `id` of the registry (D196) with the JSON text `input`
/// and returns the result envelope (`foehn.tool-result`) as JSON text.
/// Never raises for a failing tool: its diagnostic is in the envelope.
/// `observer` and `log` as for [`run`].
#[pyfunction]
#[pyo3(signature = (id, input = "{}", observer = None, log = None))]
fn call(
    py: Python<'_>,
    id: &str,
    input: &str,
    observer: Option<PyObject>,
    log: Option<PyObject>,
) -> PyResult<String> {
    let value: serde_json::Value = serde_json::from_str(input).map_err(|e| {
        FoehncraftError::new_err((
            diag::TOOL_INPUT.to_string(),
            format!("the input is not JSON: {e}"),
        ))
    })?;
    let mut c = ctx();
    let failed = callbacks(&mut c, observer, log);
    let id = id.to_string();
    let r = py.allow_threads(move || tools::call(&id, value, &mut c));
    if let Some(e) = failed.lock().unwrap().take() {
        return Err(e);
    }
    serde_json::to_string(&r).map_err(|e| py_err(e.into()))
}

/// The catalogue of every tool (`foehn.tool-catalogue`) as JSON text.
#[pyfunction]
fn commands() -> PyResult<String> {
    serde_json::to_string(&tools::catalogue()).map_err(|e| py_err(e.into()))
}

/// Checks a case file (`case.check`): validation, boundary conditions on
/// the first mesh and the estimate. Returns the check report
/// (`foehncraft.check`) as JSON.
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new(), out = None))]
fn check(
    py: Python<'_>,
    path: PathBuf,
    sets: Vec<String>,
    out: Option<PathBuf>,
) -> PyResult<String> {
    let input = tools::case::FlowInput {
        set: sets,
        out,
        ..tools::case::FlowInput::new(path)
    };
    py.allow_threads(|| typed::<tools::case::Check>(input, &mut ctx()))
}

/// Runs a case file in this process (`flow.run`) and returns the result
/// summary (`foehncraft.result`, also written to `result.json`) as JSON.
///
/// `observer(row)` gets every row of the run record as JSON text (the rows
/// of `run.jsonl`: header, iteration, step, balance, event, end); a true
/// return value stops the run as a stop request does (restart file, final
/// output, exit code 3). `log(line)` gets the progress lines. An exception
/// in a callback stops the run and is raised after it. Ctrl-C stops the
/// run the same way.
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new(), out = None, observer = None, log = None))]
fn run(
    py: Python<'_>,
    path: PathBuf,
    sets: Vec<String>,
    out: Option<PathBuf>,
    observer: Option<PyObject>,
    log: Option<PyObject>,
) -> PyResult<String> {
    let input = tools::case::FlowInput {
        set: sets,
        out,
        ..tools::case::FlowInput::new(path)
    };
    let mut c = ctx();
    let failed = callbacks(&mut c, observer, log);
    let res = py.allow_threads(move || tools::run::<tools::flow::Run>(input, &mut c));
    if let Some(e) = failed.lock().unwrap().take() {
        return Err(e);
    }
    let (summary, _) = res.map_err(py_err)?;
    serde_json::to_string(&summary).map_err(|e| py_err(e.into()))
}

/// Compiles the expression `source` with the `[constants]` and `[tables]`
/// of the TOML text `env` (relative table files from `base`) and the
/// variables it may use (a string of `x`, `y`, `z`, `t`), and evaluates it
/// at the `points` (x, y, z, t) (`expression.eval`). Returns the values and
/// the variables the expression uses. Errors carry the codes EXPR001–EXPR005.
#[pyfunction]
#[pyo3(signature = (source, env = "", allowed = "xyzt", points = Vec::new(), base = None))]
fn expression(
    source: &str,
    env: &str,
    allowed: &str,
    points: Vec<[f64; 4]>,
    base: Option<PathBuf>,
) -> PyResult<(Vec<f64>, String)> {
    let input = tools::expression::EvalInput {
        source: source.into(),
        env: env.into(),
        allowed: allowed.into(),
        points,
        base,
    };
    let (o, _) = tools::run::<tools::expression::Eval>(input, &mut ctx()).map_err(py_err)?;
    Ok((o.values, o.uses))
}

/// The material database (`foehncraft.materials`) as JSON, every entry
/// evaluated at `temperature` (K) and `pressure` (Pa); `name` picks one
/// (`materials.list`, `materials.show`). The files of `FOEHN_MATERIALS`
/// are searched as by the command line.
#[pyfunction]
#[pyo3(signature = (name = None, temperature = 293.15, pressure = 101325.0))]
fn material_database(name: Option<String>, temperature: f64, pressure: f64) -> PyResult<String> {
    if !(temperature > 0.0 && pressure > 0.0) {
        return Err(FoehncraftError::new_err((
            "MAT005".to_string(),
            "temperature and pressure must be > 0".to_string(),
        )));
    }
    use tools::materials as m;
    match name {
        Some(name) => typed::<m::Show>(
            m::ShowInput {
                name,
                temperature,
                pressure,
                include: Vec::new(),
            },
            &mut ctx(),
        ),
        None => typed::<m::List>(
            m::ListInput {
                temperature,
                pressure,
                include: Vec::new(),
            },
            &mut ctx(),
        ),
    }
}

/// The JSON Schema of the flow case file (D77), as `foehn schema` writes
/// it (`case.schema`).
#[pyfunction]
fn case_schema() -> PyResult<String> {
    typed::<tools::case::Schema>(NoInput {}, &mut ctx())
}

/// The diagnostic codes (D78) with their severity and meaning, as a JSON
/// array (`diag.codes`).
#[pyfunction]
fn codes() -> PyResult<String> {
    let (o, _) = tools::run::<tools::diag_codes::Codes>(NoInput {}, &mut ctx()).map_err(py_err)?;
    serde_json::to_string(&o.codes).map_err(|e| py_err(e.into()))
}

/// Reads a case file the way a run does and returns it, with the defaults
/// filled in and the overrides applied, as TOML text (the
/// `case.resolved.toml` of a run; `case.resolve`).
#[pyfunction]
#[pyo3(signature = (path, sets = Vec::new()))]
fn resolved(path: PathBuf, sets: Vec<String>) -> PyResult<String> {
    let input = tools::case::FlowInput {
        set: sets,
        ..tools::case::FlowInput::new(path)
    };
    let (o, _) = tools::run::<tools::case::Resolve>(input, &mut ctx()).map_err(py_err)?;
    Ok(o.toml)
}

/// The native module.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("FoehncraftError", m.py().get_type::<FoehncraftError>())?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(call, m)?)?;
    m.add_function(wrap_pyfunction!(commands, m)?)?;
    m.add_function(wrap_pyfunction!(check, m)?)?;
    m.add_function(wrap_pyfunction!(run, m)?)?;
    m.add_function(wrap_pyfunction!(expression, m)?)?;
    m.add_function(wrap_pyfunction!(material_database, m)?)?;
    m.add_function(wrap_pyfunction!(case_schema, m)?)?;
    m.add_function(wrap_pyfunction!(codes, m)?)?;
    m.add_function(wrap_pyfunction!(resolved, m)?)?;
    Ok(())
}
