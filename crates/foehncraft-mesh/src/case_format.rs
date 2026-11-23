//! The format version of case files.
//!
//! A case file may state `case_format = "0.6"`: the Foehncraft version that
//! introduced the key layout it is written in. A file without the key is in
//! format 0.5, the layout before versions were recorded. The rules in
//! [`RULES`] say which keys were renamed or moved in which format, so that
//! an older file keeps running: [`migrate`] maps its keys when the file is
//! read (the log notes each of them) and [`upgrade_text`] rewrites the text
//! of the file in the current format with its comments kept
//! (`foehn case upgrade`). A file in a newer format than the code is
//! refused; a file in the current format that still uses an old key gets the
//! usual "unknown key" error with the nearest key as a hint.
//!
//! The rules live in one table, here, so that the history of the case file
//! is in one place instead of in every file. `foehn schema` lists them.

use serde::Serialize;
use toml_edit::{DocumentMut, Item, Table};

/// The key of the format version.
pub const KEY: &str = "case_format";

/// The format written by this version of the code.
pub const CURRENT: &str = "0.7";

/// The format of a file without `case_format`.
pub const UNVERSIONED: &str = "0.5";

/// A key that was renamed or moved (dotted paths from the top of the file).
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Rename {
    /// The format that introduced the new name.
    pub since: &'static str,
    /// The old path.
    pub from: &'static str,
    /// The new path.
    pub to: &'static str,
}

const fn r(since: &'static str, from: &'static str, to: &'static str) -> Rename {
    Rename { since, from, to }
}

/// Every rename since format 0.5, oldest first.
pub const RULES: &[Rename] = &[
    // "boundary" instead of "patch" in every user-facing name
    r("0.6", "mesh.patches", "mesh.boundaries"),
    r(
        "0.6",
        "output.coefficients.patches",
        "output.coefficients.boundaries",
    ),
    r(
        "0.6",
        "verification.benchmark_patch",
        "verification.benchmark_boundary",
    ),
    // one table per equation (the "solution flags" of the log)
    r("0.6", "solver.scheme", "solver.momentum.scheme"),
    r("0.6", "solver.limiter", "solver.momentum.limiter"),
    r("0.6", "solver.limiter_k", "solver.momentum.limiter_k"),
    r(
        "0.6",
        "solver.central_blending",
        "solver.momentum.central_blending",
    ),
    r(
        "0.6",
        "solver.central_bounding",
        "solver.momentum.central_bounding",
    ),
    r(
        "0.6",
        "solver.relaxation_velocity",
        "solver.momentum.relaxation",
    ),
    r(
        "0.6",
        "solver.momentum_tolerance",
        "solver.momentum.tolerance",
    ),
    r(
        "0.6",
        "solver.momentum_max_iter",
        "solver.momentum.max_iter",
    ),
    r(
        "0.6",
        "solver.relaxation_pressure",
        "solver.pressure.relaxation",
    ),
    r(
        "0.6",
        "solver.pressure_tolerance",
        "solver.pressure.tolerance",
    ),
    r(
        "0.6",
        "solver.pressure_max_iter",
        "solver.pressure.max_iter",
    ),
    r(
        "0.6",
        "solver.pressure_preconditioner",
        "solver.pressure.preconditioner",
    ),
    r("0.6", "solver.amg", "solver.pressure.amg"),
];

/// A default that changed: a file of an older format that does not set the
/// key meant the old value, so reading it puts the old value in (the log
/// notes it) and `foehn case upgrade` writes it into the file. Every
/// existing case keeps its meaning and its results.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct ChangedDefault {
    /// The format that introduced the new default.
    pub since: &'static str,
    /// The key (dotted path from the top of the file).
    pub key: &'static str,
    /// The value of older files (a TOML string).
    pub old: &'static str,
    /// What the change means, for the note.
    pub why: &'static str,
}

/// Every changed default since format 0.5, oldest first. They apply to
/// flow cases only (see [`is_flow_case`]).
pub const DEFAULTS: &[ChangedDefault] = &[ChangedDefault {
    since: "0.7",
    key: "fluid.pressure_values",
    old: "gauge",
    why: "pressure values are absolute by default from format 0.7 (D239); \
          this file's values are gauge values",
}];

/// True for the table of a flow case: not a scalar case (`[problem]`), not
/// a Poisson case (`[linear]`) and not a mesh case (only `[mesh]`,
/// `[partition]`, `[check]`, `[output]`).
pub fn is_flow_case(has_key: impl Fn(&str) -> bool, keys: &[&str]) -> bool {
    if has_key("problem") || has_key("linear") {
        return false;
    }
    keys.iter()
        .any(|k| !matches!(*k, KEY | "mesh" | "partition" | "check" | "output"))
}

fn default_note(file: &str, d: &ChangedDefault) -> String {
    format!(
        "case format {file}: {} (`{} = \"{}\"`; `foehn case upgrade` writes it into the file)",
        d.why, d.key, d.old
    )
}

/// `major.minor` of a format string.
fn parse(v: &str) -> Option<(u32, u32)> {
    let (a, b) = v.trim().split_once('.')?;
    Some((a.parse().ok()?, b.parse().ok()?))
}

fn older(a: &str, b: &str) -> bool {
    matches!((parse(a), parse(b)), (Some(x), Some(y)) if x < y)
}

/// The rules that apply to a file of format `file`.
fn pending(file: &str) -> impl Iterator<Item = &'static Rename> + '_ {
    RULES.iter().filter(move |r| older(file, r.since))
}

fn check_version(file: &str) -> Result<(), String> {
    if parse(file).is_none() {
        return Err(format!(
            "{KEY} = \"{file}\" is not a format version (like \"{CURRENT}\")"
        ));
    }
    if older(CURRENT, file) {
        return Err(format!(
            "the case file is in format {file}, this Foehncraft reads formats up to {CURRENT}: \
             use a newer Foehncraft"
        ));
    }
    Ok(())
}

fn note(file: &str, r: &Rename) -> String {
    format!(
        "case format {file}: the key `{}` is now `{}` (`foehn case upgrade` rewrites the file)",
        r.from, r.to
    )
}

// ---------------------------------------------------------------------------
// the table in memory (reading a case)

fn take(t: &mut toml::Table, path: &str) -> Option<toml::Value> {
    let (head, rest) = match path.split_once('.') {
        Some((h, r)) => (h, Some(r)),
        None => (path, None),
    };
    match rest {
        None => t.remove(head),
        Some(rest) => {
            let toml::Value::Table(sub) = t.get_mut(head)? else {
                return None;
            };
            let v = take(sub, rest);
            if v.is_some() && sub.is_empty() {
                t.remove(head);
            }
            v
        }
    }
}

fn put(t: &mut toml::Table, path: &str, v: toml::Value) -> Result<(), String> {
    match path.split_once('.') {
        None => {
            t.insert(path.to_string(), v);
            Ok(())
        }
        Some((head, rest)) => {
            let e = t
                .entry(head.to_string())
                .or_insert_with(|| toml::Value::Table(toml::Table::new()));
            match e {
                toml::Value::Table(sub) => put(sub, rest, v),
                _ => Err(format!("`{head}` is a value, not a table")),
            }
        }
    }
}

fn has(t: &toml::Table, path: &str) -> bool {
    match path.split_once('.') {
        None => t.contains_key(path),
        Some((h, r)) => match t.get(h) {
            Some(toml::Value::Table(s)) => has(s, r),
            _ => false,
        },
    }
}

/// A rule applies to a case of the kind it was made for: the `solver.*`
/// renames are those of the flow case; the scalar case (`[problem]`) and
/// the Poisson case (`[linear]`) have their own `[solver]` keys.
fn applies(r: &Rename, has_key: impl Fn(&str) -> bool) -> bool {
    !(r.from.starts_with("solver.") && (has_key("problem") || has_key("linear")))
}

/// Reads the format of a case table, maps its old keys to the current ones
/// and states the current format in it. Returns one note per mapped key.
pub fn migrate(t: &mut toml::Table) -> Result<Vec<String>, String> {
    let file = match t.get(KEY) {
        None => UNVERSIONED.to_string(),
        Some(toml::Value::String(s)) => s.clone(),
        Some(_) => return Err(format!("{KEY} must be a string like \"{CURRENT}\"")),
    };
    check_version(&file)?;
    let mut notes = Vec::new();
    for r in pending(&file) {
        if !applies(r, |k| t.contains_key(k)) || !has(t, r.from) {
            continue;
        }
        if has(t, r.to) {
            return Err(format!(
                "both `{}` and `{}` are set: `{}` is the old name of `{}`, keep one",
                r.from, r.to, r.from, r.to
            ));
        }
        if let Some(v) = take(t, r.from) {
            put(t, r.to, v)?;
            notes.push(note(&file, r));
        }
    }
    let keys: Vec<String> = t.keys().cloned().collect();
    let keys: Vec<&str> = keys.iter().map(String::as_str).collect();
    if is_flow_case(|k| t.contains_key(k), &keys) {
        for d in DEFAULTS.iter().filter(|d| older(&file, d.since)) {
            if !has(t, d.key) {
                put(t, d.key, toml::Value::String(d.old.to_string()))?;
                notes.push(default_note(&file, d));
            }
        }
    }
    t.insert(KEY.to_string(), toml::Value::String(CURRENT.to_string()));
    Ok(notes)
}

// ---------------------------------------------------------------------------
// the text (foehn case upgrade)

fn take_item(t: &mut Table, path: &str) -> Option<Item> {
    match path.split_once('.') {
        None => t.remove(path),
        Some((head, rest)) => {
            let sub = t.get_mut(head)?.as_table_mut()?;
            let v = take_item(sub, rest);
            if v.is_some() && sub.is_empty() {
                t.remove(head);
            }
            v
        }
    }
}

fn put_item(t: &mut Table, path: &str, v: Item) -> Result<(), String> {
    match path.split_once('.') {
        None => {
            t.insert(path, v);
            Ok(())
        }
        Some((head, rest)) => {
            if !t.contains_key(head) {
                let mut sub = Table::new();
                sub.set_implicit(true);
                t.insert(head, Item::Table(sub));
            }
            match t.get_mut(head).and_then(Item::as_table_mut) {
                Some(sub) => put_item(sub, rest, v),
                None => Err(format!("`{head}` is a value, not a table")),
            }
        }
    }
}

fn has_item(t: &Table, path: &str) -> bool {
    match path.split_once('.') {
        None => t.contains_key(path),
        Some((h, r)) => t
            .get(h)
            .and_then(Item::as_table)
            .is_some_and(|s| has_item(s, r)),
    }
}

/// Puts a changed default into the text: into its table when the file has
/// one, else into a new table at the end.
fn put_default(root: &mut Table, path: &str, v: &str) -> Result<(), String> {
    let (head, key) = path.rsplit_once('.').unwrap_or(("", path));
    if head.is_empty() {
        root.insert(key, toml_edit::value(v));
        return Ok(());
    }
    put_item(root, &format!("{head}.{key}"), toml_edit::value(v))?;
    // a table created here is explicit, so that it is written
    let mut t = &mut *root;
    for part in head.split('.') {
        let Some(sub) = t.get_mut(part).and_then(Item::as_table_mut) else {
            break;
        };
        sub.set_implicit(false);
        t = sub;
    }
    Ok(())
}

/// The text of a case file in the current format, with its comments, and
/// one note per key that was renamed. The key `case_format` goes to the top.
pub fn upgrade_text(text: &str) -> Result<(String, Vec<String>), String> {
    let mut doc: DocumentMut = text.parse().map_err(|e| format!("not valid TOML: {e}"))?;
    let file = match doc.get(KEY) {
        None => UNVERSIONED.to_string(),
        Some(i) => i
            .as_str()
            .ok_or_else(|| format!("{KEY} must be a string like \"{CURRENT}\""))?
            .to_string(),
    };
    check_version(&file)?;
    let mut notes = Vec::new();
    for r in pending(&file) {
        let root = doc.as_table_mut();
        if !applies(r, |k| root.contains_key(k)) || !has_item(root, r.from) {
            continue;
        }
        if has_item(root, r.to) {
            return Err(format!(
                "both `{}` and `{}` are set: keep one",
                r.from, r.to
            ));
        }
        if let Some(v) = take_item(root, r.from) {
            put_item(root, r.to, v)?;
            notes.push(note(&file, r).replace(" (`foehn case upgrade` rewrites the file)", ""));
        }
    }
    {
        let root = doc.as_table_mut();
        let keys: Vec<String> = root.iter().map(|(k, _)| k.to_string()).collect();
        let keys: Vec<&str> = keys.iter().map(String::as_str).collect();
        if is_flow_case(|k| root.contains_key(k), &keys) {
            for d in DEFAULTS.iter().filter(|d| older(&file, d.since)) {
                if !has_item(root, d.key) {
                    put_default(root, d.key, d.old)?;
                    notes.push(format!(
                        "case format {file}: {} (`{} = \"{}\"` written)",
                        d.why, d.key, d.old
                    ));
                }
            }
        }
    }
    let mut out = doc.to_string();
    if file == CURRENT && notes.is_empty() {
        return Ok((out, notes));
    }
    // the format line goes to the top, after the leading comment block
    out = out
        .lines()
        .filter(|l| {
            !l.trim_start().starts_with(&format!("{KEY} "))
                && !l.trim_start().starts_with(&format!("{KEY}="))
        })
        .collect::<Vec<_>>()
        .join("\n")
        + if out.ends_with('\n') { "\n" } else { "" };
    let mut lines: Vec<&str> = out.lines().collect();
    let first_code = lines
        .iter()
        .position(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .unwrap_or(lines.len());
    let stated = format!("{KEY} = \"{CURRENT}\"");
    let mut at = first_code;
    // before the first key or table, after the header comment block
    while at > 0 && lines[at - 1].trim().is_empty() {
        at -= 1;
    }
    lines.insert(at, "");
    lines.insert(at + 1, &stated);
    if at + 2 < lines.len()
        && !lines[at + 2].trim().is_empty()
        && !lines[at + 2].trim().starts_with('[')
    {
        // keep a blank line before the next key
    }
    let mut s = lines.join("\n");
    s.push('\n');
    Ok((s, notes))
}

/// Reads a case of type `T` from text: the keys of an older format are
/// mapped first (the notes say which). The error of a file that needs no
/// mapping carries the line of the TOML parser.
pub fn read<T: serde::de::DeserializeOwned>(text: &str) -> anyhow::Result<(T, Vec<String>)> {
    let mut t: toml::Table = toml::from_str(text)?;
    let notes = migrate(&mut t).map_err(anyhow::Error::msg)?;
    let v = if notes.is_empty() {
        toml::from_str(text)?
    } else {
        T::deserialize(toml::Value::Table(t))?
    };
    Ok((v, notes))
}

/// The renames as JSON for `foehn schema`.
pub fn history() -> serde_json::Value {
    serde_json::json!({
        "current": CURRENT,
        "unversioned": UNVERSIONED,
        "renames": RULES,
        "changed_defaults": DEFAULTS,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn old_keys_are_mapped_and_noted() {
        let mut t: toml::Table =
            toml::from_str("[solver]\nrelaxation_velocity = 0.5\nmax_iterations = 10\n").unwrap();
        let notes = migrate(&mut t).unwrap();
        assert_eq!(notes.len(), 1);
        assert_eq!(t["solver"]["momentum"]["relaxation"].as_float(), Some(0.5));
        assert_eq!(t[KEY].as_str(), Some(CURRENT));
        assert!(
            t["solver"]
                .as_table()
                .unwrap()
                .get("relaxation_velocity")
                .is_none()
        );
    }

    #[test]
    fn a_current_file_keeps_its_old_keys_for_the_unknown_key_error() {
        let mut t: toml::Table =
            toml::from_str("case_format = \"0.7\"\n[solver]\nrelaxation_velocity = 0.5\n").unwrap();
        assert!(migrate(&mut t).unwrap().is_empty());
        assert!(
            t["solver"]
                .as_table()
                .unwrap()
                .contains_key("relaxation_velocity")
        );
    }

    #[test]
    fn both_names_are_an_error_and_newer_formats_are_refused() {
        let mut t: toml::Table = toml::from_str(
            "[solver]\nrelaxation_velocity = 0.5\n[solver.momentum]\nrelaxation = 0.6\n",
        )
        .unwrap();
        assert!(migrate(&mut t).is_err());
        let mut t: toml::Table = toml::from_str("case_format = \"9.0\"\n").unwrap();
        assert!(migrate(&mut t).unwrap_err().contains("newer"));
    }

    #[test]
    fn the_text_keeps_its_comments() {
        let text = "# A case.\n# More.\n\n[mesh]\nfile = \"a.msh\"\n\n[mesh.patches]  # rename\nwall = \"w\"\n\n[solver]\nrelaxation_velocity = 0.5  # alpha\nlog_every = 5\n";
        let (out, notes) = upgrade_text(text).unwrap();
        // two renames and the changed default of D239
        assert_eq!(notes.len(), 3);
        assert!(out.starts_with("# A case.\n# More.\n"));
        assert!(out.contains("case_format = \"0.7\""));
        assert!(out.contains("pressure_values = \"gauge\""), "{out}");
        assert!(out.contains("# alpha"));
        assert!(out.contains("[mesh.boundaries]"));
        assert!(out.contains("[solver.momentum]"));
        // the result is current: nothing more to do
        let (again, n2) = upgrade_text(&out).unwrap();
        assert!(n2.is_empty());
        assert_eq!(again, out);
        let t: toml::Table = toml::from_str(&out).unwrap();
        assert_eq!(t["solver"]["momentum"]["relaxation"].as_float(), Some(0.5));
    }

    #[test]
    fn solver_renames_leave_scalar_and_poisson_cases_alone() {
        for text in [
            "[problem]\n[solver]\nscheme = \"upwind\"\nrelaxation = 0.9\n",
            "[linear]\n[solver]\nscheme = \"upwind\"\n",
        ] {
            let (out, _) = upgrade_text(text).unwrap();
            assert!(!out.contains("[solver.momentum]"), "{out}");
            let mut t: toml::Table = toml::from_str(text).unwrap();
            migrate(&mut t).unwrap();
            assert!(t["solver"].as_table().unwrap().contains_key("scheme"));
        }
        let (out, notes) = upgrade_text("[fluid]\n[solver]\nscheme = \"upwind\"\n").unwrap();
        assert!(out.contains("momentum") && !notes.is_empty(), "{out}");
    }

    #[test]
    fn a_changed_default_keeps_the_meaning_of_older_files() {
        // an older flow case means gauge pressures: reading puts it in
        let mut t: toml::Table = toml::from_str(
            "case_format = \"0.6\"\n[boundary.out]\ntype = \"static-pressure-outlet\"\n",
        )
        .unwrap();
        let notes = migrate(&mut t).unwrap();
        assert_eq!(notes.len(), 1, "{notes:?}");
        assert_eq!(t["fluid"]["pressure_values"].as_str(), Some("gauge"));
        // a file that states it keeps its value; a current file gets nothing
        let mut t: toml::Table =
            toml::from_str("[fluid]\npressure_values = \"absolute\"\n").unwrap();
        assert!(migrate(&mut t).unwrap().is_empty());
        assert_eq!(t["fluid"]["pressure_values"].as_str(), Some("absolute"));
        let mut t: toml::Table =
            toml::from_str("case_format = \"0.7\"\n[fluid]\ndensity = 1.0\n").unwrap();
        assert!(migrate(&mut t).unwrap().is_empty());
        assert!(t["fluid"].get("pressure_values").is_none());
        // mesh, scalar and Poisson cases have no pressure values
        for text in ["[mesh]\nfile = \"a.msh\"\n", "[problem]\n", "[linear]\n"] {
            let mut t: toml::Table = toml::from_str(text).unwrap();
            migrate(&mut t).unwrap();
            assert!(t.get("fluid").is_none(), "{text}");
        }
        // the text: into the existing [fluid] table, comments kept
        let (out, _) =
            upgrade_text("# c\n[fluid]  # air\ndensity = 1.2\n\n[boundary.out]\ntype = \"wall\"\n")
                .unwrap();
        let t: toml::Table = toml::from_str(&out).unwrap();
        assert_eq!(t["fluid"]["pressure_values"].as_str(), Some("gauge"));
        assert!(out.contains("# air"));
        let (again, n2) = upgrade_text(&out).unwrap();
        assert!(n2.is_empty());
        assert_eq!(again, out);
    }
}
