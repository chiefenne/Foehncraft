//! The console and `flow.log`: frames, tables and numbers (D43).
//!
//! One formatter feeds the terminal and the log file, so the file shows
//! what the terminal showed, without colour. A log that circulates alone
//! must say what was set: the header is a series of framed blocks (the
//! version and host, the files, the run control, the mesh, the fluid, the
//! solution settings and the boundary conditions), the iterations or time
//! steps are a table with a header row, and the end of the run is a
//! summary, the warnings and a closing frame with the exit code.
//!
//! Every frame is a box as wide as its widest line, centred on a page as
//! wide as the banner. [`Style::Asterisk`] (the default) draws the box
//! with asterisks; [`Style::Box`] draws light box lines. Nothing here changes a result: only text is produced.

use serde::{Deserialize, Serialize};
use std::fmt::Display;

/// Width of the page on which the frames are centred: the widest frame,
/// the banner with the name of the code.
pub const PAGE: usize = 133;

/// Spaces between the side of a frame and its text.
const PAD: usize = 3;

/// Width of the label column of a frame.
const LABEL: usize = 35;

/// Longest value on one line of a frame; a longer one is wrapped.
const VALUE: usize = 60;

/// How frames and rules are drawn (`[log] style`).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Style {
    /// A box of asterisks.
    #[default]
    Asterisk,
    /// Light box lines (`┌─┐│└┘`).
    Box,
}

/// How much of the run the console and `flow.log` show (`[log] detail`,
/// D214). `run.jsonl` keeps every row at every level.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Detail {
    /// Everything: the header blocks, every table row, the box and the
    /// outer-iteration table of every logged time step, the summaries.
    #[default]
    Full,
    /// As `full`, but a logged time step shows its box and its table row,
    /// not the table of its outer iterations.
    Steps,
    /// The header blocks, the messages, the summaries and the closing frame:
    /// no table rows, no step boxes, no balance lines during the run.
    Summary,
}

impl Detail {
    /// Whether the table rows (and the balance lines) are shown.
    pub fn rows(self) -> bool {
        self != Detail::Summary
    }

    /// Whether the table of the outer iterations of a time step is shown.
    pub fn outer(self) -> bool {
        self == Detail::Full
    }
}

/// `[log]`: how the console and `flow.log` look.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct LogSection {
    /// Frame style: `"asterisk"` (rules of asterisks) or `"box"` (light
    /// box lines).
    pub style: Style,
    /// Repeat the header row of the iteration or time-step table every N
    /// rows (0: only once).
    pub header_every: usize,
    /// Columns of the iteration or time-step table, as
    /// `<kind>.<target>.<quantity>` (see `foehn schema`): residuals, linear
    /// iterations, values at probes or in cells, extremes with their cells,
    /// boundary mass flows, time quantities. Without it the default set
    /// is shown. May change during a run (`[[schedule]]`, `run.set`, a
    /// re-read): the header row is printed again (D214).
    pub columns: Option<Vec<String>>,
    /// Significant digits of the numbers in the table (1 to 8).
    pub digits: usize,
    /// How much the log shows: `"full"` (every row, every outer iteration
    /// of a logged time step), `"steps"` (a logged time step shows its box
    /// and its row, not its outer iterations) or `"summary"` (the header
    /// blocks, the summaries and the closing frame). `run.jsonl` keeps
    /// everything at every level.
    pub detail: Detail,
}

impl Default for LogSection {
    fn default() -> Self {
        Self {
            style: Style::Asterisk,
            header_every: 40,
            columns: None,
            digits: 3,
            detail: Detail::Full,
        }
    }
}

/// Whether the locale of the terminal is UTF-8: the first of `LC_ALL`,
/// `LC_CTYPE` and `LANG` that is set names `UTF-8` (or `utf8`). With none
/// set the terminal is taken to be UTF-8.
pub fn utf8_locale() -> bool {
    utf8_of(
        ["LC_ALL", "LC_CTYPE", "LANG"]
            .iter()
            .find_map(|k| std::env::var(k).ok().filter(|v| !v.is_empty()))
            .as_deref(),
    )
}

/// [`utf8_locale`] for a locale string.
pub fn utf8_of(locale: Option<&str>) -> bool {
    locale.is_none_or(|v| {
        let v = v.to_ascii_lowercase();
        v.contains("utf-8") || v.contains("utf8")
    })
}

/// A number as `9.88E-01`: `decimals` digits after the point, an exponent
/// with a sign and two digits, `NaN` and `Inf` as words.
pub fn sci(v: f64, decimals: usize) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf" } else { "-Inf" }.into();
    }
    let s = format!("{v:.decimals$e}");
    match s.split_once('e') {
        Some((m, e)) => {
            let e: i32 = e.parse().unwrap_or(0);
            format!("{m}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
        }
        None => s,
    }
}

/// [`sci`] with three significant digits (the residual tables).
pub fn sci3(v: f64) -> String {
    sci(v, 2)
}

/// [`sci`] with five significant digits (the values of the summaries).
pub fn sci5(v: f64) -> String {
    sci(v, 4)
}

/// A row of a [`Frame`].
#[derive(Debug, Clone)]
pub enum Row {
    /// A label and its value.
    Kv(String, String),
    /// A line of text.
    Text(String),
    /// A horizontal rule inside the frame.
    Rule,
}

/// A framed block: a title and rows.
#[derive(Debug, Clone, Default)]
pub struct Frame {
    title: String,
    /// Lines drawn, centred, in place of the title line (the banner).
    art: Vec<String>,
    /// A summary box: blank rows above and below, the title as a heading
    /// with a colon, dashed rules.
    summary: bool,
    /// No rule between the title lines and the rows.
    plain: bool,
    rows: Vec<Row>,
}

impl Frame {
    /// A frame with a title.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            art: Vec::new(),
            summary: false,
            plain: false,
            rows: Vec::new(),
        }
    }

    /// Draws the frame as a summary: blank rows above and below, the
    /// title as a heading with a colon and dashed rules between groups.
    #[must_use]
    pub fn summary(mut self) -> Self {
        self.summary = true;
        self
    }

    /// Draws the frame without the rule between the title and the rows
    /// (the box of a time step).
    #[must_use]
    pub fn plain(mut self) -> Self {
        self.plain = true;
        self
    }

    /// Draws `art`, centred between blank lines, in place of the title
    /// line; the title stays the name of the frame in the run record.
    #[must_use]
    pub fn art(mut self, art: &[&str]) -> Self {
        self.art = art.iter().map(|l| l.trim_end().to_string()).collect();
        self
    }

    /// The title as printed: a summary heading ends with a colon.
    fn heading_text(&self) -> String {
        if self.summary && !self.title.contains(':') {
            format!("{}:", self.title)
        } else {
            self.title.clone()
        }
    }

    /// A label and a value.
    #[must_use]
    pub fn kv(mut self, label: &str, value: impl Display) -> Self {
        self.rows.push(Row::Kv(label.into(), value.to_string()));
        self
    }

    /// A line of text.
    #[must_use]
    pub fn text(mut self, t: impl Into<String>) -> Self {
        self.rows.push(Row::Text(t.into()));
        self
    }

    /// A rule between two groups of rows.
    #[must_use]
    pub fn rule(mut self) -> Self {
        self.rows.push(Row::Rule);
        self
    }

    /// Adds rows in place.
    pub fn push_kv(&mut self, label: &str, value: impl Display) {
        self.rows.push(Row::Kv(label.into(), value.to_string()));
    }

    /// Adds a text row in place.
    pub fn push_text(&mut self, t: impl Into<String>) {
        self.rows.push(Row::Text(t.into()));
    }

    /// Adds a rule in place.
    pub fn push_rule(&mut self) {
        self.rows.push(Row::Rule);
    }

    /// Whether the frame has no rows.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// The title.
    pub fn title(&self) -> &str {
        &self.title
    }

    /// The rows as (label, value) pairs; a text row has an empty label.
    pub fn pairs(&self) -> Vec<(String, String)> {
        self.rows
            .iter()
            .filter_map(|r| match r {
                Row::Kv(k, v) => Some((k.clone(), v.clone())),
                Row::Text(t) => Some((String::new(), t.clone())),
                Row::Rule => None,
            })
            .collect()
    }

    /// The frame as a row of the run record: its title and the rows with
    /// label, value and unit (`section`, run record 1.1).
    pub fn record(&self) -> serde_json::Value {
        let rows: Vec<serde_json::Value> = self
            .pairs()
            .into_iter()
            .map(|(label, value)| {
                let (value, unit) = split_unit(&value);
                serde_json::json!({"label": label, "value": value, "unit": unit})
            })
            .collect();
        serde_json::json!({"type": "section", "title": self.title, "rows": rows})
    }

    fn body(&self, r: &Row) -> Vec<Option<String>> {
        match r {
            Row::Kv(k, v) => {
                let mut out = Vec::new();
                let wrapped: Vec<String> = v.lines().flat_map(|l| wrap(l, VALUE)).collect();
                let mut lines = wrapped.iter().map(String::as_str);
                let first = lines.next().unwrap_or("");
                if k.chars().count() >= LABEL {
                    // a label longer than its column stands on its own line
                    out.push(Some(k.clone()));
                    out.push(Some(format!("{:<LABEL$}{first}", "")));
                } else {
                    out.push(Some(format!("{k:<LABEL$}{first}")));
                }
                for extra in lines {
                    out.push(Some(format!("{:<LABEL$}{extra}", "")));
                }
                out
            }
            Row::Text(t) => vec![Some(t.clone())],
            Row::Rule => vec![None],
        }
    }

    /// The lines of the frame, without the blank line before it: a box of
    /// the width of its widest line, centred on the page.
    pub fn render(&self, style: Style) -> Vec<String> {
        let rows: Vec<Option<String>> = self.rows.iter().flat_map(|r| self.body(r)).collect();
        let art_w = self
            .art
            .iter()
            .map(|l| l.chars().count())
            .max()
            .unwrap_or(0);
        let inner = rows
            .iter()
            .flatten()
            .map(|l| l.chars().count())
            .chain([art_w])
            .chain(
                self.art
                    .is_empty()
                    .then(|| self.heading_text().chars().count()),
            )
            .max()
            .unwrap_or(0);
        let (sep, top, bottom, side) = match style {
            Style::Asterisk => {
                let r = "*".repeat(inner + 2 * PAD + 2);
                (r.clone(), r.clone(), r, '*')
            }
            Style::Box => {
                let bar = "─".repeat(inner + 2 * PAD);
                (
                    format!("├{bar}┤"),
                    format!("┌{bar}┐"),
                    format!("└{bar}┘"),
                    '│',
                )
            }
        };
        let pad = " ".repeat(PAD);
        let line = |l: &str| {
            let n = l.chars().count();
            format!("{side}{pad}{l}{}{pad}{side}", " ".repeat(inner - n))
        };
        let mut out = vec![top];
        let dashes = match style {
            Style::Asterisk => format!(" {} {side}", "-".repeat(inner + 2 * PAD - 2)),
            Style::Box => sep.clone(),
        };
        if self.art.is_empty() {
            if self.summary {
                out.push(line(""));
                out.push(line(&self.heading_text()));
                out.push(line(""));
            } else {
                out.push(line(&self.title));
            }
        } else {
            let left = " ".repeat((inner - art_w) / 2);
            out.push(line(""));
            for a in &self.art {
                out.push(line(&format!("{left}{a}")));
            }
            out.push(line(""));
        }
        if !rows.is_empty() {
            if !self.summary && !self.plain {
                out.push(sep.clone());
            }
            for r in rows {
                match r {
                    Some(l) => out.push(line(&l)),
                    None if self.summary => out.push(dashes.clone()),
                    None => out.push(sep.clone()),
                }
            }
            if self.summary {
                out.push(line(""));
            }
        } else if self.summary {
            out.push(line(""));
        }
        out.push(bottom);
        let width = inner + 2 * PAD + 2;
        let indent = " ".repeat(PAGE.saturating_sub(width) / 2);
        out.into_iter().map(|l| format!("{indent}{l}")).collect()
    }
}

/// A column of a [`Table`].
#[derive(Debug, Clone)]
pub struct Col {
    /// Short name in the header row.
    pub name: String,
    /// Width of the column.
    pub width: usize,
}

impl Col {
    /// A column of numbers (`1.23E-04`).
    pub fn num(name: &str) -> Self {
        Self {
            name: name.into(),
            width: 8,
        }
    }

    /// A column of integers of the given width.
    pub fn int(name: &str, width: usize) -> Self {
        Self {
            name: name.into(),
            width,
        }
    }
}

/// A group of columns under one title.
#[derive(Debug, Clone)]
pub struct Group {
    /// Title above the columns.
    pub title: String,
    /// The columns.
    pub cols: Vec<Col>,
}

impl Group {
    /// A group.
    pub fn new(title: &str, cols: Vec<Col>) -> Self {
        Self {
            title: title.into(),
            cols,
        }
    }
}

/// A table of iterations or time steps: a counter column, then groups of
/// columns; the header repeats every `header_every` rows.
#[derive(Debug, Clone)]
pub struct Table {
    first: Col,
    groups: Vec<Group>,
    style: Style,
    header_every: usize,
    rows: usize,
}

impl Table {
    /// A table whose first column is `first` (`Iter`, `Step`).
    pub fn new(first: Col, groups: Vec<Group>, style: Style, header_every: usize) -> Self {
        Self {
            first,
            groups,
            style,
            header_every,
            rows: 0,
        }
    }

    /// Prints the header again before the next row.
    pub fn restart(&mut self) {
        self.rows = 0;
    }

    /// The header: the group titles with rules, then the column names.
    pub fn header(&self) -> Vec<String> {
        let rule = match self.style {
            Style::Asterisk => '-',
            Style::Box => '─',
        };
        let mut g = " ".repeat(self.first.width + 2);
        let mut h = format!("  {:>w$}", self.first.name, w = self.first.width);
        for grp in &self.groups {
            let w: usize = grp.cols.iter().map(|c| c.width + 1).sum::<usize>() + 1;
            let t = format!("{rule}{rule} {} ", grp.title);
            let mut t: String = t.chars().take(w).collect();
            while t.chars().count() < w {
                t.push(rule);
            }
            g.push(' ');
            g.push_str(&t);
            h.push(' ');
            for c in &grp.cols {
                h.push_str(&format!(" {:>w$}", c.name, w = c.width));
            }
            h.push(' ');
        }
        vec![g.trim_end().to_string(), h]
    }

    /// A row: the counter and the cells (as many as there are columns). The
    /// header comes first on the first row and every `header_every` rows.
    pub fn row(&mut self, first: impl Display, cells: &[String]) -> Vec<String> {
        let mut out = Vec::new();
        if self.rows == 0 || (self.header_every > 0 && self.rows.is_multiple_of(self.header_every))
        {
            if self.rows > 0 {
                out.push(String::new());
            }
            out.extend(self.header());
        }
        self.rows += 1;
        let mut l = format!("  {:>w$}", first.to_string(), w = self.first.width);
        let mut k = 0;
        for grp in &self.groups {
            l.push(' ');
            for c in &grp.cols {
                let v = cells.get(k).map_or("", String::as_str);
                l.push_str(&format!(" {:>w$}", v, w = c.width));
                k += 1;
            }
            l.push(' ');
        }
        out.push(l.trim_end().to_string());
        out
    }
}

/// `1.000E+00 [kg/m3]` as the value and its unit (empty when there is
/// none at the end).
fn split_unit(v: &str) -> (String, String) {
    if let Some(stripped) = v.strip_suffix(']')
        && let Some((val, unit)) = stripped.rsplit_once(" [")
        && !unit.contains(['[', ']', ' '])
    {
        return (val.trim_end().to_string(), unit.to_string());
    }
    (v.to_string(), String::new())
}

/// Breaks a line at spaces into lines of at most `width` characters (a
/// longer word stays whole).
pub fn wrap(text: &str, width: usize) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut cur = String::new();
    for w in text.split(' ') {
        if cur.is_empty() {
            cur = w.to_string();
        } else if cur.chars().count() + 1 + w.chars().count() <= width {
            cur.push(' ');
            cur.push_str(w);
        } else {
            out.push(std::mem::take(&mut cur));
            cur = w.to_string();
        }
    }
    out.push(cur);
    out
}

/// A message line ` %FOEHN-I-CASE, text`: the severity (`I`, `W`, `E`) and
/// the name of the source file that writes the message, in capitals. The
/// file is the caller's, so each module tags its own messages.
#[track_caller]
pub fn tagged(severity: &str, msg: &str) -> String {
    let t = match severity {
        "error" => 'E',
        "warning" => 'W',
        _ => 'I',
    };
    let file = std::panic::Location::caller().file();
    let name = std::path::Path::new(file)
        .file_stem()
        .map_or(String::new(), |s| s.to_string_lossy().to_uppercase());
    format!(" %FOEHN-{t}-{name}, {msg}")
}

/// A date and time `2026-11-16 17:02:11 UTC` from the ISO form.
pub fn pretty_date(iso: &str) -> String {
    iso.replace('T', " ").replace('Z', " UTC")
}

/// The name of the host, as far as the standard library and the system
/// say it.
pub fn hostname() -> String {
    for k in ["HOSTNAME", "HOST", "COMPUTERNAME"] {
        if let Ok(v) = std::env::var(k)
            && !v.trim().is_empty()
        {
            return v.trim().to_string();
        }
    }
    for f in ["/proc/sys/kernel/hostname", "/etc/hostname"] {
        if let Ok(v) = std::fs::read_to_string(f)
            && !v.trim().is_empty()
        {
            return v.trim().to_string();
        }
    }
    std::process::Command::new("hostname")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".into())
}

/// `1234567` as `1 234 567` (a thin grouping for counts).
pub fn count(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(' ');
        }
        out.push(ch);
    }
    out
}

/// A size in bytes as `271 kB` or `8 MB`.
pub fn size(bytes: usize) -> String {
    if bytes >= 1_000_000 {
        format!("{} MB", bytes / 1_000_000 + 1)
    } else {
        format!("{} kB", bytes / 1000 + 1)
    }
}

/// Seconds as `0:03:58` (hours, minutes, seconds).
pub fn clock(secs: f64) -> String {
    let s = secs.max(0.0).round() as u64;
    format!("{}:{:02}:{:02}", s / 3600, (s / 60) % 60, s % 60)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_lines_need_a_utf8_locale() {
        assert!(utf8_of(None));
        assert!(utf8_of(Some("en_US.UTF-8")));
        assert!(utf8_of(Some("C.utf8")));
        assert!(!utf8_of(Some("C")));
        assert!(!utf8_of(Some("POSIX")));
        assert!(!utf8_of(Some("de_AT.ISO-8859-1")));
    }

    #[test]
    fn detail_levels() {
        assert!(Detail::Full.rows() && Detail::Full.outer());
        assert!(Detail::Steps.rows() && !Detail::Steps.outer());
        assert!(!Detail::Summary.rows() && !Detail::Summary.outer());
        assert_eq!(Detail::default(), Detail::Full);
    }

    #[test]
    fn numbers_have_a_signed_two_digit_exponent() {
        assert_eq!(sci3(9.877e-1), "9.88E-01");
        assert_eq!(sci3(0.0), "0.00E+00");
        assert_eq!(sci3(1.5e12), "1.50E+12");
        assert_eq!(sci3(-3.54e3), "-3.54E+03");
        assert_eq!(sci3(f64::NAN), "NaN");
        assert_eq!(sci5(1.0 / 3.0), "3.3333E-01");
        assert_eq!(sci3(1.0e-100), "1.00E-100");
    }

    #[test]
    fn art_replaces_the_title_line_and_fits_the_page() {
        let f = Frame::new("Foehncraft")
            .art(&crate::logblocks::BANNER_ART)
            .kv("Version", "0.5.0");
        assert_eq!(f.title(), "Foehncraft");
        for s in [Style::Asterisk, Style::Box] {
            let l = f.render(s);
            assert!(l.iter().all(|x| x.chars().count() == PAGE), "{l:?}");
            assert!(!l.iter().any(|x| x.contains("Foehncraft")), "{l:?}");
            assert!(l.iter().any(|x| x.contains("/:::/__\\:::\\")), "{l:?}");
        }
    }

    #[test]
    fn frames_are_boxes_centred_on_the_page() {
        let f = Frame::new("Fluid properties")
            .kv("Density", "1.000E+00 [kg/m3]")
            .rule()
            .text("free text");
        for s in [Style::Asterisk, Style::Box] {
            let l = f.render(s);
            let w = l[0].chars().count();
            assert!(l.iter().all(|x| x.chars().count() == w), "{l:?}");
            let indent = l[0].len() - l[0].trim_start().len();
            assert_eq!(indent, (PAGE - (w - indent)) / 2, "{l:?}");
        }
        let l = f.render(Style::Asterisk);
        assert!(l[1].trim().starts_with("*   Fluid properties"), "{l:?}");
        assert!(l[1].trim_end().ends_with('*'), "{l:?}");
    }

    #[test]
    fn the_table_repeats_its_header() {
        let mut t = Table::new(
            Col::int("Iter", 6),
            vec![Group::new(
                "scaled residuals",
                vec![Col::num("ResU"), Col::num("ResM")],
            )],
            Style::Asterisk,
            2,
        );
        let cells = vec![sci3(1e-3), sci3(2e-4)];
        let mut all = Vec::new();
        for i in 1..=3 {
            all.extend(t.row(i, &cells));
        }
        let heads = all.iter().filter(|l| l.contains("Iter")).count();
        assert_eq!(heads, 2, "{all:#?}");
        assert!(
            all.iter().any(|l| l.ends_with("1.00E-03 2.00E-04")),
            "{all:#?}"
        );
    }
}
