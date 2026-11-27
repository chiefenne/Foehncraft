//! The manifest `series.json` of the files a transient run writes (D268):
//! one entry per step with the field file and the parcel file of that
//! step, so that a viewer (MurmelCraft's time player) steps through a run
//! without guessing from file names.
//!
//! ```json
//! {"format": "foehncraft.series", "format_version": 1, "case": "spray_duct",
//!  "entries": [{"step": 50, "time": 0.25,
//!               "files": {"fields": "fields_000050.vtkhdf",
//!                         "parcels": "spray/parcels_000050.vtkhdf"}}]}
//! ```
//!
//! Paths are relative to the folder of `series.json` (the run's output
//! folder); entries are sorted by step; an entry holds the files written at
//! that step (fields, parcels or both). Files written once for the whole
//! run are under `static` (`{"devices": "spray/devices.vtu"}`, D345). The file is rewritten after each
//! output through a temporary file and a rename, so a running case can be
//! followed. A restart keeps the entries before its first step. Readers
//! ignore keys they do not know (D80).

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Format name of `series.json`.
pub const FORMAT: &str = "foehncraft.series";

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct Entry {
    step: usize,
    time: f64,
    files: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct File {
    format: String,
    format_version: u32,
    case: String,
    /// Files written once for the whole run (D345), for example
    /// `"devices": "spray/devices.vtu"`.
    #[serde(default, rename = "static", skip_serializing_if = "BTreeMap::is_empty")]
    fixed: BTreeMap<String, String>,
    entries: Vec<Entry>,
}

/// The manifest of a run.
pub struct Series {
    dir: PathBuf,
    case: String,
    fixed: BTreeMap<String, String>,
    entries: BTreeMap<usize, Entry>,
}

impl Series {
    /// The manifest of the run writing into `dir`; with `continue_before`
    /// (a restart) the entries of the existing file before that step are
    /// kept.
    pub fn new(dir: &Path, case: &str, continue_before: Option<usize>) -> Self {
        let mut entries = BTreeMap::new();
        if let Some(first) = continue_before
            && let Ok(t) = std::fs::read_to_string(dir.join("series.json"))
            && let Ok(f) = serde_json::from_str::<File>(&t)
        {
            for e in f.entries {
                if e.step < first {
                    entries.insert(e.step, e);
                }
            }
        }
        Self {
            dir: dir.to_path_buf(),
            case: case.into(),
            fixed: BTreeMap::new(),
            entries,
        }
    }

    /// Adds a file of `kind` written once for the whole run (the key
    /// `static` of the manifest, D345) and rewrites the manifest.
    pub fn set_static(&mut self, kind: &str, path: &Path) -> Result<()> {
        let rel = self.relative(path);
        self.fixed.insert(kind.into(), rel);
        self.write()
    }

    fn relative(&self, path: &Path) -> String {
        let rel = path.strip_prefix(&self.dir).unwrap_or(path);
        rel.components()
            .map(|c| c.as_os_str().to_string_lossy())
            .collect::<Vec<_>>()
            .join("/")
    }

    /// Adds file `path` (inside the output folder) of `kind` (`fields`,
    /// `parcels`) written at `step`, `time`, and rewrites the manifest.
    pub fn add(&mut self, step: usize, time: f64, kind: &str, path: &Path) -> Result<()> {
        let rel = self.relative(path);
        let e = self.entries.entry(step).or_insert_with(|| Entry {
            step,
            time,
            files: BTreeMap::new(),
        });
        e.time = time;
        e.files.insert(kind.into(), rel);
        self.write()
    }

    fn write(&self) -> Result<()> {
        let f = File {
            format: FORMAT.into(),
            format_version: 1,
            case: self.case.clone(),
            fixed: self.fixed.clone(),
            entries: self.entries.values().cloned().collect(),
        };
        let tmp = self.dir.join("series.json.tmp");
        std::fs::write(&tmp, serde_json::to_vec_pretty(&f)?)?;
        std::fs::rename(&tmp, self.dir.join("series.json"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entries_are_sorted_and_merged_by_step() {
        let d = std::env::temp_dir().join(format!("series-test-{}", std::process::id()));
        std::fs::create_dir_all(d.join("spray")).unwrap();
        let mut s = Series::new(&d, "c", None);
        s.add(20, 0.2, "fields", &d.join("fields_000020.vtkhdf"))
            .unwrap();
        s.add(10, 0.1, "parcels", &d.join("spray/parcels_000010.vtkhdf"))
            .unwrap();
        s.add(20, 0.2, "parcels", &d.join("spray/parcels_000020.vtkhdf"))
            .unwrap();
        let f: File =
            serde_json::from_str(&std::fs::read_to_string(d.join("series.json")).unwrap()).unwrap();
        assert_eq!(f.format, FORMAT);
        assert_eq!(
            f.entries.iter().map(|e| e.step).collect::<Vec<_>>(),
            vec![10, 20]
        );
        assert_eq!(f.entries[1].files["parcels"], "spray/parcels_000020.vtkhdf");
        assert_eq!(f.entries[1].files["fields"], "fields_000020.vtkhdf");
        // a restart from step 15 keeps step 10 only
        let s2 = Series::new(&d, "c", Some(15));
        assert_eq!(s2.entries.len(), 1);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// D345: the devices of the Spray G device file are written once, with
    /// the standard tip, the 8 hole axes and the cone of the one injector,
    /// and listed under `static`.
    #[test]
    fn devices_are_written_and_listed_as_static() {
        let d = std::env::temp_dir().join(format!("series-devices-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let file = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples/cases/flow/spray_g_device_injectors.toml");
        let sp = foehncraft_particles::SpraySection::from_file(&file).unwrap();
        let p = crate::spray::write_devices(&sp, &d).unwrap().unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        let types = text
            .split("Name=\"types\" format=\"ascii\">")
            .nth(1)
            .unwrap()
            .split("</DataArray>")
            .next()
            .unwrap()
            .split_whitespace()
            .map(|t| t.parse::<u8>().unwrap())
            .collect::<Vec<_>>();
        let n = |t: u8| types.iter().filter(|&&x| x == t).count();
        let tris = sp.devices[0].surface().unwrap().len();
        assert!(tris > 1000);
        assert_eq!((n(5), n(3), n(4)), (tris, 8, 1));
        let mut s = Series::new(&d, "c", None);
        s.set_static("devices", &p).unwrap();
        let f: File =
            serde_json::from_str(&std::fs::read_to_string(d.join("series.json")).unwrap()).unwrap();
        assert_eq!(f.fixed["devices"], "spray/devices.vtu");
        let _ = std::fs::remove_dir_all(&d);
    }
}
