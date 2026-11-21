//! Broken-input corpus of the mesh readers (D202): every reader gets
//! truncated files, damaged numbers, missing and repeated lines, flipped
//! bytes, garbage and the files of the other formats. Each input must end in
//! an error message (or in a mesh, when the damage happens to leave a valid
//! file), never in a panic, an abort or a hang. The case layer turns the
//! message into the diagnostic code `MESH003` (tested in foehncraft-flow,
//! `tests/broken_inputs.rs`).
//!
//! The corpus is generated from the valid fixtures of `tests/data`, so it
//! follows them; `tests/data/broken` holds hand-written inputs for the
//! failures found on the way (allocation from a declared count, indices out
//! of range, loops on a missing end marker).

use foehncraft_mesh::Geometry;
use foehncraft_mesh::format::{self, MeshFile, MeshFormat, ReadOptions};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::time::Duration;

/// Time one input may take (the largest fixture reads in milliseconds).
const LIMIT: Duration = Duration::from_secs(20);

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// What became of one input.
enum Outcome {
    Mesh,
    Error,
}

/// Runs `read` on its own thread: a panic or a run past [`LIMIT`] fails.
fn guarded(
    label: &str,
    read: impl FnOnce() -> Result<MeshFile, String> + Send + 'static,
) -> Outcome {
    let (tx, rx) = mpsc::channel();
    let label_owned = label.to_string();
    std::thread::Builder::new()
        .name(label_owned)
        .stack_size(64 << 20)
        .spawn(move || {
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| match read() {
                Ok(f) => {
                    // a mesh that a reader returns is checked like any other
                    if f.mesh.validate().is_ok() && f.mesh.n_cells > 0 {
                        let _ = Geometry::compute(&f.mesh);
                    }
                    Ok(())
                }
                Err(e) => Err(e),
            }));
            let _ = tx.send(r);
        })
        .expect("thread");
    match rx.recv_timeout(LIMIT) {
        Ok(Ok(Ok(()))) => Outcome::Mesh,
        Ok(Ok(Err(e))) => {
            assert!(!e.trim().is_empty(), "{label}: an empty error message");
            Outcome::Error
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

/// A small deterministic generator (no seed from the clock).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Positions spread evenly over `0..n`, at most `k` of them.
fn spread(n: usize, k: usize) -> Vec<usize> {
    if n <= k {
        return (0..n).collect();
    }
    (0..k).map(|i| i * n / k).collect()
}

/// Spans of the number tokens of a text, in order.
fn numbers(b: &[u8]) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let start = i;
        if b[i] == b'-' && i + 1 < b.len() && b[i + 1].is_ascii_digit() {
            i += 1;
        }
        if b[i].is_ascii_digit() && (start == 0 || !b[start - 1].is_ascii_alphanumeric()) {
            while i < b.len()
                && (b[i].is_ascii_digit() || matches!(b[i], b'.' | b'e' | b'E' | b'+' | b'-'))
            {
                i += 1;
            }
            v.push((start, i));
        } else {
            i = start + 1;
        }
    }
    v
}

/// Byte ranges of the lines (without their end).
fn lines(b: &[u8]) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    let mut s = 0;
    for (i, &c) in b.iter().enumerate() {
        if c == b'\n' {
            v.push((s, i));
            s = i + 1;
        }
    }
    if s < b.len() {
        v.push((s, b.len()));
    }
    v
}

fn splice(b: &[u8], r: (usize, usize), with: &[u8]) -> Vec<u8> {
    let mut v = Vec::with_capacity(b.len() + with.len());
    v.extend_from_slice(&b[..r.0]);
    v.extend_from_slice(with);
    v.extend_from_slice(&b[r.1..]);
    v
}

/// The damaged versions of one valid input.
fn mutations(b: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = vec![
        ("empty".into(), Vec::new()),
        ("one space".into(), b" ".to_vec()),
        (
            "first line only".into(),
            b[..lines(b).first().map_or(0, |l| l.1)].to_vec(),
        ),
    ];
    // truncated
    for i in 1..32 {
        let n = b.len() * i / 32;
        out.push((format!("truncated at {n}"), b[..n].to_vec()));
    }
    // numbers replaced: the first tokens (versions, counts), the first token
    // of a line (counts, tags, section sizes) and tokens spread over the file
    let nums = numbers(b);
    let ls = lines(b);
    let mut picks: Vec<usize> = (0..nums.len().min(24)).collect();
    let line_first: Vec<usize> = ls
        .iter()
        .filter_map(|&(s, e)| nums.iter().position(|&(a, _)| a >= s && a < e))
        .collect();
    let mut lf = line_first.clone();
    lf.dedup();
    picks.extend(spread(lf.len(), 48).into_iter().map(|i| lf[i]));
    picks.extend(spread(nums.len(), 24));
    picks.sort_unstable();
    picks.dedup();
    const BAD: [&str; 9] = [
        "-1",
        "0",
        "3",
        "4294967295",
        "18446744073709551616",
        "999999999999",
        "1e308",
        "nan",
        "x",
    ];
    for &p in &picks {
        for bad in BAD {
            out.push((
                format!(
                    "number {p} ({}) -> {bad}",
                    String::from_utf8_lossy(&b[nums[p].0..nums[p].1])
                ),
                splice(b, nums[p], bad.as_bytes()),
            ));
        }
    }
    // lines dropped and repeated: keyword lines (with a letter) and spread lines
    let mut lpicks: Vec<usize> = ls
        .iter()
        .enumerate()
        .filter(|&(_, &(s, e))| b[s..e].iter().any(|c| c.is_ascii_alphabetic()))
        .map(|(i, _)| i)
        .take(64)
        .collect();
    lpicks.extend(spread(ls.len(), 24));
    lpicks.sort_unstable();
    lpicks.dedup();
    for &i in &lpicks {
        let (s, e) = ls[i];
        let end = (e + 1).min(b.len());
        out.push((format!("line {i} dropped"), splice(b, (s, end), b"")));
        let mut twice = b[s..end].to_vec();
        if twice.last() != Some(&b'\n') {
            twice.push(b'\n');
        }
        twice.extend_from_slice(&b[s..end]);
        out.push((format!("line {i} twice"), splice(b, (s, end), &twice)));
    }
    // bytes flipped and words set to all ones (binary sizes, offsets, indices)
    let mut g = Lcg(0x5eed ^ b.len() as u64);
    for p in spread(b.len(), 48) {
        let mut v = b.to_vec();
        v[p] ^= 0xff;
        out.push((format!("byte {p} flipped"), v));
    }
    for _ in 0..48 {
        if b.len() < 8 {
            break;
        }
        let p = g.next() as usize % (b.len() - 4);
        let mut v = b.to_vec();
        v[p..p + 4].copy_from_slice(&[0xff; 4]);
        out.push((format!("bytes {p}..{} set to 0xff", p + 4), v));
    }
    // garbage
    let mut v = Vec::with_capacity(4096);
    for _ in 0..4096 {
        v.push(g.next() as u8);
    }
    out.push(("random bytes".into(), v));
    out.push(("zeros".into(), vec![0; 4096]));
    out.push(("header then zeros".into(), {
        let mut v = b[..b.len().min(64)].to_vec();
        v.extend(std::iter::repeat_n(0u8, 4096));
        v
    }));
    out
}

/// Damage that no reader may take for a mesh. (Other damage, a changed
/// coordinate say, can leave a valid file.)
const MUST_FAIL: &[&str] = &["empty", "one space", "random bytes", "zeros"];

/// The valid single-file fixtures and the reader that owns them.
const FILES: &[(&str, &str)] = &[
    ("box_tet_ascii.msh", "gmsh"),
    ("box_tet_binary.msh", "gmsh"),
    ("channel_mixed_ascii.msh", "gmsh"),
    ("channel_mixed_binary.msh", "gmsh"),
    ("pyramids_ascii.msh", "gmsh"),
    ("channel_2d.su2", "su2"),
    ("cubes_mixed_3d.su2", "su2"),
    ("multizone/rotor_stator.su2", "su2"),
    ("industry_msh/channel_hex.msh", "industry_msh"),
    ("industry_msh/dual_polyhedra.msh", "industry_msh"),
    ("vtk/channel_mixed_ascii.vtk", "vtk"),
    ("vtk/channel_mixed_binary.vtk", "vtk"),
    ("vtk/channel_mixed_ascii.vtu", "vtk"),
    ("vtk/channel_mixed_binary.vtu", "vtk"),
    ("vtk/channel_mixed_zlib.vtu", "vtk"),
];

fn reader(name: &str) -> &'static dyn MeshFormat {
    format::by_name(name).unwrap_or_else(|e| panic!("{e}"))
}

/// Reads every damaged version of `file` with its reader; returns how many
/// inputs ended in a mesh and how many in an error.
fn corpus_of(file: &str, name: &'static str, options: ReadOptions) -> (usize, usize) {
    let bytes = std::fs::read(data(file)).unwrap_or_else(|e| panic!("{file}: {e}"));
    let (mut meshes, mut errors) = (0, 0);
    for (what, m) in mutations(&bytes) {
        let label = format!("{name} reader, {file}, {what}");
        let must_fail = MUST_FAIL.contains(&what.as_str());
        match guarded(&label, move || reader(name).parse(&m, &options)) {
            Outcome::Mesh if must_fail => panic!("{label}: read as a mesh"),
            Outcome::Mesh => meshes += 1,
            Outcome::Error => errors += 1,
        }
    }
    (meshes, errors)
}

#[test]
fn gmsh_reader_survives_the_corpus() {
    for &(f, r) in FILES.iter().filter(|x| x.1 == "gmsh") {
        let (m, e) = corpus_of(f, r, ReadOptions::default());
        assert!(e > m, "{f}: {m} meshes, {e} errors: the damage is not seen");
    }
}

#[test]
fn su2_reader_survives_the_corpus() {
    for &(f, r) in FILES.iter().filter(|x| x.1 == "su2") {
        corpus_of(f, r, ReadOptions::default());
        corpus_of(
            f,
            r,
            ReadOptions {
                thickness: Some(0.1),
                merge_nodes: true,
            },
        );
    }
}

#[test]
fn industry_msh_reader_survives_the_corpus() {
    for &(f, r) in FILES.iter().filter(|x| x.1 == "industry_msh") {
        corpus_of(f, r, ReadOptions::default());
    }
}

#[test]
fn vtk_reader_survives_the_corpus() {
    for &(f, r) in FILES.iter().filter(|x| x.1 == "vtk") {
        corpus_of(f, r, ReadOptions::default());
    }
}

/// Every file reader is given the fixtures of the other formats.
#[test]
fn readers_refuse_the_other_formats() {
    for &(f, owner) in FILES {
        let bytes = std::fs::read(data(f)).unwrap();
        for r in format::FORMATS {
            let name = r.name();
            if name == owner || matches!(name, "openfoam" | "cgns") {
                continue;
            }
            let b = bytes.clone();
            let label = format!("{name} reader, {f}");
            let opts = ReadOptions::default();
            if let Outcome::Mesh = guarded(&label, move || reader(name).parse(&b, &opts)) {
                panic!("{label}: a file of another format was read as a mesh");
            }
        }
    }
}

/// Copies a polyMesh directory into a fresh directory of the target dir.
fn scratch(tag: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("broken-{tag}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let t = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_dir(&e.path(), &t);
        } else {
            std::fs::copy(e.path(), t).unwrap();
        }
    }
}

#[test]
fn openfoam_reader_survives_the_corpus() {
    for case in ["channel_ascii", "channel_binary", "dual_polyhedra"] {
        let src = data(&format!("openfoam/{case}"));
        let poly = Path::new("constant/polyMesh");
        let files: Vec<String> = std::fs::read_dir(src.join(poly))
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        for file in &files {
            let bytes = std::fs::read(src.join(poly).join(file)).unwrap();
            for (k, (what, m)) in mutations(&bytes).into_iter().enumerate() {
                let dir = scratch(&format!("openfoam-{case}-{file}-{k}"));
                copy_dir(&src, &dir);
                std::fs::write(dir.join(poly).join(file), &m).unwrap();
                let label = format!("openfoam reader, {case}/{file}, {what}");
                let d = dir.clone();
                guarded(&label, move || {
                    reader("openfoam").read_path(&d, &ReadOptions::default())
                });
                let _ = std::fs::remove_dir_all(&dir);
            }
            // the file missing
            let dir = scratch(&format!("openfoam-{case}-{file}-missing"));
            copy_dir(&src, &dir);
            std::fs::remove_file(dir.join(poly).join(file)).unwrap();
            let d = dir.clone();
            let label = format!("openfoam reader, {case} without {file}");
            let o = guarded(&label, move || {
                reader("openfoam").read_path(&d, &ReadOptions::default())
            });
            if file != "cellZones" && matches!(o, Outcome::Mesh) {
                panic!("{label}: read as a mesh");
            }
            let _ = std::fs::remove_dir_all(&dir);
        }
    }
}

/// CGNS files are HDF5 files: the damage is done to the bytes (the HDF5
/// library must refuse them) and, through HDF5, to the CGNS nodes.
#[cfg(feature = "hdf5")]
#[test]
fn cgns_reader_survives_the_corpus() {
    for file in [
        "cgns/box_hex.cgns",
        "cgns/box_tet.cgns",
        "cgns/dual_ngon.cgns",
        "multizone/rotor_stator.cgns",
    ] {
        let bytes = std::fs::read(data(file)).unwrap();
        let all = mutations(&bytes);
        // byte damage: every truncation and flip, a sample of the rest
        for (k, (what, m)) in all.into_iter().enumerate() {
            if !(what.starts_with("truncated")
                || what.starts_with("byte")
                || what.starts_with("bytes")
                || k % 16 == 0)
            {
                continue;
            }
            let dir = scratch(&format!("cgns-{}-{k}", file.replace('/', "-")));
            let p = dir.join("mesh.cgns");
            std::fs::write(&p, &m).unwrap();
            let label = format!("cgns reader, {file}, {what}");
            let q = p.clone();
            guarded(&label, move || {
                reader("cgns").read_path(&q, &ReadOptions::default())
            });
            let _ = std::fs::remove_dir_all(&dir);
        }
    }
}

/// Paths of the integer ` data` datasets of a CGNS file (the sizes,
/// ranges, connectivity and offsets of its nodes), in file order.
#[cfg(feature = "hdf5")]
fn cgns_int_nodes(g: &hdf5_metno::Group, at: &str, out: &mut Vec<String>) {
    use hdf5_metno::types::TypeDescriptor;
    for name in g.member_names().unwrap_or_default() {
        let path = format!("{at}/{name}");
        if name == " data" {
            if let Ok(d) = g.dataset(&name)
                && matches!(
                    d.dtype().and_then(|t| t.to_descriptor()),
                    Ok(TypeDescriptor::Integer(_))
                )
            {
                out.push(at.to_string());
            }
        } else if let Ok(c) = g.group(&name) {
            cgns_int_nodes(&c, &path, out);
        }
    }
}

/// A change of the `i`-th value `v` of a node.
#[cfg(feature = "hdf5")]
type Damage = fn(usize, i64) -> i64;

/// The CGNS nodes themselves damaged: every integer node (sizes, ranges,
/// connectivity, offsets, BC lists) set to values out of range.
#[cfg(feature = "hdf5")]
#[test]
fn cgns_reader_survives_damaged_nodes() {
    for file in [
        "cgns/box_hex.cgns",
        "cgns/box_prism.cgns",
        "cgns/dual_ngon.cgns",
        "multizone/rotor_stator.cgns",
    ] {
        let nodes = {
            let f = hdf5_metno::File::open(data(file)).unwrap();
            let mut v = Vec::new();
            cgns_int_nodes(&f, "", &mut v);
            v
        };
        assert!(
            nodes.len() >= 4,
            "{file}: only {} integer nodes",
            nodes.len()
        );
        let picks = spread(nodes.len(), 48);
        let damage: [(&str, Damage); 6] = [
            ("all -1", |_, _| -1),
            ("all 0", |_, _| 0),
            ("all i32::MAX", |_, _| i32::MAX as i64),
            (
                "first i32::MAX",
                |i, v| if i == 0 { i32::MAX as i64 } else { v },
            ),
            (
                "first -2^31",
                |i, v| if i == 0 { i32::MIN as i64 } else { v },
            ),
            ("every value + 7", |_, v| v.saturating_add(7)),
        ];
        for &k in &picks {
            let node = &nodes[k];
            for (what, f) in damage {
                let dir = scratch(&format!("cgns-node-{}-{k}", file.replace('/', "-")));
                let p = dir.join("mesh.cgns");
                std::fs::copy(data(file), &p).unwrap();
                {
                    let h = hdf5_metno::File::open_rw(&p).unwrap();
                    let d = h.dataset(&format!("{node}/ data")).unwrap();
                    let v: Vec<i64> = d.read_raw().unwrap();
                    let w: Vec<i64> = v.iter().enumerate().map(|(i, &x)| f(i, x)).collect();
                    // values that do not fit an int32 node stay as they are
                    if d.dtype().unwrap().size() == 4 {
                        let w: Vec<i32> = w
                            .iter()
                            .zip(&v)
                            .map(|(&a, &b)| i32::try_from(a).unwrap_or(b as i32))
                            .collect();
                        d.write_raw(&w).unwrap();
                    } else {
                        d.write_raw(&w).unwrap();
                    }
                }
                let label = format!("cgns reader, {file}, node {node}, {what}");
                let q = p.clone();
                guarded(&label, move || {
                    reader("cgns").read_path(&q, &ReadOptions::default())
                });
                let _ = std::fs::remove_dir_all(&dir);
            }
        }
    }
}

/// The hand-written inputs of `tests/data/broken`: `<name>.<ext>` read
/// with the reader the extension (or the `reader` prefix of the name) names;
/// each must end in an error.
#[test]
fn hand_written_broken_inputs_are_refused() {
    let dir = data("broken");
    let mut n = 0;
    let mut entries: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    entries.sort();
    for p in entries {
        let fname = p.file_name().unwrap().to_string_lossy().into_owned();
        if fname.ends_with(".md") {
            continue;
        }
        let name: &'static str = match fname.split('_').next().unwrap() {
            "gmsh" => "gmsh",
            "su2" => "su2",
            "industry" => "industry_msh",
            "vtk" => "vtk",
            "openfoam" => "openfoam",
            other => panic!("{fname}: unknown reader prefix {other}"),
        };
        let q = p.clone();
        let label = format!("{name} reader, broken/{fname}");
        let o = guarded(&label, move || {
            reader(name).read_path(&q, &ReadOptions::default())
        });
        assert!(matches!(o, Outcome::Error), "{label}: read as a mesh");
        n += 1;
    }
    assert!(n >= 5, "only {n} hand-written inputs");
}
