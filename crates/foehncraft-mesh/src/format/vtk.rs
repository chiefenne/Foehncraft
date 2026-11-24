//! Reader and writer for VTK unstructured grids (D192): legacy `.vtk` and
//! XML `.vtu`.
//!
//! Written from the public descriptions of the two file formats (the VTK file
//! formats guide).
//!
//! **Reading.** Legacy files in ASCII or BINARY (big-endian) form, also the
//! 5.1 layout with `OFFSETS` and `CONNECTIVITY`; XML files with the data in
//! ASCII, inline base64 or appended (raw or base64) form, little- or
//! big-endian, 32- or 64-bit headers, uncompressed or zlib-compressed (the
//! default of ParaView; inflated by the decoder of this crate). Cell types
//! tetrahedron (10), hexahedron (12), voxel (11), wedge (13), pyramid (14) and
//! polyhedron (42) are the cells; triangles (5), polygons (7) and
//! quadrilaterals (9) are boundary faces; vertices and lines are ignored; a
//! higher-order or 2D cell is refused (a 2D VTK mesh: extrude it in the
//! mesher).
//!
//! **Streaming.** A file is read once from the front, through a buffer of
//! 1 MiB; neither the file nor the text of an array is ever held whole.
//! The numbers go straight from the text, the base64 decoder or the
//! inflated block into compact arrays: node ids as `u32`, the face streams
//! of the polyhedra as one flat array that moves into the
//! [`ElementMesh`](super::ElementMesh) without a copy. Appended data are
//! read array by array at their offsets, and only the arrays the mesh needs
//! are decoded. A mesh so takes about the memory of its arrays: a KoflCraft
//! mesh of 10 M polyhedra is some 12 GB of ASCII text, which earlier took
//! five times that in memory.
//!
//! **Names.** A grid has numbers where Foehncraft has names. The region
//! numbers are one scalar cell-data array, the first found of `region`,
//! `gmsh:physical`, `CellEntityIds`, `PatchID`, `patch`, `boundary`, `zone`,
//! `MaterialIds` (also in other capitals): the volume cells give the cell
//! zones, the surface cells the boundaries. Names for the numbers come from
//! the string arrays `zone_names` and `boundary_names` of the field data
//! (entries `id=name`, as [`write_vtu`] writes them); else a zone is `fluid`
//! (one zone) or `zone<id>`, a boundary `boundary<id>`. Surface cells without
//! a region array form the boundary `default`; boundary faces that no surface
//! cell names go to `default` as for every format.
//!
//! **Writing.** [`write_vtu`] writes the mesh as a `.vtu` file: every cell a
//! polyhedron, every boundary face a polygon, the region array and the names
//! as above, ASCII or appended binary. The file reads back as the same mesh
//! (also in ParaView).

use super::{
    CellKind, ElementMesh, MeshFile, MeshFormat, Polyhedra, ReadOptions, check_face_stream,
    declared, inflate,
};
use crate::mesh::Mesh;
use std::collections::{BTreeMap, HashMap};
use std::io::{ErrorKind, Read, Seek, SeekFrom};
use std::path::Path;

/// The VTK format in the registry of [`crate::format`].
#[derive(Debug, Clone, Copy)]
pub struct Vtk;

impl MeshFormat for Vtk {
    fn name(&self) -> &'static str {
        "vtk"
    }

    fn title(&self) -> &'static str {
        "VTK unstructured grid, legacy .vtk and XML .vtu"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["vtk", "vtu"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        let start = head
            .iter()
            .position(|c| !c.is_ascii_whitespace())
            .unwrap_or(0);
        let h = &head[start..];
        h.starts_with(b"# vtk DataFile") || h.starts_with(b"<?xml") || h.starts_with(b"<VTKFile")
    }

    fn parse(&self, bytes: &[u8], _options: &ReadOptions) -> Result<MeshFile, String> {
        read(Stream::new(std::io::Cursor::new(bytes), bytes.len() as u64))
    }

    /// Streams the file (see the module docs).
    fn read_path(&self, path: &Path, _options: &ReadOptions) -> Result<MeshFile, String> {
        let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let len = file
            .metadata()
            .map_err(|e| format!("{}: {e}", path.display()))?
            .len();
        read(Stream::new(file, len))
    }
}

fn read<R: Read + Seek>(mut s: Stream<R>) -> Result<MeshFile, String> {
    let (grid, how) = if s.ahead(5)?.starts_with(b"# vtk") {
        parse_legacy(&mut s)?
    } else {
        parse_xml(&mut s)?
    };
    drop(s);
    grid.into_mesh_file(how)
}

/// What a file holds, before it becomes a mesh.
#[derive(Default)]
struct Grid {
    points: Vec<[f64; 3]>,
    /// VTK cell types.
    types: Vec<u8>,
    /// Where the entry of each cell in `conn` ends; the first starts at
    /// `first`.
    offsets: Vec<i64>,
    first: i64,
    /// Node lists of the cells (of a polyhedron in a legacy file: its face
    /// stream).
    conn: Vec<u32>,
    /// XML: the face streams of the polyhedra, and where each cell's
    /// stream ends (`faceoffsets`, -1 for the other cells).
    faces: Vec<u32>,
    faceoffsets: Vec<i64>,
    /// A legacy file (the face streams are in `conn`).
    legacy: bool,
    /// The scalar cell-data arrays that can hold region numbers.
    arrays: Vec<(String, Vec<i64>)>,
    /// String arrays of the field data.
    names: HashMap<String, Vec<String>>,
}

const REGION_NAMES: [&str; 12] = [
    "region",
    "gmsh:physical",
    "cellentityids",
    "regionid",
    "regionids",
    "patchid",
    "patch",
    "boundary",
    "boundaryid",
    "zone",
    "zoneid",
    "materialids",
];

fn is_region(name: &str) -> bool {
    REGION_NAMES.iter().any(|r| name.eq_ignore_ascii_case(r))
}

fn id_names(list: Option<&Vec<String>>) -> HashMap<i64, String> {
    list.into_iter()
        .flatten()
        .filter_map(|e| {
            let (i, n) = e.split_once('=')?;
            Some((i.trim().parse().ok()?, n.trim().to_string()))
        })
        .collect()
}

impl Grid {
    fn into_mesh_file(self, how: String) -> Result<MeshFile, String> {
        let Grid {
            points,
            types,
            offsets,
            first,
            conn,
            faces,
            faceoffsets,
            legacy,
            arrays,
            names,
        } = self;
        let region: Option<&Vec<i64>> = REGION_NAMES.iter().find_map(|want| {
            arrays
                .iter()
                .find(|(n, _)| n.eq_ignore_ascii_case(want))
                .map(|(_, v)| v)
        });
        let id = |i: usize| region.and_then(|r| r.get(i)).copied();
        let zone_names = id_names(names.get("zone_names"));
        let boundary_names = id_names(names.get("boundary_names"));

        // the zones and boundaries, in the order of their numbers
        let mut zone_ids: BTreeMap<i64, ()> = BTreeMap::new();
        let mut bound_ids: BTreeMap<i64, ()> = BTreeMap::new();
        for (i, t) in types.iter().enumerate() {
            match t {
                10 | 11 | 12 | 13 | 14 | 42 => {
                    zone_ids.insert(id(i).unwrap_or(0), ());
                }
                5 | 7 | 9 => {
                    bound_ids.insert(id(i).unwrap_or(i64::MIN), ());
                }
                _ => {}
            }
        }
        let n = points.len() as u32;
        let mut em = ElementMesh {
            points,
            face_word: "surface cells",
            group_word: "region",
            ..Default::default()
        };
        let single = zone_ids.len() == 1;
        let zone_of: HashMap<i64, u32> = zone_ids
            .keys()
            .map(|&k| {
                let name = zone_names.get(&k).cloned().unwrap_or_else(|| {
                    if single {
                        "fluid".to_string()
                    } else {
                        format!("zone{k}")
                    }
                });
                (k, em.zone(&name))
            })
            .collect();
        let bound_of: HashMap<i64, u32> = bound_ids
            .keys()
            .map(|&k| {
                let name = boundary_names.get(&k).cloned().unwrap_or_else(|| {
                    if k == i64::MIN {
                        "default".to_string()
                    } else {
                        format!("boundary{k}")
                    }
                });
                (k, em.boundary(&name))
            })
            .collect();

        // the face streams of the polyhedra stay where they are: in `faces`
        // (XML) or in `conn` (legacy); only their starts are noted
        let stream: &[u32] = if legacy { &conn } else { &faces };
        let mut starts: Vec<usize> = Vec::new();
        let mut n_volume = 0usize;
        let mut start = first;
        let mut fstart = 0usize;
        for (i, (&t, &end)) in types.iter().zip(&offsets).enumerate() {
            let range = usize::try_from(start)
                .ok()
                .zip(usize::try_from(end).ok())
                .filter(|&(a, b)| a <= b && b <= conn.len())
                .ok_or("the offsets go beyond the connectivity")?;
            start = end;
            let nodes = &conn[range.0..range.1];
            if t != 42
                && let Some(&v) = nodes.iter().find(|&&v| v >= n)
            {
                return Err(format!("cell {i} refers to the point {v}, beyond {n}"));
            }
            let need = |k: usize| -> Result<(), String> {
                if nodes.len() == k {
                    Ok(())
                } else {
                    Err(format!(
                        "cell {i} of type {t} has {} points, expected {k}",
                        nodes.len()
                    ))
                }
            };
            let z = *zone_of.get(&id(i).unwrap_or(0)).unwrap_or(&0);
            match t {
                10 => {
                    need(4)?;
                    em.add_cell(CellKind::Tetrahedron, nodes, z);
                }
                12 => {
                    need(8)?;
                    em.add_cell(CellKind::Hexahedron, nodes, z);
                }
                11 => {
                    need(8)?;
                    let v = [0, 1, 3, 2, 4, 5, 7, 6].map(|k| nodes[k]);
                    em.add_cell(CellKind::Hexahedron, &v, z);
                }
                13 => {
                    need(6)?;
                    em.add_cell(CellKind::Prism, nodes, z);
                }
                14 => {
                    need(5)?;
                    em.add_cell(CellKind::Pyramid, nodes, z);
                }
                42 => {
                    let (a, b) = if legacy {
                        range
                    } else {
                        if faceoffsets.is_empty() {
                            return Err("polyhedra need the faces and faceoffsets arrays".into());
                        }
                        let fend = *faceoffsets.get(i).ok_or("faceoffsets ends too early")?;
                        let fend = usize::try_from(fend)
                            .ok()
                            .filter(|&e| e >= fstart && e <= faces.len())
                            .ok_or("faceoffsets go beyond the faces")?;
                        let r = (fstart, fend);
                        fstart = fend;
                        r
                    };
                    check_face_stream(&stream[a..b], i)?;
                    starts.push(a);
                    em.cell_kind.push(CellKind::Polyhedron);
                    em.cell_zone.push(z);
                }
                5 | 7 | 9 => {
                    if nodes.len() < 3 {
                        return Err(format!("surface cell {i} has {} points", nodes.len()));
                    }
                    let b = *bound_of.get(&id(i).unwrap_or(i64::MIN)).unwrap_or(&0);
                    em.add_face(nodes, b);
                    continue;
                }
                1 | 2 | 3 | 4 | 6 | 8 => continue, // points, lines, strips, pixels
                other => {
                    return Err(format!(
                        "cell {i} has the type {other}: only linear tetrahedra, hexahedra, wedges, \
                         pyramids and polyhedra (and triangles and quadrilaterals as boundary \
                         faces) are read"
                    ));
                }
            }
            n_volume += 1;
        }
        if n_volume == 0 {
            return Err("the grid has no 3D cells (a 2D mesh: extrude it in the mesher)".into());
        }
        drop((types, offsets, faceoffsets, arrays));
        em.poly = Polyhedra::from_parts(if legacy { conn } else { faces }, starts);
        let built = em.into_mesh()?;
        Ok(MeshFile {
            mesh: built.mesh,
            reader: "vtk",
            format: how,
            warnings: built.warnings,
        })
    }
}

// ---------------------------------------------------------------------
// the stream

/// Size of the read buffer and of the pieces of binary data.
const CHUNK: usize = 1 << 20;

/// Longest tag read (a longer one is a broken file).
const MAX_TAG: usize = 1 << 16;

/// A file read through a buffer of its own: a lookahead of a few bytes for
/// the keywords of legacy files, and seeks to the appended arrays.
struct Stream<R> {
    inner: R,
    buf: Vec<u8>,
    /// The unread bytes are `buf[at..end]`.
    at: usize,
    end: usize,
    /// Position of `buf[at]` in the file.
    pos: u64,
    /// Size of the file: the bound of the counts it declares (D202).
    len: u64,
}

fn io_err(e: std::io::Error) -> String {
    format!("reading the file: {e}")
}

impl<R: Read + Seek> Stream<R> {
    fn new(inner: R, len: u64) -> Self {
        Self::with_buffer(inner, len, CHUNK)
    }

    /// A stream with a buffer of `cap` bytes (at least 64: the lookahead
    /// of the keywords).
    fn with_buffer(inner: R, len: u64, cap: usize) -> Self {
        Self {
            inner,
            buf: vec![0; cap.max(64)],
            at: 0,
            end: 0,
            pos: 0,
            len,
        }
    }

    /// Reads more bytes into the buffer; false at the end of the file.
    fn more(&mut self) -> Result<bool, String> {
        if self.at > 0 {
            self.buf.copy_within(self.at..self.end, 0);
            self.end -= self.at;
            self.at = 0;
        }
        loop {
            match self.inner.read(&mut self.buf[self.end..]) {
                Ok(0) => return Ok(false),
                Ok(n) => {
                    self.end += n;
                    return Ok(true);
                }
                Err(e) if e.kind() == ErrorKind::Interrupted => {}
                Err(e) => return Err(io_err(e)),
            }
        }
    }

    /// The unread bytes of the buffer, at least `k` of them unless the file
    /// ends first (`k` at most [`CHUNK`]).
    fn ahead(&mut self, k: usize) -> Result<&[u8], String> {
        while self.end - self.at < k && self.more()? {}
        Ok(&self.buf[self.at..self.end])
    }

    /// The unread bytes of the buffer (empty at the end of the file).
    fn data(&mut self) -> Result<&[u8], String> {
        self.ahead(1)
    }

    fn consume(&mut self, n: usize) {
        self.at += n;
        self.pos += n as u64;
    }

    fn peek(&mut self) -> Result<Option<u8>, String> {
        Ok(self.data()?.first().copied())
    }

    /// Skips bytes while `keep` holds.
    fn skip_while(&mut self, keep: impl Fn(u8) -> bool) -> Result<(), String> {
        loop {
            let d = self.data()?;
            if d.is_empty() {
                return Ok(());
            }
            match d.iter().position(|&c| !keep(c)) {
                Some(k) => {
                    self.consume(k);
                    return Ok(());
                }
                None => {
                    let k = d.len();
                    self.consume(k);
                }
            }
        }
    }

    /// Skips past the first `end` (`err` if the file ends before).
    fn skip_past(&mut self, end: &[u8], err: &str) -> Result<(), String> {
        loop {
            let d = self.ahead(end.len())?;
            if d.len() < end.len() {
                return Err(err.into());
            }
            match find(d, end) {
                Some(k) => {
                    self.consume(k + end.len());
                    return Ok(());
                }
                None => {
                    let k = d.len() + 1 - end.len();
                    self.consume(k);
                }
            }
        }
    }

    /// Moves to the byte `pos` of the file.
    fn seek(&mut self, pos: u64) -> Result<(), String> {
        if pos >= self.pos && pos - self.pos <= (self.end - self.at) as u64 {
            self.consume((pos - self.pos) as usize);
            return Ok(());
        }
        self.inner.seek(SeekFrom::Start(pos)).map_err(io_err)?;
        self.at = 0;
        self.end = 0;
        self.pos = pos;
        Ok(())
    }

    /// A word of a legacy file (bytes up to white space, at most `max`);
    /// false at the end of the file.
    fn word(&mut self, out: &mut Vec<u8>) -> Result<bool, String> {
        self.skip_while(|c| c.is_ascii_whitespace())?;
        out.clear();
        loop {
            let d = self.data()?;
            if d.is_empty() {
                break;
            }
            let k = d
                .iter()
                .position(|c| c.is_ascii_whitespace())
                .unwrap_or(d.len());
            out.extend_from_slice(&d[..k]);
            let more = k == d.len();
            self.consume(k);
            if out.len() > MAX_TAG {
                return Err("a word of more than 64 KiB (not a VTK file?)".into());
            }
            if !more {
                break;
            }
        }
        Ok(!out.is_empty())
    }

    /// A line of a legacy file, without its line break.
    fn line(&mut self) -> Result<String, String> {
        let mut out = Vec::new();
        loop {
            let d = self.data()?;
            if d.is_empty() {
                break;
            }
            match d.iter().position(|&c| c == b'\n') {
                Some(k) => {
                    out.extend_from_slice(&d[..k]);
                    self.consume(k + 1);
                    break;
                }
                None => {
                    out.extend_from_slice(d);
                    let k = d.len();
                    self.consume(k);
                }
            }
            if out.len() > MAX_TAG {
                return Err("a line of more than 64 KiB (not a VTK file?)".into());
            }
        }
        Ok(String::from_utf8_lossy(&out).trim().to_string())
    }

    /// Numbers in ASCII: `count` of them (legacy files), or up to the `<`
    /// that ends the text of an XML element (`None`). Returns how many.
    fn numbers<S: Sink>(
        &mut self,
        count: Option<u64>,
        sink: &mut S,
        what: &str,
    ) -> Result<u64, String> {
        let mut tok: Vec<u8> = Vec::with_capacity(32);
        let mut got = 0u64;
        loop {
            if count == Some(got) {
                return Ok(got);
            }
            let d = self.data()?;
            if d.is_empty() {
                if !tok.is_empty() {
                    number(&tok, sink, what)?;
                    got += 1;
                }
                if count.is_some_and(|c| got < c) {
                    return Err(format!("the file ends where {what} is expected"));
                }
                return Ok(got);
            }
            let mut i = 0;
            let mut stop = false;
            while i < d.len() {
                let c = d[i];
                i += 1;
                if c.is_ascii_whitespace() {
                    if !tok.is_empty() {
                        number(&tok, sink, what)?;
                        tok.clear();
                        got += 1;
                        if count == Some(got) {
                            break;
                        }
                    }
                } else if c == b'<' && count.is_none() {
                    i -= 1;
                    stop = true;
                    break;
                } else {
                    tok.push(c);
                    if tok.len() > 64 {
                        return Err(format!(
                            "{what}: '{}…' is not a number",
                            String::from_utf8_lossy(&tok[..16])
                        ));
                    }
                }
            }
            self.consume(i);
            if stop {
                if !tok.is_empty() {
                    number(&tok, sink, what)?;
                    got += 1;
                }
                return Ok(got);
            }
        }
    }
}

impl<R: Read + Seek> Read for Stream<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if self.at == self.end {
            if out.len() >= self.buf.len() {
                // large reads go past the buffer
                let n = self.inner.read(out)?;
                self.pos += n as u64;
                return Ok(n);
            }
            self.more().map_err(std::io::Error::other)?;
        }
        let n = out.len().min(self.end - self.at);
        out[..n].copy_from_slice(&self.buf[self.at..self.at + n]);
        self.consume(n);
        Ok(n)
    }
}

/// One number of an ASCII text: whole numbers without a detour through
/// `f64`, anything else as a real number.
fn number<S: Sink>(t: &[u8], sink: &mut S, what: &str) -> Result<(), String> {
    let (neg, digits) = match t.first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if !digits.is_empty() && digits.len() <= 18 && digits.iter().all(u8::is_ascii_digit) {
        let v = digits
            .iter()
            .fold(0i64, |v, &d| v * 10 + i64::from(d - b'0'));
        return sink.int(if neg { -v } else { v });
    }
    let v = std::str::from_utf8(t)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .ok_or_else(|| format!("{what}: '{}' is not a number", String::from_utf8_lossy(t)))?;
    sink.real(v)
}

/// Fills `buf` from `rd`, or says that the binary data end too early.
fn fill<Rd: Read>(rd: &mut Rd, buf: &mut [u8]) -> Result<(), String> {
    rd.read_exact(buf).map_err(|e| match e.kind() {
        ErrorKind::UnexpectedEof => "the binary data ends too early".to_string(),
        _ => e.to_string(),
    })
}

/// Base64 text, decoded as it is read. A `=` ends a group of four
/// characters, so pieces encoded one after the other (a header, then the
/// data) decode as one stream of bytes; a `<` (the end of the text of an
/// element) or the end of the file ends it.
struct B64<'s, R> {
    s: &'s mut Stream<R>,
    acc: u32,
    n: u32,
    /// Decoded bytes not yet handed out: `out[at..len]`.
    out: [u8; 3],
    at: usize,
    len: usize,
    done: bool,
}

impl<'s, R> B64<'s, R> {
    fn new(s: &'s mut Stream<R>) -> Self {
        Self {
            s,
            acc: 0,
            n: 0,
            out: [0; 3],
            at: 0,
            len: 0,
            done: false,
        }
    }
}

/// Ends a group of base64 characters (at padding or at the end of the
/// text): two characters give one byte, three give two. Returns how many
/// bytes are in `out`.
fn b64_end(acc: &mut u32, n: &mut u32, out: &mut [u8; 3]) -> usize {
    let a = *acc;
    let len = match *n {
        2 => {
            out[0] = (a >> 4) as u8;
            1
        }
        3 => {
            out[0] = (a >> 10) as u8;
            out[1] = (a >> 2) as u8;
            2
        }
        _ => 0,
    };
    *acc = 0;
    *n = 0;
    len
}

impl<R: Read + Seek> Read for B64<'_, R> {
    fn read(&mut self, dst: &mut [u8]) -> std::io::Result<usize> {
        let mut w = 0;
        loop {
            while self.at < self.len && w < dst.len() {
                dst[w] = self.out[self.at];
                self.at += 1;
                w += 1;
            }
            if w == dst.len() || self.done {
                return Ok(w);
            }
            let d = self.s.data().map_err(std::io::Error::other)?;
            if d.is_empty() {
                self.done = true;
                self.len = b64_end(&mut self.acc, &mut self.n, &mut self.out);
                self.at = 0;
                continue;
            }
            let mut used = 0;
            for &c in d {
                let v = match c {
                    b'A'..=b'Z' => c - b'A',
                    b'a'..=b'z' => c - b'a' + 26,
                    b'0'..=b'9' => c - b'0' + 52,
                    b'+' => 62,
                    b'/' => 63,
                    b'<' => {
                        self.done = true;
                        self.len = b64_end(&mut self.acc, &mut self.n, &mut self.out);
                        self.at = 0;
                        break;
                    }
                    b'=' => {
                        used += 1;
                        self.len = b64_end(&mut self.acc, &mut self.n, &mut self.out);
                        self.at = 0;
                        if self.len > 0 {
                            break;
                        }
                        continue;
                    }
                    c if c.is_ascii_whitespace() => {
                        used += 1;
                        continue;
                    }
                    _ => {
                        return Err(std::io::Error::new(
                            ErrorKind::InvalidData,
                            "bad base64 data",
                        ));
                    }
                };
                used += 1;
                self.acc = (self.acc << 6) | u32::from(v);
                self.n += 1;
                if self.n == 4 {
                    let a = self.acc;
                    let b = [(a >> 16) as u8, (a >> 8) as u8, a as u8];
                    self.acc = 0;
                    self.n = 0;
                    if w + 3 <= dst.len() {
                        dst[w..w + 3].copy_from_slice(&b);
                        w += 3;
                        if w == dst.len() {
                            break;
                        }
                    } else {
                        self.out = b;
                        self.at = 0;
                        self.len = 3;
                        break;
                    }
                }
            }
            self.s.consume(used);
        }
    }
}

// ---------------------------------------------------------------------
// where the numbers go

/// Where the numbers of an array go.
trait Sink {
    /// A whole number of the file.
    fn int(&mut self, v: i64) -> Result<(), String>;
    /// A real number of the file.
    fn real(&mut self, v: f64) -> Result<(), String>;
    /// The number of values to come, when the data say it (bounded by the
    /// size of the file).
    fn reserve(&mut self, _n: usize) {}
}

/// Node ids and counts.
struct Ids<'a> {
    out: &'a mut Vec<u32>,
    what: &'static str,
}

impl Sink for Ids<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        let v = u32::try_from(v).map_err(|_| format!("{}: {v} is not a valid index", self.what))?;
        self.out.push(v);
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.out.push(to_u32(v, self.what)?);
        Ok(())
    }

    fn reserve(&mut self, n: usize) {
        self.out.reserve(n);
    }
}

/// Offsets and region numbers (a real number is rounded).
struct Ints<'a>(&'a mut Vec<i64>);

impl Sink for Ints<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        self.0.push(v);
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.0.push(v.round() as i64);
        Ok(())
    }

    fn reserve(&mut self, n: usize) {
        self.0.reserve(n);
    }
}

/// Offsets where the first value is the start of the first cell (legacy
/// files of version 5).
struct Starts<'a> {
    first: &'a mut i64,
    out: &'a mut Vec<i64>,
    seen: bool,
}

impl Sink for Starts<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        if self.seen {
            self.out.push(v);
        } else {
            *self.first = v;
            self.seen = true;
        }
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.int(v.round() as i64)
    }
}

/// Cell types.
struct Types<'a>(&'a mut Vec<u8>);

impl Sink for Types<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        let t = u8::try_from(v).map_err(|_| format!("types: {v} is not a cell type"))?;
        self.0.push(t);
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.int(v as i64)
    }

    fn reserve(&mut self, n: usize) {
        self.0.reserve(n);
    }
}

/// Points, three numbers each (an incomplete last point is dropped).
struct Points<'a> {
    out: &'a mut Vec<[f64; 3]>,
    p: [f64; 3],
    k: usize,
}

impl<'a> Points<'a> {
    fn new(out: &'a mut Vec<[f64; 3]>) -> Self {
        Self {
            out,
            p: [0.0; 3],
            k: 0,
        }
    }
}

impl Sink for Points<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        self.real(v as f64)
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.p[self.k] = v;
        self.k += 1;
        if self.k == 3 {
            self.out.push(self.p);
            self.k = 0;
        }
        Ok(())
    }

    fn reserve(&mut self, n: usize) {
        self.out.reserve(n / 3);
    }
}

/// Character codes of a string array.
struct Codes<'a>(&'a mut Vec<u8>);

impl Sink for Codes<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        self.0.push(v as u8);
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.0.push(v as u8);
        Ok(())
    }
}

/// An array the mesh does not need.
struct Skip;

impl Sink for Skip {
    fn int(&mut self, _: i64) -> Result<(), String> {
        Ok(())
    }

    fn real(&mut self, _: f64) -> Result<(), String> {
        Ok(())
    }
}

/// The cell list of a legacy file before version 5: for each cell its
/// number of entries, then the entries. Values after `want` cells are
/// ignored.
struct CellList<'a> {
    conn: &'a mut Vec<u32>,
    ends: &'a mut Vec<i64>,
    left: u64,
    cells: usize,
    want: usize,
}

impl Sink for CellList<'_> {
    fn int(&mut self, v: i64) -> Result<(), String> {
        if self.left == 0 {
            if self.cells == self.want {
                return Ok(());
            }
            self.left = u64::try_from(v).map_err(|_| format!("CELLS: {v} entries"))?;
        } else {
            let v = u32::try_from(v).map_err(|_| format!("CELLS: {v} is not a valid index"))?;
            self.conn.push(v);
            self.left -= 1;
        }
        if self.left == 0 {
            self.ends.push(self.conn.len() as i64);
            self.cells += 1;
        }
        Ok(())
    }

    fn real(&mut self, v: f64) -> Result<(), String> {
        self.int(i64::from(to_u32(v, "CELLS")?))
    }
}

// ---------------------------------------------------------------------
// numbers in binary form

#[derive(Clone, Copy)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    I64,
    U64,
    F32,
    F64,
}

impl Ty {
    fn legacy(name: &str) -> Result<Self, String> {
        Ok(match name.to_ascii_lowercase().as_str() {
            "char" => Ty::I8,
            "unsigned_char" => Ty::U8,
            "short" => Ty::I16,
            "unsigned_short" => Ty::U16,
            "int" => Ty::I32,
            "unsigned_int" => Ty::U32,
            "long" | "vtktypeint64" | "vtkidtype" => Ty::I64,
            "unsigned_long" | "vtktypeuint64" => Ty::U64,
            "float" => Ty::F32,
            "double" => Ty::F64,
            other => return Err(format!("unknown data type '{other}'")),
        })
    }

    fn xml(name: &str) -> Result<Self, String> {
        Ok(match name {
            "Int8" => Ty::I8,
            "UInt8" => Ty::U8,
            "Int16" => Ty::I16,
            "UInt16" => Ty::U16,
            "Int32" => Ty::I32,
            "UInt32" => Ty::U32,
            "Int64" => Ty::I64,
            "UInt64" => Ty::U64,
            "Float32" => Ty::F32,
            "Float64" => Ty::F64,
            other => return Err(format!("unknown data type '{other}'")),
        })
    }

    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::I64 | Ty::U64 | Ty::F64 => 8,
        }
    }
}

/// Binary numbers of one type, given in pieces that may split a number.
struct Feed {
    ty: Ty,
    big: bool,
    carry: [u8; 8],
    nc: usize,
}

impl Feed {
    fn new(ty: Ty, big: bool) -> Self {
        Self {
            ty,
            big,
            carry: [0; 8],
            nc: 0,
        }
    }

    fn push<S: Sink>(&mut self, mut b: &[u8], sink: &mut S) -> Result<(), String> {
        let k = self.ty.size();
        if self.nc > 0 {
            let take = (k - self.nc).min(b.len());
            self.carry[self.nc..self.nc + take].copy_from_slice(&b[..take]);
            self.nc += take;
            b = &b[take..];
            if self.nc < k {
                return Ok(());
            }
            self.nc = 0;
            let c = self.carry;
            self.values(&c[..k], sink)?;
        }
        let whole = b.len() / k * k;
        self.values(&b[..whole], sink)?;
        let rest = &b[whole..];
        self.carry[..rest.len()].copy_from_slice(rest);
        self.nc = rest.len();
        Ok(())
    }

    fn values<S: Sink>(&self, b: &[u8], sink: &mut S) -> Result<(), String> {
        let big = self.big;
        macro_rules! each {
            ($t:ty, $n:literal, $f:ident) => {
                for &c in b.as_chunks::<$n>().0 {
                    let v = if big {
                        <$t>::from_be_bytes(c)
                    } else {
                        <$t>::from_le_bytes(c)
                    };
                    sink.$f(v.into())?;
                }
            };
        }
        match self.ty {
            Ty::I8 => each!(i8, 1, int),
            Ty::U8 => each!(u8, 1, int),
            Ty::I16 => each!(i16, 2, int),
            Ty::U16 => each!(u16, 2, int),
            Ty::I32 => each!(i32, 4, int),
            Ty::U32 => each!(u32, 4, int),
            Ty::I64 => each!(i64, 8, int),
            Ty::F32 => each!(f32, 4, real),
            Ty::F64 => each!(f64, 8, real),
            Ty::U64 => {
                for &c in b.as_chunks::<8>().0 {
                    let v = if big {
                        u64::from_be_bytes(c)
                    } else {
                        u64::from_le_bytes(c)
                    };
                    match i64::try_from(v) {
                        Ok(v) => sink.int(v)?,
                        Err(_) => sink.real(v as f64)?,
                    }
                }
            }
        }
        Ok(())
    }
}

fn to_u32(v: f64, what: &str) -> Result<u32, String> {
    if v >= 0.0 && v <= u32::MAX as f64 && v.fract() == 0.0 {
        Ok(v as u32)
    } else {
        Err(format!("{what}: {v} is not a valid index"))
    }
}

/// `n` bytes of binary numbers from `rd` into `sink`, a piece at a time.
fn binary<Rd: Read, S: Sink>(
    rd: &mut Rd,
    n: u64,
    feed: &mut Feed,
    sink: &mut S,
) -> Result<(), String> {
    let mut buf = vec![0u8; (n as usize).min(CHUNK)];
    let mut left = n;
    while left > 0 {
        let k = (left as usize).min(CHUNK);
        fill(rd, &mut buf[..k])?;
        feed.push(&buf[..k], sink)?;
        left -= k as u64;
    }
    Ok(())
}

// ---------------------------------------------------------------------
// legacy files

/// The data of a legacy keyword: `n` numbers of type `ty`.
fn legacy_values<R: Read + Seek, S: Sink>(
    s: &mut Stream<R>,
    n: u64,
    ty: Ty,
    bin: bool,
    what: &str,
    sink: &mut S,
) -> Result<(), String> {
    let n = declared(n as usize, s.len as usize, what)? as u64;
    if bin {
        // one line break follows the header of binary data
        if s.peek()? == Some(b'\r') {
            s.consume(1);
        }
        if s.peek()? == Some(b'\n') {
            s.consume(1);
        }
        let bytes = n * ty.size() as u64;
        if bytes > s.len {
            return Err(format!("{what}: the binary data ends too early"));
        }
        sink.reserve(n as usize);
        binary(s, bytes, &mut Feed::new(ty, true), sink).map_err(|e| format!("{what}: {e}"))
    } else {
        s.numbers(Some(n), sink, what).map(|_| ())
    }
}

fn need_word<R: Read + Seek>(
    s: &mut Stream<R>,
    w: &mut Vec<u8>,
    what: &str,
) -> Result<String, String> {
    if !s.word(w)? {
        return Err(format!("the file ends where {what} is expected"));
    }
    Ok(String::from_utf8_lossy(w).into_owned())
}

fn need_number<R: Read + Seek>(
    s: &mut Stream<R>,
    w: &mut Vec<u8>,
    what: &str,
) -> Result<u64, String> {
    let t = need_word(s, w, what)?;
    t.parse()
        .map_err(|_| format!("{what}: '{t}' is not a number"))
}

fn parse_legacy<R: Read + Seek>(s: &mut Stream<R>) -> Result<(Grid, String), String> {
    let first = s.line()?;
    let version = first
        .rsplit(' ')
        .next()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.0);
    s.line()?; // title
    let mode = s.line()?.to_ascii_uppercase();
    let bin = match mode.as_str() {
        "ASCII" => false,
        "BINARY" => true,
        other => return Err(format!("the file is '{other}', expected ASCII or BINARY")),
    };
    let mut g = Grid {
        legacy: true,
        ..Default::default()
    };
    let mut cell_data = false;
    let mut n_cells = 0usize;
    let mut count_data = 0u64;
    let mut w = Vec::new();
    while s.word(&mut w)? {
        let key = String::from_utf8_lossy(&w).to_ascii_uppercase();
        match key.as_str() {
            "DATASET" => {
                let k = need_word(s, &mut w, "the dataset type")?;
                if !k.eq_ignore_ascii_case("UNSTRUCTURED_GRID") {
                    return Err(format!(
                        "the dataset is {k}: only UNSTRUCTURED_GRID is read"
                    ));
                }
            }
            "POINTS" => {
                let n = need_number(s, &mut w, "the number of points")?;
                let ty = Ty::legacy(&need_word(s, &mut w, "the point type")?)?;
                g.points.clear();
                let n3 = n.checked_mul(3).ok_or("POINTS: too many points")?;
                legacy_values(s, n3, ty, bin, "POINTS", &mut Points::new(&mut g.points))?;
            }
            "CELLS" => {
                let n1 = need_number(s, &mut w, "the number of cells")?;
                let n2 = need_number(s, &mut w, "the size of the cell list")?;
                g.conn.clear();
                g.offsets.clear();
                g.first = 0;
                if version >= 5.0 && next_is(s, b"OFFSETS")? {
                    need_word(s, &mut w, "OFFSETS")?;
                    let ty = Ty::legacy(&need_word(s, &mut w, "the offsets type")?)?;
                    let mut sink = Starts {
                        first: &mut g.first,
                        out: &mut g.offsets,
                        seen: false,
                    };
                    legacy_values(s, n1, ty, bin, "OFFSETS", &mut sink)?;
                    need_word(s, &mut w, "CONNECTIVITY")?;
                    let ty = Ty::legacy(&need_word(s, &mut w, "the connectivity type")?)?;
                    let mut sink = Ids {
                        out: &mut g.conn,
                        what: "CELLS",
                    };
                    legacy_values(s, n2, ty, bin, "CONNECTIVITY", &mut sink)?;
                } else {
                    let want = declared(n1 as usize, s.len as usize, "CELLS")?;
                    let mut sink = CellList {
                        conn: &mut g.conn,
                        ends: &mut g.offsets,
                        left: 0,
                        cells: 0,
                        want,
                    };
                    legacy_values(s, n2, Ty::I32, bin, "CELLS", &mut sink)?;
                    if sink.cells < want {
                        return Err("CELLS ends too early".into());
                    }
                }
                n_cells = g.offsets.len();
            }
            "CELL_TYPES" => {
                let n = need_number(s, &mut w, "the number of cell types")?;
                g.types.clear();
                legacy_values(s, n, Ty::I32, bin, "CELL_TYPES", &mut Types(&mut g.types))?;
            }
            "CELL_DATA" => {
                count_data = need_number(s, &mut w, "the number of cell data")?;
                cell_data = true;
            }
            "POINT_DATA" => {
                count_data = need_number(s, &mut w, "the number of point data")?;
                cell_data = false;
            }
            "SCALARS" => {
                let name = need_word(s, &mut w, "the array name")?;
                let ty = Ty::legacy(&need_word(s, &mut w, "the array type")?)?;
                // an optional number of components, then LOOKUP_TABLE
                let mut comps = 1u64;
                let mut next = need_word(s, &mut w, "LOOKUP_TABLE")?;
                if let Ok(k) = next.parse::<u64>() {
                    comps = k;
                    next = need_word(s, &mut w, "LOOKUP_TABLE")?;
                }
                if !next.eq_ignore_ascii_case("LOOKUP_TABLE") {
                    return Err(format!(
                        "SCALARS {name}: expected LOOKUP_TABLE, found {next}"
                    ));
                }
                need_word(s, &mut w, "the table name")?;
                let n = count_data
                    .checked_mul(comps)
                    .ok_or("SCALARS: too many values")?;
                legacy_array(s, &mut g, name, cell_data && comps == 1, n, ty, bin)?;
            }
            "VECTORS" | "NORMALS" => {
                need_word(s, &mut w, "the array name")?;
                let ty = Ty::legacy(&need_word(s, &mut w, "the array type")?)?;
                let n = count_data
                    .checked_mul(3)
                    .ok_or("VECTORS: too many values")?;
                legacy_values(s, n, ty, bin, "VECTORS", &mut Skip)?;
            }
            "FIELD" => {
                need_word(s, &mut w, "the field name")?;
                let k = need_number(s, &mut w, "the number of arrays")?;
                for _ in 0..k {
                    let name = need_word(s, &mut w, "the array name")?;
                    let comps = need_number(s, &mut w, "the number of components")?;
                    let tuples = need_number(s, &mut w, "the number of tuples")?;
                    let ty = Ty::legacy(&need_word(s, &mut w, "the array type")?)?;
                    let n = comps.checked_mul(tuples).ok_or("FIELD: too many values")?;
                    legacy_array(s, &mut g, name, cell_data && comps == 1, n, ty, bin)?;
                }
            }
            other => {
                return Err(format!(
                    "the keyword {other} is not supported (read: POINTS, CELLS, CELL_TYPES and \
                     scalar cell data)"
                ));
            }
        }
    }
    if g.points.is_empty() {
        return Err("the file has no POINTS".into());
    }
    if g.types.len() != n_cells {
        return Err(format!("{n_cells} cells but {} cell types", g.types.len()));
    }
    let how = format!("VTK legacy {}", if bin { "binary" } else { "ASCII" });
    Ok((g, how))
}

/// A data array of a legacy file: kept when it is a scalar cell array
/// that can hold region numbers.
fn legacy_array<R: Read + Seek>(
    s: &mut Stream<R>,
    g: &mut Grid,
    name: String,
    scalar_cells: bool,
    n: u64,
    ty: Ty,
    bin: bool,
) -> Result<(), String> {
    if scalar_cells && is_region(&name) {
        let mut v = Vec::new();
        legacy_values(s, n, ty, bin, &name, &mut Ints(&mut v))?;
        g.arrays.push((name, v));
        Ok(())
    } else {
        legacy_values(s, n, ty, bin, &name, &mut Skip)
    }
}

/// Whether the next word is `w` (without reading it).
fn next_is<R: Read + Seek>(s: &mut Stream<R>, w: &[u8]) -> Result<bool, String> {
    let d = s.ahead(256)?;
    let k = d
        .iter()
        .position(|c| !c.is_ascii_whitespace())
        .unwrap_or(d.len());
    let rest = &d[k..];
    Ok(rest.len() >= w.len()
        && rest[..w.len()].eq_ignore_ascii_case(w)
        && rest.get(w.len()).is_none_or(|c| c.is_ascii_whitespace()))
}

// ---------------------------------------------------------------------
// XML files

/// A tag of an XML file.
enum Tag {
    Start {
        name: String,
        attrs: HashMap<String, String>,
        empty: bool,
    },
    End(String),
}

fn parse_attrs(s: &str) -> Result<HashMap<String, String>, String> {
    let mut m = HashMap::new();
    let mut rest = s.trim();
    while !rest.is_empty() {
        let eq = rest
            .find('=')
            .ok_or_else(|| format!("a bad attribute near '{rest}'"))?;
        let key = rest[..eq].trim().to_string();
        let after = rest[eq + 1..].trim_start();
        let q = after
            .chars()
            .next()
            .filter(|&c| c == '"' || c == '\'')
            .ok_or_else(|| format!("the attribute {key} has no quotes"))?;
        let end = after[1..]
            .find(q)
            .ok_or_else(|| format!("the attribute {key} is not closed"))?;
        m.insert(key, after[1..1 + end].to_string());
        rest = after[end + 2..].trim_start();
    }
    Ok(m)
}

fn parse_tag(inner: &str) -> Result<Tag, String> {
    if let Some(name) = inner.strip_prefix('/') {
        return Ok(Tag::End(name.trim().to_string()));
    }
    let empty = inner.ends_with('/');
    let inner = inner.trim_end_matches('/');
    let (name, attrs) = match inner.find(|c: char| c.is_ascii_whitespace()) {
        Some(k) => (&inner[..k], &inner[k..]),
        None => (inner, ""),
    };
    Ok(Tag::Start {
        name: name.to_string(),
        attrs: parse_attrs(attrs)?,
        empty,
    })
}

impl<R: Read + Seek> Stream<R> {
    /// The next tag, past text, comments and declarations; `None` at the
    /// end of the file.
    fn tag(&mut self) -> Result<Option<Tag>, String> {
        loop {
            self.skip_while(|c| c != b'<')?;
            let head = self.ahead(4)?;
            if head.is_empty() {
                return Ok(None);
            }
            if head.starts_with(b"<!--") {
                self.consume(4);
                self.skip_past(b"-->", "a comment is not closed")?;
                continue;
            }
            if head.starts_with(b"<?") {
                self.consume(2);
                self.skip_past(b"?>", "a declaration is not closed")?;
                continue;
            }
            if head.starts_with(b"<!") {
                self.consume(2);
                self.skip_past(b">", "a declaration is not closed")?;
                continue;
            }
            self.consume(1);
            // up to the end of the tag, outside quotes
            let mut t: Vec<u8> = Vec::new();
            let mut q: Option<u8> = None;
            loop {
                let d = self.data()?;
                if d.is_empty() {
                    return Err("a tag is not closed".into());
                }
                let mut k = 0;
                let mut closed = false;
                for &c in d {
                    match (q, c) {
                        (None, b'"' | b'\'') => q = Some(c),
                        (Some(x), y) if x == y => q = None,
                        (None, b'>') => {
                            closed = true;
                            break;
                        }
                        _ => {}
                    }
                    k += 1;
                }
                t.extend_from_slice(&d[..k]);
                self.consume(k + usize::from(closed));
                if closed {
                    break;
                }
                if t.len() > MAX_TAG {
                    return Err("a tag is not closed".into());
                }
            }
            let s = std::str::from_utf8(&t).map_err(|_| "the XML is not UTF-8 text")?;
            return parse_tag(s).map(Some);
        }
    }
}

/// The settings of the binary data of an XML file.
struct Xml {
    /// Header integer size in bytes (4 or 8).
    hl: usize,
    compressed: bool,
    big: bool,
}

impl Xml {
    /// `k` header integers.
    fn ints<Rd: Read>(&self, rd: &mut Rd, k: usize) -> Result<Vec<u64>, String> {
        let mut b = vec![0u8; k * self.hl];
        fill(rd, &mut b)?;
        Ok(b.chunks_exact(self.hl)
            .map(|c| {
                let mut a = [0u8; 8];
                a[..c.len()].copy_from_slice(c);
                if self.big {
                    a[..c.len()].reverse();
                }
                u64::from_le_bytes(a)
            })
            .collect())
    }

    /// One binary array: `rd` starts at its header; `len` is the size of
    /// the file.
    fn payload<Rd: Read, S: Sink>(
        &self,
        rd: &mut Rd,
        ty: Ty,
        sink: &mut S,
        len: u64,
    ) -> Result<(), String> {
        let mut feed = Feed::new(ty, self.big);
        if !self.compressed {
            let n = self.ints(rd, 1)?[0];
            if n > len {
                return Err("the binary data ends too early".into());
            }
            sink.reserve((n / ty.size() as u64) as usize);
            return binary(rd, n, &mut feed, sink);
        }
        // blocks: number, block size, size of the last block, the
        // compressed sizes, then the compressed blocks
        let h = self.ints(rd, 3)?;
        let (nb, bs, last) = (h[0], h[1], h[2]);
        let nb = declared(nb as usize, len as usize, "compressed blocks")?;
        let sizes = self.ints(rd, nb)?;
        for (k, &sz) in sizes.iter().enumerate() {
            let sz = declared(sz as usize, len as usize, "a compressed block")?;
            let mut blk = vec![0u8; sz];
            fill(rd, &mut blk).map_err(|_| "the compressed data ends too early".to_string())?;
            let d = inflate::zlib(&blk)?;
            let want = if k + 1 == nb && last != 0 { last } else { bs };
            if d.len() as u64 != want {
                return Err(format!(
                    "a compressed block has {} bytes, expected {want}",
                    d.len()
                ));
            }
            feed.push(&d, sink)?;
        }
        Ok(())
    }
}

/// What an array of the file is for.
#[derive(Clone)]
enum Role {
    Points,
    Conn,
    Offsets,
    Types,
    Faces,
    FaceOffsets,
    Region(String),
    Names(String),
    Skip,
}

/// Which roles are taken (the first array of a name counts).
#[derive(Default)]
struct Seen {
    points: bool,
    conn: bool,
    offsets: bool,
    types: bool,
    faces: bool,
    faceoffsets: bool,
}

fn role(stack: &[String], attrs: &HashMap<String, String>, seen: &mut Seen) -> Role {
    let at = |k: usize| {
        stack
            .len()
            .checked_sub(k)
            .and_then(|i| stack.get(i))
            .map(String::as_str)
    };
    let name = attrs.get("Name").map(String::as_str).unwrap_or("");
    let take = |flag: &mut bool, r: Role| {
        if *flag {
            Role::Skip
        } else {
            *flag = true;
            r
        }
    };
    match (at(2), at(1)) {
        (Some("Piece"), Some("Points")) => take(&mut seen.points, Role::Points),
        (Some("Piece"), Some("Cells")) => match name {
            "connectivity" => take(&mut seen.conn, Role::Conn),
            "offsets" => take(&mut seen.offsets, Role::Offsets),
            "types" => take(&mut seen.types, Role::Types),
            "faces" => take(&mut seen.faces, Role::Faces),
            "faceoffsets" => take(&mut seen.faceoffsets, Role::FaceOffsets),
            _ => Role::Skip,
        },
        (Some("Piece"), Some("CellData"))
            if attrs
                .get("NumberOfComponents")
                .is_none_or(|c| c.as_str() == "1")
                && is_region(name) =>
        {
            Role::Region(name.to_string())
        }
        (_, Some("FieldData")) if attrs.get("type").map(String::as_str) == Some("String") => {
            Role::Names(name.to_string())
        }
        _ => Role::Skip,
    }
}

/// The source of the numbers of one array.
enum Src<'a, R> {
    /// ASCII text up to the next `<`.
    Ascii(&'a mut Stream<R>),
    /// Binary data, from the stream as it stands (raw) or through base64.
    Binary {
        s: &'a mut Stream<R>,
        base64: bool,
        x: &'a Xml,
        ty: Ty,
    },
}

impl<R: Read + Seek> Src<'_, R> {
    fn run<S: Sink>(self, sink: &mut S, what: &str) -> Result<(), String> {
        match self {
            Src::Ascii(s) => s.numbers(None, sink, what).map(|_| ()),
            Src::Binary { s, base64, x, ty } => {
                let len = s.len;
                if base64 {
                    x.payload(&mut B64::new(s), ty, sink, len)
                } else {
                    x.payload(s, ty, sink, len)
                }
            }
        }
        .map_err(|e| {
            if e.starts_with(what) {
                e
            } else {
                format!("{what}: {e}")
            }
        })
    }
}

/// Reads one array into its place in the grid.
fn fill_role<R: Read + Seek>(
    g: &mut Grid,
    role: &Role,
    src: Src<'_, R>,
    what: &str,
) -> Result<(), String> {
    match role {
        Role::Points => src.run(&mut Points::new(&mut g.points), what),
        Role::Conn => src.run(
            &mut Ids {
                out: &mut g.conn,
                what: "connectivity",
            },
            what,
        ),
        Role::Offsets => src.run(&mut Ints(&mut g.offsets), what),
        Role::Types => src.run(&mut Types(&mut g.types), what),
        Role::Faces => src.run(
            &mut Ids {
                out: &mut g.faces,
                what: "faces",
            },
            what,
        ),
        Role::FaceOffsets => src.run(&mut Ints(&mut g.faceoffsets), what),
        Role::Region(name) => {
            let mut v = Vec::new();
            src.run(&mut Ints(&mut v), what)?;
            g.arrays.push((name.clone(), v));
            Ok(())
        }
        Role::Names(name) => {
            let mut codes = Vec::new();
            src.run(&mut Codes(&mut codes), what)?;
            let mut list = Vec::new();
            let mut cur = String::new();
            for c in codes {
                if c == 0 {
                    list.push(std::mem::take(&mut cur));
                } else {
                    cur.push(c as char);
                }
            }
            g.names.insert(name.clone(), list);
            Ok(())
        }
        Role::Skip => src.run(&mut Skip, what),
    }
}

/// An array of the appended data, read after the XML.
struct Later {
    offset: u64,
    role: Role,
    ty: Ty,
    name: String,
}

/// The checks of the root element; the settings of its binary data.
fn root(name: &str, attrs: &HashMap<String, String>) -> Result<Xml, String> {
    if name != "VTKFile" {
        return Err("not a VTKFile".into());
    }
    let ty = attrs.get("type").map(String::as_str).unwrap_or("");
    if ty != "UnstructuredGrid" {
        return Err(format!(
            "the file is a {ty}: only UnstructuredGrid is read (a parallel .pvtu: read its pieces)"
        ));
    }
    let compressor = attrs.get("compressor");
    if let Some(c) = compressor
        && c != "vtkZLibDataCompressor"
    {
        return Err(format!(
            "the compressor {c} is not read (zlib is; save the file without compression or with zlib)"
        ));
    }
    Ok(Xml {
        hl: if attrs.get("header_type").map(String::as_str) == Some("UInt64") {
            8
        } else {
            4
        },
        compressed: compressor.is_some(),
        big: attrs.get("byte_order").map(String::as_str) == Some("BigEndian"),
    })
}

fn parse_xml<R: Read + Seek>(s: &mut Stream<R>) -> Result<(Grid, String), String> {
    let mut g = Grid::default();
    let mut xml: Option<Xml> = None;
    let mut stack: Vec<String> = Vec::new();
    let mut later: Vec<Later> = Vec::new();
    let mut seen = Seen::default();
    let (mut grids, mut pieces, mut points_el, mut cells_el) = (0usize, 0usize, false, false);
    let mut points_format: Option<String> = None;
    let mut appended: Option<(u64, bool)> = None;
    while let Some(tag) = s.tag()? {
        let (name, attrs, empty) = match tag {
            Tag::End(name) => {
                let top = stack.pop().ok_or("an end tag without a start")?;
                if top != name {
                    return Err(format!("the end tag </{name}> does not match"));
                }
                if stack.is_empty() {
                    // only the first root element is read
                    break;
                }
                continue;
            }
            Tag::Start { name, attrs, empty } => (name, attrs, empty),
        };
        let Some(x) = &xml else {
            xml = Some(root(&name, &attrs)?);
            if empty {
                break;
            }
            stack.push(name);
            continue;
        };
        match name.as_str() {
            "AppendedData" => {
                // binary data follow, read below at the offsets of the arrays
                let b64 = attrs.get("encoding").is_some_and(|e| e.contains("base64"));
                s.skip_while(|c| c.is_ascii_whitespace())?;
                if s.peek()? == Some(b'_') {
                    s.consume(1);
                }
                appended = Some((s.pos, b64));
                break;
            }
            "UnstructuredGrid" if stack.len() == 1 => grids += 1,
            "Piece" => {
                pieces += 1;
                if pieces > 1 {
                    return Err(format!("{pieces} pieces: only a single piece is read"));
                }
                // room for what the piece declares, within the size of the file
                let count = |k: &str| {
                    attrs
                        .get(k)
                        .and_then(|v| v.parse::<usize>().ok())
                        .filter(|&n| n as u64 <= s.len)
                        .unwrap_or(0)
                };
                g.points.reserve(count("NumberOfPoints"));
                let nc = count("NumberOfCells");
                g.types.reserve(nc);
                g.offsets.reserve(nc);
            }
            "Points" => points_el = true,
            "Cells" => cells_el = true,
            "DataArray" => {
                let role = role(&stack, &attrs, &mut seen);
                if matches!(role, Role::Points) {
                    points_format = attrs.get("format").cloned();
                }
                if !matches!(role, Role::Skip) {
                    let what = format!(
                        "array {}",
                        attrs.get("Name").map(String::as_str).unwrap_or("?")
                    );
                    // string arrays hold UInt8 codes
                    let ty = if matches!(role, Role::Names(_)) {
                        Ty::U8
                    } else {
                        Ty::xml(attrs.get("type").map(String::as_str).unwrap_or("Float32"))?
                    };
                    match attrs.get("format").map(String::as_str).unwrap_or("ascii") {
                        "ascii" if !empty => fill_role(&mut g, &role, Src::Ascii(s), &what)?,
                        "binary" if !empty => {
                            let src = Src::Binary {
                                s,
                                base64: true,
                                x,
                                ty,
                            };
                            fill_role(&mut g, &role, src, &what)?;
                        }
                        "ascii" => {}
                        "binary" => return Err(format!("{what}: the binary data ends too early")),
                        "appended" => {
                            let offset = attrs
                                .get("offset")
                                .and_then(|o| o.trim().parse().ok())
                                .ok_or_else(|| format!("{what}: no offset"))?;
                            later.push(Later {
                                offset,
                                role,
                                ty,
                                name: what,
                            });
                        }
                        other => return Err(format!("{what}: the format '{other}' is not read")),
                    }
                }
            }
            _ => {}
        }
        if !empty {
            stack.push(name);
        }
    }
    let x = xml.ok_or("the file has no elements")?;
    match appended {
        Some((base, b64)) => {
            later.sort_by_key(|l| l.offset);
            for l in later {
                let at = base
                    .checked_add(l.offset)
                    .filter(|&a| a < s.len)
                    .ok_or_else(|| format!("{}: the offset is beyond the data", l.name))?;
                s.seek(at)?;
                let src = Src::Binary {
                    s: &mut *s,
                    base64: b64,
                    x: &x,
                    ty: l.ty,
                };
                fill_role(&mut g, &l.role, src, &l.name)?;
            }
        }
        None => {
            if !stack.is_empty() {
                return Err("an element is not closed".into());
            }
            if let Some(l) = later.first() {
                return Err(format!("{}: the offset is beyond the data", l.name));
            }
        }
    }
    if grids == 0 {
        return Err("no UnstructuredGrid".into());
    }
    if pieces != 1 {
        return Err(format!("{pieces} pieces: only a single piece is read"));
    }
    if !points_el || !seen.points {
        return Err("the piece has no Points".into());
    }
    if !cells_el {
        return Err("the piece has no Cells".into());
    }
    for (ok, n) in [
        (seen.conn, "connectivity"),
        (seen.offsets, "offsets"),
        (seen.types, "types"),
    ] {
        if !ok {
            return Err(format!("Cells has no {n}"));
        }
    }
    if g.offsets.len() != g.types.len() {
        return Err(format!(
            "{} offsets for {} cell types",
            g.offsets.len(),
            g.types.len()
        ));
    }
    let enc = match (x.compressed, points_format) {
        (true, Some(f)) => format!("{f}, zlib"),
        (false, Some(f)) => f,
        _ => String::new(),
    };
    Ok((g, format!("VTK XML ({enc})")))
}

fn find(h: &[u8], n: &[u8]) -> Option<usize> {
    h.windows(n.len()).position(|w| w == n)
}

// ---------------------------------------------------------------------
// writing

/// Writes `mesh` as a VTK XML unstructured grid (`.vtu`): every cell a
/// polyhedron, every boundary face a polygon, the cell-data array `region`
/// (cell zone of a cell, boundary of a face) and the names of the zones and
/// boundaries as field data. `binary` selects appended binary data (little
/// endian, uncompressed) instead of ASCII.
pub fn write_vtu(mesh: &Mesh, path: &Path, binary: bool) -> Result<(), String> {
    let n_cells = mesh.n_cells;
    // faces of each cell, outward
    let mut cell_faces: Vec<Vec<Vec<u32>>> = vec![Vec::new(); n_cells];
    let nf = mesh.n_faces();
    let ni = mesh.n_interior_faces();
    let nodes =
        |f: usize| &mesh.face_nodes[mesh.face_ptr[f] as usize..mesh.face_ptr[f + 1] as usize];
    for f in 0..nf {
        let n = nodes(f).to_vec();
        if f < ni {
            let mut r = n.clone();
            r.reverse();
            cell_faces[mesh.neighbour[f] as usize].push(r);
        }
        cell_faces[mesh.owner[f] as usize].push(n);
    }
    let mut conn: Vec<i64> = Vec::new();
    let mut offs: Vec<i64> = Vec::new();
    let mut types: Vec<u8> = Vec::new();
    let mut faces: Vec<i64> = Vec::new();
    let mut faceoffs: Vec<i64> = Vec::new();
    let mut region: Vec<i32> = Vec::new();
    for (c, fs) in cell_faces.iter().enumerate() {
        let mut uniq: Vec<u32> = fs.iter().flatten().copied().collect();
        uniq.sort_unstable();
        uniq.dedup();
        conn.extend(uniq.iter().map(|&v| v as i64));
        offs.push(conn.len() as i64);
        types.push(42);
        faces.push(fs.len() as i64);
        for f in fs {
            faces.push(f.len() as i64);
            faces.extend(f.iter().map(|&v| v as i64));
        }
        faceoffs.push(faces.len() as i64);
        region.push(mesh.cell_zone[c] as i32);
    }
    for (bi, b) in mesh.boundaries.iter().enumerate() {
        for f in b.range() {
            let n = nodes(f);
            conn.extend(n.iter().map(|&v| v as i64));
            offs.push(conn.len() as i64);
            types.push(match n.len() {
                3 => 5,
                4 => 9,
                _ => 7,
            });
            faceoffs.push(-1);
            region.push(bi as i32);
        }
    }
    let n_out = types.len();
    let mut pts: Vec<f64> = Vec::with_capacity(3 * mesh.points.len());
    for p in &mesh.points {
        pts.extend_from_slice(p);
    }
    let names = |list: Vec<String>| -> Vec<u8> {
        let mut v = Vec::new();
        for (i, n) in list.iter().enumerate() {
            v.extend(format!("{i}={n}").bytes());
            v.push(0);
        }
        v
    };
    let zone_names = names(mesh.zones.clone());
    let bound_names = names(mesh.boundaries.iter().map(|b| b.name.clone()).collect());

    enum Data<'a> {
        F64(&'a [f64]),
        I64(&'a [i64]),
        I32(&'a [i32]),
        U8(&'a [u8]),
    }
    impl Data<'_> {
        fn ty(&self) -> &'static str {
            match self {
                Data::F64(_) => "Float64",
                Data::I64(_) => "Int64",
                Data::I32(_) => "Int32",
                Data::U8(_) => "UInt8",
            }
        }
        fn len(&self) -> usize {
            match self {
                Data::F64(v) => v.len(),
                Data::I64(v) => v.len(),
                Data::I32(v) => v.len(),
                Data::U8(v) => v.len(),
            }
        }
        fn bytes(&self) -> Vec<u8> {
            match self {
                Data::F64(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::I64(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::I32(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
                Data::U8(v) => v.to_vec(),
            }
        }
        fn ascii(&self) -> String {
            let mut s = String::new();
            let mut put = |x: String| {
                s.push_str(&x);
                s.push(' ');
            };
            match self {
                Data::F64(v) => v.iter().for_each(|x| put(format!("{x:e}"))),
                Data::I64(v) => v.iter().for_each(|x| put(x.to_string())),
                Data::I32(v) => v.iter().for_each(|x| put(x.to_string())),
                Data::U8(v) => v.iter().for_each(|x| put(x.to_string())),
            }
            s
        }
    }
    let mut appended: Vec<u8> = Vec::new();
    let mut array = |name: &str, comps: usize, d: Data| -> String {
        if binary {
            let off = appended.len();
            let b = d.bytes();
            appended.extend((b.len() as u64).to_le_bytes());
            appended.extend(b);
            format!(
                "<DataArray type=\"{}\" Name=\"{name}\" NumberOfComponents=\"{comps}\" \
                 format=\"appended\" offset=\"{off}\"/>\n",
                d.ty()
            )
        } else {
            format!(
                "<DataArray type=\"{}\" Name=\"{name}\" NumberOfComponents=\"{comps}\" \
                 format=\"ascii\">\n{}\n</DataArray>\n",
                d.ty(),
                d.ascii()
            )
        }
    };
    let string_array = |name: &str, bytes: &[u8]| -> String {
        let codes = Data::U8(bytes).ascii();
        format!(
            "<DataArray type=\"String\" Name=\"{name}\" NumberOfTuples=\"1\" format=\"ascii\">\
             {codes}</DataArray>\n"
        )
    };
    let mut x = String::new();
    x.push_str(
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" \
         byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<UnstructuredGrid>\n",
    );
    x.push_str("<FieldData>\n");
    x.push_str(&string_array("zone_names", &zone_names));
    x.push_str(&string_array("boundary_names", &bound_names));
    x.push_str("</FieldData>\n");
    x.push_str(&format!(
        "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{n_out}\">\n<Points>\n",
        mesh.points.len()
    ));
    x.push_str(&array("Points", 3, Data::F64(&pts)));
    x.push_str("</Points>\n<Cells>\n");
    x.push_str(&array("connectivity", 1, Data::I64(&conn)));
    x.push_str(&array("offsets", 1, Data::I64(&offs)));
    x.push_str(&array("types", 1, Data::U8(&types)));
    x.push_str(&array("faces", 1, Data::I64(&faces)));
    x.push_str(&array("faceoffsets", 1, Data::I64(&faceoffs)));
    x.push_str("</Cells>\n<CellData>\n");
    x.push_str(&array("region", 1, Data::I32(&region)));
    x.push_str("</CellData>\n</Piece>\n</UnstructuredGrid>\n");
    let _ = Data::len;
    let mut out = x.into_bytes();
    if binary {
        out.extend_from_slice(b"<AppendedData encoding=\"raw\">\n_");
        out.extend(appended);
        out.extend_from_slice(b"\n</AppendedData>\n");
    }
    out.extend_from_slice(b"</VTKFile>\n");
    std::fs::write(path, out).map_err(|e| format!("{}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data(name: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    fn same(a: &Mesh, b: &Mesh, what: &str) {
        assert_eq!(a.n_cells, b.n_cells, "{what}");
        assert_eq!(a.points, b.points, "{what}");
        assert_eq!(a.face_ptr, b.face_ptr, "{what}");
        assert_eq!(a.face_nodes, b.face_nodes, "{what}");
        assert_eq!(a.owner, b.owner, "{what}");
        assert_eq!(a.neighbour, b.neighbour, "{what}");
        assert_eq!(a.cell_zone, b.cell_zone, "{what}");
        assert_eq!(a.zones, b.zones, "{what}");
        assert_eq!(a.boundaries, b.boundaries, "{what}");
    }

    /// Every encoding read through buffers of odd sizes, so that numbers,
    /// base64 groups, header integers, compressed blocks and tags fall
    /// across the ends of the buffer, gives the mesh of the whole file.
    #[test]
    fn small_buffers_read_the_same_mesh() {
        let dir = std::env::temp_dir().join(format!("vtk-stream-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut files: Vec<std::path::PathBuf> = [
            "channel_mixed_ascii.vtk",
            "channel_mixed_binary.vtk",
            "channel_mixed_ascii.vtu",
            "channel_mixed_binary.vtu",
            "channel_mixed_zlib.vtu",
            "channel_mixed_appended.vtu",
        ]
        .iter()
        .map(|f| data(&format!("vtk/{f}")))
        .collect();
        // appended raw data, polyhedra and a legacy file of polyhedra
        let poly = crate::format::read(data("openfoam/dual_polyhedra"), None)
            .unwrap()
            .mesh;
        for binary in [false, true] {
            let p = dir.join(format!("poly-{binary}.vtu"));
            write_vtu(&poly, &p, binary).unwrap();
            files.push(p);
        }
        for f in &files {
            let bytes = std::fs::read(f).unwrap();
            let whole = Vtk.parse(&bytes, &ReadOptions::default()).unwrap();
            let from_file = Vtk.read_path(f, &ReadOptions::default()).unwrap();
            same(&whole.mesh, &from_file.mesh, &format!("{}", f.display()));
            for cap in [64, 97, 1000, 4099] {
                let s = Stream::with_buffer(std::io::Cursor::new(&bytes), bytes.len() as u64, cap);
                let r = read(s).unwrap_or_else(|e| panic!("{} ({cap}): {e}", f.display()));
                same(&whole.mesh, &r.mesh, &format!("{} ({cap})", f.display()));
                assert_eq!(whole.format, r.format);
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn base64_pieces_decode_as_one_stream() {
        // "AQID" = 1 2 3, "BA==" = 4 (padded), then "BQY=" = 5 6, then
        // white space and the end of the element
        let text = b"AQID BA==\nBQY=  <";
        let mut s = Stream::with_buffer(std::io::Cursor::new(&text[..]), text.len() as u64, 64);
        let mut out = Vec::new();
        B64::new(&mut s).read_to_end(&mut out).unwrap();
        assert_eq!(out, [1, 2, 3, 4, 5, 6]);
        assert_eq!(s.peek().unwrap(), Some(b'<'));
    }

    #[test]
    fn numbers_split_across_buffers() {
        let text: String = (0..500)
            .map(|i| format!("{} ", i * 1001))
            .collect::<String>()
            + "<";
        let b = text.as_bytes();
        let mut s = Stream::with_buffer(std::io::Cursor::new(b), b.len() as u64, 64);
        let mut v = Vec::new();
        let n = s.numbers(None, &mut Ints(&mut v), "test").unwrap();
        assert_eq!(n, 500);
        assert_eq!(v, (0..500).map(|i| i * 1001).collect::<Vec<i64>>());
    }
}
