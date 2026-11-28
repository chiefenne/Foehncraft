//! STL surfaces (triangles), for the datums of reports (D262) and later
//! the surface queries of D108.
//!
//! Both forms of the format are read: ASCII (`solid … facet normal …
//! outer loop vertex x y z … endloop endfacet … endsolid`) and binary (an
//! 80-byte header, a 32-bit little-endian triangle count, then 50 bytes per
//! triangle: normal, three corners as 32-bit floats, a 16-bit attribute).
//! A file is binary when its size is exactly 84 + 50 n for the count n it
//! declares (an ASCII file may begin with `solid` as well, the test of the
//! size decides). The stored normals are not used: a triangle's direction
//! is that of its corners (right-hand rule), as most writers keep both
//! consistent and the corners are what the file is about.
//!
//! Every broken input ends in an error message, never a panic (D202):
//! counts are checked against the size of the file before allocating,
//! numbers must parse and be finite.

/// A triangle: three corners (m, after the scale of the caller).
pub type Triangle = [[f64; 3]; 3];

/// The triangles of an STL file with the tag of every triangle: the index
/// of its `solid` in an ASCII file (named in `solids`), the 16-bit
/// attribute of a binary file (`solids` empty). The tags name the surface
/// patches of an immersed body (L13 stage 2e).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Tagged {
    /// The triangles, in file order.
    pub triangles: Vec<Triangle>,
    /// The tag of every triangle.
    pub tags: Vec<u32>,
    /// The names of the solids of an ASCII file, in file order (an unnamed
    /// solid has an empty name); empty for a binary file.
    pub solids: Vec<String>,
}

/// Reads the triangles of an STL file's bytes.
pub fn parse(bytes: &[u8]) -> Result<Vec<Triangle>, String> {
    parse_tagged(bytes).map(|t| t.triangles)
}

/// Reads the triangles of an STL file's bytes with their tags.
pub fn parse_tagged(bytes: &[u8]) -> Result<Tagged, String> {
    if bytes.len() >= 84 {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        let size = n.checked_mul(50).and_then(|x| x.checked_add(84));
        if size == Some(bytes.len()) {
            return binary(bytes, n);
        }
    }
    let text = std::str::from_utf8(bytes).map_err(|_| {
        "STL: neither a binary file of the declared size nor ASCII text".to_string()
    })?;
    ascii(text)
}

/// Reads an STL file.
pub fn read(path: &std::path::Path) -> Result<Vec<Triangle>, String> {
    read_tagged(path).map(|t| t.triangles)
}

/// Reads an STL file with the tags of its triangles.
pub fn read_tagged(path: &std::path::Path) -> Result<Tagged, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse_tagged(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

fn binary(bytes: &[u8], n: usize) -> Result<Tagged, String> {
    let n = super::declared(n, bytes.len(), "STL triangles")?;
    if n == 0 {
        return Err("STL: the file holds no triangle".into());
    }
    let mut out = Vec::with_capacity(n);
    let mut tags = Vec::with_capacity(n);
    for i in 0..n {
        let rec = 84 + 50 * i;
        let f = |k: usize| -> Result<f64, String> {
            let o = rec + 12 + 4 * k;
            let b = bytes
                .get(o..o + 4)
                .ok_or_else(|| format!("STL: triangle {i} is cut off"))?;
            let v = f32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64;
            if v.is_finite() {
                Ok(v)
            } else {
                Err(format!(
                    "STL: triangle {i} has a coordinate that is not finite"
                ))
            }
        };
        let mut t = [[0.0; 3]; 3];
        for (c, corner) in t.iter_mut().enumerate() {
            for (k, x) in corner.iter_mut().enumerate() {
                *x = f(3 * c + k)?;
            }
        }
        out.push(t);
        let a = bytes
            .get(rec + 48..rec + 50)
            .ok_or_else(|| format!("STL: triangle {i} is cut off"))?;
        tags.push(u16::from_le_bytes([a[0], a[1]]) as u32);
    }
    Ok(Tagged {
        triangles: out,
        tags,
        solids: Vec::new(),
    })
}

/// The name of a solid: the words after `solid` up to its first facet (or
/// its end).
fn solid_name<'a>(words: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>) -> String {
    let mut name: Vec<&str> = Vec::new();
    while let Some(w) = words.peek() {
        if w.eq_ignore_ascii_case("facet") || w.eq_ignore_ascii_case("endsolid") {
            break;
        }
        name.push(w);
        words.next();
    }
    name.join(" ")
}

fn ascii(text: &str) -> Result<Tagged, String> {
    let mut words = text.split_ascii_whitespace().peekable();
    match words.next() {
        Some(w) if w.eq_ignore_ascii_case("solid") => {}
        _ => return Err("STL: an ASCII file must begin with `solid`".into()),
    }
    let mut solids = vec![solid_name(&mut words)];
    let mut out = Vec::new();
    let mut tags = Vec::new();
    let mut corners: Vec<[f64; 3]> = Vec::with_capacity(3);
    let mut in_loop = false;
    let mut ended = false;
    while let Some(w) = words.next() {
        match w.to_ascii_lowercase().as_str() {
            "facet" | "normal" => {}
            "endsolid" => {
                if in_loop {
                    return Err("STL: `endsolid` inside a loop".into());
                }
                ended = true;
                // the name after endsolid, up to a further solid
                while let Some(n) = words.peek() {
                    if n.eq_ignore_ascii_case("solid") {
                        break;
                    }
                    words.next();
                }
            }
            "solid" if ended && !in_loop => {
                ended = false;
                solids.push(solid_name(&mut words));
            }
            "outer" => {
                if words.next().map(|x| x.to_ascii_lowercase()) != Some("loop".into()) {
                    return Err("STL: `outer` without `loop`".into());
                }
                if in_loop {
                    return Err("STL: a loop inside a loop".into());
                }
                in_loop = true;
                corners.clear();
            }
            "vertex" => {
                if !in_loop {
                    return Err("STL: a vertex outside a loop".into());
                }
                let mut v = [0.0; 3];
                for x in &mut v {
                    let s = words
                        .next()
                        .ok_or_else(|| "STL: a vertex with fewer than three numbers".to_string())?;
                    *x = s
                        .parse::<f64>()
                        .ok()
                        .filter(|y| y.is_finite())
                        .ok_or_else(|| format!("STL: `{s}` is not a finite number"))?;
                }
                if corners.len() == 3 {
                    return Err("STL: a facet with more than three vertices".into());
                }
                corners.push(v);
            }
            "endloop" => {
                if !in_loop || corners.len() != 3 {
                    return Err("STL: a loop without exactly three vertices".into());
                }
                in_loop = false;
                out.push([corners[0], corners[1], corners[2]]);
                tags.push((solids.len() - 1) as u32);
            }
            "endfacet" => {
                if in_loop {
                    return Err("STL: `endfacet` inside a loop".into());
                }
            }
            _ => {
                // the numbers of `facet normal` and the names of solids
                if !(w.parse::<f64>().is_ok() || !in_loop) {
                    return Err(format!("STL: unexpected `{w}`"));
                }
            }
        }
    }
    if in_loop {
        return Err("STL: the file ends inside a loop".into());
    }
    if !ended {
        return Err("STL: the file ends without `endsolid`".into());
    }
    if out.is_empty() {
        return Err("STL: the file holds no triangle".into());
    }
    Ok(Tagged {
        triangles: out,
        tags,
        solids,
    })
}

/// The ASCII STL text of `tris` (for tests and examples).
pub fn to_ascii(name: &str, tris: &[Triangle]) -> String {
    let mut s = format!("solid {name}\n");
    for t in tris {
        let a = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
        let b = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
        let n = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
        s.push_str(&format!(
            "  facet normal {:e} {:e} {:e}\n    outer loop\n",
            n[0] / l,
            n[1] / l,
            n[2] / l
        ));
        for v in t {
            s.push_str(&format!("      vertex {:e} {:e} {:e}\n", v[0], v[1], v[2]));
        }
        s.push_str("    endloop\n  endfacet\n");
    }
    s.push_str(&format!("endsolid {name}\n"));
    s
}

/// The binary STL bytes of `tris` (for tests).
pub fn to_binary(tris: &[Triangle]) -> Vec<u8> {
    let mut b = vec![0u8; 80];
    b.extend_from_slice(&(tris.len() as u32).to_le_bytes());
    for t in tris {
        b.extend_from_slice(&[0u8; 12]);
        for v in t {
            for x in v {
                b.extend_from_slice(&(*x as f32).to_le_bytes());
            }
        }
        b.extend_from_slice(&[0u8; 2]);
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Vec<Triangle> {
        vec![
            [[0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [1.0, 1.0, 0.5]],
            [[0.0, 0.0, 0.5], [1.0, 1.0, 0.5], [0.0, 1.0, 0.5]],
        ]
    }

    #[test]
    fn both_forms_read_back() {
        let t = square();
        assert_eq!(parse(to_ascii("sq", &t).as_bytes()).unwrap(), t);
        assert_eq!(parse(&to_binary(&t)).unwrap(), t);
    }

    #[test]
    fn broken_files_are_refused() {
        let a = to_ascii("sq", &square());
        for bad in [
            "",
            "solid x\n",
            &a[..a.len() / 2],
            &a.replace("vertex 1e0", "vertex nan"),
            &a.replace("endloop", ""),
        ] {
            assert!(parse(bad.as_bytes()).is_err(), "{bad:?}");
        }
        let mut b = to_binary(&square());
        b[80] = 200; // a count the file cannot hold
        assert!(parse(&b).is_err());
    }

    #[test]
    fn solids_and_attributes_tag_the_triangles() {
        let t = square();
        // three solids, the last one unnamed, with and without names after
        // endsolid
        let text = to_ascii("inlet face", &t[..1])
            + &to_ascii("wall", &t[1..]).replace("endsolid wall", "endsolid")
            + &to_ascii("", &t);
        let g = parse_tagged(text.as_bytes()).unwrap();
        assert_eq!(g.solids, ["inlet face", "wall", ""]);
        assert_eq!(g.tags, [0, 1, 2, 2]);
        assert_eq!(g.triangles.len(), 4);
        // the 16-bit attribute of a binary file
        let mut b = to_binary(&t);
        b[84 + 48] = 7;
        b[84 + 50 + 49] = 1;
        let g = parse_tagged(&b).unwrap();
        assert!(g.solids.is_empty());
        assert_eq!(g.tags, [7, 256]);
    }
}
