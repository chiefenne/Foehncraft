//! Output of a radiation call (R1a): the fields of a [`Solution`] as a
//! VTKHDF file for ParaView (feature `hdf5`, D41) and as CSV tables for
//! scripts.
//!
//! The cell fields are the net source (absorbed − emitted, W/m³) and its
//! standard error; the boundary blocks carry the incident and the net flux
//! into the wall (W/m²) and the error of the net flux. Both writers take
//! the cells and boundary faces of the mesh they are given, so on ranks
//! each rank writes its own part (the parallel writer into one file comes
//! with the solver stages, as for the flow).

use crate::Solution;
use foehncraft_mesh::geometry::norm;
use foehncraft_mesh::{Geometry, Mesh};
use std::fmt::Write;

/// The cells as CSV: centroid, volume, net source and its error (W/m³),
/// absorbed and emitted power (W), then the columns `inputs` (such as the
/// temperature of the call), one row per cell.
pub fn cell_csv(mesh: &Mesh, geo: &Geometry, sol: &Solution, inputs: &[(&str, &[f64])]) -> String {
    let mut s = String::from("cell,x,y,z,volume,source,source_error,absorbed,emitted");
    for (name, _) in inputs {
        let _ = write!(s, ",{name}");
    }
    s.push('\n');
    for c in 0..mesh.n_cells {
        let [x, y, z] = geo.cell_centroid[c];
        let _ = write!(
            s,
            "{c},{x:e},{y:e},{z:e},{:e},{:e},{:e},{:e},{:e}",
            geo.cell_volume[c],
            sol.source[c],
            sol.source_error[c],
            sol.absorbed[c] as f64 * sol.quantum,
            sol.emitted[c] as f64 * sol.quantum,
        );
        for (_, v) in inputs {
            let _ = write!(s, ",{:e}", v[c]);
        }
        s.push('\n');
    }
    s
}

/// The boundary faces as CSV: boundary, centroid, area, incident flux, net
/// flux into the wall and its error (W/m²), one row per boundary face.
pub fn wall_csv(mesh: &Mesh, geo: &Geometry, sol: &Solution) -> String {
    let mut s = String::from("boundary,face,x,y,z,area,incident,net,net_error\n");
    let first = mesh.n_interior_faces();
    for b in &mesh.boundaries {
        for f in b.range() {
            let k = f - first;
            let [x, y, z] = geo.face_centroid[f];
            let _ = writeln!(
                s,
                "{},{},{x:e},{y:e},{z:e},{:e},{:e},{:e},{:e}",
                b.name,
                mesh.face_id[f],
                norm(geo.face_area[f]),
                sol.wall_incident[k],
                sol.wall_net[k],
                sol.wall_net_error[k],
            );
        }
    }
    s
}

/// Writes the mesh with the fields of `sol` (and the cell fields `inputs`,
/// such as the temperature and the absorption coefficient of the call) to
/// a VTKHDF file, with the energy account as attributes of the file.
#[cfg(feature = "hdf5")]
pub fn write_vtkhdf(
    path: impl AsRef<std::path::Path>,
    mesh: &Mesh,
    sol: &Solution,
    inputs: &[foehncraft_mesh::vtkhdf::CellField],
) -> anyhow::Result<()> {
    use foehncraft_mesh::vtkhdf::{Attr, CellField, FaceField, write_all};
    let mut cells = vec![
        CellField {
            name: "radiative source",
            values: &sol.source,
        },
        CellField {
            name: "radiative source error",
            values: &sol.source_error,
        },
    ];
    cells.extend(inputs.iter().map(|f| CellField {
        name: f.name,
        values: f.values,
    }));
    let faces = [
        FaceField {
            name: "incident radiative flux",
            values: &sol.wall_incident,
        },
        FaceField {
            name: "net radiative flux",
            values: &sol.wall_net,
        },
        FaceField {
            name: "net radiative flux error",
            values: &sol.wall_net_error,
        },
    ];
    let l = &sol.ledger;
    let q = sol.quantum;
    write_all(
        path,
        mesh,
        &cells,
        &faces,
        &[
            ("format", Attr::Str("foehncraft-radiation")),
            (
                "units",
                Attr::Str("SI: W/m3 in the cells, W/m2 on the walls"),
            ),
            ("bundles", Attr::U64(sol.bundles as u64)),
            ("quantum", Attr::F64(q)),
            ("emitted", Attr::F64(l.emitted as f64 * q)),
            ("absorbed_gas", Attr::F64(l.absorbed_gas as f64 * q)),
            ("absorbed_walls", Attr::F64(l.absorbed_walls as f64 * q)),
            ("escaped", Attr::F64(l.escaped as f64 * q)),
            ("roulette", Attr::F64(l.roulette as f64 * q)),
        ],
    )
}
