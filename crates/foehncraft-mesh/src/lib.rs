//! The mesh of Foehncraft (stage 1).
//!
//! - [`mesh`]: the polyhedral face-based mesh (D23) in the layout the
//!   kernels of `foehncraft-kernels` use;
//! - [`builder`]: assembling a mesh from faces or from cells;
//! - [`geometry`]: face area vectors and centroids, cell volumes and
//!   centroids;
//! - [`quality`]: the mesh quality report;
//! - [`generate`]: the built-in generator for canonical test meshes (D34);
//! - [`cylinder`]: the built-in channel with a cylinder (Schäfer & Turek);
//! - [`annulus`]: the built-in annulus (Taylor–Couette, rotor–stator);
//! - [`revolved`]: the built-in mesh of revolution (rotor and stator with an
//!   interface of cylinders, discs and cones, D186);
//! - [`imprint`]: the sliding interface of a rotor and a stator with different
//!   pitch, imprinted at the current rotor angle (D58, D182);
//! - [`wedge`]: axisymmetric wedge meshes revolved from a 2D mesh;
//! - [`surface`]: triangle surfaces of immersed bodies (L13 stage 0, D225):
//!   check, winding numbers, signed distance, crossings;
//! - [`format`](mod@format): the mesh file formats behind one reader interface (D33,
//!   D187): [`gmsh`] (Gmsh `.msh` 4.1) and [`format::su2`] (SU2 `.su2`),
//!   [`format::fluent_msh`] (the face-based `.msh`), [`format::vtk`] (VTK
//!   `.vtk` and `.vtu`, also written), [`format::openfoam`] (polyMesh directories) and, with the `hdf5` feature,
//!   `format::cgns` (CGNS HDF5);
//! - [`periodic`]: translational periodic boundary pairs (D29, D109);
//! - [`renumber`]: reverse Cuthill–McKee renumbering of the cells;
//! - [`locate`]: point location and cell-to-cell particle tracking (D47);
//! - [`partition`]: domain decomposition (D32);
//! - [`decompose`]: local meshes with halo cells and the halo exchange plan;
//! - [`vtkhdf`]: VTKHDF output (D41), with the feature `hdf5`.

pub mod annulus;
pub mod builder;
pub mod case;
pub mod case_format;
pub mod cylinder;
pub mod decompose;
pub mod format;
pub mod generate;
pub mod geometry;
pub mod handover;
pub mod imprint;
pub mod isosurface;
pub mod locate;
pub mod mesh;
pub mod partition;
pub mod periodic;
pub mod quality;
pub mod query;
pub mod renumber;
pub mod revolved;
pub mod surface;
pub mod wedge;

pub use format::gmsh;
pub use geometry::Geometry;
pub use mesh::{Boundary, CellShape, Mesh, Periodic, PeriodicPair};
#[cfg(feature = "hdf5")]
pub mod vtkhdf;
