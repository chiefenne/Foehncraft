//! Stage 5: the operations of the flow loop as PTX entry points `flow_*`,
//! one GPU thread per index. Each entry takes the number of indices and the
//! argument block of `foehncraft_kernels::entry` by value, and runs the
//! body that the CPU loops of `foehncraft-flow` (`device::kernel`) also
//! run. The host side is `crates/foehncraft-cuda/src/flow.rs`; the argument
//! blocks are one definition for both, so the argument lists cannot differ.

use crate::tid;
use foehncraft_kernels::entry::{self as e, Entry};

macro_rules! entries {
    ($($name:ident => $op:ty),* $(,)?) => {$(
        #[unsafe(no_mangle)]
        pub unsafe extern "ptx-kernel" fn $name(n: u64, a: <$op as Entry>::Args) {
            let i = tid();
            if (i as u64) < n {
                // SAFETY: i < n; the host built every array of `a` for the
                // mesh and keeps outputs distinct from inputs (see the
                // module docs of `foehncraft_kernels::entry`).
                unsafe { <$op as Entry>::run(i, &a) }
            }
        }
    )*};
}

entries! {
    flow_lincomb => e::Lincomb,
    flow_scaled_difference => e::ScaledDifference,
    flow_squares_add => e::SquaresAdd,
    flow_scale_vector => e::ScaleVector,
    flow_schur_diagonal => e::SchurDiagonal,
    flow_product => e::Product,
    flow_vertex_bounds => e::VertexBounds,
    flow_mlp_limiter => e::MlpLimiter,
    flow_pseudo_coefficient => e::PseudoCoefficient,
    flow_fill => e::Fill,
    flow_subtract => e::Subtract,
    flow_scale_entry => e::ScaleEntry,
    flow_block_sum_sq => e::BlockSumSq,
    flow_block_dot => e::BlockDot,
    flow_block_min_max => e::BlockMinMax,
    flow_boundary_constants => e::BoundaryConstants,
    flow_boundary_values => e::BoundaryValues,
    flow_gradient => e::Gradient,
    flow_limiter => e::Limiter,
    flow_mass_flux => e::MassFlux,
    flow_momentum_faces => e::MomentumFaces,
    flow_momentum_rows => e::MomentumRows,
    flow_momentum_operator => e::MomentumOperator,
    flow_residual_scale => e::ResidualScale,
    flow_set_diagonal => e::SetDiagonal,
    flow_diagonal_inverse => e::DiagonalInverse,
    flow_net_outflow => e::NetOutflow,
    flow_throughflow => e::Throughflow,
    flow_courant => e::Courant,
    flow_pressure_faces => e::PressureFaces,
    flow_close_backflow => e::CloseBackflow,
    flow_pressure_rows => e::PressureRows,
    flow_correct_fluxes => e::CorrectFluxes,
    flow_correct_velocity => e::CorrectVelocity,
    flow_neighbour_correction => e::NeighbourCorrection,
    flow_delta_flux => e::DeltaFlux,
    flow_transpose => e::Transpose,
    flow_matrix_product => e::MatrixProduct,
    flow_transposed_product => e::TransposedProduct,
    flow_scalar_faces => e::ScalarFaces,
    flow_scalar_rows => e::ScalarRows,
    flow_scalar_scale => e::ScalarScale,
    flow_sst_cells => e::SstCells,
    flow_strain_cells => e::StrainCells,
    flow_convective_derivative => e::ConvectiveDerivative,
    flow_curvature_cells => e::CurvatureCells,
    flow_qcr_cells => e::QcrCells,
    flow_explicit_stress_faces => e::ExplicitStressFaces,
    flow_quadratic_gradient => e::QuadraticGradient,
    flow_keps_cells => e::KepsCells,
    flow_zeta_f_cells => e::ZetaFCells,
    flow_turb_faces => e::TurbFaces,
    flow_scalar_turb_faces => e::ScalarTurbFaces,
    flow_scalar_wall_faces => e::ScalarWallFaces,
    flow_wall_faces => e::WallFaces,
    flow_wall_cells => e::WallCells,
    flow_floor => e::Floor,
    flow_quotient => e::Quotient,
    flow_vof_faces => e::VofFaces,
    flow_vof_update => e::VofUpdate,
    flow_vof_outflow_rate => e::VofOutflowRate,
    flow_vof_smooth => e::VofSmooth,
    flow_vof_gradient => e::VofGradient,
    flow_vof_curvature => e::VofCurvature,
    flow_vof_tension => e::VofTension,
    flow_vof_reconstruct => e::VofReconstruct,
    flow_two_fluid_cells => e::TwoFluidCells,
    flow_two_fluid_faces => e::TwoFluidFaces,
    flow_two_fluid_eliminate => e::TwoFluidEliminate,
}
