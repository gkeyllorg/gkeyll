#pragma once

// Private header for the fem_poisson_perp_adiabatic updater, not for direct use in user code.

#include <gkyl_fem_poisson_perp.h>
#include <gkyl_fem_poisson_perp_adiabatic.h>
#include <gkyl_mat.h>

struct gkyl_fem_poisson_perp_adiabatic {
  struct gkyl_basis basis; // Configuration-space basis.
  struct gkyl_range local; // Local range in which to solve.
  struct gkyl_range local_x; // 1D range in x of the flux-surface averages.
  struct gkyl_flux_surf_avg *fsa; // Flux-surface average and extension (not owned).
  struct gkyl_array *kJ; // Adiabatic response coefficient K.
  struct gkyl_fem_poisson_perp *helmholtz; // Solver for H = L + K.
  // In-place parallel projection of a right-hand side before the perpendicular solve, and
  // its context.
  void (*fem_parproj_op)(struct gkyl_array *rhs, void *ctx);
  void *fem_parproj_ctx;

  struct gkyl_array *phi_buff, *rhs_buff; // Buffers for the zonal response solve.
  struct gkyl_array *phi_rhs0; // Solution with rhs = 0 (Dirichlet values and biased lines).
  struct gkyl_array *psi; // Auxiliary 1D (x) potential: <phi1>, then the zonal solution psi.

  // Zonal system (I - G) psi = <phi1> of size m = Nx*num_basis_x, solved by applying the
  // inverse computed at init (with LAPACK) as a matrix-vector product (BLAS/cuBLAS).
  struct gkyl_array *zonal_rhs, *zonal_sol; // <phi1> and psi as single-cell m-vectors.
  gkyl_mat_mm_array_mem *zonal_inv; // Holds (I - G)^{-1} and applies it.
  double zonal_cond; // 1-norm condition number of I - G.
};
