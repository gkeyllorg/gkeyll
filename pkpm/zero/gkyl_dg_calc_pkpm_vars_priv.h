// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_vars_kernels.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*pkpm_grad_u_t)(
  const double *dxv, const double *jacob_pos, const double *u_l, const double *u_c, const double *u_r,
  double *GKYL_RESTRICT grad_u
);

typedef struct {
  pkpm_grad_u_t kernels[3];
} gkyl_dg_pkpm_grad_u_kern_list;

typedef void (*pkpm_diag_t)(
  const double *dxv, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *bvar,
  const double *grad_b, const double *u, const double *grad_u, const double *theta,
  const double *alpha_0, const double *xi, const double *c_h, double *GKYL_RESTRICT fluid_io,
  double *GKYL_RESTRICT pkpm_vars_io
);

typedef int (*pkpm_closure_t)(
  double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms, const double *bvar,
  const double *wn, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up, double *GKYL_RESTRICT c_0, double *GKYL_RESTRICT alpha_0,
  double *GKYL_RESTRICT c_h, double *GKYL_RESTRICT a_perp, double *GKYL_RESTRICT flux_lo,
  double *GKYL_RESTRICT flux_up, double *GKYL_RESTRICT cfl, double *GKYL_RESTRICT cond
);
typedef int (*pkpm_prim_t)(const double *fluid, double *GKYL_RESTRICT u, double *GKYL_RESTRICT xi);
typedef void (*pkpm_int_t)(
  const double *fluid, const double *vlasov_pkpm_moms, const double *u, double *GKYL_RESTRICT int_pkpm_vars
);

typedef int (*pkpm_fluid_positivity_t)(double eps, double *GKYL_RESTRICT fluid);
typedef int (*pkpm_grad_ln_n_t)(
  const double *rho, const double *grad_rho, double *GKYL_RESTRICT grad_ln_n
);

typedef int (*pkpm_normalize_t)(const double *n, const double *f, double *GKYL_RESTRICT fbar);
typedef void (*pkpm_transverse_t)(const double *bvar, double *GKYL_RESTRICT f);
typedef void (*pkpm_source_t)(const double *qmem, const double *fluid, double *GKYL_RESTRICT out);

struct gkyl_dg_calc_pkpm_vars {
  struct gkyl_rect_grid conf_grid; // Configuration space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int poly_order; // polynomial order (p=1 serendipity only)
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy)
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device)

  pkpm_grad_u_t pkpm_grad_u[3]; // weak (surface-averaged) gradient of a 3-vector in each direction
  pkpm_closure_t pkpm_closure; // closures c_0, a_P/xi and c from the moments of the stage-A rhs, with the velocity edges
  pkpm_prim_t pkpm_prim; // nodal flow and thermal speed [xi, xi^2, p] of the fluid state
  pkpm_int_t pkpm_int; // cell integrals of the diagnostics of the [rho, M, E] state
  pkpm_grad_ln_n_t pkpm_grad_ln_n; // grad(ln n) by nodal quotient of the weak gradient of rho
  pkpm_fluid_positivity_t pkpm_fluid_positivity; // vertex positivity limiter of the fluid state
  pkpm_normalize_t pkpm_normalize; // nodal constraint enforcement of the initial distributions
  pkpm_transverse_t pkpm_transverse; // transversality projection of the first-harmonic blocks
  pkpm_source_t pkpm_source; // explicit electromagnetic source of [rho, M, E]
  pkpm_diag_t pkpm_diag; // diagnostic arrays
  int num_null_rho; // number of nodes with rho <= 0 in the last nodal computation (host copy)
  int *num_null_rho_dev; // device counter for num_null_rho (0 on host)
  double closure_cond_min; // smallest det/det0 of the closure solves of the last call (host copy)
  unsigned long long *closure_cond_min_dev; // device minimum (bit pattern of the positive double)
  int *num_pos_dev; // device counters [rho, e, xi, avg] of the positivity limiter (0 on host)

  uint32_t flags;
  struct gkyl_dg_calc_pkpm_vars *on_dev; // pointer to itself or device data
};


GKYL_CU_D static const pkpm_diag_t ser_pkpm_diag_kernels[] = {
  pkpm_vars_diag_1x_ser_p1,
  pkpm_vars_diag_2x_ser_p1,
  pkpm_vars_diag_3x_ser_p1,
};

// Weak gradient of a 3-vector (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_grad_u_kern_list ser_pkpm_grad_u_kernels[] = {
  {pkpm_vars_grad_u_x_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_vars_grad_u_x_2x_ser_p1, pkpm_vars_grad_u_y_2x_ser_p1, NULL}, // 1
  {pkpm_vars_grad_u_x_3x_ser_p1, pkpm_vars_grad_u_y_3x_ser_p1, pkpm_vars_grad_u_z_3x_ser_p1}, // 2
};

GKYL_CU_D static pkpm_grad_u_t
choose_pkpm_grad_u_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_grad_u_kernels[cdim - 1].kernels[dir];
  }
  return 0;
}

GKYL_CU_D static pkpm_diag_t
choose_pkpm_diag_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_diag_kernels[cdim - 1];
  }
  return 0;
}

// Closures of the factored, normalized update: c_0, a_P/xi and the heating rate c (p=1 serendipity)
GKYL_CU_D static const pkpm_closure_t ser_pkpm_closure_kernels[] = {
  pkpm_vars_closure_1x_ser_p1,
  pkpm_vars_closure_2x_ser_p1,
  pkpm_vars_closure_3x_ser_p1,
};
GKYL_CU_D
static pkpm_closure_t
choose_pkpm_closure_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_closure_kernels[cdim - 1];
  }
  return 0;
}
// Nodal thermal speed [xi, xi^2] of the fluid state (p=1 serendipity)
GKYL_CU_D static const pkpm_prim_t ser_pkpm_prim_kernels[] = {
  pkpm_vars_prim_1x_ser_p1,
  pkpm_vars_prim_2x_ser_p1,
  pkpm_vars_prim_3x_ser_p1,
};
GKYL_CU_D
static pkpm_prim_t
choose_pkpm_prim_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_prim_kernels[cdim - 1];
  }
  return 0;
}
GKYL_CU_D static const pkpm_int_t ser_pkpm_int_kernels[] = {
  pkpm_vars_int_1x_ser_p1,
  pkpm_vars_int_2x_ser_p1,
  pkpm_vars_int_3x_ser_p1,
};
GKYL_CU_D
static pkpm_int_t
choose_pkpm_int_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_int_kernels[cdim - 1];
  }
  return 0;
}
// grad(ln n) from the weak gradient of rho over the cell average (p=1 serendipity)
GKYL_CU_D static const pkpm_fluid_positivity_t ser_pkpm_fluid_positivity_kernels[] = {
  pkpm_fluid_positivity_1x_ser_p1,
  pkpm_fluid_positivity_2x_ser_p1,
  pkpm_fluid_positivity_3x_ser_p1,
};
GKYL_CU_D
static pkpm_fluid_positivity_t
choose_pkpm_fluid_positivity_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_fluid_positivity_kernels[cdim - 1];
  }
  return 0;
}

GKYL_CU_D static const pkpm_grad_ln_n_t ser_pkpm_grad_ln_n_kernels[] = {
  pkpm_vars_grad_ln_n_1x_ser_p1,
  pkpm_vars_grad_ln_n_2x_ser_p1,
  pkpm_vars_grad_ln_n_3x_ser_p1,
};
GKYL_CU_D
static pkpm_grad_ln_n_t
choose_pkpm_grad_ln_n_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_grad_ln_n_kernels[cdim - 1];
  }
  return 0;
}
// Nodal constraint enforcement of the initial normalized distributions (p=1 serendipity)
GKYL_CU_D static const pkpm_normalize_t ser_pkpm_normalize_kernels[] = {
  pkpm_dist_normalize_1x1v_ser_p1,
  pkpm_dist_normalize_2x1v_ser_p1,
  pkpm_dist_normalize_3x1v_ser_p1,
};
GKYL_CU_D
static pkpm_normalize_t
choose_pkpm_normalize_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_normalize_kernels[cdim - 1];
  }
  return 0;
}
// Transversality projection of the first-harmonic blocks (p=1 serendipity)
GKYL_CU_D static const pkpm_transverse_t ser_pkpm_transverse_kernels[] = {
  pkpm_dist_transverse_1x1v_ser_p1,
  pkpm_dist_transverse_2x1v_ser_p1,
  pkpm_dist_transverse_3x1v_ser_p1,
};
GKYL_CU_D
static pkpm_transverse_t
choose_pkpm_transverse_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_transverse_kernels[cdim - 1];
  }
  return 0;
}
// Explicit electromagnetic source of [rho, M, E] (p=1 serendipity)
GKYL_CU_D static const pkpm_source_t ser_pkpm_source_kernels[] = {
  pkpm_vars_source_1x_ser_p1,
  pkpm_vars_source_2x_ser_p1,
  pkpm_vars_source_3x_ser_p1,
};
GKYL_CU_D
static pkpm_source_t
choose_pkpm_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    return ser_pkpm_source_kernels[cdim - 1];
  }
  return 0;
}
