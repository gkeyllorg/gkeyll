// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_euler_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <assert.h>

typedef int (*pkpm_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double *p_ij, const double *pkpm_div_ppar);

typedef void (*pkpm_copy_t)(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim, double* GKYL_RESTRICT prim_surf);

typedef int (*pkpm_u_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm);

typedef void (*pkpm_u_copy_t)(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT pkpm_u);

typedef void (*pkpm_pressure_t)(const double *bvar, const double *vlasov_pkpm_moms, 
  double* GKYL_RESTRICT p_ij);

typedef void (*pkpm_p_force_t)(const double *prim_c, const double *div_b, 
  double* GKYL_RESTRICT pkpm_accel);

typedef void (*pkpm_int_t)(const double *vlasov_pkpm_moms, 
  const double *euler_pkpm, const double* prim, 
  double* GKYL_RESTRICT int_pkpm_vars); 

typedef void (*pkpm_source_t)(const double* qmem, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  double* GKYL_RESTRICT out);

typedef void (*pkpm_io_t)(const double *vlasov_pkpm_moms, 
  const double *euler_pkpm, const double* p_ij, 
  const double* prim, const double* pkpm_accel, 
  double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io); 

typedef void (*pkpm_accel_t)(const double *dxv, 
  const double *prim_surf_l, const double *prim_surf_c, const double *prim_surf_r, 
  const double *prim_c, const double *bvar_c, const double *nu_c, 
  double* GKYL_RESTRICT pkpm_accel); 

typedef void (*pkpm_penalization_t)(double tol, bool force_lax, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_r,
  const double *p_ij_l, const double *p_ij_r,
  const double *prim_l, const double *prim_r, 
  const double *euler_pkpm_l, const double *euler_pkpm_r,
  double* GKYL_RESTRICT pkpm_lax, double* GKYL_RESTRICT pkpm_penalization); 

typedef void (*pkpm_limiter_t)(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom, const double *prim_c, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *p_ij_l, const double *p_ij_c, const double *p_ij_r, 
  double *euler_pkpm_l, double *euler_pkpm_c, double *euler_pkpm_r);

// for use in kernel tables
typedef struct { pkpm_set_t kernels[3]; } gkyl_dg_pkpm_set_kern_list;
typedef struct { pkpm_copy_t kernels[3]; } gkyl_dg_pkpm_copy_kern_list;
typedef struct { pkpm_u_set_t kernels[3]; } gkyl_dg_pkpm_u_set_kern_list;
typedef struct { pkpm_u_copy_t kernels[3]; } gkyl_dg_pkpm_u_copy_kern_list;
typedef struct { pkpm_pressure_t kernels[3]; } gkyl_dg_pkpm_pressure_kern_list;
typedef struct { pkpm_p_force_t kernels[3]; } gkyl_dg_pkpm_p_force_kern_list;

typedef struct { pkpm_int_t kernels[3]; } gkyl_dg_pkpm_int_kern_list;
typedef struct { pkpm_source_t kernels[3]; } gkyl_dg_pkpm_source_kern_list;
typedef struct { pkpm_io_t kernels[3]; } gkyl_dg_pkpm_io_kern_list;

typedef struct { pkpm_accel_t kernels[3]; } gkyl_dg_pkpm_accel_kern_list;
typedef struct { pkpm_penalization_t kernels[3]; } gkyl_dg_pkpm_penalization_kern_list;
typedef struct { pkpm_limiter_t kernels[4]; } gkyl_dg_pkpm_limiter_kern_list;

struct gkyl_dg_calc_pkpm_vars {
  struct gkyl_rect_grid conf_grid; // Configuration space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int poly_order; // polynomial order (determines whether we solve linear system or use basis_inv method)
  struct gkyl_range mem_range; // Configuration space range for linear solve

  const struct gkyl_wv_eqn *wv_eqn; // Wave equation for characteristic limiting of solution
  const struct gkyl_wave_geom *geom; // Wave geometry for rotating solution
  double limiter_fac; // Factor for relationship between cell slopes and cell average differences (default: 1/sqrt(3))

  double tol; // Tolerance in mass density and average normal velocity at the interface
              // for switching to Lax fluxes in computing penalization of the momentum solve (default: 1.0e-12)
  bool force_lax; // Boolean for forcing the system to use Lax fluxes for the momentum (default: false)

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int Ncomp; // number of components in the linear solve (6 variables being solved for)

  struct gkyl_nmat *As_u, *xs_u; // matrices for LHS and RHS for flow velocity solve
  gkyl_nmat_mem *mem_u; // memory for use in batched linear solve for velocity

  pkpm_set_t pkpm_set;  // kernel for setting matrices for linear solve
  pkpm_copy_t pkpm_copy; // kernel for copying solution to output; also computed needed surface expansions
  pkpm_u_set_t pkpm_u_set;  // kernel for setting matrices for linear solve for flow velocity
  pkpm_u_copy_t pkpm_u_copy; // kernel for copying solution for flow velocity to output
  pkpm_pressure_t pkpm_pressure; // kernel for computing pressure (Volume and surface expansion)
  pkpm_p_force_t pkpm_p_force; // kernel for computing pressure force p_force = 1/rho div(p_par b) - T_perp/m div(b)

  pkpm_int_t pkpm_int; // kernel for computing integrated pkpm variables
  pkpm_source_t pkpm_source; // kernel for computing pkpm source update
  pkpm_io_t pkpm_io; // kernel for constructing I/O arrays for pkpm diagnostics

  pkpm_accel_t pkpm_accel[3]; // kernel for computing pkpm acceleration and Lax variables
  pkpm_penalization_t pkpm_penalization[3]; // kernel for computing pkpm acceleration and Lax variables
  pkpm_limiter_t pkpm_limiter[3]; // kernel for limiting slopes of fluid variables

  uint32_t flags;
  struct gkyl_dg_calc_pkpm_vars *on_dev; // pointer to itself or device data
};

// Set matrices for computing pkpm primitive vars, e.g., ux,uy,uz (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_set_kern_list ser_pkpm_set_kernels[] = {
  { NULL, pkpm_vars_set_1x_ser_p1, pkpm_vars_set_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_set_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_set_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing pkpm primitive vars, e.g., ux,uy,uz (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_set_kern_list ten_pkpm_set_kernels[] = {
  { NULL, pkpm_vars_set_1x_ser_p1, pkpm_vars_set_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_set_2x_ser_p1, pkpm_vars_set_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_set_3x_ser_p1, NULL }, // 2
};

// Copy solution for pkpm primitive vars, e.g., ux,uy,uz (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_copy_kern_list ser_pkpm_copy_kernels[] = {
  { NULL, pkpm_vars_copy_1x_ser_p1, pkpm_vars_copy_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_copy_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_copy_3x_ser_p1, NULL }, // 2
};

// Copy solution for pkpm primitive vars, e.g., ux,uy,uz (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_copy_kern_list ten_pkpm_copy_kernels[] = {
  { NULL, pkpm_vars_copy_1x_ser_p1, pkpm_vars_copy_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_copy_2x_ser_p1, pkpm_vars_copy_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_copy_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_u_set_kern_list ser_pkpm_u_set_kernels[] = {
  { NULL, pkpm_vars_u_set_1x_ser_p1, pkpm_vars_u_set_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_u_set_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_u_set_3x_ser_p1, NULL }, // 2
};

// Set matrices for computing flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_u_set_kern_list ten_pkpm_u_set_kernels[] = {
  { NULL, pkpm_vars_u_set_1x_ser_p1, pkpm_vars_u_set_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_u_set_2x_ser_p1, pkpm_vars_u_set_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_u_set_3x_ser_p1, NULL }, // 2
};

// Copy solution for flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_u_copy_kern_list ser_pkpm_u_copy_kernels[] = {
  { NULL, pkpm_vars_u_copy_1x_ser_p1, pkpm_vars_u_copy_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_u_copy_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_u_copy_3x_ser_p1, NULL }, // 2
};

// Copy solution for flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_u_copy_kern_list ten_pkpm_u_copy_kernels[] = {
  { NULL, pkpm_vars_u_copy_1x_ser_p1, pkpm_vars_u_copy_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_u_copy_2x_ser_p1, pkpm_vars_u_copy_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_u_copy_3x_ser_p1, NULL }, // 2
};

// PKPM Pressure (p_ij = (p_par - p_perp)b_i b_j + p_perp g_ij) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ser_pkpm_pressure_kernels[] = {
  { NULL, pkpm_vars_pressure_1x_ser_p1, pkpm_vars_pressure_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_pressure_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_pressure_3x_ser_p1, NULL }, // 2
};

// PKPM Pressure (p_ij = (p_ij = (p_par - p_perp)b_i b_j + p_perp g_ij) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_pressure_kern_list ten_pkpm_pressure_kernels[] = {
  { NULL, pkpm_vars_pressure_1x_ser_p1, pkpm_vars_pressure_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_pressure_2x_ser_p1, pkpm_vars_pressure_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_pressure_3x_ser_p1, NULL }, // 2
};

// PKPM Pressure force p_force = 1/rho div(p_par b) - T_perp/m div(b) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_p_force_kern_list ser_pkpm_p_force_kernels[] = {
  { NULL, pkpm_vars_p_force_1x_ser_p1, pkpm_vars_p_force_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_p_force_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_p_force_3x_ser_p1, NULL }, // 2
};

// PKPM Pressure force p_force = 1/rho div(p_par b) - T_perp/m div(b) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_p_force_kern_list ten_pkpm_p_force_kernels[] = {
  { NULL, pkpm_vars_p_force_1x_ser_p1, pkpm_vars_p_force_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_p_force_2x_ser_p1, pkpm_vars_p_force_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_p_force_3x_ser_p1, NULL }, // 2
};

// PKPM integrated variables integral (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_int_kern_list ser_pkpm_int_kernels[] = {
  { NULL, pkpm_vars_integrated_1x_ser_p1, pkpm_vars_integrated_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_integrated_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_integrated_3x_ser_p1, NULL }, // 2
};

// PKPM integrated variables integral (rho, rhoux, rhouy, rhouz, rhoux^2, rhouy^2, rhouz^2, p_parallel, p_perp) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_int_kern_list ten_pkpm_int_kernels[] = {
  { NULL, pkpm_vars_integrated_1x_ser_p1, pkpm_vars_integrated_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_integrated_2x_ser_p1, pkpm_vars_integrated_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_integrated_3x_ser_p1, NULL }, // 2
};

// PKPM explicit source solve (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_source_kern_list ser_pkpm_source_kernels[] = {
  { NULL, euler_pkpm_source_1x_ser_p1, euler_pkpm_source_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_source_2x_ser_p1, NULL }, // 1
  { NULL, euler_pkpm_source_3x_ser_p1, NULL }, // 2
};

// PKPM explicit source solve (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_source_kern_list ten_pkpm_source_kernels[] = {
  { NULL, euler_pkpm_source_1x_ser_p1, euler_pkpm_source_1x_ser_p2 }, // 0
  { NULL, euler_pkpm_source_2x_ser_p1, euler_pkpm_source_2x_tensor_p2 }, // 1
  { NULL, euler_pkpm_source_3x_ser_p1, NULL }, // 2
};

// PKPM io variables (Serendipity kernels)
// Conserved fluid variables: [rho, rho ux, rho uy, rho uz, Pxx + rho ux^2, Pxy + rho ux uy, Pxz + rho ux uz, Pyy + rho uy^2, Pyz + rho uy uz, Pzz + rho uz^2]
// PKPM primitive and acceleration variables:  
// [ux, uy, uz, T_perp/m, m/T_perp, 1/rho div(p_par b), T_perp/m div(b), bb : grad(u)]
GKYL_CU_D
static const gkyl_dg_pkpm_io_kern_list ser_pkpm_io_kernels[] = {
  { NULL, pkpm_vars_io_1x_ser_p1, pkpm_vars_io_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_io_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_io_3x_ser_p1, NULL }, // 2
};

// PKPM io variables (Tensor kernels)
// Conserved fluid variables: [rho, rho ux, rho uy, rho uz, Pxx + rho ux^2, Pxy + rho ux uy, Pxz + rho ux uz, Pyy + rho uy^2, Pyz + rho uy uz, Pzz + rho uz^2]
// PKPM primitive and acceleration variables:  
// [ux, uy, uz, T_perp/m, m/T_perp, 1/rho div(p_par b), T_perp/m div(b), bb : grad(u)]
GKYL_CU_D
static const gkyl_dg_pkpm_io_kern_list ten_pkpm_io_kernels[] = {
  { NULL, pkpm_vars_io_1x_ser_p1, pkpm_vars_io_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_io_2x_ser_p1, pkpm_vars_io_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_io_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (in x) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ser_pkpm_accel_x_kernels[] = {
  { NULL, pkpm_vars_accel_x_1x_ser_p1, pkpm_vars_accel_x_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_accel_x_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_accel_x_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (in y) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ser_pkpm_accel_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_vars_accel_y_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_accel_y_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (in z) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ser_pkpm_accel_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_vars_accel_z_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ten_pkpm_accel_x_kernels[] = {
  { NULL, pkpm_vars_accel_x_1x_ser_p1, pkpm_vars_accel_x_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_accel_x_2x_ser_p1, pkpm_vars_accel_x_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_accel_x_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (in y) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ten_pkpm_accel_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_vars_accel_y_2x_ser_p1, pkpm_vars_accel_y_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_accel_y_3x_ser_p1, NULL }, // 2
};

// PKPM acceleration variables, e.g., bb:grad(u), (in z) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_accel_kern_list ten_pkpm_accel_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_vars_accel_z_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in x) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ser_pkpm_penalization_x_kernels[] = {
  { NULL, pkpm_vars_penalization_x_1x_ser_p1, pkpm_vars_penalization_x_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_penalization_x_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_penalization_x_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and 
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in y) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ser_pkpm_penalization_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_vars_penalization_y_2x_ser_p1, NULL }, // 1
  { NULL, pkpm_vars_penalization_y_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and 
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in z) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ser_pkpm_penalization_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_vars_penalization_z_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and 
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in x) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ten_pkpm_penalization_x_kernels[] = {
  { NULL, pkpm_vars_penalization_x_1x_ser_p1, pkpm_vars_penalization_x_1x_ser_p2 }, // 0
  { NULL, pkpm_vars_penalization_x_2x_ser_p1, pkpm_vars_penalization_x_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_penalization_x_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and 
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in y) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ten_pkpm_penalization_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, pkpm_vars_penalization_y_2x_ser_p1, pkpm_vars_penalization_y_2x_tensor_p2 }, // 1
  { NULL, pkpm_vars_penalization_y_3x_ser_p1, NULL }, // 2
};

// PKPM penalization variables, e.g., total momentum penalization and 
// Lax penalization (lambda_i = |u_i| + sqrt(3*T_ii/m)) (in z) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_penalization_kern_list ten_pkpm_penalization_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, pkpm_vars_penalization_z_3x_ser_p1, NULL }, // 2
};

// Characteristic limiter in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ser_pkpm_limiter_x_kernels[] = {
  { NULL, euler_pkpm_limiter_x_1x_ser_p1, euler_pkpm_limiter_x_1x_ser_p2, euler_pkpm_limiter_x_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_limiter_x_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_pkpm_limiter_x_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ser_pkpm_limiter_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_limiter_y_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_pkpm_limiter_y_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ser_pkpm_limiter_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_limiter_z_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ten_pkpm_limiter_x_kernels[] = {
  { NULL, euler_pkpm_limiter_x_1x_ser_p1, euler_pkpm_limiter_x_1x_ser_p2, euler_pkpm_limiter_x_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_limiter_x_2x_ser_p1, euler_pkpm_limiter_x_2x_tensor_p2, NULL }, // 1
  { NULL, euler_pkpm_limiter_x_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ten_pkpm_limiter_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, euler_pkpm_limiter_y_2x_ser_p1, euler_pkpm_limiter_y_2x_tensor_p2, NULL }, // 1
  { NULL, euler_pkpm_limiter_y_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in z (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_limiter_kern_list ten_pkpm_limiter_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, euler_pkpm_limiter_z_3x_ser_p1, NULL, NULL }, // 2
};

GKYL_CU_D
static pkpm_set_t
choose_pkpm_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_copy_t
choose_pkpm_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_u_set_t
choose_pkpm_u_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_u_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_u_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_u_copy_t
choose_pkpm_u_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_u_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_u_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_pressure_t
choose_pkpm_pressure_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_p_force_t
choose_pkpm_p_force_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_p_force_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_p_force_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_int_t
choose_pkpm_int_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_int_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_int_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_source_t
choose_pkpm_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_source_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_source_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_io_t
choose_pkpm_io_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_io_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_io_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_accel_t
choose_pkpm_accel_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_pkpm_accel_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_pkpm_accel_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ser_pkpm_accel_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return ten_pkpm_accel_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ten_pkpm_accel_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ten_pkpm_accel_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_penalization_t
choose_pkpm_penalization_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_pkpm_penalization_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_pkpm_penalization_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ser_pkpm_penalization_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return ten_pkpm_penalization_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ten_pkpm_penalization_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ten_pkpm_penalization_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_limiter_t
choose_pkpm_limiter_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_pkpm_limiter_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_pkpm_limiter_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ser_pkpm_limiter_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return ten_pkpm_limiter_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ten_pkpm_limiter_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ten_pkpm_limiter_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}
