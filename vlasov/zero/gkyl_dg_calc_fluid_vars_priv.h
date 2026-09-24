// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <assert.h>

typedef int (*fluid_set_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *fluid);

typedef void (*fluid_copy_t)(int count, struct gkyl_nmat *x, 
  double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf);

typedef void (*fluid_pressure_t)(double gas_gamma, const double *fluid, const double *u, 
  double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf);

typedef void (*fluid_ke_t)(const double *fluid, const double *u, 
  double* GKYL_RESTRICT ke);

typedef void (*fluid_limiter_t)(double limiter_fac, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
  double *fluid_l, double *fluid_c, double *fluid_r);

typedef void (*fluid_int_t)(const double *fluid, 
  const double* u_i, const double* p_ij, 
  double* GKYL_RESTRICT int_fluid_vars); 

typedef void (*fluid_source_t)(const double* app_accel, const double* fluid, 
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { fluid_set_t kernels[4]; } gkyl_dg_fluid_set_kern_list;
typedef struct { fluid_copy_t kernels[4]; } gkyl_dg_fluid_copy_kern_list;
typedef struct { fluid_pressure_t kernels[4]; } gkyl_dg_fluid_pressure_kern_list;
typedef struct { fluid_ke_t kernels[4]; } gkyl_dg_fluid_ke_kern_list;
typedef struct { fluid_limiter_t kernels[4]; } gkyl_dg_fluid_limiter_kern_list;
typedef struct { fluid_int_t kernels[4]; } gkyl_dg_fluid_int_kern_list;
typedef struct { fluid_source_t kernels[4]; } gkyl_dg_fluid_source_kern_list;

struct gkyl_dg_calc_fluid_vars {
  enum gkyl_eqn_type eqn_type; // Equation type
  const struct gkyl_wv_eqn *wv_eqn; // Wave equation for characteristic limiting of solution
  const struct gkyl_wave_geom *geom; // Wave geometry for rotating solution
  double param; // parameter for computing primitive moments/limiting solution (vt for isothermal Euler, gas_gammas for Euler)

  int cdim; // Configuration space dimensionality
  int poly_order; // polynomial order (determines whether we solve linear system or use basis_inv method)
  struct gkyl_range mem_range; // Configuration space range for linear solve

  double limiter_fac; // Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3))

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int Ncomp; // number of components in the linear solve (6 variables being solved for)

  fluid_set_t fluid_set;  // kernel for setting matrices for linear solve
  fluid_copy_t fluid_copy; // kernel for copying solution to output; also computed needed surface expansions
  fluid_pressure_t fluid_pressure; // kernel for computing pressure (Volume and surface expansion)
  fluid_ke_t fluid_ke; // kernel for computing kinetic energy (Volume expansion)
  fluid_limiter_t fluid_limiter[3]; // kernel for limiting slopes of fluid variables
  fluid_int_t fluid_int; // kernel for computing integrated fluid variables
  fluid_source_t fluid_source; // kernel for computing fluid source update

  uint32_t flags;
  struct gkyl_dg_calc_fluid_vars *on_dev; // pointer to itself or device data
};

// Set matrices for computing fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_set_kern_list ser_fluid_set_kernels[] = {
  { NULL, fluid_vars_u_set_1x_ser_p1, fluid_vars_u_set_1x_ser_p2, fluid_vars_u_set_1x_ser_p3 }, // 0
  { NULL, fluid_vars_u_set_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_u_set_3x_ser_p1, NULL, NULL }, // 2
};

// Set matrices for computing fluid flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_set_kern_list ten_fluid_set_kernels[] = {
  { NULL, fluid_vars_u_set_1x_ser_p1, fluid_vars_u_set_1x_ser_p2, fluid_vars_u_set_1x_ser_p3 }, // 0
  { NULL, fluid_vars_u_set_2x_ser_p1, fluid_vars_u_set_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_u_set_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_copy_kern_list ser_fluid_copy_kernels[] = {
  { NULL, fluid_vars_u_copy_1x_ser_p1, fluid_vars_u_copy_1x_ser_p2, fluid_vars_u_copy_1x_ser_p3 }, // 0
  { NULL, fluid_vars_u_copy_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_u_copy_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for fluid flow velocity (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_copy_kern_list ten_fluid_copy_kernels[] = {
  { NULL, fluid_vars_u_copy_1x_ser_p1, fluid_vars_u_copy_1x_ser_p2, fluid_vars_u_copy_1x_ser_p3 }, // 0
  { NULL, fluid_vars_u_copy_2x_ser_p1, fluid_vars_u_copy_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_u_copy_3x_ser_p1, NULL, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_pressure_kern_list ser_fluid_pressure_kernels[] = {
  { NULL, fluid_vars_pressure_1x_ser_p1, fluid_vars_pressure_1x_ser_p2, fluid_vars_pressure_1x_ser_p3 }, // 0
  { NULL, fluid_vars_pressure_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_pressure_3x_ser_p1, NULL, NULL }, // 2
};

// Scalar pressure Isothermal Euler -> p = vth*rho; Euler -> p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_pressure_kern_list ten_fluid_pressure_kernels[] = {
  { NULL, fluid_vars_pressure_1x_ser_p1, fluid_vars_pressure_1x_ser_p2, fluid_vars_pressure_1x_ser_p3 }, // 0
  { NULL, fluid_vars_pressure_2x_ser_p1, fluid_vars_pressure_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_pressure_3x_ser_p1, NULL, NULL }, // 2
};

// Kinetic energy = 1/2 (rho ux^2 + rho uy^2 + rho uz^2) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_ke_kern_list ser_fluid_ke_kernels[] = {
  { NULL, fluid_vars_ke_1x_ser_p1, fluid_vars_ke_1x_ser_p2, fluid_vars_ke_1x_ser_p3 }, // 0
  { NULL, fluid_vars_ke_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_ke_3x_ser_p1, NULL, NULL }, // 2
};

// Kinetic energy = 1/2 (rho ux^2 + rho uy^2 + rho uz^2) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_ke_kern_list ten_fluid_ke_kernels[] = {
  { NULL, fluid_vars_ke_1x_ser_p1, fluid_vars_ke_1x_ser_p2, fluid_vars_ke_1x_ser_p3 }, // 0
  { NULL, fluid_vars_ke_2x_ser_p1, fluid_vars_ke_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_ke_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ser_fluid_limiter_x_kernels[] = {
  { NULL, fluid_vars_limiterx_1x_ser_p1, fluid_vars_limiterx_1x_ser_p2, fluid_vars_limiterx_1x_ser_p3 }, // 0
  { NULL, fluid_vars_limiterx_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_limiterx_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ser_fluid_limiter_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, fluid_vars_limitery_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_limitery_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ser_fluid_limiter_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, fluid_vars_limiterz_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ten_fluid_limiter_x_kernels[] = {
  { NULL, fluid_vars_limiterx_1x_ser_p1, fluid_vars_limiterx_1x_ser_p2, fluid_vars_limiterx_1x_ser_p3 }, // 0
  { NULL, fluid_vars_limiterx_2x_ser_p1, fluid_vars_limiterx_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_limiterx_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ten_fluid_limiter_y_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, fluid_vars_limitery_2x_ser_p1, fluid_vars_limitery_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_limitery_3x_ser_p1, NULL, NULL }, // 2
};

// Characteristic limiter in z (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_limiter_kern_list ten_fluid_limiter_z_kernels[] = {
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, fluid_vars_limiterz_3x_ser_p1, NULL, NULL }, // 2
};

// Fluid integrated variables integral (rho, rhoux, rhouy, rhouz, rhou^2, p) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_int_kern_list ser_fluid_int_kernels[] = {
  { NULL, fluid_vars_integrated_1x_ser_p1, fluid_vars_integrated_1x_ser_p2, fluid_vars_integrated_1x_ser_p3 }, // 0
  { NULL, fluid_vars_integrated_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_integrated_3x_ser_p1, NULL, NULL }, // 2
};

// Fluid integrated variables integral (rho, rhoux, rhouy, rhouz, rhou^2, p) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_int_kern_list ten_fluid_int_kernels[] = {
  { NULL, fluid_vars_integrated_1x_ser_p1, fluid_vars_integrated_1x_ser_p2, fluid_vars_integrated_1x_ser_p3 }, // 0
  { NULL, fluid_vars_integrated_2x_ser_p1, fluid_vars_integrated_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_integrated_3x_ser_p1, NULL, NULL }, // 2
};

// Fluid explicit source solve (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_source_kern_list ser_fluid_source_kernels[] = {
  { NULL, fluid_vars_source_1x_ser_p1, fluid_vars_source_1x_ser_p2, fluid_vars_source_1x_ser_p3 }, // 0
  { NULL, fluid_vars_source_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_vars_source_3x_ser_p1, NULL, NULL }, // 2
};

// Fluid explicit source solve (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_source_kern_list ten_fluid_source_kernels[] = {
  { NULL, fluid_vars_source_1x_ser_p1, fluid_vars_source_1x_ser_p2, fluid_vars_source_1x_ser_p3 }, // 0
  { NULL, fluid_vars_source_2x_ser_p1, fluid_vars_source_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_vars_source_3x_ser_p1, NULL, NULL }, // 2
};

GKYL_CU_D
static fluid_set_t
choose_fluid_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_copy_t
choose_fluid_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_pressure_t
choose_fluid_pressure_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_pressure_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_ke_t
choose_fluid_ke_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_ke_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_ke_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_limiter_t
choose_fluid_limiter_kern(int dir, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_fluid_limiter_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_fluid_limiter_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ser_fluid_limiter_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return ten_fluid_limiter_x_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ten_fluid_limiter_y_kernels[cdim-1].kernels[poly_order];
      else if (dir == 2)
        return ten_fluid_limiter_z_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_int_t
choose_fluid_int_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_int_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_int_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_source_t
choose_fluid_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_source_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_source_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}
