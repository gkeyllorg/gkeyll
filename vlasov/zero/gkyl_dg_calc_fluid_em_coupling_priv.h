// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*fluid_em_coupling_set_t)(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em);

typedef void (*fluid_em_coupling_copy_t)(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em);

typedef void (*fluid_em_coupling_energy_t)(const double* ke_old, const double* ke_new, 
  double* GKYL_RESTRICT fluid);

// for use in kernel tables
typedef struct { fluid_em_coupling_set_t kernels[4]; } gkyl_dg_fluid_em_coupling_set_kern_list;
typedef struct { fluid_em_coupling_copy_t kernels[4]; } gkyl_dg_fluid_em_coupling_copy_kern_list;
typedef struct { fluid_em_coupling_energy_t kernels[4]; } gkyl_dg_fluid_em_coupling_energy_kern_list;

struct gkyl_dg_calc_fluid_em_coupling {
  int cdim; // Configuration space dimensionality
  int poly_order; // polynomial order (determines whether we solve linear system or use basis_inv method)
  struct gkyl_range mem_range; // Configuration space range for linear solve

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve

  fluid_em_coupling_set_t fluid_em_coupling_set;  // kernel for setting matrices for linear solve
  fluid_em_coupling_copy_t fluid_em_coupling_copy; // kernel for copying solution to output
  fluid_em_coupling_energy_t fluid_em_coupling_energy; // kernel for computing energy from updated solution (for Euler/5-moment)

  int num_fluids; // number of fluids being implicitly solved for
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  double epsilon0; // permitivvity of free space

  uint32_t flags;
  struct gkyl_dg_calc_fluid_em_coupling *on_dev; // pointer to itself or device data
};

// Set matrices for computing implicit source solve for fluid-em coupling (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_set_kern_list ser_fluid_em_coupling_set_kernels[] = {
  { NULL, fluid_em_coupling_set_1x_ser_p1, fluid_em_coupling_set_1x_ser_p2, fluid_em_coupling_set_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_set_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_em_coupling_set_3x_ser_p1, NULL, NULL }, // 2
};

// Set matrices for computing implicit source solve for fluid-em coupling (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_set_kern_list ten_fluid_em_coupling_set_kernels[] = {
  { NULL, fluid_em_coupling_set_1x_ser_p1, fluid_em_coupling_set_1x_ser_p2, fluid_em_coupling_set_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_set_2x_ser_p1, fluid_em_coupling_set_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_em_coupling_set_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for implicit source solve for fluid-em coupling (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_copy_kern_list ser_fluid_em_coupling_copy_kernels[] = {
  { NULL, fluid_em_coupling_copy_1x_ser_p1, fluid_em_coupling_copy_1x_ser_p2, fluid_em_coupling_copy_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_copy_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_em_coupling_copy_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for implicit source solve for fluid-em coupling (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_copy_kern_list ten_fluid_em_coupling_copy_kernels[] = {
  { NULL, fluid_em_coupling_copy_1x_ser_p1, fluid_em_coupling_copy_1x_ser_p2, fluid_em_coupling_copy_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_copy_2x_ser_p1, fluid_em_coupling_copy_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_em_coupling_copy_3x_ser_p1, NULL, NULL }, // 2
};

// Compute energy from updated kinetic energy and old kinetic energy (Euler/5-moment) (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_energy_kern_list ser_fluid_em_coupling_energy_kernels[] = {
  { NULL, fluid_em_coupling_energy_1x_ser_p1, fluid_em_coupling_energy_1x_ser_p2, fluid_em_coupling_energy_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_energy_2x_ser_p1, NULL, NULL }, // 1
  { NULL, fluid_em_coupling_energy_3x_ser_p1, NULL, NULL }, // 2
};

// Compute energy from updated kinetic energy and old kinetic energy (Euler/5-moment) (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_fluid_em_coupling_energy_kern_list ten_fluid_em_coupling_energy_kernels[] = {
  { NULL, fluid_em_coupling_energy_1x_ser_p1, fluid_em_coupling_energy_1x_ser_p2, fluid_em_coupling_energy_1x_ser_p3 }, // 0
  { NULL, fluid_em_coupling_energy_2x_ser_p1, fluid_em_coupling_energy_2x_tensor_p2, NULL }, // 1
  { NULL, fluid_em_coupling_energy_3x_ser_p1, NULL, NULL }, // 2
};

GKYL_CU_D
static fluid_em_coupling_set_t
choose_fluid_em_coupling_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_em_coupling_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_em_coupling_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_em_coupling_copy_t
choose_fluid_em_coupling_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_em_coupling_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_em_coupling_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static fluid_em_coupling_energy_t
choose_fluid_em_coupling_energy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_fluid_em_coupling_energy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_fluid_em_coupling_energy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}
