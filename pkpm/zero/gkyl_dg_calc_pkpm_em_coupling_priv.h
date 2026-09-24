// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_euler_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*pkpm_em_coupling_set_t)(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool pkpm_field_static, double dt, 
  struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double *pkpm_u[GKYL_MAX_SPECIES],
  double* GKYL_RESTRICT em);

typedef void (*pkpm_em_coupling_copy_t)(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double *pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em);

// for use in kernel tables
typedef struct { pkpm_em_coupling_set_t kernels[4]; } gkyl_dg_pkpm_em_coupling_set_kern_list;
typedef struct { pkpm_em_coupling_copy_t kernels[4]; } gkyl_dg_pkpm_em_coupling_copy_kern_list;

struct gkyl_dg_calc_pkpm_em_coupling {
  struct gkyl_range mem_range; // Configuration space range for linear solve

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve

  pkpm_em_coupling_set_t pkpm_em_coupling_set;  // kernel for setting matrices for linear solve
  pkpm_em_coupling_copy_t pkpm_em_coupling_copy; // kernel for copying solution to output

  bool pkpm_field_static; // bool to determine if we are updating the self-consistent EM fields (dE/dt)

  int num_species; // number of species being implicitly solved for
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  double epsilon0; // permitivvity of free space

  uint32_t flags;
  struct gkyl_dg_calc_pkpm_em_coupling *on_dev; // pointer to itself or device data
};

// Set matrices for computing implicit source solve for fluid-em coupling in the PKPM system. (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_em_coupling_set_kern_list ser_pkpm_em_coupling_set_kernels[] = {
  { NULL, euler_pkpm_em_coupling_set_1x_ser_p1, euler_pkpm_em_coupling_set_1x_ser_p2, euler_pkpm_em_coupling_set_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_em_coupling_set_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_pkpm_em_coupling_set_3x_ser_p1, NULL, NULL }, // 2
};

// Set matrices for computing implicit source solve for fluid-em coupling in the PKPM system. (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_em_coupling_set_kern_list ten_pkpm_em_coupling_set_kernels[] = {
  { NULL, euler_pkpm_em_coupling_set_1x_ser_p1, euler_pkpm_em_coupling_set_1x_ser_p2, euler_pkpm_em_coupling_set_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_em_coupling_set_2x_ser_p1, euler_pkpm_em_coupling_set_2x_tensor_p2, NULL }, // 1
  { NULL, euler_pkpm_em_coupling_set_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for implicit source solve for fluid-em coupling in the PKPM system. (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_em_coupling_copy_kern_list ser_pkpm_em_coupling_copy_kernels[] = {
  { NULL, euler_pkpm_em_coupling_copy_1x_ser_p1, euler_pkpm_em_coupling_copy_1x_ser_p2, euler_pkpm_em_coupling_copy_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_em_coupling_copy_2x_ser_p1, NULL, NULL }, // 1
  { NULL, euler_pkpm_em_coupling_copy_3x_ser_p1, NULL, NULL }, // 2
};

// Copy solution for implicit source solve for fluid-em coupling in the PKPM system. (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_pkpm_em_coupling_copy_kern_list ten_pkpm_em_coupling_copy_kernels[] = {
  { NULL, euler_pkpm_em_coupling_copy_1x_ser_p1, euler_pkpm_em_coupling_copy_1x_ser_p2, euler_pkpm_em_coupling_copy_1x_ser_p3 }, // 0
  { NULL, euler_pkpm_em_coupling_copy_2x_ser_p1, euler_pkpm_em_coupling_copy_2x_tensor_p2, NULL }, // 1
  { NULL, euler_pkpm_em_coupling_copy_3x_ser_p1, NULL, NULL }, // 2
};

GKYL_CU_D
static pkpm_em_coupling_set_t
choose_pkpm_em_coupling_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_em_coupling_set_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_em_coupling_set_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static pkpm_em_coupling_copy_t
choose_pkpm_em_coupling_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_em_coupling_copy_kernels[cdim-1].kernels[poly_order];
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_pkpm_em_coupling_copy_kernels[cdim-1].kernels[poly_order];
      break;
    default:
      assert(false);
      break;  
  }
}
