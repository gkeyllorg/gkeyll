// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fluid_em_coupling_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*fluid_em_coupling_set_t)(
  int count, int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, bool is_field_static,
  double dt, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *app_accel[GKYL_MAX_SPECIES],
  const double *ext_em, const double *app_current, const double *sigma,
  const double *fluid[GKYL_MAX_SPECIES], const double *em
);
typedef void (*fluid_em_coupling_copy_t)(
  int count, int num_species, double qbym[GKYL_MAX_SPECIES], bool has_energy[GKYL_MAX_SPECIES],
  double epsilon0, double dt, const double *sigma, struct gkyl_nmat *x,
  double *GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double *GKYL_RESTRICT em
);

struct gkyl_dg_calc_fluid_em_coupling {
  struct gkyl_range mem_range; // Configuration space range for linear solve

  struct gkyl_nmat *As, *xs; // LHS and RHS of the linear system of each node of each cell
  gkyl_nmat_mem *mem; // memory for use in batched linear solve

  fluid_em_coupling_set_t fluid_em_coupling_set; // kernel for setting matrices for linear solve
  fluid_em_coupling_copy_t fluid_em_coupling_copy; // kernel for copying solution to output

  bool is_field_static; // true if the EM fields are static: the species currents do not advance E

  int num_species; // number of species being implicitly solved for
  double qbym[GKYL_MAX_SPECIES]; // charge/mass ratio for each species
  bool has_energy[GKYL_MAX_SPECIES]; // whether the state of each species carries an energy
  double epsilon0; // permittivity of free space

  uint32_t flags;
  struct gkyl_dg_calc_fluid_em_coupling *on_dev; // pointer to itself or device data
};

// Nodal locally implicit source kernels (p=1 serendipity only; indexed [cdim-1])
GKYL_CU_D static const fluid_em_coupling_set_t ser_fluid_em_coupling_set_kernels[] = {
  fluid_em_coupling_set_1x_ser_p1, fluid_em_coupling_set_2x_ser_p1, fluid_em_coupling_set_3x_ser_p1
};
GKYL_CU_D static const fluid_em_coupling_copy_t ser_fluid_em_coupling_copy_kernels[] = {
  fluid_em_coupling_copy_1x_ser_p1, fluid_em_coupling_copy_2x_ser_p1,
  fluid_em_coupling_copy_3x_ser_p1
};

GKYL_CU_D static fluid_em_coupling_set_t
choose_fluid_em_coupling_set_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // The tensor and serendipity p = 1 bases are the same.
  assert(
    poly_order == 1 && (b_type == GKYL_BASIS_MODAL_SERENDIPITY || b_type == GKYL_BASIS_MODAL_TENSOR)
  );
  return ser_fluid_em_coupling_set_kernels[cdim - 1];
}

GKYL_CU_D static fluid_em_coupling_copy_t
choose_fluid_em_coupling_copy_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // The tensor and serendipity p = 1 bases are the same.
  assert(
    poly_order == 1 && (b_type == GKYL_BASIS_MODAL_SERENDIPITY || b_type == GKYL_BASIS_MODAL_TENSOR)
  );
  return ser_fluid_em_coupling_copy_kernels[cdim - 1];
}
