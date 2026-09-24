#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_cx_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

typedef double (*dg_cx_react_ratef_t)(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, 
  const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion,
  double* GKYL_RESTRICT v_sigma_cx) ;

// for use in kernel tables
typedef struct { dg_cx_react_ratef_t kernels[3]; } gkyl_cx_react_rate_kern_list;

// CX reaction rate kernel list

//
// Serendipity basis kernels
// 
GKYL_CU_D
static const gkyl_cx_react_rate_kern_list ser_cx_react_rate_kernels[] = {
  { sigma_cx_1x_ser_p1, sigma_cx_1x_ser_p2 }, // 0
  { sigma_cx_2x_ser_p1, sigma_cx_2x_ser_p2 }, // 4
  { sigma_cx_3x_ser_p1, NULL }, // 5
};

struct gkyl_dg_cx {
  struct gkyl_basis *cbasis;
  const struct gkyl_range *conf_rng;
  
  double a; // Fitting function coefficient.
  double b; // Fitting function coefficient.
  double vt_sq_ion_min;
  double vt_sq_neut_min; 
  
  enum gkyl_ion_type type_ion;

  uint32_t flags;
  dg_cx_react_ratef_t react_rate; // pointer to reaction rate kernel

  struct gkyl_dg_cx *on_dev; // pointer to itself or device data

};

GKYL_CU_D
static dg_cx_react_ratef_t
choose_kern(struct gkyl_basis cbasis)
{
  int cdim = cbasis.ndim;
  int poly_order = cbasis.poly_order;
  enum gkyl_basis_type b_type = cbasis.b_type;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_cx_react_rate_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break; 
  }
  return 0;
}

static void
fit_param(enum gkyl_ion_type type_ion, double *a, double *b)
{
  // These values are from E. Meier's PhD Thesis
  if (type_ion == GKYL_ION_H) {
    a[0] = 1.12e-18;
    b[0] = 7.15e-20;
  }
  else if (type_ion == GKYL_ION_D) {
    a[0] = 1.09e-18;
    b[0] = 7.15e-20;
  }
  else if (type_ion == GKYL_ION_HE) {
    a[0] = 6.484e-19;
    b[0] = 4.350e-20;
  } 
  else if (type_ion == GKYL_ION_NE) {
    a[0] = 7.95e-19;
    b[0] = 5.65e-20;
  }
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new charge exchange updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_cx* gkyl_dg_cx_cu_dev_new(struct gkyl_dg_cx_inp *inp);

/**
 * Compute CX reaction rate coefficient for use in neutral reactions
 * on the NVIDIA GPU. The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param cx charge exchange object.
 * @param maxwellian_moms_ion  Ion Maxwellian moments (n, upar, T/m).
 * @param maxwellian_moms_neut Neutral Maxwellian moments (n, ux, uy, uz, T/m).
 * @param upar_b_i Ion drift velocity vector (upar b_x, upar b_y, upar b_z).
 * @param coef_cx Output reaction rate coefficient
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T]) 
 */
void gkyl_dg_cx_coll_cu(const struct gkyl_dg_cx *up, 
  struct gkyl_array *maxwellian_moms_ion, struct gkyl_array *maxwellian_moms_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate);
#endif
