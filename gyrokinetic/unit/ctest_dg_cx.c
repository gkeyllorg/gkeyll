#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_cx.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double d_ion_mass = GKYL_PROTON_MASS*2.01410177811;
double B0 = 0.5;
double check_fac = 1.0e10;

void eval_n(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_T_over_m_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.*echarge/d_ion_mass;  
}
void eval_T_over_m_neut(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.*echarge/d_ion_mass;  
}

static inline
void proj_on_basis_copy(const struct gkyl_proj_on_basis *proj_op, double tm, struct gkyl_range *rng,
  struct gkyl_array *arr, bool use_gpu)
{
  struct gkyl_array *arr_ho = use_gpu? mkarr(false, arr->ncomp, arr->size)
                                     : gkyl_array_acquire(arr);
  gkyl_proj_on_basis_advance(proj_op, tm, rng, arr_ho);
  gkyl_array_copy(arr, arr_ho);
  gkyl_array_release(arr_ho);
}

// test 2x2v / 2x3v
void
test_coll_cx_d(bool use_gpu)
{
  int charge_state = 0;
  int poly_order = 1;
  int cdim = 2, vdim_vl = 3;

  // Grids 
  double lower_ion[] = {-2.0,-2.0}, upper_ion[] = {2.0,2.0};
  int ghost_gk[] = {0, 0};
  int cells_gk[] = {2, 2};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_ion, upper_ion, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange_ext, &confRange);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // Projection updater for moments.
  gkyl_proj_on_basis *proj_n = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_n, NULL);
  gkyl_proj_on_basis *proj_T_over_m_ion = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_ion, NULL);
  gkyl_proj_on_basis *proj_T_over_m_neut = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_neut, NULL);

  double vt_sq_ion_min = 1*echarge/d_ion_mass;
  double vt_sq_neut_min = 1*echarge/d_ion_mass;
  struct gkyl_dg_cx_inp cx_inp = {
    .cbasis = &basis,
    .conf_rng = &confRange,
    .vt_sq_ion_min = vt_sq_ion_min, 
    .vt_sq_neut_min = vt_sq_neut_min, 
    .type_ion = GKYL_ION_D,
  };
  
  // Coll struct.
  struct gkyl_dg_cx *coll_cx_up = gkyl_dg_cx_new(&cx_inp, use_gpu);;

  struct gkyl_array *n_ion = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_ion = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = mkarr(use_gpu, (1+2)*basis.num_basis, confRange.volume);
  struct gkyl_array *n_neut = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_neut = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = mkarr(use_gpu, (3+2)*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_cx = mkarr(use_gpu, basis.num_basis, confRange.volume);
  
  // Arrays necessary for prim_vars.
  struct gkyl_array *u_par_ion = mkarr(use_gpu, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_i = mkarr(use_gpu, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *b_x = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *b_y = mkarr(use_gpu, basis.num_basis, confRange.volume);
  struct gkyl_array *b_z = mkarr(use_gpu, basis.num_basis, confRange.volume);
  
  // Project moments on basis.
  proj_on_basis_copy(proj_n, 0.0, &confRange, n_ion, use_gpu);
  proj_on_basis_copy(proj_T_over_m_ion, 0.0, &confRange, T_over_m_ion, use_gpu);
  gkyl_array_set_offset(moms_ion, 1.0, n_ion, 0);
  gkyl_array_set_offset(moms_ion, 1.0, T_over_m_ion, (1+1)*basis.num_basis);

  proj_on_basis_copy(proj_n, 0.0, &confRange, n_neut, use_gpu);
  proj_on_basis_copy(proj_T_over_m_neut, 0.0, &confRange, T_over_m_neut, use_gpu);
  gkyl_array_set_offset(moms_neut, 1.0, n_neut, 0);
  gkyl_array_set_offset(moms_neut, 1.0, T_over_m_neut, (1+vdim_vl)*basis.num_basis);

  // Project b_i.
  gkyl_array_set_offset(b_i, 1.0, b_x, 0);
  gkyl_array_set_offset(b_i, 1.0, b_y, basis.num_basis);
  gkyl_array_set_offset(b_i, 1.0, b_z, 2*basis.num_basis);

  // Compute the CX reaction rate.
  gkyl_dg_cx_coll(coll_cx_up, moms_ion, moms_neut, u_par_ion, coef_cx, 0);    

  struct gkyl_array *coef_cx_ho = use_gpu? mkarr(false, coef_cx->ncomp, coef_cx->size)
                                         : gkyl_array_acquire(coef_cx);
  gkyl_array_copy(coef_cx_ho, coef_cx);

  const double *cv_cx = gkyl_array_cfetch(coef_cx_ho, gkyl_range_idx(&confRange, (int[2]) { 1, 1}));

  // Test against predicted value.
  double p1_vals[] = {3.242709205939892e-14, 0.000000000000000e+00,
		      0.000000000000000e+00, 0.000000000000000e+00};
  for (int i=0; i<basis.num_basis; ++i) {
    TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_cx[i]*check_fac, 1e-12) );
    TEST_MSG( "i:%d | Expected: %.9e | Got:%.9e\n", i, p1_vals[i]*check_fac, cv_cx[i]*check_fac );
  }

  gkyl_array_release(n_ion); 
  gkyl_array_release(T_over_m_ion); 
  gkyl_array_release(n_neut); 
  gkyl_array_release(T_over_m_neut);
  gkyl_array_release(moms_neut); 
  gkyl_array_release(moms_ion);
  gkyl_array_release(u_par_ion);
  gkyl_array_release(b_i); 
  gkyl_array_release(b_x); 
  gkyl_array_release(b_y); 
  gkyl_array_release(b_z);
  gkyl_array_release(coef_cx);
  gkyl_array_release(coef_cx_ho);
  gkyl_proj_on_basis_release(proj_n); 
  gkyl_proj_on_basis_release(proj_T_over_m_ion); 
  gkyl_proj_on_basis_release(proj_T_over_m_neut); 
  gkyl_dg_cx_release(coll_cx_up); 
}


void coll_cx_d() { test_coll_cx_d(false); }

#ifdef GKYL_HAVE_CUDA
void coll_cx_d_gpu() { test_coll_cx_d(true); }
#endif

TEST_LIST = {
  { "coll_cx_d", coll_cx_d },
#ifdef GKYL_HAVE_CUDA
  { "coll_cx_d_gpu", coll_cx_d_gpu },
#endif
  { NULL, NULL },
};
