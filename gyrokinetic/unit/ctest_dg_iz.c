#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_iz.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;
double check_fac = 1.e10;
double B0 = 0.5;

void eval_n_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_T_over_m_elc_40ev(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40.0*echarge/emass;  
}
void eval_T_over_m_elc_100ev(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 100.0*echarge/emass;  
}

void
test_coll_iz(bool use_gpu, enum gkyl_ion_type type_ion)
{
  int charge_state;
  if (type_ion == GKYL_ION_H) {  
    charge_state = 0;
  }
  else if (type_ion == GKYL_ION_LI) {
    charge_state = 1;
  }
  // use vt = 40 eV for all grids
  double vmax_elc = 4.*sqrt(40.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*40.*echarge/(2.*B0);
  int poly_order = 1;
  const int cdim = 3, vdim_gk = 2;
  int pdim_gk = cdim + vdim_gk;

  // for gk grids 
  double lower_elc[] = {-2.0,-2.0,-2.0,vmin_elc,0.0}, upper_elc[] = {2.0,2.0,2.0,vmax_elc,mumax_elc};
  int ghost_gk[] = {0, 0, 0, 0, 0};
  int cells_gk[] = {16, 16, 16, 8, 4};
  
  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&confGrid, ghost_gk, &confRange_ext, &confRange);

  // elc phase grid
  struct gkyl_rect_grid phaseGrid_elc;
  struct gkyl_range phaseRange_elc, phaseRange_ext_elc;
  gkyl_rect_grid_init(&phaseGrid_elc, pdim_gk, lower_elc, upper_elc, cells_gk);
  gkyl_create_grid_ranges(&phaseGrid_elc, ghost_gk, &phaseRange_ext_elc, &phaseRange_elc);

  // basis functions
  struct gkyl_basis phaseBasis_vl, phaseBasis_gk, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_gkhybrid(&phaseBasis_gk, cdim, vdim_gk);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *proj_n_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_n_elc, NULL);
  gkyl_proj_on_basis *proj_T_over_m_elc_40ev = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_elc_40ev, NULL);
  gkyl_proj_on_basis *proj_T_over_m_elc_100ev = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_elc_100ev, NULL);

  struct gkyl_dg_iz_inp iz_inp = {
    .cbasis = &basis,
    .conf_rng = &confRange,
    .type_ion = type_ion,
    .charge_state = charge_state,
    .type_self = GKYL_SELF_ELC,
  };
  
  // coll struct.
  struct gkyl_dg_iz *coll_iz_up = gkyl_dg_iz_new(&iz_inp, use_gpu);
  
  struct gkyl_array *n_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_iz1 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *vtSq_iz2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *coef_iz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(proj_n_elc, 0.0, &confRange, n_elc);
  if (type_ion == GKYL_ION_H) {
    gkyl_proj_on_basis_advance(proj_T_over_m_elc_40ev, 0.0, &confRange, T_over_m_elc);
  }
  else {
    gkyl_proj_on_basis_advance(proj_T_over_m_elc_100ev, 0.0, &confRange, T_over_m_elc);
  }
 
  gkyl_array_set_offset(moms_elc, 1.0, n_elc, 0);
  gkyl_array_set_offset(moms_elc, 1.0, T_over_m_elc, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *vtSq_iz1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *vtSq_iz2_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
    struct gkyl_array *coef_iz_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_elc_cu, moms_elc);
	
    gkyl_dg_iz_coll(coll_iz_up, moms_elc_cu, vtSq_iz1_cu, vtSq_iz2_cu, coef_iz_cu, 0);
    gkyl_array_copy(coef_iz, coef_iz_cu);

    gkyl_array_release(moms_elc_cu); 
    gkyl_array_release(vtSq_iz1_cu); 
    gkyl_array_release(vtSq_iz2_cu);
    gkyl_array_release(coef_iz_cu);
  }
  else {
    gkyl_dg_iz_coll(coll_iz_up, moms_elc, vtSq_iz1, vtSq_iz2, coef_iz, 0);
  }
  const double *cv_iz = gkyl_array_cfetch(coef_iz, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  // test against predicted value
  if (type_ion == GKYL_ION_H) {
    double p1_vals[] = {8.3873695794923268e-14, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00};
    for (int i=0; i<basis.num_basis; ++i) { 
      TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_iz[i]*check_fac, 1e-12) );
    }    
  }
  else if (type_ion == GKYL_ION_LI) {
    double p1_vals[] = {4.7212386637934834e-15, 0.0000000000000000e+00, 0.0000000000000000e+00,
              0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
              0.0000000000000000e+00, 0.0000000000000000e+00};
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_iz[i]*check_fac, 1e-12) );
    }
  }

  gkyl_array_release(n_elc); 
  gkyl_array_release(T_over_m_elc); 
  gkyl_array_release(moms_elc); 
  gkyl_array_release(vtSq_iz1); 
  gkyl_array_release(vtSq_iz2); 
  gkyl_array_release(coef_iz); 
  gkyl_proj_on_basis_release(proj_n_elc); 
  gkyl_proj_on_basis_release(proj_T_over_m_elc_40ev); 
  gkyl_proj_on_basis_release(proj_T_over_m_elc_100ev); 
  gkyl_dg_iz_release(coll_iz_up);
}

void coll_iz_h() { test_coll_iz(false, GKYL_ION_H); }
void coll_iz_li() { test_coll_iz(false, GKYL_ION_LI); }
void coll_iz_ar() { test_coll_iz(false, GKYL_ION_AR); }
void coll_iz_he() { test_coll_iz(false, GKYL_ION_HE); }
void coll_iz_be() { test_coll_iz(false, GKYL_ION_BE); }
void coll_iz_b() { test_coll_iz(false, GKYL_ION_B); }
void coll_iz_c() { test_coll_iz(false, GKYL_ION_C); }
void coll_iz_n() { test_coll_iz(false, GKYL_ION_N); }
void coll_iz_o() { test_coll_iz(false, GKYL_ION_O); }

#ifdef GKYL_HAVE_CUDA
void coll_iz_h_gpu() { test_coll_iz(true, GKYL_ION_H); }
void coll_iz_li_gpu() { test_coll_iz(true, GKYL_ION_LI); }
void coll_iz_ar_gpu() { test_coll_iz(true, GKYL_ION_AR); }
void coll_iz_he_gpu() { test_coll_iz(true, GKYL_ION_HE); }
void coll_iz_be_gpu() { test_coll_iz(true, GKYL_ION_BE); }
void coll_iz_b_gpu() { test_coll_iz(true, GKYL_ION_B); }
void coll_iz_c_gpu() { test_coll_iz(true, GKYL_ION_C); }
void coll_iz_n_gpu() { test_coll_iz(true, GKYL_ION_N); }
void coll_iz_o_gpu() { test_coll_iz(true, GKYL_ION_O); }
#endif

TEST_LIST = {
  { "coll_iz_h", coll_iz_h },
  { "coll_iz_li", coll_iz_li },
#ifdef GKYL_HAVE_CUDA
  { "coll_iz_h_gpu", coll_iz_h_gpu },
  { "coll_iz_li_gpu", coll_iz_li_gpu },
#endif
  { NULL, NULL },
};
