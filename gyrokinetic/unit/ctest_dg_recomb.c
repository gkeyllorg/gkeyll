#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_recomb.h>
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
void eval_T_over_m_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4.0*echarge/emass;  
}

void
test_coll_recomb(bool use_gpu, enum gkyl_ion_type type_ion)
{
  int charge_state;
  if (type_ion == GKYL_ION_H) {  
    charge_state = 0;
  }
  else if (type_ion == GKYL_ION_LI || type_ion == GKYL_ION_AR) {
    charge_state = 1;
  }
  // use vt = 4 eV for all grids
  double vmax_elc = 4.*sqrt(4.*echarge/emass);
  double vmin_elc = -vmax_elc;
  double mumax_elc = 12.*4.*echarge/(2.*B0);
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
  gkyl_proj_on_basis *proj_T_over_m_elc = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_T_over_m_elc, NULL);

  struct gkyl_dg_recomb_inp rec_inp = {
    .grid = &phaseGrid_elc,
    .cbasis = &basis,
    .pbasis = &phaseBasis_gk,
    .conf_rng = &confRange,
    .conf_rng_ext = &confRange_ext,
    .phase_rng = &phaseRange_elc,
    .type_ion = type_ion,
    .charge_state = charge_state,
  };
  
  // coll struct.
  struct gkyl_dg_recomb *coll_recomb_up = gkyl_dg_recomb_new(&rec_inp, use_gpu);
  
  struct gkyl_array *n_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *T_over_m_elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *coef_recomb = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(proj_n_elc, 0.0, &confRange, n_elc);
  gkyl_proj_on_basis_advance(proj_T_over_m_elc, 0.0, &confRange, T_over_m_elc);
 
  gkyl_array_set_offset(moms_elc, 1.0, n_elc, 0);
  gkyl_array_set_offset(moms_elc, 1.0, T_over_m_elc, 2*basis.num_basis);

  // cuda stuff
  if (use_gpu) {
    struct gkyl_array *moms_elc_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
    struct gkyl_array *coef_recomb_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

    gkyl_array_copy(moms_elc_cu, moms_elc);
	
    gkyl_dg_recomb_coll(coll_recomb_up, moms_elc_cu, coef_recomb_cu, 0);
    gkyl_array_copy(coef_recomb, coef_recomb_cu);

    gkyl_array_release(moms_elc_cu); 
    gkyl_array_release(coef_recomb_cu);
  }
  else {
    gkyl_dg_recomb_coll(coll_recomb_up, moms_elc, coef_recomb, 0);
  }
  const double *cv_r = gkyl_array_cfetch(coef_recomb, gkyl_range_idx(&confRange, (int[3]) { 1, 1, 1}));

  // test against predicted value
  if (type_ion == GKYL_ION_H) {
    double p1_vals[] = {4.0651315620487753e-19, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00};
    for (int i=0; i<basis.num_basis; ++i) { 
      TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
    }    
  }
  else if (type_ion == GKYL_ION_LI) {
    double p1_vals[] = {1.4761368114720401e-18, 0.0000000000000000e+00, 0.0000000000000000e+00,
              0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
              0.0000000000000000e+00, 0.0000000000000000e+00};
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
    }
  }
  else if (type_ion == GKYL_ION_AR) {
    double p1_vals[] = {2.6716460249415115e-18, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00};
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_vals[i]*check_fac, cv_r[i]*check_fac, 1e-12) );
    }
  }

  gkyl_array_release(n_elc); 
  gkyl_array_release(T_over_m_elc); 
  gkyl_array_release(moms_elc); 
  gkyl_array_release(coef_recomb); 
  gkyl_proj_on_basis_release(proj_n_elc); 
  gkyl_proj_on_basis_release(proj_T_over_m_elc); 
  gkyl_dg_recomb_release(coll_recomb_up);
}

void coll_recomb_h() { test_coll_recomb(false, GKYL_ION_H); }
void coll_recomb_li() { test_coll_recomb(false, GKYL_ION_LI); }
void coll_recomb_ar() { test_coll_recomb(false, GKYL_ION_AR); }
void coll_recomb_he() { test_coll_recomb(false, GKYL_ION_HE); }
void coll_recomb_be() { test_coll_recomb(false, GKYL_ION_BE); }
void coll_recomb_b() { test_coll_recomb(false, GKYL_ION_B); }
void coll_recomb_c() { test_coll_recomb(false, GKYL_ION_C); }
void coll_recomb_n() { test_coll_recomb(false, GKYL_ION_N); }
void coll_recomb_o() { test_coll_recomb(false, GKYL_ION_O); }

#ifdef GKYL_HAVE_CUDA
void coll_recomb_h_gpu() { test_coll_recomb(true, GKYL_ION_H); }
void coll_recomb_li_gpu() { test_coll_recomb(true, GKYL_ION_LI); }
void coll_recomb_ar_gpu() { test_coll_recomb(true, GKYL_ION_AR); }
void coll_recomb_he_gpu() { test_coll_recomb(true, GKYL_ION_HE); }
void coll_recomb_be_gpu() { test_coll_recomb(true, GKYL_ION_BE); }
void coll_recomb_b_gpu() { test_coll_recomb(true, GKYL_ION_B); }
void coll_recomb_c_gpu() { test_coll_recomb(true, GKYL_ION_C); }
void coll_recomb_n_gpu() { test_coll_recomb(true, GKYL_ION_N); }
void coll_recomb_o_gpu() { test_coll_recomb(true, GKYL_ION_O); }
#endif

TEST_LIST = {
  { "coll_recomb_h", coll_recomb_h },
  { "coll_recomb_li", coll_recomb_li },
  { "coll_recomb_ar", coll_recomb_ar },
#ifdef GKYL_HAVE_CUDA
  { "coll_recomb_h_gpu", coll_recomb_h_gpu },
  { "coll_recomb_li_gpu", coll_recomb_li_gpu },
  { "coll_recomb_ar_gpu", coll_recomb_ar_gpu },
#endif
  { NULL, NULL },
};
