// A test of calculation of moments of a distribution function.
//
#include <acutest.h>
#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>

static inline double
maxwellian1D(double n, double vx, double ux, double vth)
{
  double v2 = (vx-ux)*(vx-ux);
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

static inline double
maxwellian2D(double n, double vx, double vy, double ux, double uy, double vth)
{
  double v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy);
  return n/(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc1x1v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  
  fout[0] = maxwellian1D(1.0, vx, 0.0, 1.0);
}

void
evalDistFunc1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], vx = xn[1], vy = xn[2];
  
  fout[0] = maxwellian2D(1.0, vx, vy, 0.0, 0.0, 1.0);
}

void nu_prof_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  fout[0] = 1.0;
}

void nu_prof_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double vy = xn[2];
  fout[0] = 1.0;
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void
test_func(int cdim, int vdim, int poly_order, 
  evalf_t evalDistFunc, double f_check[], double vf_check[], 
  double u_check[], double vth_check[], double ucross_check[], double vthcross_check[])
{
  int pdim = cdim + vdim;  
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM], confLower[GKYL_MAX_DIM], confUpper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], confCells[GKYL_MAX_DIM];
  
  double v_bounds[2*vdim];
  
  for (int i=0; i<pdim; ++i) {
    if (i<cdim) {
      cells[i] = 4;
      lower[i] = 0.0;
      upper[i] = 1.0;
      confLower[i] = lower[i];
      confUpper[i] = upper[i];
      confCells[i] = cells[i];
    } else {
      cells[i] = 24;
      lower[i] = -2.0;
      upper[i] = 2.0;
      v_bounds[i - cdim] = lower[i];
      v_bounds[i - cdim + vdim] = upper[i];
    }
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[GKYL_MAX_DIM];
  for (int i=0; i<pdim; ++i) {
    if (i<cdim) {
      ghost[i] = confGhost[i];
    } else {
      ghost[i] = 0;
    }
  }
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalDistFunc, NULL);
  
  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // projection updater for collision frequency
  gkyl_proj_on_basis *projNu;
  if (cdim==1 && vdim==1) {
    projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis,
      poly_order+1, 1, nu_prof_1x1v, NULL);
  } else if (cdim==1 && vdim==2) {
    projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis,
      poly_order+1, 1, nu_prof_1x2v, NULL);
  }

  // create collision frequency array
  struct gkyl_array *nu;
  nu = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // project collision frequency on basis
  gkyl_proj_on_basis_advance(projNu, 0.0, &confLocal_ext, nu);
  
  struct gkyl_mom_type *vm_moms_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0M1M2, false);
  gkyl_mom_calc *moms_calc = gkyl_mom_calc_new(&grid, vm_moms_t, false);

  // create moment arrays
  struct gkyl_array *moms;
  moms = mkarr(confBasis.num_basis*(vdim+2), confLocal_ext.volume);
  
  // compute the moments
  gkyl_mom_calc_advance(moms_calc, &local, &confLocal, distf, moms);

  gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(&grid, &confBasis, &basis, v_bounds, false);
  
  // create moment arrays
  struct gkyl_array *boundary_corrections;
  boundary_corrections = mkarr((vdim+1)*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &local, &confLocal, distf, boundary_corrections);
  
  // Check boundary corrections of momentum and energy.
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *fptr = gkyl_array_fetch(boundary_corrections, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( f_check[k], fptr[k], 1e-12) );
      TEST_CHECK( gkyl_compare( vf_check[k], fptr[k+vdim*confBasis.num_basis], 1e-12) );
    }
  }

  gkyl_prim_lbo_calc *primcalc = gkyl_prim_lbo_vlasov_calc_new(&grid, &confBasis, &basis, &confLocal, false);
  gkyl_prim_lbo_cross_calc *crossprimcalc = gkyl_prim_lbo_vlasov_cross_calc_new(&grid, &confBasis, &basis, &confLocal, false);

  const struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_calc_get_prim(primcalc);

  TEST_CHECK( prim->cdim == cdim );
  TEST_CHECK( prim->pdim == pdim );
  TEST_CHECK( prim->poly_order == poly_order );
  TEST_CHECK( prim->num_config == confBasis.num_basis );
  TEST_CHECK( prim->num_phase == basis.num_basis );
  
  // create moment arrays
  struct gkyl_array *u, *vth, *prim_moms;
  u = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  vth = mkarr(confBasis.num_basis, confLocal_ext.volume);
  prim_moms = mkarr((vdim+1)*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_prim_lbo_calc_advance(primcalc, &confLocal, moms, boundary_corrections, nu, prim_moms);

  gkyl_array_set_offset(u, 1., prim_moms, 0);
  gkyl_array_set_offset(vth, 1., prim_moms, vdim*confBasis.num_basis);

  // Check u
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( u_check[k], uptr[k], 1e-12) );
  }}

  // Check vtSq.
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vth, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( vth_check[k], vthptr[k], 1e-12) );
  }}

  struct gkyl_array *u_out = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_out = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *prim_moms_out = mkarr((vdim+1)*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *greene = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_clear(greene, 1.0);
  struct gkyl_array *cross_prim_moms = prim_moms;
  double self_m = 1.;
  double cross_m = self_m;

  // MF 2022/09/13: the second moms here should be cross_moms, but we pass moms
  // for simplicity in this (infrastructure) test.
  gkyl_prim_lbo_cross_calc_advance(crossprimcalc, &confLocal, greene,
    self_m, moms, prim_moms, cross_m, moms, cross_prim_moms,
    boundary_corrections, nu, prim_moms_out);
  
  gkyl_array_set_offset(u_out, 1., prim_moms_out, 0);
  gkyl_array_set_offset(vtsq_out, 1., prim_moms_out, vdim*confBasis.num_basis);

  // Check cross u
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u_out, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( ucross_check[k], uptr[k], 1e-12) );
  }}

  // Check cross vtsq
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vtsq_out, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( vthcross_check[k], vthptr[k], 1e-12) );
      TEST_MSG("Expected: %.13e in cell (%d)", vthcross_check[k], i);
      TEST_MSG("Produced: %.13e", vthptr[k]);
  }}

  gkyl_array_release(moms);
  gkyl_mom_calc_release(moms_calc);
  gkyl_mom_type_release(vm_moms_t);

  gkyl_array_release(boundary_corrections); 
  gkyl_mom_calc_bcorr_release(bcorr_calc); 

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);

  gkyl_proj_on_basis_release(projNu);
  gkyl_array_release(nu);

  gkyl_array_release(u); gkyl_array_release(vth);
  gkyl_array_release(prim_moms);
  gkyl_prim_lbo_calc_release(primcalc);

  gkyl_array_release(u_out); gkyl_array_release(vtsq_out);
  gkyl_array_release(prim_moms_out);
  gkyl_array_release(greene);
  gkyl_prim_lbo_cross_calc_release(crossprimcalc);

}

#ifdef GKYL_HAVE_CUDA
void
test_func_cu(int cdim, int vdim, int poly_order, 
  evalf_t evalDistFunc, double f_check[], double vf_check[], 
  double u_check[], double vth_check[], double ucross_check[], double vthcross_check[])
{
  int pdim = cdim + vdim;  
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM], confLower[GKYL_MAX_DIM], confUpper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], confCells[GKYL_MAX_DIM];
  
  double v_bounds[2*vdim];
  
  for (int i=0; i<pdim; ++i) {
    if (i<cdim) {
      cells[i] = 4;
      lower[i] = 0.0;
      upper[i] = 1.0;
      confLower[i] = lower[i];
      confUpper[i] = upper[i];
      confCells[i] = cells[i];
    } else {
      cells[i] = 24;
      lower[i] = -2.0;
      upper[i] = 2.0;
      v_bounds[i - cdim] = lower[i];
      v_bounds[i - cdim + vdim] = upper[i];
    }
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[GKYL_MAX_DIM];
  for (int i=0; i<pdim; ++i) {
    if (i<cdim) {
      ghost[i] = confGhost[i];
    } else {
      ghost[i] = 0;
    }
  }
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalDistFunc, NULL);

  // create distribution function array
  struct gkyl_array *distf, *distf_cu;
  distf = mkarr(basis.num_basis, local_ext.volume);
  distf_cu = mkarr_cu(basis.num_basis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  gkyl_array_copy(distf_cu, distf);

  // projection updater for collision frequency
  gkyl_proj_on_basis *projNu;
  if (cdim==1 && vdim==1) {
    projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis,
      poly_order+1, 1, nu_prof_1x1v, NULL);
  } else if (cdim==1 && vdim==2) {
    projNu = gkyl_proj_on_basis_new(&confGrid, &confBasis,
      poly_order+1, 1, nu_prof_1x2v, NULL);
  }

  // create collision frequency array
  struct gkyl_array *nu, *nu_cu;
  nu = mkarr(confBasis.num_basis, confLocal_ext.volume);
  nu_cu = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // project collision frequency on basis
  gkyl_proj_on_basis_advance(projNu, 0.0, &confLocal_ext, nu);
  gkyl_array_copy(nu_cu, nu);
  
  struct gkyl_mom_type *vm_moms_t = gkyl_mom_vlasov_new(&confBasis, &basis, GKYL_F_MOMENT_M0M1M2, true);
  gkyl_mom_calc *moms_calc = gkyl_mom_calc_new(&grid, vm_moms_t, true);

  // create moment arrays
  struct gkyl_array *moms_cu;
  moms_cu = mkarr_cu(confBasis.num_basis*(vdim+2), confLocal_ext.volume);
  
  // compute the moments
  gkyl_mom_calc_advance_cu(moms_calc, &local, &confLocal, distf_cu, moms_cu);

  gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_lbo_vlasov_new(&grid, &confBasis, &basis, v_bounds, true);
  
  // create moment arrays
  struct gkyl_array *boundary_corrections_cu;
  boundary_corrections_cu = mkarr_cu((vdim+1)*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &local, &confLocal, distf_cu, boundary_corrections_cu);

  gkyl_prim_lbo_calc *primcalc = gkyl_prim_lbo_vlasov_calc_new(&grid, &confBasis, &basis, &confLocal, true);
  gkyl_prim_lbo_cross_calc *crossprimcalc = gkyl_prim_lbo_vlasov_cross_calc_new(&grid, &confBasis, &basis, &confLocal, true);

  // create moment arrays
  struct gkyl_array *u, *vth, *u_cu, *vth_cu, *prim_moms_cu;
  u = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  vth = mkarr(confBasis.num_basis, confLocal_ext.volume);
  u_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  vth_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  prim_moms_cu = mkarr_cu((vdim+1)*confBasis.num_basis, confLocal_ext.volume);

  // compute the moment corrections
  gkyl_prim_lbo_calc_advance(primcalc, &confLocal, moms_cu, boundary_corrections_cu, nu_cu, prim_moms_cu);

  gkyl_array_set_offset(u_cu, 1., prim_moms_cu, 0);
  gkyl_array_set_offset(vth_cu, 1., prim_moms_cu, vdim*confBasis.num_basis);

  gkyl_array_copy(u, u_cu);
  gkyl_array_copy(vth, vth_cu);

  // Check u
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( u_check[k], uptr[k], 1e-12) );
  }}

  // Check vtSq.
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vth, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( vth_check[k], vthptr[k], 1e-12) );
  }}

  struct gkyl_array *u_out = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_out = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *u_out_cu = mkarr_cu(vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq_out_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *prim_moms_out_cu = mkarr_cu((vdim+1)*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *greene_cu = mkarr_cu(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *cross_prim_moms = prim_moms_cu;
  double self_m = 1.;
  double cross_m = self_m;
  gkyl_array_clear(greene_cu, 1.0);

  // MF 2022/09/13: the second moms here should be cross_moms, but we pass moms
  // for simplicity in this (infrastructure) test.
  gkyl_prim_lbo_cross_calc_advance(crossprimcalc, &confLocal, greene_cu,
    self_m, moms_cu, prim_moms_cu, cross_m, moms_cu, cross_prim_moms,
    boundary_corrections_cu, nu_cu, prim_moms_out_cu);

  gkyl_array_set_offset(u_out_cu, 1., prim_moms_out_cu, 0);
  gkyl_array_set_offset(vtsq_out_cu, 1., prim_moms_out_cu, vdim*confBasis.num_basis);

  gkyl_array_copy(u_out, u_out_cu);
  gkyl_array_copy(vtsq_out, vtsq_out_cu);
  
  // Check cross u
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *uptr = gkyl_array_fetch(u_out, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( ucross_check[k], uptr[k], 1e-12) );
  }}

  // Check cross vtsq
  // 1-indexed for interfacing with G2 Lua layer
  for (unsigned int i=1; i<cells[0]+1; ++i) {
    int cidx[] = {i};
    long linc = gkyl_range_idx(&confLocal, cidx);
    double *vthptr = gkyl_array_fetch(vtsq_out, linc);
    for (unsigned int k=0; k<confBasis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( vthcross_check[k], vthptr[k], 1e-12) );
  }}

  // release memory for objects
  gkyl_array_release(moms_cu);
  gkyl_mom_calc_release(moms_calc);
  gkyl_mom_type_release(vm_moms_t);

  gkyl_array_release(boundary_corrections_cu); 
  gkyl_mom_calc_bcorr_release(bcorr_calc); 

  gkyl_proj_on_basis_release(projNu);
  gkyl_array_release(nu); gkyl_array_release(nu_cu);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf); gkyl_array_release(distf_cu);

  gkyl_array_release(u); gkyl_array_release(vth);
  gkyl_array_release(u_cu); gkyl_array_release(vth_cu);
  gkyl_array_release(prim_moms_cu);
  gkyl_prim_lbo_calc_release(primcalc);

  gkyl_array_release(u_out); gkyl_array_release(vtsq_out);
  gkyl_array_release(u_out_cu); gkyl_array_release(vtsq_out_cu);
  gkyl_array_release(prim_moms_out_cu);
  gkyl_array_release(greene_cu);
  gkyl_prim_lbo_cross_calc_release(crossprimcalc);
}
#endif

void
test_prim_vlasov_1x1v_p2_ho()
{
  int poly_order = 2;
  int vdim = 1, cdim = 1;

  double f_check[] = { 0.0, 0.0, 0.0 };
  double vf_check[] = { 0.30543841971927, 0.0, 0.0 };
  double u_check[] = { 0.0, 0.0, 0.0 };
  double vth_check[] = { 1.4142398195471544, 0.0, 0.0 };
  double ucross_check[] = { 0.0, 0.0, 0.0 };
  double vthcross_check[] = { 1.4142398195471544, 0.0, 0.0 };

  test_func(cdim, vdim, poly_order, evalDistFunc1x1v, f_check, vf_check, u_check, vth_check, ucross_check, vthcross_check);
}

void
test_prim_vlasov_1x2v_p2_ho()
{
  int poly_order = 2;
  int vdim = 2, cdim = 1;

  double f_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vf_check[] = { 0.583081782023233, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double u_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vth_check[] = { 1.4142398195471586, 0.0, 0.0 };
  double ucross_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vthcross_check[] = { 1.4142398195471586, 0.0, 0.0 };

  test_func(cdim, vdim, poly_order, evalDistFunc1x2v, f_check, vf_check, u_check, vth_check, ucross_check, vthcross_check);
}

#ifdef GKYL_HAVE_CUDA
void
test_prim_vlasov_1x1v_p2_dev()
{
  int poly_order = 2;
  int vdim = 1, cdim = 1;

  double f_check[] = { 0.0, 0.0, 0.0 };
  double vf_check[] = { 0.30543841971927, 0.0, 0.0 };
  double u_check[] = { 0.0, 0.0, 0.0 };
  double vth_check[] = { 1.4142398195471544, 0.0, 0.0 };
  double ucross_check[] = { 0.0, 0.0, 0.0 };
  double vthcross_check[] = { 1.4142398195471544, 0.0, 0.0 };

  test_func_cu(cdim, vdim, poly_order, evalDistFunc1x1v, f_check, vf_check, u_check, vth_check, ucross_check, vthcross_check);
}

void
test_prim_vlasov_1x2v_p2_dev()
{
  int poly_order = 2;
  int vdim = 2, cdim = 1;

  double f_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vf_check[] = { 0.583081782023233, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double u_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vth_check[] = { 1.4142398195471586, 0.0, 0.0 };
  double ucross_check[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double vthcross_check[] = { 1.4142398195471586, 0.0, 0.0 };


  test_func_cu(cdim, vdim, poly_order, evalDistFunc1x2v, f_check, vf_check, u_check, vth_check, ucross_check, vthcross_check);
}
#endif

TEST_LIST = {
  { "test_prim_vlasov_1x1v_p2_ho", test_prim_vlasov_1x1v_p2_ho },
  { "test_prim_vlasov_1x2v_p2_ho", test_prim_vlasov_1x2v_p2_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_prim_vlasov_1x1v_p2_dev", test_prim_vlasov_1x1v_p2_dev },
  { "test_prim_vlasov_1x2v_p2_dev", test_prim_vlasov_1x2v_p2_dev },
#endif
  { NULL, NULL },
};
