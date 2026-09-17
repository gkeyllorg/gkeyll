// Tests for the Vlasov FPO routines to compute moments and boundary corrections
// and resulting corrections to drag and diffusion coefficients.
//
// Included are tests for p=1 hybrid and p=2 serendipity.
// As shown in Rodman 2025 PhD Thesis, the p=2 computation enforces
// conservation of momentum and energy, so this result is taken as the
// "correct" result for the hybrid basis.

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>

#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_fpo_vlasov_coeff_correct.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_fpo_vlasov.h>
#include <gkyl_dg_updater_moment.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

// Allocate array
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  if (use_gpu) {
    struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
    return a;
  }
  else {
    struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
    return a;
  }
}

static inline double
bump_maxwellian(double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA, double bUx, double bUy, double bUz, double bS, double bVt)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  double bv2 = (vx - bUx)*(vx - bUx) + (vy - bUy)*(vy - bUy) + (vz - bUz)*(vz - bUz);
  return n/pow(sqrt(2*M_PI*vt*vt), 3)*exp(-v2/(2*vt*vt)) + n/pow(sqrt(2*M_PI*bVt*bVt), 3)*exp(-bv2/(2*bVt*bVt))*(bA*bA)/(bv2 + bS*bS);
}

void
eval_distf_square(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];

  double width = 2.0; // corresponds to a final vth of 2/3
  if(vx>-width && vx<width && vy>-width && vy<width && vz>-width && vz<width) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_distf_bump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];

  fout[0] = bump_maxwellian(1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vth0, 
    sqrt(0.15), 4.0*app->vth0, 0.0, 0.0, 0.14, 3.0*app->vth0);
}

void
eval_drag_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vy*vz;
  fout[1] = vx*vz;
  fout[2] = vy*vz;
}

void
eval_diff_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vx*vx;
  fout[1] = vx*vy;
  fout[2] = vx*vz;
  fout[3] = vy*vx;
  fout[4] = vy*vy;
  fout[5] = vz*vz;
  fout[6] = vz*vx;
  fout[7] = vz*vy;
  fout[8] = vz*vz;
}

struct fpo_ctx
create_ctx() {
  struct fpo_ctx ctx = {
    .n0 = 1.0,
    .ux0 = 0.0,
    .uy0 = 0.0,
    .uz0 = 0.0,
    .vth0 = 1.0,
    .gamma0 = 10.0,
  };
return ctx;
}

void test_1x3v_square(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_square = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_square, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  struct gkyl_array *fpo_moms_ho, *boundary_corrections_ho, *drag_diff_coeff_corrs_ho;
  fpo_moms_ho = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);
  boundary_corrections_ho = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);
  drag_diff_coeff_corrs_ho = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, use_gpu);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, use_gpu);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, use_gpu);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, use_gpu);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, use_gpu);

  // Project distribution function and coefficients and copy to device if using GPU
  if (use_gpu) {
    struct gkyl_array *distf_ho, *drag_coeff_ho, *diff_coeff_ho;
    distf_ho = mkarr(phase_basis.num_basis, phase_range_ext.volume, false);
    drag_coeff_ho = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, false);
    diff_coeff_ho = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, false);

    gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf_ho);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff_ho);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff_ho);

    gkyl_array_copy(distf, distf_ho);
    gkyl_array_copy(drag_coeff, drag_coeff_ho);
    gkyl_array_copy(diff_coeff, diff_coeff_ho);

    gkyl_array_release(distf_ho);
    gkyl_array_release(drag_coeff_ho);
    gkyl_array_release(diff_coeff_ho);
  }
  else {
    gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  }

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  else {
    gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, use_gpu);

  const double *fpo_moms_c, *bcorr_c, *coeff_corrs_c;
  if (use_gpu) {
    gkyl_array_copy(fpo_moms_ho, fpo_moms);
    gkyl_array_copy(boundary_corrections_ho, boundary_corrections);
    gkyl_array_copy(drag_diff_coeff_corrs_ho, drag_diff_coeff_corrs);

    fpo_moms_c = gkyl_array_cfetch(fpo_moms_ho, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections_ho, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs_ho, 1);
  }
  else {
    fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);
  }

  if ((poly_order == 1) && (NV == 4)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 2.8421709430e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 1.9919077569e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 2.1316282073e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 1.0361524557e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -8.6053253606e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 7.5424723327e+02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -1.0515868383e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -8.8817841970e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -6.2247117404e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -6.6613381478e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -3.2379764239e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 1.9844966563e-34, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 2.6891641752e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -7.8567420132e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 5.1560651283e-17, 1e-10) );
  }

  if ((poly_order == 1) && (NV == 8)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -2.0428103653e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -4.3478795219e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 3.5527136788e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -5.0633927005e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 8.8817841970e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -3.4525448362e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 8.2212948426e+02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -7.6594836702e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 6.3837823916e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.3587123506e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.1102230246e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 1.5823102189e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], -2.7755575616e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 1.0789202613e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -8.5638487944e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 1.6588475860e-17, 1e-10) );
  }

  if ((poly_order == 2) && (NV == 4)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -1.7763568394e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -6.3528162176e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 1.7215963393e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -7.1054273576e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.5719302273e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 3.6965525345e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -1.4210854715e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -3.3883142039e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 5.5511151231e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.9852550680e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -5.3799885604e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.2204460493e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 4.9122819603e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -1.1551726670e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 4.4408920985e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 1.0588481887e-16, 1e-10) );
  }

  if ((poly_order == 2) && (NV == 8)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -6.4837024638e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 1.2754784269e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 1.9415716862e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.6645352591e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.6729638855e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -3.2907815638e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -6.2172489379e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -1.6409199701e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 2.0261570199e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -3.9858700841e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -6.0674115194e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 8.3266726847e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 5.2280121423e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 1.0283692387e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 1.9428902931e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 5.1278749067e-17, 1e-10) );
  }

  // Release arrays
  gkyl_array_release(distf);
  gkyl_array_release(moms);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);
  gkyl_array_release(fpo_moms);
  gkyl_array_release(boundary_corrections);
  gkyl_array_release(drag_diff_coeff_corrs);
  gkyl_array_release(fpo_moms_ho);
  gkyl_array_release(boundary_corrections_ho);
  gkyl_array_release(drag_diff_coeff_corrs_ho);

  gkyl_proj_on_basis_release(proj_distf_square);
  gkyl_proj_on_basis_release(proj_drag_coeff);
  gkyl_proj_on_basis_release(proj_diff_coeff);
  gkyl_dg_updater_moment_release(mcalc);
  gkyl_mom_calc_release(fpo_mom_calc);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_fpo_vlasov_coeff_correct_release(coeff_correct_calc);
}

void test_1x3v_bump(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_bump = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_bump, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  struct gkyl_array *fpo_moms_ho, *boundary_corrections_ho, *drag_diff_coeff_corrs_ho;
  fpo_moms_ho = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);
  boundary_corrections_ho = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);
  drag_diff_coeff_corrs_ho = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume, false);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, use_gpu);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, use_gpu);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, use_gpu);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, use_gpu);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, use_gpu);

  // Project distribution function and coefficients and copy to device if using GPU
  if (use_gpu) {
    struct gkyl_array *distf_ho, *drag_coeff_ho, *diff_coeff_ho;
    distf_ho = mkarr(phase_basis.num_basis, phase_range_ext.volume, false);
    drag_coeff_ho = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume, false);
    diff_coeff_ho = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume, false);

    gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf_ho);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff_ho);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff_ho);

    gkyl_array_copy(distf, distf_ho);
    gkyl_array_copy(drag_coeff, drag_coeff_ho);
    gkyl_array_copy(diff_coeff, diff_coeff_ho);

    gkyl_array_release(distf_ho);
    gkyl_array_release(drag_coeff_ho);
    gkyl_array_release(diff_coeff_ho);
  }
  else {
    gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  }

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  else {
    gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, use_gpu);

  const double *fpo_moms_c, *bcorr_c, *coeff_corrs_c;
  if (use_gpu) {
    gkyl_array_copy(fpo_moms_ho, fpo_moms); 
    gkyl_array_copy(boundary_corrections_ho, boundary_corrections);
    gkyl_array_copy(drag_diff_coeff_corrs_ho, drag_diff_coeff_corrs);

    fpo_moms_c = gkyl_array_cfetch(fpo_moms_ho, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections_ho, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs_ho, 1);
  }
  else {
    fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);
  }

  if ((poly_order == 1) && (NV == 4)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1176082568e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -3.2040844934e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 9.4108748572e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], 2.6033758658e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 2.6020852140e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 4.4811268787e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.0529922609e+01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -5.1691319351e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 8.0716726476e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 5.4016554182e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 2.8528069664e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 3.0706526675e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -3.1035287187e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 1.2264291798e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.2286690590e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -3.8794949062e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.2727055669e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 1.4581989058e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 8.2819493451e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 7.7057824695e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 5.5063917835e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 1.1040832381e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -9.7549059286e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -8.9733795918e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -2.4692492062e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 4.3645535983e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -8.6714557454e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], -1.9740596199e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 6.7558826319e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -3.9655873472e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -7.0705797339e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -3.0865658754e-17, 1e-10) );
  }

  if ((poly_order == 1) && (NV == 8)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1442657996e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 2.5642265864e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 1.5666721392e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -6.9945533700e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.6032639626e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -6.5251993915e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.0434307818e+01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 6.8669677532e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.2895217299e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -5.2888124503e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -5.5353603103e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 6.6439114953e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -7.2590723580e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 4.7713078789e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.7158086920e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 1.6589162051e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.5290115501e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -1.1774855031e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 2.4733362060e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 3.0970016228e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 2.4733362060e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 3.0352862958e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -5.2148267350e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -1.2814567031e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -6.6464423797e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -2.4644981739e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.0428069876e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 7.2207936984e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 2.1880516673e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 6.4864002898e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -6.8762791995e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -7.9471729274e-17, 1e-10) );
  }

  if ((poly_order == 2) && (NV == 4)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1685375344e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -1.7479585815e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 4.7796917537e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -1.1115240672e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -3.1577874431e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 5.7255385992e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -6.3013830265e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 3.6390683829e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.7187609281e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -4.5431692243e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -3.3687946259e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -2.7985968577e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 1.4315267919e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], -5.2830870936e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 3.9437854024e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 1.0419175869e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], -7.2134147934e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.8875043712e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 8.2788993820e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 1.8473779258e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.6228429244e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 1.3764516147e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -8.5130270446e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -9.1181402706e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -1.2348316748e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 3.3139414810e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -1.1473213591e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 1.0264285073e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.2487628719e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -8.0844209384e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 6.5410810414e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -2.5874348491e-17, 1e-10) );
  }

  if ((poly_order == 2) && (NV == 8)) { 
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1500173245e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -4.2297996907e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], 4.8700768466e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -7.8116766528e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -2.8956004986e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -3.5526076250e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -8.2358707527e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 1.4196831814e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.9720462549e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 2.3204166111e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -8.9252038711e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -1.8465318250e-19, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -3.0004963329e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 2.2141983873e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 3.7269449679e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -9.5682677541e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.7471425182e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.9888185020e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], -8.7789158807e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[11], 3.4116279122e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.6589403275e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 8.5636755422e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[14], -1.9309248406e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -6.2612675461e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 2.4172048298e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 4.9019585586e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], -5.5783373589e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 7.2244614701e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 2.8893338044e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 4.0965708878e-18, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 7.7108655006e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -1.5502340858e-17, 1e-10) );
  }

  // Release arrays
  gkyl_array_release(distf);
  gkyl_array_release(moms);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);
  gkyl_array_release(fpo_moms);
  gkyl_array_release(boundary_corrections);
  gkyl_array_release(drag_diff_coeff_corrs);
  gkyl_array_release(fpo_moms_ho);
  gkyl_array_release(boundary_corrections_ho);
  gkyl_array_release(drag_diff_coeff_corrs_ho);

  gkyl_proj_on_basis_release(proj_distf_bump);
  gkyl_proj_on_basis_release(proj_drag_coeff);
  gkyl_proj_on_basis_release(proj_diff_coeff);
  gkyl_dg_updater_moment_release(mcalc);
  gkyl_mom_calc_release(fpo_mom_calc);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_fpo_vlasov_coeff_correct_release(coeff_correct_calc);
}

void test_1x3v_square_p1_4() { test_1x3v_square(1, 4, false); }
void test_1x3v_square_p1_8() { test_1x3v_square(1, 8, false); }
void test_1x3v_bump_p1_4() { test_1x3v_bump(1, 4, false); }
void test_1x3v_bump_p1_8() { test_1x3v_bump(1, 8, false); }

void test_1x3v_square_p2_4() { test_1x3v_square(2, 4, false); }
void test_1x3v_square_p2_8() { test_1x3v_square(2, 8, false); }
void test_1x3v_bump_p2_4() { test_1x3v_bump(2, 4, false); }
void test_1x3v_bump_p2_8() { test_1x3v_bump(2, 8, false); }

void test_1x3v_square_p1_4_cu() { test_1x3v_square(1, 4, true); }
void test_1x3v_square_p1_8_cu() { test_1x3v_square(1, 8, true); }
void test_1x3v_bump_p1_4_cu() { test_1x3v_bump(1, 4, true); }
void test_1x3v_bump_p1_8_cu() { test_1x3v_bump(1, 8, true); }

void test_1x3v_square_p2_4_cu() { test_1x3v_square(2, 4, true); }
void test_1x3v_square_p2_8_cu() { test_1x3v_square(2, 8, true); }
void test_1x3v_bump_p2_4_cu() { test_1x3v_bump(2, 4, true); }
void test_1x3v_bump_p2_8_cu() { test_1x3v_bump(2, 8, true); }

TEST_LIST = {
  { "test_1x3v_square_p1_4", test_1x3v_square_p1_4 },
  { "test_1x3v_square_p1_8", test_1x3v_square_p1_8 },
  { "test_1x3v_bump_p1_4", test_1x3v_bump_p1_4 },
  { "test_1x3v_bump_p1_8", test_1x3v_bump_p1_8 },
  { "test_1x3v_square_p2_4", test_1x3v_square_p2_4 },
  { "test_1x3v_square_p2_8", test_1x3v_square_p2_8 },
  { "test_1x3v_bump_p2_4", test_1x3v_bump_p2_4 },
  { "test_1x3v_bump_p2_8", test_1x3v_bump_p2_8 },

  #ifdef GKYL_HAVE_CUDA
  { "test_1x3v_square_p1_4_cu", test_1x3v_square_p1_4_cu },
  { "test_1x3v_square_p1_8_cu", test_1x3v_square_p1_8_cu },
  { "test_1x3v_bump_p1_4_cu", test_1x3v_bump_p1_4_cu },
  { "test_1x3v_bump_p1_8_cu", test_1x3v_bump_p1_8_cu },
  { "test_1x3v_square_p2_4_cu", test_1x3v_square_p2_4_cu },
  { "test_1x3v_square_p2_8_cu", test_1x3v_square_p2_8_cu },
  { "test_1x3v_bump_p2_4_cu", test_1x3v_bump_p2_4_cu },
  { "test_1x3v_bump_p2_8_cu", test_1x3v_bump_p2_8_cu },
  #endif
  { NULL, NULL }
};
