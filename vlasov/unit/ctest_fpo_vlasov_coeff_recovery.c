// Tests for the Vlasov FPO Maxwellian Rosenbluth potential projection and
// drag and diffusion coefficient computation

#include <math.h>

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
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

// Compare a drag or diffusion coefficient with its reference value. The coefficients are
// O(1) and the diffusion tensor takes two derivatives of the recovered potential, so
// round-off is ~1e-14 in absolute terms; a purely relative tolerance would be too strict
// for the small (1e-6) high-order coefficients, hence the absolute floor.
static int
compare_coeff(double val, double ref)
{
  return gkyl_compare_double(val, ref, 1e-10) || fabs(val - ref) < 1e-12;
}

// Allocate array
static struct gkyl_array *
mkarr(long nc, long size, bool use_gpu)
{
  if (use_gpu) {
    struct gkyl_array *a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
    return a;
  } else {
    struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
    return a;
  }
}

void
eval_gamma(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->gamma0;
}

void
eval_lte_moms(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->n0;
  fout[1] = app->ux0;
  fout[2] = app->uy0;
  fout[3] = app->uz0;
  fout[4] = app->vth0 * app->vth0;
}

struct fpo_ctx
create_ctx()
{
  struct fpo_ctx ctx = {.n0 = 1.0, .ux0 = 0.0, .uy0 = 0.0, .uz0 = 0.0, .vth0 = 1.0, .gamma0 = 1.0};
  return ctx;
}

void
test_1x3v(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim + vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {1, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 5.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  if (poly_order == 1) {
    // p=1 is the tensor hybrid basis: p=1 in configuration space, p=2 in velocity space.
    gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
    gkyl_cart_modal_tensor(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim - 1);
  } else {
    gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
    gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_serendip(&surf_basis, pdim - 1, poly_order);
  }

  gkyl_proj_on_basis *proj_gamma =
    gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order + 1, 1, eval_gamma, &ctx);
  gkyl_proj_on_basis *proj_lte =
    gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order + 1, 5, eval_lte_moms, &ctx);

  struct gkyl_array *lte_moms, *gamma, *h, *g, *h_surf, *g_surf, *dhdv_surf, *dgdv_surf,
    *d2gdv2_surf;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  lte_moms = mkarr(5 * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  gamma = mkarr(conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  h = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  g = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  h_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  g_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  dhdv_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  dgdv_surf = mkarr(2 * vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  d2gdv2_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2 * vdim * vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);

  struct gkyl_array *drag_coeff_ho, *diff_coeff_ho;
  drag_coeff_ho = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, false);
  diff_coeff_ho = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, false);

  // Array of relative offsets for recovery stencils
  long offsets[36];

  // Initialize updaters for potentials and coeffs
  gkyl_proj_maxwellian_pots_on_basis *pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
    &phase_grid, &conf_range, &phase_range, &conf_basis, &phase_basis, poly_order + 1, use_gpu
  );

  gkyl_fpo_vlasov_coeff_recovery *coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(
    &phase_grid, &phase_basis, &phase_range_ext, offsets, use_gpu
  );

  // Project moments and compute potentials
  if (use_gpu) {
    struct gkyl_array *lte_moms_ho, *gamma_ho;
    lte_moms_ho = mkarr(5 * conf_basis.num_basis, conf_range_ext.volume, false);
    gamma_ho = mkarr(conf_basis.num_basis, conf_range_ext.volume, false);
    gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma_ho);
    gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms_ho);

    gkyl_array_copy(gamma, gamma_ho);
    gkyl_array_copy(lte_moms, lte_moms_ho);

    gkyl_array_release(lte_moms_ho);
    gkyl_array_release(gamma_ho);
  } else {
    gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma);
    gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms);
  }

  gkyl_proj_maxwellian_pots_on_basis_advance(
    pot_slvr, &phase_range, &conf_range, lte_moms, h, g, h_surf, g_surf, dhdv_surf, dgdv_surf,
    d2gdv2_surf
  );

  // Compute drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(
    coeff_recovery, &phase_grid, phase_basis, &phase_range, &conf_range, gamma, h, dhdv_surf,
    drag_coeff, drag_coeff_surf, use_gpu
  );

  gkyl_calc_fpo_diff_coeff_recovery(
    coeff_recovery, &phase_grid, phase_basis, &phase_range, &conf_range, gamma, g, g_surf,
    dgdv_surf, d2gdv2_surf, diff_coeff, diff_coeff_surf, use_gpu
  );

  // Copy from device if using GPU
  if (use_gpu) {
    gkyl_array_copy(drag_coeff_ho, drag_coeff);
    gkyl_array_copy(diff_coeff_ho, diff_coeff);
  }

  // Write out drag and diffusion coefficients
  const char *fmt_drag = "ctest_fpo_drag_coeff_p%d_%d.gkyl";
  const char *fmt_diff = "ctest_fpo_diff_coeff_p%d_%d.gkyl";
  if (use_gpu) {
    int sz = gkyl_calc_strlen(fmt_drag, NV);
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof(fileNm), fmt_drag, poly_order, NV);
    gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, drag_coeff_ho, fileNm);
    snprintf(fileNm, sizeof(fileNm), fmt_diff, poly_order, NV);
    gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, diff_coeff_ho, fileNm);
  } else {
    int sz = gkyl_calc_strlen(fmt_drag, NV);
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof(fileNm), fmt_drag, poly_order, NV);
    gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, drag_coeff, fileNm);
    snprintf(fileNm, sizeof(fileNm), fmt_diff, poly_order, NV);
    gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, diff_coeff, fileNm);
  }

  // Compare values in corner cell and interior cell (of velocity space)
  // Checking first five components of a_x, a_y, D_xx, D_xy, and D_yx
  int idx_corner[] = {1, 1, 1, 1};
  int idx_int[] = {1, 2, 2, 2};
  long lin_corner = gkyl_range_idx(&phase_range, idx_corner);
  long lin_int = gkyl_range_idx(&phase_range, idx_int);

  double *drag_corner, *drag_int;
  double *diff_corner, *diff_int;
  if (use_gpu) {
    drag_corner = gkyl_array_fetch(drag_coeff_ho, lin_corner);
    drag_int = gkyl_array_fetch(drag_coeff_ho, lin_int);
    diff_corner = gkyl_array_fetch(diff_coeff_ho, lin_corner);
    diff_int = gkyl_array_fetch(diff_coeff_ho, lin_int);
  } else {
    drag_corner = gkyl_array_fetch(drag_coeff, lin_corner);
    drag_int = gkyl_array_fetch(drag_coeff, lin_int);
    diff_corner = gkyl_array_fetch(diff_coeff, lin_corner);
    diff_int = gkyl_array_fetch(diff_coeff, lin_int);
  }

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK(compare_coeff(drag_corner[0], 1.0957802248e-01));
    TEST_CHECK(compare_coeff(drag_corner[1], -4.9415567867e-17));
    TEST_CHECK(compare_coeff(drag_corner[2], 6.2671381068e-04));
    TEST_CHECK(compare_coeff(drag_corner[3], 2.1289992907e-02));
    TEST_CHECK(compare_coeff(drag_corner[4], 2.1289992907e-02));
    TEST_CHECK(compare_coeff(drag_corner[40], -2.5274433594e-18));
    TEST_CHECK(compare_coeff(drag_corner[41], -6.5106217797e-18));
    TEST_CHECK(compare_coeff(drag_corner[42], -2.1733654521e-18));
    TEST_CHECK(compare_coeff(drag_corner[43], -2.1535061821e-17));
    TEST_CHECK(compare_coeff(drag_corner[44], 5.4078567739e-05));

    TEST_CHECK(compare_coeff(drag_int[0], 5.9970334621e-01));
    TEST_CHECK(compare_coeff(drag_int[1], -3.4835699528e-17));
    TEST_CHECK(compare_coeff(drag_int[2], -1.7386662731e-01));
    TEST_CHECK(compare_coeff(drag_int[3], 1.8927389074e-01));
    TEST_CHECK(compare_coeff(drag_int[4], 1.8927389074e-01));
    TEST_CHECK(compare_coeff(drag_int[40], -2.7806857173e-17));
    TEST_CHECK(compare_coeff(drag_int[41], -1.0373791980e-17));
    TEST_CHECK(compare_coeff(drag_int[42], -6.1880002881e-18));
    TEST_CHECK(compare_coeff(drag_int[43], 6.2008886577e-19));
    TEST_CHECK(compare_coeff(drag_int[44], 5.6049908838e-05));

    TEST_CHECK(compare_coeff(diff_corner[0], 4.1009998893e-01));
    TEST_CHECK(compare_coeff(diff_corner[1], 8.5718139186e-16));
    TEST_CHECK(compare_coeff(diff_corner[2], 7.5275719022e-02));
    TEST_CHECK(compare_coeff(diff_corner[3], 2.7322422345e-03));
    TEST_CHECK(compare_coeff(diff_corner[4], 2.7322422345e-03));
    TEST_CHECK(compare_coeff(diff_corner[40], -1.8271977312e-16));
    TEST_CHECK(compare_coeff(diff_corner[41], -3.2156216911e-16));
    TEST_CHECK(compare_coeff(diff_corner[42], -6.3742957636e-16));
    TEST_CHECK(compare_coeff(diff_corner[43], -1.0314396180e-16));
    TEST_CHECK(compare_coeff(diff_corner[44], 9.6804526908e-05));
    TEST_CHECK(compare_coeff(diff_corner[80], -6.9332740679e-04));
    TEST_CHECK(compare_coeff(diff_corner[81], 8.1005109601e-16));
    TEST_CHECK(compare_coeff(diff_corner[82], -6.1891790027e-15));
    TEST_CHECK(compare_coeff(diff_corner[83], -1.9759325823e-15));
    TEST_CHECK(compare_coeff(diff_corner[84], 1.2465676662e-18));

    TEST_CHECK(compare_coeff(diff_int[0], 1.1361360918e+00));
    TEST_CHECK(compare_coeff(diff_int[1], 2.6118467600e-15));
    TEST_CHECK(compare_coeff(diff_int[2], 3.2279185522e-01));
    TEST_CHECK(compare_coeff(diff_int[3], 7.7070123529e-02));
    TEST_CHECK(compare_coeff(diff_int[4], 7.7070123529e-02));
    TEST_CHECK(compare_coeff(diff_int[40], -1.3684099487e-16));
    TEST_CHECK(compare_coeff(diff_int[41], 1.5441006644e-16));
    TEST_CHECK(compare_coeff(diff_int[42], -9.4337127140e-17));
    TEST_CHECK(compare_coeff(diff_int[43], 2.7825987072e-16));
    TEST_CHECK(compare_coeff(diff_int[44], 1.1400909128e-04));
    TEST_CHECK(compare_coeff(diff_int[80], -1.8500866926e-03));
    TEST_CHECK(compare_coeff(diff_int[81], 1.8955640785e-16));
    TEST_CHECK(compare_coeff(diff_int[82], 2.6236734572e-16));
    TEST_CHECK(compare_coeff(diff_int[83], 1.0337485959e-16));
    TEST_CHECK(compare_coeff(diff_int[84], -4.5232032317e-17));
  }

  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK(compare_coeff(drag_corner[0], 8.0437908204e-02));
    TEST_CHECK(compare_coeff(drag_corner[1], 1.6494310649e-16));
    TEST_CHECK(compare_coeff(drag_corner[2], 3.4791321315e-05));
    TEST_CHECK(compare_coeff(drag_corner[3], 6.6446451202e-03));
    TEST_CHECK(compare_coeff(drag_corner[4], 6.6446451202e-03));
    TEST_CHECK(compare_coeff(drag_corner[40], -3.1505816345e-18));
    TEST_CHECK(compare_coeff(drag_corner[41], 1.7127799509e-18));
    TEST_CHECK(compare_coeff(drag_corner[42], -3.1539068265e-18));
    TEST_CHECK(compare_coeff(drag_corner[43], 3.4832020723e-17));
    TEST_CHECK(compare_coeff(drag_corner[44], 4.8921997222e-07));

    TEST_CHECK(compare_coeff(drag_int[0], 1.5767743298e-01));
    TEST_CHECK(compare_coeff(drag_int[1], 1.3173841120e-16));
    TEST_CHECK(compare_coeff(drag_int[2], 1.7245670964e-04));
    TEST_CHECK(compare_coeff(drag_int[3], 1.8263418950e-02));
    TEST_CHECK(compare_coeff(drag_int[4], 1.8263418950e-02));
    TEST_CHECK(compare_coeff(drag_int[40], 9.0319273399e-18));
    TEST_CHECK(compare_coeff(drag_int[41], -4.6590896088e-18));
    TEST_CHECK(compare_coeff(drag_int[42], 8.3731123776e-18));
    TEST_CHECK(compare_coeff(drag_int[43], -1.7310625109e-17));
    TEST_CHECK(compare_coeff(drag_int[44], 5.1679262922e-06));

    TEST_CHECK(compare_coeff(diff_corner[0], 3.5183645246e-01));
    TEST_CHECK(compare_coeff(diff_corner[1], 8.6129172451e-15));
    TEST_CHECK(compare_coeff(diff_corner[2], 2.8049423283e-02));
    TEST_CHECK(compare_coeff(diff_corner[3], 5.6189623384e-04));
    TEST_CHECK(compare_coeff(diff_corner[4], 5.6189623385e-04));
    TEST_CHECK(compare_coeff(diff_corner[40], 5.7739416586e-16));
    TEST_CHECK(compare_coeff(diff_corner[41], 6.8560736815e-16));
    TEST_CHECK(compare_coeff(diff_corner[42], 4.9965168395e-17));
    TEST_CHECK(compare_coeff(diff_corner[43], 2.2994597736e-16));
    TEST_CHECK(compare_coeff(diff_corner[44], 1.1523233322e-06));
    TEST_CHECK(compare_coeff(diff_corner[80], -5.4094739355e-05));
    TEST_CHECK(compare_coeff(diff_corner[81], 1.7612760239e-15));
    TEST_CHECK(compare_coeff(diff_corner[82], -1.1817749284e-14));
    TEST_CHECK(compare_coeff(diff_corner[83], 1.1837281187e-14));
    TEST_CHECK(compare_coeff(diff_corner[84], 2.2845315175e-17));

    TEST_CHECK(compare_coeff(diff_int[0], 4.9271652157e-01));
    TEST_CHECK(compare_coeff(diff_int[1], -1.3150363908e-14));
    TEST_CHECK(compare_coeff(diff_int[2], 5.2728454299e-02));
    TEST_CHECK(compare_coeff(diff_int[3], 2.2095876647e-03));
    TEST_CHECK(compare_coeff(diff_int[4], 2.2095876647e-03));
    TEST_CHECK(compare_coeff(diff_int[40], 1.4651386003e-15));
    TEST_CHECK(compare_coeff(diff_int[41], -8.4127363233e-17));
    TEST_CHECK(compare_coeff(diff_int[42], -1.8993865843e-16));
    TEST_CHECK(compare_coeff(diff_int[43], -6.0384558258e-16));
    TEST_CHECK(compare_coeff(diff_int[44], 7.7354501762e-06));
    TEST_CHECK(compare_coeff(diff_int[80], -1.6365418417e-04));
    TEST_CHECK(compare_coeff(diff_int[81], -5.3378517424e-16));
    TEST_CHECK(compare_coeff(diff_int[82], -9.4709407825e-15));
    TEST_CHECK(compare_coeff(diff_int[83], 3.6437216418e-15));
    TEST_CHECK(compare_coeff(diff_int[84], 1.8168112347e-16));
  }

  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK(compare_coeff(drag_corner[0], 1.095780224753593e-01));
    TEST_CHECK(compare_coeff(drag_corner[1], 5.120107105299084e-18));
    TEST_CHECK(compare_coeff(drag_corner[2], 6.267138106830090e-04));
    TEST_CHECK(compare_coeff(drag_corner[3], 2.128999290678767e-02));
    TEST_CHECK(compare_coeff(drag_corner[4], 2.128999290678767e-02));
    TEST_CHECK(compare_coeff(drag_corner[48], 1.095780224753592e-01));
    TEST_CHECK(compare_coeff(drag_corner[49], 9.574526289061817e-18));
    TEST_CHECK(compare_coeff(drag_corner[50], 2.128999290678771e-02));
    TEST_CHECK(compare_coeff(drag_corner[51], 6.267138106832865e-04));
    TEST_CHECK(compare_coeff(drag_corner[52], 2.128999290678768e-02));

    TEST_CHECK(compare_coeff(drag_int[0], 5.997033462078909e-01));
    TEST_CHECK(compare_coeff(drag_int[1], 1.303726025098014e-16));
    TEST_CHECK(compare_coeff(drag_int[2], -1.738666273128696e-01));
    TEST_CHECK(compare_coeff(drag_int[3], 1.892738907386429e-01));
    TEST_CHECK(compare_coeff(drag_int[4], 1.892738907386429e-01));
    TEST_CHECK(compare_coeff(drag_int[48], 5.997033462078908e-01));
    TEST_CHECK(compare_coeff(drag_int[49], 3.189525375466083e-17));
    TEST_CHECK(compare_coeff(drag_int[50], 1.892738907386430e-01));
    TEST_CHECK(compare_coeff(drag_int[51], -1.738666273128699e-01));
    TEST_CHECK(compare_coeff(drag_int[52], 1.892738907386429e-01));

    TEST_CHECK(compare_coeff(diff_corner[0], 4.100999889255090e-01));
    TEST_CHECK(compare_coeff(diff_corner[1], -2.896597844153597e-15));
    TEST_CHECK(compare_coeff(diff_corner[2], 7.527571902224374e-02));
    TEST_CHECK(compare_coeff(diff_corner[3], 2.732242234506259e-03));
    TEST_CHECK(compare_coeff(diff_corner[4], 2.732242234505176e-03));
    TEST_CHECK(compare_coeff(diff_corner[48], -1.883902242384383e-01));
    TEST_CHECK(compare_coeff(diff_corner[49], 6.841182943492474e-17));
    TEST_CHECK(compare_coeff(diff_corner[50], 1.546275427632733e-03));
    TEST_CHECK(compare_coeff(diff_corner[51], -6.946834285633851e-03));
    TEST_CHECK(compare_coeff(diff_corner[52], -3.346156723907082e-02));
    TEST_CHECK(compare_coeff(diff_corner[96], -1.883902242384383e-01));
    TEST_CHECK(compare_coeff(diff_corner[97], 4.119538750209884e-17));
    TEST_CHECK(compare_coeff(diff_corner[98], 1.546275427632331e-03));
    TEST_CHECK(compare_coeff(diff_corner[99], -3.346156723907068e-02));
    TEST_CHECK(compare_coeff(diff_corner[100], -6.946834285633850e-03));

    TEST_CHECK(compare_coeff(diff_int[0], 1.136136091845040e+00));
    TEST_CHECK(compare_coeff(diff_int[1], 1.129764761964678e-16));
    TEST_CHECK(compare_coeff(diff_int[2], 3.227918552176493e-01));
    TEST_CHECK(compare_coeff(diff_int[3], 7.707012352914944e-02));
    TEST_CHECK(compare_coeff(diff_int[4], 7.707012352914944e-02));
    TEST_CHECK(compare_coeff(diff_int[48], -1.888164518268589e-01));
    TEST_CHECK(compare_coeff(diff_int[49], -7.715770094895800e-17));
    TEST_CHECK(compare_coeff(diff_int[50], 7.351883249525898e-02));
    TEST_CHECK(compare_coeff(diff_int[51], 7.351883249525724e-02));
    TEST_CHECK(compare_coeff(diff_int[52], -4.145240662121676e-02));
    TEST_CHECK(compare_coeff(diff_int[96], -1.888164518268586e-01));
    TEST_CHECK(compare_coeff(diff_int[97], -1.021219910658699e-16));
    TEST_CHECK(compare_coeff(diff_int[98], 7.351883249525877e-02));
    TEST_CHECK(compare_coeff(diff_int[99], -4.145240662121677e-02));
    TEST_CHECK(compare_coeff(diff_int[100], 7.351883249525794e-02));
  }

  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK(compare_coeff(drag_corner[0], 8.043790820400039e-02));
    TEST_CHECK(compare_coeff(drag_corner[1], 7.367084060495940e-17));
    TEST_CHECK(compare_coeff(drag_corner[2], 3.479132131522081e-05));
    TEST_CHECK(compare_coeff(drag_corner[3], 6.644645120171883e-03));
    TEST_CHECK(compare_coeff(drag_corner[4], 6.644645120171879e-03));
    TEST_CHECK(compare_coeff(drag_corner[48], 8.043790820400080e-02));
    TEST_CHECK(compare_coeff(drag_corner[49], 8.051495686328235e-17));
    TEST_CHECK(compare_coeff(drag_corner[50], 6.644645120171901e-03));
    TEST_CHECK(compare_coeff(drag_corner[51], 3.479132131554971e-05));
    TEST_CHECK(compare_coeff(drag_corner[52], 6.644645120171889e-03));

    TEST_CHECK(compare_coeff(drag_int[0], 1.576774329802842e-01));
    TEST_CHECK(compare_coeff(drag_int[1], 3.284827258687446e-17));
    TEST_CHECK(compare_coeff(drag_int[2], 1.724567096380162e-04));
    TEST_CHECK(compare_coeff(drag_int[3], 1.826341894981168e-02));
    TEST_CHECK(compare_coeff(drag_int[4], 1.826341894981168e-02));
    TEST_CHECK(compare_coeff(drag_int[48], 1.576774329802842e-01));
    TEST_CHECK(compare_coeff(drag_int[49], 7.624730483962438e-17));
    TEST_CHECK(compare_coeff(drag_int[50], 1.826341894981167e-02));
    TEST_CHECK(compare_coeff(drag_int[51], 1.724567096387267e-04));
    TEST_CHECK(compare_coeff(drag_int[52], 1.826341894981169e-02));

    TEST_CHECK(compare_coeff(diff_corner[0], 3.518364524562381e-01));
    TEST_CHECK(compare_coeff(diff_corner[1], -7.032330922844327e-16));
    TEST_CHECK(compare_coeff(diff_corner[2], 2.804942328325896e-02));
    TEST_CHECK(compare_coeff(diff_corner[3], 5.618962338477527e-04));
    TEST_CHECK(compare_coeff(diff_corner[4], 5.618962338465006e-04));
    TEST_CHECK(compare_coeff(diff_corner[48], -1.673341378406727e-01));
    TEST_CHECK(compare_coeff(diff_corner[49], 7.207343536595653e-16));
    TEST_CHECK(compare_coeff(diff_corner[50], 5.142107433236290e-04));
    TEST_CHECK(compare_coeff(diff_corner[51], -2.534916638362508e-03));
    TEST_CHECK(compare_coeff(diff_corner[52], -1.320569383053611e-02));
    TEST_CHECK(compare_coeff(diff_corner[96], -1.673341378406691e-01));
    TEST_CHECK(compare_coeff(diff_corner[97], 1.862004531478780e-16));
    TEST_CHECK(compare_coeff(diff_corner[98], 5.142107433227497e-04));
    TEST_CHECK(compare_coeff(diff_corner[99], -1.320569383053591e-02));
    TEST_CHECK(compare_coeff(diff_corner[100], -2.534916638351369e-03));

    TEST_CHECK(compare_coeff(diff_int[0], 4.927165215706736e-01));
    TEST_CHECK(compare_coeff(diff_int[1], 1.699450009781498e-14));
    TEST_CHECK(compare_coeff(diff_int[2], 5.272845429890367e-02));
    TEST_CHECK(compare_coeff(diff_int[3], 2.209587664694937e-03));
    TEST_CHECK(compare_coeff(diff_int[4], 2.209587664696642e-03));
    TEST_CHECK(compare_coeff(diff_int[48], -2.178235155393483e-01));
    TEST_CHECK(compare_coeff(diff_int[49], 1.541939283232677e-16));
    TEST_CHECK(compare_coeff(diff_int[50], 1.994932779555455e-03));
    TEST_CHECK(compare_coeff(diff_int[51], 1.994932779572272e-03));
    TEST_CHECK(compare_coeff(diff_int[52], -2.304571525855550e-02));
    TEST_CHECK(compare_coeff(diff_int[96], -2.178235155393505e-01));
    TEST_CHECK(compare_coeff(diff_int[97], 9.904047191736779e-17));
    TEST_CHECK(compare_coeff(diff_int[98], 1.994932779558042e-03));
    TEST_CHECK(compare_coeff(diff_int[99], -2.304571525855535e-02));
    TEST_CHECK(compare_coeff(diff_int[100], 1.994932779577841e-03));
  }

  // Release memory
  gkyl_array_release(lte_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(h);
  gkyl_array_release(g);
  gkyl_array_release(h_surf);
  gkyl_array_release(g_surf);
  gkyl_array_release(dhdv_surf);
  gkyl_array_release(dgdv_surf);
  gkyl_array_release(d2gdv2_surf);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);
  gkyl_array_release(drag_coeff_ho);
  gkyl_array_release(diff_coeff_ho);

  gkyl_proj_maxwellian_pots_on_basis_release(pot_slvr);
  gkyl_fpo_vlasov_coeff_recovery_release(coeff_recovery);
}

void
test_1x3v_p1_4()
{
  test_1x3v(1, 4, false);
}
void
test_1x3v_p1_8()
{
  test_1x3v(1, 8, false);
}
void
test_1x3v_p2_4()
{
  test_1x3v(2, 4, false);
}
void
test_1x3v_p2_8()
{
  test_1x3v(2, 8, false);
}
void
test_1x3v_p1_4_cu()
{
  test_1x3v(1, 4, true);
}
void
test_1x3v_p1_8_cu()
{
  test_1x3v(1, 8, true);
}
void
test_1x3v_p2_4_cu()
{
  test_1x3v(2, 4, true);
}
void
test_1x3v_p2_8_cu()
{
  test_1x3v(2, 8, true);
}

TEST_LIST = {
  {"test_1x3v_p1_4", test_1x3v_p1_4},
  {"test_1x3v_p1_8", test_1x3v_p1_8},
  {"test_1x3v_p2_4", test_1x3v_p2_4},
  {"test_1x3v_p2_8", test_1x3v_p2_8},
#ifdef GKYL_HAVE_CUDA
  {"test_1x3v_p1_4_cu", test_1x3v_p1_4_cu},
  {"test_1x3v_p1_8_cu", test_1x3v_p1_8_cu},
  {"test_1x3v_p2_4_cu", test_1x3v_p2_4_cu},
  {"test_1x3v_p2_8_cu", test_1x3v_p2_8_cu},
#endif
  {NULL, NULL}
};
