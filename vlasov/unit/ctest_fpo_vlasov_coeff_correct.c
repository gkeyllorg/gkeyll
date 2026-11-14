// Tests for the Vlasov FPO routines to compute moments and boundary corrections
// and resulting corrections to drag and diffusion coefficients.
//
// Included are tests for the p=1 tensor hybrid basis (p=1 in configuration space, p=2 in
// velocity space) and the p=2 Serendipity basis. The reference values are regression values:
// those of p=2 were generated on Stellar, those of p=1 with the CPU build that first adopted
// the tensor hybrid kernels (after checking conservation with the regression tests). The
// energy boundary corrections of the bump tests, and the coefficient corrections that depend
// on them, were regenerated with a CPU build for both bases after correcting the upper
// velocity bounds given to the boundary correction updater (they were set to the lower bounds).

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
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_vlasov_velocity_map.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

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

// Calculator of the M0, M1, M2 moments on the (uniform) velocity grid of a phase-space grid.
struct fpo_test_moms {
  struct gkyl_rect_grid vel_grid; // Velocity-space grid.
  struct gkyl_range vel_range, vel_range_ext; // Velocity-space ranges.
  struct gkyl_basis vel_basis; // Velocity-space basis.
  struct gkyl_vlasov_velocity_map *vel_map; // Identity velocity map.
  struct gkyl_array *hamil, *hamil_inv; // Velocity-space Hamiltonian, v^2/2.
  struct gkyl_mom_calc *calc; // Moment calculator.
};

static void
fpo_test_moms_init(
  struct fpo_test_moms *tm, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
)
{
  int cdim = conf_basis->ndim, vdim = phase_basis->ndim - cdim;
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < vdim; ++d) {
    lower[d] = phase_grid->lower[cdim + d];
    upper[d] = phase_grid->upper[cdim + d];
    cells[d] = phase_grid->cells[cdim + d];
  }
  gkyl_rect_grid_init(&tm->vel_grid, vdim, lower, upper, cells);
  gkyl_create_grid_ranges(&tm->vel_grid, ghost, &tm->vel_range_ext, &tm->vel_range);
  if (phase_basis->b_type == GKYL_BASIS_MODAL_HYBRID) {
    gkyl_cart_modal_tensor(&tm->vel_basis, vdim, 2);
  } else {
    gkyl_cart_modal_serendip(&tm->vel_basis, vdim, phase_basis->poly_order);
  }

  tm->hamil = mkarr(tm->vel_basis.num_basis, tm->vel_range.volume, use_gpu);
  tm->hamil_inv = mkarr(tm->vel_basis.num_basis, tm->vel_range.volume, use_gpu);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = {0};
  tm->vel_map = gkyl_vlasov_velocity_map_new(
    &tm->vel_grid, &tm->vel_range, &tm->vel_basis, inp_vmap, false, use_gpu
  );
  gkyl_dg_vlasov_calc_hamil(
    &tm->vel_grid, &tm->vel_basis, &tm->vel_range, GKYL_MODEL_DEFAULT, tm->vel_map, tm->hamil,
    tm->hamil_inv, use_gpu
  );

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = conf_basis,
    .phase_basis = phase_basis,
    .vel_range = &tm->vel_range,
    .hamil_range = &tm->vel_range,
    .hamil = tm->hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_DEFAULT),
    .mom_type = GKYL_F_MOMENT_M0M1M2,
    .use_gpu = use_gpu,
    .vel_map = tm->vel_map,
  };
  struct gkyl_mom_type *mom_type = gkyl_mom_vlasov_inew(&inp_mom);
  tm->calc = gkyl_mom_calc_new(phase_grid, mom_type, use_gpu);
  gkyl_mom_type_release(mom_type);
}

static void
fpo_test_moms_advance(
  const struct fpo_test_moms *tm, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *distf, struct gkyl_array *moms,
  bool use_gpu
)
{
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(tm->calc, phase_range, conf_range, distf, moms);
  } else {
    gkyl_mom_calc_advance(tm->calc, phase_range, conf_range, distf, moms);
  }
}

static void
fpo_test_moms_release(struct fpo_test_moms *tm)
{
  gkyl_mom_calc_release(tm->calc);
  gkyl_array_release(tm->hamil);
  gkyl_array_release(tm->hamil_inv);
  gkyl_vlasov_velocity_map_release(tm->vel_map);
}

static inline double
bump_maxwellian(
  double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA,
  double bUx, double bUy, double bUz, double bS, double bVt
)
{
  double v2 = (vx - ux) * (vx - ux) + (vy - uy) * (vy - uy) + (vz - uz) * (vz - uz);
  double bv2 = (vx - bUx) * (vx - bUx) + (vy - bUy) * (vy - bUy) + (vz - bUz) * (vz - bUz);
  return n / pow(sqrt(2 * M_PI * vt * vt), 3) * exp(-v2 / (2 * vt * vt)) +
         n / pow(sqrt(2 * M_PI * bVt * bVt), 3) * exp(-bv2 / (2 * bVt * bVt)) * (bA * bA) /
           (bv2 + bS * bS);
}

void
eval_distf_square(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];

  double width = 2.0; // corresponds to a final vth of 2/3
  if (vx > -width && vx < width && vy > -width && vy < width && vz > -width && vz < width) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_distf_bump(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];

  fout[0] = bump_maxwellian(
    1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vth0, sqrt(0.15), 4.0 * app->vth0, 0.0, 0.0, 0.14,
    3.0 * app->vth0
  );
}

void
eval_drag_coeff(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vy * vz;
  fout[1] = vx * vz;
  fout[2] = vy * vz;
}

void
eval_diff_coeff(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vx * vx;
  fout[1] = vx * vy;
  fout[2] = vx * vz;
  fout[3] = vy * vx;
  fout[4] = vy * vy;
  fout[5] = vz * vz;
  fout[6] = vz * vx;
  fout[7] = vz * vy;
  fout[8] = vz * vz;
}

struct fpo_ctx
create_ctx()
{
  struct fpo_ctx ctx = {.n0 = 1.0, .ux0 = 0.0, .uy0 = 0.0, .uz0 = 0.0, .vth0 = 1.0, .gamma0 = 10.0};
  return ctx;
}

void
test_1x3v_square(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim + vdim;

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

  int num_quad = poly_order + 2;
  gkyl_proj_on_basis *proj_distf_square =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, 1, eval_distf_square, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, vdim * vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5 * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2 * vdim * vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections =
    mkarr(2 * (vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  struct gkyl_array *fpo_moms_ho, *boundary_corrections_ho, *drag_diff_coeff_corrs_ho;
  fpo_moms_ho = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);
  boundary_corrections_ho =
    mkarr(2 * (vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);
  drag_diff_coeff_corrs_ho = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);

  // Initialize updater to compute Five Moments
  struct fpo_test_moms mcalc;
  fpo_test_moms_init(&mcalc, &phase_grid, &conf_basis, &phase_basis, use_gpu);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type *fpo_mom_type =
    gkyl_mom_fpo_vlasov_new(&conf_basis, &phase_basis, &phase_range, use_gpu);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {.a = drag_coeff, .D = diff_coeff};
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, use_gpu);

  // Initialize updater to compute boundary corrections
  double v_bounds[2 * GKYL_MAX_DIM];
  for (int d = 0; d < vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(
    &phase_grid, &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, use_gpu
  );

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc =
    gkyl_fpo_coeff_correct_new(&phase_grid, &conf_basis, &conf_range, use_gpu);

  // Project distribution function and coefficients and copy to device if using GPU
  if (use_gpu) {
    struct gkyl_array *distf_ho, *drag_coeff_ho, *diff_coeff_ho;
    distf_ho = mkarr(phase_basis.num_basis, phase_range_ext.volume, false);
    drag_coeff_ho = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, false);
    diff_coeff_ho = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, false);

    gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf_ho);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff_ho);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff_ho);

    gkyl_array_copy(distf, distf_ho);
    gkyl_array_copy(drag_coeff, drag_coeff_ho);
    gkyl_array_copy(diff_coeff, diff_coeff_ho);

    gkyl_array_release(distf_ho);
    gkyl_array_release(drag_coeff_ho);
    gkyl_array_release(diff_coeff_ho);
  } else {
    gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  }

  // Compute Five Moments
  fpo_test_moms_advance(&mcalc, &phase_range, &conf_range, distf, moms, use_gpu);

  // Compute moments, boundary corrections, and coefficient corrections
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  } else {
    gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(
    coeff_correct_calc, &conf_range, &phase_range, fpo_moms, boundary_corrections, moms,
    drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf, diff_coeff, diff_coeff_surf, use_gpu
  );

  const double *fpo_moms_c, *bcorr_c, *coeff_corrs_c;
  if (use_gpu) {
    gkyl_array_copy(fpo_moms_ho, fpo_moms);
    gkyl_array_copy(boundary_corrections_ho, boundary_corrections);
    gkyl_array_copy(drag_diff_coeff_corrs_ho, drag_diff_coeff_corrs);

    fpo_moms_c = gkyl_array_cfetch(fpo_moms_ho, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections_ho, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs_ho, 1);
  } else {
    fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);
  }

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], 3.9079850467e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], 6.4651657510e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 7.1054273576e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], 1.3564427193e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], -5.9045901882e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], 8.4475690126e+02, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], -1.0272395136e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], -1.2212453271e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], -2.0203642972e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -2.2204460493e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], -4.2388834978e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 1.3616730332e-34, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], 1.8451844338e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], -8.7995510548e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], 4.2066914300e-17, 1e-10));
  }

  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], -2.0428103653e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], -3.6213625561e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 1.1546319456e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -1.9084651579e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], 3.5527136788e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], -3.9018151609e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], 8.4475690126e+02, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], -7.8034380868e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 6.3837823916e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], 1.1316757988e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -3.6082248300e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 5.9639536183e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], -1.1102230246e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], 1.2193172378e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], -8.7995510548e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], -5.6808620362e-17, 1e-10));
  }

  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], -1.7763568394e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], -6.3528162176e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 1.7215963393e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -7.1054273576e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -1.5719302273e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], 3.6965525345e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], -1.4210854715e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], -3.3883142039e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 5.5511151231e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], 1.9852550680e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -5.3799885604e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 2.2204460493e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 4.9122819603e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], -1.1551726670e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], 4.4408920985e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], 1.0588481887e-16, 1e-10));
  }

  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], -6.4837024638e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], 1.2754784269e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 1.9415716862e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -2.6645352591e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -1.6729638855e-14, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], -3.2907815638e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], -6.2172489379e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], -1.6409199701e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 0.0000000000e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 2.0261570199e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], -3.9858700841e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -6.0674115194e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 8.3266726847e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 5.2280121423e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], 1.0283692387e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], 1.9428902931e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], 5.1278749067e-17, 1e-10));
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
  fpo_test_moms_release(&mcalc);
  gkyl_mom_calc_release(fpo_mom_calc);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_fpo_vlasov_coeff_correct_release(coeff_correct_calc);
}

void
test_1x3v_bump(int poly_order, int NV, bool use_gpu)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim + vdim;

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

  int num_quad = poly_order + 2;
  gkyl_proj_on_basis *proj_distf_bump =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, 1, eval_distf_bump, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff =
    gkyl_proj_on_basis_new(&phase_grid, &phase_basis, num_quad, vdim * vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  moms = mkarr(5 * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_coeff = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  drag_coeff_surf = mkarr(vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, use_gpu);
  diff_coeff_surf = mkarr(2 * vdim * vdim * surf_basis.num_basis, phase_range_ext.volume, use_gpu);
  fpo_moms = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  boundary_corrections =
    mkarr(2 * (vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);
  drag_diff_coeff_corrs = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, use_gpu);

  struct gkyl_array *fpo_moms_ho, *boundary_corrections_ho, *drag_diff_coeff_corrs_ho;
  fpo_moms_ho = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);
  boundary_corrections_ho =
    mkarr(2 * (vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);
  drag_diff_coeff_corrs_ho = mkarr((vdim + 1) * conf_basis.num_basis, conf_range_ext.volume, false);

  // Initialize updater to compute Five Moments
  struct fpo_test_moms mcalc;
  fpo_test_moms_init(&mcalc, &phase_grid, &conf_basis, &phase_basis, use_gpu);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type *fpo_mom_type =
    gkyl_mom_fpo_vlasov_new(&conf_basis, &phase_basis, &phase_range, use_gpu);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {.a = drag_coeff, .D = diff_coeff};
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, use_gpu);

  // Initialize updater to compute boundary corrections
  double v_bounds[2 * GKYL_MAX_DIM];
  for (int d = 0; d < vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(
    &phase_grid, &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, use_gpu
  );

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc =
    gkyl_fpo_coeff_correct_new(&phase_grid, &conf_basis, &conf_range, use_gpu);

  // Project distribution function and coefficients and copy to device if using GPU
  if (use_gpu) {
    struct gkyl_array *distf_ho, *drag_coeff_ho, *diff_coeff_ho;
    distf_ho = mkarr(phase_basis.num_basis, phase_range_ext.volume, false);
    drag_coeff_ho = mkarr(vdim * phase_basis.num_basis, phase_range_ext.volume, false);
    diff_coeff_ho = mkarr(vdim * vdim * phase_basis.num_basis, phase_range_ext.volume, false);

    gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf_ho);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff_ho);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff_ho);

    gkyl_array_copy(distf, distf_ho);
    gkyl_array_copy(drag_coeff, drag_coeff_ho);
    gkyl_array_copy(diff_coeff, diff_coeff_ho);

    gkyl_array_release(distf_ho);
    gkyl_array_release(drag_coeff_ho);
    gkyl_array_release(diff_coeff_ho);
  } else {
    gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf);
    gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
    gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  }

  // Compute Five Moments
  fpo_test_moms_advance(&mcalc, &phase_range, &conf_range, distf, moms, use_gpu);

  // Compute moments, boundary corrections, and coefficient corrections
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  } else {
    gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  }
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(
    coeff_correct_calc, &conf_range, &phase_range, fpo_moms, boundary_corrections, moms,
    drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf, diff_coeff, diff_coeff_surf, use_gpu
  );

  const double *fpo_moms_c, *bcorr_c, *coeff_corrs_c;
  if (use_gpu) {
    gkyl_array_copy(fpo_moms_ho, fpo_moms);
    gkyl_array_copy(boundary_corrections_ho, boundary_corrections);
    gkyl_array_copy(drag_diff_coeff_corrs_ho, drag_diff_coeff_corrs);

    fpo_moms_c = gkyl_array_cfetch(fpo_moms_ho, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections_ho, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs_ho, 1);
  } else {
    fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
    bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
    coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);
  }

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], 1.1791515235e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], -2.2051457391e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 8.4546085410e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], 1.6229760761e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -9.3891908137e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], 1.1734758077e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], 2.0494508020e+01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], -1.1314797668e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 9.6637375780e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], -2.4250250528e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 2.2632720351e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 2.2560327305e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 2.2022856629e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 3.3651310126e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 1.7138908839e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], -6.7008311594e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 1.6135052317e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], -5.8732794233e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 5.5402731014e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 3.0562771597e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], -1.7780915629e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], -5.1105779546e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 3.0624441269e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 1.9020439616e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 2.3093235528e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], 1.7515169410e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -7.9603892455e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], -1.4174423659e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 7.4584750752e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], -1.2331661761e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], -6.0163813831e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], 4.0478582324e-16, 1e-10));
  }

  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], 1.1610542468e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], 1.6168758246e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], -2.0873602326e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -6.2709164510e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -1.1890987327e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], 3.4231539280e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], 2.0442234527e+01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], 2.5829209294e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 1.0336666910e-02, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], -9.6296026857e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], 1.8126505071e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], 2.7912168714e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], -4.3622196784e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 1.3724039117e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 5.5150718490e-02, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 2.1863927610e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 1.7172613954e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], -1.4051844703e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], -1.1451885447e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 4.4614987849e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], -2.3920210430e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 1.7068533809e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], 1.0009841486e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], 3.2786821285e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 7.9162200171e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], -1.7023697821e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], 2.0525496613e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 9.3087033744e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 1.1590889755e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], -3.2958785618e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], -6.5190290172e+00, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], -5.9422194769e-17, 1e-10));
  }

  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], 1.1685375344e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], -1.7479585815e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 4.7796917537e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -1.1115240672e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -3.1577874431e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], 5.7255385992e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], -6.3013830265e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], 3.6390683829e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 9.7187609281e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], -4.5431692243e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], -3.3687946259e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], -2.7985968577e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], 1.4315267919e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], -5.2830870936e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 3.9437854024e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], 1.0419175869e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], -7.2134147934e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 1.9884276654e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], 8.2788993820e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 1.8473779258e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 1.6228429244e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 1.3764516147e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], -8.5130270446e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], -9.1181402706e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 5.0693737855e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], 3.3139414810e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -1.1473213591e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 1.0264285073e-15, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 3.2487628719e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], -8.0844209384e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], 6.5410810414e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], -2.5874348491e-17, 1e-10));
  }

  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[0], 1.1500173245e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[1], -4.2297996907e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[2], 4.8700768466e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[3], -7.8116766528e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[4], -2.8956004986e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[5], -3.5526076250e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[6], -8.2358707527e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(fpo_moms_c[7], 1.4196831814e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[0], 9.9720462549e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[1], 2.3204166111e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[2], -8.9252038711e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[3], -1.8465318250e-19, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[4], -3.0004963329e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[5], 2.2141983873e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[6], 3.7269449679e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[7], -9.5682677541e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[8], 1.7471425182e-20, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[9], 5.4683720399e-02, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[10], -8.7789158807e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[11], 3.4116279122e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[12], 1.6589403275e-01, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[13], 8.5636755422e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[14], -1.9309248406e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(bcorr_c[15], -6.2612675461e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[0], 4.8906611506e-03, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[1], 4.9019585586e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[2], -5.5783373589e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[3], 7.2244614701e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[4], 2.8893338044e-16, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[5], 4.0965708878e-18, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[6], 7.7108655006e-17, 1e-10));
    TEST_CHECK(gkyl_compare_double(coeff_corrs_c[7], -1.5502340858e-17, 1e-10));
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
  fpo_test_moms_release(&mcalc);
  gkyl_mom_calc_release(fpo_mom_calc);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_fpo_vlasov_coeff_correct_release(coeff_correct_calc);
}

void
test_1x3v_square_p1_4()
{
  test_1x3v_square(1, 4, false);
}
void
test_1x3v_square_p1_8()
{
  test_1x3v_square(1, 8, false);
}
void
test_1x3v_bump_p1_4()
{
  test_1x3v_bump(1, 4, false);
}
void
test_1x3v_bump_p1_8()
{
  test_1x3v_bump(1, 8, false);
}

void
test_1x3v_square_p2_4()
{
  test_1x3v_square(2, 4, false);
}
void
test_1x3v_square_p2_8()
{
  test_1x3v_square(2, 8, false);
}
void
test_1x3v_bump_p2_4()
{
  test_1x3v_bump(2, 4, false);
}
void
test_1x3v_bump_p2_8()
{
  test_1x3v_bump(2, 8, false);
}

void
test_1x3v_square_p1_4_cu()
{
  test_1x3v_square(1, 4, true);
}
void
test_1x3v_square_p1_8_cu()
{
  test_1x3v_square(1, 8, true);
}
void
test_1x3v_bump_p1_4_cu()
{
  test_1x3v_bump(1, 4, true);
}
void
test_1x3v_bump_p1_8_cu()
{
  test_1x3v_bump(1, 8, true);
}

void
test_1x3v_square_p2_4_cu()
{
  test_1x3v_square(2, 4, true);
}
void
test_1x3v_square_p2_8_cu()
{
  test_1x3v_square(2, 8, true);
}
void
test_1x3v_bump_p2_4_cu()
{
  test_1x3v_bump(2, 4, true);
}
void
test_1x3v_bump_p2_8_cu()
{
  test_1x3v_bump(2, 8, true);
}

TEST_LIST = {
  {"test_1x3v_square_p1_4", test_1x3v_square_p1_4},
  {"test_1x3v_square_p1_8", test_1x3v_square_p1_8},
  {"test_1x3v_bump_p1_4", test_1x3v_bump_p1_4},
  {"test_1x3v_bump_p1_8", test_1x3v_bump_p1_8},
  {"test_1x3v_square_p2_4", test_1x3v_square_p2_4},
  {"test_1x3v_square_p2_8", test_1x3v_square_p2_8},
  {"test_1x3v_bump_p2_4", test_1x3v_bump_p2_4},
  {"test_1x3v_bump_p2_8", test_1x3v_bump_p2_8},

#ifdef GKYL_HAVE_CUDA
  {"test_1x3v_square_p1_4_cu", test_1x3v_square_p1_4_cu},
  {"test_1x3v_square_p1_8_cu", test_1x3v_square_p1_8_cu},
  {"test_1x3v_bump_p1_4_cu", test_1x3v_bump_p1_4_cu},
  {"test_1x3v_bump_p1_8_cu", test_1x3v_bump_p1_8_cu},
  {"test_1x3v_square_p2_4_cu", test_1x3v_square_p2_4_cu},
  {"test_1x3v_square_p2_8_cu", test_1x3v_square_p2_8_cu},
  {"test_1x3v_bump_p2_4_cu", test_1x3v_bump_p2_4_cu},
  {"test_1x3v_bump_p2_8_cu", test_1x3v_bump_p2_8_cu},
#endif
  {NULL, NULL}
};
