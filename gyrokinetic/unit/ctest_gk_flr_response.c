// Single-mode validation of the Pade FLR chain (see DR #797) at the app level.
//
// Straight slab, uniform B, Dirichlet walls in x, periodic in z, and a single
// ion density mode
//   n_i = n0*(1 + A*sin(kx*x)),  n_e = n0,  kx = pi/Lx,
// which is an eigenfunction of every operator in the chain. The initial field
// solve is linear and diagonal in kx, so with b = (kx*rho_i)^2 the potential
// with FLR effects on and off satisfies
//   phi_flr/phi_noflr = (1 + b)/(1 + b/2),
// from the density gyroaverage (Gamma_1 = 1/(1+b/2)) and the FLR inversion of
// phi (A = 1 + b). Only ICs and the initial field solve are needed. We sweep
// kx*rho_i to check the rho_i->0 limit (ratio -> 1) and the large-b Pade
// saturation (ratio -> 2), and that the field-level implementations (local
// term and FEM operator) agree.
//
#include <acutest.h>

#include <glob.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

struct flr_ctx {
  int cdim, vdim;
  double me, qe, mi, qi;
  double n0, Te0, Ti0, B0;
  double rho_i; // Ion gyroradius, sets the FLR scale.
  double kx; // Perturbation wavenumber.
  double pert_amp; // Relative amplitude of the density perturbation.
  double Lx, Lz;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  int cells[GKYL_MAX_DIM];
};

static struct flr_ctx
create_ctx(double kx_rho)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double me = GKYL_ELECTRON_MASS, mi = GKYL_PROTON_MASS;
  double qe = -eV, qi = eV;

  double n0 = 1.0e19;
  double Te0 = 100.0*eV, Ti0 = 100.0*eV;
  double B0 = 1.0;

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi);
  double rho_i = sqrt(Ti0*mi)/(qi*B0);

  double kx = kx_rho/rho_i;
  double Lx = M_PI/kx;
  double Lz = 1.0;

  struct flr_ctx ctx = {
    .cdim = 2, .vdim = 2,
    .me = me, .qe = qe, .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0, .B0 = B0,
    .rho_i = rho_i,
    .kx = kx,
    .pert_amp = 1.0e-2,
    .Lx = Lx, .Lz = Lz,
    .vpar_max_elc = 5.0*vte,
    .mu_max_elc = me*pow(5.0*vte,2)/(2.0*B0),
    .vpar_max_ion = 5.0*vti,
    .mu_max_ion = mi*pow(5.0*vti,2)/(2.0*B0),
    .cells = {32, 4, 12, 8},
  };
  return ctx;
}

static void
density_elc(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->n0;
}

static void
density_ion(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  double x = xn[0];
  fout[0] = app->n0*(1.0 + app->pert_amp*sin(app->kx*x));
}

static void
temp_elc(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->Te0;
}

static void
temp_ion(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->Ti0;
}

static void
zero_func(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

static void
bfield_func(double t, const double *xc, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = app->B0;
}

// Run the initial field solve with the given field FLR options and write
// frame 0. Only the ions provide a reference gyroradius.
static void
run_case(struct flr_ctx *ctx, struct gkyl_comm *comm, bool use_gpu,
  struct gkyl_gyrokinetic_field_flr flr, const char *name)
{
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx->qe, .mass = ctx->me,
    .vdim = ctx->vdim,
    .lower = { -ctx->vpar_max_elc, 0.0},
    .upper = {  ctx->vpar_max_elc, ctx->mu_max_elc},
    .cells = { ctx->cells[2], ctx->cells[3] },
    .polarization_density = ctx->n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = ctx,
      .ctx_upar = ctx,
      .ctx_temp = ctx,
      .density = density_elc,
      .upar = zero_func,
      .temp = temp_elc,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx->qi, .mass = ctx->mi,
    .vdim = ctx->vdim,
    .lower = { -ctx->vpar_max_ion, 0.0},
    .upper = {  ctx->vpar_max_ion, ctx->mu_max_ion},
    .cells = { ctx->cells[2], ctx->cells[3] },
    .polarization_density = ctx->n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = ctx,
      .ctx_upar = ctx,
      .ctx_temp = ctx,
      .density = density_ion,
      .upar = zero_func,
      .temp = temp_ion,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},

    .flr = { .Tperp = ctx->Ti0 },
  };

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx->B0,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .flr = flr,
  };

  struct gkyl_gk gk = {
    .cdim = ctx->cdim,
    .lower = { 0.0, -ctx->Lz/2.0 },
    .upper = { ctx->Lx, ctx->Lz/2.0 },
    .cells = { ctx->cells[0], ctx->cells[1] },
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {0.0},
      .mapc2p = mapc2p,
      .c2p_ctx = ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .comm = comm,
      .cuts = {1, 1},
      .use_gpu = use_gpu,
    },
  };
  strcpy(gk.name, name);

  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  gkyl_gyrokinetic_app_write_field(app, 0.0, 0);
  gkyl_gyrokinetic_app_release(app);
}

// L2 norm of a DG field, sqrt(int f^2 dx).
static double
calc_l2(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_range range_ext,
  struct gkyl_basis basis, struct gkyl_array *f)
{
  struct gkyl_array *l2 = gkyl_array_new(GKYL_DOUBLE, 1, range_ext.volume);
  gkyl_dg_calc_l2_range(&basis, 0, l2, 0, f, range);
  gkyl_array_scale_range(l2, grid.cellVolume, &range);
  double l2red[1];
  gkyl_array_reduce_range(l2red, l2, GKYL_SUM, &range);
  gkyl_array_release(l2);
  return sqrt(l2red[0]);
}

static void
read_field0(const char *name, struct gkyl_range *local, struct gkyl_array *phi)
{
  char fname[128];
  snprintf(fname, sizeof(fname), "%s-field_0.gkyl", name);
  struct gkyl_rect_grid grid_read;
  gkyl_grid_sub_array_read(&grid_read, local, phi, fname);
}

static void
remove_outputs(const char *name)
{
  char pattern[128];
  snprintf(pattern, sizeof(pattern), "%s-*.gkyl", name);
  glob_t g;
  if (glob(pattern, 0, NULL, &g) == 0) {
    for (size_t i=0; i<g.gl_pathc; i++) remove(g.gl_pathv[i]);
  }
  globfree(&g);
}

static void
check_case(double kx_rho, bool use_gpu)
{
  struct flr_ctx ctx = create_ctx(kx_rho);
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(false, use_gpu, stderr);

  const char *name_off = "ctest_gk_flr_response_flroff";
  const char *name_on = "ctest_gk_flr_response_flron";
  const char *name_op = "ctest_gk_flr_response_flron_op";

  run_case(&ctx, comm, use_gpu, (struct gkyl_gyrokinetic_field_flr) { .type = GKYL_GK_FLR_NONE }, name_off);
  run_case(&ctx, comm, use_gpu, (struct gkyl_gyrokinetic_field_flr) { .type = GKYL_GK_FLR_PADE }, name_on);
  run_case(&ctx, comm, use_gpu, (struct gkyl_gyrokinetic_field_flr) { .type = GKYL_GK_FLR_PADE,
    .use_fem_operator = true }, name_op);

  // Read back the potentials.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ctx.cdim, (double[]) { 0.0, -ctx.Lz/2.0 },
    (double[]) { ctx.Lx, ctx.Lz/2.0 }, (int[]) { ctx.cells[0], ctx.cells[1] });
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ctx.cdim, 1);
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_array *phi_off = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *phi_on = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *phi_op = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  read_field0(name_off, &local, phi_off);
  read_field0(name_on, &local, phi_on);
  read_field0(name_op, &local, phi_op);

  // Measured amplitudes and expected factors.
  double l2_off = calc_l2(grid, local, local_ext, basis, phi_off);
  double l2_on = calc_l2(grid, local, local_ext, basis, phi_on);

  double b = pow(ctx.kx*ctx.rho_i, 2.0);
  double ratio_expected = (1.0 + b)/(1.0 + b/2.0);
  double ratio = l2_on/l2_off;

  // Analytic amplitude of the FLR-off solve: phi = qi*A*n0/(eps*kx^2)*sin(kx*x).
  double eps_pol = ctx.n0*(ctx.mi + ctx.me)/pow(ctx.B0, 2.0);
  double phi_amp = ctx.qi*ctx.pert_amp*ctx.n0/(eps_pol*pow(ctx.kx, 2.0));
  double l2_off_expected = phi_amp*sqrt(ctx.Lx*ctx.Lz/2.0);

  // Mode shape: phi_on must be the same mode scaled by the expected ratio.
  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_set(diff, 1.0, phi_on);
  gkyl_array_accumulate(diff, -ratio_expected, phi_off);
  double l2_shape = calc_l2(grid, local, local_ext, basis, diff);

  // The two field-level paths (local term vs FEM operator) must agree.
  gkyl_array_set(diff, 1.0, phi_on);
  gkyl_array_accumulate(diff, -1.0, phi_op);
  double err_op = calc_l2(grid, local, local_ext, basis, diff)/l2_on;

  double err_ratio = fabs(ratio/ratio_expected - 1.0);
  double err_abs = fabs(l2_off/l2_off_expected - 1.0);
  double err_shape = l2_shape/l2_on;

  TEST_CHECK(err_abs < 0.05);
  TEST_MSG("kx*rho_i=%g: |phi| FLR-off measured/expected = %e / %e", kx_rho, l2_off, l2_off_expected);
  TEST_CHECK(err_ratio < 0.02);
  TEST_MSG("kx*rho_i=%g: amplitude ratio measured/expected = %f / %f", kx_rho, ratio, ratio_expected);
  TEST_CHECK(err_shape < 0.02);
  TEST_MSG("kx*rho_i=%g: mode-shape error = %e", kx_rho, err_shape);
  TEST_CHECK(err_op < 0.02);
  TEST_MSG("kx*rho_i=%g: local-term vs FEM-operator difference = %e", kx_rho, err_op);

  gkyl_array_release(phi_off);
  gkyl_array_release(phi_on);
  gkyl_array_release(phi_op);
  gkyl_array_release(diff);
  remove_outputs(name_off);
  remove_outputs(name_on);
  remove_outputs(name_op);
  gkyl_gyrokinetic_comms_release(comm);
}

// Sweep kx*rho_i from the rho_i->0 limit (ratio -> 1) to the large-b Pade
// saturation (ratio -> 2).
static const double kx_rho_list[] = { 0.1, 1.0, 4.0 };

static void
test_flr_response(bool use_gpu)
{
  for (int i=0; i<sizeof(kx_rho_list)/sizeof(kx_rho_list[0]); i++)
    check_case(kx_rho_list[i], use_gpu);
}

void test_flr_response_cpu(void) { test_flr_response(false); }

#ifdef GKYL_HAVE_CUDA
void test_flr_response_gpu(void) { test_flr_response(true); }
#endif

TEST_LIST = {
  { "test_flr_response_cpu", test_flr_response_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_flr_response_gpu", test_flr_response_gpu },
#endif
  { NULL, NULL },
};
