// DIII-D negative-triangularity inner-wall-limited 3x2v regression with a
// numerical Miller equilibrium and scaled fluid neutrals.
//
// The app generates a Miller GEQDSK with the shaping, axis, toroidal field and
// q profile of the analytic (mapc2p) 3x2v DIII-D input, and builds tokamak
// geometry, the LCFS split and the core twist shift from it. The computational
// coordinates are psi (poloidal flux / 2 pi), alpha (field-line label, set from
// a binormal length) and the normalized poloidal arc length in [-pi, pi]. The
// outer surface is the analytic input's; the LCFS is placed at 2/3 of the
// uniform psi grid, which moves the inner surface slightly, so Nx must be a
// multiple of 3. Profiles, the core sources and the neutral initial condition
// are evaluated in the analytic input's radial coordinate
// x = r - r_LCFS + x_inner and Miller angle.
//
// The fluid neutrals (vdim = 0) are not advanced; their density is rescaled
// each step so that ionization balances the recycled ion flux to the targets.

#include <math.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

struct gk_app_ctx {
  int cdim, vdim;
  // Miller geometry parameters.
  double a_shift, Zvert_axis, Rmajor_axis, R0, rminor_LCFS, x_inner, r0, B0, kappa, delta;
  double q0, Cy, Bref, x_LCFS;
  struct gkyl_gen_miller_eqdsk_inp miller;
  gkyl_gen_miller_eqdsk *equilibrium; // Flux table for profile callbacks; owned by main.
  // Plasma parameters.
  double me, qe, mi, qi, n0, Te0, Ti0, c_s, rho_s;
  double rec_frac, D0_gamma, nu_frac;
  // Neutral temperature: T_n = w*Tn_cold + (1-w)*Tn_ion_frac*T_i, with the
  // cold weight w decaying from the SOL targets over Tn_cold_width (in z).
  double Tn_cold, Tn_ion_frac, Tn_cold_width;
  // Adaptive core source parameters.
  bool adapt_energy_srcCORE, adapt_particle_srcCORE;
  double center_srcCORE[3], sigma_srcCORE[3];
  double energy_srcCORE, particle_srcCORE, floor_srcCORE;
  // Grid parameters (x in outboard major radius and y in binormal length
  // before the app converts them to psi and alpha).
  double Lx, Ly, Lz;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters.
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
};

static void
miller_kappa(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct gk_app_ctx *ctx = vctx;
  out[0] = ctx->kappa;
}

static void
miller_delta(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct gk_app_ctx *ctx = vctx;
  out[0] = ctx->delta;
}

static void
miller_a_shift(double t, const double *xc, double *GKYL_RESTRICT out, void *vctx)
{
  const struct gk_app_ctx *ctx = vctx;
  out[0] = ctx->a_shift;
}

static double
qsafety(double R, void *ctx)
{
  // Cubic fit to the DIII-D q profile of the analytic input, used for
  // R >= R_match, which covers the simulated domain. The fit turns negative
  // near the axis, so inside R_match q continues as the exponential with the
  // same value and slope; it only sets psi offsets inside the domain.
  const double a[] = {
    154.51071835546747, -921.8584472748003, 1842.1077075366113, -1231.619813170522
  };
  const double R_match = 2.0;
  const double Rq = fmax(R, R_match);
  const double q = ((a[0] * Rq + a[1]) * Rq + a[2]) * Rq + a[3];
  if (R >= R_match) {
    return q;
  }
  const double dqdR = (3.0 * a[0] * Rq + 2.0 * a[1]) * Rq + a[2];
  return q * exp(dqdR / q * (R - R_match));
}

/**
 * Convert the analytic input's radial coordinate x = r - r_LCFS + x_inner to psi.
 *
 * @param ctx Simulation context with an initialized flux table.
 * @param x Radial coordinate [m].
 * @return Poloidal flux divided by 2*pi [Wb/rad].
 */
static double
psi_at_x(const struct gk_app_ctx *ctx, double x)
{
  return gkyl_gen_miller_eqdsk_psi_at_r(ctx->equilibrium, x + ctx->rminor_LCFS - ctx->x_inner);
}

/**
 * Recover the analytic input's radial coordinate x = r - r_LCFS + x_inner from psi.
 *
 * @param ctx Simulation context with an initialized flux table.
 * @param psi Poloidal flux divided by 2*pi [Wb/rad].
 * @return Radial coordinate [m].
 */
static double
x_at_psi(const struct gk_app_ctx *ctx, double psi)
{
  return gkyl_gen_miller_eqdsk_r_at_psi(ctx->equilibrium, psi) - ctx->rminor_LCFS + ctx->x_inner;
}

// Miller angle at the normalized poloidal arc length z.
static double
theta_at_z(const struct gk_app_ctx *ctx, double z)
{
  return gkyl_gen_miller_eqdsk_theta_at_theta_tok(ctx->equilibrium, z);
}

struct source_profile_ctx {
  const struct gk_app_ctx *app;
  double center[3], sigma[3], floor;
};

// Core source Gaussian in the analytic radial coordinate and Miller angle,
// uniform in the binormal direction.
static void
source_profile(double t, const double *xn, double *fout, void *vctx)
{
  const struct source_profile_ctx *src = vctx;
  const double dx = (x_at_psi(src->app, xn[0]) - src->center[0]) / src->sigma[0];
  const double dtheta = (theta_at_z(src->app, xn[2]) - src->center[2]) / src->sigma[2];
  fout[0] = exp(-0.5 * (dx * dx + dtheta * dtheta)) + src->floor;
}

// Radial tanh profile, a*tanh(b*(x - c)) + d, fitted to experimental data.
static double
tanh_profile(double x, double a, double b, double c, double d)
{
  return a * tanh(b * (x - c)) + d;
}

static void
zero_func(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
density_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  const double x = x_at_psi(ctx, xn[0]);
  fout[0] = tanh_profile(x, -8.8e18, 34.0, 0.093, 1.1e19);
}

static void
temp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  const double eV = GKYL_ELEMENTARY_CHARGE;
  const double x = x_at_psi(ctx, xn[0]);
  fout[0] = tanh_profile(x, -240.0 * eV, 32.5, 0.08, 260.0 * eV);
}

static void
temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  const double eV = GKYL_ELEMENTARY_CHARGE;
  const double x = x_at_psi(ctx, xn[0]);
  fout[0] = tanh_profile(x, -1.6e3 * eV, 3.2, 0.16, 14.0 * eV);
}

// Neutral density peaked symmetrically at both targets and suppressed inside the LCFS.
static void
neut_density_init(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  const double x = x_at_psi(ctx, xn[0]), theta = theta_at_z(ctx, xn[2]);
  const double x_LCFS = ctx->x_inner; // LCFS position in the analytic radial coordinate.
  const double alpha_neut = 0.1 * ctx->n0;
  const double radial_factor = 1.0 - tanh((x_LCFS - x) / 0.05);
  // Symmetric in theta: the same peak at both targets.
  const double dist_to_target = 0.5 * ctx->Lz - fabs(theta);
  fout[0] = alpha_neut * radial_factor * (pow(1.0 / cosh(dist_to_target / 0.5), 2) + 1.0e-6);
}

static void
temp_neut(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  const double x = x_at_psi(ctx, xn[0]);
  fout[0] = 0.6 * ctx->Ti0 * (0.5 * (1.0 + tanh(2.0 - 25.0 * x)) + 0.1);
}

// Weight of the cold (recycled) neutral temperature: 1 at the SOL targets,
// decaying over Tn_cold_width along z, and 0 in the core.
static void
neut_cold_weight(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  const double psi = xn[0], z = xn[2];
  if (psi < ctx->x_LCFS) {
    fout[0] = 0.0;
    return;
  }
  const double w = exp(-(z - ctx->z_min) / ctx->Tn_cold_width) +
                   exp(-(ctx->z_max - z) / ctx->Tn_cold_width);
  fout[0] = fmin(w, 1.0);
}

static void
udrift(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = fout[1] = fout[2] = 0.0;
}

// Velocity maps: linear up to vpar_max/2, then quadratic; quadratic in mu.
static void
mapc2p_vel_elc(double t, const double *vc, double *GKYL_RESTRICT vp, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  const double cvpar = vc[0], cmu = vc[1];
  if (fabs(cvpar) <= 0.5) {
    vp[0] = ctx->vpar_max_elc * cvpar;
  } else {
    vp[0] = copysign(ctx->vpar_max_elc * 2.0 * cvpar * cvpar, cvpar);
  }
  vp[1] = ctx->mu_max_elc * cmu * cmu;
}

static void
mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  const double cvpar = vc[0], cmu = vc[1];
  if (fabs(cvpar) <= 0.5) {
    vp[0] = ctx->vpar_max_ion * cvpar;
  } else {
    vp[0] = copysign(ctx->vpar_max_ion * 2.0 * cvpar * cvpar, cvpar);
  }
  vp[1] = ctx->mu_max_ion * cmu * cmu;
}

/**
 * Build the Miller equilibrium and derive the simulation parameters from it.
 *
 * @param out Simulation context, which must remain alive while geometry callbacks run.
 * @return True when the equilibrium was created successfully.
 */
static bool
create_ctx(struct gk_app_ctx *out)
{
  int cdim = 3, vdim = 2;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV, qe = -eV;

  // Geometry and magnetic field.
  double Zvert_axis = -0.0014645315; // Magnetic axis height [m].
  double Rmajor_axis = 1.7074685; // Magnetic axis major radius [m].
  double bmag_axis = 2.0; // Toroidal field at the axis [T]; positive B_phi as in the analytic input.
  double R_LCFSmid = 2.17; // Major radius of the LCFS at the outboard midplane [m].
  double x_inner = 0.10; // Radial extent inside the LCFS [m].
  double x_outer = 0.05; // Radial extent outside the LCFS [m].
  double R0 = R_LCFSmid + 0.5 * (x_outer - x_inner); // Outboard reference radius [m].
  double a_shift = 1.0; // Shafranov shift parameter.
  double kappa = 1.35; // Elongation.
  double delta = -0.4; // Triangularity.
  double Ly_length_over_rho_s = 150.0; // Binormal domain length divided by rho_s.

  // Miller shape functions need the final context address.
  *out = (struct gk_app_ctx){
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .kappa = kappa,
    .delta = delta,
  };
  out->miller = (struct gkyl_gen_miller_eqdsk_inp){
    .input_coord_limits_type =
      {GKYL_GK_GEOMETRY_COORD_RADIAL_R_OMP, GKYL_GK_GEOMETRY_COORD_BINORMAL_LENGTH, 0},
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .bmag_axis = bmag_axis,
    .Rmajor_LCFS = R_LCFSmid,
    .kappa = miller_kappa,
    .kappa_ctx = out,
    .delta = miller_delta,
    .delta_ctx = out,
    .a_shift = miller_a_shift,
    .a_shift_ctx = out,
    .qsafety = qsafety,
    .has_Rmajor_omp_bounds = true,
  };
  double rminor_LCFS = gkyl_gen_miller_eqdsk_rminor_at_R_omp(&out->miller, R_LCFSmid);
  double outer_R = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&out->miller, rminor_LCFS + x_outer);
  out->miller.Rmajor_omp_bounds[0] = Rmajor_axis;
  out->miller.Rmajor_omp_bounds[1] = outer_R;
  out->equilibrium = gkyl_gen_miller_eqdsk_inew(&out->miller);
  if (!out->equilibrium) {
    return false;
  }
  const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(out->equilibrium);
  rminor_LCFS = info->rminor_LCFS;
  double B0 = gkyl_gen_miller_eqdsk_bmag_at_RZ(out->equilibrium, R0, Zvert_axis);
  if (!isfinite(B0) || B0 <= 0) {
    gkyl_gen_miller_eqdsk_release(out->equilibrium);
    out->equilibrium = NULL;
    return false;
  }

  // Plasma parameters.
  double mi = 2.01410177811 * mp; // Deuterium ions.
  double Te0 = 300.0 * eV;
  double Ti0 = 300.0 * eV;
  double n0 = 2.0e19; // [1/m^3]
  double vte = sqrt(Te0 / me), vti = sqrt(Ti0 / mi);
  double c_s = sqrt(Te0 / mi);
  double rho_s = c_s / fabs(qi * B0 / mi);

  // Radial domain: keep the analytic outer surface and put the LCFS at 2/3 of
  // the psi grid.
  double Lx = x_inner + x_outer;
  out->rminor_LCFS = rminor_LCFS;
  out->x_inner = x_inner;
  double x_LCFS = gkyl_round_to_e16_9(info->psi_boundary);
  double psi_upper = psi_at_x(out, Lx);
  double psi_lower = (x_LCFS - (2.0 / 3.0) * psi_upper) / (1.0 / 3.0);
  if (!(psi_lower > info->psi_axis)) {
    gkyl_gen_miller_eqdsk_release(out->equilibrium);
    out->equilibrium = NULL;
    return false;
  }
  double x_min = gkyl_gen_miller_eqdsk_R_omp_at_psi(out->equilibrium, psi_lower);

  // Binormal domain, quantized to an integer number of toroidal wedges.
  double psi_bounds[2] = {psi_lower, psi_upper};
  double requested_Ly = Ly_length_over_rho_s * rho_s;
  double y_bounds[2] = {-0.5 * requested_Ly, 0.5 * requested_Ly};
  struct gkyl_gen_miller_eqdsk_binormal_extent binormal;
  if (!gkyl_gen_miller_eqdsk_set_binormal_extents(
        out->equilibrium, psi_bounds, y_bounds, &binormal
      )) {
    gkyl_gen_miller_eqdsk_release(out->equilibrium);
    out->equilibrium = NULL;
    return false;
  }

  double Lz = 2.0 * M_PI - 1e-10; // Domain size along the field line.

  // Core source: injects 0.75 MW per species over the torus (scaled to the
  // simulated toroidal wedge) and replaces particles lost through the inner
  // radial boundary.
  double energy_srcCORE = 0.75e6 / binormal.n_tor;

  // Grid (reduced resolution for the regression test).
  int Nx = 6; // The LCFS is at 2/3 of the domain, so Nx must be a multiple of 3.
  int Ny = 4;
  int Nz = 8;
  int Nvpar = 6;
  int Nmu = 4;

  struct gk_app_ctx ctx_tmp = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .Zvert_axis = Zvert_axis,
    .Rmajor_axis = Rmajor_axis,
    .R0 = R0,
    .rminor_LCFS = rminor_LCFS,
    .x_inner = x_inner,
    .r0 = binormal.r0,
    .B0 = B0,
    .kappa = kappa,
    .delta = delta,
    .q0 = binormal.q0,
    .Cy = binormal.Cy,
    .Bref = B0,
    .x_LCFS = x_LCFS,
    .miller = out->miller,
    .equilibrium = out->equilibrium,
    .me = me,
    .qe = qe,
    .mi = mi,
    .qi = qi,
    .n0 = n0,
    .Te0 = Te0,
    .Ti0 = Ti0,
    .c_s = c_s,
    .rho_s = rho_s,
    .rec_frac = 0.98,
    .D0_gamma = 5.0 / 3.0,
    .nu_frac = 0.5,
    .Tn_cold = 3.0 * eV, // Franck-Condon recycled atoms.
    .Tn_ion_frac = 0.5,
    .Tn_cold_width = 0.1, // About a 3 eV atom's CX mean free path [rad].
    .adapt_energy_srcCORE = true,
    .adapt_particle_srcCORE = true,
    // Centered on the inner radial surface; the ion source is at +Lz/4 and the
    // electron source at -Lz/4 (Miller angle).
    .center_srcCORE = {x_inner + gkyl_gen_miller_eqdsk_r_at_psi(out->equilibrium, psi_lower) -
                         rminor_LCFS,
                       0.0, 0.25 * Lz},
    .sigma_srcCORE = {0.03 * Lx, 0.0, Lz / 6.0},
    .energy_srcCORE = energy_srcCORE,
    .particle_srcCORE = 0.0,
    .floor_srcCORE = 1e-10,
    .Lx = Lx,
    .Ly = binormal.Ly,
    .Lz = Lz,
    .x_min = x_min,
    .x_max = outer_R,
    .y_min = y_bounds[0],
    .y_max = y_bounds[1],
    .z_min = -0.5 * Lz,
    .z_max = 0.5 * Lz,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = 1,
    .vpar_max_elc = 5.0 * vte,
    .mu_max_elc = me * pow(4.0 * vte, 2) / (2.0 * B0),
    .vpar_max_ion = 5.0 * vti,
    .mu_max_ion = mi * pow(4.0 * vti, 2) / (2.0 * B0),
    .t_end = 1.0e-7,
    .num_frames = 1,
    .write_phase_freq = 1.0,
    .int_diag_calc_num = 10,
    .dt_failure_tol = 1.0e-3,
    .num_failures_max = 20,
  };
  *out = ctx_tmp;
  return true;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = {0}; // Stable context for Miller callbacks and profiles.
  bool initialized = create_ctx(&ctx);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  int failed = !initialized || cells_x[0] < 3 || cells_x[0] % 3 != 0;
  int any_failed;
  gkyl_comm_allreduce_host(comm, GKYL_INT, GKYL_MAX, 1, &failed, &any_failed);
  if (any_failed) {
    fprintf(stderr, "DIII-D Miller input requires a valid equilibrium and Nx divisible by 3.\n");
    gkyl_gen_miller_eqdsk_release(ctx.equilibrium);
    gkyl_gyrokinetic_comms_release(comm);
#ifdef GKYL_HAVE_MPI
    if (app_args.use_mpi) {
      MPI_Finalize();
    }
#endif
    return 1;
  }

  struct source_profile_ctx core_elc = {
    .app = &ctx,
    .center = {ctx.center_srcCORE[0], 0.0, -ctx.center_srcCORE[2]},
    .sigma = {ctx.sigma_srcCORE[0], 0.0, ctx.sigma_srcCORE[2]},
    .floor = ctx.floor_srcCORE,
  };
  struct source_profile_ctx core_ion = core_elc;
  core_ion.center[2] = ctx.center_srcCORE[2];

  // Electron core source.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_elc,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_min = 0.5 * ctx.Te0,
    .temp_max = 5.0 * ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e = {
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .mapc2p = {.mapping = mapc2p_vel_elc, .ctx = &ctx},
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = density_init,
        .ctx_density = &ctx,
        .upar = zero_func,
        .ctx_upar = &ctx,
        .temp = temp_elc,
        .ctx_temp = &ctx,
      },
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .num_cross_collisions = 1,
        .collide_with = {"ion"},
        .den_ref = ctx.n0,
        .temp_ref = ctx.Te0,
        .nu_frac = ctx.nu_frac,
      },
    .react_neut =
      {
        .num_react = 2,
        .write_diagnostics = true,
        .react_type =
          {{
             .react_id = GKYL_REACT_IZ,
             .type_self = GKYL_SELF_ELC,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .donor_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_RECOMB,
             .type_self = GKYL_SELF_ELC,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .recvr_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           }},
      },
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .num_adapt_sources = 1,
        .projection[0] = proj_srcCORE_e,
        .adapt[0] = adapt_srcCORE_e,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },
    .bcs =
      {
        {.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
        {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
        {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
        {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
      },
    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics =
      {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      },
  };

  // Ion core source.
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .configuration_profile = source_profile,
    .configuration_profile_ctx = &core_ion,
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0 * ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .mapc2p = {.mapping = mapc2p_vel_ion, .ctx = &ctx},
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = density_init,
        .ctx_density = &ctx,
        .upar = zero_func,
        .ctx_upar = &ctx,
        .temp = temp_ion,
        .ctx_temp = &ctx,
      },
    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},
    .collisions =
      {
        .collision_id = GKYL_LBO_COLLISIONS,
        .num_cross_collisions = 1,
        .collide_with = {"elc"},
        .den_ref = ctx.n0,
        .temp_ref = ctx.Ti0,
        .nu_frac = ctx.nu_frac,
      },
    .react_neut =
      {
        .num_react = 3,
        .write_diagnostics = true,
        .react_type =
          {{
             .react_id = GKYL_REACT_IZ,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .donor_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_RECOMB,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .elc_nm = "elc",
             .ion_nm = "ion",
             .recvr_nm = "D0",
             .charge_state = 0,
             .ion_mass = ctx.mi,
             .elc_mass = ctx.me,
           },
           {
             .react_id = GKYL_REACT_CX,
             .type_self = GKYL_SELF_ION,
             .ion_id = GKYL_ION_H,
             .ion_nm = "ion",
             .partner_nm = "D0",
             .ion_mass = ctx.mi,
             .partner_mass = ctx.mi,
           }},
      },
    .source =
      {
        .source_id = GKYL_PROJ_SOURCE,
        .num_sources = 1,
        .num_adapt_sources = 1,
        .projection[0] = proj_srcCORE_i,
        .adapt[0] = adapt_srcCORE_i,
        .diagnostics =
          {
            .num_diag_moments = 1,
            .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
            .num_integrated_diag_moments = 1,
            .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
          },
      },
    .bcs =
      {
        {.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
        {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
        {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
        {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING},
      },
    .num_diag_moments = 9,
    .diag_moments =
      {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
       GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR,
       GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics =
      {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      },
  };

  // Fluid neutrals, rescaled each step so that ionization balances recycling
  // of the ion flux to the targets.
  struct gkyl_gyrokinetic_neut_species D0 = {
    .name = "D0",
    .mass = ctx.mi,
    .vdim = 0,
    .gas_gamma = ctx.D0_gamma,
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = neut_density_init,
        .ctx_density = &ctx,
        .udrift = udrift,
        .ctx_udrift = &ctx,
        .temp = temp_neut,
        .ctx_temp = &ctx,
      },
    // Recycled atoms are cold near the targets and heat toward a fraction of
    // T_i by charge exchange away from them.
    .temp_model =
      {
        .temp_cold = ctx.Tn_cold,
        .ion_temp_frac = ctx.Tn_ion_frac,
        .ion_nm = "ion",
        .cold_weight = neut_cold_weight,
        .cold_weight_ctx = &ctx,
      },
    .scaling =
      {
        .type = GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE,
        .impacting_ion_name = "ion",
        .impacting_ion_id = GKYL_ION_H,
        .electron_name = "elc",
        .recycling_coeff = ctx.rec_frac,
        .num_boundaries = 2,
        .boundaries_dir = {2, 2},
        .boundaries_edge = {GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
      },
  };

  // Zero potential where the LCFS meets the targets.
  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {.perp_dirs = {0, 2}, .perp_coords = {ctx.x_LCFS, ctx.z_min}, .val = 0.0},
    {.perp_dirs = {0, 2}, .perp_coords = {ctx.x_LCFS, ctx.z_max}, .val = 0.0},
  };
  struct gkyl_poisson_bias_line_list bias_line_list = {.num_bias_line = 2, .bl = target_corner_bcs};

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx.Bref,
    .poisson_bcs =
      {
        {.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
        {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
      },
    .bias_line_list = &bias_line_list,
    .time_rate_diagnostics = true,
  };

  // The app writes <app-name>-miller.geqdsk and builds tokamak geometry from
  // it, including the LCFS split and the core twist shift.
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MILLER,
    .miller_grid_info = ctx.miller,
  };

  // Binormal angle alpha uses y = Cy*alpha; rho_s is the reference sound gyroradius.
  struct gkyl_msgpack_map_elem additional_metadata[] = {
    {.key = "Cy", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Cy},
    {.key = "rho_s", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.rho_s},
    {.key = "Rmajor_axis", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Rmajor_axis},
    {.key = "Te0", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Te0},
    {.key = "c_s", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.c_s},
  };

  struct gkyl_gk app_inp = {
    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 0.5,
    .cdim = ctx.cdim,
    .lower = {ctx.x_min, ctx.y_min, ctx.z_min},
    .upper = {ctx.x_max, ctx.y_max, ctx.z_max},
    .cells = {cells_x[0], cells_x[1], cells_x[2]},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = geometry,
    .num_periodic_dir = 1,
    .periodic_dirs = {1},
    .num_species = 2,
    .species = {elc, ion},
    .num_neut_species = 1,
    .neut_species = {D0},
    .field = field,
    .parallelism =
      {
        .use_gpu = app_args.use_gpu,
        .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
        .comm = comm,
      },
    .metadata =
      {
        .num_attributes = sizeof(additional_metadata) / sizeof(additional_metadata[0]),
        .attributes = additional_metadata,
      },
  };
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping =
      {
        .t_end = ctx.t_end,
        .num_frames = ctx.num_frames,
        .write_phase_freq = ctx.write_phase_freq,
        .int_diag_calc_num = ctx.int_diag_calc_num,
        .dt_failure_tol = ctx.dt_failure_tol,
        .num_failures_max = ctx.num_failures_max,
        .is_restart = app_args.is_restart,
        .restart_frame = app_args.restart_frame,
        .num_steps = app_args.num_steps,
      },
  };
  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gen_miller_eqdsk_release(ctx.equilibrium);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
