#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// State of the pseudo orbit-averaged integrator.
enum gk_poa_state {
  GK_POA_NONE = 0, // Haven't started.
  GK_POA_OAP, // Orbit averaged phase.
  GK_POA_FDP, // Full dynamics phase.
  GK_POA_COMPLETED, // Finished simulation.
};

struct gk_poa_phase_params {
  enum gk_poa_state phase; // Type of phase.
  int num_frames; // Number of frames.
  double duration; // Duration.
  double alpha; // Factor multiplying collisionless terms.
  bool is_static_field; // Whether to evolve the field.
  bool is_positivity_enabled; // Whether positivity is enabled.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type; // Type of df/dt multipler.
};

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx {
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p;
  double beta;
  double tau;
  double Ti0;
  double nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  double vti;
  double RatZeq0; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_max;
  double psi_min;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double mu_max_ion;
  int Npsi;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  // Source parameters
  double ion_source_amplitude;
  double ion_source_sigma;
  double ion_source_temp;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  int num_phases; // Number of phases.
  struct gk_poa_phase_params *poa_phases; // Phases to run.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Geometry parameters for Lorentzian mirror
  double mcB;     // Magnetic field parameter
  double gamma;   // Width parameter for Lorentzian profile
  double Z_m;     // Mirror throat location
  double Z_min;   // Minimum Z coordinate
  double Z_max;   // Maximum Z coordinate
  double psi_in;  // Working variable for psi integration
  double z_in;    // Working variable for z integration
};

double
psi_RZ(double RIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double psi = 0.5 * pow(RIn, 2.) * mcB *
    (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
    1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))));
  return psi;
}

double
R_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rout = sqrt(2. * psiIn / (mcB *
    (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
    1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))))));
  return Rout;
}

void
Bfield_psiZ(double psiIn, double ZIn, void *ctx, double *BRad, double *BZ, double *Bmag)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rcoord = R_psiZ(psiIn, ZIn, ctx);

  BRad[0] = -(1. / 2.) * Rcoord * mcB *
    (-2. * (ZIn - Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn - Z_m) / gamma, 2.), 2.))) -
    2. * (ZIn + Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn + Z_m) / gamma, 2.), 2.))));

  BZ[0] = mcB *
    (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
    1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))) );

  Bmag[0] = sqrt(pow(BRad[0], 2) + pow(BZ[0], 2));
}

double
integrand_z_psiZ(double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_in;
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, ZIn, ctx, &BRad, &BZ, &Bmag);
  return Bmag / BZ;
}

double
z_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double eps = 0.0;
  app->psi_in = psiIn;
  struct gkyl_qr_res integral;
  if (eps <= ZIn) {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, eps, ZIn, 7, 1e-14);
  }
  else {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, ZIn, eps, 7, 1e-14);
    integral.res = -integral.res;
  }
  return integral.res;
}

// Invert z(Z) via root-finding.
double
root_Z_psiz(double Z, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  return app->z_in - z_psiZ(app->psi_in, Z, ctx);
}

double
Z_psiz(double psiIn, double zIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double maxL = app->Z_max - app->Z_min;
  double eps = maxL / app->Nz;   // Interestingly using a smaller eps yields larger errors in some geo quantities.
  app->psi_in = psiIn;
  app->z_in = zIn;
  struct gkyl_qr_res Zout;
  if (0.0 <= zIn) {
    double fl = root_Z_psiz(-eps, ctx);
    double fr = root_Z_psiz(app->Z_max + eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, -eps, app->Z_max + eps, fl, fr, 1000, 1e-14);
  }
  else {
    double fl = root_Z_psiz(app->Z_min - eps, ctx);
    double fr = root_Z_psiz(eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, app->Z_min - eps, eps, fl, fr, 1000, 1e-14);
  }
  return Zout.res;
}

// Geometry evaluation functions for the gk app
void
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  double psi = xc[0], theta = xc[1], z = xc[2];

  double Z = Z_psiz(psi, z, ctx);
  double R = R_psiZ(psi, Z, ctx);

  // Cartesian coordinates on plane perpendicular to Z axis.
  double x = R * cos(theta);
  double y = R * sin(theta);

  xp[0] = x; xp[1] = y; xp[2] = Z;
}

void
bfield_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xc[2];
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx);
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);

  double phi = xc[1];
  // zc are computational coords.
  // Set Cartesian components of magnetic field.
  fout[0] = BRad * cos(phi);
  fout[1] = BRad * sin(phi);
  fout[2] = BZ;
}

// Evaluate collision frequencies
void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[1];
  fout[0] = 1e17 * exp(-2 * pow(fabs(z), 2));
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = 0.0;
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->Ti0;
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout,
  void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[1];
  double src_amp = app->ion_source_amplitude;
  double z_src = 0.0;
  double src_sigma = app->ion_source_sigma;
  double src_amp_floor = src_amp * 1e-2;
  if (fabs(z) <= 0.98) {
    fout[0] = src_amp * (1 - pow(fabs(z), 6) / 0.98);
  }
  else {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout,
  void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout,
  void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[1];
  double TSrc0 = app->ion_source_temp;
  double Tfloor = TSrc0 * 1e-2;
  if (fabs(z) <= 0.98) {
    fout[0] = TSrc0;
  }
  else {
    fout[0] = Tfloor;
  }
}

void mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.4;
  vp[0] = vpar_max_ion * tan(cvpar * b) / tan(b);
  vp[1] = mu_max_ion * pow(cmu, 3);
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.
  int poly_order = 1;

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS;
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;

  // Ion-ion collision freq.
  double nuFrac = 1.0;
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
    (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);

  // Grid parameters
  double vpar_max_ion = 16 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nz = 64;
  int Npsi = 4;
  int Nvpar = 32; // 96 uniform
  int Nmu = 16;  // 192 uniform

  // Source parameters
  double ion_source_amplitude = 1.e20;
  double ion_source_sigma = 0.5;
  double ion_source_temp = 5000. * eV;

  // Geometry parameters.
  double RatZeq0 = 0.10; // Radius of the field line at Z=0.
  double Z_min = -2.5;
  double Z_max = 2.5;
  double mcB = 3.691260;
  double gamma = 0.226381;
  double Z_m = 0.98;

  // POA parameters
  double alpha_oap = 5e-6;  // Factor multiplying collisionless terms.
  double alpha_fdp = 1.0;
  double tau_oap = 0.001;  // Duration of each phase.
  double tau_fdp = 7e-9;
  double tau_fdp_extra = 2e-9;
  int num_cycles = 2; // Number of OAP+FDP cycles to run.

  // Frame counts for each phase type (specified independently)
  int num_frames_oap = 1;        // Frames per OAP phase
  int num_frames_fdp = 1;        // Frames per FDP phase
  int num_frames_fdp_extra = 2;  // Frames for the extra FDP phase

  // Whether to evolve the field.
  bool is_static_field_oap = true;
  bool is_static_field_fdp = false;

  // Whether positivity is enabled.
  bool is_positivity_enabled_oap = false;
  bool is_positivity_enabled_fdp = false;

  // Type of df/dt multipler.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_oap = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE;
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_fdp = GKYL_GK_FDOT_MULTIPLIER_NONE;

  // Calculate phase structure
  double t_end = (tau_oap + tau_fdp) * num_cycles + tau_fdp_extra;
  double tau_pair = tau_oap + tau_fdp; // Duration of an OAP+FDP pair.
  int num_phases = 2 * num_cycles + 1;
  int num_frames = num_cycles * (num_frames_oap + num_frames_fdp) + num_frames_fdp_extra;

  struct gk_poa_phase_params *poa_phases = gkyl_calloc(num_phases,
    sizeof(struct gk_poa_phase_params));
  for (int i = 0; i < (num_phases - 1) / 2; i++) {
    // OAPs.
    poa_phases[2 * i].phase = GK_POA_OAP;
    poa_phases[2 * i].num_frames = num_frames_oap;
    poa_phases[2 * i].duration = tau_oap;
    poa_phases[2 * i].alpha = alpha_oap;
    poa_phases[2 * i].is_static_field = is_static_field_oap;
    poa_phases[2 * i].fdot_mult_type = fdot_mult_type_oap;
    poa_phases[2 * i].is_positivity_enabled = is_positivity_enabled_oap;

    // FDPs.
    poa_phases[2 * i + 1].phase = GK_POA_FDP;
    poa_phases[2 * i + 1].num_frames = num_frames_fdp;
    poa_phases[2 * i + 1].duration = tau_fdp;
    poa_phases[2 * i + 1].alpha = alpha_fdp;
    poa_phases[2 * i + 1].is_static_field = is_static_field_fdp;
    poa_phases[2 * i + 1].fdot_mult_type = fdot_mult_type_fdp;
    poa_phases[2 * i + 1].is_positivity_enabled = is_positivity_enabled_fdp;
  }
  // The final stage is an extra, longer FDP.
  poa_phases[num_phases - 1].phase = GK_POA_FDP;
  poa_phases[num_phases - 1].num_frames = num_frames_fdp_extra;
  poa_phases[num_phases - 1].duration = tau_fdp_extra;
  poa_phases[num_phases - 1].alpha = alpha_fdp;
  poa_phases[num_phases - 1].is_static_field = is_static_field_fdp;
  poa_phases[num_phases - 1].fdot_mult_type = fdot_mult_type_fdp;
  poa_phases[num_phases - 1].is_positivity_enabled = is_positivity_enabled_fdp;

  double write_phase_freq = 1; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq = 5; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te0 = Te0,
    .n0 = n0,
    .B_p = B_p,
    .beta = beta,
    .tau = tau,
    .Ti0 = Ti0,
    .nuFrac = nuFrac,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .RatZeq0 = RatZeq0,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .Npsi = Npsi,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = { Npsi, Nz, Nvpar, Nmu },
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .num_phases = num_phases,
    .poa_phases = poa_phases,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_freq = int_diag_calc_freq,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,

    .ion_source_amplitude = ion_source_amplitude,
    .ion_source_sigma = ion_source_sigma,
    .ion_source_temp = ion_source_temp,

    .mcB = mcB,
    .gamma = gamma,
    .Z_m = Z_m,
    .Z_min = Z_min,
    .Z_max = Z_max,
  };

  // Populate a couple more values in the context.
  ctx.psi_max = psi_RZ(ctx.RatZeq0, 0., &ctx);
  ctx.psi_min = psi_RZ(ctx.RatZeq0 / 10, 0., &ctx);
  ctx.z_min = z_psiZ(ctx.psi_max, ctx.Z_min, &ctx);
  ctx.z_max = z_psiZ(ctx.psi_max, ctx.Z_max, &ctx);

  return ctx;
}

void
release_ctx(struct gk_mirror_ctx *ctx)
{
  gkyl_free(ctx->poa_phases);
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app,
  double t_curr, bool force_calc, double dt)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if (!(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger *iot_conf, struct gkyl_tm_trigger *iot_phase,
  gkyl_gyrokinetic_app *app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr - 1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
    gkyl_gyrokinetic_app_write_dt(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr - 1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

struct time_frame_state {
  double t_curr; // Current simulation time.
  double t_end; // End time of current phase.
  int frame_curr; // Current frame.
  int num_frames; // Number of frames at the end of current phase.
};

void reset_io_triggers(struct gk_mirror_ctx *ctx, struct time_frame_state *tfs,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag)
{
  // Reset I/O triggers:
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  int frame_curr = tfs->frame_curr;
  int num_frames = tfs->num_frames;
  int num_int_diag_calc = ctx->int_diag_calc_freq * num_frames;

  // Prevent division by zero when frame_curr equals num_frames
  int frames_remaining = num_frames - frame_curr;
  double time_remaining = t_end - t_curr;

  trig_write_conf->dt = time_remaining / frames_remaining;
  trig_write_conf->tcurr = t_curr;
  trig_write_conf->curr = frame_curr;

  trig_write_phase->dt = time_remaining / (ctx->write_phase_freq * frames_remaining);
  trig_write_phase->tcurr = t_curr;
  trig_write_phase->curr = frame_curr;

  int diag_frames = GKYL_MAX2(frames_remaining,
    (num_int_diag_calc / num_frames) * frames_remaining);
  trig_calc_intdiag->dt = time_remaining / diag_frames;
  trig_calc_intdiag->tcurr = t_curr;
  trig_calc_intdiag->curr = frame_curr;
}

void run_phase(gkyl_gyrokinetic_app *app, struct gk_mirror_ctx *ctx, double num_steps,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag, struct time_frame_state *tfs,
  struct gk_poa_phase_params *pparams)
{
  tfs->t_end = tfs->t_curr + pparams->duration;
  tfs->num_frames = tfs->frame_curr + pparams->num_frames;

  // Run an OAP or FDP.
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;

  // Reset I/O triggers:
  reset_io_triggers(ctx, tfs, trig_write_conf, trig_write_phase, trig_calc_intdiag);

  // Reset simulation parameters and function pointers.
  struct gkyl_gyrokinetic_collisionless collisionless_inp = {
    .type = GKYL_GK_COLLISIONLESS_ES,
    .scale_factor = pparams->alpha,
  };
  struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp = {
    .num_multipliers = 1,
    .multiplier[0] = {
      .type = pparams->fdot_mult_type,
      .cellwise_const = true,
      .write_diagnostics = true,
    }
  };
  struct gkyl_gyrokinetic_field field_inp = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx->me,
    .electron_charge = ctx->qe,
    .electron_temp = ctx->Te0,
    .polarization_bmag = ctx->B_p,
    .is_static = pparams->is_static_field,
  };
  struct gkyl_gyrokinetic_positivity positivity_inp = {
    .type = pparams->is_positivity_enabled? GKYL_GK_POSITIVITY_SHIFT : GKYL_GK_POSITIVITY_NONE,
    .write_diagnostics = pparams->is_positivity_enabled,
  };

  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "ion", collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "ion", fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "ion", positivity_inp);
  gkyl_gyrokinetic_app_reset_field(app, t_curr, field_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    dt = t_end - t_curr; // Ensure we don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, t_curr >= t_end, status.dt_actual);
    write_data(trig_write_conf, trig_write_phase, app, t_curr, t_curr >= t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ",
          dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n",
          num_failures_max);
        calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, true, status.dt_actual);
        write_data(trig_write_conf, trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  tfs->t_curr = t_curr;
  tfs->frame_curr = tfs->frame_curr + pparams->num_frames;
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mirror_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .density = eval_density_ion,
      .ctx_density = &ctx,
      .upar = eval_upar_ion,
      .ctx_upar = &ctx,
      .temp = eval_temp_ion,
      .ctx_temp = &ctx,
    },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
      .scale_factor = 1.0, // Will be replaced below.
      .write_diagnostics = true,
    },

    .time_rate_multiplier = {
      .num_multipliers = 1,
      .multiplier[0] = {
        .type = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE,
        .cellwise_const = true,
        .write_diagnostics = true,
      }
    },

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te0,
      .write_diagnostics = true,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar = eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,
      },
      .diagnostics = {
        .num_diag_moments = 6,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR,
                          GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = { GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
                      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP,
                      GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
    },
  };
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .is_static = false,
  };

  struct gkyl_gk app_inp = {  // GK app
    .name = "gk_mirror_boltz_elc_poa_2x2v_p1",
    .cdim = ctx.cdim,
    .upper = { ctx.psi_max, ctx.Z_max },
    .lower = { ctx.psi_min, ctx.Z_min },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = { 0.0 },
      .mapc2p = mapc2p, // Mapping of computational to physical space.
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // Magnetic field.
      .bfield_ctx = &ctx,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = { ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  
  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Triggers for IO.
  struct gkyl_tm_trigger trig_write_conf, trig_write_phase, trig_calc_intdiag;

  struct time_frame_state tfs = {
    .t_curr = 0.0, // Initial simulation time.
    .frame_curr = 0, // Initial frame.
    .t_end = ctx.poa_phases[0].duration, // Final time of 1st phase.
    .num_frames = ctx.poa_phases[0].num_frames, // Number of frames in 1st phase.
  };

  int phase_idx_init = 0, phase_idx_end = ctx.num_phases; // Initial and final phase index.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app,
      app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    tfs.frame_curr = status.frame;
    tfs.t_curr = status.stime;

    // Find out what phase we are in.
    double time_count = 0.0;
    int frame_count = 0;
    int pit_curr = 0;
    for (int pit = 0; pit < ctx.num_phases; pit++) {
      time_count += ctx.poa_phases[pit].duration;
      frame_count += ctx.poa_phases[pit].num_frames;
      if ((tfs.t_curr <= time_count) && (tfs.frame_curr <= frame_count)) {
        pit_curr = pit;
        break;
      }
    }
    ;
    phase_idx_init = pit_curr;

    // Change the duration and number frames so this phase reaches the expected
    // time and number of frames and not beyond.
    struct gk_poa_phase_params *pparams = &ctx.poa_phases[phase_idx_init];
    pparams->num_frames = frame_count - tfs.frame_curr;
    pparams->duration = time_count - tfs.t_curr;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", tfs.frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", tfs.t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, tfs.t_curr);

    // Write out ICs.
    reset_io_triggers(&ctx, &tfs, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag);

    calc_integrated_diagnostics(&trig_calc_intdiag, app, tfs.t_curr, true, -1.0);
    write_data(&trig_write_conf, &trig_write_phase, app, tfs.t_curr, true);
  }

  if (app_args.num_steps != INT_MAX)
    phase_idx_end = 1;

  // Loop over number of number of phases;
  for (int pit = phase_idx_init; pit < phase_idx_end; pit++) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\nRunning phase %d @ t = %.9e ... \n", pit, tfs.t_curr);
    struct gk_poa_phase_params *phase_params = &ctx.poa_phases[pit];
    run_phase(app, &ctx, app_args.num_steps, &trig_write_conf, &trig_write_phase,
      &trig_calc_intdiag, &tfs, phase_params);
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n",
      stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n",
      stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

freeresources:
  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);
  release_ctx(&ctx);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  return 0;
}
