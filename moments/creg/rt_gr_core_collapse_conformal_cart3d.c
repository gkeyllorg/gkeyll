// Mixed-scheme coupled TOV star in Cartesian coordinates: GR Euler (gamma-law) fluid evolved with the MP flux-form scheme, vacuum Einstein (Bona-Masso) spacetime evolved with the
//  wave-propagation scheme, coupled each step.
//  Coupling once per step:
//  - spacetime -> fluid: the live metric is synced into the fluid metric slots at the start of the step (eqns 52/53 from the Multifluids paper);
//  - fluid -> spacetime: the fluid stress-energy adds the standard ADM matter sources
//    -8*pi*alpha*(T_ij - 1/2 gamma_ij T) to K_ij and +8*pi*alpha^2*T^0_k to V_k, folded into the spacetime's Strang-split half-step source integration.
// The fluid still has the frozen-discrete well-balancing (toggle with gr_euler_disable_well_balanced).
// The mixed scheme always applies the atmosphere floor in the primitive recovery (fluxes and sources only).
// Both species share the same analytic TOV spacetime object for their initial data only.

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_moment.h>
#include <gkyl_moment_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_gr_tov.h>
#include <gkyl_gr_einstein_excision_priv.h>
#include <gkyl_gr_blackhole.h>
#include "tov_solver.h"

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>
#include "rt_kdrop.h"

// Compile-time diagnostic comparison only; no runtime limiter override.
#ifndef GKYL_COLLAPSE_ST_LIMITER
#define GKYL_COLLAPSE_ST_LIMITER GKYL_MONOTONIZED_CENTERED
#endif
#ifndef GKYL_COLLAPSE_INITIAL_KS
#define GKYL_COLLAPSE_INITIAL_KS 0
#endif

// Run defaults for THIS collapse problem, not global solver defaults. Existing
// environment overrides still win. Physical/grid parameters are in create_ctx.
static void
collapse_run_defaults(void)
{
  const char *defaults[][2] = {
    { "GKYL_WB", "1" }, // Hold the fluid until the scheduled trigger.
    { "GKYL_DYN_COLLAPSE", "1" }, // Set to 0 for hold-only; default trigger time is 1.
    { "GKYL_KDROP_ETK_EULERIAN", "1" }, // Preserve Eulerian E and all S_i at the drop.
    { "GKYL_EXPERIMENTAL", "1" }, // Enables the separately selectable Gamma driver.
    { "GKYL_GAMMA_DRIVER", "0.75" }, // Integrated, nonadvective reference-Gamma coefficient.
    { "GKYL_GAMMA_DRIVER_ETA", "0.05" }, // Shift damping rate.
    { "GKYL_GAMMA_DRIVER_ALPHA_FLOOR", "0" }, // No positive-lapse cutoff.
    { "GKYL_PCP", "1" }, // Admissible reconstructed fluid faces.
    { "GKYL_SPACETIME_CELL_SPD_REPAIR", "1" }, // Existing cell metric protection.
    { "GKYL_GR_EULER_SURFACE_U1_OFF", "1" }, // Do not replace surface reconstruction with U1.
    { "GKYL_PHYS_DIAG", "1" }, // Density, pressure, velocity and enclosed rest mass.
    { "GKYL_DIAG_TRACK_EVERY", "4" }, // Physical diagnostics cadence in accepted steps.
    { "GKYL_CONSTRAINT_DIAG_EVERY", "4" }, // Constraint diagnostics cadence.
    { "GKYL_CONSTRAINT_DIAG_FINAL", "1" }, // Include the final constraint measurement.
    { "GKYL_CONSTRAINT_DELTA", "1" }, // Report change from initial discrete constraints.
    { "GKYL_AH_CERT", "1" }, // Diagnostic only; never sets an excision mask.
    { "GKYL_AH_CERT_EVERY", "8" }, // Horizon scan cadence.
    { "GKYL_AH_CERT_RMIN", "0.5" }, // Smallest sampled coordinate radius.
    { "GKYL_AH_CERT_RMAX", "9.5" }, // Largest sampled coordinate radius.
    { "GKYL_AH_CERT_NR", "80" }, // Number of sampled spheres.
    { "GKYL_AH_CERT_NTHETA", "25" }, // Polar quadrature samples.
    { "GKYL_AH_CERT_NPHI", "40" }, // Azimuthal quadrature samples.
    { "GKYL_AH_CERT_MAX_BAD_FRAC", "0" }, // Require valid geometry at every sampled angle.
  };
  for (size_t i=0; i<sizeof defaults/sizeof defaults[0]; ++i)
    if (setenv(defaults[i][0], defaults[i][1], 0) != 0) {
      perror(defaults[i][0]);
      exit(EXIT_FAILURE);
    }
}

static bool
collapse_flag_enabled(const char *name)
{
  const char *v = getenv(name);
  return v && strcmp(v, "0") && strcmp(v, "false") && strcmp(v, "False")
    && strcmp(v, "off") && strcmp(v, "OFF");
}

struct gr_tov_coupled_ctx
{
  double gas_gamma; // Adiabatic index (evolution EOS).
  double K_poly; // Polytropic constant of the cold-IC EOS (p = K rho^gamma).
  double rho_c; // Central rest-mass density.
  double dr_tov; // TOV solver radial step.
  double rho_atm; // Atmosphere rest-mass density.
  double p_atm; // Atmosphere pressure floor.

  struct gkyl_tov *tov; // Frozen TOV table.
  struct gkyl_gr_spacetime *spacetime; // TOV spacetime object (shared for initial data only).
  struct gkyl_wv_eqn *fluid_eqn; // GR-Euler eqn handle (set in main; used by the mid-run WB toggle).
  double M_star;
  double R_star;

  enum gkyl_spacetime_gauge spacetime_gauge; // Fluid-side gauge tag (GR Euler eqn).
  int reinit_freq;

  // Vacuum Einstein evolution parameters.
  double excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing;   // Lapse slicing throughout the run; unchanged by the trigger.
  enum gkyl_spacetime_evolution spacetime_evolution;

  double perturb_frac; // x-momentum kick fraction in the inner star (0 => static test).

  // Collapse trigger (core-localized K-drop, same physics as rt_gr_core_collapse_cart.c).
  double R_core;      // core radius (r < R_core is softened / watched).
  double kdrop_delta; // core pressure-softening fraction: p -> (1-delta) p for r < R_core (0 => untriggered).
  bool kdrop_smooth;  // smooth raised-cosine taper (else sharp top-hat).

  int Nx, Ny, Nz;
  double Lx, Ly, Lz; // Domain [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [-Lz/2,Lz/2], star centered at origin.
  double cfl_frac;

  double t_end;
  int num_frames;
  int field_energy_calcs;
  int integrated_mom_calcs;
  double dt_failure_tol;
  int num_failures_max;
};

struct gr_tov_coupled_ctx
create_ctx(void)
{
  // Neutron-star polytrope. This is intentionally ordinary driver setup, not an environment override layer:
  // choose paper scenarios by editing these local parameters or by moving them into the eventual par file.
  // Stable reference: Gamma=2, K=100, rho_c=1.28e-3. Turning point (M_max=1.637) at rho_c=3.18e-3
  // (stage1_diagnostics/tov_sequence). rho_c=1.0e-2 is past it (unstable, spontaneous collapse).
  double gas_gamma = 2.0;
  double K_poly = 100.0;
  double rho_c = 3.15e-3; // Near the mass maximum: WB hold -> E/S-preserving K-drop.
  // Diagnostic/run-control override while selecting the Stage-1 collapse model.
  // Final production runs should move this into the simulation parameter file.
  const char *rho_c_env = getenv("GKYL_RHO_C");
  if (rho_c_env != NULL) {
    double v = atof(rho_c_env);
    if (v > 0.0) rho_c = v;
  }
  double dr_tov = 0.01;

  double rho_atm = 1e-7 * rho_c;
  double p_atm = K_poly * pow(rho_atm, gas_gamma);

  double L = 24.0;
  const char *L_env = getenv("GKYL_LBOX");
  if (L_env != NULL) {
    L = atof(L_env);
  }

  double t_end = 200.0; // Formation pilot; longer-time collapse is not yet certified.
  const char *t_end_env = getenv("GKYL_T_END");
  if (t_end_env != NULL) {
    t_end = atof(t_end_env);
  }

  int num_frames = 200; // Standard .gkyl frames, written in the working directory.
  const char *num_frames_env = getenv("GKYL_NUM_FRAMES");
  if (num_frames_env != NULL) {
    num_frames = atoi(num_frames_env);
  }

  // GKYL_CFL overrides the CFL fraction (default 0.8). Lowering it (e.g. 0.2) shrinks dt so temporal
  // error is sub-dominant -- used to cleanly separate spatial from temporal error in a convergence study.
  double cfl_frac = 0.8;
  const char *cfl_env = getenv("GKYL_CFL");
  if (cfl_env != NULL) {
    cfl_frac = atof(cfl_env);
  }

  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);
  double M_star = gkyl_tov_star_mass(tov);
  double R_star = gkyl_tov_star_radius(tov);
  printf("M_star = %e \n", M_star);
  printf("R_star = %e \n", R_star);
  printf("Compactness (2M_star / R_star) = %e \n", 2.0 * M_star / R_star);

  // Collapse trigger: core-localized K-drop (soften core pressure -> loss of support -> collapse).
  // GKYL_KDROP=delta (0..1), GKYL_RCORE=fraction of R_star (default 0.4), GKYL_KDROP_SMOOTH=1 for a
  // raised-cosine taper. delta=0 (default) => untriggered (static) star.
  double rcore_frac = 1.0; // Deplete the whole initial star, not a smaller core.
  const char *rcore_env = getenv("GKYL_RCORE");
  if (rcore_env) rcore_frac = atof(rcore_env);
  double R_core = rcore_frac * R_star;
  double kdrop_delta = 0.02; // Fractional depletion of live effective K at the event.
  const char *kdrop_env = getenv("GKYL_KDROP");
  if (kdrop_env) kdrop_delta = atof(kdrop_env);
  bool kdrop_smooth = (getenv("GKYL_KDROP_SMOOTH") != NULL);
  printf("K-drop delta = %e (R_core = %.3f = %.2f R_star, smooth=%d)\n",
    kdrop_delta, R_core, rcore_frac, (int) kdrop_smooth);
  printf("IC EOS: gamma_ini=%g K_ini=%g rho_c=%g \n", gas_gamma, K_poly, rho_c);

  // The stationary reference is regular-centre areal TOV. Its metric/shift then
  // evolve freely; this does not enforce areal coordinates during collapse.
  bool use_kerr_schild = GKYL_COLLAPSE_INITIAL_KS; // Compile-time initial-slice control; not an evolution-gauge switch.

  // GKYL_SLICING = harmonic | 1pluslog | geodesic. Default: 1+log.
  enum gkyl_spacetime_slicing slicing = GKYL_1PLUSLOG_SLICING;
  const char *slc = getenv("GKYL_SLICING");
  if (slc) {
    if (strcmp(slc, "harmonic") == 0) slicing = GKYL_HARMONIC_SLICING;
    else if (strcmp(slc, "geodesic") == 0) slicing = GKYL_GEODESIC_SLICING;
    else if (strcmp(slc, "1pluslog") == 0) slicing = GKYL_1PLUSLOG_SLICING;
  }
  // Both reductions have matching matter sources. This checkpoint uses Einstein.
  enum gkyl_spacetime_evolution evolution = GKYL_EINSTEIN_EVOLUTION;
  const char *evo = getenv("GKYL_EVOLUTION");
  if (evo) {
    if (strcmp(evo, "einstein") == 0) evolution = GKYL_EINSTEIN_EVOLUTION;
    else if (strcmp(evo, "ricci") == 0) evolution = GKYL_RICCI_EVOLUTION;
  }
  // Lapse-threshold excision. NOTE this driver previously IGNORED GKYL_EXCISION_THRESHOLD, so every stage-1 run
  // excised at alpha < 0.3 while the job script set 0.001. Lapse-based excision is in tension with singularity-
  // avoiding slicing (1+log drives alpha down exactly where curvature grows), so the value matters: honour the env.
  // A small lapse is a gauge condition, not proof of trapping. Legacy lapse
  // excision remains opt-in; the default keeps the star until a horizon is validated.
  double excision_threshold = 0.0;
  { const char *e = getenv("GKYL_EXCISION_THRESHOLD"); if (e) { double v = atof(e); if (v >= 0.0) excision_threshold = v; } }

  printf("Gauge: use_kerr_schild=%d, slicing=%d, evolution=%d, excision_threshold=%g\n",
    (int) use_kerr_schild, (int) slicing, (int) evolution, excision_threshold);
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, 0.0, 0.0, 0.0, use_kerr_schild);

  struct gr_tov_coupled_ctx ctx = {
    .gas_gamma = gas_gamma,
    .K_poly = K_poly,
    .rho_c = rho_c,
    .dr_tov = dr_tov,
    .rho_atm = rho_atm,
    .p_atm = p_atm,
    .tov = tov,
    .spacetime = spacetime,
    .M_star = M_star,
    .R_star = R_star,
    .spacetime_gauge = GKYL_STATIC_GAUGE,
    .reinit_freq = 100,
    .excision_threshold = excision_threshold,
    .spacetime_slicing = slicing, // env-selectable (default 1+log); see gauge sweep above
    // The n=1 geometric constraint addition is paired with +4*pi*alpha*gamma_ij*E.
    .spacetime_evolution = evolution,
    .perturb_frac = 0.0, // set > 0 (e.g. 1.0e-3) for the momentum-perturbation test.
    .R_core = R_core,
    .kdrop_delta = kdrop_delta,
    .kdrop_smooth = kdrop_smooth,
    .Nx = 32,
    .Ny = 32,
    .Nz = 32,
    .Lx = L,
    .Ly = L,
    .Lz = L,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = INT_MAX,
    .integrated_mom_calcs = INT_MAX,
    .dt_failure_tol = 1.0e-4,
    .num_failures_max = 20,
  };

  return ctx;
}

void
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct gr_tov_coupled_ctx *app = ctx;

  double gas_gamma = app->gas_gamma;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double r = sqrt((x * x) + (y * y) + (z * z));
  struct tov_eval_bl bl = { 0 };
  gkyl_tov_eval_bl(app->tov, r, &bl);

  double p = fmax(bl.P, app->p_atm);
  double rho = pow(p / app->K_poly, 1.0 / gas_gamma); // non-hybrid atmosphere remains exactly K rho^Gamma
  double gamma_h = gas_gamma;
  double depletion_velocity_scale = 1.0;

  // Core-localized K-drop collapse trigger: soften the pressure in r < R_core while KEEPING the density
  // profile, so the core is under-pressured and collapses while the envelope stays near equilibrium.
  // delta=0 => static (bit-identical to the untriggered star).
  // When GKYL_DYN_COLLAPSE is set, the depletion is applied MID-RUN (dyn_collapse_trigger), not at t=0,
  // so the initial data is the clean (held) star. Skip the t=0 depletion in that case.
  if (app->kdrop_delta > 0.0 && r < app->R_core && rho > app->rho_atm
      && p > app->p_atm && getenv("GKYL_DYN_COLLAPSE") == NULL) {
    double local_delta = app->kdrop_delta;
    if (app->kdrop_smooth) {
      local_delta = app->kdrop_delta * 0.5 * (1.0 + cos(M_PI * r / app->R_core));
    }
    double Kp_etk = (1.0 - local_delta) * app->K_poly;
    double G_etk = gas_gamma;
    if (getenv("GKYL_KDROP_ETK_EULERIAN") != NULL) {
      // Preserve E AND S_i of the unperturbed analytic state. The changed
      // pressure requires changed velocity even though the metric stays fixed.
      double lapse_e; double *shift_e = gkyl_malloc(sizeof(double[3]));
      double **g_e = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) g_e[i] = gkyl_malloc(sizeof(double[3]));
      spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse_e);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift_e);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &g_e);
      double ve[3] = { shift_e[0]/lapse_e, shift_e[1]/lapse_e, shift_e[2]/lapse_e };
      double vsq = 0.0;
      for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) vsq += g_e[i][j]*ve[i]*ve[j];
      double W2 = 1.0/(1.0 - vsq);
      double h0 = 1.0 + app->K_poly*pow(rho, G_etk - 1.0)*G_etk/(G_etk - 1.0);
      double E0 = rho*h0*W2 - app->K_poly*pow(rho, G_etk); // original Eulerian energy density
      double oldp=p, rp, pp, Wnew;
      bool ok=rt_kdrop_es_primitives(G_etk,Kp_etk,E0,(E0+oldp)*(E0+oldp)*vsq,&rp,&pp,&Wnew);
      if (!ok || rp<app->rho_atm || pp<app->p_atm) {
        fprintf(stderr,"Invalid energy-momentum-preserving initial K-drop at (%g,%g,%g)\n",x,y,z);
        exit(EXIT_FAILURE);
      }
      depletion_velocity_scale=(E0+oldp)/(E0+pp);
      rho=rp; p=pp;
      for (int i = 0; i < 3; i++) gkyl_free(g_e[i]);
      gkyl_free(g_e); gkyl_free(shift_e);
    }
    else if (getenv("GKYL_KDROP_ETK") != NULL) {
      // ETK constraint-preserving pressure depletion (arXiv:1111.3344 eq. 91): lower the polytropic
      // constant K -> (1-delta) K, then rescale rho so the (fluid-frame, v=0) energy density
      // T = rho + K rho^Gamma/(Gamma-1) is UNCHANGED. Exact for v=0/areal; approximate in Kerr-Schild
      // where W>1 (use GKYL_KDROP_ETK_EULERIAN for the exact Eulerian target). Newton-solve, general Gamma.
      double T_energy = rho + app->K_poly * pow(rho, G_etk) / (G_etk - 1.0); // original energy density
      double rp = rho;
      for (int it = 0; it < 40; it++) {
        double f  = rp + Kp_etk * pow(rp, G_etk) / (G_etk - 1.0) - T_energy;
        double df = 1.0 + Kp_etk * G_etk * pow(rp, G_etk - 1.0) / (G_etk - 1.0);
        double step = f / df;
        rp -= step;
        if (rp < app->rho_atm) rp = app->rho_atm;
        if (fabs(step) <= 1.0e-14 * rp) break;
      }
      rho = rp;
      p = fmax(Kp_etk * pow(rp, G_etk), app->p_atm);
    }
    else {
      p = fmax((1.0 - local_delta) * p, app->p_atm);
    }
  }

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  double dl = pow(10.0, -8.0);
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &spatial_metric_der);

  // Once and for all with the velocity (I'll make a note for myself for later): In the Kerr-Schild slicing the static star is not at rest w.r.t. the Eulerian observer: 
  // a fluid element following the static Killing vector (zero coordinate velocity) has Eulerian 3-velocity v^i = beta^i / alpha (physical speed = 2m/r, subluminal for a regular star). 
  // v=0 would initialize a moving fluid. (For the areal gauge beta^i=0, so this reduces to v=0.)
  double vel[3] = { depletion_velocity_scale*shift[0] / lapse,
    depletion_velocity_scale*shift[1] / lapse, depletion_velocity_scale*shift[2] / lapse };
  double v_sq = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];

  double W = 1.0 / sqrt(1.0 - v_sq);
  if (v_sq > 1.0 - pow(10.0, -8.0)) W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));

  double h = 1.0 + ((p / rho) * (gamma_h / (gamma_h - 1.0)));
  double sqrt_det = sqrt(spatial_det);
  double cov_vel[3] = { 0.0 };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cov_vel[i] += spatial_metric[i][j] * vel[j];
    }
  }

  double rho_rel = sqrt_det * rho * W;
  double mom_x = sqrt_det * rho * h * (W * W) * cov_vel[0];
  double mom_y = sqrt_det * rho * h * (W * W) * cov_vel[1];
  double mom_z = sqrt_det * rho * h * (W * W) * cov_vel[2];
  double Etot = sqrt_det * ((rho * h * (W * W)) - p - (rho * W));

  fout[0] = rho_rel;
  fout[1] = mom_x; fout[2] = mom_y; fout[3] = mom_z;
  fout[4] = Etot;

  fout[5] = lapse;
  fout[6] = shift[0]; fout[7] = shift[1]; fout[8] = shift[2];

  fout[9]  = spatial_metric[0][0]; fout[10] = spatial_metric[0][1]; fout[11] = spatial_metric[0][2];
  fout[12] = spatial_metric[1][0]; fout[13] = spatial_metric[1][1]; fout[14] = spatial_metric[1][2];
  fout[15] = spatial_metric[2][0]; fout[16] = spatial_metric[2][1]; fout[17] = spatial_metric[2][2];

  fout[18] = extrinsic_curvature[0][0]; fout[19] = extrinsic_curvature[0][1]; fout[20] = extrinsic_curvature[0][2];
  fout[21] = extrinsic_curvature[1][0]; fout[22] = extrinsic_curvature[1][1]; fout[23] = extrinsic_curvature[1][2];
  fout[24] = extrinsic_curvature[2][0]; fout[25] = extrinsic_curvature[2][1]; fout[26] = extrinsic_curvature[2][2];

  fout[27] = in_excision_region ? -1.0 : 1.0;

  fout[28] = lapse_der[0]; fout[29] = lapse_der[1]; fout[30] = lapse_der[2];

  fout[31] = shift_der[0][0]; fout[32] = shift_der[0][1]; fout[33] = shift_der[0][2];
  fout[34] = shift_der[1][0]; fout[35] = shift_der[1][1]; fout[36] = shift_der[1][2];
  fout[37] = shift_der[2][0]; fout[38] = shift_der[2][1]; fout[39] = shift_der[2][2];

  fout[40] = spatial_metric_der[0][0][0]; fout[41] = spatial_metric_der[0][0][1]; fout[42] = spatial_metric_der[0][0][2];
  fout[43] = spatial_metric_der[0][1][0]; fout[44] = spatial_metric_der[0][1][1]; fout[45] = spatial_metric_der[0][1][2];
  fout[46] = spatial_metric_der[0][2][0]; fout[47] = spatial_metric_der[0][2][1]; fout[48] = spatial_metric_der[0][2][2];

  fout[49] = spatial_metric_der[1][0][0]; fout[50] = spatial_metric_der[1][0][1]; fout[51] = spatial_metric_der[1][0][2];
  fout[52] = spatial_metric_der[1][1][0]; fout[53] = spatial_metric_der[1][1][1]; fout[54] = spatial_metric_der[1][1][2];
  fout[55] = spatial_metric_der[1][2][0]; fout[56] = spatial_metric_der[1][2][1]; fout[57] = spatial_metric_der[1][2][2];

  fout[58] = spatial_metric_der[2][0][0]; fout[59] = spatial_metric_der[2][0][1]; fout[60] = spatial_metric_der[2][0][2];
  fout[61] = spatial_metric_der[2][1][0]; fout[62] = spatial_metric_der[2][1][1]; fout[63] = spatial_metric_der[2][1][2];
  fout[64] = spatial_metric_der[2][2][0]; fout[65] = spatial_metric_der[2][2][1]; fout[66] = spatial_metric_der[2][2][2];

  fout[67] = 0.0;
  fout[68] = x; fout[69] = y; fout[70] = z;

  // Frozen-discrete well-balancing reference: t=0 conserved equilibrium fluid.
  fout[71] = rho_rel;
  fout[72] = Etot;
  fout[73] = mom_x;
  fout[74] = mom_y;
  fout[75] = mom_z;

  // Optional momentum perturbation in the inner star
  if (app->perturb_frac > 0.0 && r < 0.5 * app->R_star) {
    fout[1] += app->perturb_frac * rho_rel;
  }

  if (in_excision_region) {
    for (int i = 0; i < 68; i++) fout[i] = 0.0;
    fout[27] = -1.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) { gkyl_free(spatial_metric[i]); gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]); }
  gkyl_free(spatial_metric); gkyl_free(extrinsic_curvature); gkyl_free(lapse_der); gkyl_free(shift_der);
  for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]); gkyl_free(spatial_metric_der[i]); }
  gkyl_free(spatial_metric_der);
}

void
evalVacuumEinsteinInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct gr_tov_coupled_ctx *app = ctx;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  double dl = pow(10.0, -8.0);
  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &spatial_metric_der);

  // Bona-Masso conventions: 
  // D_kij = 1/2 d_k gamma_ij, (1/2) d_i beta^j, A_i = d_i ln(alpha)
  // (So what I calculate fits the vacuum einstein later)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k];
      }
      shift_der[i][j] = 0.5 * shift_der[i][j];
    }
  }
  for (int i = 0; i < 3; i++) {
    lapse_der[i] = lapse_der[i] / lapse;
  }

  double spatial_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised1[k][i][j] = 0.0;
        for (int l = 0; l < 3; l++)
          spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
      }

  double spatial_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised3[i][j][k] = 0.0;
        for (int l = 0; l < 3; l++)
          spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
      }

  double aux_vect[3];
  for (int i = 0; i < 3; i++) {
    aux_vect[i] = 0.0;
    for (int s = 0; s < 3; s++) {
      aux_vect[i] += spatial_metric_der_raised3[i][s][s];
      aux_vect[i] -= spatial_metric_der_raised1[s][s][i];
    }
  }

  fout[0] = spatial_metric[0][0]; fout[1] = spatial_metric[0][1]; fout[2] = spatial_metric[0][2];
  fout[3] = spatial_metric[1][0]; fout[4] = spatial_metric[1][1]; fout[5] = spatial_metric[1][2];
  fout[6] = spatial_metric[2][0]; fout[7] = spatial_metric[2][1]; fout[8] = spatial_metric[2][2];

  fout[9] = lapse;

  fout[10] = extrinsic_curvature[0][0]; fout[11] = extrinsic_curvature[0][1]; fout[12] = extrinsic_curvature[0][2];
  fout[13] = extrinsic_curvature[1][0]; fout[14] = extrinsic_curvature[1][1]; fout[15] = extrinsic_curvature[1][2];
  fout[16] = extrinsic_curvature[2][0]; fout[17] = extrinsic_curvature[2][1]; fout[18] = extrinsic_curvature[2][2];

  fout[19] = spatial_metric_der[0][0][0]; fout[20] = spatial_metric_der[0][0][1]; fout[21] = spatial_metric_der[0][0][2];
  fout[22] = spatial_metric_der[0][1][0]; fout[23] = spatial_metric_der[0][1][1]; fout[24] = spatial_metric_der[0][1][2];
  fout[25] = spatial_metric_der[0][2][0]; fout[26] = spatial_metric_der[0][2][1]; fout[27] = spatial_metric_der[0][2][2];

  fout[28] = spatial_metric_der[1][0][0]; fout[29] = spatial_metric_der[1][0][1]; fout[30] = spatial_metric_der[1][0][2];
  fout[31] = spatial_metric_der[1][1][0]; fout[32] = spatial_metric_der[1][1][1]; fout[33] = spatial_metric_der[1][1][2];
  fout[34] = spatial_metric_der[1][2][0]; fout[35] = spatial_metric_der[1][2][1]; fout[36] = spatial_metric_der[1][2][2];

  fout[37] = spatial_metric_der[2][0][0]; fout[38] = spatial_metric_der[2][0][1]; fout[39] = spatial_metric_der[2][0][2];
  fout[40] = spatial_metric_der[2][1][0]; fout[41] = spatial_metric_der[2][1][1]; fout[42] = spatial_metric_der[2][1][2];
  fout[43] = spatial_metric_der[2][2][0]; fout[44] = spatial_metric_der[2][2][1]; fout[45] = spatial_metric_der[2][2][2];

  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];

  fout[49] = aux_vect[0]; fout[50] = aux_vect[1]; fout[51] = aux_vect[2];

  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];

  fout[55] = shift_der[0][0]; fout[56] = shift_der[0][1]; fout[57] = shift_der[0][2];
  fout[58] = shift_der[1][0]; fout[59] = shift_der[1][1]; fout[60] = shift_der[1][2];
  fout[61] = shift_der[2][0]; fout[62] = shift_der[2][1]; fout[63] = shift_der[2][2];

  if (in_excision_region) {
    for (int i = 0; i < 64; i++) fout[i] = 0.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]); gkyl_free(inv_spatial_metric[i]);
    gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]);
    for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]);
    gkyl_free(spatial_metric_der[i]);
  }
  gkyl_free(spatial_metric); gkyl_free(inv_spatial_metric); gkyl_free(extrinsic_curvature);
  gkyl_free(lapse_der); gkyl_free(shift_der); gkyl_free(spatial_metric_der);
}

// Conformal (77-comp, chi = psi^-2) vacuum-Einstein initial data, 3D. Ported from
// rt_gr_tov_conformal_mixed_cart.c with the third dimension carried through
// (z passed to every spacetime callback; 3D radius used by the fluid IC).
// Selected with GKYL_USE_CONFORMAL=1.
void
evalVacuumEinsteinConformalInit3d(double t, const double* GKYL_RESTRICT xn,
  double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct gr_tov_coupled_ctx *app = ctx;
  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double lapse, psi, chi;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **conformal_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) conformal_metric[i] = gkyl_malloc(sizeof(double[3]));
  double **inv_conformal_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) inv_conformal_metric[i] = gkyl_malloc(sizeof(double[3]));
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));

  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));

  double ***conformal_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    conformal_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) conformal_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  double *psi_der = gkyl_malloc(sizeof(double[3]));
  double *chi_der = gkyl_malloc(sizeof(double[3]));
  double **chi_der2 = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) chi_der2[i] = gkyl_malloc(sizeof(double[3]));

  double dl = 1.0e-8;
  double dl2 = 1.0e-6;
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &conformal_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_conformal_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, dl, dl, dl, &extrinsic_curvature);
  spacetime->conformal_factor_func(spacetime, 0.0, x, y, z, &psi);
  spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, y, z, &chi);
  spacetime->conformal_factor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &psi_der);
  spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &chi_der);
  spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, y, z, dl2, dl2, dl2, &chi_der2);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, dl, dl, dl, &conformal_metric_der);

  double psi4 = psi * psi * psi * psi;
  for (int k = 0; k < 3; k++) psi_der[k] /= psi;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_metric[i][j] /= psi4;
      inv_conformal_metric[i][j] *= psi4;
      for (int k = 0; k < 3; k++) {
        conformal_metric_der[k][i][j] =
          (0.5 * conformal_metric_der[k][i][j] / psi4)
          - (2.0 * psi_der[k] * conformal_metric[i][j]);
      }
      shift_der[i][j] = 0.5 * shift_der[i][j];
    }
  }
  for (int i = 0; i < 3; i++) lapse_der[i] /= lapse;

  double cmd_raised1[3][3][3], cmd_raised3[3][3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        cmd_raised1[k][i][j] = 0.0;
        cmd_raised3[i][j][k] = 0.0;
        for (int l = 0; l < 3; l++) {
          cmd_raised1[k][i][j] += inv_conformal_metric[k][l] * conformal_metric_der[l][i][j];
          cmd_raised3[i][j][k] += inv_conformal_metric[l][k] * conformal_metric_der[i][j][l];
        }
      }

  double conformal_aux_vect[3];
  for (int i = 0; i < 3; i++) {
    conformal_aux_vect[i] = 0.0; // D below is already the conformal metric derivative.
    for (int s = 0; s < 3; s++) {
      conformal_aux_vect[i] += cmd_raised3[i][s][s];
      conformal_aux_vect[i] -= cmd_raised1[s][s][i];
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fout[3 * i + j] = conformal_metric[i][j];
      fout[10 + 3 * i + j] = extrinsic_curvature[i][j];
    }
  }
  fout[9] = lapse;
  for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        fout[19 + (9 * k) + (3 * i) + j] = conformal_metric_der[k][i][j];

  fout[46] = lapse_der[0]; fout[47] = lapse_der[1]; fout[48] = lapse_der[2];
  fout[49] = conformal_aux_vect[0]; fout[50] = conformal_aux_vect[1]; fout[51] = conformal_aux_vect[2];
  fout[52] = shift[0]; fout[53] = shift[1]; fout[54] = shift[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      fout[55 + 3 * i + j] = shift_der[i][j];

  fout[64] = chi;
  fout[65] = chi_der[0]; fout[66] = chi_der[1]; fout[67] = chi_der[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      fout[68 + 3 * i + j] = chi_der2[i][j];

  if (in_excision_region) {
    for (int i = 0; i < 77; i++) fout[i] = 0.0;
  }

  gkyl_free(shift);
  for (int i = 0; i < 3; i++) {
    gkyl_free(conformal_metric[i]); gkyl_free(inv_conformal_metric[i]);
    gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]);
    gkyl_free(chi_der2[i]);
    for (int j = 0; j < 3; j++) gkyl_free(conformal_metric_der[i][j]);
    gkyl_free(conformal_metric_der[i]);
  }
  gkyl_free(conformal_metric); gkyl_free(inv_conformal_metric); gkyl_free(extrinsic_curvature);
  gkyl_free(lapse_der); gkyl_free(shift_der); gkyl_free(conformal_metric_der);
  gkyl_free(psi_der); gkyl_free(chi_der); gkyl_free(chi_der2);
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = (force_write) ? iot->curr : iot->curr - 1;
    gkyl_moment_app_write(app, t_curr, frame);
    gkyl_moment_app_write_field_energy(app);
    gkyl_moment_app_write_integrated_mom(app);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_moment_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_moment_app_calc_integrated_mom(app, t_curr);
  }
}

// [DIAG] Per-step mixed TOV tracker. Fluid quantities come from species 0; lapse comes from
// the live Einstein species 1, matching the collapse diagnostic.
static double diag_SumD0 = -1.0;
static void
diag_read_einstein_state(const double *q, double g[3][3], double gi[3][3],
  double K[3][3], double D[3][3][3], double A[3], double V[3]);

static void
diag_force_balance(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  if (getenv("GKYL_FORCE_DIAG") == NULL || ctx->fluid_eqn == NULL) {
    return;
  }

  static long call_count = 0;
  call_count += 1;
  int every = 1;
  const char *every_env = getenv("GKYL_FORCE_DIAG_EVERY");
  if (every_env != NULL) {
    every = atoi(every_env);
    if (every < 1) every = 1;
  }
  if (((call_count - 1) % every) != 0) {
    return;
  }

  const struct gkyl_array *fluid = gkyl_moment_app_get_write_array_species(app, 0);
  const double dx = app->grid.dx[0];
  const double radii[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0.5 * ctx->R_star, ctx->R_core, ctx->R_star };
  const int nr = (int) (sizeof(radii) / sizeof(radii[0]));

  int rank = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    printf("[FORCEDIAG] t=%.6e  positive-x axis proxy; inward means rhs_Sx<0; WB_disabled=%d\n",
      t_curr, (int) gkyl_gr_euler_wb_disabled(ctx->fluid_eqn));
  }

  for (int ir = 0; ir < nr; ir++) {
    double r0 = radii[ir];
    double best_score_local = DBL_MAX;
    double vals_local[12] = { 0.0 };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[0] <= app->local.lower[0] || iter.idx[0] >= app->local.upper[0]) {
        continue;
      }

      double xc[3];
      gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
      if (xc[0] <= 0.0) {
        continue;
      }

      double score = fabs(xc[0] - r0) + fabs(xc[1]) + fabs(xc[2]);
      if (score < best_score_local) {
        int idxp[3] = { iter.idx[0] + 1, iter.idx[1], iter.idx[2] };
        int idxm[3] = { iter.idx[0] - 1, iter.idx[1], iter.idx[2] };
        long loc0 = gkyl_range_idx(&app->local_ext, iter.idx);
        long locp = gkyl_range_idx(&app->local_ext, idxp);
        long locm = gkyl_range_idx(&app->local_ext, idxm);
        const double *q0 = gkyl_array_cfetch(fluid, loc0);
        const double *qp = gkyl_array_cfetch(fluid, locp);
        const double *qm = gkyl_array_cfetch(fluid, locm);

        double prim0[71] = { 0.0 }, primp[71] = { 0.0 }, primm[71] = { 0.0 };
        double Fp[71] = { 0.0 }, Fm[71] = { 0.0 }, src[ctx->fluid_eqn->num_equations];
        for (int m = 0; m < ctx->fluid_eqn->num_equations; m++) src[m] = 0.0;
        gkyl_gr_euler_prim_vars(ctx->gas_gamma, q0, prim0);
        gkyl_gr_euler_prim_vars(ctx->gas_gamma, qp, primp);
        gkyl_gr_euler_prim_vars(ctx->gas_gamma, qm, primm);
        gkyl_gr_euler_flux(ctx->gas_gamma, qp, Fp);
        gkyl_gr_euler_flux(ctx->gas_gamma, qm, Fm);
        gkyl_wv_eqn_source(ctx->fluid_eqn, q0, src);

        double flux_div_sx = (Fp[1] - Fm[1]) / (2.0 * dx);
        double rhs_sx = -flux_div_sx + src[1];
        double dpdr = (primp[4] - primm[4]) / (2.0 * dx);
        double accel_proxy = rhs_sx / (fabs(q0[0]) + 1.0e-300);

        best_score_local = score;
        vals_local[0] = xc[0];
        vals_local[1] = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);
        vals_local[2] = prim0[0];
        vals_local[3] = prim0[4];
        vals_local[4] = prim0[1];
        vals_local[5] = dpdr;
        vals_local[6] = -flux_div_sx;
        vals_local[7] = src[1];
        vals_local[8] = rhs_sx;
        vals_local[9] = accel_proxy;
        vals_local[10] = q0[0];
        vals_local[11] = q0[1];
      }
    }

    double best_score = DBL_MAX;
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &best_score_local, &best_score);

    double vals_candidate[12] = { 0.0 };
    if (fabs(best_score_local - best_score) < 1.0e-12 * fmax(1.0, best_score)) {
      for (int m = 0; m < 12; m++) vals_candidate[m] = vals_local[m];
    }
    double vals[12] = { 0.0 };
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 12, vals_candidate, vals);

    if (rank == 0) {
      printf("[FORCEDIAG] r0=%.6e x=%.6e r=%.6e rho=%.6e p=%.6e vx=%.6e"
        " dpdx=%.6e -dFxx_dx=%.6e srcSx=%.6e rhsSx=%.6e rhsSx/D=%.6e D=%.6e Sx=%.6e\n",
        r0, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6],
        vals[7], vals[8], vals[9], vals[10], vals[11]);
    }
  }
  if (rank == 0) fflush(stdout);
}

static bool diag_spacetime_in_excision_region(double excision_threshold, const double *q);

static void
diag_track(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  static long call_count = 0;
  call_count += 1;
  int every = 1;
  const char *every_env = getenv("GKYL_DIAG_TRACK_EVERY");
  if (every_env) {
    every = atoi(every_env);
    if (every < 1) {
      every = 1;
    }
  }
  if ((call_count - 1) % every != 0) {
    return;
  }
  const char *tmin_env = getenv("GKYL_DIAG_TRACK_TMIN");
  if (tmin_env && t_curr < atof(tmin_env)) {
    return;
  }

  const struct gkyl_array *fluid = gkyl_moment_app_get_write_array_species(app, 0);
  const struct gkyl_array *einstein = gkyl_moment_app_get_write_array_species(app, 1);

  double SumD_local = 0.0;
  double Dmax_local = 0.0, r_Dmax_local = 0.0;
  double Dc_local = 0.0, lapse_c_local = 0.0, r_c_local = 1.0e300;
  double core_mom_local = 0.0, surf_mom_local = 0.0;
  double glob_mom_local = 0.0, r_glob_mom_local = 0.0;
  double lapse_min_local = 1.0e300, lapse_max_local = -1.0e300;
  double rho_max_local = 0.0, p_max_local = 0.0, v2_max_local = 0.0;
  double rho_c_local = 0.0, p_c_local = 0.0;
  double M0_r2_local = 0.0, M0_r3_local = 0.0, M0_r4_local = 0.0, M0_r5_local = 0.0;
  int n_bad_local = 0, n_floor_local = 0;
  bool phys_diag = getenv("GKYL_PHYS_DIAG") != NULL;

  double R_core = 0.5 * ctx->R_star;
  double R_star = ctx->R_star;
  double surf_lo = R_star - 0.15 * R_star, surf_hi = R_star + 0.15 * R_star;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(fluid, loc);
    const double *es = gkyl_array_cfetch(einstein, loc);

    double xc[3];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);

    if (diag_spacetime_in_excision_region(ctx->excision_threshold, es)) {
      continue;
    }

    double D = fl[0];
    double lapse = es[9];
    double mom = sqrt(fl[1] * fl[1] + fl[2] * fl[2] + fl[3] * fl[3]);

    if (!isfinite(D) || !isfinite(lapse) || !isfinite(mom)) n_bad_local += 1;

    SumD_local += D;
    if (r < 2.0) M0_r2_local += D;
    if (r < 3.0) M0_r3_local += D;
    if (r < 4.0) M0_r4_local += D;
    if (r < 5.0) M0_r5_local += D;
    if (D > Dmax_local) { Dmax_local = D; r_Dmax_local = r; }

    bool is_center_cell = r < r_c_local;
    if (is_center_cell) { r_c_local = r; Dc_local = D; lapse_c_local = lapse; }

    if (phys_diag) {
      double prim[71] = { 0.0 };
      gkyl_gr_euler_prim_vars(ctx->gas_gamma, fl, prim);
      double rho = prim[0], p = prim[4];
      double vel[3] = { prim[1], prim[2], prim[3] };
      double g[3][3], gi[3][3], Ktmp[3][3], Dtmp[3][3][3], Atmp[3], Vtmp[3];
      diag_read_einstein_state(es, g, gi, Ktmp, Dtmp, Atmp, Vtmp);
      double v2 = 0.0;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          v2 += g[i][j] * vel[i] * vel[j];
      if (isfinite(rho) && rho > rho_max_local) rho_max_local = rho;
      if (isfinite(p) && p > p_max_local) p_max_local = p;
      // maxv2 only where the C2P recovered a real state: a cell at the pressure or density floor, or in the
      // atmosphere band (rho <= 100 rho_atm, the validator's band), has no meaningful velocity (AUDIT_FLUID_STATE_WRITES,
      // End-to-end result, Q4); such cells are counted (floored=) instead.
      if (!(rho > 100.0 * ctx->rho_atm) || !(p > ctx->p_atm)) n_floor_local += 1;
      else if (isfinite(v2) && v2 > v2_max_local) v2_max_local = v2;
      if (is_center_cell) { rho_c_local = rho; p_c_local = p; }
    }

    if (r < R_core && mom > core_mom_local) core_mom_local = mom;
    if (r > surf_lo && r < surf_hi && mom > surf_mom_local) surf_mom_local = mom;
    if (mom > glob_mom_local) { glob_mom_local = mom; r_glob_mom_local = r; }

    if (r < 1.5 * R_star) {
      if (lapse < lapse_min_local) lapse_min_local = lapse;
      if (lapse > lapse_max_local) lapse_max_local = lapse;
    }
  }

  SumD_local *= app->grid.cellVolume;

  double SumD = 0.0;
  double Dmax = 0.0, core_mom = 0.0, surf_mom = 0.0, glob_mom = 0.0;
  double lapse_min = 0.0, lapse_max = 0.0, r_c = 0.0;
  double rho_max = 0.0, p_max = 0.0, v2_max = 0.0;
  double M0_r2 = 0.0, M0_r3 = 0.0, M0_r4 = 0.0, M0_r5 = 0.0;
  int n_bad = 0, n_floor = 0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &SumD_local, &SumD);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Dmax_local, &Dmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &core_mom_local, &core_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &surf_mom_local, &surf_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &glob_mom_local, &glob_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &lapse_min_local, &lapse_min);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &lapse_max_local, &lapse_max);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_c_local, &r_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &rho_max_local, &rho_max);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &p_max_local, &p_max);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &v2_max_local, &v2_max);
  M0_r2_local *= app->grid.cellVolume;
  M0_r3_local *= app->grid.cellVolume;
  M0_r4_local *= app->grid.cellVolume;
  M0_r5_local *= app->grid.cellVolume;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &M0_r2_local, &M0_r2);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &M0_r3_local, &M0_r3);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &M0_r4_local, &M0_r4);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &M0_r5_local, &M0_r5);
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_SUM, 1, &n_bad_local, &n_bad);
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_SUM, 1, &n_floor_local, &n_floor);

  double r_Dmax_candidate = 1.0e300;
  if (Dmax_local == Dmax) r_Dmax_candidate = r_Dmax_local;
  double r_glob_mom_candidate = 1.0e300;
  if (glob_mom_local == glob_mom) r_glob_mom_candidate = r_glob_mom_local;
  double Dc_candidate = 0.0, lapse_c_candidate = 0.0;
  double rho_c_candidate = 0.0, p_c_candidate = 0.0;
  if (r_c_local == r_c) {
    Dc_candidate = Dc_local;
    lapse_c_candidate = lapse_c_local;
    rho_c_candidate = rho_c_local;
    p_c_candidate = p_c_local;
  }

  double r_Dmax = 0.0, r_glob_mom = 0.0, Dc = 0.0, lapse_c = 0.0;
  double rho_c = 0.0, p_c = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_Dmax_candidate, &r_Dmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_glob_mom_candidate, &r_glob_mom);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Dc_candidate, &Dc);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &lapse_c_candidate, &lapse_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &rho_c_candidate, &rho_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &p_c_candidate, &p_c);

  if (diag_SumD0 < 0.0) diag_SumD0 = SumD;
  double dSumD = SumD - diag_SumD0;

  int rank = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    printf("[DIAG] t=%.4f  Dc=%.6e lapse_c=%.6e  Dmax=%.6e @r=%.3f  core|mom|=%.3e surf|mom|=%.3e glob|mom|=%.3e @r=%.2f  lapse[%.4f,%.4f]  SumD=%.8e dSumD=%.3e bad=%d\n",
      t_curr, Dc, lapse_c, Dmax, r_Dmax, core_mom, surf_mom, glob_mom, r_glob_mom, lapse_min, lapse_max, SumD, dSumD, n_bad);
    if (phys_diag) {
      printf("[PHYSDIAG] t=%.4f  rho_c=%.6e p_c=%.6e rho_max=%.6e p_max=%.6e maxv2=%.6e  M0(<2,3,4,5)=(%.6e %.6e %.6e %.6e) floored=%d\n",
        t_curr, rho_c, p_c, rho_max, p_max, v2_max, M0_r2, M0_r3, M0_r4, M0_r5, n_floor);
    }
    fflush(stdout);
  }

  diag_force_balance(app, ctx, t_curr);
}

static double
diag_sym3_inverse(const double g[3][3], double gi[3][3])
{
  double c00 = g[1][1] * g[2][2] - g[1][2] * g[2][1];
  double c01 = g[1][2] * g[2][0] - g[1][0] * g[2][2];
  double c02 = g[1][0] * g[2][1] - g[1][1] * g[2][0];
  double det = g[0][0] * c00 + g[0][1] * c01 + g[0][2] * c02;
  double idet = 0.0;
  if (det != 0.0) {
    idet = 1.0 / det;
  }

  gi[0][0] = c00 * idet;
  gi[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) * idet;
  gi[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) * idet;
  gi[0][1] = gi[1][0] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) * idet;
  gi[0][2] = gi[2][0] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) * idet;
  gi[1][2] = gi[2][1] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) * idet;

  return det;
}

static bool diag_is_conformal = false;
static bool
diag_spacetime_in_excision_region(double excision_threshold, const double *q)
{
  return gkyl_gr_einstein_is_excised(excision_threshold, q, diag_is_conformal ? 77 : 64);
}

static bool
diag_constraint_stencil_touches_excision(gkyl_moment_app *app,
  const struct gkyl_array *einstein, const int *idx, double excision_threshold)
{
  const double *q0 = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idx));
  if (diag_spacetime_in_excision_region(excision_threshold, q0)) {
    return true;
  }

  for (int d = 0; d < app->ndim; d++) {
    int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
    for (int e = 0; e < app->ndim; e++) {
      idxp[e] = idx[e];
      idxm[e] = idx[e];
    }
    idxp[d] += 1;
    idxm[d] -= 1;

    const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
    const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
    if (diag_spacetime_in_excision_region(excision_threshold, qp) ||
        diag_spacetime_in_excision_region(excision_threshold, qm)) {
      return true;
    }
  }

  return false;
}

static void diag_read_physical_metric(const double *q, double g[3][3], double gi[3][3]);

// [TRES] sqrt(-g) * T^i_j for one cell: the flux whose divergence enters nabla_mu T^mu_j.
// T^i_j = g_{j nu} T^{i nu} = beta_j T^{i0} + gamma_{jk} T^{ik}, and sqrt(-g) = alpha sqrt(gamma).
static void
diag_sqrtg_Tmixed(double gas_gamma, const double *qf, const double *qe, double out[3][3])
{
  double **T = gkyl_malloc(sizeof(double*[4]));
  for (int a = 0; a < 4; a++) {
    T[a] = gkyl_malloc(sizeof(double[4]));
  }
  gkyl_gr_euler_stress_energy_tensor(gas_gamma, qf, &T);

  double g[3][3], gi[3][3];
  diag_read_physical_metric(qe, g, gi);
  double det = diag_sym3_inverse(g, gi);
  double sdet = 0.0;
  if (det > 0.0) {
    sdet = sqrt(det);
  }
  double lapse = qe[9];
  double shift[3] = { qe[52], qe[53], qe[54] };

  double shift_low[3] = { 0.0, 0.0, 0.0 };
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      shift_low[j] += g[j][k] * shift[k];
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double t = shift_low[j] * T[i + 1][0];
      for (int k = 0; k < 3; k++) {
        t += g[j][k] * T[i + 1][k + 1];
      }
      out[i][j] = lapse * sdet * t;
    }
  }

  for (int a = 0; a < 4; a++) {
    gkyl_free(T[a]);
  }
  gkyl_free(T);
}

// Set once at startup: the monitor reads the PHYSICAL metric in both formulations
// so the constraint diagnostics are directly comparable. Conformal state stores
// chi = psi^-2 at slot 64 (d_k chi at 65-67) and tilde(gamma)_ij at 0-8, with
// gamma_ij = chi^-2 tilde(gamma)_ij; K_ij is already physical at 10-18. The
// derivative reduction obeys d_k gamma_ij = chi^-2 (2 tilde(D)_kij - 2 (d_k chi/chi)
// tilde(gamma)_ij), i.e. D_kij = chi^-2 (tilde(D)_kij - (d_k chi/chi) tilde(gamma)_ij).
// (Matches sync_fluid_metric_from_einstein in mom_update_ssp_rk.c.)
static void
diag_read_physical_metric(const double *q, double g[3][3], double gi[3][3])
{
  double psi4 = 1.0;
  if (diag_is_conformal) {
    double chi = q[64];
    if (!(chi > 0.0)) {
      chi = 1.0;
    }
    double inv_chi = 1.0 / chi;
    psi4 = inv_chi * inv_chi;
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      g[i][j] = psi4 * q[3 * i + j];
    }
  }
  diag_sym3_inverse(g, gi);
}

static void
diag_read_einstein_state(const double *q, double g[3][3], double gi[3][3],
  double K[3][3], double D[3][3][3], double A[3], double V[3])
{
  double psi4 = 1.0, dchi_over_chi[3] = { 0.0, 0.0, 0.0 };
  if (diag_is_conformal) {
    double chi = q[64];
    if (!(chi > 0.0)) chi = 1.0;
    double inv_chi = 1.0 / chi;
    psi4 = inv_chi * inv_chi;
    for (int k = 0; k < 3; k++) dchi_over_chi[k] = q[65 + k] * inv_chi;
  }

  double gtil[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gtil[i][j] = q[3 * i + j];
      g[i][j] = psi4 * gtil[i][j];
      K[i][j] = q[10 + 3 * i + j];
    }
  }
  diag_sym3_inverse(g, gi);

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double Dt = q[19 + 9 * k + 3 * i + j];
        D[k][i][j] = diag_is_conformal
          ? psi4 * (Dt - dchi_over_chi[k] * gtil[i][j])
          : Dt;
      }
    }
    A[k] = q[46 + k];
    V[k] = q[49 + k]; // C^V below uses the matching physical or conformal D contraction.
  }
}


static double
diag_aux_constraint(const double *q, const double gi[3][3], const double D[3][3][3], const double V[3])
{
  double Cv = 0.0;
  if (diag_is_conformal) {
    double gtilde[3][3], gtilde_inv[3][3], Dtilde[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        gtilde[i][j] = q[3 * i + j];
      }
    }
    diag_sym3_inverse(gtilde, gtilde_inv);

    for (int k = 0; k < 3; k++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          Dtilde[k][i][j] = q[19 + 9 * k + 3 * i + j];
        }
      }
    }

    double Dtilde_raised1[3][3][3], Dtilde_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          Dtilde_raised1[k][i][j] = 0.0;
          Dtilde_raised3[i][j][k] = 0.0;
          for (int l = 0; l < 3; l++) {
            Dtilde_raised1[k][i][j] += gtilde_inv[k][l] * Dtilde[l][i][j];
            Dtilde_raised3[i][j][k] += gtilde_inv[l][k] * Dtilde[i][j][l];
          }
        }
      }
    }

    for (int i = 0; i < 3; i++) {
      double aux_i = 0.0;
      for (int s = 0; s < 3; s++) {
        aux_i += Dtilde_raised3[i][s][s];
        aux_i -= Dtilde_raised1[s][s][i];
      }
      double cvi = V[i] - aux_i;
      if (fabs(cvi) > Cv) {
        Cv = fabs(cvi);
      }
    }
  }
  else {
    double D_raised1[3][3][3], D_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          D_raised1[k][i][j] = 0.0;
          D_raised3[i][j][k] = 0.0;
          for (int l = 0; l < 3; l++) {
            D_raised1[k][i][j] += gi[k][l] * D[l][i][j];
            D_raised3[i][j][k] += gi[l][k] * D[i][j][l];
          }
        }
      }
    }

    for (int i = 0; i < 3; i++) {
      double aux_i = 0.0;
      for (int s = 0; s < 3; s++) {
        aux_i += D_raised3[i][s][s];
        aux_i -= D_raised1[s][s][i];
      }
      double cvi = V[i] - aux_i;
      if (fabs(cvi) > Cv) {
        Cv = fabs(cvi);
      }
    }
  }

  return Cv;
}

// [PROJ] Christoffel at grid index `idx`, built from D_proj = 1/2 * central-difference of the
// physical metric gamma_ij at the neighbours, instead of the evolved reduction variable D_kij
// (or conformal Dtilde_kij plus chi derivatives). Depends only on primary metric variables, so
// comparing H_proj to H isolates the part driven by reduction drift C_gamma.
static void
diag_christoffel_proj(const struct gkyl_array *ef, const struct gkyl_range *rng,
  const int *idx, const double *dx, double Gamma[3][3][3])
{
  const double *q = gkyl_array_cfetch(ef, gkyl_range_idx(rng, idx));
  double g[3][3], gi[3][3];
  diag_read_physical_metric(q, g, gi);

  double D[3][3][3];
  for (int k = 0; k < 3; k++) {
    int ip[GKYL_MAX_DIM], im[GKYL_MAX_DIM];
    for (int d = 0; d < 3; d++) {
      ip[d] = idx[d];
      im[d] = idx[d];
    }
    ip[k] = idx[k] + 1;
    im[k] = idx[k] - 1;
    const double *qp = gkyl_array_cfetch(ef, gkyl_range_idx(rng, ip));
    const double *qm = gkyl_array_cfetch(ef, gkyl_range_idx(rng, im));
    double gp[3][3], gip[3][3], gm[3][3], gim[3][3];
    diag_read_physical_metric(qp, gp, gip);
    diag_read_physical_metric(qm, gm, gim);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        D[k][i][j] = (gp[i][j] - gm[i][j]) / (4.0 * dx[k]); // 1/2 * (gamma_p - gamma_m)/(2 dx)
      }
    }
  }

  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      for (int c = 0; c < 3; c++) {
        Gamma[a][b][c] = 0.0;
        for (int m = 0; m < 3; m++) {
          Gamma[a][b][c] += gi[a][m] * (D[b][m][c] + D[c][m][b] - D[m][b][c]);
        }
      }
    }
  }
}

static void
diag_christoffel_from_state(const double *q, double Gamma[3][3][3])
{
  double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
  diag_read_einstein_state(q, g, gi, K, D, A, V);

  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      for (int c = 0; c < 3; c++) {
        Gamma[a][b][c] = 0.0;
        for (int m = 0; m < 3; m++) {
          Gamma[a][b][c] += gi[a][m] * (D[b][m][c] + D[c][m][b] - D[m][b][c]);
        }
      }
    }
  }
}

static void
diag_compute_Q(const double *q, double Q[3][3])
{
  double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
  diag_read_einstein_state(q, g, gi, K, D, A, V);

  double trK = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      trK += gi[i][j] * K[i][j];
    }
  }

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      Q[j][i] = 0.0;
      for (int a = 0; a < 3; a++) {
        Q[j][i] += gi[j][a] * K[a][i];
      }
      if (j == i) {
        Q[j][i] -= trK;
      }
    }
  }
}

static void
diag_matter_projections(double gas_gamma, const double *qf, const double g[3][3],
  double *E, double S_cov[3])
{
  double prim[71] = { 0.0 };
  gkyl_gr_euler_prim_vars(gas_gamma, qf, prim);
  double rho = prim[0];
  double p = prim[4];
  double vel[3] = { prim[1], prim[2], prim[3] };

  double v_sq = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += g[i][j] * vel[i] * vel[j];
    }
  }
  if (v_sq > 1.0 - 1.0e-12) {
    v_sq = 1.0 - 1.0e-12;
  }

  double W = 1.0 / sqrt(1.0 - v_sq);
  double h = gkyl_gr_euler_specific_enthalpy(gas_gamma, rho, p);
  *E = rho * h * W * W - p;

  for (int i = 0; i < 3; i++) {
    S_cov[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      S_cov[i] += g[i][j] * rho * h * W * W * vel[j];
    }
  }
}

static bool
sync_initial_fluid_metric(gkyl_moment_app *app, int fluid_idx, int einstein_idx)
{
  bool is_conformal = app->species[einstein_idx].has_vacuum_einstein_conformal;
  double excision_threshold = is_conformal
    ? app->species[einstein_idx].vacuum_einstein_conformal_excision_threshold
    : app->species[einstein_idx].vacuum_einstein_excision_threshold;

  const struct gkyl_array *einstein_f = app->species[einstein_idx].fcurr;
  struct gkyl_array *fluid_f = app->species[fluid_idx].fcurr;

  struct gkyl_range_iter iter;
  int64_t repaired=0, invalid=0;
  double max_momentum_change=0.0;
  gkyl_range_iter_init(&iter, &app->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *qe = gkyl_array_cfetch(einstein_f, loc);
    double *qf = gkyl_array_fetch(fluid_f, loc);
    sync_fluid_metric_from_einstein(qe, qf, excision_threshold, is_conformal);
    // Ghost conserved/reference slots are populated by the species BC below.
    if (!gkyl_range_contains_idx(&app->local,iter.idx)) continue;
    if (qf[27] < 1e-8) continue;
    double zero[5]={0};
    if (gkyl_mp_scheme_gr_euler_pcp_theta(qf,zero,&qf[9]) == 1.0) continue;
    // Separate FV projections of geometry and matter need not lie in the
    // same metric-dependent cone. Project only the initial momentum, keeping
    // D and tau exactly. Freeze the identical admissible state in the WB slots.
    double safe[5];
    gkyl_mp_scheme_gr_euler_pcp_safe_anchor(qf,&qf[9],safe);
    if (safe[0]!=qf[0] || safe[4]!=qf[4] ||
        gkyl_mp_scheme_gr_euler_pcp_theta(safe,zero,&qf[9]) != 1.0) {
      invalid++;
      continue;
    }
    double norm=0.0, change=0.0;
    for (int i=0; i<3; ++i) {
      norm+=qf[1+i]*qf[1+i];
      change+=(qf[1+i]-safe[1+i])*(qf[1+i]-safe[1+i]);
      qf[1+i]=safe[1+i]; qf[73+i]=safe[1+i];
    }
    max_momentum_change=fmax(max_momentum_change,sqrt(change/fmax(norm,DBL_MIN)));
    repaired++;
  }
  int64_t repaired_all=0, invalid_all=0;
  double max_change_all=0;
  gkyl_comm_allreduce(app->comm,GKYL_INT_64,GKYL_SUM,1,&repaired,&repaired_all);
  gkyl_comm_allreduce(app->comm,GKYL_INT_64,GKYL_SUM,1,&invalid,&invalid_all);
  gkyl_comm_allreduce(app->comm,GKYL_DOUBLE,GKYL_MAX,1,&max_momentum_change,&max_change_all);
  gkyl_moment_app_cout(app,stdout,
    "[IC-ADMISSIBILITY] momentum_projected=%ld invalid=%ld max_relative_change=%.17g dD=0 dTau=0\n",
    (long)repaired_all,(long)invalid_all,max_change_all);
  moment_species_apply_bc(app,0.0,&app->species[fluid_idx],fluid_f);
  return invalid_all==0;
}

// [CONS] ADM + first-order reduction monitor for the mixed 3D TOV run.
// This is diagnostic-only: it does not feed back into the update. It checks the
// quantities that would have caught the old Cartesian TOV metric bug immediately:
// C_kij = D_kij - 1/2 d_k gamma_ij, C_i = A_i - d_i log(alpha), and
// C^V_i = V_i - (D_is{}^s - D^s{}_si). It also reports ADM Hamiltonian and
// momentum residuals using the evolved first-order metric derivatives.
static void
constraint_monitor(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  static long call_count = 0;
  static int delta_seeded = 0;
  static struct gkyl_array *delta_H0 = NULL;
  static struct gkyl_array *delta_T0 = NULL;
  bool delta_diag = (getenv("GKYL_CONSTRAINT_DELTA") != NULL);
  bool need_delta_seed = delta_diag && !delta_seeded;
  bool force_final = false;
  const char *final_env = getenv("GKYL_CONSTRAINT_DIAG_FINAL");
  if (final_env && t_curr >= ctx->t_end - 1.0e-12) {
    force_final = true;
  }
  int every = 1;
  const char *every_env = getenv("GKYL_CONSTRAINT_DIAG_EVERY");
  if (every_env) {
    every = atoi(every_env);
    if (every < 1) {
      every = 1;
    }
  }
  call_count += 1;
  bool emit_diag = true;
  if (!force_final && (call_count - 1) % every != 0) emit_diag = false;
  const char *tmin_env = getenv("GKYL_CONSTRAINT_DIAG_TMIN");
  if (!force_final && tmin_env && t_curr < atof(tmin_env)) emit_diag = false;
  if (!emit_diag && !need_delta_seed) {
    return;
  }

  const struct gkyl_array *fluid = gkyl_moment_app_get_write_array_species(app, 0);
  const struct gkyl_array *einstein = gkyl_moment_app_get_write_array_species(app, 1);
  if (delta_diag && delta_H0 == NULL) {
    delta_H0 = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
    delta_T0 = gkyl_array_new(GKYL_DOUBLE, 3, app->local_ext.volume);
    gkyl_array_clear(delta_H0, 0.0);
    gkyl_array_clear(delta_T0, 0.0);
  }

  double Hmax_local = 0.0, Mmax_local = 0.0;
  double Cgmax_local = 0.0, Camax_local = 0.0, Cvmax_local = 0.0;
  double trKmax_local = 0.0, trKc_local = 0.0, r_c_local = 1.0e300;
  double r_H_local = 0.0, r_M_local = 0.0;
  // [DETGT] Raw determinant monitor for the rescaled spatial metric in slots 0-8.
  // This is an algebraic unimodular constraint only when chi is evolved as a BSSN
  // volume variable. The current production convention freezes chi, so det(gamma-tilde)
  // may carry physical/gauge volume evolution and must not be interpreted alone as an error.
  double detgt_maxdev_local = 0.0, detgt_at_maxdev_local = 1.0, detgt_r_local = 0.0;
  double R_at_H_local = 0.0, trK2_at_H_local = 0.0, K2_at_H_local = 0.0, E_at_H_local = 0.0;
  double H_at_c_local = 0.0, R_at_c_local = 0.0, trK2_at_c_local = 0.0, K2_at_c_local = 0.0, E_at_c_local = 0.0;
  int n_bad_local = 0;

  // [LINE] H sampled cell-by-cell along the +x axis row, to distinguish a GRID-SCALE (2*dx
  // zig-zag) instability from a smooth long-wavelength one. A growth rate scaling like 1/h
  // implies the unstable mode sits at the Nyquist wavelength, which shows up here as sign
  // alternation between neighbouring cells. Serial diagnostic (rank 0 / single block).
  // Env GKYL_CONSTRAINT_LINE.
  #define NLINE 256
  double lineH[NLINE], lineX[NLINE];
  int nline = 0;
  bool line_diag = (getenv("GKYL_CONSTRAINT_LINE") != NULL);

  // [TRES] covariant matter-conservation residual (env GKYL_TRES; costs 6 extra C2P calls per cell)
  bool tres_diag = (getenv("GKYL_TRES") != NULL);
  bool hterm_diag = (getenv("GKYL_HTERM") != NULL);
  // [REGION5] see below: 0 core | 1 surface | 2 near-ext | 3 far-ext | 4 outer-boundary shell.
  #define NSPLIT 5
  double tresmax_local = 0.0, tres2_local = 0.0, tresvol_local = 0.0;
  double tresmax_split_local[NSPLIT] = { 0.0 };
  double tres2_split_local[NSPLIT] = { 0.0 };
  double tresvol_split_local[NSPLIT] = { 0.0 };
  double dHmax_int_local = 0.0, dH2_local = 0.0, dHvol_local = 0.0;
  double dHmax_split_local[NSPLIT] = { 0.0 };
  double dH2_split_local[NSPLIT] = { 0.0 };
  double dHvol_split_local[NSPLIT] = { 0.0 };
  double dTresmax_int_local = 0.0, dTres2_local = 0.0, dTresvol_local = 0.0;
  double dTresmax_split_local[NSPLIT] = { 0.0 };
  double dTres2_split_local[NSPLIT] = { 0.0 };
  double dTresvol_split_local[NSPLIT] = { 0.0 };

  // Interior / exterior split: the global max|H| is typically dominated by the
  // star-surface kink and the outer-boundary fill error, hiding whether the CORE
  // is converging. Report the interior (r <= r_int, default just inside the star
  // surface; env GKYL_CONSTRAINT_RINT) maxima separately.
  double r_int = ctx->R_star - 1.0;
  const char *rint_env = getenv("GKYL_CONSTRAINT_RINT");
  if (rint_env) r_int = atof(rint_env);
  double Hint_local = 0.0, Mint_local = 0.0, Cgint_local = 0.0;
  // L2 (volume-weighted) interior norm of H. A max-norm is a poor convergence
  // diagnostic (it is grid-alignment sensitive and samples a single cell); the
  // L2 norm is the standard measure for a self-convergence test.
  double H2_local = 0.0, vol_local = 0.0;
  double Hnorm_int_local = 0.0, H2_noctr_local = 0.0, vol_noctr_local = 0.0;
  // [CAUSAL] L2|H| inside the shrinking radius r_causal = box_half - speed*t, i.e. the region the
  // outer boundary cannot yet have influenced. If this stays clean and 2nd-order at ALL times, the
  // interior scheme is perfect and ONLY the boundary is broken ⇒ a good outgoing BC (Sommerfeld) is
  // guaranteed to help, without implementing it. speed = fastest signal (gauge speed up to alpha*sqrt(f),
  // ~2-3 for 1+log); GKYL_CAUSAL_SPEED default 3.0 is conservative (over-estimates ⇒ smaller safe region).
  double H2_causal_local = 0.0, vol_causal_local = 0.0, Hmax_causal_local = 0.0;
  double causal_speed = 3.0;
  const char *cs_env = getenv("GKYL_CAUSAL_SPEED");
  if (cs_env) {
    causal_speed = atof(cs_env);
  }
  double r_causal = app->grid.upper[0] - causal_speed * t_curr;
  // [REGION5] Finer radial/boundary split to separate a BC problem (bad at the outer-boundary
  // shell) from a surface/metric-matching problem (bad just outside the star) from a genuine
  // reduction drift (bad H co-located with bad Cgamma/CA/CV). Regions defined above (NSPLIT=5).
  // The outer-boundary shell (cells within ~2.5 dx of the cubic domain face) takes priority over the
  // radial near/far split, so it isolates ghost-fill/BC error from bulk-exterior error.
  double Hmax_split_local[NSPLIT] = { 0.0 };
  double H2_split_local[NSPLIT] = { 0.0 };
  double Hvol_split_local[NSPLIT] = { 0.0 };
  double Cgmax_split_local[NSPLIT] = { 0.0 };
  double Camax_split_local[NSPLIT] = { 0.0 };
  double Cvmax_split_local[NSPLIT] = { 0.0 };
  // [PROJ] H recomputed with the reduction variables projected to their definitions (ADM H).
  bool hproj_diag = (getenv("GKYL_HPROJ") != NULL);
  double Hproj2_split_local[NSPLIT] = { 0.0 };
  double Hprojmax_split_local[NSPLIT] = { 0.0 };
  // [FEC] Boundary shell sub-split by codimension: how many coords are in the boundary layer.
  // 1 = face, 2 = edge, 3 = corner. The falloff BC mirrors across ONE face, so an edge/corner ghost
  // reads a skin cell that is itself a ghost (ordering-dependent, inconsistently filled). If the
  // boundary badness concentrates in edge/corner, that is the mechanism.
  double Hfec2_local[3] = { 0.0 }, Hfecvol_local[3] = { 0.0 }, Hfecmax_local[3] = { 0.0 };

  // [CONS-RAD] Radial profile of the constraint violation. This is the diagnostic that
  // separates the candidate causes: an ORIGIN/Cartesian-n_i representation problem loads
  // the innermost shell; a STAR-SURFACE problem loads the shell at R_star (~9.6); a
  // frozen-WB-vs-evolving-geometry (matter) problem loads the matter interior (E != 0)
  // and stops at the surface; a genuine formulation/ADM problem is spread everywhere
  // INCLUDING the vacuum exterior. Also bin the matter density so "where E lives" is
  // directly comparable to "where H lives".
  #define NRAD 8
  double rad_hi[NRAD] = { 2.0, 4.0, 6.0, 8.0, 9.6, 11.0, 14.0, 1.0e30 };
  double radH2[NRAD] = { 0.0 }, radVol[NRAD] = { 0.0 }, radE[NRAD] = { 0.0 };
  double radDH2[NRAD] = { 0.0 }, radDHVol[NRAD] = { 0.0 }, radDHMax[NRAD] = { 0.0 };
  double radDTres2[NRAD] = { 0.0 }, radDTresVol[NRAD] = { 0.0 }, radDTresMax[NRAD] = { 0.0 };

  struct gkyl_range_iter iter;
  double density_surface_frac = 0.0, density_surface_r = 0.0, density_surface_Dmax = 0.0;
  double density_surface_abs_floor = 0.0, density_surface_threshold = 0.0;
  int density_surface_dynamic = 0;
  // Density-based interior region (polytrope convergence): r_int = radius where the spherically-binned
  // rest-mass density D drops to frac*D_max (~ rho > frac*rho_c, frac=GKYL_RHO_CORE_FRAC default 0.1).
  // Self-scaling with the star and tracks the compressing core in collapse. Overridden by an explicit
  // GKYL_CONSTRAINT_RINT. (Physical NS core = rho>rho_nuc is reserved for the hybrid EOS, not this toy polytrope.)
  if (getenv("GKYL_CONSTRAINT_RINT") == NULL) {
    double frac = 0.1;
    const char *fe = getenv("GKYL_RHO_CORE_FRAC"); if (fe) frac = atof(fe);
    const char *sf = getenv("GKYL_SURFACE_D_ABS_FLOOR"); if (sf) density_surface_abs_floor = atof(sf);
    if (density_surface_abs_floor < 0.0) density_surface_abs_floor = 0.0;
    density_surface_frac = frac;
    density_surface_dynamic = 1;
    int NBr = 200; double rmx = app->grid.upper[0]; double dbr = rmx / NBr;
    double *Dbin = calloc(NBr, sizeof(double));
    double Dmax_mon = 0.0;
    struct gkyl_range_iter it0;
    gkyl_range_iter_init(&it0, &app->local);
    while (gkyl_range_iter_next(&it0)) {
      long loc0 = gkyl_range_idx(&app->local_ext, it0.idx);
      const double *q0 = gkyl_array_cfetch(einstein, loc0);
      if (diag_spacetime_in_excision_region(ctx->excision_threshold, q0)) {
        continue;
      }
      double D0 = ((const double *) gkyl_array_cfetch(fluid, loc0))[0];
      double x0[3]; gkyl_rect_grid_cell_center(&app->grid, it0.idx, x0);
      double r0 = sqrt(x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]);
      if (D0 > Dmax_mon) Dmax_mon = D0;
      int b0 = (int)(r0 / dbr);
      if (b0 >= 0 && b0 < NBr && D0 > Dbin[b0]) Dbin[b0] = D0;
    }
    double r_rho = 0.0;
    density_surface_threshold = fmax(frac * Dmax_mon, density_surface_abs_floor);
    for (int b0 = 0; b0 < NBr; b0++) if (Dbin[b0] > density_surface_threshold) r_rho = (b0 + 0.5) * dbr;
    if (r_rho > 0.0) r_int = r_rho;
    density_surface_r = r_rho;
    density_surface_Dmax = Dmax_mon;
    free(Dbin);
  }

  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    const double *q = gkyl_array_cfetch(einstein, loc);
    const double *qf = gkyl_array_cfetch(fluid, loc);

    double xc[3];
    gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double r = sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);

    if (diag_constraint_stencil_touches_excision(app, einstein, iter.idx,
        ctx->excision_threshold)) {
      continue;
    }

    double g[3][3], gi[3][3], K[3][3], D[3][3][3], A[3], V[3];
    diag_read_einstein_state(q, g, gi, K, D, A, V);

    // [DETGT] Raw rescaled metric in slots 0-8, before the chi^-2 physical scale.
    // det via the rule of Sarrus. See the declaration above for interpretation.
    {
      double gt00 = q[0], gt01 = q[1], gt02 = q[2];
      double gt10 = q[3], gt11 = q[4], gt12 = q[5];
      double gt20 = q[6], gt21 = q[7], gt22 = q[8];
      double detgt = gt00 * (gt11 * gt22 - gt12 * gt21)
                   - gt01 * (gt10 * gt22 - gt12 * gt20)
                   + gt02 * (gt10 * gt21 - gt11 * gt20);
      if (isfinite(detgt)) {
        double dev = fabs(detgt - 1.0);
        if (dev > detgt_maxdev_local) {
          detgt_maxdev_local = dev;
          detgt_at_maxdev_local = detgt;
          detgt_r_local = r;
        }
      }
    }

    bool bad_cell = false;
    if (!isfinite(q[9])) {
      bad_cell = true;
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (!isfinite(g[i][j]) || !isfinite(K[i][j])) {
          bad_cell = true;
        }
      }
    }
    if (bad_cell) {
      n_bad_local += 1;
      continue;
    }

    double Cg = 0.0, Ca = 0.0;
    for (int k = 0; k < 3; k++) {
      int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
      for (int d = 0; d < 3; d++) {
        idxp[d] = iter.idx[d];
        idxm[d] = iter.idx[d];
      }
      idxp[k] = iter.idx[k] + 1;
      idxm[k] = iter.idx[k] - 1;
      const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
      const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
      double dx = app->grid.dx[k];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          double gp[3][3], gip[3][3], gm[3][3], gim[3][3];
          diag_read_physical_metric(qp, gp, gip);
          diag_read_physical_metric(qm, gm, gim);
          double half_dg = (gp[i][j] - gm[i][j]) / (4.0 * dx);
          double ckij = D[k][i][j] - half_dg;
          if (fabs(ckij) > Cg) {
            Cg = fabs(ckij);
          }
        }
      }

      double cp = qp[9], cm = qm[9];
      if (cp > 0.0 && cm > 0.0) {
        double dlog_alpha = (log(cp) - log(cm)) / (2.0 * dx);
        double ci = A[k] - dlog_alpha;
        if (fabs(ci) > Ca) {
          Ca = fabs(ci);
        }
      }
      else {
        n_bad_local += 1;
      }
    }

    double Cv = diag_aux_constraint(q, gi, D, V);

    double Gamma[3][3][3];
    diag_christoffel_from_state(q, Gamma);

    double Ricci[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Ricci[i][j] = 0.0;
        for (int k = 0; k < 3; k++) {
          int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
          for (int d = 0; d < 3; d++) {
            idxp[d] = iter.idx[d];
            idxm[d] = iter.idx[d];
          }
          idxp[k] = iter.idx[k] + 1;
          idxm[k] = iter.idx[k] - 1;
          const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
          const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
          double Gp[3][3][3], Gm[3][3][3];
          diag_christoffel_from_state(qp, Gp);
          diag_christoffel_from_state(qm, Gm);
          Ricci[i][j] += (Gp[k][i][j] - Gm[k][i][j]) / (2.0 * app->grid.dx[k]);
        }
        for (int k = 0; k < 3; k++) {
          int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
          for (int d = 0; d < 3; d++) {
            idxp[d] = iter.idx[d];
            idxm[d] = iter.idx[d];
          }
          idxp[j] = iter.idx[j] + 1;
          idxm[j] = iter.idx[j] - 1;
          const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
          const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
          double Gp[3][3][3], Gm[3][3][3];
          diag_christoffel_from_state(qp, Gp);
          diag_christoffel_from_state(qm, Gm);
          Ricci[i][j] -= (Gp[k][i][k] - Gm[k][i][k]) / (2.0 * app->grid.dx[j]);
        }
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            Ricci[i][j] += Gamma[k][i][j] * Gamma[l][k][l];
            Ricci[i][j] -= Gamma[l][i][k] * Gamma[k][j][l];
          }
        }
      }
    }

    double R = 0.0, trK = 0.0, KijKij = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        R += gi[i][j] * Ricci[i][j];
        trK += gi[i][j] * K[i][j];
      }
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double Kup = 0.0;
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            Kup += gi[i][a] * gi[j][b] * K[a][b];
          }
        }
        KijKij += K[i][j] * Kup;
      }
    }

    double E = 0.0, S_cov[3] = { 0.0, 0.0, 0.0 };
    diag_matter_projections(ctx->gas_gamma, qf, g, &E, S_cov);
    double H = R + trK * trK - KijKij - 16.0 * M_PI * E;

    // [PROJ] ADM Hamiltonian: same H but with R rebuilt from D_proj = 1/2 d_k gamma_ij (projected
    // reduction). If H_proj is small/convergent where H is large/non-convergent, the badness is
    // reduction drift (D_kij off its manifold), i.e. a reduction/BC problem, not real curvature.
    double H_proj = 0.0;
    if (hproj_diag) {
      double Gproj_c[3][3][3];
      diag_christoffel_proj(einstein, &app->local_ext, iter.idx, app->grid.dx, Gproj_c);
      double Rproj = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          double Ric = 0.0;
          for (int k = 0; k < 3; k++) {
            int ip[GKYL_MAX_DIM], im[GKYL_MAX_DIM];
            for (int d = 0; d < 3; d++) { ip[d] = iter.idx[d]; im[d] = iter.idx[d]; }
            ip[k] = iter.idx[k] + 1; im[k] = iter.idx[k] - 1;
            double Gp[3][3][3], Gm[3][3][3];
            diag_christoffel_proj(einstein, &app->local_ext, ip, app->grid.dx, Gp);
            diag_christoffel_proj(einstein, &app->local_ext, im, app->grid.dx, Gm);
            Ric += (Gp[k][i][j] - Gm[k][i][j]) / (2.0 * app->grid.dx[k]);
            int jp[GKYL_MAX_DIM], jm[GKYL_MAX_DIM];
            for (int d = 0; d < 3; d++) { jp[d] = iter.idx[d]; jm[d] = iter.idx[d]; }
            jp[j] = iter.idx[j] + 1; jm[j] = iter.idx[j] - 1;
            double Gpj[3][3][3], Gmj[3][3][3];
            diag_christoffel_proj(einstein, &app->local_ext, jp, app->grid.dx, Gpj);
            diag_christoffel_proj(einstein, &app->local_ext, jm, app->grid.dx, Gmj);
            Ric -= (Gpj[k][i][k] - Gmj[k][i][k]) / (2.0 * app->grid.dx[j]);
          }
          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              Ric += Gproj_c[k][i][j] * Gproj_c[l][k][l];
              Ric -= Gproj_c[l][i][k] * Gproj_c[k][j][l];
            }
          }
          Rproj += gi[i][j] * Ric;
        }
      }
      H_proj = Rproj + trK * trK - KijKij - 16.0 * M_PI * E;
    }
    // Normalized Hamiltonian: |H| divided by the sum of the magnitudes of its
    // own terms. This says whether H is a genuine O(1) imbalance or a small
    // relative residual (the standard "is the constraint violation actually big"
    // measure). And an interior L2 that EXCLUDES the innermost shell (r < 1.5 dx)
    // to separate a real distributed violation from a single center-cell FD spike.
    double Hscale = fabs(R) + trK*trK + KijKij + 16.0*M_PI*fabs(E) + 1.0e-30;
    double Hnorm = fabs(H)/Hscale;
    double dxmin = app->grid.dx[0];
    double surface_outer = ctx->R_star + 4.0 * dxmin;
    double box_half = app->grid.upper[0];
    double cheb = fmax(fabs(xc[0]), fmax(fabs(xc[1]), fabs(xc[2]))); // Chebyshev (cube face) distance
    double ext_mid = 0.5 * (surface_outer + box_half); // near/far exterior split
    // 0 core | 1 surface | 4 outer-boundary shell (priority) | 2 near-ext | 3 far-ext
    int split_region;
    if (r <= 0.8 * ctx->R_star) {
      split_region = 0;
    }
    else if (r <= surface_outer) {
      split_region = 1;
    }
    else if (cheb >= box_half - 2.5 * dxmin) {
      split_region = 4;
    }
    else if (r <= ext_mid) {
      split_region = 2;
    }
    else {
      split_region = 3;
    }
    double cellvol = app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2];
    int rad_bin = -1;
    for (int b = 0; b < NRAD; ++b) {
      if (r <= rad_hi[b]) {
        rad_bin = b;
        break;
      }
    }
    // [CAUSAL] region the boundary cannot have reached yet (exclude the star interior kink r<0.8R so
    // this measures the exterior scheme, not the surface; it is the "is only the boundary broken" test).
    if (r < r_causal && r > 0.8 * ctx->R_star) {
      H2_causal_local += H * H * cellvol;
      vol_causal_local += cellvol;
      if (fabs(H) > Hmax_causal_local) {
        Hmax_causal_local = fabs(H);
      }
    }
    H2_split_local[split_region] += H * H * cellvol;
    Hvol_split_local[split_region] += cellvol;
    if (fabs(H) > Hmax_split_local[split_region]) {
      Hmax_split_local[split_region] = fabs(H);
    }
    if (delta_diag) {
      double *h0 = gkyl_array_fetch(delta_H0, loc);
      if (!delta_seeded) {
        h0[0] = H;
      }
      double dH = H - h0[0];
      dH2_split_local[split_region] += dH * dH * cellvol;
      dHvol_split_local[split_region] += cellvol;
      if (fabs(dH) > dHmax_split_local[split_region]) {
        dHmax_split_local[split_region] = fabs(dH);
      }
      if (rad_bin >= 0) {
        radDH2[rad_bin] += dH * dH * cellvol;
        radDHVol[rad_bin] += cellvol;
        if (fabs(dH) > radDHMax[rad_bin]) {
          radDHMax[rad_bin] = fabs(dH);
        }
      }
      if (r <= r_int) {
        dH2_local += dH * dH * cellvol;
        dHvol_local += cellvol;
        if (fabs(dH) > dHmax_int_local) {
          dHmax_int_local = fabs(dH);
        }
      }
    }
    // Reduction constraints in the SAME bins, so bad H can be co-located with bad Cgamma/CA/CV.
    if (Cg > Cgmax_split_local[split_region]) {
      Cgmax_split_local[split_region] = Cg;
    }
    if (Ca > Camax_split_local[split_region]) {
      Camax_split_local[split_region] = Ca;
    }
    if (Cv > Cvmax_split_local[split_region]) {
      Cvmax_split_local[split_region] = Cv;
    }
    if (hproj_diag) {
      Hproj2_split_local[split_region] += H_proj * H_proj * cellvol;
      if (fabs(H_proj) > Hprojmax_split_local[split_region]) {
        Hprojmax_split_local[split_region] = fabs(H_proj);
      }
    }
    // [FEC] face/edge/corner sub-split of the boundary shell (uses H_proj if available, else H).
    if (split_region == 4) {
      int nbnd = 0;
      if (fabs(xc[0]) >= box_half - 2.5 * dxmin) {
        nbnd++;
      }
      if (fabs(xc[1]) >= box_half - 2.5 * dxmin) {
        nbnd++;
      }
      if (fabs(xc[2]) >= box_half - 2.5 * dxmin) {
        nbnd++;
      }
      if (nbnd >= 1 && nbnd <= 3) {
        int fec = nbnd - 1;
        double hval = hproj_diag ? H_proj : H;
        Hfec2_local[fec] += hval * hval * cellvol;
        Hfecvol_local[fec] += cellvol;
        if (fabs(hval) > Hfecmax_local[fec]) {
          Hfecmax_local[fec] = fabs(hval);
        }
      }
    }

    { // radial binning (see NRAD comment above)
      if (rad_bin >= 0) {
        radH2[rad_bin] += H*H*cellvol;
        radVol[rad_bin] += cellvol;
        radE[rad_bin] += fabs(E)*cellvol;
      }
    }

    if (r <= r_int) {
      if (Hnorm > Hnorm_int_local) Hnorm_int_local = Hnorm;
      if (r > 1.5*dxmin) {
        H2_noctr_local += H*H*cellvol;
        vol_noctr_local += cellvol;
      }
    }

    double Mmax_cell = 0.0;
    for (int i = 0; i < 3; i++) {
      double divQ = 0.0;
      for (int j = 0; j < 3; j++) {
        int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
        for (int d = 0; d < 3; d++) {
          idxp[d] = iter.idx[d];
          idxm[d] = iter.idx[d];
        }
        idxp[j] = iter.idx[j] + 1;
        idxm[j] = iter.idx[j] - 1;
        const double *qp = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxp));
        const double *qm = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idxm));
        double Qp[3][3], Qm[3][3];
        diag_compute_Q(qp, Qp);
        diag_compute_Q(qm, Qm);
        divQ += (Qp[j][i] - Qm[j][i]) / (2.0 * app->grid.dx[j]);
      }

      double Q[3][3];
      diag_compute_Q(q, Q);
      for (int j = 0; j < 3; j++) {
        for (int l = 0; l < 3; l++) {
          divQ += Gamma[j][j][l] * Q[l][i];
          divQ -= Gamma[l][j][i] * Q[j][l];
        }
      }

      double Mi = divQ - 8.0 * M_PI * S_cov[i];
      if (fabs(Mi) > Mmax_cell) {
        Mmax_cell = fabs(Mi);
      }
    }

    if (fabs(H) > Hmax_local) {
      Hmax_local = fabs(H);
      r_H_local = r;
      R_at_H_local = R;
      trK2_at_H_local = trK * trK;
      K2_at_H_local = KijKij;
      E_at_H_local = E;
    }
    if (Mmax_cell > Mmax_local) {
      Mmax_local = Mmax_cell;
      r_M_local = r;
    }
    if (r <= r_int) {
      if (fabs(H) > Hint_local) Hint_local = fabs(H);
      if (Mmax_cell > Mint_local) Mint_local = Mmax_cell;
      if (Cg > Cgint_local) Cgint_local = Cg;
      H2_local += H*H*cellvol;
      vol_local += cellvol;
    }
    if (Cg > Cgmax_local) Cgmax_local = Cg;
    if (Ca > Camax_local) Camax_local = Ca;
    // Pick the single cell nearest the +x axis on each side (y, z in the first positive cell),
    // otherwise the four cells straddling the axis are degenerate by symmetry and the sequence
    // becomes a staircase, which fakes a zig-zag signal.
    // [TRES] Covariant conservation residual of the MATTER, which is what Bianchi needs for the
    // constraints to propagate:
    //   nabla_mu T^mu_j = (1/sqrt(-g)) d_mu (sqrt(-g) T^mu_j) - (1/2) T^{ab} d_j g_{ab} = 0.
    // The fluid is held with S_j == 0 (well-balanced / frozen), so T^t_j == 0 and the time
    // derivative drops out -- the whole residual is computable from a single time level:
    //   R_j = (1/(alpha sqrt(gamma))) d_i(alpha sqrt(gamma) T^i_j) + alpha^2 A_j T^{00} - T^{kl} D_jkl
    // (using beta = 0, so d_j g_00 = -2 alpha^2 A_j, d_j g_0i = 0, d_j g_kl = 2 D_jkl).
    // At t=0 this vanishes to truncation error -- that IS the TOV equilibrium. If it fails to
    // converge with h at later times, the frozen matter is NOT a solution in the evolved metric,
    // i.e. nabla_mu T^{mu nu} != 0 and Bianchi no longer protects H. dSumD cannot see this.
    if (tres_diag) {
      double sgT_p[3][3], sgT_m[3][3];
      double divT[3] = { 0.0, 0.0, 0.0 };
      for (int i = 0; i < 3; i++) {
        int idxp[GKYL_MAX_DIM], idxm[GKYL_MAX_DIM];
        for (int d = 0; d < 3; d++) {
          idxp[d] = iter.idx[d];
          idxm[d] = iter.idx[d];
        }
        idxp[i] = iter.idx[i] + 1;
        idxm[i] = iter.idx[i] - 1;
        long lp = gkyl_range_idx(&app->local_ext, idxp);
        long lm = gkyl_range_idx(&app->local_ext, idxm);
        diag_sqrtg_Tmixed(ctx->gas_gamma, gkyl_array_cfetch(fluid, lp), gkyl_array_cfetch(einstein, lp), sgT_p);
        diag_sqrtg_Tmixed(ctx->gas_gamma, gkyl_array_cfetch(fluid, lm), gkyl_array_cfetch(einstein, lm), sgT_m);
        for (int j = 0; j < 3; j++) {
          divT[j] += (sgT_p[i][j] - sgT_m[i][j]) / (2.0 * app->grid.dx[i]);
        }
      }

      double **Tc = gkyl_malloc(sizeof(double*[4]));
      for (int a = 0; a < 4; a++) {
        Tc[a] = gkyl_malloc(sizeof(double[4]));
      }
      gkyl_gr_euler_stress_energy_tensor(ctx->gas_gamma, qf, &Tc);

      double gg[3][3], ggi[3][3];
      double detc = 0.0;
      diag_read_physical_metric(q, gg, ggi);
      detc = diag_sym3_inverse(gg, ggi);
      double sdetc = 0.0;
      if (detc > 0.0) {
        sdetc = sqrt(detc);
      }
      double lapsec = q[9];
      double sqrtmg = lapsec * sdetc;

      // Full -(1/2) T^{ab} d_j g_{ab} with beta != 0 (Kerr-Schild): the beta=0 form (alpha^2 A_j T^00 -
      // T^{kl}D_jkl) OMITS the shift terms d_j g_0i = d_j beta_i (times momentum T^{0i}) and the beta.beta
      // part of d_j g_00. Those are identically zero in the areal (beta=0) gauge but O(1) in KS, which is
      // why div-T looked non-convergent (~4e-5 flat) there -- a diagnostic incompleteness, not the scheme.
      double beta[3] = { q[52], q[53], q[54] };
      double dbeta[3][3]; // dbeta[j][k] = d_j beta^k. Slots 55-63 store HALF the derivative (the reduction
      // convention, like D_kij = 1/2 d_k gamma_ij), so multiply by 2 to recover d_j beta^k.
      for (int j2 = 0; j2 < 3; j2++) for (int k2 = 0; k2 < 3; k2++) dbeta[j2][k2] = 2.0 * q[55 + 3*j2 + k2];
      double beta_low[3] = { 0.0, 0.0, 0.0 };
      for (int i2 = 0; i2 < 3; i2++) for (int k2 = 0; k2 < 3; k2++) beta_low[i2] += gg[i2][k2] * beta[k2];
      double dbeta_low[3][3], dbeta2[3]; // d_j beta_i = 2 D_jik beta^k + gamma_ik d_j beta^k ; d_j(beta_m beta^m)
      for (int j2 = 0; j2 < 3; j2++) {
        double db2 = 0.0;
        for (int i2 = 0; i2 < 3; i2++) {
          double s2 = 0.0;
          for (int k2 = 0; k2 < 3; k2++) s2 += 2.0 * D[j2][i2][k2] * beta[k2] + gg[i2][k2] * dbeta[j2][k2];
          dbeta_low[j2][i2] = s2;
        }
        for (int m2 = 0; m2 < 3; m2++) {
          for (int n2 = 0; n2 < 3; n2++) db2 += 2.0 * D[j2][m2][n2] * beta[m2] * beta[n2];
          db2 += 2.0 * beta_low[m2] * dbeta[j2][m2];
        }
        dbeta2[j2] = db2;
      }

      double *t0 = delta_diag ? gkyl_array_fetch(delta_T0, loc) : NULL;
      for (int j = 0; j < 3; j++) {
        double src = (lapsec * lapsec) * A[j] * Tc[0][0];
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            src -= Tc[k + 1][l + 1] * D[j][k][l];
          }
        }
        // beta != 0 additions: -(1/2) T^{00} d_j(beta.beta) - T^{0i} d_j beta_i.
        src -= 0.5 * Tc[0][0] * dbeta2[j];
        for (int i2 = 0; i2 < 3; i2++) src -= Tc[0][i2 + 1] * dbeta_low[j][i2];
        double Rj = src;
        if (sqrtmg > 0.0) {
          Rj += divT[j] / sqrtmg;
        }
        if (fabs(Rj) > tresmax_local) {
          tresmax_local = fabs(Rj);
        }
        if (r <= r_int) {
          tres2_local += Rj * Rj * cellvol;
        }
        tres2_split_local[split_region] += Rj * Rj * cellvol;
        if (fabs(Rj) > tresmax_split_local[split_region]) {
          tresmax_split_local[split_region] = fabs(Rj);
        }
        if (delta_diag) {
          if (!delta_seeded) {
            t0[j] = Rj;
          }
          double dRj = Rj - t0[j];
          dTres2_split_local[split_region] += dRj * dRj * cellvol;
          if (fabs(dRj) > dTresmax_split_local[split_region]) {
            dTresmax_split_local[split_region] = fabs(dRj);
          }
          if (rad_bin >= 0) {
            radDTres2[rad_bin] += dRj * dRj * cellvol;
            if (fabs(dRj) > radDTresMax[rad_bin]) {
              radDTresMax[rad_bin] = fabs(dRj);
            }
          }
          if (r <= r_int) {
            dTres2_local += dRj * dRj * cellvol;
            if (fabs(dRj) > dTresmax_int_local) {
              dTresmax_int_local = fabs(dRj);
            }
          }
        }
      }
      if (r <= r_int) {
        tresvol_local += cellvol;
        if (delta_diag) {
          dTresvol_local += cellvol;
        }
      }
      tresvol_split_local[split_region] += cellvol;
      if (delta_diag) {
        dTresvol_split_local[split_region] += cellvol;
        if (rad_bin >= 0) {
          radDTresVol[rad_bin] += cellvol;
        }
      }

      for (int a = 0; a < 4; a++) {
        gkyl_free(Tc[a]);
      }
      gkyl_free(Tc);
    }

    if (line_diag && nline < NLINE
      && xc[1] > 0.0 && xc[1] < app->grid.dx[1] && xc[2] > 0.0 && xc[2] < app->grid.dx[2] && xc[0] > 0.0) {
      lineX[nline] = xc[0];
      lineH[nline] = H;
      nline += 1;
    }

    if (Cv > Cvmax_local) Cvmax_local = Cv;
    if (fabs(trK) > fabs(trKmax_local)) trKmax_local = trK;
    if (r < r_c_local) {
      r_c_local = r;
      trKc_local = trK;
      H_at_c_local = H;
      R_at_c_local = R;
      trK2_at_c_local = trK * trK;
      K2_at_c_local = KijKij;
      E_at_c_local = E;
    }

    if (!isfinite(H) || !isfinite(Mmax_cell) || !isfinite(Cg) || !isfinite(Ca) || !isfinite(Cv)) {
      n_bad_local += 1;
    }
  }

  double Hmax = 0.0, Mmax = 0.0, Cgmax = 0.0, Camax = 0.0, Cvmax = 0.0, abs_trKmax = 0.0, r_c = 0.0;
  int n_bad = 0;
  double abs_trKmax_local = fabs(trKmax_local);
  double Hint = 0.0, Mint = 0.0, Cgint = 0.0, H2sum = 0.0, volsum = 0.0;
  double Hnorm_int = 0.0, H2noctr = 0.0, volnoctr = 0.0;
  double Hmax_split[NSPLIT] = { 0.0 };
  double H2_split[NSPLIT] = { 0.0 };
  double Hvol_split[NSPLIT] = { 0.0 };
  double Cgmax_split[NSPLIT] = { 0.0 };
  double Camax_split[NSPLIT] = { 0.0 };
  double Cvmax_split[NSPLIT] = { 0.0 };
  double dHmax_int = 0.0, dH2 = 0.0, dHvol = 0.0;
  double dHmax_split[NSPLIT] = { 0.0 };
  double dH2_split[NSPLIT] = { 0.0 };
  double dHvol_split[NSPLIT] = { 0.0 };
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hnorm_int_local, &Hnorm_int);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H2_noctr_local, &H2noctr);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &vol_noctr_local, &volnoctr);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, Hmax_split_local, Hmax_split);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, H2_split_local, H2_split);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, Hvol_split_local, Hvol_split);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, Cgmax_split_local, Cgmax_split);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, Camax_split_local, Camax_split);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, Cvmax_split_local, Cvmax_split);
  if (delta_diag) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &dHmax_int_local, &dHmax_int);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &dH2_local, &dH2);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &dHvol_local, &dHvol);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, dHmax_split_local, dHmax_split);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, dH2_split_local, dH2_split);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, dHvol_split_local, dHvol_split);
  }
  double Hproj2_split[NSPLIT] = { 0.0 };
  double Hprojmax_split[NSPLIT] = { 0.0 };
  if (hproj_diag) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, Hproj2_split_local, Hproj2_split);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, Hprojmax_split_local, Hprojmax_split);
  }
  double Hfec2[3] = { 0.0 }, Hfecvol[3] = { 0.0 }, Hfecmax[3] = { 0.0 };
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 3, Hfec2_local, Hfec2);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 3, Hfecvol_local, Hfecvol);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 3, Hfecmax_local, Hfecmax);
  double radDH2_all[NRAD] = { 0.0 }, radDHVol_all[NRAD] = { 0.0 }, radDHMax_all[NRAD] = { 0.0 };
  double radDTres2_all[NRAD] = { 0.0 }, radDTresVol_all[NRAD] = { 0.0 }, radDTresMax_all[NRAD] = { 0.0 };
  if (delta_diag) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NRAD, radDH2, radDH2_all);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NRAD, radDHVol, radDHVol_all);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NRAD, radDHMax, radDHMax_all);
    if (tres_diag) {
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NRAD, radDTres2, radDTres2_all);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NRAD, radDTresVol, radDTresVol_all);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NRAD, radDTresMax, radDTresMax_all);
    }
  }
  double H2_causal = 0.0, vol_causal = 0.0, Hmax_causal = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H2_causal_local, &H2_causal);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &vol_causal_local, &vol_causal);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hmax_causal_local, &Hmax_causal);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H2_local, &H2sum);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &vol_local, &volsum);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hint_local, &Hint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Mint_local, &Mint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cgint_local, &Cgint);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Hmax_local, &Hmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Mmax_local, &Mmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cgmax_local, &Cgmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Camax_local, &Camax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &Cvmax_local, &Cvmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &abs_trKmax_local, &abs_trKmax);
  double detgt_maxdev = detgt_maxdev_local;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &detgt_maxdev_local, &detgt_maxdev);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_c_local, &r_c);
  gkyl_comm_allreduce(app->comm, GKYL_INT, GKYL_SUM, 1, &n_bad_local, &n_bad);

  double r_H_candidate = 1.0e300;
  if (Hmax_local == Hmax) r_H_candidate = r_H_local;
  double r_M_candidate = 1.0e300;
  if (Mmax_local == Mmax) r_M_candidate = r_M_local;
  double trKmax_candidate = -1.0e300;
  if (fabs(trKmax_local) == abs_trKmax) trKmax_candidate = trKmax_local;
  double trKc_candidate = -1.0e300;
  if (r_c_local == r_c) trKc_candidate = trKc_local;
  double R_at_H_candidate = 0.0, trK2_at_H_candidate = 0.0, K2_at_H_candidate = 0.0, E_at_H_candidate = 0.0;
  if (Hmax_local == Hmax) {
    R_at_H_candidate = R_at_H_local;
    trK2_at_H_candidate = trK2_at_H_local;
    K2_at_H_candidate = K2_at_H_local;
    E_at_H_candidate = E_at_H_local;
  }
  double H_at_c_candidate = 0.0, R_at_c_candidate = 0.0, trK2_at_c_candidate = 0.0, K2_at_c_candidate = 0.0, E_at_c_candidate = 0.0;
  if (r_c_local == r_c) {
    H_at_c_candidate = H_at_c_local;
    R_at_c_candidate = R_at_c_local;
    trK2_at_c_candidate = trK2_at_c_local;
    K2_at_c_candidate = K2_at_c_local;
    E_at_c_candidate = E_at_c_local;
  }

  double r_H = 0.0, r_M = 0.0, trKmax = 0.0, trKc = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_H_candidate, &r_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &r_M_candidate, &r_M);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &trKmax_candidate, &trKmax);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &trKc_candidate, &trKc);
  double R_at_H = 0.0, trK2_at_H = 0.0, K2_at_H = 0.0, E_at_H = 0.0;
  double H_at_c = 0.0, R_at_c = 0.0, trK2_at_c = 0.0, K2_at_c = 0.0, E_at_c = 0.0;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &R_at_H_candidate, &R_at_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &trK2_at_H_candidate, &trK2_at_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &K2_at_H_candidate, &K2_at_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &E_at_H_candidate, &E_at_H);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &H_at_c_candidate, &H_at_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &R_at_c_candidate, &R_at_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &trK2_at_c_candidate, &trK2_at_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &K2_at_c_candidate, &K2_at_c);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &E_at_c_candidate, &E_at_c);

  double tresmax = 0.0, tres2 = 0.0, tresvol = 0.0;
  double tresmax_split[NSPLIT] = { 0.0 };
  double tres2_split[NSPLIT] = { 0.0 };
  double tresvol_split[NSPLIT] = { 0.0 };
  double dTresmax_int = 0.0, dTres2 = 0.0, dTresvol = 0.0;
  double dTresmax_split[NSPLIT] = { 0.0 };
  double dTres2_split[NSPLIT] = { 0.0 };
  double dTresvol_split[NSPLIT] = { 0.0 };
  if (tres_diag) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &tresmax_local, &tresmax);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &tres2_local, &tres2);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &tresvol_local, &tresvol);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, tresmax_split_local, tresmax_split);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, tres2_split_local, tres2_split);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, tresvol_split_local, tresvol_split);
    if (delta_diag) {
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &dTresmax_int_local, &dTresmax_int);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &dTres2_local, &dTres2);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &dTresvol_local, &dTresvol);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, NSPLIT, dTresmax_split_local, dTresmax_split);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, dTres2_split_local, dTres2_split);
      gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, NSPLIT, dTresvol_split_local, dTresvol_split);
    }
  }

  int rank = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    if (density_surface_dynamic) {
      printf("[SURFACE] t=%.4f r_rho=%.4f frac=%.3g Dmax=%.6e Dfloor=%.6e Dthr=%.6e  (outer radius with shell max D > max(frac*Dmax,Dfloor))\n",
        t_curr, density_surface_r, density_surface_frac, density_surface_Dmax,
        density_surface_abs_floor, density_surface_threshold);
    }
    printf("[CONS] t=%.4f  max|H|=%.6e @r=%.3f  max|M|=%.6e @r=%.3f  max|Cgamma|=%.6e max|CA|=%.6e max|CV|=%.6e  trK_c=%.6e max|trK|=%.6e  bad=%d\n",
      t_curr, Hmax, r_H, Mmax, r_M, Cgmax, Camax, Cvmax, trKc, trKmax, n_bad);
    printf("[DETGT] t=%.4f  max|det(gamma-tilde)-1|=%.6e  det=%.6e @r=%.3f  (raw conformal-metric determinant)\n",
      t_curr, detgt_maxdev, detgt_at_maxdev_local, detgt_r_local);
    if (hterm_diag) {
      printf("[CONS-HTERM] t=%.4f  at max|H| r=%.3f: R=%+.16e trK2=%+.16e KijKij=%+.16e 16piE=%+.16e H=%+.16e\n",
        t_curr, r_H, R_at_H, trK2_at_H, K2_at_H, 16.0 * M_PI * E_at_H,
        R_at_H + trK2_at_H - K2_at_H - 16.0 * M_PI * E_at_H);
      printf("[CONS-HTERM] t=%.4f  at center r=%.3f: H=%+.16e R=%+.16e trK2=%+.16e KijKij=%+.16e 16piE=%+.16e\n",
        t_curr, r_c, H_at_c, R_at_c, trK2_at_c, K2_at_c, 16.0 * M_PI * E_at_c);
    }
    double H_l2 = (volsum > 0.0) ? sqrt(H2sum/volsum) : 0.0;
    double H_l2_noctr = (volnoctr > 0.0) ? sqrt(H2noctr/volnoctr) : 0.0;
    printf("[CONS-INT] t=%.4f  r<=%.1f: max|H|=%.6e L2|H|=%.6e L2|H|(no-ctr)=%.6e maxHnorm=%.4f max|M|=%.6e max|Cgamma|=%.6e\n",
      t_curr, r_int, Hint, H_l2, H_l2_noctr, Hnorm_int, Mint, Cgint);
    printf("[CONS-SPLIT] t=%.4f  core(r<=0.8R): max|H|=%.6e L2|H|=%.6e  surface(0.8R<r<=R+4dx): max|H|=%.6e L2|H|=%.6e  atm(r>R+4dx): max|H|=%.6e L2|H|=%.6e\n",
      t_curr,
      Hmax_split[0], (Hvol_split[0] > 0.0) ? sqrt(H2_split[0] / Hvol_split[0]) : 0.0,
      Hmax_split[1], (Hvol_split[1] > 0.0) ? sqrt(H2_split[1] / Hvol_split[1]) : 0.0,
      Hmax_split[2] > Hmax_split[3] ? Hmax_split[2] : Hmax_split[3],
      (Hvol_split[2] + Hvol_split[3] > 0.0) ? sqrt((H2_split[2] + H2_split[3]) / (Hvol_split[2] + Hvol_split[3])) : 0.0);
    // [REGION5] the finer split: near-ext / far-ext / outer-boundary, each with H AND the reduction
    // constraints, so a BC problem (bad at bnd, reductions bad) is separable from a bulk-exterior or
    // surface-matching one. L2|H| per region + max Cgamma/CA/CV per region.
    const char *rn[NSPLIT] = { "core", "surf", "nearE", "farE", "BND" };
    for (int rg = 0; rg < NSPLIT; rg++) {
      double l2 = (Hvol_split[rg] > 0.0) ? sqrt(H2_split[rg] / Hvol_split[rg]) : 0.0;
      printf("[CONS-REG] t=%.4f  %-5s: L2|H|=%.6e max|H|=%.6e  max|Cg|=%.6e max|CA|=%.6e max|CV|=%.6e\n",
        t_curr, rn[rg], l2, Hmax_split[rg], Cgmax_split[rg], Camax_split[rg], Cvmax_split[rg]);
    }
    if (hproj_diag) {
      for (int rg = 0; rg < NSPLIT; rg++) {
        double l2p = (Hvol_split[rg] > 0.0) ? sqrt(Hproj2_split[rg] / Hvol_split[rg]) : 0.0;
        printf("[CONS-HPROJ] t=%.4f  %-5s: L2|H_proj|=%.6e max|H_proj|=%.6e  (vs L2|H|=%.6e; small ⇒ badness is reduction drift)\n",
          t_curr, rn[rg], l2p, Hprojmax_split[rg], (Hvol_split[rg] > 0.0) ? sqrt(H2_split[rg] / Hvol_split[rg]) : 0.0);
      }
    }
    {
      double l2c = (vol_causal > 0.0) ? sqrt(H2_causal / vol_causal) : 0.0;
      printf("[CONS-CAUSAL] t=%.4f  r<%.2f (0.8R<r, boundary-disconnected): L2|H|=%.6e max|H|=%.6e ncells~%.0f\n",
        t_curr, r_causal, l2c, Hmax_causal,
        vol_causal / (app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2]));
    }
    const char *fecn[3] = { "face", "edge", "corner" };
    for (int fec = 0; fec < 3; fec++) {
      double l2 = (Hfecvol[fec] > 0.0) ? sqrt(Hfec2[fec] / Hfecvol[fec]) : 0.0;
      printf("[CONS-FEC] t=%.4f  BND-%-6s: L2=%.6e max=%.6e ncells~%.0f  (%s)\n",
        t_curr, fecn[fec], l2, Hfecmax[fec],
        (Hfecvol[fec] / (app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2])),
        hproj_diag ? "H_proj" : "H");
    }

    // [LINE] Signed H along +x, plus the discrete curvature ratio
    //   Q = sum (H_{i+1} - 2 H_i + H_{i-1})^2 / sum H_i^2,
    // which is O(h^4) for a smooth profile but O(1) for a 2*dx zig-zag.
    if (line_diag && nline > 2) {
      // insertion sort by x (range iteration order is not guaranteed monotonic)
      for (int i = 1; i < nline; i++) {
        double kx = lineX[i], kh = lineH[i];
        int j = i - 1;
        while (j >= 0 && lineX[j] > kx) {
          lineX[j + 1] = lineX[j];
          lineH[j + 1] = lineH[j];
          j -= 1;
        }
        lineX[j + 1] = kx;
        lineH[j + 1] = kh;
      }

      double num = 0.0, den = 0.0;
      for (int i = 1; i < nline - 1; i++) {
        double d2 = lineH[i + 1] - (2.0 * lineH[i]) + lineH[i - 1];
        num += d2 * d2;
        den += lineH[i] * lineH[i];
      }
      printf("[CONS-LINE] t=%.4f  n=%d  zigzag Q=%.4e  H(x):", t_curr, nline, (den > 0.0) ? num / den : 0.0);
      for (int i = 0; i < nline; i++) {
        printf(" %+.2e", lineH[i]);
      }
      printf("\n");
    }
    if (tres_diag) {
      double tl2 = 0.0;
      if (tresvol > 0.0) {
        tl2 = sqrt(tres2 / tresvol);
      }
      printf("[TRES] t=%.4f  r<=%.1f: max|div T|=%.6e  L2|div T|=%.6e\n", t_curr, r_int, tresmax, tl2);
      const char *trn[NSPLIT] = { "core", "surf", "nearE", "farE", "BND" };
      for (int rg = 0; rg < NSPLIT; rg++) {
        double l2 = (tresvol_split[rg] > 0.0) ? sqrt(tres2_split[rg] / tresvol_split[rg]) : 0.0;
        printf("[TRES-REG] t=%.4f  %-5s: max|div T|=%.6e L2|div T|=%.6e\n", t_curr, trn[rg], tresmax_split[rg], l2);
      }
    }
    if (delta_diag) {
      double dHl2 = (dHvol > 0.0) ? sqrt(dH2 / dHvol) : 0.0;
      printf("[CONS-DELTA] t=%.4f  r<=%.1f: max|dH|=%.6e L2|dH|=%.6e seeded=%d\n",
        t_curr, r_int, dHmax_int, dHl2, delta_seeded);
      const char *drn[NSPLIT] = { "core", "surf", "nearE", "farE", "BND" };
      for (int rg = 0; rg < NSPLIT; rg++) {
        double l2 = (dHvol_split[rg] > 0.0) ? sqrt(dH2_split[rg] / dHvol_split[rg]) : 0.0;
        printf("[CONS-DELTA-REG] t=%.4f  %-5s: max|dH|=%.6e L2|dH|=%.6e\n",
          t_curr, drn[rg], dHmax_split[rg], l2);
      }
      printf("[CONS-DELTA-RAD] t=%.4f  L2|dH| per shell:", t_curr);
      double rlo = 0.0;
      for (int b = 0; b < NRAD; ++b) {
        double dhb = (radDHVol_all[b] > 0.0) ? sqrt(radDH2_all[b] / radDHVol_all[b]) : 0.0;
        double rhi = (rad_hi[b] > 1.0e29) ? 99.0 : rad_hi[b];
        printf("  [%.0f-%.0f]=%.3e(max=%.1e)", rlo, rhi, dhb, radDHMax_all[b]);
        rlo = rhi;
      }
      printf("\n");
      if (tres_diag) {
        double dtl2 = (dTresvol > 0.0) ? sqrt(dTres2 / dTresvol) : 0.0;
        printf("[TRES-DELTA] t=%.4f  r<=%.1f: max|dDivT|=%.6e L2|dDivT|=%.6e\n",
          t_curr, r_int, dTresmax_int, dtl2);
        for (int rg = 0; rg < NSPLIT; rg++) {
          double l2 = (dTresvol_split[rg] > 0.0) ? sqrt(dTres2_split[rg] / dTresvol_split[rg]) : 0.0;
          printf("[TRES-DELTA-REG] t=%.4f  %-5s: max|dDivT|=%.6e L2|dDivT|=%.6e\n",
            t_curr, drn[rg], dTresmax_split[rg], l2);
        }
        printf("[TRES-DELTA-RAD] t=%.4f  L2|dDivT| per shell:", t_curr);
        rlo = 0.0;
        for (int b = 0; b < NRAD; ++b) {
          double dtb = (radDTresVol_all[b] > 0.0) ? sqrt(radDTres2_all[b] / radDTresVol_all[b]) : 0.0;
          double rhi = (rad_hi[b] > 1.0e29) ? 99.0 : rad_hi[b];
          printf("  [%.0f-%.0f]=%.3e(max=%.1e)", rlo, rhi, dtb, radDTresMax_all[b]);
          rlo = rhi;
        }
        printf("\n");
      }
    }
    printf("[CONS-RAD] t=%.4f  L2|H| per shell:", t_curr);
    double rlo = 0.0;
    for (int b = 0; b < NRAD; ++b) {
      double hb = (radVol[b] > 0.0) ? sqrt(radH2[b]/radVol[b]) : 0.0;
      double eb = (radVol[b] > 0.0) ? radE[b]/radVol[b] : 0.0;
      double rhi = (rad_hi[b] > 1.0e29) ? 99.0 : rad_hi[b];
      printf("  [%.0f-%.0f]=%.3e(E=%.1e)", rlo, rhi, hb, eb);
      rlo = rhi;
    }
    printf("\n");
    fflush(stdout);
  }
  if (need_delta_seed) {
    delta_seeded = 1;
  }
}

static bool
ahcert_interp_einstein_state(gkyl_moment_app *app, const struct gkyl_array *einstein,
  const double x[3], double qout[77])
{
  int idx0[3], idx1[3];
  double w1[3];

  for (int d = 0; d < 3; ++d) {
    double pos = (x[d] - app->grid.lower[d]) / app->grid.dx[d] + 0.5;
    if (!(isfinite(pos) && pos >= app->local.lower[d] && pos <= app->local.upper[d])) {
      return false;
    }

    int i0 = (int) floor(pos);
    double f = pos - i0;
    if (i0 < app->local.lower[d]) {
      i0 = app->local.lower[d];
      f = 0.0;
    }
    if (i0 >= app->local.upper[d]) {
      i0 = app->local.upper[d];
      f = 0.0;
    }
    idx0[d] = i0;
    idx1[d] = (i0 < app->local.upper[d]) ? i0 + 1 : i0;
    w1[d] = f;
  }

  int ncomp = diag_is_conformal ? 77 : 64;
  for (int m = 0; m < 77; ++m) {
    qout[m] = 0.0;
  }

  for (int bx = 0; bx <= 1; ++bx) {
    for (int by = 0; by <= 1; ++by) {
      for (int bz = 0; bz <= 1; ++bz) {
        int idx[3] = {
          bx ? idx1[0] : idx0[0],
          by ? idx1[1] : idx0[1],
          bz ? idx1[2] : idx0[2],
        };
        double wt = (bx ? w1[0] : 1.0 - w1[0])
          * (by ? w1[1] : 1.0 - w1[1])
          * (bz ? w1[2] : 1.0 - w1[2]);
        const double *q = gkyl_array_cfetch(einstein, gkyl_range_idx(&app->local_ext, idx));
        for (int m = 0; m < ncomp; ++m) {
          qout[m] += wt * q[m];
        }
      }
    }
  }

  return true;
}

static bool
ahcert_spherical_area_density(gkyl_moment_app *app, const struct gkyl_array *einstein,
  struct gr_tov_coupled_ctx *ctx, const double x[3], double R0, double theta, double phi,
  double *dA_out)
{
  double q[77];
  if (!ahcert_interp_einstein_state(app, einstein, x, q)) {
    return false;
  }
  if (diag_spacetime_in_excision_region(ctx->excision_threshold, q)) {
    return false;
  }

  double g[3][3], gi[3][3];
  diag_read_physical_metric(q, g, gi);

  double st = sin(theta), ct = cos(theta);
  double sp = sin(phi), cp = cos(phi);
  double e_theta[3] = { R0*ct*cp, R0*ct*sp, -R0*st };
  double e_phi[3] = { -R0*st*sp, R0*st*cp, 0.0 };
  double E = 0.0, F = 0.0, G = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      E += g[i][j] * e_theta[i] * e_theta[j];
      F += g[i][j] * e_theta[i] * e_phi[j];
      G += g[i][j] * e_phi[i] * e_phi[j];
    }
  }
  double area2 = E*G - F*F;
  if (!(isfinite(area2) && area2 > 0.0)) {
    return false;
  }

  *dA_out = sqrt(area2);
  return isfinite(*dA_out);
}

struct ahcert_scan_cfg {
  double Rmin, Rmax, theta_thresh, max_bad_frac;
  int NR, NT, NP;
};

struct ahcert_scan_result {
  int ran;
  int NR, NT, NP;
  int best_valid, trap_valid, bad, total;
  double best_R, best_max, best_min, best_R_areal;
  double trap_R, trap_max, trap_min, trap_R_areal;
  double Rmin, Rmax, theta_thresh, max_bad_frac;
  // Horizon location: zero of the sphere-averaged Theta_+ (outermost + -> <=0 crossing, linear in R).
  int root_valid;
  double root_R, root_R_areal;
};

static struct ahcert_scan_cfg
ahcert_scan_cfg_from_env(gkyl_moment_app *app, const char *prefix, double default_Rmax,
  int default_NR, int default_NT, int default_NP)
{
  double dxmax = fmax(fmax(app->grid.dx[0], app->grid.dx[1]), app->grid.dx[2]);
  struct ahcert_scan_cfg cfg = {
    .Rmin = 1.5*dxmax,
    .Rmax = default_Rmax,
    .theta_thresh = 0.0,
    .max_bad_frac = 0.0,
    .NR = default_NR,
    .NT = default_NT,
    .NP = default_NP,
  };
  if (!(cfg.Rmax > cfg.Rmin)) cfg.Rmax = 0.82*app->grid.upper[0];

  char key[128];
  const char *ev = NULL;
  snprintf(key, sizeof(key), "%s_RMIN", prefix);
  ev = getenv(key); if (ev) cfg.Rmin = atof(ev);
  snprintf(key, sizeof(key), "%s_RMAX", prefix);
  ev = getenv(key); if (ev) cfg.Rmax = atof(ev);
  if (!(cfg.Rmax > cfg.Rmin)) cfg.Rmax = cfg.Rmin + dxmax;

  snprintf(key, sizeof(key), "%s_NR", prefix);
  ev = getenv(key); if (ev) cfg.NR = atoi(ev);
  snprintf(key, sizeof(key), "%s_NTHETA", prefix);
  ev = getenv(key); if (ev) cfg.NT = atoi(ev);
  snprintf(key, sizeof(key), "%s_NPHI", prefix);
  ev = getenv(key); if (ev) cfg.NP = atoi(ev);
  if (cfg.NR < 4) cfg.NR = 4;
  if (cfg.NT < 5) cfg.NT = 5;
  if (cfg.NP < 4) cfg.NP = 4;

  snprintf(key, sizeof(key), "%s_THETA_THRESH", prefix);
  ev = getenv(key); if (ev) cfg.theta_thresh = atof(ev);

  snprintf(key, sizeof(key), "%s_MAX_BAD_FRAC", prefix);
  ev = getenv(key);
  if (ev) {
    cfg.max_bad_frac = atof(ev);
    if (cfg.max_bad_frac < 0.0) cfg.max_bad_frac = 0.0;
    if (cfg.max_bad_frac > 1.0) cfg.max_bad_frac = 1.0;
  }

  return cfg;
}

// [AHCERT THETA FIELD] Theta_+ = D_i s^i + K_ij s^i s^j - trK on the GRID (s^i = unit outward coordinate-radial normal
// about the origin), with D_i s^i = d_i(sqrt(gamma) s^i)/sqrt(gamma) from second-order central differences of cell values
// (one-sided next to excised cells). The scan interpolates this field to sphere points. The previous pointwise version
// differentiated a trilinear INTERPOLANT of the state (first order, noisy across cell faces) and placed the horizon of exact
// Kerr-Schild Schwarzschild 5-10% inside 2M (GKYL_AH_SELFTEST). Component 0 = Theta_+, component 1 = valid flag.
static struct gkyl_array *
ahcert_theta_field(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, const struct gkyl_array *einstein)
{
  double exc = ctx->excision_threshold;
  struct gkyl_array *sg = gkyl_array_new(GKYL_DOUBLE, 4, app->local_ext.volume);
  struct gkyl_array *th = gkyl_array_new(GKYL_DOUBLE, 2, app->local_ext.volume);
  gkyl_array_clear(sg, 0.0);
  gkyl_array_clear(th, 0.0);

  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    long loc = gkyl_range_idx(&app->local_ext, it.idx);
    const double *es = gkyl_array_cfetch(einstein, loc);
    double *o = gkyl_array_fetch(sg, loc);
    double xc[3]; gkyl_rect_grid_cell_center(&app->grid, it.idx, xc);
    double r = sqrt(xc[0]*xc[0] + xc[1]*xc[1] + xc[2]*xc[2]);
    if (diag_spacetime_in_excision_region(exc, es) || r < 1.0e-9) { o[3] = 0.0; continue; }
    double g[3][3], gi[3][3];
    diag_read_physical_metric(es, g, gi);
    double det = diag_sym3_inverse(g, gi);
    if (!(det > 0.0)) { o[3] = 0.0; continue; }
    double v[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) v[i] += gi[i][j] * xc[j];
    double nn = xc[0]*v[0] + xc[1]*v[1] + xc[2]*v[2];
    if (!(nn > 0.0)) { o[3] = 0.0; continue; }
    double sdet = sqrt(det), invn = 1.0/sqrt(nn);
    for (int i = 0; i < 3; i++) o[i] = sdet * v[i] * invn; // sqrt(gamma) s^i
    o[3] = 1.0;
  }

  gkyl_range_iter_init(&it, &app->local);
  while (gkyl_range_iter_next(&it)) {
    long loc = gkyl_range_idx(&app->local_ext, it.idx);
    const double *o0 = gkyl_array_cfetch(sg, loc);
    if (o0[3] < 0.5) continue;
    const double *es = gkyl_array_cfetch(einstein, loc);
    double g[3][3], gi[3][3], K[3][3], Dk[3][3][3], A[3], V[3];
    diag_read_einstein_state(es, g, gi, K, Dk, A, V);
    double det = diag_sym3_inverse(g, gi);
    if (!(det > 0.0)) continue;
    double sdet = sqrt(det), divs = 0.0;
    bool ok = true;
    for (int d = 0; d < 3; d++) {
      int ip[3] = { it.idx[0], it.idx[1], it.idx[2] };
      int im[3] = { it.idx[0], it.idx[1], it.idx[2] };
      ip[d] += 1; im[d] -= 1;
      const double *op = gkyl_array_cfetch(sg, gkyl_range_idx(&app->local_ext, ip));
      const double *om = gkyl_array_cfetch(sg, gkyl_range_idx(&app->local_ext, im));
      bool pv = op[3] >= 0.5, mv = om[3] >= 0.5;
      if (pv && mv)      divs += (op[d] - om[d]) / (2.0 * app->grid.dx[d]);
      else if (pv)       divs += (op[d] - o0[d]) / app->grid.dx[d];
      else if (mv)       divs += (o0[d] - om[d]) / app->grid.dx[d];
      else { ok = false; break; }
    }
    if (!ok) continue;
    double s[3] = { o0[0]/sdet, o0[1]/sdet, o0[2]/sdet };
    double trK = 0.0, Kss = 0.0;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { trK += gi[i][j]*K[i][j]; Kss += K[i][j]*s[i]*s[j]; }
    double *t = gkyl_array_fetch(th, loc);
    t[0] = divs/sdet + Kss - trK;
    t[1] = 1.0;
  }
  gkyl_array_release(sg);
  return th;
}

// Trilinear interpolation of the Theta_+ grid field; false if the point is outside the interior cell centers or any of
// the 8 surrounding cells is invalid (excised / no stencil).
static bool
ahcert_interp_theta(gkyl_moment_app *app, const struct gkyl_array *th, const double x[3], double *val)
{
  int idx0[3];
  double w1[3];
  for (int d = 0; d < 3; ++d) {
    double pos = (x[d] - app->grid.lower[d]) / app->grid.dx[d] + 0.5;
    if (!isfinite(pos)) return false;
    int i0 = (int) floor(pos);
    if (i0 < app->local.lower[d] || i0 + 1 > app->local.upper[d]) return false;
    idx0[d] = i0;
    w1[d] = pos - i0;
  }
  double sum = 0.0;
  for (int bx = 0; bx <= 1; ++bx)
    for (int by = 0; by <= 1; ++by)
      for (int bz = 0; bz <= 1; ++bz) {
        int idx[3] = { idx0[0] + bx, idx0[1] + by, idx0[2] + bz };
        const double *t = gkyl_array_cfetch(th, gkyl_range_idx(&app->local_ext, idx));
        if (t[1] < 0.5) return false;
        double wt = (bx ? w1[0] : 1.0 - w1[0]) * (by ? w1[1] : 1.0 - w1[1]) * (bz ? w1[2] : 1.0 - w1[2]);
        sum += wt * t[0];
      }
  *val = sum;
  return true;
}

static struct ahcert_scan_result
ahcert_scan(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, struct ahcert_scan_cfg cfg)
{
  struct ahcert_scan_result res;
  memset(&res, 0, sizeof(res));

  int comm_sz = 1;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz > 1) {
    return res;
  }
  res.ran = 1;

  const struct gkyl_array *einstein = gkyl_moment_app_get_write_array_species(app, 1);
  double Rmin = cfg.Rmin;
  double Rmax = cfg.Rmax;
  int NR = cfg.NR, NT = cfg.NT, NP = cfg.NP;
  res.NR = NR; res.NT = NT; res.NP = NP;
  res.Rmin = Rmin; res.Rmax = Rmax;
  res.theta_thresh = cfg.theta_thresh;
  res.max_bad_frac = cfg.max_bad_frac;

  struct gkyl_array *th = ahcert_theta_field(app, ctx, einstein);
  double *Rs = calloc(NR, sizeof(double)), *means = calloc(NR, sizeof(double)), *areals = calloc(NR, sizeof(double));
  int *oks = calloc(NR, sizeof(int));

  double best_R = 0.0, best_max = DBL_MAX, best_min = 0.0, best_R_areal = 0.0;
  double trap_R = 0.0, trap_max = 0.0, trap_min = 0.0, trap_R_areal = 0.0;
  int best_valid = 0, trap_valid = 0, bad = 0, total = 0;

  double dtheta = M_PI / (double) NT;
  double dphi = 2.0 * M_PI / (double) NP;

  for (int ir = 0; ir < NR; ++ir) {
    double R0 = (NR == 1) ? Rmin : Rmin + (Rmax - Rmin) * ir / (double)(NR - 1);
    double max_exp = -DBL_MAX, min_exp = DBL_MAX, sum_exp = 0.0;
    double area_sum = 0.0;
    int valid = 0, area_valid = 0;
    int surf_bad = 0;
    int nsurf = NP * NT;

    for (int itheta = 0; itheta < NT; ++itheta) {
      double theta = M_PI * (itheta + 0.5) / (double) NT;
      double s = sin(theta), c = cos(theta);
      for (int iphi = 0; iphi < NP; ++iphi) {
        double phi = 2.0 * M_PI * iphi / (double) NP;
        double x[3] = { R0*s*cos(phi), R0*s*sin(phi), R0*c };
        double expansion = 0.0;
        total += 1;
        if (!ahcert_interp_theta(app, th, x, &expansion)) {
          bad += 1;
          surf_bad += 1;
          continue;
        }
        if (expansion > max_exp) max_exp = expansion;
        if (expansion < min_exp) min_exp = expansion;
        sum_exp += expansion;
        valid += 1;

        double dA = 0.0;
        if (ahcert_spherical_area_density(app, einstein, ctx, x, R0, theta, phi, &dA)) {
          area_sum += dA * dtheta * dphi;
          area_valid += 1;
        }
      }
    }

    if (valid == 0 || (double) surf_bad > cfg.max_bad_frac * (double) nsurf) continue;

    double R_areal = 0.0;
    if (area_valid > 0) {
      double area_est = area_sum * (double) nsurf / (double) area_valid;
      R_areal = sqrt(area_est / (4.0 * M_PI));
    }
    Rs[ir] = R0; means[ir] = sum_exp / (double) valid; areals[ir] = R_areal; oks[ir] = (R_areal > 0.0);

    if (fabs(max_exp) < fabs(best_max)) {
      best_R = R0;
      best_max = max_exp;
      best_min = min_exp;
      best_valid = valid;
      best_R_areal = R_areal;
    }
    if (max_exp <= cfg.theta_thresh && R0 > trap_R) {
      trap_R = R0;
      trap_max = max_exp;
      trap_min = min_exp;
      trap_valid = valid;
      trap_R_areal = R_areal;
    }
  }

  // Outermost + -> <=0 crossing of the sphere-averaged Theta_+ (scanning inward), linear in R.
  for (int ir = NR - 1; ir >= 1; --ir) {
    if (!(oks[ir] && oks[ir - 1])) continue;
    if (means[ir] > 0.0 && means[ir - 1] <= 0.0) {
      double w = means[ir - 1] / (means[ir - 1] - means[ir]);
      res.root_R = Rs[ir - 1] + w * (Rs[ir] - Rs[ir - 1]);
      res.root_R_areal = areals[ir - 1] + w * (areals[ir] - areals[ir - 1]);
      res.root_valid = 1;
      break;
    }
  }
  free(Rs); free(means); free(areals); free(oks);
  gkyl_array_release(th);

  res.best_R = best_R; res.best_max = best_max; res.best_min = best_min;
  res.best_R_areal = best_R_areal; res.best_valid = best_valid;
  res.trap_R = trap_R; res.trap_max = trap_max; res.trap_min = trap_min;
  res.trap_R_areal = trap_R_areal; res.trap_valid = trap_valid;
  res.bad = bad; res.total = total;
  return res;
}

// Spherical trapped-surface certificate for the nearly spherical TOV collapse.
// It reports trapped=1 only when every sampled point on a closed coordinate
// sphere has Theta_+ <= threshold. This is stricter than the shell-averaged AH
// monitor and is meant as a pass/fail check, not a cheap always-on diagnostic.
static void
ah_spherical_cert_diag(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  if (getenv("GKYL_AH_CERT") == NULL) return;

  int comm_sz = 1;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz > 1) {
    return;
  }

  static long calls = 0;
  calls += 1;
  double tmin = -DBL_MAX;
  const char *tmin_env = getenv("GKYL_AH_CERT_TMIN");
  if (tmin_env) tmin = atof(tmin_env);
  if (t_curr < tmin) return;

  int every = 1;
  const char *ev = getenv("GKYL_AH_CERT_EVERY");
  if (ev) every = atoi(ev);
  else {
    ev = getenv("GKYL_AH_EVERY");
    if (ev) every = atoi(ev);
  }
  if (every < 1) every = 1;
  if ((calls - 1) % every != 0) return;

  // [AHCERT bad-sample handling] By default a surface is vetoed if ANY angular sample is bad.
  // GKYL_AH_CERT_MAX_BAD_FRAC=f (>0) tolerates a fraction f of excised-interior samples.
  struct ahcert_scan_cfg cert_cfg = ahcert_scan_cfg_from_env(app, "GKYL_AH_CERT",
    0.82*app->grid.upper[0], 32, 21, 32);
  struct ahcert_scan_result cert = ahcert_scan(app, ctx, cert_cfg);

  printf("[AHCERT] t=%.4f trapped=%d r=%.4f thetaMax=%+.6e thetaMin=%+.6e R_areal=%.4f M_AH=%.4f "
    "bestR=%.4f bestR_areal=%.4f bestMax=%+.6e bestMin=%+.6e valid=%d trap_valid=%d bad=%d total=%d scan=(%d,%d,%d) Rrange=[%.4f,%.4f]\n",
    t_curr, cert.trap_R > 0.0 ? 1 : 0, cert.trap_R, cert.trap_max, cert.trap_min,
    cert.trap_R_areal, 0.5 * cert.trap_R_areal,
    cert.best_R, cert.best_R_areal, cert.best_max, cert.best_min,
    cert.best_valid, cert.trap_valid, cert.bad, cert.total, cert.NR, cert.NT, cert.NP,
    cert.Rmin, cert.Rmax);
  printf("[AHCERT-ROOT] t=%.4f root_valid=%d R_AH=%.4f R_AH_areal=%.4f M_AH=%.4f (zero of sphere-averaged Theta_+)\n",
    t_curr, cert.root_valid, cert.root_R, cert.root_R_areal, 0.5 * cert.root_R_areal);
  fflush(stdout);
  { extern void route1_fx_ah(double t, int valid, double R); route1_fx_ah(t_curr, cert.root_valid, cert.root_R); } // [FX]
}





// [DYN] Mid-run dynamical-collapse trigger: at t >= GKYL_DYN_TRIGGER_TIME (default 2), apply a
// Energy-preserving pressure depletion to the LIVE fluid state in r < R_core, once.
// This turns a held/marginal star into an under-pressured one that then collapses -- the
// "keep it in equilibrium, then trigger the collapse" story in a SINGLE run. Gated by GKYL_DYN_COLLAPSE.
// Depletion fraction = GKYL_KDROP (default zero, a release-only control). WB is allowed during the hold phase and is
// switched off at this same event before the released star evolves. Recovers primitives (C2P) at the live
// metric, rescales rho to preserve either the fluid-frame ETK-91 energy density or, with
// GKYL_KDROP_ETK_EULERIAN, the normal-observer energy density used by the KS initial-data path, and
// rebuilds the conserved variables at the live metric. The Eulerian variant
// preserves BOTH E and S_i by changing velocity, not merely rescaling density.
static bool
dyn_collapse_enabled(void)
{
  return collapse_flag_enabled("GKYL_DYN_COLLAPSE");
}

static double
dyn_collapse_trigger_time(void)
{
  double t_trig = 1.0;
  const char *tt = getenv("GKYL_DYN_TRIGGER_TIME");
  if (tt != NULL) {
    t_trig = atof(tt);
  }
  return t_trig;
}

static bool dyn_fired = false; // [FX] file scope so a restart can mark the trigger fired

static bool
dyn_collapse_trigger(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double t_curr)
{
  if (dyn_fired) return true;
  if (!dyn_collapse_enabled()) return true;
  double t_trig = dyn_collapse_trigger_time();
  if (t_curr < t_trig) return true;

  double G = ctx->gas_gamma;
  double Korig = ctx->K_poly;
  double delta = ctx->kdrop_delta;
  // Commit all cells together, or leave the original state and WB untouched.
  struct gkyl_array *fluid = gkyl_array_clone(app->species[0].fcurr);
  bool eulerian = getenv("GKYL_KDROP_ETK_EULERIAN") != NULL;
  double audit[5]={0}, global_audit[5], maxerr=0, global_maxerr;
  int invalid=0, global_invalid;
  int ncells = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local_ext, iter.idx);
    double xc[3]; gkyl_rect_grid_cell_center(&app->grid, iter.idx, xc);
    double r = sqrt(xc[0]*xc[0] + xc[1]*xc[1] + xc[2]*xc[2]);
    if (delta == 0) continue; // release-only control: exactly no state remap
    if (r >= ctx->R_core) continue;

    double *q = gkyl_array_fetch(fluid, loc);
    if (!(isfinite(q[0]) && q[0] > 0.0)) { invalid++; continue; }

    double v[71];
    gkyl_gr_euler_prim_vars(G, q, v);
    double rho = v[0], p = v[4];
    if (!(rho > ctx->rho_atm) || !(p > 0.0)) continue;

    double local_delta = delta;
    if (ctx->kdrop_smooth) {
      local_delta = delta * 0.5 * (1.0 + cos(M_PI * r / ctx->R_core));
    }
    // Deplete the live effective K. C2P(cell average) need not obey the
    // pointwise initial polytrope, especially across the surface. Resetting
    // it to Korig would impose a finite pressure jump even as delta -> 0.
    double Kp = (1.0 - local_delta) * p / pow(rho, G);

    // P2C at the LIVE metric (fluid slots 9-17), same convention as the IC.
    double g[3][3] = {
      { q[9],  q[10], q[11] },
      { q[12], q[13], q[14] },
      { q[15], q[16], q[17] } };
    double det = g[0][0]*(g[1][1]*g[2][2] - g[1][2]*g[2][1])
               - g[0][1]*(g[1][0]*g[2][2] - g[1][2]*g[2][0])
               + g[0][2]*(g[1][0]*g[2][1] - g[1][1]*g[2][0]);
    if (!(isfinite(det) && det>0)) { invalid++; continue; }
    double sdet = sqrt(det);
    double vel[3] = { v[1], v[2], v[3] };
    double vsq = 0.0, cov[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        vsq += g[i][j] * vel[i] * vel[j];
        cov[i] += g[i][j] * vel[j];
      }
    if (!(isfinite(vsq) && vsq>=0 && vsq<1)) { invalid++; continue; }
    double W = 1.0 / sqrt(1.0 - vsq);

    double oldE=q[0]+q[4], oldD=q[0], oldS[3]={q[1],q[2],q[3]};
    double pp;
    if (eulerian) {
      double inverse[3][3], *rows[3]={inverse[0],inverse[1],inverse[2]}, **ip=rows;
      gkyl_gr_euler_inv_spatial_metric(q,&ip);
      if (!rt_kdrop_es_conserved(G,Kp,sdet,inverse,ctx->rho_atm,ctx->p_atm,q,&pp)) {
        invalid++; continue;
      }
    }
    else {
      double rp;
      if (!rt_kdrop_density(G,Kp,1,rho+p/(G-1),&rp)) { invalid++; continue; }
      pp=Kp*pow(rp,G);
      if (!(rp>=ctx->rho_atm && pp>=ctx->p_atm && isfinite(pp))) { invalid++; continue; }
      double h=1+(pp/rp)*G/(G-1);
      q[0]=sdet*rp*W;
      for (int i=0;i<3;++i) q[1+i]=sdet*rp*h*W*W*cov[i];
      q[4]=sdet*(rp*h*W*W-pp-rp*W);
    }
    maxerr=fmax(maxerr,fabs((q[0]+q[4])-oldE)/fabs(oldE));
    audit[0]+=q[0]-oldD;
    audit[1]+=(q[0]+q[4])-oldE;
    double ds=0; for (int i=0;i<3;++i) ds+=(q[1+i]-oldS[i])*(q[1+i]-oldS[i]);
    audit[2]+=sqrt(ds);
    audit[3]+=pp>p ? 1 : 0;
    ncells++;
  }
  audit[4]=ncells;
  gkyl_comm_allreduce(app->comm,GKYL_INT,GKYL_MAX,1,&invalid,&global_invalid);
  if (global_invalid) {
    gkyl_array_release(fluid);
    fprintf(stderr,"[DYN] rejected invalid remap; state and WB unchanged\n");
    return false;
  }
  gkyl_comm_allreduce(app->comm,GKYL_DOUBLE,GKYL_SUM,5,audit,global_audit);
  gkyl_comm_allreduce(app->comm,GKYL_DOUBLE,GKYL_MAX,1,&maxerr,&global_maxerr);
  if (delta>0 && global_audit[4]==0) {
    gkyl_array_release(fluid);
    fprintf(stderr,"[DYN] rejected empty depletion region; state and WB unchanged\n");
    return false;
  }
  gkyl_array_copy(app->species[0].fcurr,fluid);
  gkyl_array_release(fluid);
  dyn_fired=true;

  // Switch the well-balancer OFF from here on, so the depleted (out-of-equilibrium) star is free to
  // collapse instead of being pinned back to the frozen t=0 equilibrium. This is the "WB off + deplete"
  // instant: WB holds the marginal star for t < t_trigger, then releases it here.
  bool wb_was_on = false;
  if (ctx->fluid_eqn != NULL && !gkyl_gr_euler_wb_disabled(ctx->fluid_eqn)) {
    gkyl_gr_euler_set_wb_disabled(ctx->fluid_eqn, true);
    wb_was_on = true;
  }


  int rank = 0; gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    double dv=app->grid.dx[0]*app->grid.dx[1]*app->grid.dx[2];
    printf("[DYN-BUDGET] t=%.17g cells=%.0f deltaD=%.17g deltaE=%.17g integralAbsDeltaS=%.17g maxRelEnergyRemap=%.17g pressureRaised=%.0f\n",
      t_curr,global_audit[4],dv*global_audit[0],dv*global_audit[1],dv*global_audit[2],global_maxerr,global_audit[3]);
    printf("[DYN] t=%.4f FIRED %s depletion delta=%.3f (live K scaled; initial K=%.4g) in %d core cells (r<%.3f, smooth=%d); WB %s; slicing %s\n",
      t_curr, getenv("GKYL_KDROP_ETK_EULERIAN") != NULL ? "ETK-Eulerian" : "ETK-91",
      delta, Korig, (int)global_audit[4], ctx->R_core, (int) ctx->kdrop_smooth, wb_was_on ? "switched OFF" : "already off",
      "unchanged");
    fflush(stdout);
  }
  moment_species_apply_bc(app,t_curr,&app->species[0],app->species[0].fcurr);
  return true;
}

// [AH SELF-TEST] Known-answer check of this driver's trapped-sphere certificate
// on exact Kerr-Schild Schwarzschild data built with
// this driver's own conformal initial-data routine (dev data convention). gkyl_gr_blackhole returns K in its own tested
// convention 2 alpha K = -(D_i beta_j + D_j beta_i) (ctest_gr_spacetime.c), so K is negated into the solver's MTW
// convention. Expected: outermost trapped coordinate sphere R -> 2M, areal R_AH = 2M, M_AH = M. Pass 2 flips K back
// (wrong convention for this solver): a correct finder must report NO trapped sphere.
// GKYL_AH_SELFTEST=<M>, black hole at the origin, excision only for r < M/2; exits after printing.
static double ah_selftest_rexc = 0.0;

static void
ah_selftest_excision(const struct gkyl_gr_spacetime *spacetime, const double t, const double x, const double y,
  const double z, bool *in_excision_region)
{
  *in_excision_region = sqrt(x*x + y*y + z*z) < ah_selftest_rexc;
}

static void
ah_selftest(gkyl_moment_app *app, struct gr_tov_coupled_ctx *ctx, double mass)
{
  struct gkyl_gr_spacetime *bh = gkyl_gr_blackhole_new(false, mass, 0.0, 0.0, 0.0, 0.0);
  ah_selftest_rexc = 0.5*mass;
  bh->excision_region_func = ah_selftest_excision;
  struct gkyl_gr_spacetime *saved = ctx->spacetime;
  ctx->spacetime = bh;
  struct gkyl_array *es = (struct gkyl_array *) gkyl_moment_app_get_write_array_species(app, 1);
  struct gkyl_range_iter it;
  gkyl_range_iter_init(&it, &app->local_ext);
  while (gkyl_range_iter_next(&it)) {
    double xc[3]; gkyl_rect_grid_cell_center(&app->grid, it.idx, xc);
    double *q = gkyl_array_fetch(es, gkyl_range_idx(&app->local_ext, it.idx));
    evalVacuumEinsteinConformalInit3d(0.0, xc, q, ctx);
    for (int m = 10; m < 19; ++m) q[m] = -q[m];
  }
  ctx->spacetime = saved;

  for (int pass = 0; pass < 2; ++pass) {
    struct ahcert_scan_cfg cfg = ahcert_scan_cfg_from_env(app, "GKYL_AH_SELFTEST", 4.0*mass, 66, 25, 40);
    cfg.Rmin = 0.75*mass;
    cfg.Rmax = 4.0*mass;
    struct ahcert_scan_result r = ahcert_scan(app, ctx, cfg);
    printf("[AH-SELFTEST] pass=%d (%s K) M=%g dx=%.4f dR=%.4f: HORIZON root_valid=%d R_AH=%.5f (err %+.2e) areal=%.5f"
      " M_AH=%.5f (err %+.2e) | certificate trapped=%d trap_R=%.4f | bad=%d/%d\n",
      pass, pass == 0 ? "solver-convention" : "FLIPPED", mass, app->grid.dx[0], (cfg.Rmax - cfg.Rmin)/(cfg.NR - 1),
      r.root_valid, r.root_R, r.root_R - 2.0*mass, r.root_R_areal, 0.5*r.root_R_areal, 0.5*r.root_R_areal - mass,
      r.trap_R > 0.0 ? 1 : 0, r.trap_R, r.bad, r.total);
    fflush(stdout);
    gkyl_range_iter_init(&it, &app->local_ext);
    while (gkyl_range_iter_next(&it)) {
      double *q = gkyl_array_fetch(es, gkyl_range_idx(&app->local_ext, it.idx));
      for (int m = 10; m < 19; ++m) q[m] = -q[m];
    }
  }
  gkyl_gr_spacetime_release(bh);
}

int
main(int argc, char **argv)
{
  collapse_run_defaults();
  // Keep complete diagnostic lines together when stdout and stderr share a log.
  setvbuf(stdout, NULL, _IOLBF, 0);
#ifdef GKYL_HAVE_MPI
  if (argc > 1) {
    struct gkyl_app_args app_args_pre = parse_app_args(argc, argv);
    if (app_args_pre.use_mpi) MPI_Init(&argc, &argv);
  }
#endif

  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gr_tov_coupled_ctx ctx = create_ctx();

  // A static KS TOV is not stationary under 1+log. Do not subtract its genuine
  // gauge evolution as an equilibrium residual. Excision needs separate tests.
  if (GKYL_COLLAPSE_INITIAL_KS || ctx.excision_threshold != 0.0) {
    fprintf(stderr, "This equilibrium-corrected driver requires areal initial data and zero excision threshold.\n");
    gkyl_gr_spacetime_release(ctx.spacetime);
    gkyl_tov_solution_release(ctx.tov);
    return 1;
  }

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ny);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], ctx.Nz);

  if (dyn_collapse_enabled() && !(isfinite(ctx.kdrop_delta) && ctx.kdrop_delta >= 0.0 && ctx.kdrop_delta < 1.0
      && isfinite(dyn_collapse_trigger_time()) && dyn_collapse_trigger_time() >= 0.0)) {
    fprintf(stderr, "Scheduled release requires 0 <= KDROP < 1 and a finite nonnegative trigger time.\n");
    return 1;
  }

  // This problem starts held; the scheduled trigger disables only fluid WB.
  bool wb_on = collapse_flag_enabled("GKYL_WB");
  bool wb_off = !wb_on;
  printf("WB: disable_well_balanced=%d\n", (int) wb_off);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_inew(&(struct gkyl_wv_gr_euler_inp) {
    .gas_gamma = ctx.gas_gamma,
    .spacetime_gauge = ctx.spacetime_gauge,
    .reinit_freq = ctx.reinit_freq,
    .spacetime = ctx.spacetime,
    .rp_type = WV_GR_EULER_RP_LAX,
    .use_gpu = app_args.use_gpu,
    .tov_eq = ctx.tov,
    .disable_well_balanced = wb_off,
    .p_atm = ctx.p_atm,
    .rho_atm = ctx.rho_atm,
  });
  ctx.fluid_eqn = gr_euler; // handle for the mid-run WB toggle (dyn_collapse_trigger)

  // Vacuum Einstein (Bona-Masso), evolved with wave propagation and coupled to the fluid.
  // Spacetime-sector limiter is a driver choice. Keep runtime limiter sweeps in a diagnostic patch.
  enum gkyl_wave_limiter st_limiter = GKYL_COLLAPSE_ST_LIMITER;
  printf("Spacetime limiter: %d (1=none,2=minmod,5=MC,8=EP)\n", (int) st_limiter);

  // This is the conformal collapse driver. Non-conformal comparisons should use a separate driver.
  bool use_conformal = true;
  diag_is_conformal = use_conformal;

  printf("Formulation: CONFORMAL (77)\n");

  struct gkyl_wv_eqn *vacuum_einstein;
  vacuum_einstein = gkyl_wv_vacuum_einstein_conformal_inew(
    &(struct gkyl_wv_vacuum_einstein_conformal_inp) {
      .excision_threshold = ctx.excision_threshold,
      .spacetime_slicing = ctx.spacetime_slicing,
      .spacetime_evolution = ctx.spacetime_evolution,
      .rp_type = WV_VACUUM_EINSTEIN_CONFORMAL_RP_HLL,
      .use_gpu = app_args.use_gpu,
    });
  printf("Spacetime solver: HLL; t_end=%g frames=%d; trigger=%d at t=%g\n",
    ctx.t_end, ctx.num_frames, dyn_collapse_enabled(), dyn_collapse_trigger_time());
  printf("Shift: reference-Gamma mu=%s eta=%s (experimental=%s); fluid WB holds it until release\n",
    getenv("GKYL_GAMMA_DRIVER"), getenv("GKYL_GAMMA_DRIVER_ETA"), getenv("GKYL_EXPERIMENTAL"));

  enum gkyl_species_bc_type einstein_bc = GKYL_SPECIES_RADIAL_FALLOFF;
  const char *sommerfeld_env = getenv("GKYL_USE_FALLOFF_BC");
  if (sommerfeld_env != NULL) {
    if (strcmp(sommerfeld_env, "0") == 0) {
      einstein_bc = GKYL_SPECIES_COPY;
    }
  }
  const char *einstein_bc_name = "Sommerfeld";
  if (einstein_bc == GKYL_SPECIES_COPY) {
    einstein_bc_name = "Copy";
  }
  printf("Einstein BC: %s\n", einstein_bc_name);

  bool force_low_order_spacetime = false;
  if (force_low_order_spacetime) {
    printf("Spacetime wave-prop: force low-order flux\n");
  }

  struct gkyl_moment_species fluid = {
    .name = "gr_euler",
    .equation = gr_euler,
    .init = evalGREulerInit,
    .ctx = &ctx,

    .scheme_type = GKYL_MOMENT_MP, // fluid sector: MP + frozen-discrete WB

    .has_gr_euler = true,
    .gr_euler_gas_gamma = ctx.gas_gamma,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcy = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
    .bcz = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },
  };

  struct gkyl_moment_species einstein = {
    .name = "vacuum_einstein",
    .equation = vacuum_einstein,
    .init = use_conformal ? evalVacuumEinsteinConformalInit3d : evalVacuumEinsteinInit,
    .ctx = &ctx,

    .scheme_type = GKYL_MOMENT_WAVE_PROP, // spacetime sector: wave-propagation
    .limiter = st_limiter,
    .force_low_order_flux = force_low_order_spacetime,

    .has_vacuum_einstein = !use_conformal,
    .vacuum_einstein_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_spacetime_evolution = ctx.spacetime_evolution,

    .has_vacuum_einstein_conformal = use_conformal,
    .vacuum_einstein_conformal_excision_threshold = ctx.excision_threshold,
    .vacuum_einstein_conformal_spacetime_slicing = ctx.spacetime_slicing,
    .vacuum_einstein_conformal_spacetime_evolution = ctx.spacetime_evolution,

    .bcx = { einstein_bc, einstein_bc },
    .bcy = { einstein_bc, einstein_bc },
    .bcz = { einstein_bc, einstein_bc },
  };

  int cells[] = { NX, NY, NZ };
  int dim = sizeof(cells) / sizeof(cells[0]);
  int cuts[dim];
  for (int d = 0; d < dim; d++) cuts[d] = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) for (int d = 0; d < dim; d++) cuts[d] = app_args.cuts[d];
#endif

  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) { .mpi_comm = MPI_COMM_WORLD } );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app_args.use_gpu } );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { .use_gpu = app_args.use_gpu } );
#endif

  struct gkyl_moment app_inp = {
    .name = "gr_core_collapse_conformal",

    .ndim = 3,
    .lower = { -0.5 * ctx.Lx, -0.5 * ctx.Ly, -0.5 * ctx.Lz },
    .upper = {  0.5 * ctx.Lx,  0.5 * ctx.Ly,  0.5 * ctx.Lz },
    .cells = { NX, NY, NZ },

    .scheme_type = GKYL_MOMENT_MIXED,
    .mp_recon = app_args.mp_recon,

    .cfl_frac = ctx.cfl_frac,

    .num_species = 2,
    .species = { fluid, einstein },

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { cuts[0], cuts[1], cuts[2] },
      .comm = comm,
    },
  };

  gkyl_moment_app *app = gkyl_moment_app_new(&app_inp);

  double t_curr = 0.0, t_end = ctx.t_end;
  int frame_curr = 0;
  // Restart files do not persist the one-shot trigger/WB state. Refuse an
  // ambiguous restart rather than silently applying the pressure drop twice.
  if (app_args.is_restart) {
    fprintf(stderr,"Collapse restart needs persisted equilibrium, gauge and trigger state; start a fresh campaign.\n");
    gkyl_moment_app_release(app);
    gkyl_wv_eqn_release(gr_euler);
    gkyl_wv_eqn_release(vacuum_einstein);
    gkyl_gr_spacetime_release(ctx.spacetime);
    gkyl_tov_solution_release(ctx.tov);
    gkyl_comm_release(comm);
    return 1;
  }
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status =
      gkyl_moment_app_read_from_frame(app, app_args.restart_frame);
    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_moment_app_cout(app, stderr, "*** Failed to read restart frame %d, status=%d\n",
        app_args.restart_frame, status.io_status);
      gkyl_moment_app_release(app);
      gkyl_wv_eqn_release(gr_euler);
      gkyl_wv_eqn_release(vacuum_einstein);
      gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
      gkyl_tov_solution_release(ctx.tov);
      gkyl_comm_release(comm);
      return 1;
    }
    t_curr = status.stime;
    frame_curr = status.frame;
    gkyl_moment_app_cout(app, stdout, "Restarted from frame %d at t = %g\n", frame_curr, t_curr);
  }
  else {
    gkyl_moment_app_apply_ic(app, t_curr);
    if (!moment_spacetime_set_equilibrium(app, &app->species[1],
        app->species[0].fcurr, ctx.gas_gamma)) {
      fprintf(stderr, "Stationary spacetime reference rejected.\n");
      gkyl_moment_app_release(app);
      gkyl_wv_eqn_release(gr_euler);
      gkyl_wv_eqn_release(vacuum_einstein);
      gkyl_gr_spacetime_release(ctx.spacetime);
      gkyl_tov_solution_release(ctx.tov);
      gkyl_comm_release(comm);
      return 1;
    }
    printf("[SPACETIME-EQUILIBRIUM] shared-face perturbation flux; source reference persists after trigger\n");
    if (!sync_initial_fluid_metric(app, 0, 1)) {
      fprintf(stderr,"Initial fluid state is inadmissible after metric synchronization.\n");
      gkyl_moment_app_release(app);
      gkyl_wv_eqn_release(gr_euler);
      gkyl_wv_eqn_release(vacuum_einstein);
      gkyl_gr_spacetime_release(ctx.spacetime);
      gkyl_tov_solution_release(ctx.tov);
      gkyl_comm_release(comm);
      return 1;
    }
  }

  // [FX-RESTART] GKYL_FX_RESTART_FRAME=N (this private build): after the t = 0 initial data and references above,
  // read frame N of both species from the run directory (a copy of the original run's frame files), set the clock,
  // mark the trigger fired (WB off, no remap). The Gamma-driver reference comes from GKYL_FX_DRIVER_REF_LOAD, the
  // band speed from GKYL_ROUTE1_BAND_S; NewRad's evolved ghosts restart from the boundary fill.
  if (getenv("GKYL_FX_RESTART_FRAME")) {
    int rframe = atoi(getenv("GKYL_FX_RESTART_FRAME"));
    struct gkyl_app_restart_status s0 = gkyl_moment_app_from_frame_species(app, 0, rframe);
    // The library refuses a Gamma-driver spacetime restart because the activation reference is not in the frame;
    // this build reloads it (GKYL_FX_DRIVER_REF_LOAD, required here), so the spacetime array is read directly.
    if (!getenv("GKYL_FX_DRIVER_REF_LOAD")) { fprintf(stderr, "[FX-RESTART] needs GKYL_FX_DRIVER_REF_LOAD\n"); return 1; }
    char sfn[256];
    snprintf(sfn, sizeof sfn, "%s-%s_%d.gkyl", app->name, app->species[1].name, rframe);
    struct gkyl_app_restart_status s1 = { .io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local,
      app->species[1].fcurr, sfn), .frame = rframe, .stime = s0.stime };
    if (s0.io_status != GKYL_ARRAY_RIO_SUCCESS || s1.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      fprintf(stderr, "[FX-RESTART] failed to read frame %d (status %d, %d)\n", rframe, s0.io_status, s1.io_status);
      return 1;
    }
    t_curr = s0.stime; frame_curr = rframe; app->tcurr = t_curr;
    for (int s = 0; s < 2; ++s) moment_species_apply_bc(app, t_curr, &app->species[s], app->species[s].fcurr);
    dyn_fired = true;
    if (ctx.fluid_eqn != NULL && !gkyl_gr_euler_wb_disabled(ctx.fluid_eqn)) gkyl_gr_euler_set_wb_disabled(ctx.fluid_eqn, true);
    printf("[FX-RESTART] frame %d t=%.17g (fluid %.17g, spacetime %.17g): trigger marked fired, no remap, WB off\n",
      rframe, t_curr, s0.stime, s1.stime);
    fflush(stdout);
  }
  if (!dyn_collapse_trigger(app,&ctx,t_curr)) {
    gkyl_moment_app_release(app);
    gkyl_wv_eqn_release(gr_euler);
    gkyl_wv_eqn_release(vacuum_einstein);
    gkyl_gr_spacetime_release(ctx.spacetime);
    gkyl_tov_solution_release(ctx.tov);
    gkyl_comm_release(comm);
    return 1;
  }
  struct gkyl_tm_trigger io_trig = { .dt = t_end / ctx.num_frames, .tcurr = t_curr, .curr = frame_curr };
  diag_track(app, &ctx, t_curr);
  constraint_monitor(app, &ctx, t_curr);
  ah_spherical_cert_diag(app, &ctx, t_curr);
  if (getenv("GKYL_AH_SELFTEST") != NULL) {
    ah_selftest(app, &ctx, atof(getenv("GKYL_AH_SELFTEST")));
  }
  if (getenv("GKYL_DIAG_ONLY") != NULL || getenv("GKYL_AH_SELFTEST") != NULL) {
    gkyl_moment_app_release(app);
    gkyl_wv_eqn_release(gr_euler);
    gkyl_wv_eqn_release(vacuum_einstein);
    gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
    gkyl_tov_solution_release(ctx.tov);
    gkyl_comm_release(comm);
    return 0;
  }
  write_data(&io_trig, app, t_curr, false);

  double dt = fmin(t_end - t_curr, gkyl_moment_app_max_dt(app));
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  struct gkyl_array *snap_fluid = gkyl_array_clone(app->species[0].fcurr);
  struct gkyl_array *snap_einstein = gkyl_array_clone(app->species[1].fcurr);
  const int max_spd_retry = 12;
  double guard_dt_ceiling = DBL_MAX;
  double guard_dt_growth = 1.25;
  const char *guard_dt_growth_env = getenv("GKYL_GUARD_DT_GROWTH");
  if (guard_dt_growth_env != NULL) {
    double v = atof(guard_dt_growth_env);
    if (v > 1.0) {
      guard_dt_growth = v;
    }
  }

  bool run_failed = false;
  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    dt = fmin(dt, t_end - t_curr);
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);

    if (dyn_collapse_enabled()) {
      double dt_to_trigger = dyn_collapse_trigger_time() - t_curr;
      if (dt_to_trigger > 0.0 && dt > dt_to_trigger) {
        dt = dt_to_trigger;
      }
    }

    gkyl_array_copy(snap_fluid, app->species[0].fcurr);
    gkyl_array_copy(snap_einstein, app->species[1].fcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    int spd_retry = 0;
    while (!status.success && status.dt_suggested > 0.0 && status.dt_suggested < dt && spd_retry < max_spd_retry) {
      gkyl_array_copy(app->species[0].fcurr, snap_fluid);
      gkyl_array_copy(app->species[1].fcurr, snap_einstein);
      dt = status.dt_suggested;
      gkyl_moment_app_cout(app, stdout, " [update retry] rejected, retry dt = %g ...", dt);
      status = gkyl_moment_update(app, dt);
      spd_retry++;
    }
    bool had_guard_retry = spd_retry > 0;
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_array_copy(app->species[0].fcurr, snap_fluid);
      gkyl_array_copy(app->species[1].fcurr, snap_einstein);
      gkyl_moment_app_cout(app, stdout, "** Update method failed (retries exhausted)! Aborting simulation ....\n");
      run_failed = true;
      break;
    }

    t_curr += status.dt_actual;
    // If a hyperbolicity/admissibility guard rejected the trial step, keep the
    // next step close to the accepted retry instead of jumping straight back to
    // the large CFL estimate and repeatedly re-hitting the same rejected state.
    if (had_guard_retry) {
      guard_dt_ceiling = status.dt_actual;
    }
    else if (guard_dt_ceiling < DBL_MAX) {
      guard_dt_ceiling *= guard_dt_growth;
      if (guard_dt_ceiling >= status.dt_suggested) {
        guard_dt_ceiling = DBL_MAX;
      }
    }
    dt = fmin(status.dt_suggested, guard_dt_ceiling);

    if (!dyn_collapse_trigger(app, &ctx, t_curr)) { run_failed=true; break; }
    write_data(&io_trig, app, t_curr, false);
    diag_track(app, &ctx, t_curr);
    constraint_monitor(app, &ctx, t_curr);
    ah_spherical_cert_diag(app, &ctx, t_curr);

    if (dt_init < 0.0) dt_init = status.dt_actual;
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;
      gkyl_moment_app_cout(app, stdout, "WARNING: dt = %g below %g*dt_init; num_failures = %d\n",
        status.dt_actual, dt_failure_tol, num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: too many small steps. Aborting.\n");
        run_failed = true;
        break;
      }
    }
    else num_failures = 0;

    step += 1;
  }

  write_data(&io_trig, app, t_curr, true);
  gkyl_moment_app_stat_write(app);
  bool completed = !run_failed && t_curr >= t_end;
  gkyl_moment_app_cout(app, stdout, "[RUN-STATUS] completed=%d t=%.17g target=%.17g failed=%d\n",
    completed, t_curr, t_end, run_failed);

  gkyl_array_release(snap_fluid);
  gkyl_array_release(snap_einstein);
  gkyl_moment_app_release(app);
  gkyl_wv_eqn_release(gr_euler);
  gkyl_wv_eqn_release(vacuum_einstein);
  gkyl_gr_tov_spacetime_free(&ctx.spacetime->ref_count);
  gkyl_tov_solution_release(ctx.tov);
  gkyl_comm_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) MPI_Finalize();
#endif

  return completed ? 0 : 1;
}
