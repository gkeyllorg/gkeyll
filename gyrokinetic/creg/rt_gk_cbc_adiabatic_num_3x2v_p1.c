// Numerical Miller equivalent of gk_cbc_nl_adiabatic_c2p_3x2v_p1.c.
// Radial input bounds are outboard major radii; parallel bounds are Miller
// angles. The app converts both to computational coordinates. Callbacks receive
// (psi,alpha,theta). Profiles use discrete mc2p to recover
// rMr0=rminor-rminor0 and retain ylength=Cy*alpha. The radial mesh
// is uniform in psi, rather than in the original rMr0 coordinate.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double Rmajor_axis, Zvert_axis, Rmajor0, rminor_LCFS, rminor0, bmag_axis, bmag0, q0, Cy;
  double kappa; // Miller elongation.
  double delta; // Miller triangularity.
  double a_shift; // Shafranov-shift strength.
  // Plasma parameters
  double me, qe, mi, qi, n0, Te0, Ti0;
  double c_s, rho_s;
  // Initial condition parameters
  bool const_grad; // Whether to use profile with constant gradient.
  double Ln, LTi,
    prof_width; // Gradient lengths at rMr0=0 and tanh width of the profiles.
  double delta_n, pert_width; // Amplitude and radial width of the seed.
  int n_tor; // Toroidal mode number of the box (n_tor*L_ylength = 2 pi Cy).
  double nu_krook; // Krook buffer rate.
  double nu_krook_norm; // In units of c_s/Rmajor_axis.
  int num_cell_krook_buff;
  double Rmajor_LCFS; // Outboard major radius of the Miller LCFS [m].
  double R_omp_min,
    R_omp_max; // Radial initialization bounds [m]; app converts to psi.
  double ylength_min, ylength_max; // Requested binormal-length bounds [m].
  // Reference profile extents and numerical grid counts.
  double L_rMr0, L_ylength;
  double rMr0_min, thetaMiller_min, rMr0_max, thetaMiller_max;
  int N_psi, N_alpha, N_theta, Nvpar, Nmu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_ion, mu_max_ion, vpar_max_elc, mu_max_elc;
  double nu_star, nu_ion; // Ion collisionality nu_ii q Rmajor_axis/(eps^{3/2}
    // v_ti) and the LBO frequency; 0 = collisionless.
  // Simulation control parameters
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
  struct gkyl_gen_miller_eqdsk_inp miller;
  const struct gk_geometry *geo; // Borrowed host geometry, owned by the app.
};

/**
 * Safety-factor profile (Goerler et al. 2011, Sec. 4.1).
 *
 * @param rminor Minor radius [m].
 * @param rminor_LCFS LCFS minor radius [m].
 * @return Positive safety factor.
 */
static double
qsafety_at_rminor(double rminor, double rminor_LCFS)
{
  double normalized_radius = rminor / rminor_LCFS;
  return 0.854 + 2.373 * pow(normalized_radius, 2) - 0.466 * pow(normalized_radius, 3) +
         0.498 * pow(normalized_radius, 4);
}

/**
 * Supply the Miller generator with q as a function of outboard major radius.
 *
 * @param Rmajor_omp Outboard midplane major radius [m].
 * @param ctx CBC simulation context.
 * @return Safety factor for the Miller surface, including radial shift.
 */
static double
qsafety(double Rmajor_omp, void *ctx)
{
  const struct gk_app_ctx *app = ctx;
  double rminor = gkyl_gen_miller_eqdsk_rminor_at_R_omp(&app->miller, Rmajor_omp);
  if (!gkyl_is_dfinite(rminor)) {
    gkyl_exit("CBC: invalid outboard radius in safety-factor evaluation.");
  }
  return qsafety_at_rminor(rminor, app->rminor_LCFS);
}

static void
miller_kappa(double t, const double *xc, double *GKYL_RESTRICT out, void *ctx)
{
  const struct gk_app_ctx *app = ctx;
  out[0] = app->kappa;
}

static void
miller_delta(double t, const double *xc, double *GKYL_RESTRICT out, void *ctx)
{
  const struct gk_app_ctx *app = ctx;
  out[0] = app->delta;
}

static void
miller_a_shift(double t, const double *xc, double *GKYL_RESTRICT out, void *ctx)
{
  const struct gk_app_ctx *app = ctx;
  out[0] = app->a_shift;
}

/**
 * Save the app-owned host geometry before profile callbacks are evaluated.
 *
 * @param geo App's geometry object,
 * @param vctx CBC simulation context.
 */
static void
set_geo_ptr_in_ctx(const struct gk_geometry *geo, void *vctx)
{
  struct gk_app_ctx *ctx = vctx;
  ctx->geo = geo;
}

/**
 * Compute the Miller minor-radius offset r-r0 from the discrete spatial map.
 *
 * @param app Simulation context.
 * @param xc Gkeyll computational coordinates (psi,alpha,theta).
 * @return Minor-radius offset rminor-rminor0 [m], not a computational
 * coordinate.
 */
static double
rMr0_at_comp_coord(const struct gk_app_ctx *app, const double *xc)
{
  double miller_coords[3];
  gkyl_gk_geometry_geo_int_comp2miller_coords(app->geo, xc, &app->miller, miller_coords);
  return miller_coords[0] - app->rminor0;
}

/**
 * Reference exponential or tanh equilibrium profiles.
 * f0 exp(-(w/L) tanh(rMr0/w)): gradient 1/L at rMr0=0, flat at
 * the walls; or f0 exp(-rMr0/L) (constant gradient 1/L) if const_grad.
 *
 * @param rMr0 Miller minor-radius offset r-r0 [m].
 * @param f0 Profile value at rMr0=0.
 * @param Lgrad Gradient scale length [m].
 * @param width Tanh width [m].
 * @param const_grad Whether to use profile with constant gradient.
 * @return Equilibrium profile value.
 */
double
init_profile(double rMr0, double f0, double Lgrad, double width, bool const_grad)
{
  if (const_grad) {
    return f0 * exp(-rMr0 / Lgrad);
  }
  return f0 * exp(-(width / Lgrad) * tanh(rMr0 / width));
}

/**
 * Logarithmic radial derivative of the equilibrium profile.
 *
 * @param rMr0 Miller minor-radius offset r-r0 [m].
 * @param Lgrad Gradient scale length [m].
 * @param width Tanh width [m].
 * @param const_grad Whether to use profile with constant gradient.
 * @return Logarithmic gradient [1/m].
 */
double
dln_init_profile_drMr0(double rMr0, double Lgrad, double width, bool const_grad)
{
  if (const_grad) {
    return -1.0 / Lgrad;
  }
  return -1.0 / (Lgrad * pow(cosh(rMr0 / width), 2));
}

/**
 * Unperturbed density profile.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_density_eq(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double rMr0 = rMr0_at_comp_coord(app, xn); // Shifted Miller minor radius.
  fout[0] = init_profile(rMr0, app->n0, app->Ln, app->prof_width, app->const_grad);
}

/**
 * Perturbed ion density.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double alpha = xn[1];
  double rMr0 = rMr0_at_comp_coord(app, xn), ylength = app->Cy * alpha;
  double miller_coords[3];
  gkyl_gk_geometry_geo_int_comp2miller_coords(app->geo, xn, &app->miller, miller_coords);
  double thetaMiller = miller_coords[1];
  // Broadband seed: every ky mode of the box up to the Nyquist with
  // pseudo-random phases, localized in rMr0 and ballooning at the outboard
  // midplane (thetaMiller=0).
  double env = exp(-pow(rMr0 / app->pert_width, 2)) * exp(-pow(thetaMiller / (M_PI / 2.0), 2));
  double n_eq = init_profile(rMr0, app->n0, app->Ln, app->prof_width, app->const_grad);
  double pert = 0.0;
  for (int m = 1; m <= app->N_alpha; m++) {
    pert += cos(2.0 * M_PI * (m * ylength / app->L_ylength + fmod(0.6180339887 * m * m, 1.0)));
  }
  fout[0] = n_eq * (1.0 + app->delta_n * env * pert);
}

/**
 * Evaluate the reference canonical-Maxwellian parallel-flow correction.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  // To first order in rho_pol/L the canonical Maxwellian n(r*) F_M(T(r*)), with
  // r* = r + b_phi vpar/Omega_pol the orbit-centre radius (b_phi = -Bt/B here),
  // is a local Maxwellian with the parallel flow u = -(Ti/mi) (Bt/B)/Omega_pol
  // dln(p)/d(r-r0) (plus a heat flux we do not seed).
  double rMr0 = rMr0_at_comp_coord(app, xn);
  double Ti = init_profile(rMr0, app->Ti0, app->LTi, app->prof_width, app->const_grad);
  double dln_p_drMr0 = dln_init_profile_drMr0(rMr0, app->Ln, app->prof_width, app->const_grad) +
                       dln_init_profile_drMr0(rMr0, app->LTi, app->prof_width, app->const_grad);
  double btor[2], bmag;
  gkyl_gk_geometry_geo_int_bfieldtor_at_comp_coord(app->geo, xn, btor);
  gkyl_gk_geometry_geo_int_bmag_at_comp_coord(app->geo, xn, &bmag);
  double omega_pol = app->qi * btor[1] / app->mi;
  // The old Bt was a positive magnitude; btor[0] is the signed B_phi=-Bt.
  fout[0] = (Ti / app->mi) * (btor[0] / bmag) / omega_pol * dln_p_drMr0;
}

/**
 * Evaluate the ion temperature at the minor radius recovered from the discrete
 * map.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double rMr0 = rMr0_at_comp_coord(app, xn);
  fout[0] = init_profile(rMr0, app->Ti0, app->LTi, app->prof_width, app->const_grad);
}

/**
 * Set the adiabatic-electron parallel flow to zero.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

/**
 * Evaluate the reference electron temperature profile.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_temp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double rMr0 = rMr0_at_comp_coord(app, xn);
  fout[0] = init_profile(rMr0, app->Te0, app->LTi, app->prof_width,
                         app->const_grad); // LTe = LTi.
}

// Equilibrium Maxwellian (density, flow and temperature of the IC) for the
// Krook buffers.
/**
 * Evaluate the ion Krook Maxwellian using converted configuration coordinates.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_feq_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double rMr0 = rMr0_at_comp_coord(app, xn);
  double vpar = xn[3], mu = xn[4];
  double n = init_profile(rMr0, app->n0, app->Ln, app->prof_width, app->const_grad);
  double T = init_profile(rMr0, app->Ti0, app->LTi, app->prof_width, app->const_grad);
  double u;
  eval_upar_ion(t, xn, &u, ctx);
  double B;
  gkyl_gk_geometry_geo_int_bmag_at_comp_coord(app->geo, xn, &B);
  fout[0] =
    n / pow(2.0 * M_PI * T / app->mi, 1.5) * exp(-(0.5 * app->mi * pow(vpar - u, 2) + mu * B) / T);
}

// Krook rate: nu_krook in num_cell_krook_buff cells at each rMr0 wall, tanh
// transition (Grandgirard et al. 2016 Eq. 49).
/**
 * Evaluate Krook buffers as a function of r-r_0.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_krook_rate(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double rMr0 = rMr0_at_comp_coord(app, xn), Bs = 0.015 * app->L_rMr0,
         Bl = app->num_cell_krook_buff * app->L_rMr0 / app->N_psi;
  fout[0] = app->nu_krook * (1.0 + 0.5 * (tanh((rMr0 - app->rMr0_max + Bl) / Bs) -
                                          tanh((rMr0 - app->rMr0_min - Bl) / Bs)));
}

/**
 * Return the optional ion collision frequency.
 *
 * @param t Time [s].
 * @param xn Numerical (psi,alpha,theta) coordinates; velocity coordinates
 * follow for phase-space callbacks.
 * @param fout Output profile value.
 * @param ctx CBC simulation context.
 */
void
eval_nu_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

/**
 * Map computational velocity coordinates to parallel velocity and magnetic
 * moment.
 *
 * @param t Time [s], unused.
 * @param vc Computational velocity coordinates.
 * @param vp Output parallel velocity and magnetic moment.
 * @param ctx CBC simulation context.
 */
void
mapc2p_vel_ion(double t, const double *vc, double *GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  if (fabs(cvpar) <= 0.5) {
    vp[0] = vpar_max_ion * cvpar;
  } else if (cvpar < -0.5) {
    vp[0] = -vpar_max_ion * 2.0 * pow(cvpar, 2);
  } else {
    vp[0] = vpar_max_ion * 2.0 * pow(cvpar, 2);
  }
  vp[1] = mu_max_ion * pow(cmu, 2);
}

/**
 * Initialize CBC physics, profiles, domain sizes, and Miller inputs.
 *
 * @param ctx Caller-owned context; must remain at this address throughout the
 * simulation.
 */
void
create_ctx(struct gk_app_ctx *ctx)
{
  int cdim = 3, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field. (Taken from Greenfield et al. 1997 table II)
  double Rmajor_axis = 1.6714; // Magnetic axis major radius [m].
  double Zvert_axis = 0.0; // Magnetic axis vertical position [m].
  double kappa = 1.0; // Miller elongation.
  double delta = 0.0; // Miller triangularity.
  double a_shift = 0.0; // Shafranov-shift strength.
  double bmag_axis = 1.54; // Magnetic field at the magnetic axis [T].
  double rminor_LCFS = 0.604; // Miller surface label of the LCFS [m].
  double rminor0 = 0.5 * rminor_LCFS; // Reference Miller surface label [m].

  double mi_me_ratio = 100.0; // Reference ion-to-electron mass ratio (0 = hydrogen, mp/me).
  // rhostar_inv = rminor_LCFS/rho_s sets bmag_axis (0 = 1.54 T, i.e.
  // rminor_LCFS/rho_s = 203); box and profile widths are fixed in rho_s units.
  double rhostar_inv = 217.0; // Gives L_ylength = 120 rho_s with the 1/4 wedge (n_tor = 4).

  // Plasma parameters.
  double AMU = 1.0; // Hydrogen.
  double mi = mp * AMU;
  if (mi_me_ratio > 0.0) {
    me = mi / mi_me_ratio;
  }
  double Te0 = 2000 * eV;
  double Ti0 = 2000 * eV;
  double n0 = 4.5e19; // [1/m^3] according to Fig 5 of Greenfield et al. 1997.

  double vte = sqrt(Te0 / me), vti = sqrt(Ti0 / mi); // Thermal speeds.
  double c_s = sqrt(Te0 / mi);
  double omega_ci = fabs(qi * bmag_axis / mi);
  double rho_s = c_s / omega_ci;
  double rho_i = vti / omega_ci;
  if (rhostar_inv > 0.0) {
    rho_s = rminor_LCFS / rhostar_inv;
    bmag_axis = mi * c_s / (qi * rho_s);
    omega_ci = fabs(qi * bmag_axis / mi);
    rho_i = vti / omega_ci;
  }

  // Configuration domain parameters
  // Target L_ylength in rho_s; the Miller extent resolver quantizes the wedge.
  double L_ylength_rhos = 120.0;
  double L_rMr0_rhos = 96.0;
  double L_rMr0 = L_rMr0_rhos * rho_s; // Domain size along rMr0.
  double L_thetaMiller = 2. * M_PI - 1e-10; // Reference Miller-angle extent [rad].
  double rMr0_min = -L_rMr0 / 2;
  double rMr0_max = L_rMr0 / 2;
  double requested_L_ylength = L_ylength_rhos * rho_s;
  double ylength_min = -requested_L_ylength / 2;
  double ylength_max = requested_L_ylength / 2;
  double thetaMiller_min = -L_thetaMiller / 2;
  double thetaMiller_max = L_thetaMiller / 2;

  if (!gkyl_is_dfinite(rMr0_min) || !gkyl_is_dfinite(rMr0_max) || rminor0 + rMr0_min <= 0 ||
      rMr0_max <= rMr0_min || rminor0 + rMr0_max >= rminor_LCFS) {
    gkyl_exit("CBC: radial bounds must describe a positive-width core annulus "
              "inside the LCFS.");
  }

  // The Miller shape functions (e.g. kappa) expect the context.
  *ctx = (struct gk_app_ctx){
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .rminor_LCFS = rminor_LCFS,
    .kappa = kappa,
    .delta = delta,
    .a_shift = a_shift,
  };

  struct gkyl_gen_miller_eqdsk_inp miller = {
    .input_coord_limits_type =
      {GKYL_GK_GEOMETRY_COORD_RADIAL_R_OMP, GKYL_GK_GEOMETRY_COORD_BINORMAL_LENGTH,
       GKYL_GK_GEOMETRY_COORD_PARALLEL_THETA_MILLER},
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .bmag_axis = -bmag_axis, // Negative B on axis (ion gradB_vel down).
    .kappa = miller_kappa,
    .kappa_ctx = ctx,
    .delta = miller_delta,
    .delta_ctx = ctx,
    .a_shift = miller_a_shift,
    .a_shift_ctx = ctx,
    .qsafety = qsafety,
    .qsafety_ctx = ctx,
    .has_Rmajor_omp_bounds = true,
  };
  double Rmajor0 = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&miller, rminor0);
  double Rmajor_LCFS = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&miller, rminor_LCFS);
  double R_omp_min = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&miller, rminor0 + rMr0_min);
  double R_omp_max = gkyl_gen_miller_eqdsk_R_omp_at_rminor(&miller, rminor0 + rMr0_max);
  if (!gkyl_is_dfinite(Rmajor0) || !gkyl_is_dfinite(Rmajor_LCFS)) {
    gkyl_exit("CBC: invalid reference or LCFS Miller radius.");
  }
  if (!gkyl_is_dfinite(R_omp_min) || !gkyl_is_dfinite(R_omp_max) || R_omp_max <= R_omp_min) {
    gkyl_exit("CBC: outboard bounds are not increasing.");
  }
  miller.Rmajor_LCFS = Rmajor_LCFS;
  miller.Rmajor_omp_bounds[0] = R_omp_min;
  miller.Rmajor_omp_bounds[1] = R_omp_max;

  ctx->miller = miller;
  gkyl_gen_miller_eqdsk *equilibrium = gkyl_gen_miller_eqdsk_inew(&miller);
  if (!equilibrium) {
    gkyl_exit("CBC: failed to construct Miller equilibrium.");
  }

  double bmag0 = gkyl_gen_miller_eqdsk_bmag_at_RZ(equilibrium, Rmajor0, Zvert_axis);
  double psi_bounds[2] = {
    gkyl_gen_miller_eqdsk_psi_at_R_omp(equilibrium, R_omp_min),
    gkyl_gen_miller_eqdsk_psi_at_R_omp(equilibrium, R_omp_max)
  };
  double y_bounds[2] = {ylength_min, ylength_max};
  struct gkyl_gen_miller_eqdsk_binormal_extent binormal;
  if (!gkyl_gen_miller_eqdsk_set_binormal_extents(equilibrium, psi_bounds, y_bounds, &binormal)) {
    gkyl_gen_miller_eqdsk_release(equilibrium);
    gkyl_exit("CBC: failed to resolve binormal length bounds.");
  }

  double q0 = binormal.q0, Cy = binormal.Cy, L_ylength = binormal.Ly;
  int n_tor = binormal.n_tor;
  gkyl_gen_miller_eqdsk_release(equilibrium);

  // Constant-gradient profiles f0 exp(-rMr0/L) (R/LT, R/Ln everywhere in the
  // box) instead of tanh profiles whose gradient peaks at rMr0=0 and vanishes
  // at walls.
  bool const_grad = true;
  // Temperature gradient drive Rmajor_axis/LTi (0 = Cyclone 6.92);
  // Rmajor_axis/Ln stays at 2.22.
  double RoLT = 0.0;

  // Equilibrium gradients at rMr0=0 (Dimits et al. 2000) and ITG seed.
  // Cyclone R/LTi and R/Ln are defined with R = Rmajor_axis (the same R as eps
  // = rminor0/Rmajor_axis = 0.18).
  double kTi = RoLT > 0.0 ? RoLT : 6.92; // Rmajor_axis/LTi.
  double kn = 2.22; // Rmajor_axis/Ln.
  double LTi = Rmajor_axis / kTi;
  double Ln = Rmajor_axis / kn;
  double prof_width = L_rMr0 / 3.0; // Gradient > 0.6 of its peak over the
    // central 2/3 of the box (tanh only).
  double delta_n = 1.0e-6;
  double pert_width = L_rMr0 / 4.0;

  // Time unit. Linear rates in this setup: gamma = 0.20 c_s/Rmajor0 at ky*rho_s
  // = 0.35 (max ~0.23 near 0.45), so saturation is expected by t ~ 50-100
  // Rmajor0/c_s.
  double t_unit = Rmajor_axis / c_s; // Cyclone time unit R/c_s with R = Rmajor_axis.
  double inv_asp_ratio = rminor0 / Rmajor_axis;

  // Grid parameters
  int N_psi = 48 / 8; // L_rMr0/N_psi = 2 rho_s; the grid is uniform in psi.
  int N_alpha = 24 / 4; // Controls the alpha grid and reference seed spectrum.
  int N_theta = 16 / 4;
  int Nvpar = 8 / 2;
  int Nmu = 4 / 2;
  int poly_order = 1;

  // Velocity box dimensions
  double vpar_max_ion = 3. * vti * sqrt(2); // GENE normalizes with sqrt(2T/m)
  double mu_max_ion = 9. * Ti0 / bmag_axis; // GENE uses bmag_axis as reference field.
  double vpar_max_elc = 3. * vte * sqrt(2);
  double mu_max_elc = 9. * Te0 / bmag_axis;
  // Optional ion-ion LBO collisions, nu_star = nu_ii q Rmajor_axis/(eps^{3/2}
  // v_ti) with v_ti = sqrt(Ti/mi).
  double nu_star = 0.0;
  double nu_ion = nu_star * pow(inv_asp_ratio, 1.5) * vti / (q0 * Rmajor_axis);

  // Preserve Krook widths num_cell_buff*L_rMr0/N_psi.
  // This is a reference radial scale, not the spacing of the uniform-psi grid.
  // CLI resolution overrides preserve the buffer widths and seed spectrum.
  double nu_krook_norm = 5.0; // In units of c_s/Rmajor_axis.
  int num_cell_krook_buff = 4;
  double nu_krook = nu_krook_norm / t_unit;

  double t_end = 10. * t_unit;
  int num_frames = 100;
  double write_phase_freq = 0.5; // f every 2 frames for restarts across 12 h jobs.
  int int_diag_calc_num = num_frames * 100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  printf(
    "a/rho_s = %g, bmag_axis = %g T, bmag0 = %g T, n_tor = %d, L_ylength = %g rho_s "
    "(ky*rho_s "
    "multiples of %.3f, max %.2f), L_rMr0 = %g rho_s, Rmajor_axis/c_s = "
    "%g s\n",
    rminor_LCFS / rho_s, bmag_axis, bmag0, n_tor, L_ylength / rho_s, 2. * M_PI * rho_s / L_ylength,
    M_PI * N_alpha * rho_s / L_ylength, L_rMr0 / rho_s, t_unit
  );

  *ctx = (struct gk_app_ctx){
    .cdim = cdim,
    .vdim = vdim,
    .Rmajor_axis = Rmajor_axis,
    .Zvert_axis = Zvert_axis,
    .kappa = kappa,
    .delta = delta,
    .a_shift = a_shift,
    .Rmajor0 = Rmajor0,
    .Rmajor_LCFS = Rmajor_LCFS,
    .R_omp_min = R_omp_min,
    .R_omp_max = R_omp_max,
    .ylength_min = ylength_min,
    .ylength_max = ylength_max,
    .rminor_LCFS = rminor_LCFS,
    .rminor0 = rminor0,
    .bmag_axis = bmag_axis,
    .bmag0 = bmag0,
    .q0 = q0,
    .Cy = Cy,
    .L_rMr0 = L_rMr0,
    .L_ylength = L_ylength,
    .const_grad = const_grad,
    .Ln = Ln,
    .LTi = LTi,
    .prof_width = prof_width,
    .delta_n = delta_n,
    .pert_width = pert_width,
    .n_tor = n_tor,
    .nu_krook = nu_krook,
    .nu_krook_norm = nu_krook_norm,
    .num_cell_krook_buff = num_cell_krook_buff,
    .rMr0_min = rMr0_min,
    .rMr0_max = rMr0_max,
    .thetaMiller_min = thetaMiller_min,
    .thetaMiller_max = thetaMiller_max,
    .me = me,
    .qe = qe,
    .mi = mi,
    .qi = qi,
    .n0 = n0,
    .Te0 = Te0,
    .Ti0 = Ti0,
    .c_s = c_s,
    .rho_s = rho_s,
    .N_psi = N_psi,
    .N_alpha = N_alpha,
    .N_theta = N_theta,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {N_psi, N_alpha, N_theta, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .nu_star = nu_star,
    .nu_ion = nu_ion,
    .write_phase_freq = write_phase_freq,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .miller = miller,
  };
}

/**
 * Run the numerical Miller CBC case through the shared simulation driver.
 *
 * @param argc Number of command-line arguments.
 * @param argv Resolution, parallelism, restart, and step-limit arguments.
 * @return Zero after the shared driver returns.
 */
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

  struct gk_app_ctx ctx; // Context for initialization and profile callbacks.
  create_ctx(&ctx);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d = 0; d < ctx.cdim; d++) {
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  }
  for (int d = 0; d < ctx.vdim; d++) {
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim + d]);
  }
  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);
  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    gkyl_comm_get_rank(comm, &my_rank);
  }
#endif

  if (my_rank == 0) {
    printf(
      "CBC Miller: R_omp=[%.12g, %.12g] m, B_phi(axis)=%.12g T, requested y "
      "width=%.12g m.\n",
      ctx.R_omp_min, ctx.R_omp_max, ctx.miller.bmag_axis, ctx.ylength_max - ctx.ylength_min
    );
  }

  // ions
  struct gkyl_gyrokinetic_projection ion_bc = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .density = eval_density_eq,
    .ctx_density = &ctx,
    .upar = eval_upar_ion,
    .ctx_upar = &ctx,
    .temp = eval_temp_ion,
    .ctx_temp = &ctx,
    .correct_all_moms = true,
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-1.0 / sqrt(2.0), 0.0},
    .upper = {1.0 / sqrt(2.0), 1.0},
    .mapc2p.mapping = mapc2p_vel_ion,
    .mapc2p.ctx = &ctx,
    .cells = {cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,

    // Moment-corrected Maxwellian so the ion charge matches the electron
    // density profile.
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_density_ion,
        .ctx_density = &ctx,
        .upar = eval_upar_ion,
        .ctx_upar = &ctx,
        .temp = eval_temp_ion,
        .ctx_temp = &ctx,
        .correct_all_moms = true,
      },

    .collisionless = {.type = GKYL_GK_COLLISIONLESS_ES},

    // Fixed Maxwellian ghost cells: the walls act as a reservoir for the
    // radial magnetic drift (zero-flux and copy BCs build up boundary
    // layers).
    .bcs =
      {{
         .dir = 0,
         .edge = GKYL_LOWER_EDGE,
         .type = GKYL_BC_GK_SPECIES_FIXED_FUNC,
         .projection = ion_bc,
       },
       {
         .dir = 0,
         .edge = GKYL_UPPER_EDGE,
         .type = GKYL_BC_GK_SPECIES_FIXED_FUNC,
         .projection = ion_bc,
       },
       {.dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT},
       {.dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_TWISTSHIFT}},

    // Krook buffers at the radial boundaries.
    .source_bgk =
      {
        .source_bgk_id = GKYL_SOURCE_BGK_STATIC,
        .rate_profile = eval_krook_rate,
        .rate_profile_ctx = &ctx,
        .feq_shape = eval_feq_ion,
        .feq_shape_ctx = &ctx,
      },

    .num_diag_moments = 4,
    .diag_moments =
      {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
    .boundary_flux_diagnostics =
      {.num_integrated_diag_moments = 1, .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN}},
    .time_rate_diagnostics = true,
  };

  // Adiabatic electrons: not evolved, rescaled to n0*(1 + e*(phi-<phi>)/Te0)
  // every step.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .projection =
      {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_density_eq, // Match the ion equilibrium density.
        .ctx_density = &ctx,
        .upar = eval_upar_elc,
        .ctx_upar = &ctx,
        .temp = eval_temp_elc,
        .ctx_temp = &ctx,
        .correct_all_moms = true,
      },
    .scaling = {.type = GKYL_GK_SPECIES_SCALING_ADIABATIC},
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  if (ctx.nu_star > 0.0) {
    ion.collisions.collision_id = GKYL_LBO_COLLISIONS;
    ion.collisions.self_nu = eval_nu_ion;
    ion.collisions.self_nu_ctx = &ctx;
  }

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx.bmag_axis,
    .poisson_bcs =
      {{.dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}},
       {.dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}}},
    .time_rate_diagnostics = true,
  };

  // The default Miller full-contour mapping also covers this entirely closed
  // core annulus. Numerical TOKAMAK geometry supplies the two core twist
  // shifts.
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MILLER,
    .miller_grid_info = ctx.miller,
    .set_gk_geometry = set_geo_ptr_in_ctx,
    .set_gk_geometry_ctx = &ctx,
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1], app_args.cuts[2]},
    .use_gpu = app_args.use_gpu,
  };

  struct gkyl_msgpack_map_elem additional_metadata[] = {
    {.key = "Cy", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Cy},
    {.key = "Rmajor_axis", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Rmajor_axis},
    {.key = "rminor_LCFS", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.rminor_LCFS},
    {.key = "Te0", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.Te0},
    {.key = "c_s", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.c_s},
    {.key = "rho_s", .elem_type = GKYL_MP_DOUBLE, .dval = ctx.rho_s}
  };

  // GK app
  struct gkyl_gk app_inp = {
    .cfl_frac_omegaH = 1.7,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = {ctx.R_omp_min, ctx.ylength_min, ctx.thetaMiller_min},
    .upper = {ctx.R_omp_max, ctx.ylength_max, ctx.thetaMiller_max},
    .cells = {cells_x[0], cells_x[1], cells_x[2]},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = geometry,

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2, // Adiabatic electrons are a species (scaling ADIABATIC) on main;
    // dropping them gives an ion-only Poisson solve.
    .species = {ion, elc},

    .field = field,

    .parallelism = parallelism,

    .metadata =
      {
        .num_attributes = sizeof(additional_metadata) / sizeof(additional_metadata[0]),
        .attributes = additional_metadata,
      },
  };

  // Set app output name from the executable name (argv[0]).
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
    .print_verbosity = {.enabled = true, .frequency = 0.01},
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
