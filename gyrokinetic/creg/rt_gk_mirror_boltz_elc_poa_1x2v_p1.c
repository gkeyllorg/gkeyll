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
  enum gkyl_gyrokinetic_damping_type damping_type; // Type of damping.
};

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.

  // Plasma parameters
  double mi; // Ion mass.
  double me; // Electron mass.
  double qi; // Ion charge.
  double qe; // Electron charge.
  double Te0; // Electron temperature.
  double Ti0; // Ion temperature.
  double n0; // Density.
  double B_p; // Plasma magnetic field (mirror center).
  double beta; // Plasma beta in the center.
  double tau; // Temperature ratio.

  // Parameters controlling initial conditions.
  double alim;
  double alphaIC0;
  double alphaIC1;
  double Ti_perp0; // Reference ion perp temperature.
  double Ti_par0; // Reference ion par temperature.
  double Ti_perp_m; // Ion perp temperature at the throat.
  double Ti_par_m; // Ion par temperature at the throat.
  double cs_m; // Ion sound speed at the throat.

  double nuFrac; // Fraction multiplying collision frequency.
  double logLambdaIon; // Ion Coulomb logarithm.
  double nuIon; // Ion-ion collision freq.

  double vti; // Ion thermal speed.
  double vte; // Electron thermal speed.
  double c_s; // Ion sound speed.
  double omega_ci; // Ion gyrofrequency.
  double rho_s; // Ion sound gyroradius.

  double RatZeq0; // Radius of the field line at Z=0.
  double Z_min; // Minimum axial coordinate Z.
  double Z_max; // Maximum axial coordinate Z.
  double z_min; // Minimum value of the position along the field line.
  double z_max; // Maximum value of the position along the field line.
  double psi_eval; // Psi (poloidal flux) of the field line.
  double psi_in, z_in; // Auxiliary psi and z.

  // Magnetic equilibrium model.
  double mcB;
  double gamma;
  double Z_m; // Axial coordinate at mirror throat.
  double z_m; // Computational coordinate at mirror throat.

  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;

  // Physical velocity space limits.
  double vpar_min_ion, vpar_max_ion;
  double mu_max_ion;
  // Computational velocity space limits.
  double vpar_lin_fac_inv, mu_lin_fac_inv; // Inverse factor of where linear mapping ends.
  double vpar_pow, mu_pow; // Power of the velocity grid.
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  // Grid DOF.
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  int num_phases; // Number of phases.
  struct gk_poa_phase_params *poa_phases; // Phases to run.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
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
  double Rout = sqrt(2.0 * psiIn / (app->mcB * 
    (1.0 / (M_PI * app->gamma * (1.0 + pow((ZIn - app->Z_m) / app->gamma, 2.))) +
     1.0 / (M_PI * app->gamma * (1.0 + pow((ZIn + app->Z_m) / app->gamma, 2.))))));
  return Rout;
}

void
Bfield_psiZ(double psiIn, double ZIn, void *ctx, double *BRad, double *BZ, double *Bmag)
{
  struct gk_mirror_ctx *app = ctx;
  double Rcoord = R_psiZ(psiIn, ZIn, ctx);
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;
  *BRad = -(1.0 / 2.0) * Rcoord * mcB *
          (-2.0 * (ZIn - Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn - Z_m) / gamma, 2.), 2.))) -
            2.0 * (ZIn + Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn + Z_m) / gamma, 2.), 2.))));
  *BZ = mcB *
        (1.0 / (M_PI * gamma * (1.0 + pow((ZIn - Z_m) / gamma, 2.))) +
         1.0 / (M_PI * gamma * (1.0 + pow((ZIn + Z_m) / gamma, 2.))));
  *Bmag = sqrt(pow(*BRad, 2) + pow(*BZ, 2));
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
  app->psi_in = psiIn;
  double eps = 0.0;
  struct gkyl_qr_res integral;
  if (eps <= ZIn)
  {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, eps, ZIn, 7, 1e-14);
  }
  else
  {
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
  if (zIn >= 0.0)
  {
    double fl = root_Z_psiz(-eps, ctx);
    double fr = root_Z_psiz(app->Z_max + eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, -eps, app->Z_max + eps, fl, fr, 1000, 1e-14);
  }
  else
  {
    double fl = root_Z_psiz(app->Z_min - eps, ctx);
    double fr = root_Z_psiz(eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, app->Z_min - eps, eps, fl, fr, 1000, 1e-14);
  }
  return Zout.res;
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;

  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx); // Cylindrical axial coordinate.

  if (fabs(Z) <= app->Z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                              exp(-pow(z - zSrc, 2) / (2.0 * pow(sigSrc, 2))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double sigSrc = app->sigSrcIon;
  double TSrc0 = app->TSrc0Ion;
  double Tfloor = app->TSrcFloorIon;

  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = (fabs(z) / z) * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m));
  }
}

void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double Ti_par0 = app->Ti_par0;
  double Ti_par_m = app->Ti_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_par_m+(Ti_par0-Ti_par_m)*tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_par_m;
  }
}

void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double Ti_perp0 = app->Ti_perp0;
  double Ti_perp_m = app->Ti_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_perp_m - Ti_perp0*tanh(3.*fabs(z_m-fabs(z)));
  }
  else
  {
    fout[0] = Ti_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

// Geometry evaluation functions for the gk app
// mapc2p must assume a 3d input xc
void
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  double psi = xc[0];
  double theta = xc[1];
  double z = xc[2];

  double Z = Z_psiz(psi, z, ctx);
  double R = R_psiZ(psi, Z, ctx);

  // Cartesian coordinates on plane perpendicular to Z axis.
  double x = R * cos(theta);
  double y = R * sin(theta);
  xp[0] = x;
  xp[1] = y;
  xp[2] = Z;
}

// bmag_func must assume a 3d input xc
void
bmag_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[2];

  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx);
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);
  fout[0] = Bmag;
}

// bfield_func must assume a 3d input xc
void
bfield_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[2];

  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx);
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);

  double phi = xc[1];
  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = BRad*cos(phi);
  fout[1] = BRad*sin(phi);
  fout[2] = BZ;
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/lin_frac_inv, then a power grid.
  double vpar_lin_fac_inv = app->vpar_lin_fac_inv;
  double vpar_pow = app->vpar_pow;
  if (fabs(cvpar) <= 1.0/vpar_lin_fac_inv)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -1.0/vpar_lin_fac_inv)
    vp[0] = -vpar_max_ion*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);
  else
    vp[0] =  vpar_max_ion*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);

//  // Quadratic mu.
//  vp[1] = mu_max_ion*pow(cmu,2.0);
  // Linear map up to mu_max/lin_frac_inv, then a power grid.
  double mu_lin_fac_inv = app->mu_lin_fac_inv;
  double mu_pow = app->mu_pow;
//  if (cmu <= 1.0/mu_lin_fac_inv)
//    vp[0] = mu_max_ion*cmu;
//  else
//    vp[0] = mu_max_ion*pow(mu_lin_fac_inv,mu_pow-1)*pow(cmu,mu_pow);
  double w = 0.3;
  double f = 0.012;
  double a = mu_max_ion*(f-1.0)/(w*w-1.0);
  double b = mu_max_ion*(w*w-f)/(w*w-1.0);
  if (cmu <= w)
    vp[1] = (f*mu_max_ion/w)*cmu;
  else
    vp[1] = a*pow(cmu,2)+b;

}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

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

  // Parameters controlling initial conditions.
  double alim = 0.125;
  double alphaIC0 = 2;
  double alphaIC1 = 10;

  double nuFrac = 1.0;
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Geometry parameters.
  double RatZeq0 = 0.10; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  // the boundary of a cell (due to AD errors).
  double Z_min = -2.5;
  double Z_max =  2.5;

  // Parameters controlling the magnetic equilibrium model.
  double mcB = 6.51292;
  double gamma = 0.124904;
  double Z_m = 0.98;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0 / 40.0 / 2.0 * 1.25;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = Z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;

  // Grid parameters
  double vpar_max_ion = 16 * vti;
  double vpar_min_ion = -vpar_max_ion;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);

  // Computational velocity space limits.
  double vpar_lin_fac_inv = 4;
  double vpar_pow = 3;
  double vpar_min_ion_c = -1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double vpar_max_ion_c =  1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  double mu_lin_fac_inv = 1.0/0.012;
  double mu_pow = 2;
//  double mu_min_ion_c = 0.0;
//  double mu_max_ion_c = 1.0/pow(mu_lin_fac_inv,(mu_pow-1)/mu_pow);

  // Grid DOF:
  int Nz = 64; // Number of cells in z direction.
  int Nvpar = 32; // Number of cells in parallel velocity direction.
  int Nmu = 16;  // Number of cells in mu direction.
  int poly_order = 1;

  // Initial conditions parameter.s
  double Ti_perp0 = 10000 * eV;
  double Ti_par0 = 7500 * eV;

  // Parameters at mirror throat
  double Ti_perp_m = 15000 * eV;
  double Ti_par_m = 1000 * eV;
  double z_m = 0.982544;
  double cs_m = sqrt((Te0+3.0*Ti_par_m)/mi);

  // Factor multiplying collisionless terms.
  double alpha_oap = 0.01;
  double alpha_fdp = 1.0;
  // Duration of each phase.
  double tau_oap = 2400.0e-9;
  double tau_fdp = 24.0e-9;
  double tau_fdp_extra = 2*tau_fdp;
  int num_cycles = 2; // Number of OAP+FDP cycles to run.

  // Frame counts for each phase type (specified independently)
  int num_frames_oap = 1; // Frames per OAP phase
  int num_frames_fdp = 1; // Frames per FDP phase
  int num_frames_fdp_extra = 2*num_frames_fdp;  // Frames for the extra FDP phase

  // Whether to evolve the field.
  bool is_static_field_oap = true;
  bool is_static_field_fdp = false;
  // Whether to enable positivity.
  bool is_positivity_enabled_oap = false;
  bool is_positivity_enabled_fdp = true;
  // Type of df/dt multipler.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_oap = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE;
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_fdp = GKYL_GK_FDOT_MULTIPLIER_NONE;

  // Calculate phase structure
  double t_end = (tau_oap + tau_fdp)*num_cycles + tau_fdp_extra;
  double tau_pair = tau_oap+tau_fdp; // Duration of an OAP+FDP pair.
  int num_phases = 2*num_cycles + 1;
  int num_frames = num_cycles * (num_frames_oap + num_frames_fdp) + num_frames_fdp_extra;

  struct gk_poa_phase_params *poa_phases = gkyl_malloc(num_phases * sizeof(struct gk_poa_phase_params));
  for (int i=0; i<(num_phases-1)/2; i++) {
    // OAPs.
    poa_phases[2*i].phase = GK_POA_OAP;
    poa_phases[2*i].num_frames = num_frames_oap;
    poa_phases[2*i].duration = tau_oap;
    poa_phases[2*i].alpha = alpha_oap;
    poa_phases[2*i].is_static_field = is_static_field_oap;
    poa_phases[2*i].fdot_mult_type = fdot_mult_type_oap;
    poa_phases[2*i].is_positivity_enabled = is_positivity_enabled_oap;
    poa_phases[2*i].damping_type = GKYL_GK_DAMPING_NONE;

    // FDPs.
    poa_phases[2*i+1].phase = GK_POA_FDP;
    poa_phases[2*i+1].num_frames = num_frames_fdp;
    poa_phases[2*i+1].duration = tau_fdp;
    poa_phases[2*i+1].alpha = alpha_fdp;
    poa_phases[2*i+1].is_static_field = is_static_field_fdp;
    poa_phases[2*i+1].fdot_mult_type = fdot_mult_type_fdp;
    poa_phases[2*i+1].is_positivity_enabled = is_positivity_enabled_fdp;
    poa_phases[2*i+1].damping_type = GKYL_GK_DAMPING_LOW_PASS_FILTER;
  }
  // Add an extra, longer FDP.
  poa_phases[num_phases-1].phase = GK_POA_FDP;
  poa_phases[num_phases-1].num_frames = num_frames_fdp_extra;
  poa_phases[num_phases-1].duration = tau_fdp_extra;
  poa_phases[num_phases-1].alpha = alpha_fdp;
  poa_phases[num_phases-1].is_static_field = is_static_field_fdp;
  poa_phases[num_phases-1].fdot_mult_type = fdot_mult_type_fdp;
  poa_phases[num_phases-1].is_positivity_enabled = is_positivity_enabled_fdp;
  poa_phases[num_phases-1].damping_type = GKYL_GK_DAMPING_LOW_PASS_FILTER;

  double write_phase_freq = 0.5; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq = 5; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,  .vdim = vdim,
    .mi = mi,  .qi = qi,
    .me = me,  .qe = qe,
    .Te0 = Te0,  .Ti0 = Ti0,  .n0 = n0,
    .B_p = B_p,  .beta = beta,  .tau = tau,
    .alim = alim,
    .alphaIC0 = alphaIC0,
    .alphaIC1 = alphaIC1,
    .nuFrac = nuFrac,  .logLambdaIon = logLambdaIon,  .nuIon = nuIon,
    .vti = vti,  .vte = vte,  .c_s = c_s,
    .omega_ci = omega_ci,  .rho_s = rho_s,
    .RatZeq0 = RatZeq0,
    .Z_min = Z_min,  .Z_max = Z_max,
    // Parameters controlling the magnetic equilibrium model.
    .mcB = mcB,  .gamma = gamma,
    .Z_m = Z_m,
    .z_m = z_m,
    // Initial condition parameters.
    .Ti_perp0 = Ti_perp0,  .Ti_par0 = Ti_par0,
    .Ti_perp_m = Ti_perp_m,  .Ti_par_m = Ti_par_m,  .cs_m = cs_m,
    // Source parameters
    .NSrcIon = NSrcIon,  .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,  .TSrcFloorIon = TSrcFloorIon,
    .lineLengthSrcIon = lineLengthSrcIon,  .sigSrcIon = sigSrcIon,
    // Physical velocity space limits.
    .vpar_min_ion = vpar_min_ion,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    // Computational velocity space limits.
    .vpar_lin_fac_inv = vpar_lin_fac_inv,
    .vpar_pow = vpar_pow,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_lin_fac_inv = mu_lin_fac_inv,
    .mu_pow = mu_pow,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    // Grid DOF.
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    // Time integration and I/O parameters.
    .t_end = t_end,
    .num_frames = num_frames,
    .num_phases = num_phases,
    .poa_phases = poa_phases,
    .write_phase_freq     = write_phase_freq    , 
    .int_diag_calc_freq   = int_diag_calc_freq  , 
    .dt_failure_tol       = dt_failure_tol      , 
    .num_failures_max     = num_failures_max    , 
  };

  // Populate a couple more values in the context.
  ctx.psi_eval = psi_RZ(ctx.RatZeq0, 0., &ctx);
  ctx.z_min    = z_psiZ(ctx.psi_eval, ctx.Z_min, &ctx);
  ctx.z_max    = z_psiZ(ctx.psi_eval, ctx.Z_max, &ctx);

  return ctx;
}

void
release_ctx(struct gk_mirror_ctx *ctx)
{
  gkyl_free(ctx->poa_phases);
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool force_calc, double dt)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
    gkyl_gyrokinetic_app_write_dt(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

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
  int num_int_diag_calc = ctx->int_diag_calc_freq*num_frames;

  // Prevent division by zero when frame_curr equals num_frames
  int frames_remaining = num_frames - frame_curr;
  double time_remaining = t_end - t_curr;

  trig_write_conf->dt = time_remaining / frames_remaining;
  trig_write_conf->tcurr = t_curr;
  trig_write_conf->curr = frame_curr;

  trig_write_phase->dt = time_remaining / (ctx->write_phase_freq * frames_remaining);
  trig_write_phase->tcurr = t_curr;
  trig_write_phase->curr = frame_curr;

  int diag_frames = GKYL_MAX2(frames_remaining, (num_int_diag_calc/num_frames) * frames_remaining);
  trig_calc_intdiag->dt = time_remaining / diag_frames;
  trig_calc_intdiag->tcurr = t_curr;
  trig_calc_intdiag->curr = frame_curr;
}

void run_phase(gkyl_gyrokinetic_app* app, struct gk_mirror_ctx *ctx, double num_steps,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag,  struct time_frame_state *tfs,
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
    .type = pparams->fdot_mult_type,
    .cellwise_const = true,
    .write_diagnostics = true,
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
  struct gkyl_gyrokinetic_damping damping_inp = {
    .type = pparams->damping_type,
    .cellwise_const = false,
    .rate_const = 1/5e-6,
    .write_fbar = true,
    .write_rate = true,
    .do_not_reset_fbar = true,
  };

  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "ion", collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "ion", fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "ion", positivity_inp);
  gkyl_gyrokinetic_app_reset_species_damping(app, t_curr, "ion", damping_inp);
  gkyl_gyrokinetic_app_reset_field(app, t_curr, field_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps))
  {
    if (step == 1 || step % 20 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);

    dt = fmin(dt, t_end - t_curr); // Don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step == 1 || step % 20 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success)
    {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, t_curr > t_end, status.dt_actual);
    write_data(trig_write_conf, trig_write_phase, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
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
  tfs->frame_curr = tfs->frame_curr+pparams->num_frames;
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
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,  .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = eval_density_ion,
      .upar = eval_upar_ion,
      .temppar = eval_temp_par_ion,
      .tempperp = eval_temp_perp_ion,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temppar = &ctx,
      .ctx_tempperp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
      .scale_factor = 1.0, // Will be replaced below.
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalNuIon,
      .self_nu_ctx = &ctx,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = eval_density_ion_source,
        .upar = eval_upar_ion_source,
        .temp = eval_temp_ion_source,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
      }, 
    },

    .time_rate_multiplier = {
      .type = GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE, // So solvers are allocated.
      .cellwise_const = true,
      .write_diagnostics = true,
    },

    .positivity = {
      .type = GKYL_GK_POSITIVITY_SHIFT,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 4,
    .diag_moments = {GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN},
  };

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .is_static = false, // So solvers are allocated.
  };

  // GK app
  struct gkyl_gk app_inp = { 
    .cdim = ctx.cdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {ctx.psi_eval, 0.0},
      .mapc2p = mapc2p, // Mapping of computational to physical space.
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // Magnetic field.
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {ion},

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
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
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    tfs.frame_curr = status.frame;
    tfs.t_curr = status.stime;

    // Find out what phase we are in.
    double time_count = 0.0;
    int frame_count = 0;
    int pit_curr = 0;
    for (int pit=0; pit<ctx.num_phases; pit++) {
      time_count += ctx.poa_phases[pit].duration;
      frame_count += ctx.poa_phases[pit].num_frames;
      if ((tfs.t_curr <= time_count) && (tfs.frame_curr <= frame_count)) {
        pit_curr = pit;
        break;
      }
    };
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

  // if (app_args.num_steps != INT_MAX)
  //   phase_idx_end = 1;

  // Loop over number of number of phases;
  for (int pit=phase_idx_init; pit<phase_idx_end; pit++) {
    gkyl_gyrokinetic_app_cout(app, stdout, "\nRunning phase %d @ t = %.9e ... \n", pit, tfs.t_curr);
    struct gk_poa_phase_params *phase_params = &ctx.poa_phases[pit];
    run_phase(app, &ctx, app_args.num_steps, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag, &tfs, phase_params);
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0)
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
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
