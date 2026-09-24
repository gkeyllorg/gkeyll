#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double a_shift, Z_axis, R_axis, R0, a_mid, x_inner, r0, B0, kappa, delta, q0, Cy, Bref, x_LCFS;
  // Plasma parameters
  double me, qe, mi, qi, n0, Te0, Ti0;
  // Collision parameters
  double nuFrac, nuElc, nuIon;
  // Source parameters
  int num_sources;
  bool adapt_energy_srcCORE, adapt_particle_srcCORE; 
  double center_srcCORE[2], sigma_srcCORE[2];
  double energy_srcCORE, particle_srcCORE;
  double floor_srcCORE;
  bool adapt_energy_srcRECY, adapt_particle_srcRECY;
  double adapt_particle_fraction_srcRECY;
  double center_srcRECY[2], sigma_srcRECY[2];
  double energy_srcRECY, particle_srcRECY;
  double floor_srcRECY;
  // Grid parameters
  double Lx, Lz;
  double x_min, x_max, z_min, z_max;
  int Nx, Nz, Nvpar, Nmu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
};

// Geometry related functions 
double r_x(double x, double a_mid, double x_inner)
{
  return x+a_mid-x_inner;
}

// cubic polynomial fit to TCV NT q profile (discharge #65130)
double qprofile(double R) {
  double qfit[4] = {
    497.3420166252413, -1408.736172826569, 1331.4134861681464, -419.00692601227627
  };
  return qfit[0]*pow(R,3) + qfit[1]*pow(R,2) + qfit[2]*R + qfit[3];
}

double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return R_axis - a_shift*r*r/(2.*R_axis) + r*cos(theta + asin(delta)*sin(theta));
}

double Z_rtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Z_axis = app->Z_axis;
  double kappa = app->kappa;
  return Z_axis + kappa*r*sin(theta);
}

double dRdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double a_shift = app->a_shift;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return - a_shift*r/(R_axis) + cos(theta + asin(delta)*sin(theta));
}

double dRdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return -r*sin(theta + asin(delta)*sin(theta))*(1.+asin(delta)*cos(theta));
}

double dZdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*sin(theta);
}

double dZdtheta(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*cos(theta);
}

double Jr(double r, double theta, void *ctx)
{
  return R_rtheta(r,theta,ctx)*( dRdr(r,theta,ctx) * dZdtheta(r,theta,ctx)
		                -dRdtheta(r,theta,ctx) * dZdr(r,theta,ctx) );
}

struct integrand_ctx {
  struct gk_app_ctx *app_ctx;
  double r;
};

double integrand(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}

double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);
  double R_omp = R_rtheta(r,0.0,ctx);
  return ( R*Bt/(2.*M_PI*qprofile(R_omp)))*integral.res;
}

double alpha(double r, double theta, double phi, void *ctx)
{
  double twrap = theta;
  while (twrap < -M_PI) twrap = twrap+2.*M_PI;
  while (M_PI < twrap) twrap = twrap-2.*M_PI;

  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., twrap, 7, 1e-10);
  } else {
    integral = gkyl_dbl_exp(integrand, &tmp_ctx, twrap, 0., 7, 1e-10);
    integral.res = -integral.res;
  }

  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);

  return phi - R*Bt*integral.res/dPsidr(r,theta,ctx);
}

double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Density initial condition (like TCV exp profile)
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_app_ctx *app = ctx;
  double n0 = 5e19;
  double x0 = -0.03;
  double c1 = 0.5;
  double c2 = 8.0;
  double c3 = 0.005;
  fout[0] = n0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

// Electron temperature initial conditions
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_app_ctx *app = ctx;
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.03; // position of the transition region
  double c0 = 1.3; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 8.0; // control the width of the transition region
  double c3 = 0.1; // control the temperature at the SOL
  fout[0] = c0*T0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

// Ion temperature initial conditions
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_app_ctx *app = ctx;
  double T0 = 200 * GKYL_ELEMENTARY_CHARGE;
  double x0 = -0.04; // position of the transition region
  double c0 = 1.0; // multiplicative factor
  double c1 = 0.5; // control the temperature at the _core
  double c2 = 3.0; // control the width of the transition region
  double c3 = 0.2; // control the temperature at the SOL
  fout[0] = c0*T0*(c1*(1.+tanh(c2*(-10*(x+x0))))+c3);
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  fout[0] = 0.5; // Diffusivity [m^2/s].
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double Cy = app->Cy;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  // Map to cylindrical (R, Z, phi) coordinates.
  double R   = R_rtheta(r, z, ctx);
  double Z   = Z_rtheta(r, z, ctx);
  double phi = y/Cy + alpha(r, z, 0, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);
  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

// Taken from rt gk d3d 3x2c, is this the non uniform v grid mapping?
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double x_inner = app->x_inner;
  double r = r_x(x,a_mid,x_inner);
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);

  double drdtheta = dRdtheta(r,z,ctx);
  double dzdtheta = dZdtheta(r,z,ctx);
  double den = sqrt(pow(drdtheta,2) + pow(dzdtheta,2));
  double B_r = Bp*drdtheta/den;
  double B_z = Bp*dzdtheta/den;
  double phi = y/Cy + alpha(r, z, 0, ctx);
  double R   = R_rtheta(r, z, ctx);

  // xc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = B_r * cos(phi) + Bt * sin(phi);
  fout[1] = B_r * sin(phi) - Bt * cos(phi);
  fout[2] = B_z;
}


struct gk_app_ctx create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double a_shift = 0.5; // Parameter in Shafranov shift.
  double Z_axis = 0.1414361745; // Magnetic axis height [m].
  double R_axis = 0.8867856264; // Magnetic axis major radius [m].
  double B_axis = 1.4; // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 1.0870056099999; // Major radius of the LCFS at the outboard midplane [m
  double x_inner = 0.04; // Radial extent inside LCFS    
  double x_outer = 0.08; // Radial extent outside LCFS
  double Rmid_min = R_LCFSmid - x_inner; // Minimum midplane major radius of simulation box [m].
  double Rmid_max = R_LCFSmid + x_outer; // Maximum midplane major radius of simulation box [m].
  double R0 = 0.5*(Rmid_min+Rmid_max); // Major radius of the simulation box [m].

  // Minor radius at outboard midplane [m]. Redefine it with
  // Shafranov shift, to ensure LCFS radial location.
  double a_mid = fabs(a_shift)<1e-13? R_LCFSmid-R_axis :
    R_axis/a_shift - sqrt(R_axis*(R_axis - 2*a_shift*R_LCFSmid + 2*a_shift*R_axis))/a_shift;

  double r0 = R0-R_axis; // Minor radius of the simulation box [m].
  double B0 = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa = 1.4; // Elongation (=1 for no elongation).
  double delta = -0.38; // Triangularity (=0 for no triangularity).

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double AMU = 2.01410177811;
  double mi = mp*AMU; // Deuterium ions.
  double Te0 = 100*eV;
  double Ti0 = 100*eV;
  double n0 = 2.0e19; // [1/m^3]
  double Bref = 1.129; // Reference magnetic field [T].
  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Configuration domain parameters 
  double Lx = Rmid_max-Rmid_min; // Domain size along x.
  double x_min = 0.;
  double x_max = Lx;
  double x_LCFS = R_LCFSmid - Rmid_min; // Radial location of the last closed flux surface.
  double q0 = qprofile(R0); // Safety factor in the center of domain.
  double Cy = r0/q0; // Normalization in binormal coordinate.

  double Lz = 2.*M_PI-1e-10; // Domain size along magnetic field.
  double z_min = -Lz/2.;
  double z_max = Lz/2.;

  // Collision frequencies
  // We note that initial profile are unstable for the full collision frequency (nu=1.0).
  // This factor can be increased once a saturated state is reached.
  double nuFrac = 0.1;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4) * n0 /
    (6*sqrt(2.) * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(me) * pow(Te0,3./2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  // Source parameters
  int num_sources = 2;
  double P_exp = 0.34e6; // P_sol measured [W]
  double recycling_rate = 0.95; // High recycling regime.
  // Core source:
  // - Injects energy only in the core region (0.25MW per species).
  // - The particles injection is only the one that are lost through the inner radial boundary.
  bool adapt_energy_srcCORE = true; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcCORE = true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcCORE = P_exp; // What the source must inject in energy [W]
  double particle_srcCORE = 0.0;// What the source must inject in particle [1/s]
  double center_srcCORE[2] = {x_min, -Lz/4}; // This is the position of the ion source,
  double sigma_srcCORE[2] = {0.03*Lx, Lz/6}; //  the electron source will be at +Lz/2.
  double floor_srcCORE = 1e-10;
  // Recycling source:
  // - Reinjects particles that are absorbed by the wall.
  // - Energy is free to leave the system.
  bool adapt_energy_srcRECY = false;
  bool adapt_particle_srcRECY = true;
  double adapt_particle_fraction_srcRECY = recycling_rate;
  double energy_srcRECY = 0.0; // [W]
  double particle_srcRECY = 0.0; // [1/s]
  double center_srcRECY[2] = {0.5*x_LCFS, M_PI};
  double sigma_srcRECY[2] = {0.25*x_LCFS, 0.05*Lz};
  double floor_srcRECY = 1e-10;

  // Grid parameters
  int Nx = 9; // The LCFS is positionned at 1/3 of the domain -> the resolution must be divisible by 3.
  int Nz = 8;
  int Nvpar = 8;
  int Nmu = 8;
  int poly_order = 1;
  // Velocity box dimensions
  double vpar_max_elc = 5.*vte;
  double mu_max_elc = me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 5.*vti;
  double mu_max_ion = mi*pow(4*vti,2)/(2*B0);
  double t_end = 1.e-7; // Should take 8 time steps
  int num_frames = 1;
  double write_phase_freq = 1.0;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .R_axis = R_axis,
    .R0 = R0,
    .a_mid = a_mid,
    .x_inner = x_inner,
    .r0 = r0,
    .B0 = B0,
    .kappa = kappa,
    .delta = delta,
    .q0 = q0,
    .Cy = Cy,
    .Lx = Lx,
    .Lz = Lz,
    .x_min = x_min, .x_max = x_max,
    .z_min = z_min, .z_max = z_max,
    .Bref = Bref,
    .x_LCFS = x_LCFS,
    .me = me, .qe = qe,
    .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0,
    .nuFrac = nuFrac, .nuElc = nuElc, .nuIon = nuIon,
    .num_sources = num_sources,

    .adapt_energy_srcCORE = adapt_energy_srcCORE,
    .adapt_particle_srcCORE = adapt_particle_srcCORE,
    .center_srcCORE = {center_srcCORE[0], center_srcCORE[1]},
    .sigma_srcCORE = {sigma_srcCORE[0], sigma_srcCORE[1]},
    .energy_srcCORE = energy_srcCORE, .particle_srcCORE = particle_srcCORE,
    .floor_srcCORE = floor_srcCORE,

    .adapt_energy_srcRECY = adapt_energy_srcRECY,
    .adapt_particle_srcRECY = adapt_particle_srcRECY,
    .adapt_particle_fraction_srcRECY = adapt_particle_fraction_srcRECY,
    .center_srcRECY = {center_srcRECY[0], center_srcRECY[1]},
    .sigma_srcRECY = {sigma_srcRECY[0], sigma_srcRECY[1]},
    .energy_srcRECY = energy_srcRECY, .particle_srcRECY = particle_srcRECY,
    .floor_srcRECY = floor_srcRECY,
    
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc, .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion, .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .t_end = t_end, .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

int 
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif
  
  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Projection parameters for the sources.
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcCORE[0], -ctx.center_srcCORE[1]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcRECY[0], ctx.center_srcRECY[1]},
    .gaussian_std_dev = {ctx.sigma_srcRECY[0], ctx.sigma_srcRECY[1]},
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcRECY,
  };
  struct gkyl_gyrokinetic_projection proj_srcRECY_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcRECY[0], ctx.center_srcRECY[1]},
    .gaussian_std_dev = {ctx.sigma_srcRECY[0], ctx.sigma_srcRECY[1]},
    .total_num_particles = ctx.particle_srcRECY,
    .total_kin_energy = ctx.energy_srcRECY,
    .temp_max = 5.0*ctx.Ti0,
    .temp_min = 0.1*ctx.Ti0,
    .f_floor = ctx.floor_srcRECY,
  };
  
  // Source adaptation parameters.
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e ={
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i ={
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1, // Only the inner radial boundary.
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_e = {
    .adapt_to_species = "ion", // Adapt to ion losses to maintain ambipolarity.
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .has_adapt_particle_fraction = true,
    .adapt_particle_fraction = ctx.adapt_particle_fraction_srcRECY,
    .num_boundaries = 3, // Outer radial boundary and both z boundaries.
    .dir = {0, 1, 1},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcRECY_i = {
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcRECY,
    .adapt_energy = ctx.adapt_energy_srcRECY,
    .has_adapt_particle_fraction = true,
    .adapt_particle_fraction = ctx.adapt_particle_fraction_srcRECY,
    .num_boundaries = 3, // Outer radial boundary and both z boundaries.
    .dir = {0, 1, 1},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };

  // Fill the species structures.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_elc,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "ion"},
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_e,
      .adapt[0] = adapt_srcCORE_e,
      .projection[1] = proj_srcRECY_e,
      .adapt[1] = adapt_srcRECY_e,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
    .time_rate_diagnostics = true,
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init,
      .upar = zero_func,
      .temp = temp_ion,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nuFrac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "elc"},
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_i,
      .adapt[0] = adapt_srcCORE_i,
      .projection[1] = proj_srcRECY_i,
      .adapt[1] = adapt_srcRECY_i,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 9,
    .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_BIMAXWELLIAN, 
      GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, 
      GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP},
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
    .time_rate_diagnostics = true,
  };

  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {
     .perp_dirs = {0, 1}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_min}, // Coordinates of the line in perpendicular directions.
     .val = 0.0, // Biasing value.
    },
    {
     .perp_dirs = {0, 1}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_max}, // Coordinates of the line in perpendicular directions.
     .val = 0.0, // Biasing value.
    },
  };
  
  struct gkyl_poisson_bias_line_list bias_line_list = {
    .num_bias_line = 2,
    .bl = target_corner_bcs,
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx.Bref,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .bias_line_list = &bias_line_list,
    .time_rate_diagnostics = true,
  };

  // Geometry
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.},
    .mapc2p = mapc2p, // Mapping of computational to physical space.
    .c2p_ctx = &ctx,
    .bfield_func = bfield_func, // magnetic field
    .bfield_ctx = &ctx,
    .has_LCFS = true,
    .x_LCFS = ctx.x_LCFS,
  };

  // Parallelism
  struct gkyl_app_parallelism_inp parallelism = {
    .comm = comm,
    .cuts = {app_args.cuts[0], app_args.cuts[1]},
    .use_gpu = app_args.use_gpu,
  };

  // GK app
  struct gkyl_gk app_inp = {
    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = { ctx.x_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = geometry,

    .num_periodic_dir = 0,

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = parallelism
  };

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp.name, sizeof(app_inp.name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
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

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
