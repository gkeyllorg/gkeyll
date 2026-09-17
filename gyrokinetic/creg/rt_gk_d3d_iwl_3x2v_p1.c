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

// Define the context of the simulation. This is basically all the globals
struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.
  
  // Geometry and magnetic field.
  double a_shift; // Parameter in Shafranov shift.
  double Z_axis; // Magnetic axis height [m].
  double R_axis; // Magnetic axis major radius [m].

  double R0; // Major radius of the simulation box [m].
  double a_mid; // Minor radius at outboard midplane [m].
  double r0; // Minor radius of the simulation box [m].
  double B0; // Magnetic field magnitude in the simulation box [T].
  double kappa; // Elongation (=1 for no elongation).
  double delta; // Triangularity (=0 for no triangularity).
  double q0; // Magnetic safety factor in the center of domain.
  double Cy; // Prefactor in binormal coordinate.

  double x_LCFS; // Radial location of the last closed flux surface.
  double x_inner; // Domain size inside the separatrix.

  double side_wall_bias; // Potential of the side wall.
  // Plasma parameters.
  double me; double qe;
  double mi; double qi;
  double n0; double Te0; double Ti0; 

  // Collisions.
  double nu_frac;

  // Source parameters.
  double n_srcOMP; // Amplitude of the OMP source
  double x_srcOMP; // Radial location of the OMP source.
  double Te_srcOMP; // Te for the OMP source.
  double Ti_srcOMP; // Ti for the OMP source.
  double sigma_srcOMP; // Radial spread of the OMP source.
  double n_srcGB; // Amplitude of the grad-B source
  double x_srcGB; // Radial location of the grad-B source.
  double sigma_srcGB; // Radial spread of the grad-B source.
  double bfac_srcGB; // Field aligned spread of the grad-B source. 
  double Te_srcGB; // Te for the grad-B source.
  double Ti_srcGB; // Ti for the grad-B source.
  double floor_src; // Source floor.

  // Grid parameters.
  double Lx; // Domain size in radial direction.
  double Ly; // Domain size in binormal direction.
  double Lz; // Domain size along magnetic field.
  double x_min; double x_max;
  double y_min; double y_max;
  double z_min; double z_max;
  int Nx;
  int Ny;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_elc; double mu_max_elc;
  double vpar_max_ion; double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

double r_x(double x, double a_mid, double x_inner)
{
  return x+a_mid-x_inner;
}

double qprofile(double R) 
{
  // Magnetic safety factor as a function of minor radius r.
  double a[] = {49.46395467479657, -260.79513158768754, 458.42618139184754, -267.63441353752336};
  return a[0]*pow(R,3) + a[1]*pow(R,2) + a[2]*R + a[3];
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
  // Z (height) as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double Z_axis = app->Z_axis;
  double kappa = app->kappa;
  return Z_axis + kappa*r*sin(theta);
}

// Partial derivatives of R(r,theta) and Z(r,theta)
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

// Common source profiles.
void density_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_srcOMP = app->n_srcOMP;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double floor_src = app->floor_src;

  fout[0] = n_srcOMP*(exp(-(pow(x-x_srcOMP,2))/(2.*pow(sigma_srcOMP,2)))+floor_src);
}
void zero_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Electron source profiles.
void density_elc_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_srcGB = app->n_srcGB;
  double x_srcGB = app->x_srcGB;
  double sigma_srcGB = app->sigma_srcGB;
  double bfac_srcGB = app->bfac_srcGB;

  fout[0] = n_srcGB*exp(-pow(x-x_srcGB,2)/(2.*pow(sigma_srcGB,2)))
           *GKYL_MAX2(sin(z)*exp(-pow(fabs(z),1.5)/(2*pow(bfac_srcGB,2))),0.);
}
void temp_elc_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Te_srcOMP = app->Te_srcOMP;

  if (x < x_srcOMP + 3*sigma_srcOMP) {
    fout[0] = Te_srcOMP;
  } else {
    fout[0] = Te_srcOMP*3./8.;
  }
}
void temp_elc_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Te_srcGB = app->Te_srcGB;

  fout[0] = Te_srcGB;
}

// Ion source profiles.
void density_ion_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n_srcGB = app->n_srcGB;
  double x_srcGB = app->x_srcGB;
  double sigma_srcGB = app->sigma_srcGB;
  double bfac_srcGB = app->bfac_srcGB;

  fout[0] = n_srcGB*exp(-pow(x-x_srcGB,2)/(2.*pow(sigma_srcGB,2)))
           *GKYL_MAX2(-sin(z)*exp(-pow(fabs(z),1.5)/(2*pow(bfac_srcGB,2))),0.);
}
void temp_ion_srcOMP(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double x_srcOMP = app->x_srcOMP;
  double sigma_srcOMP = app->sigma_srcOMP;
  double Ti_srcOMP = app->Ti_srcOMP;

  if (x < x_srcOMP + 3*sigma_srcOMP) {
    fout[0] = Ti_srcOMP;
  } else {
    fout[0] = Ti_srcOMP*3./8.;
  }
}
void temp_ion_srcGB(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Ti_srcGB = app->Ti_srcGB;

  fout[0] = Ti_srcGB;
}

// Initial density.
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;

  fout[0] = 0.5*n0*(0.5*(1.+tanh(2.*(2.-25.*x)))+0.01);
}

// Initial electron temperature.
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Te0 = app->Te0;

  fout[0] = Te0*((1./3.)*(2.+tanh(2.*(2.-25.*x)))+0.01);
}

// Initial ion temperature.
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;

  fout[0] = Ti0*((1./3.)*(2.+tanh(2.*(2.-25.*x)))+0.01);
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

void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x,a_mid,x_inner);

  fout[0] = Cy*( alpha(r, z_min, 0.0, ctx) - alpha(r, z_max, 0.0, ctx) );
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double x_inner = app->x_inner;
  double Cy = app->Cy;
  double z_min = app->z_min;
  double z_max = app->z_max;

  double r = r_x(x,a_mid,x_inner);

  fout[0] = -Cy*( alpha(r, z_min, 0.0, ctx) - alpha(r, z_max, 0.0, ctx) );
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double a_shift = 0.0; // Parameter in Shafranov shift.
  double Z_axis = 0.013055028; // Magnetic axis height [m].
  double R_axisTrue = 1.6486461; // Change R_axis to fit geometry better.
  double R_axis = 1.6; // Magnetic axis major radius [m].
  double B_axis = 2.0*R_axisTrue/R_axis; // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 2.17; // Major radius of the LCFS at the outboard midplane [m].
  double x_inner = 5*0.15/8; // Radial extent inside LCFS    
  double x_outer = 3*0.15/8; // Radial extent outside LCFS
  double Rmid_min = R_LCFSmid - x_inner; // Minimum midplane major radius of simulation box [m].
  double Rmid_max = R_LCFSmid + x_outer; // Maximum midplane major radius of simulation box [m].
  double R0 = 0.5*(Rmid_min+Rmid_max); // Major radius of the simulation box [m].

  // Minor radius at outboard midplane [m]. Redefine it with
  // Shafranov shift, to ensure LCFS radial location.
  double a_mid = fabs(a_shift)<1e-13? R_LCFSmid-R_axis :
    R_axis/a_shift - sqrt(R_axis*(R_axis - 2*a_shift*R_LCFSmid + 2*a_shift*R_axis))/a_shift;

  double side_wall_bias = 0.0; // Potential of the side wall.

  double r0 = R0-R_axis; // Minor radius of the simulation box [m].
  double B0 = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa = 1.35; // Elongation (=1 for no elongation).
  double delta = 0.4; // Triangularity (=0 for no triangularity).

  double x_LCFS = R_LCFSmid - Rmid_min; // Radial location of the last closed flux surface.

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double AMU = 2.01410177811;
  double mi = mp*AMU; // Deuterium ions.
  double Te0 = 100*eV;
  double Ti0 = 100*eV;
  double n0 = 2.0e19; // [1/m^3]

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi); // Sound speed.
  double omega_ci = fabs(qi*B0/mi); // Ion cyclotron frequency.
  double rho_s = c_s/omega_ci; // Ion sound gyroradius.

  double Lx = Rmid_max-Rmid_min; // Domain size along x.
  double Ly = 150*rho_s; // Domain size along y.
  double Lz = 2.*M_PI-1e-10; // Domain size along magnetic field.
  double x_min = 0.;
  double x_max = Lx;
  double y_min = -Ly/2.;
  double y_max = Ly/2.;
  double z_min = -Lz/2.;
  double z_max = Lz/2.;

  double q0 = qprofile(R0); // Magnetic safety factor in the center of domain.
  double Cy = r0/q0; // Normalization in binormal coordinate.

  double nu_frac = 0.1;

  // Source parameters
  double n_srcOMP = 9.e22;
  double x_srcOMP = x_min;
  double Te_srcOMP = 2*Te0;
  double Ti_srcOMP = 2*Ti0;
  double sigma_srcOMP = 0.03*Lx;
  double n_srcGB = 1.1*8.092675420182799e+21;
  double x_srcGB = x_min;
  double sigma_srcGB = 10*rho_s;
  double bfac_srcGB = 1.2;
  double Te_srcGB = 100*eV;
  double Ti_srcGB = 100*eV;
  double floor_src = 1e-2;

  // Grid parameters
  int Nx = 32;
  int Ny = 4;
  int Nz = 8;
  int Nvpar = 8;
  int Nmu = 4;
  int poly_order = 1;

  double vpar_max_elc = 4.*vte;
  double mu_max_elc = me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = mi*pow(4*vti,2)/(2*B0);

  double t_end = 1.e-7; // End time, should terminate in 43 steps.
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
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
    .side_wall_bias = side_wall_bias,
    .Cy = Cy,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .x_min = x_min, .x_max = x_max,
    .y_min = y_min, .y_max = y_max,
    .z_min = z_min, .z_max = z_max,

    .x_LCFS = x_LCFS,
  
    .me = me, .qe = qe,
    .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0,
  
    .nu_frac = nu_frac,
  
    .n_srcOMP = n_srcOMP,
    .x_srcOMP = x_srcOMP,
    .Te_srcOMP = Te_srcOMP,
    .Ti_srcOMP = Ti_srcOMP,
    .sigma_srcOMP = sigma_srcOMP,
    .n_srcGB = n_srcGB,
    .x_srcGB = x_srcGB,
    .sigma_srcGB = sigma_srcGB,
    .bfac_srcGB = bfac_srcGB,
    .Te_srcGB = Te_srcGB,
    .Ti_srcGB = Ti_srcGB,
    .floor_src = floor_src,
  
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_elc = vpar_max_elc, .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion, .mu_max_ion = mu_max_ion,

    .t_end = t_end, .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
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
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // electrons
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
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Te0,
      .nu_frac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 2,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_elc_srcOMP,
      },
      .projection[1] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_elc_srcGB,
        .upar = zero_func,
        .temp = temp_elc_srcGB,
      },
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 1,
    .diag_moments = { GKYL_F_MOMENT_MAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .num_time_rate_diagnostics = 2,
    .time_rate_diagnostics = {
      GKYL_GK_TIME_RATE_DIAGNOSTIC_FDOT_INTEGRATED_MOMENTS,
      GKYL_GK_TIME_RATE_DIAGNOSTIC_FDOT_ABS_INTEGRATED_MOMENTS,
    },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
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
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .den_ref = ctx.n0,
      .temp_ref = ctx.Ti0, 
      .nu_frac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 2,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_srcOMP,
        .upar = zero_func,
        .temp = temp_ion_srcOMP,
      },
      .projection[1] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
        .density = density_ion_srcGB,
        .upar = zero_func,
        .temp = temp_ion_srcGB,
      },
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 1,
    .diag_moments = { GKYL_F_MOMENT_MAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .num_time_rate_diagnostics = 2,
    .time_rate_diagnostics = {
      GKYL_GK_TIME_RATE_DIAGNOSTIC_FDOT_INTEGRATED_MOMENTS,
      GKYL_GK_TIME_RATE_DIAGNOSTIC_FDOT_ABS_INTEGRATED_MOMENTS,
    },

    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  struct gkyl_poisson_bias_line target_corner_bcs[] = {
    {
     .perp_dirs = {0, 2}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_min}, // Coordinates of the line in perpendicular directions.
     .val = ctx.side_wall_bias, // Biasing value.
    },
    {
     .perp_dirs = {0, 2}, // Directions perpendicular to line.
     .perp_coords = {ctx.x_LCFS, ctx.z_max}, // Coordinates of the line in perpendicular directions.
     .val = ctx.side_wall_bias, // Biasing value.
    },
  };

  struct gkyl_poisson_bias_line_list bias_line_list = {
    .num_bias_line = 2,
    .bl = target_corner_bcs,
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {ctx.side_wall_bias} },
    },
    .bias_line_list = &bias_line_list,
    .time_rate_diagnostics = true,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cfl_frac_omegaH = 1.0,
    .cfl_frac = 1.0,

    .cdim = ctx.cdim,
    .lower = { ctx.x_min, ctx.y_min, ctx.z_min },
    .upper = { ctx.x_max, ctx.y_max, ctx.z_max },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // magnetic field
      .bfield_ctx = &ctx,
      .has_LCFS = true,
      .x_LCFS = ctx.x_LCFS, // Location of last closed flux surface.
      .parallel_lower_bc_shift_func = bc_shift_func_lo,
      .parallel_upper_bc_shift_func = bc_shift_func_up,
      .parallel_lower_bc_shift_ctx = &ctx,
      .parallel_upper_bc_shift_ctx = &ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .comm = comm,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
    }
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
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
