#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_math.h>
#include <rt_arg_parse.h>

// Function to perform 1D linear interpolation based on a lookup table (LUT).
double interp_1x_lut(double x, double *lut_grid, double *lut_val, int N)
{
  double x_min = lut_grid[0];
  double x_max = lut_grid[N-1];
  // Clamp values if x goes outside of the LUT range.
  if (x <= x_min || x >= x_max) {
    fprintf(stderr, "Warning: x = %g out of LUT range [%g, %g]. Clamping to boundaries.\n", x, x_min, x_max);
    if (x <= x_min) return lut_val[0];
    if (x >= x_max) return lut_val[N-1];
  }
  // Assume uniform spacing of the lookup table.
  double dx = (x_max - x_min)/(N-1);
  int idx = (int)((x - x_min)/dx);
  // Safety check.
  if (idx < 0 || idx >= N-1) {
    fprintf(stderr, "Error: idx out of bounds in interp_Psi_lut: idx = %d, x = %g\n", idx, x);
    return 0.0;
  }
  // Linear interpolation.
  double x0 = lut_grid[idx];
  double x1 = lut_grid[idx+1];
  double f0 = lut_val[idx];
  double f1 = lut_val[idx+1];
  return f0 + (f1 - f0)*(x - x0)/(x1 - x0);
}

// Define the context of the simulation. This stores global parameters.
struct gk_app_ctx {
  int cdim, vdim;
  // Geometry and magnetic field parameters
  double a_shift, Z_axis, R_axis, R0, a_mid, r0, B0, kappa, delta, q0, Cy, qaxis, qlcfs;
  // Plasma parameters
  double me, qe, mi, qi, n0, Te0, Ti0;
  // Collision parameters
  double nuFrac;
  // Initial condition parameters
  double Ln, LTe, LTi;
  bool can_max;
  // Krook and buffer parameters.
  double nu_krook;
  int num_cell_buff;
  // Grid parameters
  double Lx, Lz;
  double x_min, z_min, x_max, z_max;
  int Nx, Nz, Nvpar, Nmu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;
  double max_run_time; // Maximum run time in seconds, 0 means no limit.
  // Table for intPsi interpolation and lookup table (LUT).
  int psi_lut_size;
  double *r_lut;
  double *psi_lut;
  // Table for dPsidr.
  double *dPsidr_int_lut;
};

// Geometry related functions 
double r_x(double x, double r0)
{
  return x + r0;
} 

// quadratic q profile
double qprofile(double r, double a_mid, double qaxis, double qlcfs)
{
  // Profile from Grandgirard et al. 2008
  // return 1.0 + 2.78 * pow(r/a_mid,2.8);
  // Profile from Qu et al. 2026 https://hal.science/hal-05342386v2
  // return 0.86 - 0.16 * (r/a_mid) + 2.52 * pow(r/a_mid,2);
  return 1.4;
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
  return -a_shift*r/(R_axis) + cos(theta + asin(delta)*sin(theta));
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
  double theta;
};

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}

double integrand_JoRsq(double t, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  double r = inctx->r;
  struct gk_app_ctx *app = inctx->app_ctx;
  return Jr(r,t,app) / pow(R_rtheta(r,t,app),2);
}

double intdPsidr(double r, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand_JoRsq, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);
  return integral.res;
}

double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  double integral_val = interp_1x_lut(r, app->r_lut, app->dPsidr_int_lut, app->psi_lut_size);
  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);
  return ( R*Bt/(2.*M_PI*qprofile(r, app->a_mid, app->qaxis, app->qlcfs)))*integral_val;
}

double integrant_dpsi(double r, void *int_ctx)
{
  struct integrand_ctx *inctx = int_ctx;
  struct gk_app_ctx *app = inctx->app_ctx;
  double q = qprofile(r, app->a_mid, app->qaxis, app->qlcfs);
  double dpsi = -r/q; // This is valid in the circular limit. The minus sign comes dBtheta/r < 0
  return dpsi; // I am not sure about this sign, I think it must be from dx/dr = -1.
  // return -dPsidr(r, inctx->theta, app); // Seems to be the exact way but it makes the profile super flat (and -1 factor looks important here).
}

double intPsi(double r0, double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .theta = theta};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrant_dpsi, &tmp_ctx, r0, r, 7, 1e-10);
  return integral.res;
}

double compute_alpha_integral(double r, double twrap, void *ctx) {
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;

  if (twrap == 0.0) return 0.0;
  
  if (0. < twrap) {
    integral = gkyl_dbl_exp(integrand_JoRsq, &tmp_ctx, 0., twrap, 7, 1e-10);
    return integral.res;
  } else {
    integral = gkyl_dbl_exp(integrand_JoRsq, &tmp_ctx, twrap, 0., 7, 1e-10);
    return -integral.res;
  }
}

double alpha(double r, double theta, double phi, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double twrap = theta;
  while (twrap < -M_PI) twrap = twrap+2.*M_PI;
  while (M_PI < twrap) twrap = twrap-2.*M_PI;

  double integral_val = compute_alpha_integral(r, twrap, ctx);

  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);

  return phi - R*Bt*integral_val/dPsidr(r,theta,ctx);
}

double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;
  double r0 = app->r0;
  double r = r_x(x,r0);
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);

  double drdtheta = dRdtheta(r,z,ctx);
  double dzdtheta = dZdtheta(r,z,ctx);
  double den = sqrt(pow(drdtheta,2) + pow(dzdtheta,2));
  double B_r = Bp*drdtheta/den;
  double B_z = Bp*dzdtheta/den;
  double phi = y/Cy + alpha(r, z, 0, ctx);
  double R = R_rtheta(r, z, ctx);

  // xc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = B_r * cos(phi) + Bt * sin(phi);
  fout[1] = B_r * sin(phi) - Bt * cos(phi);
  fout[2] = B_z;
}

double Bmag(double r, double theta, void *ctx)
{
  double R = R_rtheta(r,theta,ctx);
  double Bt = Bphi(R,ctx);
  double Bp = dPsidr(r,theta,ctx)/R_rtheta(r,theta,ctx)*gradr(r,theta,ctx);

  return sqrt(pow(Bt,2) + pow(Bp,2));
}

double rbar(double m, double q, double r, double theta, double vpar, double mu, void *ctx)
{
  struct gk_app_ctx *app = ctx;

  double e = GKYL_ELEMENTARY_CHARGE;
  double R = R_rtheta(r,theta,ctx);
  double r0 = app->r0;
  double R0 = app->R0;
  double q0 = app->q0;

  double psi_r = interp_1x_lut(r, app->r_lut, app->psi_lut, app->psi_lut_size);
  double psi_r0 = interp_1x_lut(r0, app->r_lut, app->psi_lut, app->psi_lut_size);

  double rpsi = q0/r0 * (psi_r - psi_r0);
  
  double rmin = r0 - app->x_min;
  double rmax = r0 + app->x_max;
  double zmin = app->z_min;
  double zmax = app->z_max;
  double Bmax = Bmag(rmax, zmax, ctx);
  
  double B0 = app->B0;
  double B = Bmag(r,theta,ctx);
  double E = 0.5 * m * pow(vpar,2) + mu*B;  
  double sgnvpar = vpar>0 ? 1.0 : -1.0;
  double vparbar = E - mu*Bmax > 0 ? 
    sgnvpar*sqrt(2/m)*sqrt(E - mu*Bmax) : 0.0;

  double rvpar = m*q0/(q*B0*r0) * (R*vpar - R0*vparbar);

  double rbar = r0 - rpsi - rvpar;

  return app->can_max ? rbar : r;
}

// Static BGK source term to maintain eq. profile at lower x boundary. 
void double_buffer_profile(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  struct gk_app_ctx *app = ctx;
  int nx = app->Nx;
  double x_max = app->x_max;
  double x_min = app->x_min;
  double Lx = app->Lx;
  double Hbuff;
  // See Eq. 49 of V. Grandgirard et al. / Computer Physics Communications 207 (2016) 35–68
  double Bs = 0.015; // Buffer transition width as a fraction of the domain size.
  double Bl = app->num_cell_buff / (double)nx; // Buffer fraction.
  Hbuff = 1 + 0.5 * (tanh((x - x_max + Bl*Lx)/(Bs*Lx)) - tanh((x - x_min - Bl*Lx)/(Bs*Lx)));
  fout[0] = Hbuff * app->nu_krook;
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  fout[0] = 0.5; // Diffusivity [m^2/s].
}

double tanh_profile(double x, double v0, double Lgrad, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  // Profile use in Gysela (see V. Grandgirard et al. / Computer Physics Communications 207 (2016) 35–68).
  double buff_frac = (double)app->num_cell_buff / app->Nx; 
  double delta = 0.5*(1-buff_frac) * app->Lx;
  double arg = x / (app->a_mid * delta);
  double prof_factor = (app->a_mid * delta) / Lgrad;
  return v0 * exp(-prof_factor * tanh(arg));
}

// Density initial condition (like TCV exp profile)
double density_init(double x, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  return tanh_profile(x, app->n0, app->Ln, ctx);
}

// Electron temperature initial conditions
double temp_init_elc(double x, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  return tanh_profile(x, app->Te0, app->LTe, ctx);
}

// Ion temperature initial conditions
double temp_init_ion(double x, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  return tanh_profile(x, app->Ti0, app->LTi, ctx);
}

double maxellian(double m, double n, double T, double E)
{
  return n/pow(2*M_PI*T/m,1.5) * exp(-E/T);
}

void eval_canon_maxwellian_e(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0], z = xn[1], vpar = xn[2], mu = xn[3];
  double m = app->me;

  double r = r_x(x, app->r0);
  double theta = z;

  double rb = app->can_max ? rbar(app->me, app->qe, r, theta, vpar, mu, ctx) : r;
  double xb = rb - app->r0;
  double dens = density_init(xb, ctx);
  double temp = temp_init_elc(xb, ctx);

  double envelope = exp(-pow(xb / (app->Lx / 4.0), 2)) * exp(-pow((theta - M_PI) / (M_PI / 4.0), 2));
  dens = dens * envelope;

  double energy = 0.5*m*pow(vpar,2) + mu*Bmag(r,theta,ctx);

  fout[0] = maxellian(m, dens, temp, energy);
}

void eval_canon_maxwellian_i(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0], z = xn[1], vpar = xn[2], mu = xn[3];
  double m = app->mi;

  double r = r_x(x, app->r0);
  double theta = z;

  double rb = app->can_max ? rbar(app->mi, app->qi, r, theta, vpar, mu, ctx) : r;
  double xb = rb - app->r0;
  double dens = density_init(xb, ctx);
  double temp = temp_init_ion(xb, ctx);

  double energy = 0.5*m*pow(vpar,2) + mu*Bmag(r,theta,ctx);

  fout[0] = maxellian(m, dens, temp, energy);
}

void eval_dens_ic(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = density_init(xn[0], ctx);
}

void eval_temp_elc_ic(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = temp_init_elc(xn[0], ctx);
}

void eval_temp_ion_ic(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = temp_init_ion(xn[0], ctx);
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];

  struct gk_app_ctx *app = ctx;
  double a_mid = app->a_mid;
  double Cy = app->Cy;

  double r = r_x(x,app->r0);

  // Map to cylindrical (R, Z, phi) coordinates.
  double R = R_rtheta(r, z, ctx);
  double Z = Z_rtheta(r, z, ctx);
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

struct gk_app_ctx create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field. (Taken from Greenfield et al. 1997 table II)
  double a_shift   = 0.0; // Parameter in Shafranov shift.
  double Z_axis    = 0.0; // Magnetic axis height [m].
  double R_axis    = 1.6714; // Magnetic axis major radius [m].
  double B_axis    = 1.54; // Magnetic field at the magnetic axis [T].
  double a_mid     = 0.604; // Minor radius at OBMP [m].
  double R0        = R_axis + 0.5 * a_mid; // Major radius of the simulation box [m].
  double r0        = 0.5*a_mid; // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double kappa     = 1.0; // cirular geometry for comparison with Dimits et al. 2000.
  double delta     = 0.0; //
  double qaxis     = 1.2; // Safety factor at r=0.
  double qlcfs     = 2.0; // Safety factor at the LCFS.

  // Plasma parameters.
  double AMU = 1.0; // Hydrogen.
  double mi  = mp*AMU;
  double Te0 = 2000*eV;
  double Ti0 = 2000*eV;
  double n0  = 4.5e19; // [1/m^3] according to Fig 5 of Greenfield et al. 1997.
  double nuFrac = 1.0; // Collision factor.

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;
  double rho_i = vti/omega_ci;
  double q0 = qprofile(r0, a_mid, qaxis, qlcfs);
  double Cy = r0/q0; // Cylindrical coordinate shift for field-alignment.

  // Configuration domain parameters 
  double Lx = 150*rho_s;   // Domain size along x.
  double Lz = 2.*M_PI-1e-10;       // Domain size along magnetic field.

  double x_min = -Lx/2;
  double x_max =  Lx/2;
  double z_min = -Lz/2;
  double z_max =  Lz/2;

  // Initial conditions and gradients
  double kTi = 6.92; // R/LTi from Dimits et al. 2000.
  double etai = 3.114; // Lni/LTi from Dimits et al. 2000.
  double LTi = R0/kTi; // Ion temperature gradient scale length [m].
  double Ln = LTi*etai; // Density gradient scale length [m].
  double LTe = LTi; // Electron temperature gradient scale length [m].

  double t_unit = R0/c_s;

  // Check the global shear r/q*dq/dr in the middle of the domain
  double rL = r0 - Lx/2;
  double rR = r0 + Lx/2;
  double qL = qprofile(rL, a_mid, qaxis, qlcfs);
  double qR = qprofile(rR, a_mid, qaxis, qlcfs);
  double dq_dr = (qR-qL)/(rR-rL);
  double s0 = r0/q0*dq_dr;

  // Expected ITG growth rate and time scale
  double gamma_itg_norm = 0.1; // Fig. 1 Dimits et al. 2000 (gamma R0/vti)
  double ky_itg_norm = 0.1; // Fig. 1 Dimits et al. 2000 (k*rho_i)
  double gamma_itg = gamma_itg_norm*vti/Ln;
  double ky_itg = ky_itg_norm/rho_i;
  double t_itg = 1.0/gamma_itg;
  double taue = Te0/Ti0;
  // Eq. 2.9 from Sugama & Watanabe, JPP 2006.
  double wgam_sw2007 = sqrt(7.0 + 4.0*taue)/2.0 * q0 * (vti/R0/q0) * sqrt(1.0 + 2.0*(23.0 + 16.0*taue + 4.0*taue*taue)/pow(q0*(7.0 + 4.0*taue), 2));
  double t_gam = 2.*M_PI/wgam_sw2007;

  double rhostar = rho_s/a_mid;
  double inv_asp_ratio = a_mid/R0;

  // Grid parameters
  int Nx = 4;
  int Nz = 4;
  int Nvpar = 4;
  int Nmu = 2;
  int poly_order = 1;

  // IC, Krook, buffer, and integration lookup tables (LUTs) parameters.
  double nu_krook = 1.0/1.0e-7;
  int num_cell_buff = 2; // Number of cells in the buffer region on each side.
  bool can_max = false; // Whether to use the canonical maxwellian formulation for the IC.
  int psi_lut_nfact = 100*(poly_order+1); // Resolution factor for the psi lookup table.

  // Velocity box dimensions
  double vpar_max_elc = 4.*vte;
  double mu_max_elc = 7*Te0/B0;
  double vpar_max_ion = 4.*vti;
  double mu_max_ion = 7*Ti0/B0;
  double t_end = 0.01*t_itg;
  int num_frames = 1;
  double write_phase_freq = 1.0;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // printf("-- Simulation parameters and info ---\n");
  // printf("vte = %1.2e, vti = %1.2e, c_s = %1.1e [m/s]\n", vte, vti, c_s);
  // printf("Lx = %1.2g, rho_s = %1.2g [m]\n", Lx, rho_s);
  // printf("x_min = %1.2g, x_max = %1.2g [m]\n", x_min, x_max);
  // printf("Ln = %1.2g, LTe = %1.2g, LTi = %1.2g [m]\n", Ln, LTe, LTi);
  // printf("Cy = %1.2g, q0 = %1.2g, qL = %1.2g, qR = %1.2g, s0 = %1.2g\n", Cy, q0, qL, qR, s0);
  // printf("ε = %1.2g, ⍴* = 1/%2.2g\n", inv_asp_ratio, 1/rhostar);
  // printf("R0/c_s = %1.2e [s]\n", t_unit);
  // printf("t_itg = %1.2e [s], t_gam = %1.2e [s]\n", t_itg, t_gam);
  // printf("t_itg c_s/R0 = %1.2e, t_gam c_s/R0 = %1.2e\n", t_itg/t_unit, t_gam/t_unit);
  // printf("Num cell for ky rhos = 0.5: %1.2g\n", (2.*M_PI/ky_itg)/(Ly/Ny));

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .a_shift = a_shift,
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .r0     = r0    ,
    .B0     = B0    ,
    .kappa  = kappa ,
    .delta  = delta ,
    .q0     = q0    ,
    .Cy     = Cy    ,
    .qaxis  = qaxis ,
    .qlcfs  = qlcfs ,
    .Ln     = Ln   ,
    .LTe    = LTe   ,
    .LTi    = LTi   ,
    .Lx     = Lx    ,
    .Lz     = Lz    ,
    .num_cell_buff = num_cell_buff,
    .nu_krook = nu_krook,
    .can_max = can_max,
    .x_min = x_min,  .x_max = x_max,
    .z_min = z_min,  .z_max = z_max,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    .nuFrac = nuFrac,
    .Nx     = Nx,
    .Nz     = Nz,
    .Nvpar  = Nvpar,
    .Nmu    = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .poly_order   = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .t_end = t_end,  .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .psi_lut_size = psi_lut_nfact*Nx,
  };
  return ctx;
}

int 
main(int argc, char **argv)
{
  struct timespec timer_global = gkyl_wall_clock();
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
  int my_rank = 0;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    gkyl_comm_get_rank(comm, &my_rank);
#endif

  // Populate lookup tables to avoid integral redundancy in the geom and IC evaluations.
  ctx.r_lut = malloc(ctx.psi_lut_size*sizeof(double));
  ctx.dPsidr_int_lut = malloc(ctx.psi_lut_size * sizeof(double));
  double r_lut_min = 0.0;
  double r_lut_max = 2.0 * ctx.a_mid;
  for (int i=0; i<ctx.psi_lut_size; i++) {
    ctx.r_lut[i] = i * (r_lut_max - r_lut_min) / (ctx.psi_lut_size-1);
    ctx.dPsidr_int_lut[i] = intdPsidr(ctx.r_lut[i], &ctx);
  }
  ctx.psi_lut = malloc(ctx.psi_lut_size*sizeof(double));
  for (int i=0; i<ctx.psi_lut_size; i++) {
    ctx.psi_lut[i] = intPsi(1e-12, ctx.r_lut[i], 0, &ctx);
  }

  // electrons
  struct gkyl_gyrokinetic_projection elc_bc = {
    .proj_id = GKYL_PROJ_FUNC,
    .func = eval_canon_maxwellian_e,
    .ctx_func = &ctx,
  };

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
      .proj_id = GKYL_PROJ_FUNC,
      .func = eval_canon_maxwellian_e,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .source_bgk = {
      .source_bgk_id = GKYL_SOURCE_BGK_STATIC,
      .rate_profile = double_buffer_profile,
      .rate_profile_ctx = &ctx,
      .feq_shape = eval_canon_maxwellian_e,
      .feq_shape_ctx = &ctx,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bc, },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
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
  struct gkyl_gyrokinetic_projection ion_bc = {
    .proj_id = GKYL_PROJ_FUNC,
    .func = eval_canon_maxwellian_i,
    .ctx_func = &ctx,    
  };

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
      .proj_id = GKYL_PROJ_FUNC,
      .func = eval_canon_maxwellian_i,
      .ctx_func = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .source_bgk = {
      .source_bgk_id = GKYL_SOURCE_BGK_STATIC,
      .rate_profile = double_buffer_profile,
      .rate_profile_ctx = &ctx,
      .feq_shape = eval_canon_maxwellian_i,
      .feq_shape_ctx = &ctx,
      .write_diagnostics = true,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bc, },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
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

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
    .time_rate_diagnostics = true,
  };

  // Geometry
  struct gkyl_gyrokinetic_geometry geometry = {
    .geometry_id = GKYL_GEOMETRY_MAPC2P,
    .world = {0.},
    .mapc2p = mapc2p, // mapping of computational to physical space
    .c2p_ctx = &ctx,
    .bfield_func = bfield_func, // magnetic field magnitude
    .bfield_ctx = &ctx,
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

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

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
    .print_verbosity = {
      .disable_timings = true,
      .enabled = true,
      .frequency = 0.01,
    }
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  free(ctx.r_lut);
  free(ctx.psi_lut);
  free(ctx.dPsidr_int_lut);

  return 0;
}
