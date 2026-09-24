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
  double R_axis;    // Magnetic axis major radius [m].
  double R0;        // Major radius of the simulation box [m].
  double a_mid;     // Minor radius at outboard midplane [m].
  double r0;        // Minor radius of the simulation box [m].
  double B0;        // Magnetic field magnitude in the simulation box [T].
  double qSep;      // Safety factor at the separatrix.
  double sSep;      // Magnetic shear at the separatrix.
  double kappa;     // Elongation (=1 for no elongation).
  double delta;     // Triangularity (=0 for no triangularity).
  double q0;        // Magnetic safety factor in the center of domain.
  double Lz;        // Domain size along magnetic field.
  double z_min;  double z_max;

  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Collisions.
  double nuFrac;  double nuIon;

  // Source parameters.
  double n_src;  double Ti_src;

  // Grid parameters.
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double vpar_max_ion;  double mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

double qprofile(double r, double a_mid, double qSep, double sSep) 
{
  // Magnetic safety factor as a function of minor radius r.
  return qSep/(1.-sSep*((r-a_mid)/a_mid));
}

double R_rtheta(double r, double theta, void *ctx)
{
  // Major radius as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double R_axis = app->R_axis;
  double delta = app->delta;
  return R_axis + r*cos(theta + asin(delta)*sin(theta));
}

double Z_rtheta(double r, double theta, void *ctx)
{
  // Z (height) as a function of minor radius r and poloidal angle theta.
  struct gk_app_ctx *app = ctx;
  double kappa = app->kappa;
  return kappa*r*sin(theta);
}

// Partial derivatives of R(r,theta) and Z(r,theta)
double dRdr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double delta = app->delta;
  return cos(theta + asin(delta)*sin(theta));
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
double dPsidr(double r, double theta, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  struct integrand_ctx tmp_ctx = {.app_ctx = app, .r = r};
  struct gkyl_qr_res integral;
  integral = gkyl_dbl_exp(integrand, &tmp_ctx, 0., 2.*M_PI, 7, 1e-10);

  double B0 = app->B0;
  double R_axis = app->R_axis;
  double a_mid = app->a_mid;
  double qSep = app->qSep;
  double sSep = app->sSep;

  return ( B0*R_axis/(2.*M_PI*qprofile(r,a_mid,qSep,sSep)) )*integral.res;
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

  double B0 = app->B0;
  double R_axis = app->R_axis;

  return phi - B0*R_axis*integral.res/dPsidr(r,theta,ctx);
}

double Bphi(double R, void *ctx)
{
  // Toroidal magnetic field.
  struct gk_app_ctx *app = ctx;
  double B0 = app->B0;
  double R0 = app->R0;
  return B0*R0/R;
}
double gradr(double r, double theta, void *ctx)
{
  return (R_rtheta(r,theta,ctx)/Jr(r,theta,ctx))*sqrt(pow(dRdtheta(r,theta,ctx),2) + pow(dZdtheta(r,theta,ctx),2));
}

// Ion source profiles.
void density_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Lz = app->Lz;
  double n_src = app->n_src;
  double z = xn[0];

  if (fabs(z) < Lz/4.) {
    fout[0] = app->n_src;
  } else {
    fout[0] = 1e-40;
  }
}
void upar_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_ion_src(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti_src = app->Ti_src;
  double z = xn[0];
  fout[0] = Ti_src;
}

// Ion initial conditions
void density_init(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Lz = app->Lz;
  double mi = app->mi;

  double z = xn[0];
  double Ls = Lz/4.;
  double xn0[] = {0.};
  double n_src0[1];
  density_src(t,xn0,n_src0,ctx);
  double effSrc  = n_src0[0];
  double eV = GKYL_ELEMENTARY_CHARGE;
  double Te_src  = 410*eV;
  double c_ss    = sqrt(5./3.*Te_src/mi);
  double nPeak   = 89.*sqrt(5.)/3./c_ss*Ls*effSrc/2.;
  double perturb = 0.;
  if (fabs(z) <= Ls) {
    fout[0] = nPeak*(1.+sqrt(1.-pow(z/Ls,2)))/2.*(1.+perturb);
  } else {
    fout[0] = nPeak/2.*(1.+perturb);
  }
}
void upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double z = xn[0];
  fout[0] = 0.0;
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti0 = app->Ti0;

  double z = xn[0];
  fout[0] = Ti0;
}

// Collision frequencies.
void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
}

// Geometry evaluation functions for the gk app
void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;
  double q0 = app->q0;
  double kappa = app->kappa;

  double x = xc[0], y = xc[1], z = xc[2];

  double r = x+r0;

  // Map to cylindrical (R, Z, phi) coordinates.
  double R   = R_rtheta(r, z, ctx);
  double Z   = kappa*r*sin(z);
  double phi = -q0/r0*y - alpha(r, z, 0, ctx);
  // Map to Cartesian (X, Y, Z) coordinates.
  double X = R*cos(phi);
  double Y = R*sin(phi);

  xp[0] = X; xp[1] = Y; xp[2] = Z;
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double r0 = app->r0;
  double q0 = app->q0;

  double x = 0., y = 0., z = xc[2];
  double r = x+r0;
  double Bt = Bphi(R_rtheta(r,z,ctx),ctx);
  double Bp = dPsidr(r,z,ctx)/R_rtheta(r,z,ctx)*gradr(r,z,ctx);

  double drdtheta = dRdtheta(r,z,ctx);
  double dzdtheta = dZdtheta(r,z,ctx);
  double den = sqrt(pow(drdtheta,2) + pow(dzdtheta,2));
  double B_r = Bp*drdtheta/den;
  double B_z = Bp*dzdtheta/den;
  double phi = -q0/r0*y - alpha(r, z, 0, ctx);
  double R   = R_rtheta(r, z, ctx);

  // xc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = B_r * cos(phi) - Bt * sin(phi);
  fout[1] = B_r * sin(phi) + Bt * cos(phi);
  fout[2] = B_z;
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS, me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Geometry and magnetic field.
  double R_axis    = 0.4;                // Magnetic axis major radius [m].
  double B_axis    = 0.224374548;        // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 0.5948;             // Major radius of the LCFS at the outboard midplane [m].
  double R0        = R_LCFSmid+0.025;    // Major radius of the simulation box [m].
  double a_mid     = R_LCFSmid-R_axis;   // Minor radius at outboard midplane [m].
  double r0        = R0-R_axis;          // Minor radius of the simulation box [m].
  double B0        = B_axis*(R_axis/R0); // Magnetic field magnitude in the simulation box [T].
  double qSep      = 3.69546081;         // Safety factor at the separatrix.
  double sSep      = 2.27976219;         // Magnetic shear at the separatrix.
  double kappa     = 1.57;               // Elongation (=1 for no elongation).
  double delta     = 0.6;                // Triangularity (=0 for no triangularity).
  double Lz        = 0.62*2.*M_PI;    // Domain size along magnetic field.

  double q0       = qprofile(0.+r0,a_mid,qSep,sSep);    // Magnetic safety factor in the center of domain.
  double epsilon0 = r0/R0;              // Inverse aspect ratio in the center of the domain.

  // Plasma parameters. Chosen based on the value of a cubic sline
  // between the last TS data inside the LCFS and the probe data in
  // in the far SOL, near R=0.475 m.
  double mi  = mp;   // Hydrogen ions.
  double Te0 = 178*eV;
  double Ti0 = 70*eV;
  double n0  = 1.78e18;   // [1/m^3]

  double nuFrac = 1.0;
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0/1e20) + 1.5 * log(Ti0/eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4) * n0 /
    (12 * pow(M_PI,3./2.) * pow(eps0,2) * sqrt(mi) * pow(Ti0,3./2.));

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi); // Thermal speeds.

  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Source parameters
  double n_src = 1.95e22;
  double Ti_src = 40*eV;

  // Grid parameters
  int Nz = 64;
  int Nvpar = 16;
  int Nmu = 45;
  int poly_order = 1;

  double vpar_max_ion = 4.*vti;
  double mu_max_ion = mi*pow(1.5*4*vti,2)/(2*B0);

  double t_end = 8.0e-6;
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .R_axis = R_axis,
    .R0     = R0    ,
    .a_mid  = a_mid ,
    .r0     = r0    ,
    .B0     = B0    ,
    .qSep   = qSep  ,
    .sSep   = sSep  ,
    .kappa  = kappa ,
    .delta  = delta ,
    .q0     = q0    ,
    .Lz     = Lz    ,
    .z_min = -Lz/2.,  .z_max = Lz/2.,
  
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
  
    .nuFrac = nuFrac,  .nuIon = nuIon,
  
    .n_src = n_src,  .Ti_src = Ti_src,
  
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,

    .t_end = t_end,  .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
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

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = density_init,
      .ctx_upar = &ctx,
      .upar= upar_ion,
      .ctx_temp = &ctx,
      .temp = temp_ion,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
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
        .ctx_density = &ctx,
        .density = density_src,
        .ctx_upar = &ctx,
        .upar= upar_ion_src,
        .ctx_temp = &ctx,
        .temp = temp_ion_src,  
      }, 
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 5,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
  };

  // GK app
  struct gkyl_gk app_inp = {

    .cdim = ctx.cdim,
    .lower = { ctx.z_min },
    .upper = { ctx.z_max },
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {0., 0.},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func, // magnetic field
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 1,
    .species = { ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
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
