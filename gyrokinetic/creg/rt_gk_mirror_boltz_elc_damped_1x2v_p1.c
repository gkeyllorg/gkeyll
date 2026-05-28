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
struct gk_mirror_ctx
{
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
  // Parameters controlling initial conditions.
  double alim;
  double alphaIC0;
  double alphaIC1;
  double nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  // Thermal speeds.
  double vti;
  double vte;
  double c_s;
  // Gyrofrequencies and gyroradii.
  double omega_ci;
  double rho_s;
  double RatZeq0; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double Z_min;
  double Z_max;
  double z_min;
  double z_max;
  double psi_eval;
  double psi_in;
  double z_in;
  // Magnetic equilibrium model.
  double mcB;
  double gamma;
  double Z_m;
  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double B_bt;
  double R_bt;
  double Z_bt;
  double z_bt;
  double R_m;
  double B_m;
  double z_m;
  // Physics parameters at mirror throat
  double n_m;
  double Ti_m;
  double cs_m;
  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;
  double alpha; // Multirate factor.
  // Grid parameters
  double vpar_max_ion;
  double mu_max_ion;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
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

// -- Source functions.

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double z = xn[0];
  double Z = Z_psiz(psi, z, ctx); // Cylindrical axial coordinate.
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;
  if (fabs(Z) <= app->Z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2))));
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
  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double z = xn[0];
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
  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double z = xn[0];
  double Z = Z_psiz(psi, z, ctx); // Cylindrical axial coordinate.
  double R = R_psiZ(psi, Z, ctx); // Cylindrical radial coordinate.
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);
  if (fabs(Z) <= app->Z_bt)
  {
    fout[0] = app->n0 * pow(1.0 - pow((R - app->R_bt) / app->alim, 2), app->alphaIC0 / 2);
  }
  else if (fabs(Z) <= app->Z_m)
  {
    fout[0] = app->n0 * pow(1.0 - pow((R - app->R_bt) / app->alim, 2), app->alphaIC1 / 2);
  }
  else
  {
    fout[0] = app->n_m * sqrt(Bmag / app->B_m);
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double z = xn[0];
  if (fabs(z) <= app->z_m)
  {
    fout[0] = 0.0;
  }
  else if (z > app->z_m)
  {
    fout[0] = app->cs_m * (z - app->z_m); //* (z -  / app->z_m);
  }
  else
  {
    fout[0] = app->cs_m * (z + app->z_m); //* (z + app->z_m) / app->z_m;
  }
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double z = xn[0];
  double Z = Z_psiz(psi, z, ctx); // Cylindrical axial coordinate.
  double R = R_psiZ(psi, Z, ctx); // Cylindrical radial coordinate.
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);
  if (fabs(Z) <= app->Z_bt)
  {
    fout[0] = app->Ti0 * pow((1.0 - pow((R - app->R_bt) / app->alim, 2)), app->alphaIC0 / 2);
  }
  else if (fabs(Z) <= app->Z_m)
  {
    fout[0] = app->Ti0 * pow((1.0 - pow((R - app->R_bt) / app->alim, 2)), app->alphaIC1 / 2);
  }
  else
  {
    fout[0] = app->Ti_m * sqrt(Bmag / app->B_m);
  }
}

// Evaluate collision frequencies
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
  struct gk_mirror_ctx *app = ctx;
  double z = xc[2];
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
  struct gk_mirror_ctx *app = ctx;
  double z = xc[2];
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

void
loss_cone_damping_rate_profile(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0], vpar = xn[1], mu = xn[2];

  struct gk_mirror_ctx *app = ctx;

  double B_m = app->B_m;
  double charge = app->qi;
  double mass = app->mi;

  double bmag[1];
  double xc[] = {0.0, 0.0, z};
  bmag_func(t, xc, bmag, ctx);

  double Rm = B_m/bmag[0]; // Mirror ratio.

  double Delta_phi = 0.0; // Potential difference.

  double mu_bound = (0.5*mass*pow(vpar,2)+charge*Delta_phi)/(bmag[0]*(Rm-1));
  double loss_cone_mask;
  if (mu <= mu_bound)
    loss_cone_mask = 1.0;
  else
    loss_cone_mask = 0.0;

  double vt = app->vti;
  double z_m = app->z_m;
  double z_max = app->z_max;
  double Lpar;
  if (fabs(z) <= z_m)
    Lpar = z_m;
  else
    Lpar = z_max - z_m;

  double alpha = app->alpha;

  fout[0] = (1.0-alpha) * (vt/Lpar) * loss_cone_mask;
}

void
loss_cone_damping_rate_scaling(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0], vpar = xn[1], mu = xn[2];

  struct gk_mirror_ctx *app = ctx;
  double vt = app->vti;
  double z_m = app->z_m;
  double alpha = app->alpha;

  if (fabs(z) < z_m) {
    double Lpar = z_m;
    fout[0] = (1.0-alpha) * vt/Lpar;
  }
  else
    fout[0] = 0.0;
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
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
  double Z_max = 2.5;
  double z_min = -2.515312;
  double z_max = 2.515312;
  double psi_eval = 0.0026530898059565;

  // Parameters controlling the magnetic equilibrium model.
  double mcB = 6.51292;
  double gamma = 0.124904;
  double Z_m = 0.98;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = Z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;

  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double B_bt = 1.058278;
  double R_bt = 0.071022;
  double Z_bt = 0.467101;
  double z_bt = 0.468243;
  double R_m = 0.017845;
  double B_m = 16.662396;
  double z_m = 0.982544;

  // Physics parameters at mirror throat
  double n_m = 1.105617e19;
  double Te_m = 346.426583 * eV;
  double Ti_m = 3081.437703 * eV;
  double cs_m = 4.037740e5;

  double alpha = 0.01; // Multirate factor.

  // Grid parameters
  double vpar_max_ion = 20 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nz = 32;
  int Nvpar = 32; // Number of cells in the paralell velocity direction 96
  int Nmu = 48;  // Number of cells in the mu direction 192
  int poly_order = 1;

  double t_end = 6.0e-7;
  int num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
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
    .alim = alim,
    .alphaIC0 = alphaIC0,
    .alphaIC1 = alphaIC1,
    .nuFrac = nuFrac,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .RatZeq0 = RatZeq0,
    .Z_min = Z_min,
    .Z_max = Z_max,
    .z_min = z_min,
    .z_max = z_max,
    .psi_eval = psi_eval,
    .mcB = mcB,
    .gamma = gamma,
    .Z_m = Z_m,
    .B_bt = B_bt,
    .R_bt = R_bt,
    .Z_bt = Z_bt,
    .z_bt = z_bt,
    .R_m = R_m,
    .B_m = B_m,
    .z_m = z_m,
    .n_m = n_m,
    .Ti_m = Ti_m,
    .cs_m = cs_m,
    .NSrcIon = NSrcIon,
    .lineLengthSrcIon = lineLengthSrcIon,
    .sigSrcIon = sigSrcIon,
    .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,
    .TSrcFloorIon = TSrcFloorIon,
    .alpha = alpha,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_app_write_field_energy(app);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_write_dt(app);
    }
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
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
    .charge = ctx.qi,
    .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = {-ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = { cells_v[0], cells_v[1] },

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar_ion,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },

    .collisionless = {
//    .type = GKYL_GK_COLLISIONLESS_NONE,
      .type = GKYL_GK_COLLISIONLESS_ES,
      .scale_factor = ctx.alpha,
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
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },

    .damping = {
      .type = GKYL_GK_DAMPING_USER_INPUT,
      .rate_profile = loss_cone_damping_rate_scaling,
      .rate_profile_ctx = &ctx,
      .write_rate = true,
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
    .polarization_bmag = ctx.B_p,
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
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end           = ctx.t_end,
      .num_frames      = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol  = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart      = app_args.is_restart,
      .restart_frame   = app_args.restart_frame,
      .num_steps       = app_args.num_steps,
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
