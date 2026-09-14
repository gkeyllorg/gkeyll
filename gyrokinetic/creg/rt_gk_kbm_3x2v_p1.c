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
  double mu0;

  // Geometry and magnetic field.
  double R; // Major radius of the simulation box [m].
  double a; // Minor radius at outboard midplane [m].
  double r; // Minor radius of the simulation box [m].
  double B; // Magnetic field magnitude in the simulation box [T].

  double kperp;  // Perpendicular wavenumber.
  double beta_hat; // Plasma beta parameter.

  // Plasma parameters.
  double me;  double qe;
  double mi;  double qi;
  double n0;  double Te0;  double Ti0; 

  // Initial profile parameters.
  double r0, wTi, wTe, wN, kTi, kTe, kN;
  double L_n, L_Te, L_Ti;

  // Grid parameters.
  int Nx, Ny, Nz, Nvpar, Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double x_min, y_min, z_min;
  double x_max, y_max, z_max;
  double vpar_max_elc, mu_max_elc;
  double vpar_max_ion, mu_max_ion;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Ion initial conditions.
void density_init_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double kN = app->kN;
  double wN = app->wN;
  double a = app->a;
  double Lref = app->R;
  double r0 = app->r0;

  fout[0] = app->n0;// * exp(-kN * wN * a/Lref * tanh((x-r0)/(wN * a)));
}
void upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double kTi = app->kTi;
  double wTi = app->wTi;
  double a = app->a;
  double Lref = app->R;
  double r0 = app->r0;

  fout[0] = app->Ti0 ;//* exp(-kTi * wTi * a/Lref * tanh((x-r0)/(wTi * a)));
}

// Electron initial conditions.
void density_init_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double kN = app->kN;
  double wN = app->wN;
  double a = app->a;
  double Lref = app->R;
  double r0 = app->r0;

  fout[0] = app->n0;// * exp(-kN * wN * a/Lref * tanh((x-r0)/(wN * a)));
}
void upar_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}
void temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double kTe = app->kTe;
  double wTe = app->wTe;
  double a = app->a;
  double Lref = app->R;
  double r0 = app->r0;

  fout[0] = app->Te0;// * exp(-kTe * wTe * a/Lref * tanh((x-r0)/(wTe * a)));
}

// Slab magnetic equilibrium.
static inline void
mapc2p(double t, const double* xc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xc[0], y = xc[1], z = xc[2];
  xp[0] = x; xp[1] = y; xp[2] = z;
}
void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double x = xc[0], y = xc[1], z = xc[2];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B*app->R/(app->R + x); // this may be wrong!
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


void
write_ctx_to_json(const struct gk_app_ctx *ctx)
{
  FILE *fp = fopen("ctx.json", "w");
  if (!fp) {
    perror("Error opening file for writing context");
    return;
  }

  fprintf(fp, "{\n");
  fprintf(fp, "  \"cdim\": %d,\n", ctx->cdim);
  fprintf(fp, "  \"vdim\": %d,\n", ctx->vdim);
  fprintf(fp, "  \"mu0\": %e,\n", ctx->mu0);
  fprintf(fp, "  \"R\": %e,\n", ctx->R);
  fprintf(fp, "  \"a\": %e,\n", ctx->a);
  fprintf(fp, "  \"r\": %e,\n", ctx->r);
  fprintf(fp, "  \"B\": %e,\n", ctx->B);
  fprintf(fp, "  \"x_min\": %e,\n", ctx->x_min);
  fprintf(fp, "  \"x_max\": %e,\n", ctx->x_max);
  fprintf(fp, "  \"y_min\": %e,\n", ctx->y_min);
  fprintf(fp, "  \"y_max\": %e,\n", ctx->y_max);
  fprintf(fp, "  \"z_min\": %e,\n", ctx->z_min);
  fprintf(fp, "  \"z_max\": %e,\n", ctx->z_max);
  fprintf(fp, "  \"kperp\": %e,\n", ctx->kperp);
  fprintf(fp, "  \"beta_hat\": %e,\n", ctx->beta_hat);
  fprintf(fp, "  \"me\": %e,\n", ctx->me);
  fprintf(fp, "  \"qe\": %e,\n", ctx->qe);
  fprintf(fp, "  \"mi\": %e,\n", ctx->mi);
  fprintf(fp, "  \"qi\": %e,\n", ctx->qi);
  fprintf(fp, "  \"n0\": %e,\n", ctx->n0);
  fprintf(fp, "  \"Te0\": %e,\n", ctx->Te0);
  fprintf(fp, "  \"Ti0\": %e,\n", ctx->Ti0);
  fprintf(fp, "  \"r0\": %e,\n", ctx->r0);
  fprintf(fp, "  \"wTi\": %e,\n", ctx->wTi);
  fprintf(fp, "  \"wTe\": %e,\n", ctx->wTe);
  fprintf(fp, "  \"wN\": %e,\n", ctx->wN);
  fprintf(fp, "  \"kTi\": %e,\n", ctx->kTi);
  fprintf(fp, "  \"kTe\": %e,\n", ctx->kTe);
  fprintf(fp, "  \"kN\": %e,\n", ctx->kN);
  fprintf(fp, "  \"L_n\": %e,\n", ctx->L_n);
  fprintf(fp, "  \"L_Te\": %e,\n", ctx->L_Te);
  fprintf(fp, "  \"L_Ti\": %e,\n", ctx->L_Ti);
  fprintf(fp, "  \"Nx\": %d,\n", ctx->Nx);
  fprintf(fp, "  \"Ny\": %d,\n", ctx->Ny);
  fprintf(fp, "  \"Nz\": %d,\n", ctx->Nz);
  fprintf(fp, "  \"Nvpar\": %d,\n", ctx->Nvpar);
  fprintf(fp, "  \"Nmu\": %d,\n", ctx->Nmu);
  fprintf(fp, "  \"poly_order\": %d,\n", ctx->poly_order);
  fprintf(fp, "  \"vpar_max_elc\": %e,\n", ctx->vpar_max_elc);
  fprintf(fp, "  \"mu_max_elc\": %e,\n", ctx->mu_max_elc);
  fprintf(fp, "  \"vpar_max_ion\": %e,\n", ctx->vpar_max_ion);
  fprintf(fp, "  \"mu_max_ion\": %e,\n", ctx->mu_max_ion);
  fprintf(fp, "  \"t_end\": %e,\n", ctx->t_end);
  fprintf(fp, "  \"num_frames\": %d,\n", ctx->num_frames);
  fprintf(fp, "  \"write_phase_freq\": %e,\n", ctx->write_phase_freq);
  fprintf(fp, "  \"int_diag_calc_num\": %d,\n", ctx->int_diag_calc_num);
  fprintf(fp, "  \"dt_failure_tol\": %e,\n", ctx->dt_failure_tol);
  fprintf(fp, "  \"num_failures_max\": %d\n", ctx->num_failures_max);
  fprintf(fp, "}\n");

  fclose(fp);
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge
  double me = GKYL_ELECTRON_MASS;
  double mi = 2*GKYL_PROTON_MASS;
  double mu0 = GKYL_MU0;
  // Geometry, magnetic field and gradients.
  double R0 = 1.313; // Major radius of the simulation box [m].
  double a = 0.4701; // Minor radius at outboard midplane [m].
  double B0 = 1.91; // Magnetic field magnitude in the simulation box [T].
  double Ti0 = 1000*eV;
  double eps_n = 0.2; // aspect ratio (= normalized density gradient length).
  double eta_e = 2.0; // Normalized electron temperature gradient.
  double eta_i = 2.5; // Normalized ion temperature gradient.
  double ky_rhoi = 0.5; // Normalized perpendicular wavenumber.
  double kz_Ln = 0.1; // Normalized parallel wavenumber.
  double beta_i = 0.0045; // Ion beta.
  double tau = 1.0; // Temperature ratio Ti/Te.

  // Derived paramters.
  double r0 = 0.5*a; // Minor radius of the simulation box [m].
  double R = R0 + r0; // Major radius at center of simulation box [m].
  double B = B0 * R0/R; // Magnetic field at center of simulation box [T].
  double Te0 = Ti0/tau; // Electron temperature [J].
  double n0 = beta_i*B*B/(2*mu0*Ti0); // Density [m^-3].
  double vte = sqrt(Te0/me); // Electron thermal speed [m/s].
  double vti = sqrt(Ti0/mi); // Ion thermal speed [m/s].
  double c_s = sqrt(Te0/mi); // Sound speed [m/s].
  double omega_ci = fabs(qi*B/mi); // Ion cyclotron frequency [rad/s].
  double rho_s = c_s/omega_ci; // Ion sound gyroradius [m]
  double rho_e = vte/omega_ci; // Electron thermal gyroradius [m].
  double rho_i = vti/omega_ci; // Ion thermal gyroradius [m].
  double ky_min = ky_rhoi / rho_i; // Minimum ky.
  double dr = 2*M_PI/ky_min; // Cell size in x-direction [m].
  double L_n = R * eps_n;
  double L_Te = L_n/eta_e;
  double L_Ti = L_n/eta_i;
  double kz_min = kz_Ln / L_n;
  double L_parallel = 2*M_PI/kz_min;

  double kTi = R0/L_Ti;
  double kTe = R0/L_Te;
  double kN = R0/L_n;
  double wTi = 0.3;
  double wTe = 0.3;
  double wN = 0.3;

  // Grid parameters
  int Nx = 16;
  int Ny = 8;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 8;

  double x_min = r0;
  double x_max = r0 * (1+2*wN);
  double y_min = -dr/2;
  double y_max = dr/2;
  double z_min = -L_parallel/2;  // Min z-coordinate [m].
  double z_max = L_parallel/2;   // Max z-coordinate [m].
  double vpar_max_elc = 4.0*vte; // Max electron parallel velocity [m/s].
  double vpar_min_elc = -vpar_max_elc; // Min electron parallel velocity [m/s].
  double vpar_max_ion = 4.0*vti; // Max ion parallel velocity [m/s].
  double vpar_min_ion = -vpar_max_ion; // Min ion parallel
  double mu_max_elc = me*pow(4*vte,2)/(2*B); // Max electron mu [J/T].
  double mu_max_ion = mi*pow(4*vti,2)/(2*B); // Max ion mu [J/T].

  int poly_order = 1;

  double t_end = 0.01*L_n/vti; // End time [s].
  int num_frames = 100;
  double write_phase_freq = 0.1; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*1000;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 50; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mu0 = mu0,
    .R = R,
    .B = B,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,

    .kN = kN,  .kTi = kTi,  .kTe = kTe,
    .wN = wN,  .wTi = wTi,  .wTe = wTe,
    .r0 = r0,  .L_n = L_n,  .L_Te = L_Te,  .L_Ti = L_Ti,
    .a = a,

    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .x_min = x_min,
    .x_max = x_max,
    .y_min = y_min,
    .y_max = y_max,
    .z_min = z_min,
    .z_max = z_max,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,

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

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    // .mapc2p = {
    //   .mapping = mapc2p_vel_elc,
    //   .ctx = &ctx,
    // },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init_elc,
      .upar= upar_elc,
      .temp = temp_elc,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_EM,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    // .mapc2p = {
    //   .mapping = mapc2p_vel_ion,
    //   .ctx = &ctx,
    // },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
      .density = density_init_ion,
      .upar= upar_ion,
      .temp = temp_ion,      
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_EM,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },
    },

    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
    },
    .ampere_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
      { .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
    },
    .mu0 = ctx.mu0,
    .time_rate_diagnostics = true,
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "wk/kbm",
    .cfl_frac = 1.0,
    .cfl_frac_omegaH = 1.0,

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
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .comm = comm,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
    },
  };

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
      .enabled = true,
      .frequency = 0.1,
    }
  };

  write_ctx_to_json(&ctx);

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  return 0;
}
