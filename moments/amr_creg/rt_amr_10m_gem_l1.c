// Geospace Environmental Modeling (GEM) magnetic reconnection test, using static, block-structured mesh refinement with a single refinement block (2x refinement), for the 10-moment equations.
// Input parameters match the equilibrium and initial conditions in Section 2, from the article:
// J. Birn et al. (2001), "Geospace Environmental Modeling (GEM) Magnetic Reconnection Challenge",
// Journal of Geophysical Research: Space Physics, Volume 106 (A3): 3715-3719.
// https://agupubs.onlinelibrary.wiley.com/doi/10.1029/1999JA900449

#include <gkyl_amr_core.h>

struct amr_10m_gem_ctx
{
  // Mathematical constants (dimensionless).
  double pi;
  
  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double Ti_over_Te; // Ion temperature / electron temperature.
  double lambda; // Wavelength.
  double n0; // Reference number density.
  double nb_over_n0; // Background number density / reference number density.
  double B0; // Reference magnetic field strength.
  double beta; // Plasma beta.

  double k0_elc; // Electron closure parameter.
  double k0_ion; // Ion closure parameter.
  
  // Derived physical quantities (using normalized code units).
  double psi0; // Reference magnetic scalar potential.

  double Ti_frac; // Fraction of total temperature from ions.
  double Te_frac; // Fraction of total temperature from electrons.
  double T_tot; // Total temperature.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int Ny; // Coarse cell count (y-direction).
  int ref_factor; // Refinement factor.
  double Lx; // Coarse domain size (x-direction).
  double Ly; // Coarse domain size (y-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double fine_Ly; // Fine domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct amr_10m_gem_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Ion mass.
  double charge_ion = 1.0; // Ion charge.
  double mass_elc = 1.0 / 25.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double Ti_over_Te = 5.0; // Ion temperature / electron temperature.
  double lambda = 0.5; // Wavelength
  double n0 = 1.0; // Reference number density.
  double nb_over_n0 = 0.2; // Background number density / reference number density.
  double B0 = 0.1; // Reference magnetic field strength.
  double beta = 1.0; // Plasma beta.

  double k0_elc = 5.0; // Electron closure parameter.
  double k0_ion = 5.0; // Ion closure parameter.
  
  // Derived physical quantities (using normalized code units).
  double psi0 = 0.1 * B0; // Reference magnetic scalar potential.

  double Ti_frac = Ti_over_Te / (1.0 + Ti_over_Te); // Fraction of total temperature from ions.
  double Te_frac = 1.0 / (1.0 + Ti_over_Te); // Fraction of total temperature from electrons.
  double T_tot = beta * (B0 * B0) / 2.0 / n0; // Total temperature;

  // Simulation parameters.
  int Nx = 64; // Coarse cell count (x-direction).
  int Ny = 32; // Coarse cell count (y-direction).
  int ref_factor = 2; // Refinement factor.
  double Lx = 25.6; // Coarse domain size (x-direction).
  double Ly = 12.8; // Coarse domain size (y-direction).
  double fine_Lx = 16.0; // Fine domain size (x-direction).
  double fine_Ly = 6.0; // Fine domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 250.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_10m_gem_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .Ti_over_Te = Ti_over_Te,
    .lambda = lambda,
    .n0 = n0,
    .nb_over_n0 = nb_over_n0,
    .B0 = B0,
    .beta = beta,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
    .psi0 = psi0,
    .Ti_frac = Ti_frac,
    .Te_frac = Te_frac,
    .T_tot = T_tot,
    .Nx = Nx,
    .Ny = Ny,
    .ref_factor = ref_factor,
    .Lx = Lx,
    .Ly = Ly,
    .fine_Lx = fine_Lx,
    .fine_Ly = fine_Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_10m_gem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_gem_ctx *app = &new_ctx;

  double mass_elc = app->mass_elc;
  double charge_elc = app->charge_elc;
  double lambda = app->lambda;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double beta = app->beta;

  double Te_frac = app->Te_frac;
  double T_tot = app->T_tot;

  double sech_sq = (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rhoe = n * mass_elc; // Electron mass density.
  double momze = (mass_elc / charge_elc) * Jz * Te_frac; // Electron momentum density (z-direction).
  double pre = n * T_tot * Te_frac; // Electron pressure (scalar).

  // Set electron mass density.
  fout[0] = rhoe;
  // Set electron momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = momze;
  // Set electron pressure tensor.
  fout[4] = pre; fout[5] = 0.0; fout[6] = 0.0;
  fout[7] = pre; fout[8] = 0.0; fout[9] = pre + momze * momze / rhoe;
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_10m_gem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_gem_ctx *app = &new_ctx;

  double mass_ion = app->mass_ion;
  double charge_ion = app->charge_ion;
  double lambda = app->lambda;
  double n0 = app->n0;
  double nb_over_n0 = app->nb_over_n0;
  double B0 = app->B0;
  double beta = app->beta;

  double Ti_frac = app->Ti_frac;
  double T_tot = app->T_tot;

  double sech_sq = (1.0 / cosh(y / lambda)) * (1.0 / cosh(y / lambda)); // Hyperbolic secant squared.

  double n = n0 * (sech_sq + nb_over_n0); // Total number density.
  double Jz = -(B0 / lambda) * sech_sq; // Total current density (z-direction).

  double rhoi = n * mass_ion; // Ion mass density.
  double momzi = (mass_ion / charge_ion) * Jz * Ti_frac; // Ion momentum density (z-direction).
  double pri = n * T_tot * Ti_frac; // Ion pressure (scalar).

  // Set ion mass density.
  fout[0] = rhoi;
  // Set ion momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = momzi;
  // Set ion pressure tensor.
  fout[4] = pri; fout[5] = 0.0; fout[6] = 0.0;
  fout[7] = pri; fout[8] = 0.0; fout[9] = pri + momzi * momzi / rhoi;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_10m_gem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_10m_gem_ctx *app = &new_ctx;

  double pi = app->pi;

  double lambda = app->lambda;
  double B0 = app->B0;

  double psi0 = app->psi0;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double Bxb = B0 * tanh(y / lambda); // Total magnetic field strength.
  double Bx = Bxb - psi0 * (pi / Ly) * cos(2.0 * pi * x / Lx) * sin(pi * y / Ly); // Total magnetic field (x-direction).
  double By = psi0 * (2.0 * pi / Lx) * sin(2.0 * pi * x / Lx) * cos(pi * y / Ly); // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

int main(int argc, char **argv)
{
  struct amr_10m_gem_ctx ctx = create_ctx(); // Context for initialization functions.

  struct ten_moment_2d_single_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = -0.5 * ctx.Lx,
    .coarse_y1 = -0.5 * ctx.Ly,
    .coarse_x2 = 0.5 * ctx.Lx,
    .coarse_y2 = 0.5 * ctx.Ly,

    .refined_x1 = -0.5 * ctx.fine_Lx,
    .refined_y1 = -0.5 * ctx.fine_Ly,
    .refined_x2 = 0.5 * ctx.fine_Lx,
    .refined_y2 = 0.5 * ctx.fine_Ly,

    .eval_elc = evalElcInit,
    .eval_ion = evalIonInit,
    .eval_field = evalFieldInit,

    .k0_elc = ctx.k0_elc,
    .k0_ion = ctx.k0_ion,

    .light_speed = 1.0,
    .e_fact = 0.0,
    .b_fact = 1.0,

    .epsilon0 = ctx.epsilon0,
    .mass_elc = ctx.mass_elc,
    .charge_elc = ctx.charge_elc,
    .mass_ion = ctx.mass_ion,
    .charge_ion = ctx.charge_ion,

    .copy_x = true,
    .copy_y = false,

    .wall_x = false,
    .wall_y = true,

    .ten_moment_output = "amr_10m_gem_l1",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  ten_moment_2d_run_single(argc, argv, &init);
}