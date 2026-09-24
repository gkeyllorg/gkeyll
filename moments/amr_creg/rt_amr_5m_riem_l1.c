// Generalized Brio-Wu Riemann problem, using static, patch-structured mesh refinement with a single refinement patch (4x refinement), for the 5-moment equations.
// Input parameters match the initial conditions found in entry JE4 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je4/je4-twofluid-shock.html), adapted from Section 7.1 of the article:
// A. Hakim, J. Loverich and U. Shumlak (2006), "A high resolution wave propagation scheme for ideal Two-Fluid plasma equations",
// Journal of Computational Physics, Volume 219 (1): 418-442.
// https://www.sciencedirect.com/science/article/pii/S0021999106001707

#include <gkyl_amr_core.h>

struct amr_5m_riem_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double rhol_ion; // Left ion mass density.
  double rhor_ion; // Right ion mass density;
  double pl; // Left electron/ion pressure.
  double pr; // Right electron/ion pressure.

  double Bx; // Total magnetic field (x-direction).
  double Bzl; // Left total magneic field (z-direction).
  double Bzr; // Right total magnetic field (z-direction).

  bool has_collision; // Whether to include collisions.
  double nu_base_ei; // Base electron-ion collision frequency.

  double k0_elc; // Electron closure parameter.
  double k0_ion; // Ion closure parameter.

  // Derived physical quantities (using normalized code units).
  double rhol_elc; // Left electron mass density.
  double rhor_elc; // Right electron mass density.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int ref_factor; // Refinement factor.
  double Lx; // Coarse domain size (x-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct amr_5m_riem_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_ion = 1.0; // Proton mass.
  double charge_ion = 1.0; // Proton charge.
  double mass_elc = 1.0 / 1836.2; // Electron mass.
  double charge_elc = -1.0; // Electron charge.

  double rhol_ion = 1.0; // Left ion mass density.
  double rhor_ion = 0.125; // Right ion mass density;
  double pl = 5.0e-5; // Left electron/ion pressure.
  double pr = 5.0e-6; // Right electron/ion pressure.

  double Bx = 0.75e-2; // Total magnetic field (x-direction).
  double Bzl = 1.0e-2; // Left total magneic field (z-direction).
  double Bzr = -1.0e-2; // Right total magnetic field (z-direction).

  bool has_collision = false; // Whether to include collisions.
  double nu_base_ei = 0.5; // Base electron-ion collision frequency.

  double k0_elc = 0.0; // Electron closure parameter.
  double k0_ion = 0.0; // Ion closure parameter.

  // Derived physical quantities (using normalized code units).
  double rhol_elc = rhol_ion * mass_elc / mass_ion; // Left electron mass density.
  double rhor_elc = rhor_ion * mass_elc / mass_ion; // Right electron mass density.

  // Simulation parameters.
  int Nx = 32; // Coarse cell count (x-direction).
  int ref_factor = 4; // Refinement factor.
  double Lx = 1.0; // Coarse domain size (x-direction).
  double fine_Lx = 0.2; // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 10.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_5m_riem_ctx ctx = {
    .gas_gamma = gas_gamma,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .rhol_ion = rhol_ion,
    .rhor_ion = rhor_ion,
    .pl = pl,
    .pr = pr,
    .Bx = Bx,
    .Bzl = Bzl,
    .Bzr = Bzr,
    .has_collision = has_collision,
    .nu_base_ei = nu_base_ei,
    .k0_elc = k0_elc,
    .k0_ion = k0_ion,
    .rhol_elc = rhol_elc,
    .rhor_elc = rhor_elc,
    .Nx = Nx,
    .ref_factor = ref_factor,
    .Lx = Lx,
    .fine_Lx = fine_Lx,
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
  double x = xn[0];
  struct amr_5m_riem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_5m_riem_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double pl = app->pl;
  double pr = app->pr;

  double rhol_elc = app->rhol_elc;
  double rhor_elc = app->rhor_elc;

  double rho = 0.0;
  double p = 0.0;

  if (x < 0.5) {
    rho = rhol_elc; // Electron mass density (left).
    p = pl; // Electron pressure (left).
  }
  else {
    rho = rhor_elc; // Electron mass density (right).
    p = pr; // Electron pressure (right).
  }

  // Set electron mass density.
  fout[0] = rho;
  // Set electron momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set electron total energy density.
  fout[4] = p / (gas_gamma - 1.0);  
}

void
evalIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_5m_riem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_5m_riem_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rhol_ion = app->rhol_ion;
  double rhor_ion = app->rhor_ion;
  double pl = app->pl;
  double pr = app->pr;

  double rho = 0.0;
  double p = 0.0;

  if (x < 0.5) {
    rho = rhol_ion; // Ion mass density (left).
    p = pl; // Ion pressure (left).
  }
  else {
    rho = rhor_ion; // Ion mass density (right).
    p = pr; // Ion pressure (right).
  }

  // Set ion mass density.
  fout[0] = rho;
  // Set ion momentum density.
  fout[1] = 0.0; fout[2] = 0.0; fout[3] = 0.0;
  // Set ion total energy density.
  fout[4] = p / (gas_gamma - 1.0);  
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_5m_riem_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_5m_riem_ctx *app = &new_ctx;

  double Bx = app->Bx;
  double Bzl = app->Bzl;
  double Bzr = app->Bzr;

  double Bz = 0.0;

  if (x < 0.5) {
    Bz = Bzl; // Total magnetic field (z-direction, left).
  }
  else {
    Bz = Bzr; // Total magnetic field (z-direction, right).
  }

  // Set electric field.
  fout[0] = 0.0, fout[1] = 0.0; fout[2] = 0.0;
  // Set magnetic field.
  fout[3] = Bx, fout[4] = 0.0; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

int main(int argc, char **argv)
{
  struct amr_5m_riem_ctx ctx = create_ctx(); // Context for initialization functions.

  struct five_moment_1d_single_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = 0.0,
    .coarse_x2 = ctx.Lx,

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),

    .eval_elc = evalElcInit,
    .eval_ion = evalIonInit,
    .eval_field = evalFieldInit,

    .gas_gamma = ctx.gas_gamma,
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

    .five_moment_output = "amr_5m_riem_l1",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  five_moment_1d_run_single(argc, argv, &init);
}