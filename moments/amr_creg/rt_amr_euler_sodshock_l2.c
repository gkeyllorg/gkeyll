// Sod-type shock tube test, using static, patch-structured mesh refinement with doubly-nested refinement patches (4x4x refinement), for the 5-moment (Euler) equations.
// Input parameters match the initial conditions in Section 2.6.2, with the contact discontinuity placed at x = 0.75 rather than x = 0.5, from the thesis:
// A. Hakim (2006), "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations",
// PhD Thesis, University of Washington.
// https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf

#include <gkyl_amr_core.h>

struct amr_euler_sodshock_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int ref_factor1; // First refinement factor (coarse-to-intermediate).
  int ref_factor2; // Second refinement factor (intermediate-to-fine).
  double Lx; // Coarse domain size (x-direction).
  double intermediate_Lx; // Intermediate domain size (x-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct amr_euler_sodshock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 3.0; // Left fluid mass density.
  double ul = 0.0; // Left fluid velocity.
  double pl = 3.0; // Left fluid pressure.

  double rhor = 1.0; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 1.0; // Right fluid pressure.

  // Simulation parameters.
  int Nx = 8; // Coarse cell count (x-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 1.0; // Coarse domain size (x-direction).
  double intermediate_Lx = 0.4; // Intermediate domain size (x-direction).
  double fine_Lx = 0.2; // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.1; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_euler_sodshock_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .Nx = Nx,
    .ref_factor1 = ref_factor1,
    .ref_factor2 = ref_factor2,
    .Lx = Lx,
    .intermediate_Lx = intermediate_Lx,
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
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_euler_sodshock_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_euler_sodshock_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < 0.75) {
    rho = rhol; // Fluid mass density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else {
    rho = rhor; // Fluid mass density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
  }
  
  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * u * u;
}

int main(int argc, char **argv)
{
  struct amr_euler_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  struct euler1d_double_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = 0.25,
    .coarse_x2 = 0.25 + ctx.Lx,

    .intermediate_x1 = (0.25 + (0.5 * ctx.Lx)) - (0.5 * ctx.intermediate_Lx),
    .intermediate_x2 = (0.25 + (0.5 * ctx.Lx)) + (0.5 * ctx.intermediate_Lx),

    .refined_x1 = (0.25 + (0.5 * ctx.Lx)) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.25 + (0.5 * ctx.Lx)) + (0.5 * ctx.fine_Lx),

    .eval = evalEulerInit,
    .gas_gamma = ctx.gas_gamma,

    .euler_output = "amr_euler_sodshock_l2",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler1d_run_double(argc, argv, &init);
}