#include <gkyl_amr_core.h>

struct amr_euler_shock_bubble_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rho_pre; // Pre-shock fluid mass density.
  double u_pre; // Pre-shock fluid velocity (x-direction).
  double p_pre; // Pre-shock fluid pressure.

  double rho_post; // Post-shock fluid mass density.
  double u_post; // Post-shock fluid velocity (x-direction).
  double p_post; // Post-shock fluid pressure.

  double rho_bub; // Bubble fluid mass density.
  double u_bub; // Bubble fluid velocity (x-direction).
  double p_bub; // Bubble fluid pressure.

  // Simulation parameters.
  int Nx; // Coarse cell count (x-direction).
  int Ny; // Coarse cell count (y-direction).
  int ref_factor1; // First refinement factor (coarse-to-intermediate).
  int ref_factor2; // Second refinement factor (intermediate-to-fine).
  double Lx; // Coarse domain size (x-direction).
  double Ly; // Coarse domain size (y-direction).
  double intermediate_Lx; // Intermediate domain size (x-direction).
  double intermediate_Ly; // Intermediate domain size (y-direction).
  double fine_Lx; // Fine domain size (x-direction).
  double fine_Ly; // Fine domain size (y-direction).
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double x_loc; // Shock location (x-direction).
  double bub_loc; // Bubble location (x-direction).
  double bub_rad; // Bubble radius.
};

struct amr_euler_shock_bubble_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rho_pre = 1.0; // Pre-shock fluid mass density.
  double u_pre = -6.0; // Pre-shock fluid velocity (x-direction).
  double p_pre = 1.0; // Pre-shock fluid pressure.

  double rho_post = 5.799; // Post-shock fluid mass density.
  double u_post = 5.75; // Post-shock fluid velocity (x-direction).
  double p_post = 167.833; // Post-shock fluid pressure.

  double rho_bub = 0.138; // Bubble fluid mass density.
  double u_bub = -6.0; // Bubble fluid velocity (x-direction).
  double p_bub = 1.0; // Bubble fluid pressure.

  // Simulation parameters.
  int Nx = 16; // Coarse cell count (x-direction).
  int Ny = 16; // Coarse cell count (y-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 1.0; // Coarse domain size (x-direction).
  double Ly = 1.0; // Coarse domain size (y-direction).
  double intermediate_Lx = 0.75; // Intermediate domain size (x-direction).
  double intermediate_Ly = 0.75; // Intermediate domain size (y-direction).
  double fine_Lx = 0.5; // Fine domain size (x-direction).
  double fine_Ly = 0.5; // Fine domain size (y-direction).
  double cfl_frac = 0.85; // CFL coefficient.

  double t_end = 0.075; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double x_loc = 0.05; // Shock location (x-direction).
  double bub_loc = 0.25; // Bubble location (x-direction).
  double bub_rad = 0.15; // Bubble radius.

  struct amr_euler_shock_bubble_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rho_pre = rho_pre,
    .u_pre = u_pre,
    .p_pre = p_pre,
    .rho_post = rho_post,
    .u_post = u_post,
    .p_post = p_post,
    .rho_bub = rho_bub,
    .u_bub = u_bub,
    .p_bub = p_bub,
    .Nx = Nx,
    .Ny = Ny,
    .ref_factor1 = ref_factor1,
    .ref_factor2 = ref_factor2,
    .Lx = Lx,
    .Ly = Ly,
    .intermediate_Lx = intermediate_Lx,
    .intermediate_Ly = intermediate_Ly,
    .fine_Lx = fine_Lx,
    .fine_Ly = fine_Ly,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .x_loc = x_loc,
    .bub_loc = bub_loc,
    .bub_rad = bub_rad,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_euler_shock_bubble_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_euler_shock_bubble_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rho_pre = app->rho_pre;
  double u_pre = app->u_pre;
  double p_pre = app->p_pre;

  double rho_post = app->rho_post;
  double u_post = app->u_post;
  double p_post = app->p_post;

  double rho_bub = app->rho_bub;
  double u_bub = app->u_bub;
  double p_bub = app->p_bub;

  double x_loc = app->x_loc;
  double bub_loc = app->bub_loc;
  double bub_rad = app->bub_rad;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  double r = sqrt((x - bub_loc) * (x - bub_loc) + y * y);

  if (x < x_loc) {
    rho = rho_post; // Fluid mass density (post-shock).
    u = u_post; // Fluid velocity (post-shock).
    p = p_post; // Fluid pressure (post-shock).
  }
  else {
    rho = rho_pre; // Fluid mass density (pre-shock).
    u = u_pre; // Fluid velocity (pre-shock).
    p = p_pre; // Fluid pressure (pre-shock).
  }

  if (r < bub_rad) {
    rho = rho_bub; // Fluid mass density (bubble).
    u = u_bub; // Fluid velocity (bubble).
    p = p_bub; // Fluid pressure (bubble).
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
  struct amr_euler_shock_bubble_ctx ctx = create_ctx(); // Context for initialization functions.

  struct euler2d_double_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = 0.0,
    .coarse_y1 = -0.5 * ctx.Ly,
    .coarse_x2 = ctx.Lx,
    .coarse_y2 = 0.5 * ctx.Ly,

    .intermediate_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.intermediate_Lx),
    .intermediate_y1 = -0.5 * ctx.intermediate_Ly,
    .intermediate_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.intermediate_Lx),
    .intermediate_y2 = 0.5 * ctx.intermediate_Ly,

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_y1 = -0.5 * ctx.fine_Ly,
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),
    .refined_y2 = 0.5 * ctx.fine_Ly,

    .eval = evalEulerInit,
    .gas_gamma = ctx.gas_gamma,

    .copy_x = true,
    .copy_y = true,

    .wall_x = false,
    .wall_y = false,

    .euler_output = "amr_euler_shock_bubble_l2",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler2d_run_double(argc, argv, &init);
}