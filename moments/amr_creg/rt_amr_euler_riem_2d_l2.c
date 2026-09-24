// 2D Riemann (quadrant) problem, using static, block-structured mesh refinement with doubly-nested refinement blocks (4x4x refinement), for the 5-moment (Euler) equations.
// Input parameters match the initial conditions in Section 4.3, Case 3, with final time set to t = 0.8 rather than t = 0.3, from the article:
// R. Liska and B. Wendroff (2003), "Comparison of Several Difference Schemes on 1D and 2D Test Problems for the Euler Equations",
// SIAM Journal on Scientific Computing, Volume 25 (3): 995-1017.
// https://epubs.siam.org/doi/10.1137/S1064827502402120

#include <gkyl_amr_core.h>

struct amr_euler_riem_2d_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rho_ul; // Upper left fluid mass density.
  double u_ul; // Upper left fluid x-velocity.
  double v_ul; // Upper left fluid y-velocity.
  double p_ul; // Upper left fluid pressure.

  double rho_ur; // Upper right fluid mass density.
  double u_ur; // Upper right fluid x-velocity.
  double v_ur; // Upper right fluid y-velocity.
  double p_ur; // Upper left fluid pressure.
  
  double rho_ll; // Lower left fluid mass density.
  double u_ll; // Lower left fluid x-velocity.
  double v_ll; // Lower left fluid y-velocity.
  double p_ll; // Lower left fluid pressure.

  double rho_lr; // Lower right fluid mass density.
  double u_lr; // Lower right fluid x-velocity.
  double v_lr; // Lower right fluid y-velocity.
  double p_lr; // Lower right fluid pressure.

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

  double loc; // Fluid boundaries (both x and y coordinates).
};

struct amr_euler_riem_2d_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rho_ul = 0.5323; // Upper-left fluid mass density.
  double u_ul = 1.206; // Upper-left fluid x-velocity.
  double v_ul = 0.0; // Upper-left fluid y-velocity.
  double p_ul = 0.3; // Upper-left fluid pressure.

  double rho_ur = 1.5; // Upper-right fluid mass density.
  double u_ur = 0.0; // Upper-right fluid x-velocity.
  double v_ur = 0.0; // Upper-right fluid y-velocity.
  double p_ur = 1.5; // Upper-right fluid pressure.
  
  double rho_ll = 0.138; // Lower-left fluid mass density.
  double u_ll = 1.206; // Lower-left fluid x-velocity.
  double v_ll = 1.206; // Lower-left fluid y-velocity.
  double p_ll = 0.029; // Lower-left fluid pressure.

  double rho_lr = 0.5323; // Lower-right fluid mass density.
  double u_lr = 0.0; // Lower-right fluid x-velocity.
  double v_lr = 1.206; // Lower-right fluid y-velocity.
  double p_lr = 0.3; // Lower-right fluid pressure.

  // Simulation parameters.
  int Nx = 8; // Coarse cell count (x-direction).
  int Ny = 8; // Coarse cell count (y-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 1.0; // Coarse domain size (x-direction).
  double Ly = 1.0; // Coarse domain size (y-direction).
  double intermediate_Lx = 0.75; // Intermediate domain size (x-direction).
  double intermediate_Ly = 0.75; // Intermediate domain size (y-direction).
  double fine_Lx = 0.25; // Fine domain size (x-direction).
  double fine_Ly = 0.25; // Fine domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.8; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double loc = 0.8; // Fluid boundaries (both x and y coordinates).

  struct amr_euler_riem_2d_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rho_ul = rho_ul,
    .u_ul = u_ul,
    .v_ul = v_ul,
    .p_ul = p_ul,
    .rho_ur = rho_ur,
    .u_ur = u_ur,
    .v_ur = v_ur,
    .p_ur = p_ur,
    .rho_ll = rho_ll,
    .u_ll = u_ll,
    .v_ll = v_ll,
    .p_ll = p_ll,
    .rho_lr = rho_lr,
    .u_lr = u_lr,
    .v_lr = v_lr,
    .p_lr = p_lr,
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
    .loc = loc,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_euler_riem_2d_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_euler_riem_2d_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rho_ul = app->rho_ul;
  double u_ul = app->u_ul;
  double v_ul = app->v_ul;
  double p_ul = app->p_ul;

  double rho_ur = app->rho_ur;
  double u_ur = app->u_ur;
  double v_ur = app->v_ur;
  double p_ur = app->p_ur;

  double rho_ll = app->rho_ll;
  double u_ll = app->u_ll;
  double v_ll = app->v_ll;
  double p_ll = app->p_ll;

  double rho_lr = app->rho_lr;
  double u_lr = app->u_lr;
  double v_lr = app->v_lr;
  double p_lr = app->p_lr;

  double loc = app->loc;

  double rho = 0.0;
  double u = 0.0;
  double v = 0.0;
  double p = 0.0;

  if (y > loc) {
    if (x < loc) {
      rho = rho_ul; // Fluid mass density (upper-left).
      u = u_ul; // Fluid x-velocity (upper-left).
      v = v_ul; // Fluid y-velocity (upper-left).
      p = p_ul; // Fluid pressure (upper-left).
    }
    else {
      rho = rho_ur; // Fluid mass density (upper-right).
      u = u_ur; // Fluid x-velocity (upper-right).
      v = v_ur; // Fluid y-velocity (upper-right).
      p = p_ur; // Fluid pressure (upper-right).
    }
  }
  else {
    if (x < loc) {
      rho = rho_ll; // Fluid mass density (lower-left).
      u = u_ll; // Fluid x-velocity (lower-left).
      v = v_ll; // Fluid y-velocity (lower-left).
      p = p_ll; // Fluid pressure (lower-left).
    }
    else {
      rho = rho_lr; // Fluid mass density (lower-right).
      u = u_lr; // Fluid x-velocity (lower-right).
      v = v_lr; // Fluid y-velocity (lower-right).
      p = p_lr; // Fluid pressure (lower-right).
    }
  }
  
  // Set fluid mass density.
  fout[0] = rho;
  // Set fluid momentum density.
  fout[1] = rho * u; fout[2] = rho * v; fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = p / (gas_gamma - 1.0) + 0.5 * rho * (u * u + v * v);
}

int main(int argc, char **argv)
{
  struct amr_euler_riem_2d_ctx ctx = create_ctx(); // Context for initialization functions.

  struct euler2d_double_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = 0.0,
    .coarse_y1 = 0.0,
    .coarse_x2 = ctx.Lx,
    .coarse_y2 = ctx.Ly,

    .intermediate_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.intermediate_Lx),
    .intermediate_y1 = (0.5 * ctx.Ly) - (0.5 * ctx.intermediate_Ly),
    .intermediate_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.intermediate_Lx),
    .intermediate_y2 = (0.5 * ctx.Ly) + (0.5 * ctx.intermediate_Ly),

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_y1 = (0.5 * ctx.Ly) - (0.5 * ctx.fine_Ly),
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),
    .refined_y2 = (0.5 * ctx.Ly) + (0.5 * ctx.fine_Ly),

    .eval = evalEulerInit,
    .gas_gamma = ctx.gas_gamma,

    .copy_x = true,
    .copy_y = true,

    .wall_x = false,
    .wall_y = false,

    .euler_output = "amr_euler_riem_2d_l2",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler2d_run_double(argc, argv, &init);
}