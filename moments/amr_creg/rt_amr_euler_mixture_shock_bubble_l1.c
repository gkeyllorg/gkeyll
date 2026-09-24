#include <gkyl_amr_core.h>
#include <gkyl_alloc.h>

struct amr_shock_bubble_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma1; // First species adiabatic index.
  double gas_gamma2; // Second species adiabatic index.

  double rho_pre; // Pre-shock fluid mass density.
  double u_pre; // Pre-shock fluid velocity (x-direction).
  double alpha1_pre; // Pre-shock volume fraction (first species).

  double rho_post; // Post-shock fluid mass density.
  double u_post; // Post-shock fluid velocity (x-direction).
  double alpha1_post; // Post-shock volume fraction (first species).

  double rho_bub; // Bubble fluid mass density.
  double u_bub; // Bubble fluid velocity (x-direction).
  double alpha1_bub; // Bubble volume fraction (first species).

  // Derived physical quantities (using normalized code units).
  double p_pre; // Pre-shock fluid pressure.
  double p_post; // Post-shock fluid pressure.
  double p_bub; // Bubble fluid pressure.

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

  double x_loc; // Shock location (x-direction).
  double bub_loc_x; // Bubble location (x-direction).
  double bub_loc_y; // Bubble location (y-direction).
  double bub_rad; // Bubble radius.
};

struct amr_shock_bubble_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma1 = 1.4; // First species adiabatic index.
  double gas_gamma2 = 1.648; // Second species adiabatic index.

  double rho_pre = 1.0; // Pre-shock fluid mass density.
  double u_pre = 0.0; // Pre-shock fluid velocity (x-direction).
  double alpha1_pre = 0.99999; // Pre-shock volume fraction (first species).

  double rho_post = 1.3764; // Post-shock fluid mass density.
  double u_post = -0.3336; // Post-shock fluid velocity (x-direction).
  double alpha1_post = 0.99999; // Post-shock volume fraction (first species).

  double rho_bub = 0.1818; // Bubble fluid mass density.
  double u_bub = 0.0; // Bubble fluid velocity (x-direction).
  double alpha1_bub = 0.00001; // Bubble volume fraction (first species).

  // Derived physical quantities (using normalized code units).
  double p_pre = 1.0 / gas_gamma1; // Pre-shock fluid pressure.
  double p_post = 1.5698 / gas_gamma1; // Post-shock fluid pressure.
  double p_bub = 1.0 / gas_gamma1; // Bubble fluid pressure.

  // Simulation parameters.
  int Nx = 32; // Coarse cell count (x-direction).
  int Ny = 8; // Coarse cell count (y-direction).
  int ref_factor = 4; // Refinement factor.
  double Lx = 0.325; // Coarse domain size (x-direction).
  double Ly = 0.089; // Coarse domain size (y-direction).
  double fine_Lx = 0.24; // Fine domain size (x-direction).
  double fine_Ly = 0.06; // Fine domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.3; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double x_loc = 0.225; // Shock location (x-direction).
  double bub_loc_x = 0.175; // Bubble location (x-direction).
  double bub_loc_y = 0.5 * Ly; // Bubble location (y-direction).
  double bub_rad = 0.025; // Bubble radius.

  struct amr_shock_bubble_ctx ctx = {
    .gas_gamma1 = gas_gamma1,
    .gas_gamma2 = gas_gamma2,
    .rho_pre = rho_pre,
    .u_pre = u_pre,
    .alpha1_pre = alpha1_pre,
    .rho_post = rho_post,
    .u_post = u_post,
    .alpha1_post = alpha1_post,
    .rho_bub = rho_bub,
    .u_bub = u_bub,
    .alpha1_bub = alpha1_bub,
    .p_pre = p_pre,
    .p_post = p_post,
    .p_bub = p_bub,
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
    .x_loc = x_loc,
    .bub_loc_x = bub_loc_x,
    .bub_loc_y = bub_loc_y,
    .bub_rad = bub_rad,
  };

  return ctx;
}

void
evalEulerMixtureInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_shock_bubble_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_shock_bubble_ctx *app = &new_ctx;

  double gas_gamma1 = app->gas_gamma1;
  double gas_gamma2 = app->gas_gamma2;

  double rho_pre = app->rho_pre;
  double u_pre = app->u_pre;
  double alpha1_pre = app->alpha1_pre;

  double rho_post = app->rho_post;
  double u_post = app->u_post;
  double alpha1_post = app->alpha1_post;

  double rho_bub = app->rho_bub;
  double u_bub = app->u_bub;
  double alpha1_bub = app->alpha1_bub;

  double p_pre = app->p_pre;
  double p_post = app->p_post;
  double p_bub = app->p_bub;

  double x_loc = app->x_loc;
  double bub_loc_x = app->bub_loc_x;
  double bub_loc_y = app->bub_loc_y;
  double bub_rad = app->bub_rad;

  double rho1 = 0.0;
  double rho2 = 0.0;
  double alpha1 = 0.0;

  double vx_total = 0.0;
  double vy_total = 0.0;
  double vz_total = 0.0;
  double p_total = 0.0;

  double r = sqrt(((x - bub_loc_x) * (x - bub_loc_x)) + ((y - bub_loc_y) * (y - bub_loc_y)));

  if (x > x_loc) {
    rho1 = rho_post; // First species fluid mass density (post-shock).
    rho2 = rho_bub; // Second species fluid mass density (bubble).
    alpha1 = alpha1_post; // First species volume fraction (post-shock).

    vx_total = u_post; // Total mixture velocity (post-shock).
    p_total = p_post; // Total mixture pressure (post-shock).
  }
  else {
    rho1 = rho_pre; // First species fluid mass density (pre-shock).
    rho2 = rho_bub; // Second species fluid mass density (bubble).
    alpha1 = alpha1_pre; // First species volume fraction (pre-shock).

    vx_total = u_pre; // Total mixture velocity (pre-shock).
    p_total = p_pre; // Total mixture pressure (pre-shock).
  }

  if (r < bub_rad) {
    rho1 = rho_pre; // First species fluid mass density (pre-shock).
    rho2 = rho_bub; // Second species fluid mass density (bubble).
    alpha1 = alpha1_bub; // First species volume fraction (bubble).

    vx_total = u_bub; // Total mixture velocity (bubble).
    p_total = p_bub; // Total mixture pressure (bubble).
  }

  double rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2); // Total mixture density.

  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * (vx_total * vx_total)); // First species total energy.
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * (vx_total * vx_total)); // Second species total energy.
  double E_total = (alpha1 * E1) + ((1.0 - alpha1) * E2); // Total mixture energy.

  // Set fluid mixture total mass density.
  fout[0] = rho_total;
  // Set fluid mixture total momentum density.
  fout[1] = rho_total * vx_total; fout[2] = rho_total * vy_total; fout[3] = rho_total * vz_total;
  // Set fluid mixture total energy density.
  fout[4] = E_total;
  // Set fluid mixture weighted volume fraction (first species).
  fout[5] = rho_total * alpha1;
  // Set fluid mixture volume-weighted mass densities (first and second species).
  fout[6] = alpha1 * rho1; fout[7] = (1.0 - alpha1) * rho2;
}

int main(int argc, char **argv)
{
  struct amr_shock_bubble_ctx ctx = create_ctx(); // Context for initialization functions.

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = ctx.gas_gamma1;
  gas_gamma_s[1] = ctx.gas_gamma2;

  struct euler_mixture2d_single_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = 0.0,
    .coarse_y1 = 0.0,
    .coarse_x2 = ctx.Lx,
    .coarse_y2 = ctx.Ly,

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_y1 = (0.5 * ctx.Ly) - (0.5 * ctx.fine_Ly),
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),
    .refined_y2 = (0.5 * ctx.Ly) + (0.5 * ctx.fine_Ly),

    .eval = evalEulerMixtureInit,
    .num_species = 2,
    .gas_gamma_s = gas_gamma_s,

    .copy_x = true,
    .copy_y = false,

    .wall_x = false,
    .wall_y = true,

    .euler_mixture_output = "amr_euler_mixture_shock_bubble_l1",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler_mixture2d_run_single(argc, argv, &init);

  gkyl_free(gas_gamma_s);
}