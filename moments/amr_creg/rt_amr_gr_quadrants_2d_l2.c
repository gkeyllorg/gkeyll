// 2D quadrants test, using static, block-structured mesh refinement with doubly-nested refinement blocks (4x4x refinement), for the general relativistic Euler equations.
// Input parameters taken from the initial conditions in Section 4.2 (Riemann 2-D), from the article:
// L. Del Zanna and N. Bucciantini (2002), "An efficient shock-capturing central-type scheme for multdimensional flows. I. Hydrodynamics",
// Astronomy and Astrophysics, Volume 390 (3): 1177-1186.
// https://arxiv.org/abs/astro-ph/0205290

#include <gkyl_amr_core.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_alloc.h>

struct amr_gr_quadrants_2d_ctx
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

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

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

struct amr_gr_quadrants_2d_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.

  double rho_ul = 0.1; // Upper-left fluid mass density.
  double u_ul = 0.99; // Upper-left fluid x-velocity.
  double v_ul = 0.0; // Upper-left fluid y-velocity.
  double p_ul = 1.0; // Upper-left fluid pressure.

  double rho_ur = 0.1; // Upper-right fluid mass density.
  double u_ur = 0.0; // Upper-right fluid x-velocity.
  double v_ur = 0.0; // Upper-right fluid y-velocity.
  double p_ur = 0.01; // Upper-right fluid pressure.
  
  double rho_ll = 0.5; // Lower-left fluid mass density.
  double u_ll = 0.0; // Lower-left fluid x-velocity.
  double v_ll = 0.0; // Lower-left fluid y-velocity.
  double p_ll = 1.0; // Lower-left fluid pressure.

  double rho_lr = 0.1; // Lower-right fluid mass density.
  double u_lr = 0.0; // Lower-right fluid x-velocity.
  double v_lr = 0.99; // Lower-right fluid y-velocity.
  double p_lr = 1.0; // Lower-right fluid pressure.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  // Simulation parameters.
  int Nx = 8; // Coarse cell count (x-direction).
  int Ny = 8; // Coarse cell count (y-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 1.0; // Coarse domain size (x-direction).
  double Ly = 1.0; // Coarse domain size (y-direction).
  double intermediate_Lx = 0.75; // Intermediate domain size (x-direction).
  double intermediate_Ly = 0.75; // Intermediate domain size (y-direction).
  double fine_Lx = 0.5; // Fine domain size (x-direction).
  double fine_Ly = 0.5; // Fine domain size (y-direction).
  double cfl_frac = 0.9; // CFL coefficient.
  double t_end = 0.3; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double loc = 0.5; // Fluid boundaries (both x and y coordinates).

  struct amr_gr_quadrants_2d_ctx ctx = {
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
    .spacetime = spacetime,
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
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_gr_quadrants_2d_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_gr_quadrants_2d_ctx *app = &new_ctx;

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

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

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

  double spatial_det, lapse;
  double *shift = gkyl_malloc(sizeof(double[3]));
  bool in_excision_region;

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);

  double *vel = gkyl_malloc(sizeof(double[3]));
  double v_sq = 0.0;
  vel[0] = u; vel[1] = v; vel[2] = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      v_sq += spatial_metric[i][j] * vel[i] * vel[j];
    }
  }

  double W = 1.0 / (sqrt(1.0 - v_sq));
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
  }

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));
  
  // Set fluid mass density.
  fout[0] = sqrt(spatial_det) * rho * W;
  // Set fluid momentum density.
  fout[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
  fout[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
  fout[3] = 0.0;
  // Set fluid total energy density.
  fout[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

  // Set spatial metric determinant.
  fout[5] = spatial_det;
  // Set lapse gauge variable.
  fout[6] = lapse;
  // Set shift gauge variables.
  fout[7] = shift[0]; fout[8] = shift[1]; fout[9] = shift[2];

  // Set spatial metric tensor.
  fout[10] = spatial_metric[0][0]; fout[11] = spatial_metric[0][1]; fout[12] = spatial_metric[0][2];
  fout[13] = spatial_metric[1][0]; fout[14] = spatial_metric[1][1]; fout[15] = spatial_metric[1][2];
  fout[16] = spatial_metric[2][0]; fout[17] = spatial_metric[2][1]; fout[18] = spatial_metric[2][2];

  // Set inverse spatial metric tensor.
  fout[19] = inv_spatial_metric[0][0]; fout[20] = inv_spatial_metric[0][1]; fout[21] = inv_spatial_metric[0][2];
  fout[22] = inv_spatial_metric[1][0]; fout[23] = inv_spatial_metric[1][1]; fout[24] = inv_spatial_metric[1][2];
  fout[25] = inv_spatial_metric[2][0]; fout[26] = inv_spatial_metric[2][1]; fout[27] = inv_spatial_metric[2][2];

  // Set excision boundary conditions.
  if (in_excision_region) {
    for (int i = 0; i < 28; i++) {
      fout[i] = 0.0;
    }

    fout[28] = -1.0;
  }
  else {
    fout[28] = 1.0;
  }

  // Free all tensorial quantities.
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(inv_spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
  gkyl_free(inv_spatial_metric);
  gkyl_free(shift);
  gkyl_free(vel);
}

int main(int argc, char **argv)
{
  struct amr_gr_quadrants_2d_ctx ctx = create_ctx(); // Context for initialization functions.

  struct gr_euler2d_double_init init = {
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

    .eval = evalGREulerInit,
    .gas_gamma = ctx.gas_gamma,
    .spacetime = ctx.spacetime,

    .copy_x = true,
    .copy_y = true,

    .wall_x = false,
    .wall_y = false,

    .gr_euler_output = "amr_gr_quadrants_2d_l2",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  gr_euler2d_run_double(argc, argv, &init);
}