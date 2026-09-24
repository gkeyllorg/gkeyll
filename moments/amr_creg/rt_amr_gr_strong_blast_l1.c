// Strongly relativistic blast wave test, using static, patch-structured mesh refinement with a single refinement patch (4x refinement), for the general relativistic Euler equations.
// Input parameters taken from the initial conditions in Section 4.1 (blast wave 2), from the article:
// L. Del Zanna and N. Bucciantini (2002), "An efficient shock-capturing central-type scheme for multidimensional relativistic flows. I. Hydrodynamics",
// Astronomy and Astrophysics, Volume 390 (3): 1177-1186.
// https://arxiv.org/abs/astro-ph/0205290

#include <gkyl_amr_core.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_alloc.h>

struct amr_gr_strong_blast_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime;

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

struct amr_gr_strong_blast_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.

  double rhol = 1.0; // Left fluid mass density.
  double ul = 0.0; // Left fluid velocity.
  double pl = 1000.0; // Left fluid pressure.

  double rhor = 1.0; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 0.01; // Right fluid pressure.

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  // Simulation parameters.
  int Nx = 128; // Coarse cell count (x-direction).
  int ref_factor = 4; // Refinement factor.
  double Lx = 1.0; // Coarse domain size (x-direction).
  double fine_Lx = 0.4; // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.35; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_gr_strong_blast_ctx ctx = {
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .spacetime = spacetime,
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
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_gr_strong_blast_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_gr_strong_blast_ctx *app = &new_ctx;

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

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  if (x < 0.5) {
    rho = rhol; // Fluid mass density (left).
    u = ul; // Fluid velocity (left).
    p = pl; // Fluid pressure (left).
  }
  else {
    rho = rhor; // Fluid mass density (right).
    u = ur; // Fluid velocity (right).
    p = pr; // Fluid pressure (right).
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

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, 0.0, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, 0.0, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, 0.0, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, 0.0, 0.0, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, 0.0, 0.0, &inv_spatial_metric);

  double *vel = gkyl_malloc(sizeof(double[3]));
  double v_sq = 0.0;
  vel[0] = u; vel[1] = 0.0; vel[2] = 0.0;

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
  fout[2] = 0.0;
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
  struct amr_gr_strong_blast_ctx ctx = create_ctx(); // Context for initialization functions.

  struct gr_euler1d_single_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor = ctx.ref_factor,

    .coarse_x1 = 0.0,
    .coarse_x2 = ctx.Lx,

    .refined_x1 = (0.75 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.75 * ctx.Lx) + (0.5 * ctx.fine_Lx),

    .eval = evalGREulerInit,
    .gas_gamma = ctx.gas_gamma,
    .spacetime = ctx.spacetime,

    .gr_euler_output = "amr_gr_strong_blast_l1",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  gr_euler1d_run_single(argc, argv, &init);
}