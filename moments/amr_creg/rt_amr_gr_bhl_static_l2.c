// 2D Bondi-Hoyle-Lyttleton accretion problem onto a static (Schwarzschild) black hole, using static, block-structured mesh refinement with doubly-nested refinement blocks (4x4x refinement), for the general relativistic Euler equations.
// Input parameters describe wind accretion of a cold relativistic gas onto a non-rotating black hole.
// Based on the analytical solution for stiff relativistic fluids presented in the article:
// L. I. Petrich, S. L. Shapiro and S. A. Teukolsky (1988), "Accretion onto a moving black hole: An exact solution",
// Physical Review Letters, Volume 60 (18): 1781-1784.
// https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.60.1781

#include <gkyl_amr_core.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_alloc.h>

struct amr_gr_bhl_static_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.

  // Spacetime parameters (using geometric units).
  double mass; // Mass of the black hole.
  double spin; // Spin of the black hole.

  double pos_x; // Position of the black hole (x-direction).
  double pos_y; // Position of the black hole (y-direction).
  double pos_z; // Position of the black hole (z-direction).

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

  double x_loc; // Shock location (x-direction).
};

struct amr_gr_bhl_static_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 5.0 / 3.0; // Adiabatic index.

  double rhol = 3.0; // Left fluid mass density.
  double ul = 0.3; // Left fluid velocity.
  double pl = 0.05; // Left fluid pressure.

  double rhor = 0.01; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 0.01; // Right fluid pressure.

  // Spacetime parameters (using geometric units).
  double mass = 0.3; // Mass of the black hole.
  double spin = 0.0; // Spin of the black hole.

  double pos_x = 2.5; // Position of the black hole (x-direction).
  double pos_y = 2.5; // Position of the black hole (y-direction).
  double pos_z = 0.0; // Position of the black hole (z-direction).

  // Pointer to spacetime metric.
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, mass, spin, pos_x, pos_y, pos_z);

  // Simulation parameters.
  int Nx = 8; // Coarse cell count (x-direction).
  int Ny = 8; // Coarse cell count (y-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 5.0; // Coarse domain size (x-direction).
  double Ly = 5.0; // Coarse domain size (y-direction).
  double intermediate_Lx = 4.0; // Intermediate domain size (x-direction).
  double intermediate_Ly = 3.5; // Intermediate domain size (y-direction).
  double fine_Lx = 2.5; // Fine domain size (x-direction).
  double fine_Ly = 2.0; // Fine domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.
  
  double t_end = 15.0; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double x_loc = 1.0; // Shock location (x-direction).

  struct amr_gr_bhl_static_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .mass = mass,
    .spin = spin,
    .pos_x = pos_x,
    .pos_y = pos_y,
    .pos_z = pos_z,
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
    .x_loc = x_loc,
  };

  return ctx;
}

void
evalGREulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_gr_bhl_static_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_gr_bhl_static_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  struct gkyl_gr_spacetime *spacetime = app->spacetime;

  double x_loc = app->x_loc;

  double Lx = app->Lx;
  double Ly = app->Ly;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (x < x_loc) {
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

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
  
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);

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
  struct amr_gr_bhl_static_ctx ctx = create_ctx(); // Context for initialization functions.

  struct gr_euler2d_double_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = 0.0,
    .coarse_y1 = 0.0,
    .coarse_x2 = ctx.Lx,
    .coarse_y2 = ctx.Ly,

    .intermediate_x1 = (0.55 * ctx.Lx) - (0.5 * ctx.intermediate_Lx),
    .intermediate_y1 = (0.5 * ctx.Ly) - (0.5 * ctx.intermediate_Ly),
    .intermediate_x2 = (0.55 * ctx.Lx) + (0.5 * ctx.intermediate_Lx),
    .intermediate_y2 = (0.5 * ctx.Ly) + (0.5 * ctx.intermediate_Ly),

    .refined_x1 = (0.6 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_y1 = (0.5 * ctx.Ly) - (0.5 * ctx.fine_Ly),
    .refined_x2 = (0.6 * ctx.Lx) + (0.5 * ctx.fine_Lx),
    .refined_y2 = (0.5 * ctx.Ly) + (0.5 * ctx.fine_Ly),

    .eval = evalGREulerInit,
    .gas_gamma = ctx.gas_gamma,
    .spacetime = ctx.spacetime,

    .copy_x = true,
    .copy_y = true,

    .wall_x = false,
    .wall_y = false,

    .gr_euler_output = "amr_gr_bhl_static_l2",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  gr_euler2d_run_double(argc, argv, &init);
}