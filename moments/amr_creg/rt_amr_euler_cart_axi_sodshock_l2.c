// 2D Sod-type shock tube test in axial symmetry, in Cartesian coordinates, using static, block-structured mesh refinement with doubly-nested refinement blocks (4x4x refinement), for the 5-moment (Euler) equations.
// Input parameters are an axisymmetric generalization of those in Section 2.6.2, with the contact discontinuity placed at r = 0.75, from the thesis:
// A. Hakim (2006), "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations",
// PhD Thesis, University of Washington.
// https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf

#include <gkyl_amr_core.h>

struct amr_euler_cart_axi_sodshock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double gas_gamma; // Adiabatic index.

  double rhol; // Left/inner fluid mass density.
  double ul; // Left/inner fluid velocity.
  double pl; // Left/inner fluid pressure.

  double rhor; // Right/outer fluid mass density.
  double ur; // Right/outer fluid velocity.
  double pr; // Right/outer fluid pressure.

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

  double rloc; // Fluid boundary (radial coordinate).
};

struct amr_euler_cart_axi_sodshock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double gas_gamma = 1.4; // Adiabatic index.

  double rhol = 3.0; // Left/inner fluid mass density.
  double ul = 0.0; // Left/inner fluid velocity.
  double pl = 3.0; // Left/inner fluid pressure.

  double rhor = 1.0; // Right/outer fluid mass density.
  double ur = 0.0; // Right/outer fluid velocity.
  double pr = 1.0; // Right/outer fluid pressure.

  // Simulation parameters.
  int Nx = 8; // Coarse cell count (x-direction).
  int Ny = 8; // Coarse cell count (y-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 2.5; // Coarse domain size (x-direction).
  double Ly = 2.5; // Coarse domain size (y-direction).
  double intermediate_Lx = 1.8; // Intermediate domain size (x-direction).
  double intermediate_Ly = 1.8; // Intermediate domain size (y-direction).
  double fine_Lx = 1.0; // Fine domain size (x-direction).
  double fine_Ly = 1.0; // Fine domain size (y-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.2; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  double rloc = 0.5 * (0.25 + 1.25); // Fluid boundary (radial coordinate).

  struct amr_euler_cart_axi_sodshock_ctx ctx = {
    .pi = pi,
    .gas_gamma = gas_gamma,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
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
    .rloc = rloc,
  };

  return ctx;
}

void
evalEulerInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0], y = xn[1];
  struct amr_euler_cart_axi_sodshock_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_euler_cart_axi_sodshock_ctx *app = &new_ctx;

  double gas_gamma = app->gas_gamma;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;

  double rloc = app->rloc;

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  double r = sqrt(x * x + y * y);

  if (r < rloc) {
    rho = rhol; // Fluid mass density (left/inner).
    u = ul; // Fluid velocity (left/inner).
    p = pl; // Fluid pressure (left/inner).
  }
  else {
    rho = rhor; // Fluid mass density (right/outer).
    u = ur; // Fluid velocity (right/outer).
    p = pr; // Fluid pressure (right/outer).
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
  struct amr_euler_cart_axi_sodshock_ctx ctx = create_ctx(); // Context for initialization functions.

  struct euler2d_double_init init = {
    .base_Nx = ctx.Nx,
    .base_Ny = ctx.Ny,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = -0.5 * ctx.Lx,
    .coarse_y1 = -0.5 * ctx.Ly,
    .coarse_x2 = 0.5 * ctx.Lx,
    .coarse_y2 = 0.5 * ctx.Ly,

    .intermediate_x1 = -0.5 * ctx.intermediate_Lx,
    .intermediate_y1 = -0.5 * ctx.intermediate_Ly,
    .intermediate_x2 = 0.5 * ctx.intermediate_Lx,
    .intermediate_y2 = 0.5 * ctx.intermediate_Ly,

    .refined_x1 = -0.5 * ctx.fine_Lx,
    .refined_y1 = -0.5 * ctx.fine_Ly,
    .refined_x2 = 0.5 * ctx.fine_Lx,
    .refined_y2 = 0.5 * ctx.fine_Ly,

    .eval = evalEulerInit,
    .gas_gamma = ctx.gas_gamma,

    .copy_x = true,
    .copy_y = true,

    .wall_x = false,
    .wall_y = false,

    .euler_output = "amr_euler_cart_axi_sodshock_l2",

    .low_order_flux = false,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler2d_run_double(argc, argv, &init);
}