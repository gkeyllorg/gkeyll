#include <gkyl_amr_core.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_alloc.h>

struct amr_fedkiw_shock_ctx
{
  // Physical constants (using normalized code units).
  double gas_gamma1; // First species adiabatic index.
  double gas_gamma2; // Second species adiabatic index.

  double rhol; // Left fluid mass density.
  double ul; // Left fluid velocity.
  double pl; // Left fluid pressure.
  double alpha1_l; // Left fluid volume fraction (first species).

  double rhoc; // Central fluid mass density.
  double uc; // Central fluid velocity.
  double pc; // Central fluid pressure.
  double alpha1_c; // Central fluid volume fraction (first species).

  double rhor; // Right fluid mass density.
  double ur; // Right fluid velocity.
  double pr; // Right fluid pressure.
  double alpha1_r; // Central fluid volume fraction (first species).

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

struct amr_fedkiw_shock_ctx
create_ctx(void)
{
  // Physical constants (using normalized code units).
  double gas_gamma1 = 1.4; // First species adiabatic index.
  double gas_gamma2 = 1.67; // Second species adiabatic index.

  double rhol = 1.3333; // Left fluid mass density.
  double ul = 0.3535 * sqrt(pow(10.0, 5.0)); // Left fluid velocity.
  double pl = 1.5 * pow(10.0, 5.0); // Left fluid pressure.
  double alpha1_l = 0.99999; // Left fluid volume fraction (first species).

  double rhoc = 1.0; // Central fluid mass density.
  double uc = 0.0; // Central fluid velocity.
  double pc = 1.0 * pow(10.0, 5.0); // Central fluid pressure.
  double alpha1_c = 0.99999; // Central fluid volume fraction (first species).

  double rhor = 0.1379; // Right fluid mass density.
  double ur = 0.0; // Right fluid velocity.
  double pr = 1.0 * pow(10.0, 5.0); // Right fluid pressure.
  double alpha1_r = 0.00001; // Right fluid volume fraction (first species).

  // Simulation parameters.
  int Nx = 64; // Coarse cell count (x-direction).
  int ref_factor1 = 4; // First refinement factor (coarse-to-intermediate).
  int ref_factor2 = 4; // Second refinement factor (intermediate-to-fine).
  double Lx = 1.0; // Coarse domain size (x-direction).
  double intermediate_Lx = 0.5; // Intermediate domain size (x-direction).
  double fine_Lx = 0.2; // Fine domain size (x-direction).
  double cfl_frac = 0.95; // CFL coefficient.

  double t_end = 0.0012; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct amr_fedkiw_shock_ctx ctx = {
    .gas_gamma1 = gas_gamma1,
    .gas_gamma2 = gas_gamma2,
    .rhol = rhol,
    .ul = ul,
    .pl = pl,
    .alpha1_l = alpha1_l,
    .rhoc = rhoc,
    .uc = uc,
    .pc = pc,
    .alpha1_c = alpha1_c,
    .rhor = rhor,
    .ur = ur,
    .pr = pr,
    .alpha1_r = alpha1_r,
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
evalEulerMixtureInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  struct amr_fedkiw_shock_ctx new_ctx = create_ctx(); // Context for initialization functions.
  struct amr_fedkiw_shock_ctx *app = &new_ctx;

  double gas_gamma1 = app->gas_gamma1;
  double gas_gamma2 = app->gas_gamma2;

  double rhol = app->rhol;
  double ul = app->ul;
  double pl = app->pl;
  double alpha1_l = app->alpha1_l;

  double rhoc = app->rhoc;
  double uc = app->uc;
  double pc = app->pc;
  double alpha1_c = app->alpha1_c;

  double rhor = app->rhor;
  double ur = app->ur;
  double pr = app->pr;
  double alpha1_r = app->alpha1_r;

  double rho1 = 0.0;
  double rho2 = 0.0;
  double alpha1 = 0.0;

  double vx_total = 0.0;
  double vy_total = 0.0;
  double vz_total = 0.0;
  double p_total = 0.0;

  if (x < 0.05) {
    rho1 = rhol; // First species fluid mass density (left).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_l; // First species volume fraction (left).

    vx_total = ul; // Total mixture velocity (left).
    p_total = pl; // Total mixture pressure (left).
  }
  else if (x < 0.5) {
    rho1 = rhoc; // First species fluid mass density (central).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_c; // First species volume fraction (central).

    vx_total = uc; // Total mixture velocity (central).
    p_total = pc; // Total mixture pressure (central).
  }
  else {
    rho1 = rhoc; // First species fluid mass density (central).
    rho2 = rhor; // Second species fluid mass density (right).
    alpha1 = alpha1_r; // First species volume fraction (right).

    vx_total = ur; // Total mixture velocity (right).
    p_total = pr; // Total mixture pressure (right).
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
  struct amr_fedkiw_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = ctx.gas_gamma1;
  gas_gamma_s[1] = ctx.gas_gamma2;

  struct euler_mixture1d_double_init init = {
    .base_Nx = ctx.Nx,
    .ref_factor1 = ctx.ref_factor1,
    .ref_factor2 = ctx.ref_factor2,

    .coarse_x1 = 0.0,
    .coarse_x2 = ctx.Lx,

    .intermediate_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.intermediate_Lx),
    .intermediate_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.intermediate_Lx),

    .refined_x1 = (0.5 * ctx.Lx) - (0.5 * ctx.fine_Lx),
    .refined_x2 = (0.5 * ctx.Lx) + (0.5 * ctx.fine_Lx),

    .eval = evalEulerMixtureInit,
    .num_species = 2,
    .gas_gamma_s = gas_gamma_s,

    .euler_mixture_output = "amr_euler_mixture_fedkiw_shock_l2",

    .low_order_flux = true,
    .cfl_frac = ctx.cfl_frac,

    .t_end = ctx.t_end,
    .num_frames = ctx.num_frames,
    .dt_failure_tol = ctx.dt_failure_tol,
    .num_failures_max = ctx.num_failures_max,
  };

  euler_mixture1d_run_double(argc, argv, &init);

  gkyl_free(gas_gamma_s);
}