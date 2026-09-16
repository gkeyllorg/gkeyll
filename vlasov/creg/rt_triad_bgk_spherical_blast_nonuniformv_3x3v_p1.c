#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct spherical_blast_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double mass; // Neutral mass.
  double charge; // Neutral charge.

  double nl; // Left/inner number density.
  double Tl; // Left/inner temperature.
  double V_r_drift_l; // Left/inner drift velocity (radial direction).
  double V_theta_drift_l; // Left/inner drift velocity (angular direction).
  double V_phi_drift_l; // Left/inner drift velocity (phi-direction).


  double nr; // Right/outer number density.
  double Tr; // Right/outer temperature.
  double V_r_drift_r; // Right/outer drift velocity (radial direction).
  double V_theta_drift_r; // Right/outer drift velocity (angular direction).
  double V_phi_drift_r; // Right/outer drift velocity (phi-direction).

  double vt; // Thermal velocity.
  double nu; // Collision frequency.

  // Simulation parameters.
  int Nr; // Cell count (configuration space: radial direction).
  int Ntheta; // Cell count (configuration space: theta direction).
  int Nphi; // Cell count (configuration space: phi direction).
  int Nvr; // Cell count (velocity space: radial direction).
  int Nvtheta; // Cell count (velocity space: angular direction).
  int Nvphi; // Cell count (velocity space: phi-direction).
  double Lr; // Domain size (configuration space: radial direction).
  double theta_lo; // Domain lower boundary (configuration space: theta direction).
  double theta_up; // Domain upper boundary (configuration space: theta direction).
  double phi_lo; // Domain lower boundary (configuration space: phi direction).
  double phi_up; // Domain upper boundary (configuration space: phi direction).
  double vr_max; // Domain boundary (velocity space: radial direction).
  double vtheta_max; // Domain boundary (velocity space: angular direction).
  double vphi_max; // Domain boundary (velocity space: phi-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  double r0; // Radius of the high-density region

};

struct spherical_blast_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double mass = 1.0; // Neutral mass.
  double charge = 0.0; // Neutral charge.

  // Inner Material
  double nl = 1.0; // Left/inner number density.
  double Tl = 5.0; // Left/inner temperature.
  double V_r_drift_l = 0.0; // Left/inner drift velocity (radial direction).
  double V_theta_drift_l = 0.0; // Left/inner drift velocity (angular direction).
  double V_phi_drift_l = 0.0; // Left/inner drift velocity (phi-direction).

  // Outer Material
  double nr = 1.0; // Right/outer number density.
  double Tr = 1.0; // Right/outer temperature.
  double V_r_drift_r = 0.0; // Right/outer drift velocity (radial direction).
  double V_theta_drift_r = 0.0; // Right/outer drift velocity (angular direction).
  double V_phi_drift_r = 0.0; // Right/outer drift velocity (phi-direction).

  double vt = 1.0; // Thermal velocity.
  double nu = 15000.0; // Collision frequency.

  // Simulation parameters.
  int Nr = 4; // Cell count (configuration space: radial direction).
  int Ntheta = 4; // Cell count (configuration space: theta direction).
  int Nphi = 4; // Cell count (configuration space: phi direction).
  int Nvr = 6; // Cell count (velocity space: radial direction).
  int Nvtheta = 6; // Cell count (velocity space: angular direction).
  int Nvphi = 6; // Cell count (velocity space: phi-direction).
  double Lr = 1.5; // Domain size (configuration space: radial direction).
  double theta_lo = 0.25 * pi; // Domain lower boundary (configuration space: theta direction).
  double theta_up = 0.75 * pi; // Domain upper boundary (configuration space: theta direction).
  double phi_lo = 0.0; // Domain lower boundary (configuration space: phi direction).
  double phi_up = 0.5 * pi; // Domain upper boundary (configuration space: phi direction).
  double vr_max = 10.0 * vt; // Domain boundary (velocity space: radial direction).
  double vtheta_max = 10.0 * vt; // Domain boundary (velocity space: angular direction).
  double vphi_max = 10.0 * vt; // Domain boundary (velocity space: phi-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 0.001; // Final simulation time ( reference value = 0.7 ).
  int num_frames = 1; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = INT_MAX; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps. 

  double r0 = 0.2; // Radius of the high-density region

  struct spherical_blast_ctx ctx = {
    .pi = pi,
    .mass = mass,
    .charge = charge,
    .nl = nl,
    .Tl = Tl,
    .V_r_drift_l = V_r_drift_l,
    .V_theta_drift_l = V_theta_drift_l,
   .V_phi_drift_l = V_phi_drift_l,
    .nr = nr,
    .Tr = Tr,
    .V_r_drift_r = V_r_drift_r,
    .V_theta_drift_r = V_theta_drift_r,
    .V_phi_drift_r = V_phi_drift_r,
    .vt = vt,
    .nu = nu,
    .Nr = Nr,
    .Ntheta = Ntheta,
    .Nphi = Nphi,
    .Nvr = Nvr,
    .Nvtheta = Nvtheta,
    .Nvphi = Nvphi,
    .Lr = Lr,
    .theta_lo = theta_lo,
    .theta_up = theta_up,
    .phi_lo = phi_lo,
    .phi_up = phi_up,
    .vr_max = vr_max,
    .vtheta_max = vtheta_max,
    .vphi_max = vphi_max,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .r0 = r0,
  };

  return ctx;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;
  double r = xn[0];
  double theta = xn[1];
  double r0 = app->r0;
  double nl = app->nl;
  double nr = app->nr;
  double n = 0.0;

  if (r < r0) {
    n = nl; // Total number density (left/inner).
  }
  else {
    n = nr; // Total number density (right/outer).
  }

  double metric_det = r * r * sin(theta); // Jacobian for spherical coordinates 

  // Set total number density.
  fout[0] = metric_det * n;
}

void
evalTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;
  double r = xn[0];
  double r0 = app->r0;

  double Tl = app->Tl;
  double Tr = app->Tr;

  double T = 0.0;

  if (r < r0) {
    T = Tl; // Isotropic temperature (left/inner).
  }
  else {
    T = Tr; // Isotropic temperature (right/outer).
  }

  // Set isotropic temperature.
  fout[0] = T;
}

void
evalVDriftInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;
  double r = xn[0];
  double r0 = app->r0;

  double V_r_drift_l = app->V_r_drift_l;
  double V_theta_drift_l = app->V_theta_drift_l;
  double V_phi_drift_l = app->V_phi_drift_l;

  double V_r_drift_r = app->V_r_drift_r;
  double V_theta_drift_r = app->V_theta_drift_r;
  double V_phi_drift_r = app->V_phi_drift_r;

  double V_r_drift = 0.0;
  double V_theta_drift = 0.0;
  double V_phi_drift = 0.0;

  if (r < r0) {
    V_r_drift = V_r_drift_l; // Radial drift velocity (left/inner).
    V_theta_drift = V_theta_drift_l; // Angular drift velocity (left/inner).
    V_phi_drift = V_phi_drift_l; // phi drift velocity (left/inner).
  }
  else {
    V_r_drift = V_r_drift_r; // Radial drift velocity (right/outer).
    V_theta_drift = V_theta_drift_r; // Angular drift velocity (right/outer).
    V_phi_drift = V_phi_drift_r; // phi drift velocity (right/outer).
  }

  // Set total drift velocity.
  fout[0] = V_r_drift; fout[1] = V_theta_drift; fout[2] = V_phi_drift; 
}

void
evalNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;

  double nu = app->nu;

  // Set collision frequency.
  fout[0] = nu;
}

// Nonuniform velocity map: symmetric refinement toward v = 0,
// v_phys = v_c * (1 + s * (v_c/vmax)^2) / (1 + s), monotonic for s > -1/3.
struct vel_map_ctx { double vmax; double s; };

static void
mapc2p_vel_r(double t, const double *xc, double *xp, void *ctx)
{
  struct vel_map_ctx *vctx = ctx;
  double vc = xc[0], vm = vctx->vmax, s = vctx->s;
  xp[0] = vc*(1.0 + s*(vc/vm)*(vc/vm))/(1.0 + s);
}

void
evalCovTangentBasis(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;

  double q_r = xn[0];
  double q_theta = xn[1];
  double q_phi = xn[2];

  // vals_ij = e_{x} . \sigma_{x}
  double e_rx = cos(q_phi) * sin(q_theta); // Covariant Tangent Basis Coefficients (r-x coefficient).
  double e_ry = sin(q_phi) * sin(q_theta); // Covariant Tangent Basis Coefficients (r-y coefficient).
  double e_rz = cos(q_theta); // Covariant Tangent Basis Coefficients (r-z coefficient).
  double e_tx = q_r * cos(q_phi) * cos(q_theta); // Covariant Tangent Basis Coefficients (theta-x coefficient).
  double e_ty = q_r * sin(q_phi) * cos(q_theta); // Covariant Tangent Basis Coefficients (theta-y coefficient).
  double e_tz = - q_r * sin(q_theta); // Covariant Tangent Basis Coefficients (theta-z coefficient).
  double e_px = - q_r * sin(q_phi) * sin(q_theta); // Covariant Tangent Basis Coefficients (phi-x coefficient).
  double e_py = q_r * cos(q_phi) * sin(q_theta); // Covariant Tangent Basis Coefficients (phi-y coefficient).
  double e_pz = 0.0; // Covariant Tangent Basis Coefficients (phi-z coefficient).


  fout[0] = e_rx;
  fout[1] = e_ry;
  fout[2] = e_rz;
  fout[3] = e_tx;
  fout[4] = e_ty;
  fout[5] = e_tz;
  fout[6] = e_px;
  fout[7] = e_py;
  fout[8] = e_pz;
}

void
evalTriadBasis(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;

  double q_r = xn[0];
  double q_theta = xn[1];
  double q_phi = xn[2];

  // vals_ij = \sigma_{x} . \sigma_{x}
  double sigma_rx = cos(q_phi) * sin(q_theta); // Triad Basis Coefficients (r-x coefficient).
  double sigma_ry = sin(q_phi) * sin(q_theta); // Triad Basis Coefficients (r-y coefficient).
  double sigma_rz = cos(q_theta); // Triad Basis Coefficients (r-z coefficient).
  double sigma_tx = cos(q_phi) * cos(q_theta); // Triad Basis Coefficients (theta-x coefficient).
  double sigma_ty = sin(q_phi) * cos(q_theta); // Triad Basis Coefficients (theta-y coefficient).
  double sigma_tz = - sin(q_theta); // Triad Basis Coefficients (theta-z coefficient).
  double sigma_px = - sin(q_phi); // Triad Basis Coefficients (phi-x coefficient).
  double sigma_py = cos(q_phi); // Triad Basis Coefficients (phi-y coefficient).
  double sigma_pz = 0.0; // Triad Basis Coefficients (phi-z coefficient).
  
  fout[0] = sigma_rx;
  fout[1] = sigma_ry;
  fout[2] = sigma_rz;
  fout[3] = sigma_tx;
  fout[4] = sigma_ty;
  fout[5] = sigma_tz;
  fout[6] = sigma_px;
  fout[7] = sigma_py;
  fout[8] = sigma_pz;
}

void
evalTriadBasisGradient(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct spherical_blast_ctx *app = ctx;

  double q_r = xn[0];
  double q_theta = xn[1];
  double q_phi = xn[2];

  // d(vals_ij)/dx^k = (\sigma_{x} . \sigma_{x})/dx^k
  // d/dr components
  double d_sigma_rx_dr = 0.0; // Triad Basis Gradient Coefficients (r-x coefficient).
  double d_sigma_ry_dr = 0.0; // Triad Basis Gradient Coefficients (r-y coefficient).
  double d_sigma_rz_dr = 0.0; // Triad Basis Gradient Coefficients (r-z coefficient).
  double d_sigma_tx_dr = 0.0; // Triad Basis Gradient Coefficients (theta-x coefficient).
  double d_sigma_ty_dr = 0.0; // Triad Basis Gradient Coefficients (theta-y coefficient).
  double d_sigma_tz_dr = 0.0; // Triad Basis Gradient Coefficients (theta-z coefficient).
  double d_sigma_px_dr = 0.0; // Triad Basis Gradient Coefficients (phi-x coefficient).
  double d_sigma_py_dr = 0.0; // Triad Basis Gradient Coefficients (phi-y coefficient).
  double d_sigma_pz_dr = 0.0; // Triad Basis Gradient Coefficients (phi-z coefficient).

  // d/dtheta components
  double d_sigma_rx_dtheta = cos(q_phi) * cos(q_theta); // Triad Basis Gradient Coefficients (r-x coefficient).
  double d_sigma_ry_dtheta = sin(q_phi) * cos(q_theta); // Triad Basis Gradient Coefficients (r-y coefficient).
  double d_sigma_rz_dtheta = - sin(q_theta); // Triad Basis Gradient Coefficients (r-z coefficient).
  double d_sigma_tx_dtheta = - cos(q_phi) * sin(q_theta); // Triad Basis Gradient Coefficients (theta-x coefficient).
  double d_sigma_ty_dtheta = - sin(q_phi) * sin(q_theta); // Triad Basis Gradient Coefficients (theta-y coefficient).
  double d_sigma_tz_dtheta = - cos(q_theta); // Triad Basis Gradient Coefficients (theta-z coefficient).
  double d_sigma_px_dtheta = 0.0; // Triad Basis Gradient Coefficients (phi-x coefficient).
  double d_sigma_py_dtheta = 0.0; // Triad Basis Gradient Coefficients (phi-y coefficient).
  double d_sigma_pz_dtheta = 0.0; // Triad Basis Gradient Coefficients (phi-z coefficient).
  
  // d/dphi components
  double d_sigma_rx_dp = - sin(q_phi) * sin(q_theta); // Triad Basis Gradient Coefficients (r-x coefficient).
  double d_sigma_ry_dp = cos(q_phi) * sin(q_theta); // Triad Basis Gradient Coefficients (r-y coefficient).
  double d_sigma_rz_dp = 0.0; // Triad Basis Gradient Coefficients (r-z coefficient).
  double d_sigma_tx_dp = - sin(q_phi) * cos(q_theta); // Triad Basis Gradient Coefficients (theta-x coefficient).
  double d_sigma_ty_dp = cos(q_phi) * cos(q_theta); // Triad Basis Gradient Coefficients (theta-y coefficient).
  double d_sigma_tz_dp = 0.0; // Triad Basis Gradient Coefficients (theta-z coefficient).
  double d_sigma_px_dp = - cos(q_phi); // Triad Basis Gradient Coefficients (phi-x coefficient).
  double d_sigma_py_dp = - sin(q_phi); // Triad Basis Gradient Coefficients (phi-y coefficient).
  double d_sigma_pz_dp = 0.0; // Triad Basis Gradient Coefficients (phi-z coefficient).


  fout[0] = d_sigma_rx_dr;
  fout[1] = d_sigma_ry_dr;
  fout[2] = d_sigma_rz_dr;
  fout[3] = d_sigma_tx_dr;
  fout[4] = d_sigma_ty_dr;
  fout[5] = d_sigma_tz_dr;
  fout[6] = d_sigma_px_dr;
  fout[7] = d_sigma_py_dr;
  fout[8] = d_sigma_pz_dr;

  fout[9] = d_sigma_rx_dtheta;
  fout[10] = d_sigma_ry_dtheta;
  fout[11] = d_sigma_rz_dtheta;
  fout[12] = d_sigma_tx_dtheta;
  fout[13] = d_sigma_ty_dtheta;
  fout[14] = d_sigma_tz_dtheta;
  fout[15] = d_sigma_px_dtheta;
  fout[16] = d_sigma_py_dtheta;
  fout[17] = d_sigma_pz_dtheta;

  fout[18] = d_sigma_rx_dp;
  fout[19] = d_sigma_ry_dp;
  fout[20] = d_sigma_rz_dp;
  fout[21] = d_sigma_tx_dp;
  fout[22] = d_sigma_ty_dp;
  fout[23] = d_sigma_tz_dp;
  fout[24] = d_sigma_px_dp;
  fout[25] = d_sigma_py_dp;
  fout[26] = d_sigma_pz_dp;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, frame);
    gkyl_vlasov_app_write_field_energy(app);
    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_integrated_L2_f(app);
    gkyl_vlasov_app_write_mom(app, t_curr, frame);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct spherical_blast_ctx ctx = create_ctx(); // Context for initialization functions.

  int NR = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nr);
  int NTHETA = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.Ntheta);
  int NPHI = APP_ARGS_CHOOSE(app_args.xcells[2], ctx.Nphi);
  int NVR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvr);
  int NVTHETA = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvtheta);
  int NVPHI = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.Nvphi);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NR, NTHETA, NPHI };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }



  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Neutral species.
  struct vel_map_ctx vmap_ctx_r = { .vmax = ctx.vr_max, .s = 1.5 };
  struct vel_map_ctx vmap_ctx_theta = { .vmax = ctx.vtheta_max, .s = 1.5 };
  struct vel_map_ctx vmap_ctx_phi = { .vmax = ctx.vphi_max, .s = 1.5 };

  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .model_id = GKYL_MODEL_TRIAD,
    .charge = ctx.charge, .mass = ctx.mass,
    .lower = { -ctx.vr_max, -ctx.vtheta_max, -ctx.vphi_max,},
    .upper = { ctx.vr_max, ctx.vtheta_max, ctx.vphi_max },
    .cells = { NVR, NVTHETA, NVPHI },

    .mapc2p_vel = {
      { .mapc2p_vel_func = mapc2p_vel_r, .mapc2p_vel_ctx = &vmap_ctx_r },
      { .mapc2p_vel_func = mapc2p_vel_r, .mapc2p_vel_ctx = &vmap_ctx_theta },
      { .mapc2p_vel_func = mapc2p_vel_r, .mapc2p_vel_ctx = &vmap_ctx_phi },
    },
    .cov_tangent_basis = evalCovTangentBasis,
    .cov_tangent_basis_ctx = &ctx,
    .triad_basis = evalTriadBasis,
    .triad_basis_ctx = &ctx,
    .triad_basis_gradient = evalTriadBasisGradient,
    .triad_basis_gradient_ctx = &ctx,

    .num_init = 1, 
    .projection[0] = {
      .proj_id = GKYL_PROJ_VLASOV_LTE,
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .temp = evalTempInit,
      .ctx_temp = &ctx,
      .V_drift = evalVDriftInit,
      .ctx_V_drift = &ctx,
      .correct_all_moms = true,
      .iter_eps = 0.0,
      .max_iter = 0,
      .use_last_converged = false,
    },
    
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalNu,
      .self_nu_ctx = &ctx,
      .is_implicit = true,
    },

    .correct = {
      .correct_all_moms = true,
      .iter_eps = 1.0e-12,
      .max_iter = 100,
      .use_last_converged = false,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_REFLECT, },
    },
    .bcy = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_REFLECT, },
    },
    .bcz = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_REFLECT, },
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_LTE, GKYL_F_MOMENT_ENERGY },
  };

  // Vlasov-Maxwell app.
  struct gkyl_vm app_inp = {
   .name = "triad_bgk_spherical_blast_nonuniformv_3x3v_p1",

   .cdim = 3, .vdim = 3,
   .lower = { 0.05, ctx.theta_lo, ctx.phi_lo },
   .upper = { 0.05 + ctx.Lr, ctx.theta_up, ctx.phi_up },
   .cells = { NR, NTHETA, NPHI },

   .poly_order = ctx.poly_order,
   .basis_type = app_args.basis_type,
   .cfl_frac = ctx.cfl_frac,

   .num_periodic_dir = 0,
   .periodic_dirs = { },

   .num_species = 1,
   .species = { neut },

   .skip_field = true,

   .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  double io_tm =  stat.field_io_tm + stat.species_io_tm + stat.field_diag_io_tm + stat.species_diag_io_tm;
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", io_tm);


freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
