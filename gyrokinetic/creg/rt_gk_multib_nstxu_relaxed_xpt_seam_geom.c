#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_efit.h>
#include <gkyl_gk_block_geom.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

// Regression coverage for the relaxed X-point seam diagnostic optimizer and
// guarded production apply (improve-xpt-finder branch, commits "Optimize
// smooth bounded delta-s for X-point interfaces" and "Use relaxed X-point
// seams with guarded fallback").
//
// Builds a minimal, self-contained 2-block half-domain PF seam pair -- the
// smallest X-point seam pair type -- for one of three NSTX-U shots, so the
// optimizer's pair-discovery/search/apply machinery is exercised end to end
// without the unrelated complexity of a full 8-block mesh or any
// species/velocity-space setup. Geometry-only: gkyl_gyrokinetic_multib_app_new_geom
// is what already runs the optimizer internally during block construction;
// no species are needed to exercise it.
//
// Shots (geqdsk data bundled under gyrokinetic/data/eqdsk/):
//   203877 ms655 (default) -- known difficult physical-boundary case
//   203742 ms480           -- smoother control
//   203730 ms839            -- used at higher resolution for folded-cell
//                              detection
//
// Usage: -o "shot=203742" selects a shot (default 203877); -o "bound=0.001"
// overrides the delta-s bound (default 0.005 m); -o "optimize=0" disables
// relaxed_xpt_seam_optimize on both blocks (so the CSV output can be
// compared to the exact same zero-coefficient straight-map construction).
// Standard -x/-y override the psi/theta cell counts (default 8x6 per
// block); -M requests MPI; -d requests theta ("z") decomposition. Radial
// decomposition (-c) is intentionally rejected by the optimizer itself, not
// by this test -- see the unsupported_configuration guard added for
// cuts[0] != 1 in gyrokinetic_multib.c.

// ---- plate shape helpers, verbatim from the validated NSTX-U elong_tri_scan
// research harness (private-flux boundary plates near the lower divertor) ----
static void
extend_segment_eval(double R0, double R1, double Z0, double Z1,
  double alpha, double s, double *R, double *Z)
{
  double dR = R1-R0, dZ = Z1-Z0;
  double expand = (alpha-1.0)/2.0;
  double R0n = R0-expand*dR;
  double R1n = R1+expand*dR;
  double Z0n = Z0-expand*dZ;
  double Z1n = Z1+expand*dZ;
  *R = R0n+(R1n-R0n)*s;
  *Z = Z0n+(Z1n-Z0n)*s;
}

static void
shaped_pfunc_lower_outer(double s, double *RZ)
{
  double R, Z;
  if (s < 1.0/3.0)
    extend_segment_eval(0.415, 0.5715, -1.6234, -1.6234, 1.0, s*3.0, &R, &Z);
  else if (s < 2.0/3.0)
    extend_segment_eval(0.5715, 0.617, -1.6234, -1.628, 1.0,
      (s-1.0/3.0)*3.0, &R, &Z);
  else
    extend_segment_eval(0.617, 1.194, -1.628, -1.401, 1.0,
      (s-2.0/3.0)*3.0, &R, &Z);
  RZ[0] = R; RZ[1] = Z;
}

static void
shaped_pfunc_lower_inner(double s, double *RZ)
{
  double R, Z;
  if (s < 1.0/3.0)
    extend_segment_eval(0.3148, 0.3148, -0.6, -1.05, 1.0, s*3.0, &R, &Z);
  else if (s < 2.0/3.0)
    extend_segment_eval(0.3148, 0.415, -1.05, -1.27, 1.0,
      (s-1.0/3.0)*3.0, &R, &Z);
  else
    extend_segment_eval(0.415, 0.415, -1.27, -1.6234, 1.0,
      (s-2.0/3.0)*3.0, &R, &Z);
  RZ[0] = R; RZ[1] = Z;
}

// Poloidal flux at normalized radial coordinate rho (NSTX-U convention).
static double
psi_rho(double rho, double psi_axis, double psi_sep)
{
  return rho*rho*(psi_sep-psi_axis)+psi_axis;
}

struct shot_case {
  int shot, time_ms;
  const char *geqdsk;
};

static struct shot_case
select_shot_case(const char *opt_args)
{
  int shot = 0;
  if (opt_args && opt_args[0])
    sscanf(opt_args, "shot=%d", &shot);
  if (shot == 203742)
    return (struct shot_case) {
      203742, 480, "gyrokinetic/data/eqdsk/nstxu_shot203742_ms480.geqdsk"
    };
  if (shot == 203730)
    return (struct shot_case) {
      203730, 839, "gyrokinetic/data/eqdsk/nstxu_shot203730_ms839.geqdsk"
    };
  return (struct shot_case) {
    203877, 655, "gyrokinetic/data/eqdsk/nstxu_shot203877_ms655.geqdsk"
  };
}

static double
opt_get_double(const char *opt_args, const char *key, double def)
{
  if (!opt_args) return def;
  const char *p = strstr(opt_args, key);
  if (!p || p[strlen(key)] != '=') return def;
  return atof(p+strlen(key)+1);
}

static bool
opt_get_bool(const char *opt_args, const char *key, bool def)
{
  if (!opt_args) return def;
  const char *p = strstr(opt_args, key);
  if (!p || p[strlen(key)] != '=') return def;
  return atoi(p+strlen(key)+1) != 0;
}

static struct gkyl_gk_block_geom*
create_gk_block_geom(struct shot_case shot, int npsi_pf, int ntheta,
  double bound, bool optimize, int cuts_psi, int cuts_theta)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 2);

  struct gkyl_efit_inp efit_inp = {
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  strncpy(efit_inp.filepath, shot.geqdsk, sizeof(efit_inp.filepath)-1);
  efit_inp.filepath[sizeof(efit_inp.filepath)-1] = '\0';

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psimag = efit->simag, psisep = efit->psisep_cubic; // the flux of the representation traced on (C1, 2026-10-05)
  printf("Relaxed X-point seam regression: shot=%d ms=%d psisep=%.6g "
    "psimag=%.6g\n", shot.shot, shot.time_ms, psisep, psimag);
  gkyl_efit_release(efit);

  // Pseudo private-flux rho boundary, matching the validated elong_tri_scan
  // NSTX-U configuration.
  double pseudo_rho_pf = 0.992;
  double psi_lo_pf = psi_rho(pseudo_rho_pf, psimag, psisep);
  double psi_up_pf = psisep;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;
  double rright_out = 1.6, zinner = 1.22, zouter = 1.5;

  // Block 0: PF_LO_R (lower outer private-flux region). Its upper theta
  // edge is the X-point seam shared with block 1; every other edge (both
  // radial edges, which in the full 8-block mesh connect to a SOL/core
  // block, and the lower theta edge, which is the divertor plate) is
  // physical here since this reduced mesh only needs the one seam pair.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo },
      .upper = { psi_up_pf, theta_up },
      .cells = { npsi_pf, ntheta },
      .cuts = { cuts_psi, cuts_theta },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_R,
          .use_cubics = true, // bicubic (C1) psi for tracing, every block (2026-10-05)
          .half_domain = true,
          .straight_xpt_ray = true,
          .relaxed_xpt_seam = true,
          .relaxed_xpt_seam_sweep = false,
          .relaxed_xpt_seam_delta_s_coeff = 0.0,
          .relaxed_xpt_seam_delta_s_bound = bound,
          .relaxed_xpt_seam_optimize = optimize,
          .rright = rright_out,
          .rleft = 0.3,
          .rmin = 0.29,
          .rmax = 1.6,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      .connections[0] = {
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
      },
      .connections[1] = {
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE },
      }
    }
  );

  // Block 1: PF_LO_L (lower inner private-flux region), the reciprocal seam
  // partner of block 0.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo },
      .upper = { psi_up_pf, theta_up },
      .cells = { npsi_pf, ntheta },
      .cuts = { cuts_psi, cuts_theta },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_L,
          .use_cubics = true, // bicubic (C1) psi for tracing, every block (2026-10-05)
          .half_domain = true,
          .straight_xpt_ray = true,
          .relaxed_xpt_seam = true,
          .relaxed_xpt_seam_sweep = false,
          .relaxed_xpt_seam_delta_s_coeff = 0.0,
          .relaxed_xpt_seam_delta_s_bound = bound,
          .relaxed_xpt_seam_optimize = optimize,
          .rright = rright_out,
          .rleft = 0.3,
          .rmin = 0.29,
          .rmax = 1.6,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      .connections[0] = {
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL },
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL },
      },
      .connections[1] = {
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL },
      }
    }
  );

  return bgeom;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  struct shot_case shot = select_shot_case(app_args.opt_args);
  double bound = opt_get_double(app_args.opt_args, "bound", 0.005);
  bool optimize = opt_get_bool(app_args.opt_args, "optimize", true);

  int npsi_pf = APP_ARGS_CHOOSE(app_args.xcells[0], 8);
  int ntheta = APP_ARGS_CHOOSE(app_args.xcells[1], 6);

  printf("  npsi=%d ntheta=%d bound=%.6g optimize=%d cuts=(%d,%d)\n",
    npsi_pf, ntheta, bound, optimize, app_args.cuts[0], app_args.cuts[1]);

  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom(shot, npsi_pf,
    ntheta, bound, optimize, app_args.cuts[0], app_args.cuts[1]);

  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi,
    app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_multib mbinp = {
    .cdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .gk_block_geom = bgeom,
    .cfl_frac = 1.0,
    .num_periodic_dir = 0,
    .num_species = 0,
    .num_neut_species = 0,
    .comm = comm,
  };
  snprintf(mbinp.name, sizeof(mbinp.name), "nstxu_relaxed_xpt_shot%d_ms%d",
    shot.shot, shot.time_ms);

  gkyl_gyrokinetic_multib_app *app =
    gkyl_gyrokinetic_multib_app_new_geom(&mbinp);
  int status = 0;
  if (!app) {
    fprintf(stderr, "ERROR: failed to construct multiblock geometry app\n");
    status = 1;
  }
  else {
    gkyl_gyrokinetic_multib_app_write_geometry(app);
    gkyl_gyrokinetic_multib_app_release_geom(app);
  }

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  return status;
}
