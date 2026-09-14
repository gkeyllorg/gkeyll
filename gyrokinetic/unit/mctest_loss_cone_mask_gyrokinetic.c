#include <acutest.h>

#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_null_comm.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#endif

struct lc_ctx {
  double z_max;
  double b_m;
  double r_m;
  double t0;
  double eV;
  double phi_fac;
};

typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

static struct gkyl_array*
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static struct gkyl_comm*
comm_new(struct gkyl_rect_decomp *decomp)
{
#ifdef GKYL_HAVE_MPI
  return gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp) {
    .mpi_comm = MPI_COMM_WORLD,
    .decomp = decomp,
  });
#else
  return gkyl_null_comm_inew(&(struct gkyl_null_comm_inp) { .use_gpu = false });
#endif
}

static void
bmag_func_1x(double t, const double *xn, double *fout, void *ctx)
{
  struct lc_ctx *p = ctx;
  double z = xn[0];
  fout[0] = p->b_m * (1.0 - ((p->r_m - 1.0) / p->r_m) * pow(cos(z), 2.0));
}

static void
phi_func_1x(double t, const double *xn, double *fout, void *ctx)
{
  struct lc_ctx *p = ctx;
  double z = xn[0];
  fout[0] = p->phi_fac * p->t0 / p->eV * (1.0 - pow(cos(4.0 * z), 2.0));
}

static void
eval_dg_on_range(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_range *range, evalf_t eval, void *ctx, struct gkyl_array *out)
{
  struct gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_new(grid, basis, 1, eval, ctx);
  gkyl_eval_on_nodes_advance(ev, 0.0, range, out);
  gkyl_eval_on_nodes_release(ev);
}

static void
test_loss_cone_mask_parallel_4dom(void)
{
  // The four ranks each evaluate one z slab using globally gathered B and
  // phi. The independently evaluated single-domain result is the oracle; an
  // exact cell-by-cell match checks that decomposition and allgather do not
  // change the mask.
  int rank = 0, size = 1;
#ifdef GKYL_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if (size != 4) {
    if (rank == 0) {
      TEST_MSG("mctest_loss_cone_mask_gyrokinetic requires 4 MPI ranks, got %d", size);
    }
    return;
  }

  const int cdim = 1, vdim = 2, pdim = 3;
  const int nz = 64, nvpar = 16, nmu = 16;
  const int cells_conf[] = { nz };
  const int cells_phase[] = { nz, nvpar, nmu };
  const int cells_vel[] = { nvpar, nmu };
  const int cuts[] = { 4 };
  const double z_max = GKYL_PI - 0.5;
  const double lower_phase[] = { -z_max, -6.0, 0.0 };
  const double upper_phase[] = { z_max, 6.0, 18.0 };
  const double lower_conf[] = { -z_max };
  const double upper_conf[] = { z_max };
  const double lower_vel[] = { -6.0, 0.0 };
  const double upper_vel[] = { 6.0, 18.0 };

  struct gkyl_rect_grid grid_conf, grid_phase, grid_vel;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  gkyl_rect_grid_init(&grid_phase, pdim, lower_phase, upper_phase, cells_phase);
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  struct gkyl_range conf_global, phase_global, vel_global;
  gkyl_create_global_range(cdim, cells_conf, &conf_global);
  gkyl_create_global_range(pdim, cells_phase, &phase_global);
  gkyl_create_global_range(vdim, cells_vel, &vel_global);

  struct gkyl_rect_decomp *conf_decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &conf_global);
  int phase_cuts[] = { 4, 1, 1 };
  struct gkyl_rect_decomp *phase_decomp = gkyl_rect_decomp_new_from_cuts(pdim, phase_cuts, &phase_global);
  struct gkyl_comm *comm_conf = comm_new(conf_decomp);
  struct gkyl_comm *comm_phase = comm_new(phase_decomp);

  const struct gkyl_range *conf_local = &conf_decomp->ranges[rank];
  const struct gkyl_range *phase_local = &phase_decomp->ranges[rank];

  struct gkyl_basis basis_conf, basis_phase;
  gkyl_cart_modal_serendip(&basis_conf, cdim, 1);
  gkyl_cart_modal_gkhybrid(&basis_phase, cdim, vdim);

  struct lc_ctx ctx = {
    .z_max = z_max,
    .b_m = 4.0,
    .r_m = 8.0,
    .t0 = 100.0 * GKYL_ELEMENTARY_CHARGE,
    .eV = GKYL_ELEMENTARY_CHARGE,
    .phi_fac = 5.0,
  };

  struct gkyl_array *bmag_local = mkarr(basis_conf.num_basis, conf_local->volume);
  struct gkyl_array *phi_local = mkarr(basis_conf.num_basis, conf_local->volume);
  struct gkyl_array *bmag_global_gather = mkarr(basis_conf.num_basis, conf_global.volume);
  struct gkyl_array *phi_global_gather = mkarr(basis_conf.num_basis, conf_global.volume);

  eval_dg_on_range(&grid_conf, &basis_conf, conf_local, bmag_func_1x, &ctx, bmag_local);
  eval_dg_on_range(&grid_conf, &basis_conf, conf_local, phi_func_1x, &ctx, phi_local);

  gkyl_comm_array_allgather(comm_conf, conf_local, &conf_global, bmag_local, bmag_global_gather);
  gkyl_comm_array_allgather(comm_conf, conf_local, &conf_global, phi_local, phi_global_gather);

  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm_local = gkyl_velocity_map_new(c2p_in, grid_phase, grid_vel,
    *phase_local, *phase_local, vel_global, vel_global, false);

  struct gkyl_loss_cone_mask_gyrokinetic *up_local =
    gkyl_loss_cone_mask_gyrokinetic_inew(&(struct gkyl_loss_cone_mask_gyrokinetic_inp) {
      .conf_basis = &basis_conf,
      .vel_map = gvm_local,
      .mass = 2.014 * GKYL_PROTON_MASS,
      .charge = GKYL_ELEMENTARY_CHARGE,
    });

  struct gkyl_array *mask_local = mkarr(1, phase_local->volume);
  gkyl_loss_cone_mask_gyrokinetic_advance(up_local, phase_local, &conf_global,
    bmag_global_gather, phi_global_gather, 0, 0, mask_local);

  struct gkyl_array *mask_dist_global = mkarr(1, phase_global.volume);
  gkyl_comm_array_allgather(comm_phase, phase_local, &phase_global, mask_local, mask_dist_global);

  if (rank == 0) {
    struct gkyl_array *bmag_ref = mkarr(basis_conf.num_basis, conf_global.volume);
    struct gkyl_array *phi_ref = mkarr(basis_conf.num_basis, conf_global.volume);
    eval_dg_on_range(&grid_conf, &basis_conf, &conf_global, bmag_func_1x, &ctx, bmag_ref);
    eval_dg_on_range(&grid_conf, &basis_conf, &conf_global, phi_func_1x, &ctx, phi_ref);

    struct gkyl_velocity_map *gvm_global = gkyl_velocity_map_new(c2p_in, grid_phase, grid_vel,
      phase_global, phase_global, vel_global, vel_global, false);
    struct gkyl_loss_cone_mask_gyrokinetic *up_ref =
      gkyl_loss_cone_mask_gyrokinetic_inew(&(struct gkyl_loss_cone_mask_gyrokinetic_inp) {
        .conf_basis = &basis_conf,
        .vel_map = gvm_global,
        .mass = 2.014 * GKYL_PROTON_MASS,
        .charge = GKYL_ELEMENTARY_CHARGE,
      });

    struct gkyl_array *mask_ref = mkarr(1, phase_global.volume);
    gkyl_loss_cone_mask_gyrokinetic_advance(up_ref, &phase_global, &conf_global,
      bmag_ref, phi_ref, 0, 0, mask_ref);

    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &phase_global);
    while (gkyl_range_iter_next(&it)) {
      long linidx = gkyl_range_idx(&phase_global, it.idx);
      const double *md = gkyl_array_cfetch(mask_dist_global, linidx);
      const double *mr = gkyl_array_cfetch(mask_ref, linidx);
      TEST_CHECK(fabs(md[0] - mr[0]) < 1e-12);
    }

    gkyl_array_release(mask_ref);
    gkyl_loss_cone_mask_gyrokinetic_release(up_ref);
    gkyl_velocity_map_release(gvm_global);
    gkyl_array_release(phi_ref);
    gkyl_array_release(bmag_ref);
  }

  gkyl_array_release(mask_dist_global);
  gkyl_array_release(mask_local);
  gkyl_loss_cone_mask_gyrokinetic_release(up_local);
  gkyl_velocity_map_release(gvm_local);
  gkyl_array_release(phi_global_gather);
  gkyl_array_release(bmag_global_gather);
  gkyl_array_release(phi_local);
  gkyl_array_release(bmag_local);
  gkyl_comm_release(comm_phase);
  gkyl_comm_release(comm_conf);
  gkyl_rect_decomp_release(phase_decomp);
  gkyl_rect_decomp_release(conf_decomp);
}

TEST_LIST = {
  { "loss_cone_mask_parallel_4dom", test_loss_cone_mask_parallel_4dom },
  { NULL, NULL },
};
