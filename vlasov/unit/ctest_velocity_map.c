#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_velocity_map.h>

#include <gkyl_comm.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
#else
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
#endif
  return a;
}

void test_vmap_1x2v_p1_mapc2p_vel_vpar(double t, const double *zc, double* GKYL_RESTRICT vp, void *ctx)
{
  double vparc = zc[0];
  vp[0] = vparc<0.? -pow(vparc,2) : pow(vparc,2); // Quadratic mapping.
}
void test_vmap_1x2v_p1_mapc2p_vel_mu(double t, const double *zc, double* GKYL_RESTRICT vp, void *ctx)
{
  double muc = zc[0];
  vp[0] = pow(muc,2); // Quadratic mapping.
}
void test_vmap_1x2v_p1_mapc2p_vel(double t, const double *zc, double* GKYL_RESTRICT vp, void *ctx)
{
  double vparc[] = {zc[0]}, muc[] = {zc[1]};
  double vparp[1], mup[1];
  test_vmap_1x2v_p1_mapc2p_vel_vpar(t, vparc, vparp, ctx);
  test_vmap_1x2v_p1_mapc2p_vel_mu(t, muc, mup, ctx);

  vp[0] = vparp[0];
  vp[1] = mup[0];
}

void
test_vmap_1x2v_p1(bool use_gpu)
{
  int poly_order = 1;
  double lower[] = {-M_PI, -1.0, 0.0}, upper[] = {M_PI, 1.0, 1.0};
  int cells[] = {2, 12, 6};
  int vdim = 2;
  int pdim = sizeof(lower)/sizeof(lower[0]);
  int cdim = pdim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = cells[cdim+d];
  }

  // Grids.
  struct gkyl_rect_grid grid, grid_conf, grid_vel;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis, basis_conf;
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost[pdim];
  for (int d=0; d<pdim; d++) ghost[d] = 0;
  for (int d=0; d<cdim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);
  int ghost_vel[vdim];
  for (int d=0; d<vdim; d++) ghost_vel[d] = 0;
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  struct gkyl_mapc2p_inp c2p_in = {
    .mapping = test_vmap_1x2v_p1_mapc2p_vel,
    .ctx = NULL,
  };
  
  // Velocity space mapping.
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);

  // Check vmap.
  struct gkyl_array *vmap_ho = mkarr(false, gvm->vmap->ncomp, gvm->vmap->size);
  gkyl_array_copy(vmap_ho, gvm->vmap);
  struct gkyl_range_iter iter;
  double lower_ref[1], upper_ref[1];
  int cells_ref[1], ghost_ref[1];
  struct gkyl_rect_grid grid_ref;
  struct gkyl_range local_ref, local_ext_ref; // 1D local, local-ext velocity-space ranges
  struct gkyl_basis basis_ref;
  gkyl_cart_modal_serendip(&basis_ref, 1, 1);
  struct gkyl_array *vmap_ref;

  // Check vpar mapping.
  lower_ref[0] = lower[cdim+0];
  upper_ref[0] = upper[cdim+0];
  cells_ref[0] = cells[cdim+0];
  gkyl_rect_grid_init(&grid_ref, 1, lower_ref, upper_ref, cells_ref);
  ghost_ref[0] = ghost_vel[0];
  gkyl_create_grid_ranges(&grid_ref, ghost_ref, &local_ext_ref, &local_ref);
  vmap_ref = mkarr(false, basis_ref.num_basis, local_ext_ref.volume);
  struct gkyl_eval_on_nodes* evOnNod_vpar = gkyl_eval_on_nodes_new(&grid_ref, &basis_ref,
    1, test_vmap_1x2v_p1_mapc2p_vel_vpar, NULL);
  gkyl_eval_on_nodes_advance(evOnNod_vpar, 0., &local_ref, vmap_ref);
  gkyl_eval_on_nodes_release(evOnNod_vpar);

  gkyl_range_iter_init(&iter, &local_vel);
  while (gkyl_range_iter_next(&iter)) {
    int idx_ref[] = {iter.idx[0]};
    long loc = gkyl_range_idx(&local_vel, iter.idx);
    long loc_ref = gkyl_range_idx(&local_ref, idx_ref);
    double *vm = gkyl_array_fetch(vmap_ho, loc);
    double *vm_ref = gkyl_array_fetch(vmap_ref, loc_ref);
    for (int i=0; i<basis_ref.num_basis; i++) {
      TEST_CHECK( gkyl_compare(vm[i], vm_ref[i], 1e-12) );
      TEST_MSG("Produced: %.13e | Expected: %.13e\n",vm[i], vm_ref[i]);
    }
  }
  gkyl_array_release(vmap_ref);

  // Check mu mapping.
  lower_ref[0] = lower[cdim+1];
  upper_ref[0] = upper[cdim+1];
  cells_ref[0] = cells[cdim+1];
  gkyl_rect_grid_init(&grid_ref, 1, lower_ref, upper_ref, cells_ref);
  ghost_ref[0] = ghost_vel[1];
  gkyl_create_grid_ranges(&grid_ref, ghost_ref, &local_ext_ref, &local_ref);
  vmap_ref = mkarr(false, basis_ref.num_basis, local_ext_ref.volume);
  struct gkyl_eval_on_nodes* evOnNod_mu = gkyl_eval_on_nodes_new(&grid_ref, &basis_ref,
    1, test_vmap_1x2v_p1_mapc2p_vel_mu, NULL);
  gkyl_eval_on_nodes_advance(evOnNod_mu, 0., &local_ref, vmap_ref);
  gkyl_eval_on_nodes_release(evOnNod_mu);

  gkyl_range_iter_init(&iter, &local_vel);
  while (gkyl_range_iter_next(&iter)) {
    int idx_ref[] = {iter.idx[1]};
    long loc = gkyl_range_idx(&local_vel, iter.idx);
    long loc_ref = gkyl_range_idx(&local_ref, idx_ref);
    double *vm = gkyl_array_fetch(vmap_ho, loc);
    double *vm_ref = gkyl_array_fetch(vmap_ref, loc_ref);
    for (int i=0; i<basis_ref.num_basis; i++)
      TEST_CHECK( gkyl_compare(vm[basis_ref.num_basis+i], vm_ref[i], 1e-12) );
  }
  gkyl_array_release(vmap_ref);

  // Check min/max cell lengths over the whole grid.
  double dv_min[vdim], dv_max[vdim];
  gkyl_velocity_map_reduce_dv(gvm, GKYL_MIN, dv_min);
  gkyl_velocity_map_reduce_dv(gvm, GKYL_MAX, dv_max);

  double dv_comp[vdim];
  for (int d=0; d<vdim; d++) {
    dv_comp[d] = (upper_vel[d] - lower_vel[d])/cells_vel[d];
  }
  for (int d=0; d<vdim; d++) {
    double vc_lo[1], vc_up[1];
    double vp_lo[1], vp_up[1];
  
    // Check min dv.
    vc_lo[0] = 0.0;
    vc_up[0] = dv_comp[d];
    if (d==0) {
      test_vmap_1x2v_p1_mapc2p_vel_vpar(0.0, vc_lo, vp_lo, NULL);
      test_vmap_1x2v_p1_mapc2p_vel_vpar(0.0, vc_up, vp_up, NULL);
    }
    else {
      test_vmap_1x2v_p1_mapc2p_vel_mu(0.0, vc_lo, vp_lo, NULL);
      test_vmap_1x2v_p1_mapc2p_vel_mu(0.0, vc_up, vp_up, NULL);
    }
    double dv_min_ref = vp_up[0] - vp_lo[0];
    TEST_CHECK( gkyl_compare(dv_min[d], dv_min_ref, 1e-12) );

    // Check max dv.
    vc_lo[0] = upper_vel[d]-dv_comp[d];
    vc_up[0] = upper_vel[d];
    if (d==0) {
      test_vmap_1x2v_p1_mapc2p_vel_vpar(0.0, vc_lo, vp_lo, NULL);
      test_vmap_1x2v_p1_mapc2p_vel_vpar(0.0, vc_up, vp_up, NULL);
    }
    else {
      test_vmap_1x2v_p1_mapc2p_vel_mu(0.0, vc_lo, vp_lo, NULL);
      test_vmap_1x2v_p1_mapc2p_vel_mu(0.0, vc_up, vp_up, NULL);
    }
    double dv_max_ref = vp_up[0] - vp_lo[0];
    TEST_CHECK( gkyl_compare(dv_max[1], dv_max_ref, 1e-12) );
  }

  gkyl_array_release(vmap_ho);
  gkyl_velocity_map_release(gvm);
}

void
test_vmap_1x2v_p1_ho()
{
  test_vmap_1x2v_p1(false);
}

#ifdef GKYL_HAVE_CUDA
void
test_vmap_1x2v_p1_dev()
{
  test_vmap_1x2v_p1(true);
}
#endif

TEST_LIST = {
  { "test_vmap_1x2v_p1_ho", test_vmap_1x2v_p1_ho},

#ifdef GKYL_HAVE_CUDA
  { "test_vmap_1x2v_p1_dev", test_vmap_1x2v_p1_dev},
#endif
  { NULL, NULL },
};
