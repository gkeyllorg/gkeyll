#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_velocity_map.h>
#include <gkyl_velocity_map_priv.h>

#include <float.h>
#include <assert.h>

// Context for comp. coords = phys. coords (default).
struct mapc2p_vel_identity_ctx {
  int vdim;
};

// Comp. coords = phys. coords mapping (default).
static inline void
mapc2p_vel_identity(double t, const double *zc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct mapc2p_vel_identity_ctx *identity_ctx = ctx;
  int vdim = identity_ctx->vdim;
  for (int d=0; d<vdim; d++) vp[d] = zc[d];
}

void
gkyl_velocity_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_velocity_map *gvm = container_of(ref, struct gkyl_velocity_map, ref_count);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_velocity_map_is_cu_dev(gvm))
    gkyl_cart_modal_basis_release_cu(gvm->vmap_basis);
  else
    gkyl_cart_modal_basis_release(gvm->vmap_basis);
#else
  gkyl_cart_modal_basis_release(gvm->vmap_basis);
#endif
  gkyl_array_release(gvm->vmap);
  gkyl_array_release(gvm->vmap_sq);
  gkyl_array_release(gvm->vmap_prime);
  gkyl_array_release(gvm->jacobvel);
  gkyl_array_release(gvm->vmap_ho);

  if (gkyl_velocity_map_is_cu_dev(gvm))
    gkyl_cu_free(gvm->on_dev);

  gkyl_free(gvm);
}

struct gkyl_velocity_map*
gkyl_velocity_map_new(struct gkyl_mapc2p_inp mapc2p_in,
  struct gkyl_rect_grid grid, struct gkyl_rect_grid grid_vel,
  struct gkyl_range local, struct gkyl_range local_ext,
  struct gkyl_range local_vel, struct gkyl_range local_ext_vel, bool use_gpu)
{
  struct gkyl_velocity_map *gvm = gkyl_malloc(sizeof(*gvm));

  gvm->is_identity = mapc2p_in.mapping == 0;
  gvm->grid = grid;
  gvm->grid_vel = grid_vel;
  gvm->local = local;
  gvm->local_ext = local_ext;
  gvm->local_vel = local_vel;
  gvm->local_ext_vel = local_ext_vel;

  // Project velocity mappings, compute their derivative
  // and the corresponding Jacobian.
  int vmap_poly_order = 1; // Polynomial order used to discretize the map.
  int vmap_sq_poly_order = 2; // Poly order the squared map is represented with.

  int pdim = local.ndim, vdim = local_vel.ndim;
  int cdim = pdim - vdim;

  // Create the velocity space bases.
  struct gkyl_basis vmap_basis_vdim, vmap_sq_basis;
  gvm->vmap_basis = gkyl_cart_modal_serendip_new(1, vmap_poly_order);
  gkyl_cart_modal_serendip(&vmap_basis_vdim, vdim, vmap_poly_order);
  gkyl_cart_modal_serendip(&vmap_sq_basis, 1, vmap_sq_poly_order);
  // Need a host copy of vmap_basis for some IC projection options.
  gkyl_cart_modal_serendip(&gvm->vmap_basis_ho, 1, vmap_poly_order);

  gvm->vmap       = mkarr(false, vdim*gvm->vmap_basis->num_basis, gvm->local_ext_vel.volume);
  gvm->vmap_sq    = mkarr(false, vdim*vmap_sq_basis.num_basis, gvm->local_ext_vel.volume);
  gvm->vmap_prime = mkarr(false, vdim, gvm->local_ext_vel.volume);
  gvm->jacobvel   = mkarr(false, 1, gvm->local_ext.volume);
  // Need a host copy of vmap for some IC setting options.
  gvm->vmap_ho    = gkyl_array_acquire(gvm->vmap);

  // Project the velocity mapping (onto a vdim basis).
  struct gkyl_array *vmap_vdim = mkarr(false, vdim*vmap_basis_vdim.num_basis, gvm->vmap->size);
  struct gkyl_array *vmapc1 = mkarr(false, 1, gvm->vmap->size);

  gkyl_eval_on_nodes *evup;
  if (gvm->is_identity) {
    struct mapc2p_vel_identity_ctx identity_ctx = { .vdim = vdim };
    evup = gkyl_eval_on_nodes_new(&gvm->grid_vel, &vmap_basis_vdim,
      vdim, mapc2p_vel_identity, &identity_ctx);
  }
  else {
    evup = gkyl_eval_on_nodes_new(&gvm->grid_vel, &vmap_basis_vdim,
      vdim, mapc2p_in.mapping, mapc2p_in.ctx);
  }
  gkyl_eval_on_nodes_advance(evup, 0., &gvm->local_vel, vmap_vdim);

  // Extract the 1D basis expansions from the vdim basis expansion.
  for (int d=0; d<vdim; d++) {
    int numb_1d = gvm->vmap_basis->num_basis;
    int numb_vdim = vmap_basis_vdim.num_basis;
    gkyl_array_set_offset(vmapc1, 1./sqrt(vdim), vmap_vdim, d*numb_vdim);
    gkyl_array_set_offset(gvm->vmap, 1., vmapc1, d*numb_1d);
    gkyl_array_set_offset(vmapc1, 1./sqrt(vdim), vmap_vdim, d*numb_vdim+d+1);
    gkyl_array_set_offset(gvm->vmap, 1., vmapc1, d*numb_1d+1);
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(vmapc1);
  gkyl_array_release(vmap_vdim);

  struct gkyl_array *vmap1d       = mkarr(false, gvm->vmap_basis->num_basis, gvm->vmap->size);
  struct gkyl_array *vmap1d_p2    = mkarr(false, vmap_sq_basis.num_basis, gvm->vmap->size);
  struct gkyl_array *vmap_prime1d = mkarr(false, 1, gvm->vmap_prime->size);
  gkyl_array_clear(vmap1d_p2, 0.);
  for (int d=0; d<vdim; d++) {
    gkyl_array_set_offset(vmap1d, 1., gvm->vmap, d*gvm->vmap_basis->num_basis);

    // Compute the square mapping via weak multiplication.
    gkyl_array_set_offset(vmap1d_p2, 1., vmap1d, 0);
    gkyl_dg_mul_op(&vmap_sq_basis, d, gvm->vmap_sq, 0, vmap1d_p2, 0, vmap1d_p2);

    // Compute the derivative of the mapping.
    gkyl_array_set_offset(vmap_prime1d, sqrt(6.)/gvm->grid_vel.dx[d], vmap1d, 1);
    gkyl_array_set_offset(gvm->vmap_prime, 1., vmap_prime1d, d);
  }
  gkyl_array_release(vmap1d);
  gkyl_array_release(vmap_prime1d);
  gkyl_array_release(vmap1d_p2);

  // Compute the velocity space Jacobian (cell-wise constant).
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &gvm->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&gvm->local, iter.idx);
    double *jacv_d = gkyl_array_fetch(gvm->jacobvel, linidx);
    jacv_d[0] = 1.;

    int vidx[vdim];
    for (int d=0; d<vdim; d++) vidx[d] = iter.idx[cdim+d];
    long vlinidx = gkyl_range_idx(&gvm->local_vel, vidx);
    double *vprime_d = gkyl_array_fetch(gvm->vmap_prime, vlinidx);
    for (int d=0; d<vdim; d++)
      jacv_d[0] *= fabs(vprime_d[d]);
  }

  // Save the velocity at the boundaries.
  gkyl_velocity_map_get_boundary_values(gvm, gvm->vbounds);

  gvm->flags = 0;
  GKYL_CLEAR_CU_ALLOC(gvm->flags);
  gvm->ref_count = gkyl_ref_count_init(gkyl_velocity_map_free);
  gvm->on_dev = gvm; // CPU eqn gvm points to itself

  struct gkyl_velocity_map *gvm_out = gvm;
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gvm_out = gkyl_velocity_map_new_cu_dev(gvm);
    gkyl_velocity_map_release(gvm);
  }
#endif

  return gvm_out;
}

bool
gkyl_velocity_map_is_cu_dev(const struct gkyl_velocity_map* gvm)
{
  return GKYL_IS_CU_ALLOC(gvm->flags);
}

void
gkyl_velocity_map_write(const struct gkyl_velocity_map* gvm, struct gkyl_comm* species_comm,
  const char* app_name, const char* species_name)
{
  // Write out the velocity space mapping.
  struct gkyl_array *vmap_ho = gvm->vmap_ho, *jacobvel_ho = gvm->jacobvel;
  if (gkyl_velocity_map_is_cu_dev(gvm)) {
    jacobvel_ho = mkarr(false, gvm->jacobvel->ncomp, gvm->jacobvel->size);
    gkyl_array_copy(jacobvel_ho, gvm->jacobvel);
  }

  int rank, sz;
  int err = gkyl_comm_get_rank(species_comm, &rank);
  if (rank == 0) {
    // Write out the velocity mapping.
    const char *fmt0 = "%s-%s_mapc2p_vel.gkyl";
    sz = gkyl_calc_strlen(fmt0, app_name, species_name);
    char fileNm0[sz+1]; // ensures no buffer overflow
    snprintf(fileNm0, sizeof fileNm0, fmt0, app_name, species_name);
    gkyl_grid_sub_array_write(&gvm->grid_vel, &gvm->local_vel, NULL, vmap_ho, fileNm0);
  }

  // Write out the velocity space Jacobian.
  const char *fmt1 = "%s-%s_jacobvel.gkyl";
  sz = gkyl_calc_strlen(fmt1, app_name, species_name);
  char fileNm1[sz+1]; // ensures no buffer overflow
  snprintf(fileNm1, sizeof fileNm1, fmt1, app_name, species_name);
  gkyl_comm_array_write(species_comm, &gvm->grid, &gvm->local, NULL, jacobvel_ho, fileNm1);

  if (gkyl_velocity_map_is_cu_dev(gvm)) {
    gkyl_array_release(jacobvel_ho);
  }
}

void
gkyl_velocity_map_get_boundary_values(const struct gkyl_velocity_map* gvm, double *vbounds)
{
  int vdim = gvm->local_vel.ndim;
  for (int d=0; d<vdim; ++d) {
    long vlinidx;
    double vlog[1], *vmap_d;
    int off = d*gvm->vmap_basis_ho.num_basis;

    vlog[0] = -1.0;
    vlinidx = gkyl_range_idx(&gvm->local_vel, gvm->local_vel.lower);
    vmap_d = gkyl_array_fetch(gvm->vmap_ho, vlinidx);
    vbounds[d] = gvm->vmap_basis_ho.eval_expand(vlog, off+vmap_d);

    vlog[0] = 1.0;
    vlinidx = gkyl_range_idx(&gvm->local_vel, gvm->local_vel.upper);
    vmap_d = gkyl_array_fetch(gvm->vmap_ho, vlinidx);
    vbounds[d + vdim] = gvm->vmap_basis_ho.eval_expand(vlog, off+vmap_d);
  }
}

void
gkyl_velocity_map_reduce_dv_range(const struct gkyl_velocity_map* gvm, enum gkyl_array_op op,
  double *dv_m, struct gkyl_range range_vel)
{
  int vdim = gvm->local_vel.ndim;
  if (op == GKYL_MIN)
    for (int d=0; d<vdim; ++d) dv_m[d] = DBL_MAX;
  else if (op == GKYL_MAX)
    for (int d=0; d<vdim; ++d) dv_m[d] = -DBL_MAX;
  else
    assert(false);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range_vel);
  while (gkyl_range_iter_next(&iter)) {

    long linidx = gkyl_range_idx(&range_vel, iter.idx);
    double *vmap_d = gkyl_array_fetch(gvm->vmap_ho, linidx);

    for (int d=0; d<vdim; ++d) {
      int off = d*gvm->vmap_basis_ho.num_basis;
      double vlog[1];

      vlog[0] = -1.0;
      double vlo = gvm->vmap_basis_ho.eval_expand(vlog, off+vmap_d);

      vlog[0] = 1.0;
      double vup = gvm->vmap_basis_ho.eval_expand(vlog, off+vmap_d);

      double dv = vup - vlo;

      if (op == GKYL_MIN)
        dv_m[d] = GKYL_MIN2(dv_m[d], dv);
      else if (op == GKYL_MAX)
        dv_m[d] = GKYL_MAX2(dv_m[d], dv);
    }
  }
}

void
gkyl_velocity_map_reduce_dv(const struct gkyl_velocity_map* gvm, enum gkyl_array_op op, double *dv_m)
{
  gkyl_velocity_map_reduce_dv_range(gvm, op, dv_m, gvm->local_vel);
}

void
gkyl_velocity_map_eval_c2p(const struct gkyl_velocity_map* gvm, const double *zc, double *vp)
{
  // Find the index of the cell containing zc.
  int idx_zc[GKYL_MAX_VDIM];
  for(int d = 0; d<gvm->local_ext_vel.ndim; d++){
    int idx = gvm->local_ext_vel.lower[d] + (int) floor((zc[d] - (gvm->grid_vel.lower[d]) )/gvm->grid_vel.dx[d]);
    // Bound idx to the range in the grid. If it falls outside of that is due
    // to floating point arithmetic, or due to an error in the code.
    idx = GKYL_MIN2(idx, gvm->local_ext_vel.upper[d]);
    idx = GKYL_MAX2(idx, gvm->local_ext_vel.lower[d]);
    idx_zc[d] = idx;
  }

  // Fetch DG coefficients of the velocity map in idx_gvm.
  long lidx_zc = gkyl_range_idx(&gvm->local_ext_vel, idx_zc);
  const double *vmap_c = gkyl_array_cfetch(gvm->vmap_ho, lidx_zc);

  double zc_cc[GKYL_MAX_VDIM];
  gkyl_rect_grid_cell_center(&gvm->grid_vel, idx_zc, zc_cc);

  for(int d=0; d<gvm->local_ext_vel.ndim; d++){
    // Convert computational to logical coord.
    double zlog[] = {(zc[d] - zc_cc[d]) / (0.5*gvm->grid_vel.dx[d])};
    // Evaluate vmap expansion at logical coord.
    vp[d] = gvm->vmap_basis_ho.eval_expand(zlog, &vmap_c[d*gvm->vmap_basis_ho.num_basis]);
  }
}

struct gkyl_velocity_map*
gkyl_velocity_map_acquire(const struct gkyl_velocity_map* gvm)
{
  gkyl_ref_count_inc(&gvm->ref_count);
  return (struct gkyl_velocity_map*) gvm;
}

void
gkyl_velocity_map_release(const struct gkyl_velocity_map *gvm)
{
  gkyl_ref_count_dec(&gvm->ref_count);
}
