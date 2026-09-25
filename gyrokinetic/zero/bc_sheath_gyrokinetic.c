#include "gkyl_bc_sheath_gyrokinetic.h"
#include "gkyl_bc_sheath_gyrokinetic_priv.h"
#include <gkyl_alloc.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array_rio.h>
#include <assert.h>
#include <gkyl_gk_bc_type.h>

static struct gkyl_array *
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array *a;
  if (on_gpu) {
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  } else {
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  }
  return a;
}

static void
bc_sheath_gyrokinetic_diag_init(
  struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_range *phase_local
)
{
  // Decomposition is only along the parallel direction, so the (perpendicular
  // conf-space + mu) vcutsq array is owned entirely by the single rank abutting
  // this edge's boundary. That rank (non-empty skin range) writes it out
  // serially; all others do nothing.
  up->vcutsq_write = up->skin_r->volume > 0;
  if (!phase_grid || !phase_local || !up->vcutsq_write) {
    return;
  }

  int cdim = up->cdim;
  int pdim = phase_grid->ndim;
  int vdim = pdim - cdim;
  int vc_ndim = up->vcutsq_dim; // (cdim-1) + (vdim-1).

  // vcutsq is only meaningful (and the sheath kernels only exist) for vdim>1.
  if (vdim < 2) {
    return;
  }

  // Build the output grid spanning the perpendicular conf-space directions and
  // mu, mapping the (already SOL-aware) vcutsq index range to physical coords.
  double vc_lo[GKYL_MAX_DIM], vc_up[GKYL_MAX_DIM];
  int vc_cells[GKYL_MAX_DIM];

  int c = 0;
  for (int d = 0; d < cdim - 1; d++) { // perpendicular conf-space directions.
    int idx_lo = up->vcutsq_local.lower[c], idx_up = up->vcutsq_local.upper[c];
    vc_lo[c] = phase_grid->lower[d] + (idx_lo - phase_local->lower[d]) * phase_grid->dx[d];
    vc_up[c] = phase_grid->lower[d] + (idx_up - phase_local->lower[d] + 1) * phase_grid->dx[d];
    vc_cells[c] = idx_up - idx_lo + 1;
    c++;
  }
  for (int d = cdim + 1; d < pdim; d++) { // velocity directions other than vpar (mu, ...).
    vc_lo[c] = phase_grid->lower[d];
    vc_up[c] = phase_grid->upper[d];
    vc_cells[c] = phase_grid->cells[d];
    c++;
  }
  gkyl_rect_grid_init(&up->vcutsq_grid, vc_ndim, vc_lo, vc_up, vc_cells);

  // Range matching the output grid (1-based), used to write the local array.
  // Its shape matches up->vcutsq_local, so the flat data layout is identical.
  int nghost[GKYL_MAX_DIM] = {0};
  gkyl_create_grid_ranges(
    &up->vcutsq_grid, nghost, &up->vcutsq_diag_range_ext, &up->vcutsq_diag_range
  );

  // Host output buffer. On CPU we alias the working array to avoid a copy.
  up->vcutsq_ho =
    up->use_gpu ? gkyl_array_new(GKYL_DOUBLE, up->vcutsq_basis.num_basis, up->vcutsq_local.volume) :
                  gkyl_array_acquire(up->vcutsq);

  up->vcutsq_diag_on = true;
}

void
bc_gksheath_update_vcutsq_surrogate(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
  const struct gkyl_range *conf_r
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    bc_gksheath_update_vcutsq_surrogate_cu(
      up, phi, phi_wall, dens, temp, bmag, bimpact_angle, conf_r
    );
    return;
  }
#endif

  int cidx[GKYL_MAX_CDIM]; // Configuration space index.
  int csurfidx[GKYL_MAX_CDIM]; // Configuration space surface index.
  int vidx[GKYL_MAX_VDIM]; // Velocity space index.

  // Set the fixed parallel config index and vpar lower index.
  cidx[up->cdim - 1] = up->edge == GKYL_LOWER_EDGE ? up->skin_r->lower[up->cdim - 1] :
                                                     up->skin_r->upper[up->cdim - 1];
  csurfidx[up->cdim - 1] = up->edge == GKYL_LOWER_EDGE ? cidx[up->cdim - 1] :
                                                         cidx[up->cdim - 1] + 1;
  vidx[0] = up->skin_r->lower[up->cdim];

  // Loop over perpendicular config cells only to fill the KANN input vectors.
  // build_input writes n_nodes_per_cell nodes into nn_inp_out using stride dim_in:
  //   feature f of local node n -> nn_inp_out[f + dim_in * n]
  // This matches gkyl_kn_vec layout: node g data at data[g*N], feature f at data[g*N + f].
  struct gkyl_range_iter iter_xy;
  gkyl_range_iter_init(&iter_xy, &up->perp_local);
  while (gkyl_range_iter_next(&iter_xy)) {
    for (int d = 0; d < up->cdim - 1; d++) {
      cidx[d] = iter_xy.idx[d];
      csurfidx[d] = iter_xy.idx[d];
    }
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    long conf_surf_loc = gkyl_range_idx(conf_r, csurfidx);
    const double *phi_p = (const double *)gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double *)gkyl_array_cfetch(phi_wall, conf_loc);
    const double *dens_p = (const double *)gkyl_array_cfetch(dens, conf_loc);
    const double *temp_p = (const double *)gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double *)gkyl_array_cfetch(bmag, conf_loc);
    const double *bimpact_angle_p = (const double *)gkyl_array_cfetch(bimpact_angle, conf_surf_loc);

    // Contiguous block stride for the input parameter of the NN.
    long perp_k = gkyl_range_idx(&up->perp_local, iter_xy.idx);
    int stride = perp_k * up->perp_node_per_cell * up->kann_inp->N;

    up->kernels->vcutsq_input(
      phi_p, phi_wall_p, dens_p, temp_p, bmag_p, bimpact_angle_p, up->kann_inp->N,
      up->kann_inp->data + stride
    );
  }

  // Infer NN.
  gkyl_kann_net_apply(up->kann_net, up->kann_inp, up->kann_out);

  // Loop over perp config + mu cells to project NN output onto DG vcut_fact.
  // vcut_calc reads local node n's output block as nn_out + n*dim_out, matching
  // gkyl_kn_vec layout where node g data starts at data[g*N].
  struct gkyl_range_iter iter_xymu;
  gkyl_range_iter_init(&iter_xymu, &up->vcutsq_local);
  while (gkyl_range_iter_next(&iter_xymu)) {
    for (int d = 0; d < up->cdim - 1; d++) {
      cidx[d] = iter_xymu.idx[d];
    }
    vidx[1] = iter_xymu.idx[up->cdim - 1];

    long cperp_mu_loc = gkyl_range_idx(&up->vcutsq_local, iter_xymu.idx);
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    long vel_loc = gkyl_range_idx(&up->vel_map->local_vel, vidx);
    // 0-based sequential perp cell index from the perp part of the index.
    long perp_k = gkyl_range_idx(&up->perp_local, cidx);

    const double *vmap_p = (const double *)gkyl_array_cfetch(up->vel_map->vmap, vel_loc);
    const double *temp_p = (const double *)gkyl_array_cfetch(temp, conf_loc);
    const double *bmag_p = (const double *)gkyl_array_cfetch(bmag, conf_loc);
    double *vcutsq_p = (double *)gkyl_array_cfetch(up->vcutsq, cperp_mu_loc);

    int stride = perp_k * up->perp_node_per_cell * up->kann_out->N;
    // Compute the DG representation of the interpolated vcutsq values.
    up->kernels->vcutsq_surr(
      vmap_p, up->kann_out->data + stride, up->kann_out->N, temp_p, bmag_p, vcutsq_p
    );
  }
}

void
bc_gksheath_update_vcutsq_const(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
  const struct gkyl_range *conf_r
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    bc_gksheath_update_vcutsq_const_cu(up, phi, phi_wall, dens, temp, bmag, bimpact_angle, conf_r);
    return;
  }
#endif

  int cidx[GKYL_MAX_CDIM]; // Configuration space index.

  // Set the fixed parallel config index and vpar lower index.
  cidx[up->cdim - 1] = up->edge == GKYL_LOWER_EDGE ? up->skin_r->lower[up->cdim - 1] :
                                                     up->skin_r->upper[up->cdim - 1];

  // Loop over perpendicular config cells to compute vcutsq using the conducting sheath model.
  struct gkyl_range_iter iter_perpmu;
  gkyl_range_iter_init(&iter_perpmu, &up->vcutsq_local);
  while (gkyl_range_iter_next(&iter_perpmu)) {
    long perpmu_loc = gkyl_range_idx(&up->vcutsq_local, iter_perpmu.idx);
    double *vcutsq_p = (double *)gkyl_array_cfetch(up->vcutsq, perpmu_loc);

    for (int d = 0; d < up->cdim - 1; d++) {
      cidx[d] = iter_perpmu.idx[d];
    }
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    const double *phi_p = (const double *)gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double *)gkyl_array_cfetch(phi_wall, conf_loc);

    up->kernels->vcutsq_const(phi_p, phi_wall_p, up->q2Dm, vcutsq_p);
  }
}
struct gkyl_bc_sheath_gyrokinetic *
gkyl_bc_sheath_gyrokinetic_new(
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, const struct gkyl_velocity_map *vel_map, int cdim, double q2Dm,
  enum gkyl_gyrokinetic_bc_type type, const char *surrogate_model_path,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *phase_local,
  unsigned int poly_order, bool use_gpu
)
{
  // Allocate space for new updater.
  struct gkyl_bc_sheath_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->q2Dm = q2Dm;
  up->basis = basis;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;
  up->vel_map = gkyl_velocity_map_acquire(vel_map);
  int vdim = skin_r->ndim - cdim;
  up->type = type;
  up->update_vcutsq = bc_gksheath_update_vcutsq_const;
  up->vcutsq_dim = cdim - 1 + vdim - 1;

  // vcutsq diagnostic output is set up below (bc_sheath_gyrokinetic_diag_init);
  // these defaults keep release safe when the diagnostic is disabled.
  up->vcutsq_diag_on = false;
  up->vcutsq_write = false;
  up->vcutsq_ho = NULL;

  gkyl_cart_modal_serendip(&up->vcutsq_basis, up->vcutsq_dim, poly_order);

  int lower[up->vcutsq_dim], upper[up->vcutsq_dim];
  for (int d = 0; d < cdim - 1; d++) {
    lower[d] = skin_r->lower[d];
    upper[d] = skin_r->upper[d];
  }
  lower[up->vcutsq_dim - 1] = skin_r->lower[skin_r->ndim - 1];
  upper[up->vcutsq_dim - 1] = skin_r->upper[skin_r->ndim - 1];
  gkyl_range_init(&up->vcutsq_local, up->vcutsq_basis.ndim, lower, upper);

  int perp_lower[GKYL_MAX_CDIM], perp_upper[GKYL_MAX_CDIM];
  for (int d = 0; d < up->cdim - 1; d++) {
    perp_lower[d] = up->skin_r->lower[d];
    perp_upper[d] = up->skin_r->upper[d];
  }
  gkyl_range_init(&up->perp_local, up->cdim - 1, perp_lower, perp_upper);
  up->vcutsq = mkarr(use_gpu, up->vcutsq_basis.num_basis, up->vcutsq_local.volume);

  if (type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
    assert(vdim > 1);
    assert(
      q2Dm == -2 * GKYL_ELEMENTARY_CHARGE / GKYL_ELECTRON_MASS
    ); // The surrogate is only valid for electrons.
    // Assert existence of surrogate model path.
    FILE *fp = fopen(surrogate_model_path, "r");
    if (fp == NULL) {
      gkyl_exit("surrogate_model_path does not point to a valid file");
    }
    up->perp_node_per_cell = 1 << (up->cdim - 1);

    // gkyl KANN wrapper
    up->kann_net = gkyl_kann_net_load(surrogate_model_path, up->use_gpu);
    int dim_in = gkyl_kann_net_dim_in(up->kann_net);
    int dim_out = gkyl_kann_net_dim_out(up->kann_net);
    // The surrogate must be grid-augmented so dim_out is even.
    assert(dim_in == 3 && dim_out > 0 && dim_out % 2 == 0);
    int nperp_nodes = up->perp_node_per_cell * up->perp_local.volume;
    up->kann_inp = up->use_gpu ? gkyl_kn_vec_cu_dev_new(nperp_nodes, dim_in) :
                                 gkyl_kn_vec_new(nperp_nodes, dim_in);
    up->kann_out = up->use_gpu ? gkyl_kn_vec_cu_dev_new(nperp_nodes, dim_out) :
                                 gkyl_kn_vec_new(nperp_nodes, dim_out);
    up->kann_infer_xy_out = mkarr(up->use_gpu, dim_out * pow(2, up->cdim - 1), up->skin_r->volume);
    up->nn_out = up->use_gpu ? gkyl_cu_malloc(dim_out * pow(2, up->cdim - 1) * sizeof(double)) :
                               gkyl_malloc(dim_out * pow(2, up->cdim - 1) * sizeof(double));

    up->update_vcutsq = bc_gksheath_update_vcutsq_surrogate;
  }

  // Choose the kernels that does the reflection/no reflection/partial reflection,
  // and surrogate kernels if enabled.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
    gkyl_bc_gksheath_choose_reflectedf_kernel_cu(basis, edge, up->kernels_cu);
    if (type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      gkyl_bc_gksheath_choose_surrogate_kernels_cu(basis, edge, up->kernels_cu);
    } else {
      gkyl_bc_gksheath_choose_const_kernels_cu(basis, edge, up->kernels_cu);
    }
  } else {
    up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
    assert(up->kernels->reflectedf);
    if (type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      bc_gksheath_choose_vcutsq_surr_kernels(basis, edge, up->kernels);
    } else {
      bc_gksheath_choose_vcutsq_const_kernels(basis, edge, up->kernels);
    }
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
  assert(up->kernels->reflectedf);
  if (type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
    bc_gksheath_choose_vcutsq_surr_kernels(basis, edge, up->kernels);
  } else {
    bc_gksheath_choose_vcutsq_const_kernels(basis, edge, up->kernels);
  }
  up->kernels_cu = up->kernels;
#endif

  // Set up diagnostic output of the sheath velocity cutoff (vcutsq).
  bc_sheath_gyrokinetic_diag_init(up, phase_grid, phase_local);

  return up;
}

/* Modeled after gkyl_array_flip_copy_to_buffer_fn */
void
gkyl_bc_sheath_gyrokinetic_advance(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *density,
  const struct gkyl_array *temperature, const struct gkyl_array *bmag,
  const struct gkyl_array *bimpact_angle, struct gkyl_array *distf, const struct gkyl_range *conf_r
)
{
  // Update the sheath velocity cutoff.
  up->update_vcutsq(up, phi, phi_wall, density, temperature, bmag, bimpact_angle, conf_r);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_sheath_gyrokinetic_advance_cu(up, distf, conf_r);
    return;
  }
#endif

  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int vidx[2];
  int vcutsq_idx[up->vcutsq_local.ndim];

  int pdim = up->skin_r->ndim;
  int vpar_dir = up->cdim;
  int uplo = up->skin_r->upper[vpar_dir] + up->skin_r->lower[vpar_dir];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->skin_r);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, fidx);
    fidx[vpar_dir] = uplo - iter.idx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[up->dir] = up->ghost_r->lower[up->dir];

    long skin_loc = gkyl_range_idx(up->skin_r, iter.idx);
    long ghost_loc = gkyl_range_idx(up->ghost_r, fidx);

    const double *inp = (const double *)gkyl_array_cfetch(distf, skin_loc);
    double *out = (double *)gkyl_array_fetch(distf, ghost_loc);

    for (int d = up->cdim; d < pdim; d++) {
      vidx[d - up->cdim] = iter.idx[d];
    }
    long conf_loc = gkyl_range_idx(conf_r, iter.idx);
    long vel_loc = gkyl_range_idx(&up->vel_map->local_vel, vidx);

    const double *vmap_p = (const double *)gkyl_array_cfetch(up->vel_map->vmap, vel_loc);

    int vcutsq_dim = up->vcutsq_local.ndim;
    for (int d = 0; d < vcutsq_dim - 1; d++) {
      vcutsq_idx[d] = iter.idx[d]; // config space perp directions.
    }
    vcutsq_idx[vcutsq_dim - 1] = iter.idx[pdim - 1]; // mu direction.
    long vcutsq_loc = gkyl_range_idx(&up->vcutsq_local, vcutsq_idx);
    const double *vcutsq_p = (const double *)gkyl_array_cfetch(up->vcutsq, vcutsq_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->kernels->reflectedf(vmap_p, vcutsq_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void
gkyl_bc_sheath_gyrokinetic_write_vcutsq(
  struct gkyl_bc_sheath_gyrokinetic *up, struct gkyl_msgpack_data *meta, const char *fname
)
{
  if (!up->vcutsq_diag_on || !up->vcutsq_write) {
    return;
  }

  // Copy from device to the host buffer if needed (on CPU vcutsq_ho aliases vcutsq).
  if (up->use_gpu) {
    gkyl_array_copy(up->vcutsq_ho, up->vcutsq);
  }

  // Serial write: this rank owns the entire (perp conf-space + mu) array.
  gkyl_grid_sub_array_write(&up->vcutsq_grid, &up->vcutsq_diag_range, meta, up->vcutsq_ho, fname);
}

void
gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
  }
#endif
  gkyl_free(up->kernels);
  gkyl_velocity_map_release(up->vel_map);
  gkyl_array_release(up->vcutsq);

  if (up->vcutsq_diag_on && up->vcutsq_ho) {
    gkyl_array_release(up->vcutsq_ho);
  }

  if (up->type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
    gkyl_kann_net_release(up->kann_net);
    gkyl_kn_vec_release(up->kann_inp);
    gkyl_kn_vec_release(up->kann_out);
    gkyl_array_release(up->kann_infer_xy_out);
    if (up->use_gpu) {
      gkyl_cu_free(up->nn_out);
    } else {
      gkyl_free(up->nn_out);
    }
  }
  gkyl_free(up);
}