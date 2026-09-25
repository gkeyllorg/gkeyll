#include <string.h>
#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_velocity_map_priv.h>

// Computational coordinates = physical coordinates mapping (default).
static void
vlasov_vmap_identity(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = vc[0];
}

static struct gkyl_array *
mkarr(bool on_gpu, long nc, long size)
{
  return on_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                  gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Populate the (host) velocity map arrays from the input 1D mappings by
// constructing a C^1 piecewise cubic representation in each velocity
// dimension and evaluating its derivative at the Gauss-Legendre quadrature
// points needed by volume, surface, and projection operations.
static void
vlasov_velocity_map_c1_cubic(
  const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange, vmap_cubic_t vmap_op,
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM], struct gkyl_array *vmap,
  struct gkyl_array *jacob_vel, struct gkyl_array *jacob_vel_surf, struct gkyl_array *vmap_pgkyl,
  struct gkyl_array *vmap_avg_pgkyl, struct gkyl_array *jacob_vel_gauss
)
{
  int vdim = vgrid->ndim;
  struct gkyl_array *v_nodal[3];
  struct gkyl_array *v_cubic[3];
  struct gkyl_dg_basis_op_mem *mem[3];

  // Make 1D cubic basis for constructing C^1 expansion
  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  // 1D ranges for indexing 1D cubic bases
  struct gkyl_range local[3], local_ext[3];

  // Loop over number of dimensions and construct 1D mappings
  for (int i = 0; i < vdim; ++i) {
    double lower[] = {vgrid->lower[i]}, upper[] = {vgrid->upper[i]};
    int cells[] = {vgrid->cells[i]};

    struct gkyl_rect_grid grid_1d;
    gkyl_rect_grid_init(&grid_1d, 1, lower, upper, cells);

    // nodal grid used for constructing physical coordinates
    double nc_lower[] = {lower[0] - 0.5 * grid_1d.dx[0]};
    double nc_upper[] = {upper[0] + 0.5 * grid_1d.dx[0]};
    int nc_cells[] = {cells[0] + 1};
    struct gkyl_rect_grid nc_grid;
    gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

    int nghost[GKYL_MAX_CDIM] = {0, 0};
    gkyl_create_grid_ranges(&grid_1d, nghost, &local_ext[i], &local[i]);

    struct gkyl_range nc_local, nc_local_ext;
    gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

    v_nodal[i] = gkyl_array_new(GKYL_DOUBLE, 1, cells[0] + 1);
    v_cubic[i] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext[i].volume);
    mem[i] = gkyl_dg_alloc_cubic_1d(cells[0]);
    double xn[1];

    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid_1d, iter.idx, xn);

      double *pn = gkyl_array_fetch(v_nodal[i], nidx);
      inp_vmap[i].eval_vmap(0.0, xn, pn, inp_vmap[i].ctx);
    }

    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem[i], cells[0], grid_1d.dx[0], v_nodal[i], v_cubic[i]);
  }

  // initialize the mapping
  const double *v_cubic_dir[3]; // 1D cubic in each direction
  int vidx_1D[1]; // 1D index for indexing correct cubic mapping

  struct gkyl_range_iter iter_vmap;
  gkyl_range_iter_init(&iter_vmap, vrange);
  while (gkyl_range_iter_next(&iter_vmap)) {
    long loc_vel = gkyl_range_idx(vrange, iter_vmap.idx);

    for (int i = 0; i < vdim; ++i) {
      vidx_1D[0] = iter_vmap.idx[i];
      long loc_vel_1D = gkyl_range_idx(&local[i], vidx_1D);
      v_cubic_dir[i] = gkyl_array_cfetch(v_cubic[i], loc_vel_1D);
    }
    double *vmap_d = gkyl_array_fetch(vmap, loc_vel);
    double *jacob_vel_d = gkyl_array_fetch(jacob_vel, loc_vel);
    double *jacob_vel_surf_d = gkyl_array_fetch(jacob_vel_surf, loc_vel);
    double *vmap_pgkyl_d = gkyl_array_fetch(vmap_pgkyl, loc_vel);
    double *vmap_avg_pgkyl_d = gkyl_array_fetch(vmap_avg_pgkyl, loc_vel);
    double *jacob_vel_gauss_d = gkyl_array_fetch(jacob_vel_gauss, loc_vel);

    vmap_op(
      vgrid->dx, v_cubic_dir, vmap_d, jacob_vel_d, jacob_vel_surf_d, vmap_pgkyl_d, vmap_avg_pgkyl_d,
      jacob_vel_gauss_d
    );
  }

  // free temporary memory
  for (int i = 0; i < vdim; ++i) {
    gkyl_array_release(v_nodal[i]);
    gkyl_array_release(v_cubic[i]);
    gkyl_dg_basis_op_mem_release(mem[i]);
  }
}

// Populate the (host) velocity map arrays from the input 1D mappings by
// constructing a continuous (C^0) piecewise linear representation in each
// velocity dimension: an eval_on_nodes at the Gauss-Lobatto nodes (the cell
// vertices for a p=1 basis, shared between neighboring cells so the
// representation is continuous by construction). The linear map is stored as
// a degenerate cubic in the same vdim*4 layout as the C^1 cubic
// representation (quadratic and cubic coefficients identically zero), so
// consumers that evaluate the map with a 1D cubic basis work unchanged. The
// per-direction Jacobian is piecewise constant (discontinuous at cell
// interfaces), so every quadrature-point value in jacob_vel, jacob_vel_surf,
// and jacob_vel_gauss within a cell is the same constant; the arrays keep
// the cubic layout with that constant replicated.
static void
vlasov_velocity_map_c0_linear(
  const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange, int v_poly_order,
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM], struct gkyl_array *vmap,
  struct gkyl_array *jacob_vel, struct gkyl_array *jacob_vel_surf, struct gkyl_array *vmap_pgkyl,
  struct gkyl_array *vmap_avg_pgkyl, struct gkyl_array *jacob_vel_gauss
)
{
  int vdim = vgrid->ndim;
  struct gkyl_array *v_linear[3];

  // Make 1D linear basis for constructing C^0 expansion
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // 1D ranges for indexing 1D linear bases
  struct gkyl_range local[3], local_ext[3];

  // Loop over number of dimensions and construct 1D mappings
  for (int i = 0; i < vdim; ++i) {
    double lower[] = {vgrid->lower[i]}, upper[] = {vgrid->upper[i]};
    int cells[] = {vgrid->cells[i]};

    struct gkyl_rect_grid grid_1d;
    gkyl_rect_grid_init(&grid_1d, 1, lower, upper, cells);

    int nghost[GKYL_MAX_CDIM] = {0, 0};
    gkyl_create_grid_ranges(&grid_1d, nghost, &local_ext[i], &local[i]);

    v_linear[i] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext[i].volume);

    gkyl_eval_on_nodes *ev_map =
      gkyl_eval_on_nodes_new(&grid_1d, &basis, 1, inp_vmap[i].eval_vmap, inp_vmap[i].ctx);
    gkyl_eval_on_nodes_advance(ev_map, 0.0, &local[i], v_linear[i]);
    gkyl_eval_on_nodes_release(ev_map);
  }

  // Location of the four 1D cubic coefficients of direction d inside the
  // full-dimensional p=3 serendipity vmap_pgkyl representation, and the scale
  // factor relating them to the 1D coefficients (same placement the C^1 cubic
  // kernels use; the quadratic and cubic slots just hold zeros here).
  static const int pgkyl_idx[3][3][4] = {
    {{0, 1, 2, 3}, {0}, {0}}, // 1V
    {{0, 1, 4, 8}, {0, 2, 5, 9}, {0}}, // 2V
    {{0, 1, 7, 17}, {0, 2, 8, 18}, {0, 3, 9, 19}} // 3V
  };
  static const double pgkyl_scale[3] = {1.0, 1.4142135623730951, 2.0};

  int num_quad_vol = v_poly_order + 1, num_quad_surf = v_poly_order + 2;
  long pgkyl_nb = vmap_pgkyl->ncomp / vdim;
  int vidx_1D[1]; // 1D index for indexing correct linear mapping

  struct gkyl_range_iter iter_vmap;
  gkyl_range_iter_init(&iter_vmap, vrange);
  while (gkyl_range_iter_next(&iter_vmap)) {
    long loc_vel = gkyl_range_idx(vrange, iter_vmap.idx);

    double *vmap_d = gkyl_array_fetch(vmap, loc_vel);
    double *jacob_vel_d = gkyl_array_fetch(jacob_vel, loc_vel);
    double *jacob_vel_surf_d = gkyl_array_fetch(jacob_vel_surf, loc_vel);
    double *vmap_pgkyl_d = gkyl_array_fetch(vmap_pgkyl, loc_vel);
    double *vmap_avg_pgkyl_d = gkyl_array_fetch(vmap_avg_pgkyl, loc_vel);
    double *jacob_vel_gauss_d = gkyl_array_fetch(jacob_vel_gauss, loc_vel);

    double jacob_tot = 1.0;
    for (int i = 0; i < vdim; ++i) {
      vidx_1D[0] = iter_vmap.idx[i];
      long loc_vel_1D = gkyl_range_idx(&local[i], vidx_1D);
      const double *lin = gkyl_array_cfetch(v_linear[i], loc_vel_1D);

      double *p = &vmap_d[4 * i];
      p[0] = lin[0];
      p[1] = lin[1];
      p[2] = 0.0;
      p[3] = 0.0;

      // Piecewise constant Jacobian: (2/dv)*d/dz of the linear map, with
      // sqrt(3/2) from the derivative of the 1D linear basis function.
      double jac = (2.0 / vgrid->dx[i]) * 1.224744871391589 * p[1];
      jacob_tot *= jac;
      for (int q = 0; q < num_quad_vol; ++q) {
        jacob_vel_d[i * num_quad_vol + q] = jac;
      }
      for (int q = 0; q < num_quad_surf; ++q) {
        jacob_vel_surf_d[i * num_quad_surf + q] = jac;
      }

      for (int k = 0; k < 4; ++k) {
        vmap_pgkyl_d[i * pgkyl_nb + pgkyl_idx[vdim - 1][i][k]] = pgkyl_scale[vdim - 1] * p[k];
      }
      vmap_avg_pgkyl_d[i] = p[0] / 1.4142135623730951;
    }
    for (int q = 0; q < jacob_vel_gauss->ncomp; ++q) {
      jacob_vel_gauss_d[q] = jacob_tot;
    }
  }

  // free temporary memory
  for (int i = 0; i < vdim; ++i) {
    gkyl_array_release(v_linear[i]);
  }
}

static void
gkyl_vlasov_velocity_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_velocity_map *vvm =
    container_of(ref, struct gkyl_vlasov_velocity_map, ref_count);

  gkyl_array_release(vvm->vmap);
  gkyl_array_release(vvm->jacob_vel);
  gkyl_array_release(vvm->jacob_vel_surf);
  gkyl_array_release(vvm->jacob_vel_gauss);

  gkyl_array_release(vvm->vmap_host);
  gkyl_array_release(vvm->jacob_vel_host);
  gkyl_array_release(vvm->jacob_vel_surf_host);
  gkyl_array_release(vvm->jacob_vel_gauss_host);

  gkyl_array_release(vvm->vmap_pgkyl_host);
  gkyl_array_release(vvm->vmap_avg_pgkyl_host);

  gkyl_free(vvm);
}

struct gkyl_vlasov_velocity_map *
gkyl_vlasov_velocity_map_new(
  const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange,
  const struct gkyl_basis *vel_basis, struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM],
  bool use_lo, bool use_gpu
)
{
  struct gkyl_vlasov_velocity_map *vvm = gkyl_malloc(sizeof(*vvm));

  int vdim = vgrid->ndim;
  int poly_order = vel_basis->poly_order;
  // The tensor p=2 velocity basis with use_lo is the low-order variant of the
  // tensor p=1 hybrid: its surface kernels use p+1 = 3 quadrature nodes, so
  // jacob_vel_surf is stored at 3 nodes per direction instead of p+2. Every
  // other basis has a single surface node count (p+2), use_lo included: the
  // hybrid is the only basis with distinct lo/ho surface variants.
  bool hyb_lo = (vel_basis->b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 2) && use_lo;

  vvm->grid_vel = *vgrid;
  vvm->local_vel = *vrange;
  vvm->basis_vel = *vel_basis;
  // Representation is chosen by the velocity basis: tensor bases use the C^1
  // cubic, Serendipity bases use the C^0 linear whose Jacobian is piecewise
  // constant. Tensor p=1 hybrid species pass their p=2 tensor velocity basis
  // here (the hybrid is p=2 in velocity space), so p=1 is always Serendipity
  // and always the C^0 linear.
  vvm->rep = (vel_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY || poly_order == 1) ?
               GKYL_VLASOV_VMAP_C0_LINEAR :
               GKYL_VLASOV_VMAP_C1_CUBIC;

  vvm->is_identity = true;
  for (int v = 0; v < vdim; ++v) {
    if (inp_vmap[v].eval_vmap) {
      vvm->is_identity = false;
    }
  }
  // The map is always constructed and consumed, for every basis and for
  // uniform (identity) grids alike; there is a single code path.
  vvm->is_mapped = true;

  // velocity map is always a C^1 cubic representation in each direction (up to 3V; 3*4=12 components)
  vvm->vmap = mkarr(use_gpu, vdim * 4, vrange->volume);
  // velocity-space Jacobian at quadrature points and "surface" quadrature points. Used to compute
  // 1/Jvi nodally in volume and surface operations respectively, with surface operations utilizing
  // more quadrature points to eliminate aliasing errors.
  vvm->jacob_vel = mkarr(use_gpu, vdim * (poly_order + 1), vrange->volume);
  vvm->jacob_vel_surf =
    mkarr(use_gpu, vdim * (hyb_lo ? poly_order + 1 : poly_order + 2), vrange->volume);
  // need special basis sets to get the correct number of coefficients in 2V and 3V for constructing
  // the mapping in post-processing, as well as storing the velocity-space Jacobian at quadrature points.
  struct gkyl_basis vmap_pgkyl_basis, jacob_vel_basis;
  gkyl_cart_modal_serendip(&vmap_pgkyl_basis, vdim, 3);
  gkyl_cart_modal_tensor(&jacob_vel_basis, vdim, poly_order);
  // Keep the I/O basis so the write can label the file with the basis the data
  // is actually stored in (the species' velocity basis does not describe these arrays).
  vvm->basis_pgkyl = vmap_pgkyl_basis;
  // velocity-space Jacobian at Gaussian quadrature points for projecting distribution functions
  // and dividing out velocity-space Jacobian nodally.
  vvm->jacob_vel_gauss = mkarr(use_gpu, jacob_vel_basis.num_basis, vrange->volume);

  // host-side mirrors for initialization and I/O
  if (use_gpu) {
    vvm->vmap_host = mkarr(false, vvm->vmap->ncomp, vvm->vmap->size);
    vvm->jacob_vel_host = mkarr(false, vvm->jacob_vel->ncomp, vvm->jacob_vel->size);
    vvm->jacob_vel_surf_host = mkarr(false, vvm->jacob_vel_surf->ncomp, vvm->jacob_vel_surf->size);
    vvm->jacob_vel_gauss_host =
      mkarr(false, vvm->jacob_vel_gauss->ncomp, vvm->jacob_vel_gauss->size);
  } else {
    vvm->vmap_host = gkyl_array_acquire(vvm->vmap);
    vvm->jacob_vel_host = gkyl_array_acquire(vvm->jacob_vel);
    vvm->jacob_vel_surf_host = gkyl_array_acquire(vvm->jacob_vel_surf);
    vvm->jacob_vel_gauss_host = gkyl_array_acquire(vvm->jacob_vel_gauss);
  }

  // velocity map for I/O
  vvm->vmap_pgkyl_host = mkarr(false, vdim * vmap_pgkyl_basis.num_basis, vrange->volume);
  vvm->vmap_avg_pgkyl_host = mkarr(false, vdim, vrange->volume);

  // Default unspecified directions to the identity map.
  struct gkyl_vlasov_velocity_map_inp inp[GKYL_MAX_CDIM];
  for (int v = 0; v < vdim; ++v) {
    if (inp_vmap[v].eval_vmap) {
      inp[v] = inp_vmap[v];
    } else {
      inp[v].eval_vmap = vlasov_vmap_identity;
      inp[v].ctx = 0;
    }
  }

  if (vvm->rep == GKYL_VLASOV_VMAP_C0_LINEAR) {
    vlasov_velocity_map_c0_linear(
      vgrid, vrange, poly_order, inp, vvm->vmap_host, vvm->jacob_vel_host, vvm->jacob_vel_surf_host,
      vvm->vmap_pgkyl_host, vvm->vmap_avg_pgkyl_host, vvm->jacob_vel_gauss_host
    );
  } else {
    // C^1 cubic kernels exist for tensor p=2,3; the hyb_lo variants differ
    // only in writing jacob_vel_surf at 3 surface nodes instead of 4.
    vmap_cubic_t vmap_op = hyb_lo ? choose_vmap_hyb_lo_kern(vdim) :
                                    choose_vmap_kern(vdim, poly_order);
    assert(vmap_op);
    vlasov_velocity_map_c1_cubic(
      vgrid, vrange, vmap_op, inp, vvm->vmap_host, vvm->jacob_vel_host, vvm->jacob_vel_surf_host,
      vvm->vmap_pgkyl_host, vvm->vmap_avg_pgkyl_host, vvm->jacob_vel_gauss_host
    );
  }

  if (use_gpu) {
    gkyl_array_copy(vvm->vmap, vvm->vmap_host);
    gkyl_array_copy(vvm->jacob_vel, vvm->jacob_vel_host);
    gkyl_array_copy(vvm->jacob_vel_surf, vvm->jacob_vel_surf_host);
    gkyl_array_copy(vvm->jacob_vel_gauss, vvm->jacob_vel_gauss_host);
  }

  vvm->flags = 0;
  if (use_gpu) {
    GKYL_SET_CU_ALLOC(vvm->flags);
  }
  vvm->ref_count = gkyl_ref_count_init(gkyl_vlasov_velocity_map_free);

  return vvm;
}

bool
gkyl_vlasov_velocity_map_is_cu_dev(const struct gkyl_vlasov_velocity_map *vvm)
{
  return GKYL_IS_CU_ALLOC(vvm->flags);
}

void
gkyl_vlasov_velocity_map_eval_c2p(
  const struct gkyl_vlasov_velocity_map *vvm, const double *vc, double *vp
)
{
  // Evaluate the stored (DG) velocity map at the computational velocity vc. The
  // map is stored per direction in the 4-slot (cubic) layout, so a 1D p=3 modal
  // basis evaluates each direction's expansion.
  int vdim = vvm->grid_vel.ndim;

  // Identity map: computational coordinates ARE physical coordinates. Return
  // them exactly (bitwise) instead of reconstructing v through the DG
  // expansion, so uniform grids are truly unaffected by the map machinery.
  if (vvm->is_identity) {
    for (int d = 0; d < vdim; ++d) {
      vp[d] = vc[d];
    }
    return;
  }

  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // Find the index of the cell containing vc.
  int idx_vc[GKYL_MAX_VDIM];
  for (int d = 0; d < vdim; ++d) {
    int idx =
      vvm->local_vel.lower[d] + (int)floor((vc[d] - vvm->grid_vel.lower[d]) / vvm->grid_vel.dx[d]);
    idx = GKYL_MIN2(idx, vvm->local_vel.upper[d]);
    idx = GKYL_MAX2(idx, vvm->local_vel.lower[d]);
    idx_vc[d] = idx;
  }

  // Fetch DG coefficients of the velocity map in idx_vc.
  long lidx_vc = gkyl_range_idx(&vvm->local_vel, idx_vc);
  const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, lidx_vc);

  double vc_cc[GKYL_MAX_VDIM];
  gkyl_rect_grid_cell_center(&vvm->grid_vel, idx_vc, vc_cc);
  for (int d = 0; d < vdim; ++d) {
    // Convert computational to logical coord, then evaluate the expansion.
    double vlog = (vc[d] - vc_cc[d]) / (0.5 * vvm->grid_vel.dx[d]);
    vp[d] = b1.eval_expand(&vlog, &vmap_c[d * 4]);
  }
}

void
gkyl_vlasov_velocity_map_write(
  const struct gkyl_vlasov_velocity_map *vvm, struct gkyl_comm *comm, const char *app_name,
  const char *species_name
)
{
  int rank;
  gkyl_comm_get_rank(comm, &rank);
  // Since velocity space has no decomposition currently, make sure to only
  // have rank 0 write out the velocity map.
  if (rank != 0) {
    return;
  }

  // The map is static and written once, so both files carry frame-0 metadata.
  struct gkyl_msgpack_data *mt_vmap = gkyl_msgpack_create(
    4,
    (struct gkyl_msgpack_map_elem[]){
      {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
      {.key = "frame", .elem_type = GKYL_MP_INT, .ival = 0},
      {.key = "polyOrder", .elem_type = GKYL_MP_INT, .ival = vvm->basis_pgkyl.poly_order},
      {.key = "basisType", .elem_type = GKYL_MP_STRING, .cval = vvm->basis_pgkyl.id},
    }
  );

  const char *fmt_vmap = "%s-%s_vmap.gkyl";
  int sz_vmap = gkyl_calc_strlen(fmt_vmap, app_name, species_name);
  char fileNm_vmap[sz_vmap + 1]; // ensures no buffer overflow
  snprintf(fileNm_vmap, sizeof fileNm_vmap, fmt_vmap, app_name, species_name);

  gkyl_grid_sub_array_write(
    &vvm->grid_vel, &vvm->local_vel, mt_vmap, vvm->vmap_pgkyl_host, fileNm_vmap
  );
  gkyl_msgpack_data_release(mt_vmap);

  // The cell-average file holds one value per cell per direction: p=0 data,
  // regardless of the map's I/O basis order.
  struct gkyl_msgpack_data *mt_vmap_avg = gkyl_msgpack_create(
    4,
    (struct gkyl_msgpack_map_elem[]){
      {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
      {.key = "frame", .elem_type = GKYL_MP_INT, .ival = 0},
      {.key = "polyOrder", .elem_type = GKYL_MP_INT, .ival = 0},
      {.key = "basisType", .elem_type = GKYL_MP_STRING, .cval = vvm->basis_pgkyl.id},
    }
  );

  const char *fmt_vmap_avg = "%s-%s_vmap_avg.gkyl";
  int sz_vmap_avg = gkyl_calc_strlen(fmt_vmap_avg, app_name, species_name);
  char fileNm_vmap_avg[sz_vmap_avg + 1]; // ensures no buffer overflow
  snprintf(fileNm_vmap_avg, sizeof fileNm_vmap_avg, fmt_vmap_avg, app_name, species_name);

  gkyl_grid_sub_array_write(
    &vvm->grid_vel, &vvm->local_vel, mt_vmap_avg, vvm->vmap_avg_pgkyl_host, fileNm_vmap_avg
  );
  gkyl_msgpack_data_release(mt_vmap_avg);
}

void
gkyl_vlasov_velocity_map_divide_jacobvel(
  const struct gkyl_vlasov_velocity_map *vvm, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J
)
{
  // Both representations divide nodally at Gauss-Legendre quadrature points;
  // the kernels are selected by the conf basis type. For the C^0 linear
  // representation the Jacobian is constant in the cell, so the Serendipity
  // kernels reduce to division by that constant.
#ifdef GKYL_HAVE_CUDA
  if (gkyl_vlasov_velocity_map_is_cu_dev(vvm)) {
    gkyl_vlasov_velocity_map_divide_jacobvel_cu(
      conf_basis, phase_basis, &vvm->local_vel, phase_range, vvm->jacob_vel_gauss, Jf, f_no_J
    );
    return;
  }
#endif

  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  divide_Jv_t divide_Jv;
  switch (gkyl_basis_phase_kernel_type(conf_basis, phase_basis)) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      divide_Jv = choose_ser_divide_Jv_kern(cdim, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      divide_Jv = choose_tensor_divide_Jv_kern(cdim, vdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i = 0; i < vdim; ++i) {
      idx_vel[i] = iter.idx[cdim + i];
    }
    long vidx = gkyl_range_idx(&vvm->local_vel, idx_vel);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *Jf_d = gkyl_array_cfetch(Jf, pidx);
    double *f_no_J_d = gkyl_array_fetch(f_no_J, pidx);
    divide_Jv(gkyl_array_cfetch(vvm->jacob_vel_gauss, vidx), Jf_d, f_no_J_d);
  }
}

void
gkyl_vlasov_velocity_map_rescale_jacobvel(
  const struct gkyl_vlasov_velocity_map *vvm, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf
)
{
  // Both representations rescale nodally at Gauss-Legendre quadrature points;
  // the kernels are selected by the conf basis type. For the C^0 linear
  // representation the Jacobian is constant in the cell, so the Serendipity
  // kernels reduce to multiplication by that constant.
#ifdef GKYL_HAVE_CUDA
  if (gkyl_vlasov_velocity_map_is_cu_dev(vvm)) {
    gkyl_vlasov_velocity_map_rescale_jacobvel_cu(
      conf_basis, phase_basis, &vvm->local_vel, phase_range, vvm->jacob_vel_gauss, f_no_J, Jf
    );
    return;
  }
#endif

  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  rescale_Jv_t rescale_Jv;
  switch (gkyl_basis_phase_kernel_type(conf_basis, phase_basis)) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      rescale_Jv = choose_ser_rescale_Jv_kern(cdim, vdim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      rescale_Jv = choose_tensor_rescale_Jv_kern(cdim, vdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_vel[GKYL_MAX_DIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i = 0; i < vdim; ++i) {
      idx_vel[i] = iter.idx[cdim + i];
    }
    long vidx = gkyl_range_idx(&vvm->local_vel, idx_vel);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *f_no_J_d = gkyl_array_cfetch(f_no_J, pidx);
    double *Jf_d = gkyl_array_fetch(Jf, pidx);
    rescale_Jv(gkyl_array_cfetch(vvm->jacob_vel_gauss, vidx), f_no_J_d, Jf_d);
  }
}

struct gkyl_vlasov_velocity_map *
gkyl_vlasov_velocity_map_acquire(const struct gkyl_vlasov_velocity_map *vvm)
{
  gkyl_ref_count_inc(&vvm->ref_count);
  return (struct gkyl_vlasov_velocity_map *)vvm;
}

void
gkyl_vlasov_velocity_map_release(const struct gkyl_vlasov_velocity_map *vvm)
{
  gkyl_ref_count_dec(&vvm->ref_count);
}
