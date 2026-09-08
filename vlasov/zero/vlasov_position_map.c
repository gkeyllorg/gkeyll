#include <string.h>
#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_vlasov_position_map_priv.h>

// Computational coordinates = physical coordinates mapping (default).
static void
vlasov_pmap_identity(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0];
}

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  return on_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Populate the (host) position map arrays from the input 1D mappings by
// constructing a continuous (C^0) piecewise linear representation in each
// configuration dimension: an eval_on_nodes at the Gauss-Lobatto nodes (the
// cell vertices for a p=1 basis, shared between neighboring cells so the
// representation is continuous by construction). The linear map is stored as a
// degenerate cubic in the same cdim*4 layout the velocity map's C^1 cubic uses
// (quadratic and cubic coefficients identically zero), so the post-processing
// tooling is shared. The per-direction Jacobian is piecewise constant
// (discontinuous at cell interfaces), so every quadrature-point value in
// jacob_pos, jacob_pos_surf, and jacob_pos_gauss within a cell is the same
// constant; the arrays keep the layout with that constant replicated.
//
// The arrays are defined on the EXTENDED configuration range (crange_ext) so
// boundary-streaming kernels and IC projection over the extended range have a
// well-defined map in the ghost cells. Geometry is only defined up to the
// domain boundary, so ghost cells COPY the adjacent skin cell's map and
// Jacobian (the per-direction 1D index is clamped to the interior): this keeps
// periodic boundaries symmetric and avoids extrapolating the map past the edge.
static void
vlasov_position_map_c0_linear(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange_ext,
  int c_poly_order, struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM],
  struct gkyl_array *pmap, struct gkyl_array *jacob_pos, struct gkyl_array *jacob_pos_surf,
  struct gkyl_array *pmap_pgkyl, struct gkyl_array *pmap_avg_pgkyl,
  struct gkyl_array *jacob_pos_gauss)
{
  int cdim = cgrid->ndim;
  struct gkyl_array *c_linear[3];

  // Make 1D linear basis for constructing C^0 expansion
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // 1D ranges for indexing 1D linear bases
  struct gkyl_range local[3], local_ext[3];

  // Loop over number of dimensions and construct 1D mappings
  for (int i=0; i<cdim; ++i) {
    double lower[] = { cgrid->lower[i] }, upper[] = { cgrid->upper[i] };
    int cells[] = { cgrid->cells[i] };

    struct gkyl_rect_grid grid_1d;
    gkyl_rect_grid_init(&grid_1d, 1, lower, upper, cells);

    int nghost[GKYL_MAX_CDIM] = { 0, 0 };
    gkyl_create_grid_ranges(&grid_1d, nghost, &local_ext[i], &local[i]);

    c_linear[i] = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext[i].volume);

    gkyl_eval_on_nodes *ev_map = gkyl_eval_on_nodes_new(&grid_1d, &basis, 1,
      inp_pmap[i].eval_pmap, inp_pmap[i].ctx);
    gkyl_eval_on_nodes_advance(ev_map, 0.0, &local[i], c_linear[i]);
    gkyl_eval_on_nodes_release(ev_map);
  }

  // Location of the four 1D cubic coefficients of direction d inside the
  // full-dimensional p=3 serendipity pmap_pgkyl representation, and the scale
  // factor relating them to the 1D coefficients (same placement the velocity
  // C^1 cubic kernels use; the quadratic and cubic slots just hold zeros here).
  static const int pgkyl_idx[3][3][4] = {
    { {0, 1, 2, 3}, {0}, {0} }, // 1X
    { {0, 1, 4, 8}, {0, 2, 5, 9}, {0} }, // 2X
    { {0, 1, 7, 17}, {0, 2, 8, 18}, {0, 3, 9, 19} }, // 3X
  };
  static const double pgkyl_scale[3] = { 1.0, 1.4142135623730951, 2.0 };

  int num_quad_vol = c_poly_order+1, num_quad_surf = c_poly_order+2;
  long pgkyl_nb = pmap_pgkyl->ncomp/cdim;
  int cidx_1D[1]; // 1D index for indexing correct linear mapping

  struct gkyl_range_iter iter_pmap;
  gkyl_range_iter_init(&iter_pmap, crange_ext);
  while (gkyl_range_iter_next(&iter_pmap)) {
    long loc_pos = gkyl_range_idx(crange_ext, iter_pmap.idx);

    double *pmap_d = gkyl_array_fetch(pmap, loc_pos);
    double *jacob_pos_d = gkyl_array_fetch(jacob_pos, loc_pos);
    double *jacob_pos_surf_d = gkyl_array_fetch(jacob_pos_surf, loc_pos);
    double *pmap_pgkyl_d = gkyl_array_fetch(pmap_pgkyl, loc_pos);
    double *pmap_avg_pgkyl_d = gkyl_array_fetch(pmap_avg_pgkyl, loc_pos);
    double *jacob_pos_gauss_d = gkyl_array_fetch(jacob_pos_gauss, loc_pos);

    double jacob_tot = 1.0;
    for (int i=0; i<cdim; ++i) {
      // Clamp ghost indices to the adjacent skin cell so ghost cells copy the
      // skin cell's map/Jacobian (geometry defined only up to the boundary).
      int ci = iter_pmap.idx[i];
      if (ci < 1) ci = 1;
      if (ci > cgrid->cells[i]) ci = cgrid->cells[i];
      cidx_1D[0] = ci;
      long loc_pos_1D = gkyl_range_idx(&local[i], cidx_1D);
      const double *lin = gkyl_array_cfetch(c_linear[i], loc_pos_1D);

      double *p = &pmap_d[4*i];
      p[0] = lin[0];
      p[1] = lin[1];
      p[2] = 0.0;
      p[3] = 0.0;

      // Piecewise constant Jacobian: (2/dx)*d/dz of the linear map, with
      // sqrt(3/2) from the derivative of the 1D linear basis function.
      double jac = (2.0/cgrid->dx[i])*1.224744871391589*p[1];
      jacob_tot *= jac;
      for (int q=0; q<num_quad_vol; ++q)
        jacob_pos_d[i*num_quad_vol+q] = jac;
      for (int q=0; q<num_quad_surf; ++q)
        jacob_pos_surf_d[i*num_quad_surf+q] = jac;

      for (int k=0; k<4; ++k)
        pmap_pgkyl_d[i*pgkyl_nb + pgkyl_idx[cdim-1][i][k]] = pgkyl_scale[cdim-1]*p[k];
      pmap_avg_pgkyl_d[i] = p[0]/1.4142135623730951;
    }
    for (int q=0; q<jacob_pos_gauss->ncomp; ++q)
      jacob_pos_gauss_d[q] = jacob_tot;
  }

  // free temporary memory
  for (int i=0; i<cdim; ++i)
    gkyl_array_release(c_linear[i]);
}

static void
gkyl_vlasov_position_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_position_map *vpm = container_of(ref, struct gkyl_vlasov_position_map, ref_count);

  gkyl_array_release(vpm->pmap);
  gkyl_array_release(vpm->jacob_pos);
  gkyl_array_release(vpm->jacob_pos_surf);
  gkyl_array_release(vpm->jacob_pos_gauss);

  gkyl_array_release(vpm->pmap_host);
  gkyl_array_release(vpm->jacob_pos_host);
  gkyl_array_release(vpm->jacob_pos_surf_host);
  gkyl_array_release(vpm->jacob_pos_gauss_host);

  gkyl_array_release(vpm->pmap_pgkyl_host);
  gkyl_array_release(vpm->pmap_avg_pgkyl_host);

  gkyl_free(vpm);
}

struct gkyl_vlasov_position_map*
gkyl_vlasov_position_map_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_range *crange_ext, const struct gkyl_basis *conf_basis,
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM], bool use_gpu)
{
  struct gkyl_vlasov_position_map *vpm = gkyl_malloc(sizeof(*vpm));

  int cdim = cgrid->ndim;
  int poly_order = conf_basis->poly_order;

  vpm->grid_pos = *cgrid;
  // local_pos is the interior range (used for I/O sub-range writes and conf
  // indexing of the local phase range); it is a sub-range of crange_ext and
  // shares its layout, so it indexes the (extended) arrays correctly.
  vpm->local_pos = *crange;
  vpm->local_ext_pos = *crange_ext;
  vpm->basis_pos = *conf_basis;

  vpm->is_identity = true;
  for (int c=0; c<cdim; ++c) {
    if (inp_pmap[c].eval_pmap) {
      vpm->is_identity = false;
    }
  }

  // Arrays are defined on the EXTENDED conf range so ghost cells (copies of the
  // adjacent skin cell) are available to boundary-streaming kernels and IC
  // projection over the extended range.
  // position map is stored in the degenerate-cubic layout in each direction (up to 3X; 3*4=12 components)
  vpm->pmap = mkarr(use_gpu, cdim*4, crange_ext->volume);
  // configuration-space Jacobian at quadrature points and "surface" quadrature points, used by the
  // streaming kernels to divide grad_x by the (per-cell constant) conf Jacobian in volume and surface
  // operations respectively, with surface operations utilizing more quadrature points.
  vpm->jacob_pos = mkarr(use_gpu, cdim*(poly_order+1), crange_ext->volume);
  vpm->jacob_pos_surf = mkarr(use_gpu, cdim*(poly_order+2), crange_ext->volume);
  // special basis sets to get the correct number of coefficients in 2X and 3X for constructing the
  // mapping in post-processing, as well as storing the total conf-space Jacobian at quadrature points.
  struct gkyl_basis pmap_pgkyl_basis, jacob_pos_basis;
  gkyl_cart_modal_serendip(&pmap_pgkyl_basis, cdim, 3);
  gkyl_cart_modal_tensor(&jacob_pos_basis, cdim, poly_order);
  // Keep the I/O basis so the write can label the file with the basis the data
  // is actually stored in (the app's conf basis does not describe these arrays).
  vpm->basis_pgkyl = pmap_pgkyl_basis;
  // total conf-space Jacobian at Gaussian quadrature points (constant per cell), used to divide
  // out the conf-space Jacobian from the distribution.
  vpm->jacob_pos_gauss = mkarr(use_gpu, jacob_pos_basis.num_basis, crange_ext->volume);

  // host-side mirrors for initialization and I/O
  if (use_gpu) {
    vpm->pmap_host = mkarr(false, vpm->pmap->ncomp, vpm->pmap->size);
    vpm->jacob_pos_host = mkarr(false, vpm->jacob_pos->ncomp, vpm->jacob_pos->size);
    vpm->jacob_pos_surf_host = mkarr(false, vpm->jacob_pos_surf->ncomp, vpm->jacob_pos_surf->size);
    vpm->jacob_pos_gauss_host = mkarr(false, vpm->jacob_pos_gauss->ncomp, vpm->jacob_pos_gauss->size);
  }
  else {
    vpm->pmap_host = gkyl_array_acquire(vpm->pmap);
    vpm->jacob_pos_host = gkyl_array_acquire(vpm->jacob_pos);
    vpm->jacob_pos_surf_host = gkyl_array_acquire(vpm->jacob_pos_surf);
    vpm->jacob_pos_gauss_host = gkyl_array_acquire(vpm->jacob_pos_gauss);
  }

  // position map for I/O (extended-range layout; only the interior is written)
  vpm->pmap_pgkyl_host = mkarr(false, cdim*pmap_pgkyl_basis.num_basis, crange_ext->volume);
  vpm->pmap_avg_pgkyl_host = mkarr(false, cdim, crange_ext->volume);

  // Default unspecified directions to the identity map.
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM];
  for (int c=0; c<cdim; ++c) {
    if (inp_pmap[c].eval_pmap) {
      inp[c] = inp_pmap[c];
    }
    else {
      inp[c].eval_pmap = vlasov_pmap_identity;
      inp[c].ctx = 0;
    }
  }

  // The position map supports only the C^0 piecewise linear representation.
  vlasov_position_map_c0_linear(cgrid, crange_ext, poly_order, inp,
    vpm->pmap_host, vpm->jacob_pos_host, vpm->jacob_pos_surf_host,
    vpm->pmap_pgkyl_host, vpm->pmap_avg_pgkyl_host, vpm->jacob_pos_gauss_host);

  if (use_gpu) {
    gkyl_array_copy(vpm->pmap, vpm->pmap_host);
    gkyl_array_copy(vpm->jacob_pos, vpm->jacob_pos_host);
    gkyl_array_copy(vpm->jacob_pos_surf, vpm->jacob_pos_surf_host);
    gkyl_array_copy(vpm->jacob_pos_gauss, vpm->jacob_pos_gauss_host);
  }

  vpm->flags = 0;
  if (use_gpu) {
    GKYL_SET_CU_ALLOC(vpm->flags);
  }
  vpm->ref_count = gkyl_ref_count_init(gkyl_vlasov_position_map_free);

  return vpm;
}

bool
gkyl_vlasov_position_map_is_cu_dev(const struct gkyl_vlasov_position_map *vpm)
{
  return GKYL_IS_CU_ALLOC(vpm->flags);
}

void
gkyl_vlasov_position_map_eval_mc2p(const struct gkyl_vlasov_position_map *vpm,
  const double *xc, double *xp)
{
  // Evaluate the stored (DG) position map at the computational coordinate xc.
  // The map is stored per direction in the degenerate-cubic 4-slot layout, so a
  // 1D p=3 modal basis evaluates each direction's expansion (linear here: only
  // the first two modes are nonzero). Mirrors gkyl_velocity_map_eval_c2p.
  int cdim = vpm->grid_pos.ndim;

  // Identity map: computational coordinates ARE physical coordinates. Return
  // them exactly (bitwise) instead of reconstructing x through the DG
  // expansion, so uniform grids are truly unaffected by the map machinery.
  if (vpm->is_identity) {
    for (int d=0; d<cdim; ++d) xp[d] = xc[d];
    return;
  }

  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // Find the index of the (interior) cell containing xc. grid_pos.lower is the
  // interior lower edge, so the base must be the interior range's lower index
  // (mirrors gkyl_vlasov_velocity_map_eval_c2p, which uses local_vel). Using the
  // extended range's lower here would shift every cell into its lower neighbour
  // (the first interior cell would resolve to the lower ghost). Coordinates in
  // the ghost region clamp to the nearest interior cell, so they pick up the
  // skin cell's map (the ghost-copies-skin convention).
  int idx_xc[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; ++d) {
    int idx = vpm->local_pos.lower[d] + (int) floor((xc[d] - vpm->grid_pos.lower[d])/vpm->grid_pos.dx[d]);
    idx = GKYL_MIN2(idx, vpm->local_pos.upper[d]);
    idx = GKYL_MAX2(idx, vpm->local_pos.lower[d]);
    idx_xc[d] = idx;
  }

  // Fetch DG coefficients of the position map in idx_xc.
  long lidx_xc = gkyl_range_idx(&vpm->local_ext_pos, idx_xc);
  const double *pmap_c = gkyl_array_cfetch(vpm->pmap_host, lidx_xc);

  double xc_cc[GKYL_MAX_CDIM];
  gkyl_rect_grid_cell_center(&vpm->grid_pos, idx_xc, xc_cc);
  for (int d=0; d<cdim; ++d) {
    // Convert computational to logical coord, then evaluate the expansion.
    double xlog = (xc[d] - xc_cc[d])/(0.5*vpm->grid_pos.dx[d]);
    xp[d] = b1.eval_expand(&xlog, &pmap_c[d*4]);
  }
}

void
gkyl_vlasov_position_map_write(const struct gkyl_vlasov_position_map *vpm,
  struct gkyl_comm *comm, const char *app_name, const char *name)
{
  // The map is static and written once, so both files carry frame-0 metadata.
  struct gkyl_msgpack_data *mt_pmap = gkyl_msgpack_create(4,
    (struct gkyl_msgpack_map_elem []) {
      { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
      { .key = "frame", .elem_type = GKYL_MP_INT, .ival = 0 },
      { .key = "polyOrder", .elem_type = GKYL_MP_INT, .ival = vpm->basis_pgkyl.poly_order },
      { .key = "basisType", .elem_type = GKYL_MP_STRING, .cval = vpm->basis_pgkyl.id }
    }
  );

  const char *fmt_pmap = "%s-%s_pmap.gkyl";
  int sz_pmap = gkyl_calc_strlen(fmt_pmap, app_name, name);
  char fileNm_pmap[sz_pmap+1]; // ensures no buffer overflow
  snprintf(fileNm_pmap, sizeof fileNm_pmap, fmt_pmap, app_name, name);

  gkyl_comm_array_write(comm, &vpm->grid_pos, &vpm->local_pos,
    mt_pmap, vpm->pmap_pgkyl_host, fileNm_pmap);
  gkyl_msgpack_data_release(mt_pmap);

  // The cell-average file holds one value per cell per direction: p=0 data,
  // regardless of the map's I/O basis order.
  struct gkyl_msgpack_data *mt_pmap_avg = gkyl_msgpack_create(4,
    (struct gkyl_msgpack_map_elem []) {
      { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
      { .key = "frame", .elem_type = GKYL_MP_INT, .ival = 0 },
      { .key = "polyOrder", .elem_type = GKYL_MP_INT, .ival = 0 },
      { .key = "basisType", .elem_type = GKYL_MP_STRING, .cval = vpm->basis_pgkyl.id }
    }
  );

  const char *fmt_pmap_avg = "%s-%s_pmap_avg.gkyl";
  int sz_pmap_avg = gkyl_calc_strlen(fmt_pmap_avg, app_name, name);
  char fileNm_pmap_avg[sz_pmap_avg+1]; // ensures no buffer overflow
  snprintf(fileNm_pmap_avg, sizeof fileNm_pmap_avg, fmt_pmap_avg, app_name, name);

  gkyl_comm_array_write(comm, &vpm->grid_pos, &vpm->local_pos,
    mt_pmap_avg, vpm->pmap_avg_pgkyl_host, fileNm_pmap_avg);
  gkyl_msgpack_data_release(mt_pmap_avg);
}

void
gkyl_vlasov_position_map_divide_jacobpos(const struct gkyl_vlasov_position_map *vpm,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J)
{
  // The conf Jacobian is a per-cell constant (C^0 linear map), so dividing it
  // out is an exact scalar operation on every phase-space coefficient. The
  // total conf Jacobian of a cell is held (constant-replicated) in
  // jacob_pos_gauss; its first component is that constant.
  int cdim = vpm->basis_pos.ndim;
#ifdef GKYL_HAVE_CUDA
  if (gkyl_vlasov_position_map_is_cu_dev(vpm)) {
    gkyl_vlasov_position_map_divide_jacobpos_cu(cdim, phase_basis,
      &vpm->local_pos, phase_range, vpm->jacob_pos_gauss, Jf, f_no_J);
    return;
  }
#endif

  int num_basis = phase_basis->num_basis;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_conf[GKYL_MAX_CDIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<cdim; ++i) {
      idx_conf[i] = iter.idx[i];
    }
    long cidx = gkyl_range_idx(&vpm->local_pos, idx_conf);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *jacob_pos_gauss_d = gkyl_array_cfetch(vpm->jacob_pos_gauss, cidx);
    double jacob_inv = 1.0/jacob_pos_gauss_d[0];
    const double *Jf_d = gkyl_array_cfetch(Jf, pidx);
    double *f_no_J_d = gkyl_array_fetch(f_no_J, pidx);
    for (int k=0; k<num_basis; ++k)
      f_no_J_d[k] = jacob_inv*Jf_d[k];
  }
}

void
gkyl_vlasov_position_map_rescale_jacobpos(const struct gkyl_vlasov_position_map *vpm,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf)
{
  // Inverse of gkyl_vlasov_position_map_divide_jacobpos: multiply every
  // phase-space coefficient by the per-cell constant conf Jacobian.
  int cdim = vpm->basis_pos.ndim;
#ifdef GKYL_HAVE_CUDA
  if (gkyl_vlasov_position_map_is_cu_dev(vpm)) {
    gkyl_vlasov_position_map_rescale_jacobpos_cu(cdim, phase_basis,
      &vpm->local_pos, phase_range, vpm->jacob_pos_gauss, f_no_J, Jf);
    return;
  }
#endif

  int num_basis = phase_basis->num_basis;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  int idx_conf[GKYL_MAX_CDIM];

  while (gkyl_range_iter_next(&iter)) {
    for (int i=0; i<cdim; ++i) {
      idx_conf[i] = iter.idx[i];
    }
    long cidx = gkyl_range_idx(&vpm->local_pos, idx_conf);
    long pidx = gkyl_range_idx(phase_range, iter.idx);

    const double *jacob_pos_gauss_d = gkyl_array_cfetch(vpm->jacob_pos_gauss, cidx);
    double jacob = jacob_pos_gauss_d[0];
    const double *f_no_J_d = gkyl_array_cfetch(f_no_J, pidx);
    double *Jf_d = gkyl_array_fetch(Jf, pidx);
    for (int k=0; k<num_basis; ++k)
      Jf_d[k] = jacob*f_no_J_d[k];
  }
}

void
gkyl_vlasov_position_map_divide_jacobpos_conf(const struct gkyl_vlasov_position_map *vpm,
  const struct gkyl_range *conf_range, int num_coeff_divide,
  const struct gkyl_array *Jmom, struct gkyl_array *mom_no_J)
{
  // Host-side division of a configuration-space field (e.g. a velocity moment)
  // by the per-cell constant total conf Jacobian. Moments of the stored
  // J_x J_v f carry J (the conf Jacobian factors out of the velocity integral),
  // so dividing by it yields physical moments for I/O. The map is C^0 linear so
  // J is constant in the cell and the division is exact. Only the first
  // num_coeff_divide coefficients are scaled (the J-carrying components); the
  // rest are copied through (intensive/ratio components, e.g. LTE V_drift, T/m).
  // Uses the host Jacobian since moment I/O is host-side. Device arrays are
  // staged through host copies (this is an I/O / initialization path).
  int ncomp = Jmom->ncomp;

  const struct gkyl_array *Jmom_h = Jmom;
  struct gkyl_array *mom_no_J_h = mom_no_J;
  struct gkyl_array *Jmom_tmp = 0, *mom_no_J_tmp = 0;
  if (gkyl_array_is_cu_dev(Jmom)) {
    Jmom_tmp = gkyl_array_new(Jmom->type, Jmom->ncomp, Jmom->size);
    gkyl_array_copy(Jmom_tmp, Jmom);
    Jmom_h = Jmom_tmp;
  }
  if (gkyl_array_is_cu_dev(mom_no_J)) {
    mom_no_J_tmp = (Jmom == mom_no_J) ? gkyl_array_acquire(Jmom_tmp)
      : gkyl_array_new(mom_no_J->type, mom_no_J->ncomp, mom_no_J->size);
    mom_no_J_h = mom_no_J_tmp;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx_jac = gkyl_range_idx(&vpm->local_pos, iter.idx);
    long cidx_mom = gkyl_range_idx(conf_range, iter.idx);

    const double *jacob_pos_gauss_d = gkyl_array_cfetch(vpm->jacob_pos_gauss_host, cidx_jac);
    double jacob_inv = 1.0/jacob_pos_gauss_d[0];
    const double *Jmom_d = gkyl_array_cfetch(Jmom_h, cidx_mom);
    double *mom_no_J_d = gkyl_array_fetch(mom_no_J_h, cidx_mom);
    for (int k=0; k<num_coeff_divide; ++k)
      mom_no_J_d[k] = jacob_inv*Jmom_d[k];
    for (int k=num_coeff_divide; k<ncomp; ++k)
      mom_no_J_d[k] = Jmom_d[k];
  }

  if (mom_no_J_tmp) {
    gkyl_array_copy(mom_no_J, mom_no_J_tmp);
    gkyl_array_release(mom_no_J_tmp);
  }
  if (Jmom_tmp) gkyl_array_release(Jmom_tmp);
}

void
gkyl_vlasov_position_map_rescale_jacobpos_conf(const struct gkyl_vlasov_position_map *vpm,
  const struct gkyl_range *conf_range, const struct gkyl_array *a_no_J, struct gkyl_array *Ja)
{
  // Host-side multiplication of a configuration-space field (e.g. the EM field)
  // by the per-cell constant total conf Jacobian J. Inverse of
  // gkyl_vlasov_position_map_divide_jacobpos_conf with the full component count;
  // converts a physical conf field to the J-weighted form used for evolution
  // (e.g. E -> J*E at initialization/restart). Exact (J constant in cell).
  // Device arrays are staged through host copies (initialization/restart path).
  int ncomp = a_no_J->ncomp;

  const struct gkyl_array *a_no_J_h = a_no_J;
  struct gkyl_array *Ja_h = Ja;
  struct gkyl_array *a_no_J_tmp = 0, *Ja_tmp = 0;
  if (gkyl_array_is_cu_dev(a_no_J)) {
    a_no_J_tmp = gkyl_array_new(a_no_J->type, a_no_J->ncomp, a_no_J->size);
    gkyl_array_copy(a_no_J_tmp, a_no_J);
    a_no_J_h = a_no_J_tmp;
  }
  if (gkyl_array_is_cu_dev(Ja)) {
    Ja_tmp = (a_no_J == Ja) ? gkyl_array_acquire(a_no_J_tmp)
      : gkyl_array_new(Ja->type, Ja->ncomp, Ja->size);
    Ja_h = Ja_tmp;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx_jac = gkyl_range_idx(&vpm->local_pos, iter.idx);
    long cidx_a = gkyl_range_idx(conf_range, iter.idx);

    const double *jacob_pos_gauss_d = gkyl_array_cfetch(vpm->jacob_pos_gauss_host, cidx_jac);
    double jacob = jacob_pos_gauss_d[0];
    const double *a_no_J_d = gkyl_array_cfetch(a_no_J_h, cidx_a);
    double *Ja_d = gkyl_array_fetch(Ja_h, cidx_a);
    for (int k=0; k<ncomp; ++k)
      Ja_d[k] = jacob*a_no_J_d[k];
  }

  if (Ja_tmp) {
    gkyl_array_copy(Ja, Ja_tmp);
    gkyl_array_release(Ja_tmp);
  }
  if (a_no_J_tmp) gkyl_array_release(a_no_J_tmp);
}

struct gkyl_vlasov_position_map*
gkyl_vlasov_position_map_acquire(const struct gkyl_vlasov_position_map *vpm)
{
  gkyl_ref_count_inc(&vpm->ref_count);
  return (struct gkyl_vlasov_position_map*) vpm;
}

void
gkyl_vlasov_position_map_release(const struct gkyl_vlasov_position_map *vpm)
{
  gkyl_ref_count_dec(&vpm->ref_count);
}
