#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_comm.h>
#include <gkyl_array.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_vlasov_triad_geom_preset_priv.h>

void
gkyl_vlasov_triad_geom_from_basis(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor
)
{
  int num_pt_indices[3] = {1, 6, 18};
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim);
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim);
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim);
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(vdim);

  struct gkyl_array *cov_tangent_basis_at_nodes;
  struct gkyl_array *triad_basis_at_nodes;
  struct gkyl_array *triad_basis_gradient_at_nodes;
  struct gkyl_array *h_ij_at_nodes;
  struct gkyl_array *h_ij_inv_at_nodes;
  struct gkyl_array *det_h_at_nodes;
  struct gkyl_array *conf_poisson_tensor_at_nodes;

  // Convert nodal to modal using the computed nodal quantities
  gkyl_eval_on_nodes *h_ij_proj;
  gkyl_eval_on_nodes *h_ij_inv_proj;
  gkyl_eval_on_nodes *det_h_proj;
  gkyl_eval_on_nodes *conf_poisson_tensor_proj;

  h_ij_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim * (vdim + 1) / 2, NULL, NULL);
  h_ij_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim * (vdim + 1) / 2, NULL, NULL);
  det_h_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, 1, NULL, NULL);
  conf_poisson_tensor_proj =
    gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim - 1], NULL, NULL);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int num_basis = cbasis.num_basis;

  // Eval functions
  cov_tangent_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim, num_basis);
  triad_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim, num_basis);
  triad_basis_gradient_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim * vdim, num_basis);

  // Derived quantities from eval functions
  h_ij_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * (vdim + 1) / 2, num_basis);
  h_ij_inv_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * (vdim + 1) / 2, num_basis);
  det_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_basis);
  conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim - 1], num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  // Grab the local node list
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, cbasis.num_basis);
  cbasis.node_list(gkyl_array_fetch(nodes, 0));

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);

    for (int i = 0; i < num_basis; ++i) {
      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i), cgrid->dx, xc, xmu);
      // Sample the geometry at physical coordinates on non-uniform conf meshes
      // (NULL c2p => identity, i.e. a uniform mesh).
      if (inp_triad_geom.c2p_func) {
        inp_triad_geom.c2p_func(xmu, xmu, inp_triad_geom.c2p_func_ctx);
      }

      // Evaluate the functions for basis and their gradients at a nodal point xmu
      inp_triad_geom.eval_cov_tangent_basis(
        0.0, xmu, gkyl_array_fetch(cov_tangent_basis_at_nodes, i),
        inp_triad_geom.eval_cov_tangent_basis_ctx
      );
      inp_triad_geom.eval_triad_basis(
        0.0, xmu, gkyl_array_fetch(triad_basis_at_nodes, i), inp_triad_geom.eval_triad_basis_ctx
      );
      inp_triad_geom.eval_triad_basis_gradient(
        0.0, xmu, gkyl_array_fetch(triad_basis_gradient_at_nodes, i),
        inp_triad_geom.eval_triad_basis_gradient_ctx
      );

      // Fill the remaining nodal quantities:
      double *pn_cov_tangent_basis = gkyl_array_fetch(cov_tangent_basis_at_nodes, i);
      double *pn_triad_basis = gkyl_array_fetch(triad_basis_at_nodes, i);
      double *pn_h_ij = gkyl_array_fetch(h_ij_at_nodes, i);
      double *pn_h_ij_inv = gkyl_array_fetch(h_ij_inv_at_nodes, i);
      double *pn_det_h = gkyl_array_fetch(det_h_at_nodes, i);
      double *pn_conf_poisson_tensor = gkyl_array_fetch(conf_poisson_tensor_at_nodes, i);
      double *pn_triad_basis_gradient = gkyl_array_fetch(triad_basis_gradient_at_nodes, i);

      // Compute the metric at nodal points
      compute_h_ij_from_cov_tangent_basis(pn_cov_tangent_basis, pn_h_ij);

      // Compute the metric inverse at nodal points
      compute_h_ij_inv_from_h_ij(pn_h_ij, pn_h_ij_inv);

      // Compute the sqrt(det(h_ij)) at nodal points
      compute_det_h_ij(pn_h_ij, pn_det_h);

      // Compute the Poisson Tensor in configruation space components at nodal points
      compute_triad_poisson_tensor_ij(
        pn_h_ij_inv, pn_triad_basis, pn_cov_tangent_basis, pn_triad_basis_gradient,
        pn_conf_poisson_tensor
      );
    }

    long lidx = gkyl_range_idx(crange, iter.idx);
    gkyl_eval_on_nodes_nod2mod(
      conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes,
      gkyl_array_fetch(conf_poisson_tensor, lidx)
    );
  }

  gkyl_array_release(nodes);
  gkyl_eval_on_nodes_release(h_ij_proj);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(cov_tangent_basis_at_nodes);
  gkyl_array_release(triad_basis_at_nodes);
  gkyl_array_release(triad_basis_gradient_at_nodes);
  gkyl_array_release(h_ij_at_nodes);
  gkyl_array_release(h_ij_inv_at_nodes);
  gkyl_array_release(det_h_at_nodes);
  gkyl_array_release(conf_poisson_tensor_at_nodes);
}

void
gkyl_vlasov_triad_geom_from_nodal(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor
)
{
  int num_pt_indices[3] = {1, 6, 18};
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim);
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim);
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim);
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(vdim);

  int num_basis = cbasis.num_basis;

  // Nodal inputs are per conf cell, node-major in cbasis.node_list order:
  // [node0: all components][node1: all components]...
  assert(cov_tangent_basis_nodal->ncomp == num_basis * vdim * vdim);
  assert(triad_basis_nodal->ncomp == num_basis * vdim * vdim);
  assert(triad_basis_gradient_nodal->ncomp == num_basis * vdim * vdim * vdim);

  // Derived quantities at nodes
  struct gkyl_array *h_ij_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * (vdim + 1) / 2, num_basis);
  struct gkyl_array *h_ij_inv_at_nodes =
    gkyl_array_new(GKYL_DOUBLE, vdim * (vdim + 1) / 2, num_basis);
  struct gkyl_array *det_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_basis);
  struct gkyl_array *conf_poisson_tensor_at_nodes =
    gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim - 1], num_basis);

  gkyl_eval_on_nodes *conf_poisson_tensor_proj =
    gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim - 1], NULL, NULL);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(crange, iter.idx);

    const double *cell_cov_tangent_basis = gkyl_array_cfetch(cov_tangent_basis_nodal, lidx);
    const double *cell_triad_basis = gkyl_array_cfetch(triad_basis_nodal, lidx);
    const double *cell_triad_basis_gradient = gkyl_array_cfetch(triad_basis_gradient_nodal, lidx);

    for (int i = 0; i < num_basis; ++i) {
      const double *pn_cov_tangent_basis = &cell_cov_tangent_basis[i * vdim * vdim];
      const double *pn_triad_basis = &cell_triad_basis[i * vdim * vdim];
      const double *pn_triad_basis_gradient = &cell_triad_basis_gradient[i * vdim * vdim * vdim];
      double *pn_h_ij = gkyl_array_fetch(h_ij_at_nodes, i);
      double *pn_h_ij_inv = gkyl_array_fetch(h_ij_inv_at_nodes, i);
      double *pn_det_h = gkyl_array_fetch(det_h_at_nodes, i);
      double *pn_conf_poisson_tensor = gkyl_array_fetch(conf_poisson_tensor_at_nodes, i);

      // Compute the metric at nodal points
      compute_h_ij_from_cov_tangent_basis(pn_cov_tangent_basis, pn_h_ij);

      // Compute the metric inverse at nodal points
      compute_h_ij_inv_from_h_ij(pn_h_ij, pn_h_ij_inv);

      // Compute the sqrt(det(h_ij)) at nodal points
      compute_det_h_ij(pn_h_ij, pn_det_h);

      // Compute the Poisson Tensor in configruation space components at nodal points
      compute_triad_poisson_tensor_ij(
        pn_h_ij_inv, pn_triad_basis, pn_cov_tangent_basis, pn_triad_basis_gradient,
        pn_conf_poisson_tensor
      );
    }

    gkyl_eval_on_nodes_nod2mod(
      conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes,
      gkyl_array_fetch(conf_poisson_tensor, lidx)
    );
  }

  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(h_ij_at_nodes);
  gkyl_array_release(h_ij_inv_at_nodes);
  gkyl_array_release(det_h_at_nodes);
  gkyl_array_release(conf_poisson_tensor_at_nodes);
}

void
gkyl_vlasov_triad_geom_from_nodal_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor
)
{
  int num_pt_indices[3] = {1, 6, 18};
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;
  int num_pt = num_pt_indices[vdim - 1];

  // The 2-nodes-per-direction gather below (and the interior-node layout of the
  // gyrokinetic geometry this consumes) is p=1 only.
  assert(cbasis.poly_order == 1);

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim);
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim);
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim);
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(vdim);

  int num_basis = cbasis.num_basis;

  assert(cov_tangent_basis_nodal->ncomp == vdim * vdim);
  assert(triad_basis_nodal->ncomp == vdim * vdim);
  assert(triad_basis_gradient_nodal->ncomp == vdim * vdim * vdim);

  double h_ij[6], h_ij_inv[6], det_h[1];
  double pt_at_nodes[8 * 18]; // num_basis <= 8 at p1, num_pt <= 18
  double fnodal[8];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  while (gkyl_range_iter_next(&iter)) {
    // Gather this cell's Gauss-Legendre interior nodes from the global nodal
    // range: 2 nodes per direction per cell, node i decomposed lexicographically
    // with the last dimension fastest — the same mapping gkyl_nodal_ops_n2m_interior
    // uses, so the node order matches what quad_nodal_to_modal expects.
    int nidx[GKYL_MAX_CDIM];
    for (int i = 0; i < num_basis; ++i) {
      for (int j = 0; j < cdim; ++j) {
        nidx[j] = (iter.idx[j] - crange->lower[j]) * 2 + ((i >> (cdim - 1 - j)) & 1);
      }
      long lin_nidx = gkyl_range_idx(nrange, nidx);

      const double *pn_cov_tangent_basis = gkyl_array_cfetch(cov_tangent_basis_nodal, lin_nidx);
      const double *pn_triad_basis = gkyl_array_cfetch(triad_basis_nodal, lin_nidx);
      const double *pn_triad_basis_gradient =
        gkyl_array_cfetch(triad_basis_gradient_nodal, lin_nidx);

      compute_h_ij_from_cov_tangent_basis(pn_cov_tangent_basis, h_ij);
      compute_h_ij_inv_from_h_ij(h_ij, h_ij_inv);
      compute_det_h_ij(h_ij, det_h);
      compute_triad_poisson_tensor_ij(
        h_ij_inv, pn_triad_basis, pn_cov_tangent_basis, pn_triad_basis_gradient,
        &pt_at_nodes[i * num_pt]
      );
    }

    double *pt_modal = gkyl_array_fetch(conf_poisson_tensor, gkyl_range_idx(crange, iter.idx));
    for (int c = 0; c < num_pt; ++c) {
      for (int i = 0; i < num_basis; ++i) {
        fnodal[i] = pt_at_nodes[i * num_pt + c];
      }
      for (int k = 0; k < num_basis; ++k) {
        cbasis.quad_nodal_to_modal(fnodal, &pt_modal[c * num_basis], k);
      }
    }
  }
}

// Pointwise b-aligned Gram-Schmidt triad from the Cartesian covariant tangents
// in field-aligned coordinates (the magnetic field along the third tangent):
// bhat = t_3/|t_3|, e1 = normalize(t_1 - (t_1.bhat) bhat), e2 = bhat x e1
// (right-handed by construction). Rows of tri: e1, e2, bhat. Returns the
// ratio |t_1 - (t_1.bhat) bhat|/|t_1|, which vanishes when t_1 is parallel to
// bhat and the construction is degenerate.
static double
triad_gs(const double *tan, double *tri)
{
  const double *t1 = &tan[0], *t3 = &tan[6];
  double t3mag = sqrt(t3[0] * t3[0] + t3[1] * t3[1] + t3[2] * t3[2]);
  double bhat[3] = {t3[0] / t3mag, t3[1] / t3mag, t3[2] / t3mag};

  double t1b = t1[0] * bhat[0] + t1[1] * bhat[1] + t1[2] * bhat[2];
  double e1[3] = {t1[0] - t1b * bhat[0], t1[1] - t1b * bhat[1], t1[2] - t1b * bhat[2]};
  double e1mag = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  double t1mag = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
  e1[0] /= e1mag;
  e1[1] /= e1mag;
  e1[2] /= e1mag;

  tri[0] = e1[0];
  tri[1] = e1[1];
  tri[2] = e1[2];
  tri[3] = bhat[1] * e1[2] - bhat[2] * e1[1];
  tri[4] = bhat[2] * e1[0] - bhat[0] * e1[2];
  tri[5] = bhat[0] * e1[1] - bhat[1] * e1[0];
  tri[6] = bhat[0];
  tri[7] = bhat[1];
  tri[8] = bhat[2];

  return e1mag / t1mag;
}

// Derivative of the b-aligned Gram-Schmidt triad along one direction by the
// chain rule: tan and tri as in triad_gs, dtan = d(tangents) along that
// direction (9: row a = d t_a), dtri = d(triad) (9: row a = d sigma_a).
static void
triad_gs_grad(const double *tan, const double *tri, const double *dtan, double *dtri)
{
  const double *t1 = &tan[0], *t3 = &tan[6];
  const double *dt1 = &dtan[0], *dt3 = &dtan[6];
  const double *e1 = &tri[0], *b = &tri[6];
  double t3mag = sqrt(t3[0] * t3[0] + t3[1] * t3[1] + t3[2] * t3[2]);

  // d bhat = (d t_3 - bhat (bhat . d t_3))/|t_3|
  double bdt3 = b[0] * dt3[0] + b[1] * dt3[1] + b[2] * dt3[2];
  double db[3] = {
    (dt3[0] - b[0] * bdt3) / t3mag, (dt3[1] - b[1] * bdt3) / t3mag, (dt3[2] - b[2] * bdt3) / t3mag
  };

  // u = t_1 - c bhat with c = t_1 . bhat; e1 = u/|u|.
  double c = t1[0] * b[0] + t1[1] * b[1] + t1[2] * b[2];
  double dc =
    dt1[0] * b[0] + dt1[1] * b[1] + dt1[2] * b[2] + t1[0] * db[0] + t1[1] * db[1] + t1[2] * db[2];
  double u[3] = {t1[0] - c * b[0], t1[1] - c * b[1], t1[2] - c * b[2]};
  double umag = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  double du[3] = {
    dt1[0] - dc * b[0] - c * db[0], dt1[1] - dc * b[1] - c * db[1], dt1[2] - dc * b[2] - c * db[2]
  };
  double e1du = e1[0] * du[0] + e1[1] * du[1] + e1[2] * du[2];
  double de1[3] = {
    (du[0] - e1[0] * e1du) / umag, (du[1] - e1[1] * e1du) / umag, (du[2] - e1[2] * e1du) / umag
  };

  dtri[0] = de1[0];
  dtri[1] = de1[1];
  dtri[2] = de1[2];
  // d e2 = d bhat x e1 + bhat x d e1
  dtri[3] = db[1] * e1[2] - db[2] * e1[1] + b[1] * de1[2] - b[2] * de1[1];
  dtri[4] = db[2] * e1[0] - db[0] * e1[2] + b[2] * de1[0] - b[0] * de1[2];
  dtri[5] = db[0] * e1[1] - db[1] * e1[0] + b[0] * de1[1] - b[1] * de1[0];
  dtri[6] = db[0];
  dtri[7] = db[1];
  dtri[8] = db[2];
}

// Construction-time checks of the b-aligned triad at one node: Gram-Schmidt
// degeneracy (t_1 parallel to bhat) and, when bhat_check is given, agreement
// of bhat with the independently computed field direction. Failures print
// and, with exit_at_checks, assert.
static void
triad_gs_check(
  double gs_ratio, const double *tri, const double *bhat_check, const int *idx, int ndim,
  bool exit_at_checks
)
{
  int nidx[3] = {0, 0, 0};
  for (int d = 0; d < ndim; ++d) {
    nidx[d] = idx[d];
  }
  if (gs_ratio < 1e-10) {
    fprintf(
      stderr,
      "vlasov_triad_geom: first tangent is (near) parallel to bhat at node (%d,%d,%d); b-aligned triad is degenerate\n",
      nidx[0], nidx[1], nidx[2]
    );
    assert(!exit_at_checks);
  }
  if (bhat_check) {
    const double *bhat = &tri[6];
    double diff = sqrt(
      (bhat[0] - bhat_check[0]) * (bhat[0] - bhat_check[0]) +
      (bhat[1] - bhat_check[1]) * (bhat[1] - bhat_check[1]) +
      (bhat[2] - bhat_check[2]) * (bhat[2] - bhat_check[2])
    );
    if (diff > 1e-10) {
      fprintf(
        stderr,
        "vlasov_triad_geom: bhat from tangents differs from provided field direction by %.3e at node (%d,%d,%d)\n",
        diff, nidx[0], nidx[1], nidx[2]
      );
      assert(!exit_at_checks);
    }
  }
}

void
gkyl_vlasov_triad_geom_from_tangents_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *bhat_check_nodal,
  bool exit_at_checks, struct gkyl_array *conf_poisson_tensor
)
{
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // The b-aligned Gram-Schmidt construction is a genuine 3v triad, and the
  // finite-difference gradients need all three conf directions on the grid.
  // Deflated (axisymmetric 1x/2x) support requires analytic d/dalpha terms
  // and is not implemented yet.
  assert(vdim == 3);
  assert(cdim == 3);
  assert(cbasis.poly_order == 1);
  assert(cov_tangent_basis_nodal->ncomp == 9);

  struct gkyl_array *triad_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array *triad_grad_nodal = gkyl_array_new(GKYL_DOUBLE, 27, nrange->volume);

  // Pass 1: pointwise b-aligned Gram-Schmidt triad from the tangents.
  // The magnetic field is aligned with the third tangent (field-aligned
  // coordinates), so bhat = e_3/|e_3| needs no field input; the optional
  // bhat_check_nodal (e.g. gyrokinetic bcart_nodal) guards that assumption.
  struct gkyl_range_iter niter;
  gkyl_range_iter_init(&niter, nrange);
  while (gkyl_range_iter_next(&niter)) {
    long nl = gkyl_range_idx(nrange, niter.idx);
    const double *tan = gkyl_array_cfetch(cov_tangent_basis_nodal, nl);
    double *tri = gkyl_array_fetch(triad_nodal, nl);
    double gs_ratio = triad_gs(tan, tri);
    // Cross-check bhat against the independently computed field direction
    // (guards divergence of the geometry provider and this construction).
    const double *bc = bhat_check_nodal ? gkyl_array_cfetch(bhat_check_nodal, nl) : NULL;
    triad_gs_check(gs_ratio, tri, bc, niter.idx, 3, exit_at_checks);
  }

  // Pass 2: gradients d(triad)/dz^j by finite differences across the interior
  // node grid.
  gkyl_nodal_ops_grad_interior(cgrid, nrange, 9, triad_nodal, triad_grad_nodal);

  gkyl_vlasov_triad_geom_from_nodal_interior(
    cgrid, crange, cbasis, pgrid, prange, pbasis, nrange, cov_tangent_basis_nodal, triad_nodal,
    triad_grad_nodal, conf_poisson_tensor
  );

  gkyl_array_release(triad_nodal);
  gkyl_array_release(triad_grad_nodal);
}

void
gkyl_vlasov_triad_geom_from_tangents_grad_interior(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis, const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal,
  const struct gkyl_array *cov_tangent_basis_grad_nodal, const struct gkyl_array *bhat_check_nodal,
  const int *dir_map, bool exit_at_checks, struct gkyl_array *conf_poisson_tensor
)
{
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // The b-aligned Gram-Schmidt construction is a genuine 3v triad.
  assert(vdim == 3);
  assert(cdim <= 3);
  assert(cbasis.poly_order == 1);
  assert(cov_tangent_basis_nodal->ncomp == 9);
  assert(cov_tangent_basis_grad_nodal->ncomp == 27);

  int map[3] = {0, 1, 2};
  if (dir_map) {
    for (int q = 0; q < 3; ++q) {
      map[q] = dir_map[q];
    }
  }

  // Tangents, triad and triad gradients at the nodes with the coordinate index
  // (tangent row, gradient direction) in Poisson-tensor order.
  struct gkyl_array *tangent_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array *triad_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange->volume);
  struct gkyl_array *triad_grad_nodal = gkyl_array_new(GKYL_DOUBLE, 27, nrange->volume);

  struct gkyl_range_iter niter;
  gkyl_range_iter_init(&niter, nrange);
  while (gkyl_range_iter_next(&niter)) {
    long nl = gkyl_range_idx(nrange, niter.idx);
    const double *tan = gkyl_array_cfetch(cov_tangent_basis_nodal, nl);
    const double *dtan = gkyl_array_cfetch(cov_tangent_basis_grad_nodal, nl);
    double *tan_p = gkyl_array_fetch(tangent_nodal, nl);
    double *tri = gkyl_array_fetch(triad_nodal, nl);
    double *grad = gkyl_array_fetch(triad_grad_nodal, nl);

    // The triad is a physical frame (e1, e2, bhat): built in the 3D coordinate
    // order, in which the third tangent is field-aligned.
    double gs_ratio = triad_gs(tan, tri);
    const double *bc = bhat_check_nodal ? gkyl_array_cfetch(bhat_check_nodal, nl) : NULL;
    triad_gs_check(gs_ratio, tri, bc, niter.idx, nrange->ndim, exit_at_checks);

    // Coordinate indices reordered to Poisson-tensor order: the tangent row and
    // the direction of differentiation transform together, so the metric and
    // the connection terms assembled from them are simply re-indexed.
    for (int q = 0; q < 3; ++q) {
      int j = map[q];
      for (int c = 0; c < 3; ++c) {
        tan_p[3 * q + c] = tan[3 * j + c];
      }
      triad_gs_grad(tan, tri, &dtan[9 * j], &grad[9 * q]);
    }
  }

  gkyl_vlasov_triad_geom_from_nodal_interior(
    cgrid, crange, cbasis, pgrid, prange, pbasis, nrange, tangent_nodal, triad_nodal,
    triad_grad_nodal, conf_poisson_tensor
  );

  gkyl_array_release(tangent_nodal);
  gkyl_array_release(triad_nodal);
  gkyl_array_release(triad_grad_nodal);
}

// Index of interior node `q` (lexicographic, first direction slowest) of cell idx
// in a nodal range whose first cell is cell_lo.
static void
nodal_interior_idx(int cdim, const int *idx, const int *cell_lo, int q, int *nidx)
{
  for (int j = 0; j < cdim; ++j) {
    nidx[j] = (idx[j] - cell_lo[j]) * 2 + ((q >> (cdim - 1 - j)) & 1);
  }
}

void
gkyl_vlasov_triad_geom_nodal_from_tangents(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_range *global, struct gkyl_comm *comm,
  bool use_gpu, const struct gkyl_range *nrange_int, const struct gkyl_array *tan_nodal,
  const struct gkyl_array *tan_grad_nodal, const struct gkyl_array *bhat_nodal,
  const struct gkyl_range *nrange_surf, struct gkyl_array *const *tan_surf_nodal,
  const int *dir_map, bool exit_at_checks, struct gkyl_array *pt_vol, struct gkyl_array *pt_surf
)
{
  const int vdim = 3, num_pt = 18;
  int cdim = cgrid->ndim;
  int num_nodes = 1 << cdim, num_face_nodes = 1 << (cdim - 1);
  assert(cdim == 3 || tan_grad_nodal);
  assert(tan_nodal->ncomp == 9);
  assert(pt_vol->ncomp == num_nodes * num_pt);
  assert(pt_surf->ncomp == cdim * 2 * num_face_nodes * vdim);

  int map[3] = {0, 1, 2};
  if (dir_map) {
    for (int q = 0; q < 3; ++q) {
      map[q] = dir_map[q];
    }
  }

  metric_t compute_h_ij = choose_metric_kern(vdim);
  metric_inv_t compute_h_ij_inv = choose_metric_inv_kern(vdim);
  conf_poisson_tensor_t compute_pt = choose_conf_poisson_tensor_kern(vdim);

  // Pass 1: the triad at the interior nodes of the local cells, stored per cell
  // ([node][9]) so that the halo can be exchanged like any conf-space array.
  struct gkyl_array *triad_cell = gkyl_array_new(GKYL_DOUBLE, num_nodes * 9, local_ext->volume);
  gkyl_array_clear(triad_cell, 0.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    double *tri_c = gkyl_array_fetch(triad_cell, gkyl_range_idx(local_ext, iter.idx));
    for (int q = 0; q < num_nodes; ++q) {
      int nidx[GKYL_MAX_CDIM];
      nodal_interior_idx(cdim, iter.idx, local->lower, q, nidx);
      long nl = gkyl_range_idx(nrange_int, nidx);
      double gs_ratio = triad_gs(gkyl_array_cfetch(tan_nodal, nl), &tri_c[9 * q]);
      const double *bc = bhat_nodal ? gkyl_array_cfetch(bhat_nodal, nl) : NULL;
      triad_gs_check(gs_ratio, &tri_c[9 * q], bc, nidx, cdim, exit_at_checks);
    }
  }
  if (use_gpu) {
    struct gkyl_array *triad_cell_dev =
      gkyl_array_cu_dev_new(GKYL_DOUBLE, num_nodes * 9, local_ext->volume);
    gkyl_array_copy(triad_cell_dev, triad_cell);
    gkyl_comm_array_sync(comm, local, local_ext, triad_cell_dev);
    gkyl_array_copy(triad_cell, triad_cell_dev);
    gkyl_array_release(triad_cell_dev);
  } else {
    gkyl_comm_array_sync(comm, local, local_ext, triad_cell);
  }

  // Pass 2: triad gradients along the grid directions by finite differences over
  // the interior nodes of the local cells plus one exchanged cell across each
  // internal decomposition boundary. Ends at the global boundary stay one sided,
  // exactly as in an undecomposed run.
  int cell_lo[GKYL_MAX_CDIM], cell_up[GKYL_MAX_CDIM], nshape[GKYL_MAX_CDIM];
  for (int j = 0; j < cdim; ++j) {
    cell_lo[j] = local->lower[j] - (local->lower[j] > global->lower[j] ? 1 : 0);
    cell_up[j] = local->upper[j] + (local->upper[j] < global->upper[j] ? 1 : 0);
    nshape[j] = 2 * (cell_up[j] - cell_lo[j] + 1);
  }
  struct gkyl_range halo_cells, nrange_halo;
  gkyl_sub_range_init(&halo_cells, local_ext, cell_lo, cell_up);
  gkyl_range_init_from_shape(&nrange_halo, cdim, nshape);
  struct gkyl_array *triad_nodal = gkyl_array_new(GKYL_DOUBLE, 9, nrange_halo.volume);
  struct gkyl_array *triad_grad_nodal = gkyl_array_new(GKYL_DOUBLE, 9 * cdim, nrange_halo.volume);
  gkyl_range_iter_init(&iter, &halo_cells);
  while (gkyl_range_iter_next(&iter)) {
    const double *tri_c = gkyl_array_cfetch(triad_cell, gkyl_range_idx(local_ext, iter.idx));
    for (int q = 0; q < num_nodes; ++q) {
      int nidx[GKYL_MAX_CDIM];
      nodal_interior_idx(cdim, iter.idx, cell_lo, q, nidx);
      memcpy(
        gkyl_array_fetch(triad_nodal, gkyl_range_idx(&nrange_halo, nidx)), &tri_c[9 * q],
        9 * sizeof(double)
      );
    }
  }
  gkyl_nodal_ops_grad_interior(cgrid, &nrange_halo, 9, triad_nodal, triad_grad_nodal);

  // Pass 3: the Poisson tensor at the interior nodes. Tangent rows and gradient
  // directions are in Poisson-tensor order (grid directions first).
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    long lin = gkyl_range_idx(local_ext, iter.idx);
    const double *tri_c = gkyl_array_cfetch(triad_cell, lin);
    double *pt_c = gkyl_array_fetch(pt_vol, lin);
    for (int q = 0; q < num_nodes; ++q) {
      int nidx[GKYL_MAX_CDIM], hidx[GKYL_MAX_CDIM];
      nodal_interior_idx(cdim, iter.idx, local->lower, q, nidx);
      nodal_interior_idx(cdim, iter.idx, cell_lo, q, hidx);
      long nl = gkyl_range_idx(nrange_int, nidx);
      const double *tan = gkyl_array_cfetch(tan_nodal, nl);
      const double *tri = &tri_c[9 * q];
      const double *grad_grid =
        gkyl_array_cfetch(triad_grad_nodal, gkyl_range_idx(&nrange_halo, hidx));

      double tan_p[9], grad[27], h_ij[6], h_ij_inv[6];
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          tan_p[3 * r + c] = tan[3 * map[r] + c];
        }
        if (r < cdim) {
          memcpy(&grad[9 * r], &grad_grid[9 * r], 9 * sizeof(double));
        } else {
          triad_gs_grad(
            tan, tri, &((const double *)gkyl_array_cfetch(tan_grad_nodal, nl))[9 * map[r]],
            &grad[9 * r]
          );
        }
      }
      compute_h_ij(tan_p, h_ij);
      compute_h_ij_inv(h_ij, h_ij_inv);
      compute_pt(h_ij_inv, tri, tan_p, grad, &pt_c[q * num_pt]);
    }
  }

  // Faces: row dir of the conf-velocity block from the geometry at the face nodes
  // (no gradients enter it).
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_range_iter_init(&iter, local);
    while (gkyl_range_iter_next(&iter)) {
      double *ps = gkyl_array_fetch(pt_surf, gkyl_range_idx(local_ext, iter.idx));
      for (int e = 0; e < 2; ++e) {
        for (int s = 0; s < num_face_nodes; ++s) {
          int nidx[GKYL_MAX_CDIM];
          for (int j = 0, t = 0; j < cdim; ++j) {
            if (j == dir) {
              nidx[j] = iter.idx[j] - local->lower[j] + e;
            } else {
              nidx[j] = (iter.idx[j] - local->lower[j]) * 2 + ((s >> (cdim - 2 - t)) & 1);
              t += 1;
            }
          }
          const double *tan =
            gkyl_array_cfetch(tan_surf_nodal[dir], gkyl_range_idx(&nrange_surf[dir], nidx));
          double tri[9], tan_p[9], grad[27] = {0.0}, h_ij[6], h_ij_inv[6], pt[18];
          double gs_ratio = triad_gs(tan, tri);
          triad_gs_check(gs_ratio, tri, NULL, nidx, cdim, exit_at_checks);
          for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
              tan_p[3 * r + c] = tan[3 * map[r] + c];
            }
          }
          compute_h_ij(tan_p, h_ij);
          compute_h_ij_inv(h_ij, h_ij_inv);
          compute_pt(h_ij_inv, tri, tan_p, grad, pt);
          for (int b = 0; b < vdim; ++b) {
            ps[((dir * 2 + e) * num_face_nodes + s) * vdim + b] = pt[vdim * dir + b];
          }
        }
      }
    }
  }

  gkyl_array_release(triad_cell);
  gkyl_array_release(triad_nodal);
  gkyl_array_release(triad_grad_nodal);
}

// Configuration-space Poisson tensor at the logical point xlog of the cell with
// center xc, from the geometry evaluation functions.
static void
triad_pt_at_point(
  const struct gkyl_rect_grid *cgrid, const double *xc, const double *xlog, int vdim,
  const struct gkyl_vlasov_triad_geom_inp *inp, double *pt
)
{
  double xmu[GKYL_MAX_DIM];
  log_to_comp(cgrid->ndim, xlog, cgrid->dx, xc, xmu);
  if (inp->c2p_func) {
    inp->c2p_func(xmu, xmu, inp->c2p_func_ctx);
  }
  if (inp->use_vierbein || inp->use_preset_geom) {
    double vb[9], vb_inv[9], vb_grad[27];
    inp->eval_vierbein(0.0, xmu, vb, inp->eval_vierbein_ctx);
    inp->eval_vierbein_gradient(0.0, xmu, vb_grad, inp->eval_vierbein_gradient_ctx);
    choose_vierbein_inv_kern(vdim)(vb, vb_inv);
    choose_conf_poisson_tensor_vierbein_kern(vdim)(vb_inv, vb_grad, pt);
  } else {
    double tan[9], tri[9], tri_grad[27], h_ij[6], h_ij_inv[6];
    inp->eval_cov_tangent_basis(0.0, xmu, tan, inp->eval_cov_tangent_basis_ctx);
    inp->eval_triad_basis(0.0, xmu, tri, inp->eval_triad_basis_ctx);
    inp->eval_triad_basis_gradient(0.0, xmu, tri_grad, inp->eval_triad_basis_gradient_ctx);
    choose_metric_kern(vdim)(tan, h_ij);
    choose_metric_inv_kern(vdim)(h_ij, h_ij_inv);
    choose_conf_poisson_tensor_kern(vdim)(h_ij_inv, tri, tan, tri_grad, pt);
  }
}

void
gkyl_vlasov_triad_geom_nodal_new(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, int vdim,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *pt_vol,
  struct gkyl_array *pt_surf
)
{
  int num_pt_indices[3] = {1, 6, 18};
  int cdim = cgrid->ndim, num_pt = num_pt_indices[vdim - 1];
  int num_nodes = 1 << cdim, num_face_nodes = 1 << (cdim - 1);
  assert(pt_vol->ncomp == num_nodes * num_pt);
  assert(pt_surf->ncomp == cdim * 2 * num_face_nodes * vdim);
  if (inp_triad_geom.use_preset_geom) {
    inp_triad_geom.eval_vierbein =
      choose_vierbein_kern(inp_triad_geom.triad_preset_geom_type, vdim);
    inp_triad_geom.eval_vierbein_gradient =
      choose_vierbein_gradient_kern(inp_triad_geom.triad_preset_geom_type, vdim);
    assert(inp_triad_geom.eval_vierbein);
    assert(inp_triad_geom.eval_vierbein_gradient);
  }

  const double glpt = 1.0 / sqrt(3.0);
  double xc[GKYL_MAX_DIM], xlog[GKYL_MAX_DIM], pt[18];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);
    long lin = gkyl_range_idx(crange, iter.idx);
    double *pv = gkyl_array_fetch(pt_vol, lin);
    double *ps = gkyl_array_fetch(pt_surf, lin);
    for (int q = 0; q < num_nodes; ++q) {
      for (int j = 0; j < cdim; ++j) {
        xlog[j] = ((q >> (cdim - 1 - j)) & 1) ? glpt : -glpt;
      }
      triad_pt_at_point(cgrid, xc, xlog, vdim, &inp_triad_geom, &pv[q * num_pt]);
    }
    for (int dir = 0; dir < cdim; ++dir) {
      for (int e = 0; e < 2; ++e) {
        for (int s = 0; s < num_face_nodes; ++s) {
          for (int j = 0, t = 0; j < cdim; ++j) {
            if (j == dir) {
              xlog[j] = e == 0 ? -1.0 : 1.0;
            } else {
              xlog[j] = ((s >> (cdim - 2 - t)) & 1) ? glpt : -glpt;
              t += 1;
            }
          }
          triad_pt_at_point(cgrid, xc, xlog, vdim, &inp_triad_geom, pt);
          for (int b = 0; b < vdim; ++b) {
            ps[((dir * 2 + e) * num_face_nodes + s) * vdim + b] = pt[vdim * dir + b];
          }
        }
      }
    }
  }
}

void
gkyl_vlasov_triad_geom_from_vierbein(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor
)
{
  int num_pt_indices[3] = {1, 6, 18};
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // Choose functions
  vierbein_inv_t compute_vierbein_inv_from_vierbein = choose_vierbein_inv_kern(vdim);
  conf_poisson_tensor_vierbein_t compute_triad_poisson_tensor_ij =
    choose_conf_poisson_tensor_vierbein_kern(vdim);

  struct gkyl_array *vierbein_at_nodes;
  struct gkyl_array *vierbein_inverse_at_nodes;
  struct gkyl_array *vierbein_gradient_at_nodes;
  struct gkyl_array *conf_poisson_tensor_at_nodes;

  // Convert nodal to modal using the computed nodal quantities
  gkyl_eval_on_nodes *vierbein_proj;
  gkyl_eval_on_nodes *vierbein_inv_proj;
  gkyl_eval_on_nodes *vierbein_gradient_proj;
  gkyl_eval_on_nodes *conf_poisson_tensor_proj;

  vierbein_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim * vdim, NULL, NULL);
  vierbein_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim * vdim, NULL, NULL);
  vierbein_gradient_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim * vdim * vdim, NULL, NULL);
  conf_poisson_tensor_proj =
    gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim - 1], NULL, NULL);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int num_basis = cbasis.num_basis;

  // Eval functions
  vierbein_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim, num_basis);
  vierbein_gradient_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim * vdim, num_basis);

  // Derived quantities from eval functions
  vierbein_inverse_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim * vdim * vdim, num_basis);
  conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim - 1], num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  // Grab the local node list
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, cbasis.num_basis);
  cbasis.node_list(gkyl_array_fetch(nodes, 0));

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);

    for (int i = 0; i < num_basis; ++i) {
      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i), cgrid->dx, xc, xmu);
      // Sample the geometry at physical coordinates on non-uniform conf meshes
      // (NULL c2p => identity, i.e. a uniform mesh).
      if (inp_triad_geom.c2p_func) {
        inp_triad_geom.c2p_func(xmu, xmu, inp_triad_geom.c2p_func_ctx);
      }

      // Evaluate the functions for basis and their gradients at a nodal point xmu
      inp_triad_geom.eval_vierbein(
        0.0, xmu, gkyl_array_fetch(vierbein_at_nodes, i), inp_triad_geom.eval_vierbein_ctx
      );
      inp_triad_geom.eval_vierbein_gradient(
        0.0, xmu, gkyl_array_fetch(vierbein_gradient_at_nodes, i),
        inp_triad_geom.eval_vierbein_gradient_ctx
      );

      // Fill the remaining nodal quantities:
      double *pn_vierbein = gkyl_array_fetch(vierbein_at_nodes, i);
      double *pn_vierbein_gradient = gkyl_array_fetch(vierbein_gradient_at_nodes, i);
      double *pn_vierbein_inv = gkyl_array_fetch(vierbein_inverse_at_nodes, i);
      double *pn_conf_poisson_tensor = gkyl_array_fetch(conf_poisson_tensor_at_nodes, i);

      // Compute the vierbein upstairs componets via inverse at nodal points
      compute_vierbein_inv_from_vierbein(pn_vierbein, pn_vierbein_inv);

      // Compute the Poisson Tensor in configruation space components at nodal points
      compute_triad_poisson_tensor_ij(pn_vierbein_inv, pn_vierbein_gradient, pn_conf_poisson_tensor);
    }

    long lidx = gkyl_range_idx(crange, iter.idx);
    gkyl_eval_on_nodes_nod2mod(
      conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes,
      gkyl_array_fetch(conf_poisson_tensor, lidx)
    );
  }

  gkyl_array_release(nodes);
  gkyl_eval_on_nodes_release(vierbein_proj);
  gkyl_eval_on_nodes_release(vierbein_inv_proj);
  gkyl_eval_on_nodes_release(vierbein_gradient_proj);
  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(vierbein_at_nodes);
  gkyl_array_release(vierbein_inverse_at_nodes);
  gkyl_array_release(vierbein_gradient_at_nodes);
  gkyl_array_release(conf_poisson_tensor_at_nodes);
}

void
gkyl_vlasov_triad_geom_new(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_basis cbasis, const struct gkyl_rect_grid *pgrid,
  const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor
)
{
  // Choose between constructing from Vierbeins or constructing from basis vectors
  if (inp_triad_geom.use_preset_geom) {
    inp_triad_geom.eval_vierbein =
      choose_vierbein_kern(inp_triad_geom.triad_preset_geom_type, pgrid->ndim - cgrid->ndim);
    inp_triad_geom.eval_vierbein_gradient = choose_vierbein_gradient_kern(
      inp_triad_geom.triad_preset_geom_type, pgrid->ndim - cgrid->ndim
    );

    // ensure non-NULL pointers, requested triad geometry is not currently supported
    assert(inp_triad_geom.eval_vierbein);
    assert(inp_triad_geom.eval_vierbein_gradient);
  }

  if (inp_triad_geom.use_vierbein || inp_triad_geom.use_preset_geom) {
    gkyl_vlasov_triad_geom_from_vierbein(
      cgrid, crange, cbasis, pgrid, prange, pbasis, inp_triad_geom, conf_poisson_tensor
    );
  } else {
    gkyl_vlasov_triad_geom_from_basis(
      cgrid, crange, cbasis, pgrid, prange, pbasis, inp_triad_geom, conf_poisson_tensor
    );
  }
}

evalf_t
gkyl_vlasov_triad_preset_hamil(
  const int cdim, const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type
)
{
  return choose_hamil_kern(preset_geom_type, cdim, vdim);
}

evalf_t
gkyl_vlasov_triad_preset_vierbein(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_vierbein_kern(preset_geom_type, vdim);
}

evalf_t
gkyl_vlasov_triad_preset_vierbein_inv(
  const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type
)
{
  return choose_vierbein_inv_preset_kern(preset_geom_type, vdim);
}
