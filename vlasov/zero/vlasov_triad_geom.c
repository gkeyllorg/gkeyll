#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>
#include <gkyl_vlasov_triad_geom_preset_priv.h>


void
gkyl_vlasov_triad_geom_from_basis(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor)
{

  int num_pt_indices[3] = { 1 , 6, 18 };
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

  h_ij_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  h_ij_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  det_h_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, 1, NULL, NULL);
  conf_poisson_tensor_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim-1], NULL, NULL);


  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int num_basis = cbasis.num_basis;

  // Eval functions
  cov_tangent_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, num_basis);
  triad_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, num_basis);
  triad_basis_gradient_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim*vdim, num_basis);

  // Derived quantities from eval functions
  h_ij_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  h_ij_inv_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  det_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_basis);
  conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim-1], num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  
  // Grab the local node list
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, cbasis.num_basis);
  cbasis.node_list(gkyl_array_fetch(nodes, 0));

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {

      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i),
        cgrid->dx, xc, xmu);
      // Sample the geometry at physical coordinates on non-uniform conf meshes
      // (NULL c2p => identity, i.e. a uniform mesh).
      if (inp_triad_geom.c2p_func)
        inp_triad_geom.c2p_func(xmu, xmu, inp_triad_geom.c2p_func_ctx);

      // Evaluate the functions for basis and their gradients at a nodal point xmu
      inp_triad_geom.eval_cov_tangent_basis(0.0, xmu, gkyl_array_fetch(cov_tangent_basis_at_nodes, i), inp_triad_geom.eval_cov_tangent_basis_ctx);
      inp_triad_geom.eval_triad_basis(0.0, xmu, gkyl_array_fetch(triad_basis_at_nodes, i), inp_triad_geom.eval_triad_basis_ctx);
      inp_triad_geom.eval_triad_basis_gradient(0.0, xmu, gkyl_array_fetch(triad_basis_gradient_at_nodes, i), inp_triad_geom.eval_triad_basis_gradient_ctx);

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
      compute_triad_poisson_tensor_ij(pn_h_ij_inv, pn_triad_basis, pn_cov_tangent_basis,
        pn_triad_basis_gradient, pn_conf_poisson_tensor);
        
    }

    long lidx = gkyl_range_idx(crange, iter.idx);
    gkyl_eval_on_nodes_nod2mod(conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes, gkyl_array_fetch(conf_poisson_tensor, lidx));
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
gkyl_vlasov_triad_geom_from_nodal(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor)
{

  int num_pt_indices[3] = { 1 , 6, 18 };
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
  assert(cov_tangent_basis_nodal->ncomp == num_basis*vdim*vdim);
  assert(triad_basis_nodal->ncomp == num_basis*vdim*vdim);
  assert(triad_basis_gradient_nodal->ncomp == num_basis*vdim*vdim*vdim);

  // Derived quantities at nodes
  struct gkyl_array *h_ij_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  struct gkyl_array *h_ij_inv_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  struct gkyl_array *det_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_basis);
  struct gkyl_array *conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim-1], num_basis);

  gkyl_eval_on_nodes *conf_poisson_tensor_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim-1], NULL, NULL);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(crange, iter.idx);

    const double *cell_cov_tangent_basis = gkyl_array_cfetch(cov_tangent_basis_nodal, lidx);
    const double *cell_triad_basis = gkyl_array_cfetch(triad_basis_nodal, lidx);
    const double *cell_triad_basis_gradient = gkyl_array_cfetch(triad_basis_gradient_nodal, lidx);

    for (int i=0; i<num_basis; ++i) {

      const double *pn_cov_tangent_basis = &cell_cov_tangent_basis[i*vdim*vdim];
      const double *pn_triad_basis = &cell_triad_basis[i*vdim*vdim];
      const double *pn_triad_basis_gradient = &cell_triad_basis_gradient[i*vdim*vdim*vdim];
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
      compute_triad_poisson_tensor_ij(pn_h_ij_inv, pn_triad_basis, pn_cov_tangent_basis,
        pn_triad_basis_gradient, pn_conf_poisson_tensor);

    }

    gkyl_eval_on_nodes_nod2mod(conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes, gkyl_array_fetch(conf_poisson_tensor, lidx));
  }

  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(h_ij_at_nodes);
  gkyl_array_release(h_ij_inv_at_nodes);
  gkyl_array_release(det_h_at_nodes);
  gkyl_array_release(conf_poisson_tensor_at_nodes);
}

void
gkyl_vlasov_triad_geom_from_nodal_interior(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *triad_basis_nodal,
  const struct gkyl_array *triad_basis_gradient_nodal, struct gkyl_array *conf_poisson_tensor)
{

  int num_pt_indices[3] = { 1 , 6, 18 };
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;
  int num_pt = num_pt_indices[vdim-1];

  // The 2-nodes-per-direction gather below (and the interior-node layout of the
  // gyrokinetic geometry this consumes) is p=1 only.
  assert(cbasis.poly_order == 1);

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim);
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim);
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim);
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(vdim);

  int num_basis = cbasis.num_basis;

  assert(cov_tangent_basis_nodal->ncomp == vdim*vdim);
  assert(triad_basis_nodal->ncomp == vdim*vdim);
  assert(triad_basis_gradient_nodal->ncomp == vdim*vdim*vdim);

  double h_ij[6], h_ij_inv[6], det_h[1];
  double pt_at_nodes[8*18]; // num_basis <= 8 at p1, num_pt <= 18
  double fnodal[8];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);

  while (gkyl_range_iter_next(&iter)) {

    // Gather this cell's Gauss-Legendre interior nodes from the global nodal
    // range: 2 nodes per direction per cell, node i decomposed lexicographically
    // with the last dimension fastest — the same mapping gkyl_nodal_ops_n2m_interior
    // uses, so the node order matches what quad_nodal_to_modal expects.
    int nidx[GKYL_MAX_CDIM];
    for (int i=0; i<num_basis; ++i) {
      for (int j=0; j<cdim; ++j)
        nidx[j] = (iter.idx[j]-crange->lower[j])*2 + ((i >> (cdim-1-j)) & 1);
      long lin_nidx = gkyl_range_idx(nrange, nidx);

      const double *pn_cov_tangent_basis = gkyl_array_cfetch(cov_tangent_basis_nodal, lin_nidx);
      const double *pn_triad_basis = gkyl_array_cfetch(triad_basis_nodal, lin_nidx);
      const double *pn_triad_basis_gradient = gkyl_array_cfetch(triad_basis_gradient_nodal, lin_nidx);

      compute_h_ij_from_cov_tangent_basis(pn_cov_tangent_basis, h_ij);
      compute_h_ij_inv_from_h_ij(h_ij, h_ij_inv);
      compute_det_h_ij(h_ij, det_h);
      compute_triad_poisson_tensor_ij(h_ij_inv, pn_triad_basis, pn_cov_tangent_basis,
        pn_triad_basis_gradient, &pt_at_nodes[i*num_pt]);
    }

    double *pt_modal = gkyl_array_fetch(conf_poisson_tensor, gkyl_range_idx(crange, iter.idx));
    for (int c=0; c<num_pt; ++c) {
      for (int i=0; i<num_basis; ++i)
        fnodal[i] = pt_at_nodes[i*num_pt + c];
      for (int k=0; k<num_basis; ++k)
        cbasis.quad_nodal_to_modal(fnodal, &pt_modal[c*num_basis], k);
    }
  }
}

// Nonuniform-spacing finite-difference weights for d/dz at the middle node
// (spacing hl to the left neighbor, hr to the right): second-order accurate.
static inline double
fd_central_nonuniform(double fm, double f0, double fp, double hl, double hr)
{
  return (hl*hl*fp - hr*hr*fm + (hr*hr - hl*hl)*f0)/(hl*hr*(hl+hr));
}

// Second-order one-sided difference at the first node of a line: nodes at
// 0, h1, h1+h2 with values f0, f1, f2.
static inline double
fd_onesided(double f0, double f1, double f2, double h1, double h2)
{
  return -(2.0*h1+h2)/(h1*(h1+h2))*f0 + (h1+h2)/(h1*h2)*f1 - h1/(h2*(h1+h2))*f2;
}

void
gkyl_vlasov_triad_geom_from_tangents_interior(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis,
  const struct gkyl_range *nrange,
  const struct gkyl_array *cov_tangent_basis_nodal, const struct gkyl_array *bhat_check_nodal,
  bool exit_at_checks, struct gkyl_array *conf_poisson_tensor)
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

    const double *t1 = &tan[0], *t3 = &tan[6];
    double t3mag = sqrt(t3[0]*t3[0] + t3[1]*t3[1] + t3[2]*t3[2]);
    double bhat[3] = { t3[0]/t3mag, t3[1]/t3mag, t3[2]/t3mag };

    double t1b = t1[0]*bhat[0] + t1[1]*bhat[1] + t1[2]*bhat[2];
    double e1[3] = { t1[0]-t1b*bhat[0], t1[1]-t1b*bhat[1], t1[2]-t1b*bhat[2] };
    double e1mag = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    double t1mag = sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
    if (e1mag < 1e-10*t1mag) {
      fprintf(stderr, "vlasov_triad_geom: first tangent is (near) parallel to bhat at node (%d,%d,%d); b-aligned triad is degenerate\n",
        niter.idx[0], niter.idx[1], niter.idx[2]);
      assert(!exit_at_checks);
    }
    e1[0] /= e1mag; e1[1] /= e1mag; e1[2] /= e1mag;

    // e2 = bhat x e1 makes (e1, e2, bhat) right-handed by construction.
    double e2[3] = { bhat[1]*e1[2]-bhat[2]*e1[1], bhat[2]*e1[0]-bhat[0]*e1[2], bhat[0]*e1[1]-bhat[1]*e1[0] };

    tri[0] = e1[0]; tri[1] = e1[1]; tri[2] = e1[2];
    tri[3] = e2[0]; tri[4] = e2[1]; tri[5] = e2[2];
    tri[6] = bhat[0]; tri[7] = bhat[1]; tri[8] = bhat[2];

    // Cross-check bhat against the independently computed field direction
    // (guards divergence of the geometry provider and this construction).
    if (bhat_check_nodal) {
      const double *bc = gkyl_array_cfetch(bhat_check_nodal, nl);
      double diff = sqrt((bhat[0]-bc[0])*(bhat[0]-bc[0]) + (bhat[1]-bc[1])*(bhat[1]-bc[1]) + (bhat[2]-bc[2])*(bhat[2]-bc[2]));
      if (diff > 1e-10) {
        fprintf(stderr, "vlasov_triad_geom: bhat from tangents differs from provided field direction by %.3e at node (%d,%d,%d)\n",
          diff, niter.idx[0], niter.idx[1], niter.idx[2]);
        assert(!exit_at_checks);
      }
    }
  }

  // Pass 2: gradients d(triad)/dz^j by finite differences across the interior
  // node grid. Nodes sit at the Gauss-Legendre points of each cell, so the
  // spacing alternates: g*dx within a cell, (1-g)*dx across a cell boundary,
  // with g = 1/sqrt(3).
  double glpt = 1.0/sqrt(3.0);
  gkyl_range_iter_init(&niter, nrange);
  while (gkyl_range_iter_next(&niter)) {
    long nl = gkyl_range_idx(nrange, niter.idx);
    double *grad = gkyl_array_fetch(triad_grad_nodal, nl);

    for (int j=0; j<3; ++j) {
      double dx = cgrid->dx[j];
      double gap_in = glpt*dx, gap_across = (1.0-glpt)*dx;
      int nj = niter.idx[j];
      int par = (nj - nrange->lower[j]) & 1; // 0: lower GL node of its cell, 1: upper
      // spacing to the left/right neighbor
      double hl = par == 0 ? gap_across : gap_in;
      double hr = par == 0 ? gap_in : gap_across;

      int idx_n[GKYL_MAX_CDIM];
      gkyl_copy_int_arr(3, niter.idx, idx_n);

      const double *tri_0 = gkyl_array_cfetch(triad_nodal, nl);

      if (nj > nrange->lower[j] && nj < nrange->upper[j]) {
        idx_n[j] = nj-1;
        const double *tri_m = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        idx_n[j] = nj+1;
        const double *tri_p = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        for (int m=0; m<9; ++m)
          grad[j*9 + m] = fd_central_nonuniform(tri_m[m], tri_0[m], tri_p[m], hl, hr);
      }
      else if (nrange->upper[j] - nrange->lower[j] < 2) {
        // Only two nodes along this direction (single cell): two-point difference.
        int no = nj == nrange->lower[j] ? nj+1 : nj-1;
        double sgn = nj == nrange->lower[j] ? 1.0 : -1.0;
        idx_n[j] = no;
        const double *tri_o = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        for (int m=0; m<9; ++m)
          grad[j*9 + m] = sgn*(tri_o[m] - tri_0[m])/gap_in;
      }
      else if (nj == nrange->lower[j]) {
        // One-sided at the lower end: this node is the lower GL node of the
        // first cell, so the next gaps are g*dx then (1-g)*dx.
        idx_n[j] = nj+1;
        const double *tri_1 = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        idx_n[j] = nj+2;
        const double *tri_2 = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        for (int m=0; m<9; ++m)
          grad[j*9 + m] = fd_onesided(tri_0[m], tri_1[m], tri_2[m], gap_in, gap_across);
      }
      else {
        // One-sided at the upper end (mirrored, sign flipped).
        idx_n[j] = nj-1;
        const double *tri_1 = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        idx_n[j] = nj-2;
        const double *tri_2 = gkyl_array_cfetch(triad_nodal, gkyl_range_idx(nrange, idx_n));
        for (int m=0; m<9; ++m)
          grad[j*9 + m] = -fd_onesided(tri_0[m], tri_1[m], tri_2[m], gap_in, gap_across);
      }
    }
  }

  gkyl_vlasov_triad_geom_from_nodal_interior(cgrid, crange, cbasis, pgrid, prange, pbasis,
    nrange, cov_tangent_basis_nodal, triad_nodal, triad_grad_nodal, conf_poisson_tensor);

  gkyl_array_release(triad_nodal);
  gkyl_array_release(triad_grad_nodal);
}

void
gkyl_vlasov_triad_geom_from_vierbein(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis,
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor)
{

  int num_pt_indices[3] = { 1 , 6, 18 };
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // Choose functions
  vierbein_inv_t compute_vierbein_inv_from_vierbein = choose_vierbein_inv_kern(vdim); 
  conf_poisson_tensor_vierbein_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_vierbein_kern(vdim); 

  struct gkyl_array *vierbein_at_nodes;
  struct gkyl_array *vierbein_inverse_at_nodes;
  struct gkyl_array *vierbein_gradient_at_nodes;
  struct gkyl_array *conf_poisson_tensor_at_nodes;

  // Convert nodal to modal using the computed nodal quantities
  gkyl_eval_on_nodes *vierbein_proj;
  gkyl_eval_on_nodes *vierbein_inv_proj;
  gkyl_eval_on_nodes *vierbein_gradient_proj;
  gkyl_eval_on_nodes *conf_poisson_tensor_proj;

  vierbein_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*vdim, NULL, NULL);
  vierbein_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*vdim, NULL, NULL);
  vierbein_gradient_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*vdim*vdim, NULL, NULL);
  conf_poisson_tensor_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim-1], NULL, NULL);

  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int num_basis = cbasis.num_basis;

  // Eval functions
  vierbein_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, num_basis);
  vierbein_gradient_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim*vdim, num_basis);

  // Derived quantities from eval functions
  vierbein_inverse_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim*vdim, num_basis);
  conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim-1], num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  
  // Grab the local node list
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, cbasis.num_basis);
  cbasis.node_list(gkyl_array_fetch(nodes, 0));

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {

      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i),
        cgrid->dx, xc, xmu);
      // Sample the geometry at physical coordinates on non-uniform conf meshes
      // (NULL c2p => identity, i.e. a uniform mesh).
      if (inp_triad_geom.c2p_func)
        inp_triad_geom.c2p_func(xmu, xmu, inp_triad_geom.c2p_func_ctx);

      // Evaluate the functions for basis and their gradients at a nodal point xmu
      inp_triad_geom.eval_vierbein(0.0, xmu, gkyl_array_fetch(vierbein_at_nodes, i), inp_triad_geom.eval_vierbein_ctx);
      inp_triad_geom.eval_vierbein_gradient(0.0, xmu, gkyl_array_fetch(vierbein_gradient_at_nodes, i), inp_triad_geom.eval_vierbein_gradient_ctx);

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
    gkyl_eval_on_nodes_nod2mod(conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes, gkyl_array_fetch(conf_poisson_tensor, lidx));
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
gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_triad_geom, struct gkyl_array *conf_poisson_tensor)
{

  // Choose between constructing from Vierbeins or constructing from basis vectors
  if(inp_triad_geom.use_preset_geom) {
    inp_triad_geom.eval_vierbein = choose_vierbein_kern(inp_triad_geom.triad_preset_geom_type, pgrid->ndim - cgrid->ndim);
    inp_triad_geom.eval_vierbein_gradient = choose_vierbein_gradient_kern(inp_triad_geom.triad_preset_geom_type, pgrid->ndim - cgrid->ndim);

    // ensure non-NULL pointers, requested triad geometry is not currently supported
    assert(inp_triad_geom.eval_vierbein);
    assert(inp_triad_geom.eval_vierbein_gradient);
  }

  if (inp_triad_geom.use_vierbein || inp_triad_geom.use_preset_geom) {
    gkyl_vlasov_triad_geom_from_vierbein(cgrid, crange, cbasis, pgrid, prange, pbasis, 
      inp_triad_geom, conf_poisson_tensor);
  }
  else {
    gkyl_vlasov_triad_geom_from_basis(cgrid, crange, cbasis, pgrid, prange, pbasis, 
      inp_triad_geom, conf_poisson_tensor);
  }
}

evalf_t
gkyl_vlasov_triad_preset_hamil(const int cdim, const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_hamil_kern(preset_geom_type, cdim, vdim);
}

evalf_t
gkyl_vlasov_triad_preset_vierbein(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_vierbein_kern(preset_geom_type, vdim);
}

evalf_t
gkyl_vlasov_triad_preset_vierbein_inv(const int vdim, enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_vierbein_inv_preset_kern(preset_geom_type, vdim);
}
