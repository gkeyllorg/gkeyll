#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_vlasov_velocity_map.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array* a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
#else
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
#endif
  return a;
}

#define VMAX 6.0 // Mapped velocity extent (6 thermal speeds with vth = 1).

// Quadratic mapping from the nonuniform velocity grid regression tests:
// vp = sign(vc)*VMAX*vc^2 on the computational domain [-1, 1].
static void
eval_quad_vmap(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = vc[0] < 0.0 ? -VMAX*vc[0]*vc[0] : VMAX*vc[0]*vc[0];
}

// Project a Maxwellian with constant (n=1, V_drift=0, T/m=1) on a quadratically
// mapped velocity grid with the LTE projection (whose quadrature evaluates the
// mapped velocity at the same Gauss-Legendre nodes the Jacobian is known at and
// multiplies by the Jacobian, producing Jf), then divide the Jacobian out
// nodally. Because the tensor-basis Gauss-Legendre projection is a nodal
// interpolation and the division happens at those same nodes, the result must
// reproduce the Maxwellian evaluated at the mapped quadrature points up to one
// global constant (the projection's density-rescale factor).
static void
test_divide_jacobvel(int vdim, bool use_gpu)
{
  int poly_order = 2;
  int cdim = 1, pdim = cdim+vdim;

  double lower[GKYL_MAX_DIM] = { 0.0, -1.0, -1.0, -1.0 };
  double upper[GKYL_MAX_DIM] = { 1.0, 1.0, 1.0, 1.0 };
  // 32 cells per velocity direction keeps the within-cell dynamic range of the
  // Maxwellian in the outermost mapped cells small (~8.6 per direction), so
  // the nodal checks below can use tight relative tolerances at every node.
  int cells[GKYL_MAX_DIM] = { 1, 32, 32, 32 };

  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, &lower[cdim], &upper[cdim], &cells[cdim]);

  // Basis functions: tensor bases, as required by the nodal Jacobian division.
  struct gkyl_basis pbasis, cbasis, vbasis;
  gkyl_cart_modal_tensor(&pbasis, pdim, poly_order);
  gkyl_cart_modal_tensor(&cbasis, cdim, poly_order);
  gkyl_cart_modal_tensor(&vbasis, vdim, poly_order);

  int conf_ghost[] = { 1 };
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = { 0, 0, 0 };
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[GKYL_MAX_DIM] = { 1, 0, 0, 0 };
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  // Velocity map with the quadratic mapping in every direction.
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    inp_vmap[d].eval_vmap = eval_quad_vmap;
  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vel_grid,
    &vel_local, &vbasis, inp_vmap, false, use_gpu);

  // Velocity-space Hamiltonian (H = v^2/2 on the mapped grid), needed by the
  // moment computation inside the LTE projection's density rescale.
  struct gkyl_array *hamil = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &vbasis, &vel_local,
    GKYL_MODEL_DEFAULT, vvm, hamil, gamma_inv, use_gpu);

  // Constant LTE moments (n, V_drift, T/m) = (1, 0, 1).
  struct gkyl_array *moms_ho = gkyl_array_new(GKYL_DOUBLE, (vdim+2)*cbasis.num_basis, conf_local_ext.volume);
  gkyl_array_clear(moms_ho, 0.0);
  double *moms_d = gkyl_array_fetch(moms_ho, gkyl_range_idx(&conf_local, (int[]) { 1 }));
  moms_d[0] = sqrt(2.0); // n = 1 (1D p2 cell-constant coefficient normalization).
  moms_d[(vdim+1)*cbasis.num_basis] = sqrt(2.0); // T/m = 1.
  struct gkyl_array *moms = use_gpu ? mkarr(use_gpu, moms_ho->ncomp, moms_ho->size) : gkyl_array_acquire(moms_ho);
  gkyl_array_copy(moms, moms_ho);

  // Project the LTE (Maxwellian) distribution; output is Jf on the mapped grid.
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &cbasis,
    .vel_basis = &vbasis,
    .phase_basis = &pbasis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .vel_map = vvm,
    .hamil_range = &vel_local,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = use_gpu,
  };
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_proj);

  struct gkyl_array *f_lte = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &phase_local, &conf_local, moms, f_lte);

  // Divide out the velocity-space Jacobian nodally, and rescale back for the
  // round-trip check.
  struct gkyl_array *f_no_J = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *Jf2 = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  gkyl_vlasov_velocity_map_divide_jacobvel(vvm, &cbasis, &pbasis, &phase_local, f_lte, f_no_J);
  gkyl_vlasov_velocity_map_rescale_jacobvel(vvm, &cbasis, &pbasis, &phase_local, f_no_J, Jf2);

  // Host copies for checking.
  struct gkyl_array *f_lte_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *f_no_J_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *Jf2_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  gkyl_array_copy(f_lte_ho, f_lte);
  gkyl_array_copy(f_no_J_ho, f_no_J);
  gkyl_array_copy(Jf2_ho, Jf2);

  // 1D cubic basis for evaluating the per-direction map at quadrature points.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // Quadrature ranges matching the projection's node ordering.
  int pqshape[GKYL_MAX_DIM], vqshape[GKYL_MAX_DIM];
  for (int d=0; d<pdim; ++d) pqshape[d] = poly_order+1;
  for (int d=0; d<vdim; ++d) vqshape[d] = poly_order+1;
  struct gkyl_range phase_qrange, vel_qrange;
  gkyl_range_init_from_shape(&phase_qrange, pdim, pqshape);
  gkyl_range_init_from_shape(&vel_qrange, vdim, vqshape);
  const double *ord = gkyl_gauss_ordinates[poly_order+1];

  double maxwell_norm = 1.0/sqrt(pow(2.0*GKYL_PI, vdim));

  // First pass: find the global maximum nodal value of f_no_J and the
  // global ratio r0 = f_no_J/M there (the projection's density-rescale factor).
  double fN_max = 0.0, r0 = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&phase_local, iter.idx);
    long vloc = gkyl_range_idx(&vel_local, &iter.idx[cdim]);
    const double *fnoJ_c = gkyl_array_cfetch(f_no_J_ho, ploc);
    const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);

    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &phase_qrange);
    while (gkyl_range_iter_next(&qiter)) {
      double eta[GKYL_MAX_DIM];
      for (int d=0; d<pdim; ++d) eta[d] = ord[qiter.idx[d]];
      double fN = pbasis.eval_expand(eta, fnoJ_c);
      if (fN > fN_max) {
        double vsq = 0.0;
        for (int d=0; d<vdim; ++d) {
          double zv[] = { eta[cdim+d] };
          double vp = b1.eval_expand(zv, &vmap_c[4*d]);
          vsq += vp*vp;
        }
        fN_max = fN;
        r0 = fN/(maxwell_norm*exp(-0.5*vsq));
      }
    }
  }

  // The density rescale should leave the projection close to the target
  // density. The deviation of r0 from 1 is set by the Maxwellian tail mass
  // truncated beyond +/- 6 vth plus the Gauss-Legendre quadrature error of the
  // density moment on the mapped grid, both per velocity dimension: measured
  // r0-1 = -vdim*1.26e-8 at this resolution.
  TEST_CHECK( fabs(r0 - 1.0) < 1e-7 );
  TEST_MSG("vdim=%d: global ratio r0 = %.15e", vdim, r0);

  // Second pass: per-cell nodal checks.
  gkyl_range_iter_init(&iter, &phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&phase_local, iter.idx);
    long vloc = gkyl_range_idx(&vel_local, &iter.idx[cdim]);
    const double *Jf_c = gkyl_array_cfetch(f_lte_ho, ploc);
    const double *Jf2_c = gkyl_array_cfetch(Jf2_ho, ploc);
    const double *fnoJ_c = gkyl_array_cfetch(f_no_J_ho, ploc);
    const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
    const double *jacob_gauss_c = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, vloc);

    // (a) Round trip: rescale(divide(Jf)) reproduces Jf to machine precision.
    double Jf_c_max = 0.0;
    for (int k=0; k<pbasis.num_basis; ++k)
      Jf_c_max = fmax(Jf_c_max, fabs(Jf_c[k]));
    for (int k=0; k<pbasis.num_basis; ++k) {
      TEST_CHECK( fabs(Jf2_c[k] - Jf_c[k]) < 1e-14*Jf_c_max );
      TEST_MSG("vdim=%d cell=(%d,%d): round trip coeff %d: %.15e vs %.15e",
        vdim, iter.idx[0], iter.idx[1], k, Jf2_c[k], Jf_c[k]);
    }

    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &phase_qrange);
    while (gkyl_range_iter_next(&qiter)) {
      double eta[GKYL_MAX_DIM];
      for (int d=0; d<pdim; ++d) eta[d] = ord[qiter.idx[d]];

      // Mapped velocity coordinates at this quadrature node.
      double vsq = 0.0;
      for (int d=0; d<vdim; ++d) {
        double zv[] = { eta[cdim+d] };
        double vp = b1.eval_expand(zv, &vmap_c[4*d]);
        vsq += vp*vp;
      }
      double M = maxwell_norm*exp(-0.5*vsq);
      double fN = pbasis.eval_expand(eta, fnoJ_c);
      double JfN = pbasis.eval_expand(eta, Jf_c);

      // (b) The divided distribution at the quadrature nodes is exactly the
      // Maxwellian evaluated at the mapped nodes, up to the global constant
      // r0, since projection and division are both nodal at these points.
      // The roundoff of the nodal <-> modal transforms is amplified by the
      // within-cell dynamic range of the Maxwellian (largest in the outermost
      // mapped corner cells, ~640 at this resolution; measured worst-case
      // deviation 6e-12 in 1x3v), so a strict relative check holds at every
      // node, including the deep tail.
      TEST_CHECK( gkyl_compare_double(fN/(r0*M), 1.0, 1e-10) );
      TEST_MSG("vdim=%d cell=(%d,%d) node=(%d): f_no_J %.15e vs r0*M %.15e",
        vdim, iter.idx[0], iter.idx[1], qiter.idx[0], fN, r0*M);

      // (c) Jf and f_no_J at the nodes differ by exactly the total Jacobian at
      // the corresponding velocity quadrature node (ties the projection, the
      // division kernels, and the map to the same node ordering).
      long vqidx = gkyl_range_idx(&vel_qrange, &qiter.idx[cdim]);
      double J = jacob_gauss_c[vqidx];
      TEST_CHECK( gkyl_compare_double(JfN/(J*fN), 1.0, 1e-10) );
      TEST_MSG("vdim=%d cell=(%d,%d): Jf %.15e vs J*f %.15e (J=%.15e)",
        vdim, iter.idx[0], iter.idx[1], JfN, J*fN, J);
    }
  }

  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_array_release(f_lte);
  gkyl_array_release(f_no_J);
  gkyl_array_release(Jf2);
  gkyl_array_release(f_lte_ho);
  gkyl_array_release(f_no_J_ho);
  gkyl_array_release(Jf2_ho);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_vlasov_velocity_map_release(vvm);
}

// Independent Gauss-Legendre quadrature projection of the constant-moment
// (n=1, V_drift=0, T/m=1) Maxwellian on one cell: out[k] = sum_q w_q
// psi_k(eta_q) M(vmap(eta_q)). For Serendipity bases this is the same
// operation the LTE projection performs (there are more quadrature points
// than basis functions, so the projection is not a nodal interpolation and
// the modal coefficients, not nodal values, are the quantities to compare).
static void
quadproj_maxwellian(const struct gkyl_basis *pbasis, const struct gkyl_basis *b1,
  const struct gkyl_range *phase_qrange, int poly_order, int cdim, int vdim,
  const double *vmap_c, double maxwell_norm, double *out)
{
  const double *ord = gkyl_gauss_ordinates[poly_order+1];
  const double *wgt = gkyl_gauss_weights[poly_order+1];
  double bvals[160];

  for (int k=0; k<pbasis->num_basis; ++k) out[k] = 0.0;

  struct gkyl_range_iter qiter;
  gkyl_range_iter_init(&qiter, phase_qrange);
  while (gkyl_range_iter_next(&qiter)) {
    double eta[GKYL_MAX_DIM], wtot = 1.0;
    for (int d=0; d<cdim+vdim; ++d) {
      eta[d] = ord[qiter.idx[d]];
      wtot *= wgt[qiter.idx[d]];
    }
    double vsq = 0.0;
    for (int d=0; d<vdim; ++d) {
      double zv[] = { eta[cdim+d] };
      double vp = b1->eval_expand(zv, &vmap_c[4*d]);
      vsq += vp*vp;
    }
    double M = maxwell_norm*exp(-0.5*vsq);
    pbasis->eval(eta, bvals);
    for (int k=0; k<pbasis->num_basis; ++k)
      out[k] += wtot*bvals[k]*M;
  }
}

// Serendipity (C^0 linear map) version: the velocity-space Jacobian is
// piecewise constant, so dividing it out of Jf is exact coefficient-by-
// coefficient. Project the Maxwellian with the LTE projection (producing Jf
// on the mapped grid), divide the Jacobian out, and check (a) the round
// trip, (b) that the division is exactly by the cell's constant, and (c)
// that the divided modal coefficients reproduce an independently computed
// Gauss-Legendre quadrature projection of the Maxwellian, up to the
// projection's global density-rescale factor r0.
static void
test_divide_jacobvel_ser(int vdim, int poly_order, bool use_gpu)
{
  int cdim = 1, pdim = cdim+vdim;

  double lower[GKYL_MAX_DIM] = { 0.0, -1.0, -1.0, -1.0 };
  double upper[GKYL_MAX_DIM] = { 1.0, 1.0, 1.0, 1.0 };
  int cells[GKYL_MAX_DIM] = { 1, 32, 32, 32 };

  struct gkyl_rect_grid phase_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, &lower[cdim], &upper[cdim], &cells[cdim]);

  struct gkyl_basis pbasis, cbasis, vbasis;
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&vbasis, vdim, poly_order);

  int conf_ghost[] = { 1 };
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[] = { 0, 0, 0 };
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int phase_ghost[GKYL_MAX_DIM] = { 1, 0, 0, 0 };
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&phase_grid, phase_ghost, &phase_local_ext, &phase_local);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    inp_vmap[d].eval_vmap = eval_quad_vmap;
  struct gkyl_vlasov_velocity_map *vvm = gkyl_vlasov_velocity_map_new(&vel_grid,
    &vel_local, &vbasis, inp_vmap, false, use_gpu);

  TEST_CHECK( vvm->rep == GKYL_VLASOV_VMAP_C0_LINEAR );

  struct gkyl_array *hamil = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  gkyl_dg_vlasov_calc_hamil(&vel_grid, &vbasis, &vel_local,
    GKYL_MODEL_DEFAULT, vvm, hamil, gamma_inv, use_gpu);

  struct gkyl_array *moms_ho = gkyl_array_new(GKYL_DOUBLE, (vdim+2)*cbasis.num_basis, conf_local_ext.volume);
  gkyl_array_clear(moms_ho, 0.0);
  double *moms_d = gkyl_array_fetch(moms_ho, gkyl_range_idx(&conf_local, (int[]) { 1 }));
  moms_d[0] = sqrt(2.0); // n = 1.
  moms_d[(vdim+1)*cbasis.num_basis] = sqrt(2.0); // T/m = 1.
  struct gkyl_array *moms = use_gpu ? mkarr(use_gpu, moms_ho->ncomp, moms_ho->size) : gkyl_array_acquire(moms_ho);
  gkyl_array_copy(moms, moms_ho);

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &cbasis,
    .vel_basis = &vbasis,
    .phase_basis = &pbasis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_local_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .vel_map = vvm,
    .hamil_range = &vel_local,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .use_gpu = use_gpu,
  };
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_proj);

  struct gkyl_array *f_lte = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &phase_local, &conf_local, moms, f_lte);

  struct gkyl_array *f_no_J = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *Jf2 = mkarr(use_gpu, pbasis.num_basis, phase_local_ext.volume);
  gkyl_vlasov_velocity_map_divide_jacobvel(vvm, &cbasis, &pbasis, &phase_local, f_lte, f_no_J);
  gkyl_vlasov_velocity_map_rescale_jacobvel(vvm, &cbasis, &pbasis, &phase_local, f_no_J, Jf2);

  struct gkyl_array *f_lte_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *f_no_J_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  struct gkyl_array *Jf2_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  gkyl_array_copy(f_lte_ho, f_lte);
  gkyl_array_copy(f_no_J_ho, f_no_J);
  gkyl_array_copy(Jf2_ho, Jf2);

  // 1D cubic basis for evaluating the per-direction (degenerate cubic) map.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  int pqshape[GKYL_MAX_DIM];
  for (int d=0; d<pdim; ++d) pqshape[d] = poly_order+1;
  struct gkyl_range phase_qrange;
  gkyl_range_init_from_shape(&phase_qrange, pdim, pqshape);

  double maxwell_norm = 1.0/sqrt(pow(2.0*GKYL_PI, vdim));
  double f_expected[160];

  // First pass: locate the cell with the largest |c0| of f_no_J and measure
  // the global density-rescale factor r0 there from the modal coefficient
  // ratio.
  double c0_max = 0.0, r0 = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&phase_local, iter.idx);
    const double *fnoJ_c = gkyl_array_cfetch(f_no_J_ho, ploc);
    if (fabs(fnoJ_c[0]) > c0_max) {
      long vloc = gkyl_range_idx(&vel_local, &iter.idx[cdim]);
      const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
      quadproj_maxwellian(&pbasis, &b1, &phase_qrange, poly_order, cdim, vdim,
        vmap_c, maxwell_norm, f_expected);
      c0_max = fabs(fnoJ_c[0]);
      r0 = fnoJ_c[0]/f_expected[0];
    }
  }

  // The density rescale should leave the projection close to the target
  // density; deviations come from the truncated Maxwellian tail and the
  // quadrature error of the density moment on the mapped grid.
  TEST_CHECK( fabs(r0 - 1.0) < 1e-3 );
  TEST_MSG("vdim=%d poly_order=%d: global ratio r0 = %.15e", vdim, poly_order, r0);

  // Second pass: per-cell checks.
  gkyl_range_iter_init(&iter, &phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&phase_local, iter.idx);
    long vloc = gkyl_range_idx(&vel_local, &iter.idx[cdim]);
    const double *Jf_c = gkyl_array_cfetch(f_lte_ho, ploc);
    const double *Jf2_c = gkyl_array_cfetch(Jf2_ho, ploc);
    const double *fnoJ_c = gkyl_array_cfetch(f_no_J_ho, ploc);
    const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
    const double *jacob_gauss_c = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, vloc);
    double jac = jacob_gauss_c[0];

    double Jf_c_max = 0.0, fnoJ_c_max = 0.0;
    for (int k=0; k<pbasis.num_basis; ++k) {
      Jf_c_max = fmax(Jf_c_max, fabs(Jf_c[k]));
      fnoJ_c_max = fmax(fnoJ_c_max, fabs(fnoJ_c[k]));
    }

    quadproj_maxwellian(&pbasis, &b1, &phase_qrange, poly_order, cdim, vdim,
      vmap_c, maxwell_norm, f_expected);

    for (int k=0; k<pbasis.num_basis; ++k) {
      // (a) Round trip: rescale(divide(Jf)) reproduces Jf to machine precision.
      TEST_CHECK( fabs(Jf2_c[k] - Jf_c[k]) < 1e-14*Jf_c_max );
      TEST_MSG("vdim=%d p=%d cell=(%d,%d): round trip coeff %d: %.15e vs %.15e",
        vdim, poly_order, iter.idx[0], iter.idx[1], k, Jf2_c[k], Jf_c[k]);

      // (b) The division is exactly by the cell's constant Jacobian,
      // coefficient-by-coefficient.
      TEST_CHECK( fabs(jac*fnoJ_c[k] - Jf_c[k]) < 1e-14*Jf_c_max );
      TEST_MSG("vdim=%d p=%d cell=(%d,%d): J*f coeff %d: %.15e vs %.15e (J=%.15e)",
        vdim, poly_order, iter.idx[0], iter.idx[1], k, jac*fnoJ_c[k], Jf_c[k], jac);

      // (c) The divided modal coefficients reproduce the independent
      // quadrature projection of the Maxwellian, up to the global r0.
      TEST_CHECK( fabs(fnoJ_c[k] - r0*f_expected[k]) < 1e-11*fnoJ_c_max );
      TEST_MSG("vdim=%d p=%d cell=(%d,%d): modal coeff %d: %.15e vs r0*expected %.15e",
        vdim, poly_order, iter.idx[0], iter.idx[1], k, fnoJ_c[k], r0*f_expected[k]);
    }
  }

  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(moms_ho);
  gkyl_array_release(moms);
  gkyl_array_release(f_lte);
  gkyl_array_release(f_no_J);
  gkyl_array_release(Jf2);
  gkyl_array_release(f_lte_ho);
  gkyl_array_release(f_no_J_ho);
  gkyl_array_release(Jf2_ho);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_vlasov_velocity_map_release(vvm);
}

static void test_divide_jacobvel_1x1v(void) { test_divide_jacobvel(1, false); }
static void test_divide_jacobvel_1x2v(void) { test_divide_jacobvel(2, false); }
static void test_divide_jacobvel_1x3v(void) { test_divide_jacobvel(3, false); }

static void test_divide_jacobvel_ser_1x1v_p1(void) { test_divide_jacobvel_ser(1, 1, false); }
static void test_divide_jacobvel_ser_1x1v_p2(void) { test_divide_jacobvel_ser(1, 2, false); }
static void test_divide_jacobvel_ser_1x1v_p3(void) { test_divide_jacobvel_ser(1, 3, false); }
static void test_divide_jacobvel_ser_1x2v_p1(void) { test_divide_jacobvel_ser(2, 1, false); }
static void test_divide_jacobvel_ser_1x2v_p2(void) { test_divide_jacobvel_ser(2, 2, false); }
static void test_divide_jacobvel_ser_1x3v_p1(void) { test_divide_jacobvel_ser(3, 1, false); }
static void test_divide_jacobvel_ser_1x3v_p2(void) { test_divide_jacobvel_ser(3, 2, false); }

#ifdef GKYL_HAVE_CUDA
static void test_divide_jacobvel_1x1v_gpu(void) { test_divide_jacobvel(1, true); }
static void test_divide_jacobvel_1x2v_gpu(void) { test_divide_jacobvel(2, true); }
static void test_divide_jacobvel_1x3v_gpu(void) { test_divide_jacobvel(3, true); }
static void test_divide_jacobvel_ser_1x1v_p1_gpu(void) { test_divide_jacobvel_ser(1, 1, true); }
static void test_divide_jacobvel_ser_1x2v_p2_gpu(void) { test_divide_jacobvel_ser(2, 2, true); }
static void test_divide_jacobvel_ser_1x3v_p1_gpu(void) { test_divide_jacobvel_ser(3, 1, true); }
#endif

TEST_LIST = {
  { "divide_jacobvel_1x1v", test_divide_jacobvel_1x1v },
  { "divide_jacobvel_1x2v", test_divide_jacobvel_1x2v },
  { "divide_jacobvel_1x3v", test_divide_jacobvel_1x3v },
  { "divide_jacobvel_ser_1x1v_p1", test_divide_jacobvel_ser_1x1v_p1 },
  { "divide_jacobvel_ser_1x1v_p2", test_divide_jacobvel_ser_1x1v_p2 },
  { "divide_jacobvel_ser_1x1v_p3", test_divide_jacobvel_ser_1x1v_p3 },
  { "divide_jacobvel_ser_1x2v_p1", test_divide_jacobvel_ser_1x2v_p1 },
  { "divide_jacobvel_ser_1x2v_p2", test_divide_jacobvel_ser_1x2v_p2 },
  { "divide_jacobvel_ser_1x3v_p1", test_divide_jacobvel_ser_1x3v_p1 },
  { "divide_jacobvel_ser_1x3v_p2", test_divide_jacobvel_ser_1x3v_p2 },
#ifdef GKYL_HAVE_CUDA
  { "divide_jacobvel_1x1v_gpu", test_divide_jacobvel_1x1v_gpu },
  { "divide_jacobvel_1x2v_gpu", test_divide_jacobvel_1x2v_gpu },
  { "divide_jacobvel_1x3v_gpu", test_divide_jacobvel_1x3v_gpu },
  { "divide_jacobvel_ser_1x1v_p1_gpu", test_divide_jacobvel_ser_1x1v_p1_gpu },
  { "divide_jacobvel_ser_1x2v_p2_gpu", test_divide_jacobvel_ser_1x2v_p2_gpu },
  { "divide_jacobvel_ser_1x3v_p1_gpu", test_divide_jacobvel_ser_1x3v_p1_gpu },
#endif
  { NULL, NULL },
};
