// Unit tests of the nodal triad kernels (vlasov/ker/{vlasov,flux_vlasov}/*triad_nodal*): the
// tensor hybrid p=1 basis (p=1 tensor in configuration space x p=2 tensor in velocity space)
// with the Poisson tensor stored at the configuration-space Gauss-Legendre nodes, for the
// sparse and dense velocity Hamiltonians and the phase-space Hamiltonian (GR triads).
//
// The kernels are checked against their definition, evaluated independently of the kernel
// generator with the basis functions' own eval/eval_grad_expand and ordinary Gauss quadrature:
//  - volume term: sum over the conf-vel pairs (x_d, v_b) of
//      (2/dx_d)(1/J_xd)(2/dv_b) int [d_xd(psi) A_db - d_vb(psi) B_db] F_b,
//    A_db = I[Pi_db d_vb(H)], B_db = I[Pi_db d_xd(H)], and over the vel-vel pairs (v_a, v_b) of
//      (2/dv_a)(2/dv_b) int [d_va(psi) I[w_ab d_vb(H)] - d_vb(psi) I[w_ab d_va(H)]] F_ab,
//    w_ab = sum_j Omega_ab^j v_j, where I[.] is the p=1 interpolant in configuration space of
//    its argument's values at the 2^cdim interior Gauss-Legendre nodes (the exact Galerkin
//    integral of the interpolated characteristic, which the kernels evaluate with 2-point
//    quadrature), and F_b, F_ab are f divided by the velocity-map Jacobian at the 3
//    Gauss-Legendre points along v_b (and v_a), re-interpolated;
//  - surface forces: the characteristic at the surface nodes, alpha^{x_d} = sum_b Pi_db
//    (2/dv_b) d_vb(H)/J_vb on the configuration faces, alpha^{v_a} = -sum_c Pi_ca (2/dx_c)
//    (1/J_xc) d_xc(H) + sum_b Pi^{ab} (2/dv_b) d_vb(H)/J_vb on the lower velocity face, and
//    the sum over terms of the factors of their _shared (sum-factorized) forms;
//  - the representations against each other: dense == sparse for a separable H and
//    phase == dense for an H independent of the configuration-space coordinates;
//  - the canonical-PB phase-space Hamiltonian volume kernels of the tensor p=1 basis
//    (vlasov_hamil_phase_vol_*_tensor_p1, same generator): the volume term above with the
//    canonical bracket, Pi_db = delta_db, no Jacobians and no velocity-velocity block;
//  - the magnetic-field force volume kernels (vlasov_B<d>_hamil_<rep>_vol_*_tensor_p1):
//    (2/dv_a)(2/dv_b) int [d_va(psi) I[B d_vb(H)] - d_vb(psi) I[B d_va(H)]] F for B along d,
//    (a, b) the transverse velocity directions, F = f times the velocity-map Jacobian of the
//    third direction (nodally for a velocity Hamiltonian, the cell-constant J for a
//    phase-space Hamiltonian; F = f in 2V);
//  - the velocity-Jacobian division and rescaling kernels (vlasov_{divide,rescale}_Jv_*_tensor_p1):
//    at any configuration-space point and every velocity Gauss-Legendre node, J times the
//    output equals the input (output equals J times the input);
//  - the moment kernels (mom_vlasov_hamil_<rep>_{M1i,five_moments,int_five_moments}_*_tensor_p1):
//    M0 = int psi_k Jf, M1i = (2/dv_i) int psi_k d_vi(H) F_i and M2 = int psi_k H Jf over the
//    velocity cell (and the configuration-space cell for the integrated moments, psi_k = 1),
//    F_i = Jf divided by the velocity-map Jacobian at the 3 Gauss-Legendre points along v_i;
//  - the velocity-space Hamiltonian streaming kernels (vlasov_hamil_<rep>_{vol,surf<d>,boundary_surf<d>}):
//    sum_d (2/dx_d)(2/dv_d)/J_xd int d_xd(psi) d_vd(H) F_d over the cell, and on the faces x_d = -1, +1
//    (2/dx_d)(2/dv_d) int psi d_vd(H) U with U the upwinded (on the sign of w_vd) F_d/J_xd of the cell
//    on the left or right of the face; the force kernels vlasov_{E,phi,rad}_vol: (2/dv_d) int d_vd(psi) a_d F_d
//    with a_d = E_d(x), -(2/dx_d)/J_xd d_xd(phi) and the radiation drag rad_d(v).
// Velocity maps are random cubics, Poisson tensors and Hamiltonians random. With CUDA, the
// triad updaters (conf/vel surface fluxes and hyper_dg) are also run on CPU and GPU and
// compared, which also exercises the kernel selection of the tensor p=1 triads.
#include <acutest.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_conf_flux_surf.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_vlasov_triad_nodal_kernels.h>
#include <gkyl_vlasov_velocity_map.h>

typedef void (*tn_vol_t)(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *pt, const double *hamil,
  const double *f, const double *f_no_J, double *out);
typedef void (*tn_conf_t)(const double *w, const double *dxv, const int edge, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *pt, const double *hamil,
  double *alpha_quad);
typedef int (*tn_conf_sh_t)(int tid, int off, const double *w, const double *dxv, const int edge,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *pt,
  const double *hamil, double *O, double *I);
typedef void (*tn_vel_t)(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *pt, const double *hamil,
  double *alpha_quad);
typedef int (*tn_vel_sh_t)(int tid, int off, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *pt,
  const double *hamil, double *O, double *I);

// Kernels of one (Hamiltonian, cdim, vdim): volume, and the surface force producers of each
// direction in their standard and ho_ forms ([ho][dir]), plain and _shared.
struct tn_kernels {
  tn_vol_t vol;
  tn_conf_t conf[2][3];
  tn_conf_sh_t conf_sh[2][3];
  tn_vel_t vel[2][3];
  tn_vel_sh_t vel_sh[2][3];
};

#define TN_VOL(r, cv) vlasov_triad_nodal_##r##_vol_##cv##_tensor_p1
#define TN_K(r, ho, d, cv) triad_nodal_hamil_##r##_##ho##alpha_quad_##d##_##cv##_tensor_p1
#define TN_KS(r, ho, d, cv) triad_nodal_hamil_##r##_##ho##alpha_quad_##d##_##cv##_tensor_p1_shared
#define TN_C1(r, ho, cv) { TN_K(r, ho, x, cv), NULL, NULL }
#define TN_C2(r, ho, cv) { TN_K(r, ho, x, cv), TN_K(r, ho, y, cv), NULL }
#define TN_C3(r, ho, cv) { TN_K(r, ho, x, cv), TN_K(r, ho, y, cv), TN_K(r, ho, z, cv) }
#define TN_C1S(r, ho, cv) { TN_KS(r, ho, x, cv), NULL, NULL }
#define TN_C2S(r, ho, cv) { TN_KS(r, ho, x, cv), TN_KS(r, ho, y, cv), NULL }
#define TN_C3S(r, ho, cv) { TN_KS(r, ho, x, cv), TN_KS(r, ho, y, cv), TN_KS(r, ho, z, cv) }
#define TN_V1(r, ho, cv) { TN_K(r, ho, vx, cv), NULL, NULL }
#define TN_V2(r, ho, cv) { TN_K(r, ho, vx, cv), TN_K(r, ho, vy, cv), NULL }
#define TN_V3(r, ho, cv) { TN_K(r, ho, vx, cv), TN_K(r, ho, vy, cv), TN_K(r, ho, vz, cv) }
#define TN_V1S(r, ho, cv) { TN_KS(r, ho, vx, cv), NULL, NULL }
#define TN_V2S(r, ho, cv) { TN_KS(r, ho, vx, cv), TN_KS(r, ho, vy, cv), NULL }
#define TN_V3S(r, ho, cv) { TN_KS(r, ho, vx, cv), TN_KS(r, ho, vy, cv), TN_KS(r, ho, vz, cv) }
#define TN_SET(r, cv, CN, VN) { TN_VOL(r, cv), { CN(r, , cv), CN(r, ho_, cv) }, \
  { CN##S(r, , cv), CN##S(r, ho_, cv) }, { VN(r, , cv), VN(r, ho_, cv) }, \
  { VN##S(r, , cv), VN##S(r, ho_, cv) } }

enum tn_rep { TN_SPARSE, TN_DENSE, TN_PHASE };
static const char *tn_rep_name[3] = { "vel_sparse", "vel_dense", "phase" };
static const int tn_cv[6][2] = { {1, 1}, {1, 2}, {2, 2}, {1, 3}, {2, 3}, {3, 3} };

// [Hamiltonian][cv]; sparse is not generated in 1V (it aliases dense there).
static const struct tn_kernels tn_tab[3][6] = {
  { {0}, TN_SET(vel_sparse, 1x2v, TN_C1, TN_V2), TN_SET(vel_sparse, 2x2v, TN_C2, TN_V2),
    TN_SET(vel_sparse, 1x3v, TN_C1, TN_V3), TN_SET(vel_sparse, 2x3v, TN_C2, TN_V3),
    TN_SET(vel_sparse, 3x3v, TN_C3, TN_V3) },
  { TN_SET(vel_dense, 1x1v, TN_C1, TN_V1), TN_SET(vel_dense, 1x2v, TN_C1, TN_V2),
    TN_SET(vel_dense, 2x2v, TN_C2, TN_V2), TN_SET(vel_dense, 1x3v, TN_C1, TN_V3),
    TN_SET(vel_dense, 2x3v, TN_C2, TN_V3), TN_SET(vel_dense, 3x3v, TN_C3, TN_V3) },
  { TN_SET(phase, 1x1v, TN_C1, TN_V1), TN_SET(phase, 1x2v, TN_C1, TN_V2),
    TN_SET(phase, 2x2v, TN_C2, TN_V2), TN_SET(phase, 1x3v, TN_C1, TN_V3),
    TN_SET(phase, 2x3v, TN_C2, TN_V3), TN_SET(phase, 3x3v, TN_C3, TN_V3) },
};

// Number of Poisson-tensor components (conf-vel block vdim x vdim, then the vel-vel block).
static const int tn_npt[4] = { 0, 1, 6, 18 };

// Gauss-Legendre nodes and weights on [-1,1], ascending.
static const double tn_gl_x[6][5] = {
  {0}, {0}, { -0.57735026918962576, 0.57735026918962576 },
  { -0.77459666924148338, 0.0, 0.77459666924148338 },
  { -0.86113631159405258, -0.33998104358485626, 0.33998104358485626, 0.86113631159405258 },
  { -0.90617984593866399, -0.53846931010568309, 0.0, 0.53846931010568309, 0.90617984593866399 },
};
static const double tn_gl_w[6][5] = {
  {0}, {0}, { 1.0, 1.0 }, { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
  { 0.34785484513745386, 0.65214515486254614, 0.65214515486254614, 0.34785484513745386 },
  { 0.23692688505618909, 0.47862867049936647, 0.56888888888888889, 0.47862867049936647,
    0.23692688505618909 },
};

static double
tn_rnd(double lo, double hi)
{
  return lo + (hi - lo) * rand() / (double)RAND_MAX;
}

// One test case: inputs of a cell for a given Hamiltonian, cdim, vdim.
struct tn_case {
  enum tn_rep rep;
  bool canonical; // canonical bracket (canonical-PB phase-space Hamiltonian volume)
  int cdim, vdim;
  struct gkyl_basis pbasis, vbasis, cbasis; // phase (tensor hybrid), velocity (tensor p=2), conf
  double w[6], dxv[6], vmap[12], jacob_pos[6], jacob_vel[9], jvs[2][12]; // jvs[ho]
  double pt_vol[8 * 18], pt_surf[3 * 2 * 4 * 3], hamil[216], f[216];
  double qmem[6 * 8]; // E and B, configuration-space expansions
};

// Physical velocity of the cubic map of one direction (1D p=3 modal expansion) and its
// derivative with respect to the logical coordinate.
static double
tn_vmap(const double *vm, double eta)
{
  return vm[0] / sqrt(2.0) + vm[1] * sqrt(1.5) * eta
    + vm[2] * sqrt(2.5) * (3.0 * eta * eta - 1.0) / 2.0
    + vm[3] * sqrt(3.5) * (5.0 * eta * eta * eta - 3.0 * eta) / 2.0;
}

static double
tn_vmap_deriv(const double *vm, double eta)
{
  return vm[1] * sqrt(1.5) + vm[2] * sqrt(2.5) * 3.0 * eta
    + vm[3] * sqrt(3.5) * (15.0 * eta * eta - 3.0) / 2.0;
}

// Coordinates of node q of the gaussOrd(2, ndim) set (first coordinate slowest).
static void
tn_node(int ndim, int q, double *x)
{
  for (int d = 0; d < ndim; ++d)
    x[d] = tn_gl_x[2][(q >> (ndim - 1 - d)) & 1];
}

// Whether velocity basis function m depends on at most one velocity coordinate.
static bool
tn_is_sparse_mode(const struct gkyl_basis *vb, int m)
{
  double z[3] = { 0.3, -0.7, 0.45 }, bp[27], bm[27];
  int nvar = 0;
  for (int d = 0; d < vb->ndim; ++d) {
    double zp[3], zm[3];
    memcpy(zp, z, sizeof z); memcpy(zm, z, sizeof z);
    zp[d] += 0.5; zm[d] -= 0.5;
    vb->eval(zp, bp); vb->eval(zm, bm);
    if (fabs(bp[m] - bm[m]) > 1e-12) nvar += 1;
  }
  return nvar <= 1;
}

static void
tn_case_init(struct tn_case *cs, enum tn_rep rep, int cdim, int vdim)
{
  cs->rep = rep; cs->canonical = false; cs->cdim = cdim; cs->vdim = vdim;
  gkyl_cart_modal_hybrid(&cs->pbasis, cdim, vdim);
  gkyl_cart_modal_tensor(&cs->vbasis, vdim, 2);
  gkyl_cart_modal_tensor(&cs->cbasis, cdim, 1);
  int pdim = cdim + vdim, np = cs->pbasis.num_basis, nv = cs->vbasis.num_basis;
  for (int d = 0; d < pdim; ++d) {
    cs->w[d] = tn_rnd(-1.0, 1.0);
    cs->dxv[d] = tn_rnd(0.5, 2.0);
  }
  for (int d = 0; d < cdim; ++d) {
    cs->jacob_pos[2 * d] = tn_rnd(0.5, 1.5);
    cs->jacob_pos[2 * d + 1] = tn_rnd(-0.2, 0.2);
  }
  // Monotone cubic velocity maps; the Jacobians (2/dv) dv/deta at the 3 Gauss-Legendre
  // points (volume) and at the 3 or 4 surface nodes (standard or ho_ producers).
  for (int b = 0; b < vdim; ++b) {
    double *vm = &cs->vmap[4 * b], dv = cs->dxv[cdim + b];
    vm[0] = tn_rnd(-1.0, 1.0);
    vm[1] = 0.5 * dv * sqrt(2.0 / 3.0) * tn_rnd(0.8, 1.2);
    vm[2] = 0.05 * dv * tn_rnd(-1.0, 1.0);
    vm[3] = 0.05 * dv * tn_rnd(-1.0, 1.0);
    for (int i = 0; i < 3; ++i)
      cs->jacob_vel[3 * b + i] = (2.0 / dv) * tn_vmap_deriv(vm, tn_gl_x[3][i]);
    for (int i = 0; i < 3; ++i)
      cs->jvs[0][3 * b + i] = (2.0 / dv) * tn_vmap_deriv(vm, tn_gl_x[3][i]);
    for (int i = 0; i < 4; ++i)
      cs->jvs[1][4 * b + i] = (2.0 / dv) * tn_vmap_deriv(vm, tn_gl_x[4][i]);
  }
  int npt = tn_npt[vdim];
  for (int i = 0; i < (1 << cdim) * npt; ++i) cs->pt_vol[i] = tn_rnd(-1.0, 1.0);
  for (int i = 0; i < cdim * 2 * (1 << (cdim - 1)) * vdim; ++i) cs->pt_surf[i] = tn_rnd(-1.0, 1.0);
  for (int k = 0; k < np; ++k) cs->f[k] = tn_rnd(-1.0, 1.0);
  for (int k = 0; k < 6 * cs->cbasis.num_basis; ++k) cs->qmem[k] = tn_rnd(-1.0, 1.0);
  int nh = rep == TN_PHASE ? np : nv;
  for (int k = 0; k < nh; ++k) cs->hamil[k] = tn_rnd(-1.0, 1.0);
  if (rep == TN_SPARSE) {
    for (int m = 0; m < nv; ++m)
      if (!tn_is_sparse_mode(&cs->vbasis, m)) cs->hamil[m] = 0.0;
  }
}

// Derivative of H along phase-space direction dir at (x, eta).
static double
tn_ham_grad(const struct tn_case *cs, int dir, const double *x, const double *eta)
{
  int cdim = cs->cdim, vdim = cs->vdim;
  if (cs->rep == TN_PHASE) {
    double z[6];
    for (int d = 0; d < cdim; ++d) z[d] = x[d];
    for (int b = 0; b < vdim; ++b) z[cdim + b] = eta[b];
    return cs->pbasis.eval_grad_expand(dir, z, cs->hamil);
  }
  return dir < cdim ? 0.0 : cs->vbasis.eval_grad_expand(dir - cdim, eta, cs->hamil);
}

// f divided by the velocity-map Jacobian at the 3 Gauss-Legendre points along velocity
// direction a (and then b, when b >= 0), re-interpolated, at the phase-space point z.
static double
tn_f_div(const struct tn_case *cs, const double *z, int a, int b)
{
  const double *xi = tn_gl_x[3];
  int ia = cs->cdim + a, ib = cs->cdim + b;
  double zz[6], sum = 0.0;
  memcpy(zz, z, sizeof zz);
  for (int i = 0; i < 3; ++i) {
    double la = 1.0;
    for (int m = 0; m < 3; ++m)
      if (m != i) la *= (z[ia] - xi[m]) / (xi[i] - xi[m]);
    zz[ia] = xi[i];
    if (b < 0) {
      sum += la * cs->pbasis.eval_expand(zz, cs->f) / cs->jacob_vel[3 * a + i];
      continue;
    }
    for (int j = 0; j < 3; ++j) {
      double lb = 1.0;
      for (int m = 0; m < 3; ++m)
        if (m != j) lb *= (z[ib] - xi[m]) / (xi[j] - xi[m]);
      zz[ib] = xi[j];
      sum += la * lb * cs->pbasis.eval_expand(zz, cs->f)
        / (cs->jacob_vel[3 * a + i] * cs->jacob_vel[3 * b + j]);
    }
  }
  return sum;
}

// Index of Omega_ab^j (a < b) in the Poisson-tensor components.
static int
tn_om_idx(int vdim, int a, int b, int j)
{
  if (vdim == 2) return 4 + j;
  int pair = a == 0 ? b - 1 : 2;
  return 9 + 3 * j + pair;
}

// Pi^{ab} = sum_m Omega_ab^m v_m at Poisson-tensor storage pt (antisymmetric in a, b).
static double
tn_pi_vv(const struct tn_case *cs, const double *pt, int a, int b, const double *eta)
{
  int lo = a < b ? a : b, hi = a < b ? b : a;
  double s = 0.0;
  for (int m = 0; m < cs->vdim; ++m)
    s += pt[tn_om_idx(cs->vdim, lo, hi, m)] * tn_vmap(&cs->vmap[4 * m], eta[m]);
  return a < b ? s : -s;
}

// Brute-force volume term (see the file comment): 3 Gauss-Legendre points per configuration
// and 5 per velocity direction integrate every term exactly.
static void
tn_ref_vol(const struct tn_case *cs, double *out)
{
  const int cdim = cs->cdim, vdim = cs->vdim, pdim = cdim + vdim, np = cs->pbasis.num_basis;
  const int nq = 1 << cdim, npt = tn_npt[vdim];
  double dx1[3], jinv[3], dv1[3];
  for (int d = 0; d < cdim; ++d) {
    dx1[d] = 2.0 / cs->dxv[d];
    jinv[d] = cs->canonical ? 1.0 : 1.0 / cs->jacob_pos[2 * d];
  }
  for (int b = 0; b < vdim; ++b) dv1[b] = 2.0 / cs->dxv[cdim + b];
  for (int k = 0; k < np; ++k) out[k] = 0.0;

  int nquad[6];
  long ntot = 1;
  for (int d = 0; d < pdim; ++d) {
    nquad[d] = d < cdim ? 3 : 5;
    ntot *= nquad[d];
  }
  for (long p = 0; p < ntot; ++p) {
    double z[6], wgt = 1.0;
    long rem = p;
    for (int d = pdim - 1; d >= 0; --d) {
      int i = rem % nquad[d];
      rem /= nquad[d];
      z[d] = tn_gl_x[nquad[d]][i];
      wgt *= tn_gl_w[nquad[d]][i];
    }
    const double *eta = &z[cdim];

    // Interpolated characteristic coefficients from the values at the configuration nodes.
    double A[3][3] = {{0}}, B[3][3] = {{0}}, C[3] = {0}, D[3] = {0};
    int pa[3], pb[3], npair = 0;
    for (int a = 0; a < vdim && !cs->canonical; ++a)
      for (int b = a + 1; b < vdim; ++b) {
        pa[npair] = a; pb[npair] = b; npair += 1;
      }
    for (int q = 0; q < nq; ++q) {
      double xq[3], lq = 1.0;
      tn_node(cdim, q, xq);
      for (int d = 0; d < cdim; ++d)
        lq *= 0.5 * (1.0 + (xq[d] > 0 ? 1.0 : -1.0) * sqrt(3.0) * z[d]);
      const double *pt = &cs->pt_vol[q * npt];
      double hv[3], hx[3];
      for (int b = 0; b < vdim; ++b) hv[b] = tn_ham_grad(cs, cdim + b, xq, eta);
      for (int d = 0; d < cdim; ++d) hx[d] = tn_ham_grad(cs, d, xq, eta);
      for (int d = 0; d < cdim; ++d)
        for (int b = 0; b < vdim; ++b) {
          double pdb = cs->canonical ? (d == b ? 1.0 : 0.0) : pt[d * vdim + b];
          A[d][b] += lq * pdb * hv[b];
          B[d][b] += lq * pdb * hx[d];
        }
      for (int s = 0; s < npair; ++s) {
        double om = tn_pi_vv(cs, pt, pa[s], pb[s], eta);
        C[s] += lq * om * hv[pb[s]];
        D[s] += lq * om * hv[pa[s]];
      }
    }

    double Fb[3], Fab[3];
    for (int b = 0; b < vdim; ++b)
      Fb[b] = cs->canonical ? cs->pbasis.eval_expand(z, cs->f) : tn_f_div(cs, z, b, -1);
    for (int s = 0; s < npair; ++s) Fab[s] = tn_f_div(cs, z, pa[s], pb[s]);

    // Derivatives of the basis functions: central differences, exact for the at most
    // quadratic dependence on each coordinate.
    double dpsi[6][216], bp[216], bm[216];
    for (int d = 0; d < pdim; ++d) {
      double zp[6], zm[6];
      memcpy(zp, z, sizeof zp); memcpy(zm, z, sizeof zm);
      zp[d] += 0.5; zm[d] -= 0.5;
      cs->pbasis.eval(zp, bp); cs->pbasis.eval(zm, bm);
      for (int k = 0; k < np; ++k) dpsi[d][k] = bp[k] - bm[k];
    }

    for (int k = 0; k < np; ++k) {
      double g = 0.0;
      for (int d = 0; d < cdim; ++d)
        for (int b = 0; b < vdim; ++b)
          g += dx1[d] * jinv[d] * dv1[b]
            * (dpsi[d][k] * A[d][b] - dpsi[cdim + b][k] * B[d][b]) * Fb[b];
      for (int s = 0; s < npair; ++s)
        g += dv1[pa[s]] * dv1[pb[s]]
          * (dpsi[cdim + pa[s]][k] * C[s] - dpsi[cdim + pb[s]][k] * D[s]) * Fab[s];
      out[k] += wgt * g;
    }
  }
}

// Relative max-norm difference of a against the reference b.
static double
tn_rel_diff(const double *a, const double *b, int n)
{
  double err = 0.0, nrm = 0.0;
  for (int i = 0; i < n; ++i) {
    err = fmax(err, fabs(a[i] - b[i]));
    nrm = fmax(nrm, fabs(b[i]));
  }
  return err / (nrm > 0.0 ? nrm : 1.0);
}

static void
tn_call_vol(const struct tn_case *cs, tn_vol_t vol, double *out)
{
  for (int k = 0; k < cs->pbasis.num_basis; ++k) out[k] = 0.0;
  vol(cs->w, cs->dxv, cs->vmap, cs->jacob_pos, cs->jacob_vel, cs->pt_vol, cs->hamil, cs->f,
    cs->f, out);
}

static const double tn_tol = 1e-12;

void
test_triad_nodal_vol(void)
{
  srand(1170);
  for (int r = 0; r < 3; ++r)
    for (int i = 0; i < 6; ++i) {
      const struct tn_kernels *kn = &tn_tab[r][i];
      if (!kn->vol) continue;
      int cdim = tn_cv[i][0], vdim = tn_cv[i][1], ncase = cdim + vdim == 6 ? 1 : 3;
      for (int n = 0; n < ncase; ++n) {
        struct tn_case cs;
        tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
        double out[216], ref[216];
        tn_call_vol(&cs, kn->vol, out);
        tn_ref_vol(&cs, ref);
        double rel = tn_rel_diff(out, ref, cs.pbasis.num_basis);
        TEST_CHECK(rel < tn_tol);
        TEST_MSG("%s %dx%dv volume vs definition: rel diff %.2e", tn_rep_name[r], cdim, vdim, rel);
      }
    }
}

// Projection onto the phase basis of the velocity Hamiltonian of cs (constant in x).
static void
tn_project_vel_hamil(const struct tn_case *cs, double *hp)
{
  const int cdim = cs->cdim, vdim = cs->vdim, pdim = cdim + vdim, np = cs->pbasis.num_basis;
  for (int k = 0; k < np; ++k) hp[k] = 0.0;
  long ntot = 1;
  for (int d = 0; d < pdim; ++d) ntot *= 3;
  for (long p = 0; p < ntot; ++p) {
    double z[6], wgt = 1.0, bz[216];
    long rem = p;
    for (int d = pdim - 1; d >= 0; --d) {
      int i = rem % 3;
      rem /= 3;
      z[d] = tn_gl_x[3][i];
      wgt *= tn_gl_w[3][i];
    }
    double h = cs->vbasis.eval_expand(&z[cdim], cs->hamil);
    cs->pbasis.eval(z, bz);
    for (int k = 0; k < np; ++k) hp[k] += wgt * bz[k] * h;
  }
}

void
test_triad_nodal_reps(void)
{
  srand(2026);
  for (int i = 0; i < 6; ++i) {
    int cdim = tn_cv[i][0], vdim = tn_cv[i][1];
    // dense == sparse for a separable Hamiltonian (volume and surface forces)
    const struct tn_kernels *ks = &tn_tab[TN_SPARSE][i], *kd = &tn_tab[TN_DENSE][i];
    if (ks->vol && kd->vol) {
      struct tn_case cs;
      tn_case_init(&cs, TN_SPARSE, cdim, vdim);
      double o1[216], o2[216];
      tn_call_vol(&cs, ks->vol, o1);
      tn_call_vol(&cs, kd->vol, o2);
      double rel = tn_rel_diff(o1, o2, cs.pbasis.num_basis);
      TEST_CHECK(rel < tn_tol);
      TEST_MSG("%dx%dv volume, sparse vs dense: rel diff %.2e", cdim, vdim, rel);
      for (int ho = 0; ho < 2; ++ho) {
        int nn = ho ? 4 : 3;
        for (int d = 0; d < cdim; ++d) {
          int nalpha = (1 << (cdim - 1)) * (int)lround(pow(nn, vdim));
          double a1[256] = {0}, a2[256] = {0};
          ks->conf[ho][d](cs.w, cs.dxv, -1, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_surf,
            cs.hamil, a1);
          kd->conf[ho][d](cs.w, cs.dxv, -1, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_surf,
            cs.hamil, a2);
          rel = tn_rel_diff(a1, a2, nalpha);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%dx%dv conf force dir %d ho %d, sparse vs dense: rel diff %.2e", cdim, vdim,
            d, ho, rel);
        }
        for (int a = 0; a < vdim; ++a) {
          int nalpha = (1 << cdim) * (int)lround(pow(nn, vdim - 1));
          double a1[256] = {0}, a2[256] = {0};
          ks->vel[ho][a](cs.w, cs.dxv, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_vol, cs.hamil, a1);
          kd->vel[ho][a](cs.w, cs.dxv, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_vol, cs.hamil, a2);
          rel = tn_rel_diff(a1, a2, nalpha);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%dx%dv vel force dir %d ho %d, sparse vs dense: rel diff %.2e", cdim, vdim,
            a, ho, rel);
        }
      }
    }
    // phase == dense for a Hamiltonian independent of the configuration-space coordinates
    const struct tn_kernels *kp = &tn_tab[TN_PHASE][i];
    if (kp->vol && kd->vol) {
      struct tn_case cs;
      tn_case_init(&cs, TN_DENSE, cdim, vdim);
      struct tn_case csp = cs;
      csp.rep = TN_PHASE;
      tn_project_vel_hamil(&cs, csp.hamil);
      double o1[216], o2[216];
      tn_call_vol(&csp, kp->vol, o1);
      tn_call_vol(&cs, kd->vol, o2);
      double rel = tn_rel_diff(o1, o2, cs.pbasis.num_basis);
      TEST_CHECK(rel < tn_tol);
      TEST_MSG("%dx%dv volume, phase vs dense: rel diff %.2e", cdim, vdim, rel);
    }
  }
}

// alpha(i, j) = sum_t O[t*NO + i]*I[t*NI + j] from a _shared producer (thread loop on host).
static void
tn_compose(int nt, int no, int ni, const double *O, const double *I, double *alpha)
{
  for (int i = 0; i < no; ++i)
    for (int j = 0; j < ni; ++j) {
      double s = 0.0;
      for (int t = 0; t < nt; ++t) s += O[t * no + i] * I[t * ni + j];
      alpha[i * ni + j] = s;
    }
}

void
test_triad_nodal_conf_force(void)
{
  srand(4242);
  for (int r = 0; r < 3; ++r)
    for (int i = 0; i < 6; ++i) {
      const struct tn_kernels *kn = &tn_tab[r][i];
      if (!kn->vol) continue;
      int cdim = tn_cv[i][0], vdim = tn_cv[i][1];
      struct tn_case cs;
      tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
      int ns = 1 << (cdim - 1);
      for (int ho = 0; ho < 2; ++ho) {
        int nn = ho ? 4 : 3, ni = (int)lround(pow(nn, vdim));
        for (int d = 0; d < cdim; ++d)
          for (int edge = -1; edge <= 1; edge += 2) {
            double alpha[256] = {0}, ref[256], comp[256];
            kn->conf[ho][d](cs.w, cs.dxv, edge, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_surf,
              cs.hamil, alpha);
            const double *pt_face = &cs.pt_surf[(d * 2 + (edge == -1 ? 0 : 1)) * ns * vdim];
            for (int s = 0; s < ns; ++s) {
              double xt[3], x[3];
              tn_node(cdim - 1, s, xt);
              for (int dd = 0, t = 0; dd < cdim; ++dd) x[dd] = dd == d ? (double)edge : xt[t++];
              for (int j = 0; j < ni; ++j) {
                double eta[3];
                int jb[3];
                for (int b = 0, rem = j; b < vdim; ++b) {
                  int stride = (int)lround(pow(nn, vdim - 1 - b));
                  jb[b] = rem / stride;
                  rem %= stride;
                  eta[b] = tn_gl_x[nn][jb[b]];
                }
                double a = 0.0;
                for (int b = 0; b < vdim; ++b)
                  a += pt_face[s * vdim + b] * tn_ham_grad(&cs, cdim + b, x, eta)
                    * (2.0 / cs.dxv[cdim + b]) / cs.jvs[ho][nn * b + jb[b]];
                ref[s * ni + j] = a;
              }
            }
            double rel = tn_rel_diff(alpha, ref, ns * ni);
            TEST_CHECK(rel < tn_tol);
            TEST_MSG("%s %dx%dv conf force dir %d edge %d ho %d vs definition: rel diff %.2e",
              tn_rep_name[r], cdim, vdim, d, edge, ho, rel);

            int nt = kn->conf_sh[ho][d](0, 0, cs.w, cs.dxv, edge, cs.vmap, cs.jacob_pos,
              cs.jvs[ho], cs.pt_surf, cs.hamil, NULL, NULL);
            double *O = (double *)gkyl_malloc(sizeof(double) * nt * ns);
            double *I = (double *)gkyl_malloc(sizeof(double) * nt * ni);
            for (int tid = 0; tid < (ns > ni ? ns : ni); ++tid)
              kn->conf_sh[ho][d](tid, 0, cs.w, cs.dxv, edge, cs.vmap, cs.jacob_pos, cs.jvs[ho],
                cs.pt_surf, cs.hamil, O, I);
            tn_compose(nt, ns, ni, O, I, comp);
            rel = tn_rel_diff(comp, alpha, ns * ni);
            TEST_CHECK(rel < tn_tol);
            TEST_MSG("%s %dx%dv conf force dir %d edge %d ho %d, _shared vs plain: rel diff %.2e",
              tn_rep_name[r], cdim, vdim, d, edge, ho, rel);
            gkyl_free(O);
            gkyl_free(I);
          }
      }
    }
}

void
test_triad_nodal_vel_force(void)
{
  srand(777);
  for (int r = 0; r < 3; ++r)
    for (int i = 0; i < 6; ++i) {
      const struct tn_kernels *kn = &tn_tab[r][i];
      if (!kn->vol) continue;
      int cdim = tn_cv[i][0], vdim = tn_cv[i][1], nq = 1 << cdim, npt = tn_npt[vdim];
      struct tn_case cs;
      tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
      for (int ho = 0; ho < 2; ++ho) {
        int nn = ho ? 4 : 3, ni = (int)lround(pow(nn, vdim - 1));
        for (int a = 0; a < vdim; ++a) {
          double alpha[256] = {0}, ref[256], comp[256];
          kn->vel[ho][a](cs.w, cs.dxv, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_vol, cs.hamil,
            alpha);
          for (int q = 0; q < nq; ++q) {
            double xq[3];
            tn_node(cdim, q, xq);
            const double *pt = &cs.pt_vol[q * npt];
            for (int j = 0; j < ni; ++j) {
              double eta[3];
              int jb[3];
              for (int b = 0, t = 0, rem = j; b < vdim; ++b) {
                if (b == a) {
                  eta[b] = -1.0;
                  jb[b] = -1;
                  continue;
                }
                int stride = (int)lround(pow(nn, vdim - 2 - t));
                jb[b] = rem / stride;
                rem %= stride;
                eta[b] = tn_gl_x[nn][jb[b]];
                t += 1;
              }
              double al = 0.0;
              for (int c = 0; c < cdim; ++c)
                al -= pt[c * vdim + a] * (2.0 / cs.dxv[c]) / cs.jacob_pos[2 * c]
                  * tn_ham_grad(&cs, c, xq, eta);
              for (int b = 0; b < vdim; ++b) {
                if (b == a) continue;
                al += tn_pi_vv(&cs, pt, a, b, eta) * tn_ham_grad(&cs, cdim + b, xq, eta)
                  * (2.0 / cs.dxv[cdim + b]) / cs.jvs[ho][nn * b + jb[b]];
              }
              ref[q * ni + j] = al;
            }
          }
          double rel = tn_rel_diff(alpha, ref, nq * ni);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%s %dx%dv vel force dir %d ho %d vs definition: rel diff %.2e",
            tn_rep_name[r], cdim, vdim, a, ho, rel);

          int nt = kn->vel_sh[ho][a](0, 0, cs.w, cs.dxv, cs.vmap, cs.jacob_pos, cs.jvs[ho],
            cs.pt_vol, cs.hamil, NULL, NULL);
          if (nt == 0) continue; // no force in this direction (1V velocity Hamiltonian)
          double *O = (double *)gkyl_malloc(sizeof(double) * nt * nq);
          double *I = (double *)gkyl_malloc(sizeof(double) * nt * ni);
          for (int tid = 0; tid < (nq > ni ? nq : ni); ++tid)
            kn->vel_sh[ho][a](tid, 0, cs.w, cs.dxv, cs.vmap, cs.jacob_pos, cs.jvs[ho], cs.pt_vol,
              cs.hamil, O, I);
          tn_compose(nt, nq, ni, O, I, comp);
          rel = tn_rel_diff(comp, alpha, nq * ni);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%s %dx%dv vel force dir %d ho %d, _shared vs plain: rel diff %.2e",
            tn_rep_name[r], cdim, vdim, a, ho, rel);
          gkyl_free(O);
          gkyl_free(I);
        }
      }
    }
}

// Magnetic-field force volume kernels, [Hamiltonian][cv][B direction].
typedef void (*tn_bvol_t)(const double *w, const double *dxv, const double *jacob_vel,
  const double *hamil, const double *qmem, const double *f, double *out);
#define TN_B(r, d, cv) vlasov_B##d##_hamil_##r##_vol_##cv##_tensor_p1
#define TN_B2(r, cv) { NULL, NULL, TN_B(r, z, cv) }
#define TN_B3(r, cv) { TN_B(r, x, cv), TN_B(r, y, cv), TN_B(r, z, cv) }
static const tn_bvol_t tn_bvol[3][6][3] = {
  { {0}, TN_B2(vel_sparse, 1x2v), TN_B2(vel_sparse, 2x2v), TN_B3(vel_sparse, 1x3v),
    TN_B3(vel_sparse, 2x3v), TN_B3(vel_sparse, 3x3v) },
  { {0}, TN_B2(vel_dense, 1x2v), TN_B2(vel_dense, 2x2v), TN_B3(vel_dense, 1x3v),
    TN_B3(vel_dense, 2x3v), TN_B3(vel_dense, 3x3v) },
  { {0}, TN_B2(phase, 1x2v), TN_B2(phase, 2x2v), TN_B3(phase, 1x3v), TN_B3(phase, 2x3v),
    TN_B3(phase, 3x3v) },
};

// f times the velocity-map Jacobian at the 3 Gauss-Legendre points along velocity
// direction t, re-interpolated, at the phase-space point z.
static double
tn_f_mul(const struct tn_case *cs, const double *z, int t)
{
  const double *xi = tn_gl_x[3];
  int it = cs->cdim + t;
  double zz[6], sum = 0.0;
  memcpy(zz, z, sizeof zz);
  for (int i = 0; i < 3; ++i) {
    double l = 1.0;
    for (int m = 0; m < 3; ++m)
      if (m != i) l *= (z[it] - xi[m]) / (xi[i] - xi[m]);
    zz[it] = xi[i];
    sum += l * cs->pbasis.eval_expand(zz, cs->f) * cs->jacob_vel[3 * t + i];
  }
  return sum;
}

// Brute-force magnetic-field force volume term (see the file comment) for B along bdir.
static void
tn_ref_bvol(const struct tn_case *cs, int bdir, double *out)
{
  const int cdim = cs->cdim, vdim = cs->vdim, pdim = cdim + vdim, np = cs->pbasis.num_basis;
  const int nq = 1 << cdim;
  const int cyc[3][3] = { {1, 2, 0}, {2, 0, 1}, {0, 1, 2} };
  int a = vdim == 2 ? 0 : cyc[bdir][0], b = vdim == 2 ? 1 : cyc[bdir][1], t = vdim == 2 ? -1 : cyc[bdir][2];
  const double *Bc = &cs->qmem[(3 + bdir) * cs->cbasis.num_basis];
  double dv1a = 2.0 / cs->dxv[cdim + a], dv1b = 2.0 / cs->dxv[cdim + b];
  for (int k = 0; k < np; ++k) out[k] = 0.0;

  int nquad[6];
  long ntot = 1;
  for (int d = 0; d < pdim; ++d) {
    nquad[d] = d < cdim ? 3 : 5;
    ntot *= nquad[d];
  }
  for (long p = 0; p < ntot; ++p) {
    double z[6], wgt = 1.0;
    long rem = p;
    for (int d = pdim - 1; d >= 0; --d) {
      int i = rem % nquad[d];
      rem /= nquad[d];
      z[d] = tn_gl_x[nquad[d]][i];
      wgt *= tn_gl_w[nquad[d]][i];
    }
    const double *eta = &z[cdim];
    double ga = 0.0, gb = 0.0; // I[B d_vb(H)], I[B d_va(H)]
    for (int q = 0; q < nq; ++q) {
      double xq[3], lq = 1.0;
      tn_node(cdim, q, xq);
      for (int d = 0; d < cdim; ++d)
        lq *= 0.5 * (1.0 + (xq[d] > 0 ? 1.0 : -1.0) * sqrt(3.0) * z[d]);
      double Bq = cs->cbasis.eval_expand(xq, Bc);
      ga += lq * Bq * tn_ham_grad(cs, cdim + b, xq, eta);
      gb += lq * Bq * tn_ham_grad(cs, cdim + a, xq, eta);
    }
    double F = t < 0 ? cs->pbasis.eval_expand(z, cs->f) :
      (cs->rep == TN_PHASE ? cs->pbasis.eval_expand(z, cs->f) * cs->jacob_vel[3 * t] :
                             tn_f_mul(cs, z, t));
    double bp[216], bm[216], zp[6], zm[6];
    for (int s = 0; s < 2; ++s) {
      int d = cdim + (s == 0 ? a : b);
      memcpy(zp, z, sizeof zp); memcpy(zm, z, sizeof zm);
      zp[d] += 0.5; zm[d] -= 0.5;
      cs->pbasis.eval(zp, bp); cs->pbasis.eval(zm, bm);
      double g = s == 0 ? ga : -gb;
      for (int k = 0; k < np; ++k) out[k] += wgt * dv1a * dv1b * (bp[k] - bm[k]) * g * F;
    }
  }
}

void
test_triad_nodal_bvol(void)
{
  srand(8086);
  for (int r = 0; r < 3; ++r)
    for (int i = 1; i < 6; ++i)
      for (int bdir = 0; bdir < 3; ++bdir) {
        tn_bvol_t kern = tn_bvol[r][i][bdir];
        if (!kern) continue;
        int cdim = tn_cv[i][0], vdim = tn_cv[i][1], ncase = cdim + vdim == 6 ? 1 : 2;
        for (int n = 0; n < ncase; ++n) {
          struct tn_case cs;
          tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
          double out[216] = {0}, ref[216];
          kern(cs.w, cs.dxv, cs.jacob_vel, cs.hamil, cs.qmem, cs.f, out);
          tn_ref_bvol(&cs, bdir, ref);
          double rel = tn_rel_diff(out, ref, cs.pbasis.num_basis);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%s %dx%dv B%c volume vs definition: rel diff %.2e", tn_rep_name[r], cdim, vdim,
            "xyz"[bdir], rel);
        }
      }
}

// Velocity-Jacobian division and rescaling kernels, [cv] over tn_mcv.
typedef void (*tn_jv_t)(const double *jacob_vel_gauss, const double *in, double *out);
static const int tn_mcv[7][2] = { {1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 3}, {2, 3}, {3, 3} };
#define TN_JV(op, cv) vlasov_##op##_Jv_##cv##_tensor_p1
static const tn_jv_t tn_jv[2][7] = {
  { TN_JV(divide, 1x1v), TN_JV(divide, 1x2v), TN_JV(divide, 2x1v), TN_JV(divide, 2x2v),
    TN_JV(divide, 1x3v), TN_JV(divide, 2x3v), TN_JV(divide, 3x3v) },
  { TN_JV(rescale, 1x1v), TN_JV(rescale, 1x2v), TN_JV(rescale, 2x1v), TN_JV(rescale, 2x2v),
    TN_JV(rescale, 1x3v), TN_JV(rescale, 2x3v), TN_JV(rescale, 3x3v) },
};

void
test_jacob_vel_factored(void)
{
  srand(4242);
  for (int op = 0; op < 2; ++op)
    for (int i = 0; i < 7; ++i) {
      if (!tn_jv[op][i]) continue;
      int cdim = tn_mcv[i][0], vdim = tn_mcv[i][1], pdim = cdim + vdim, nqv = 1;
      for (int d = 0; d < vdim; ++d) nqv *= 3;
      struct gkyl_basis pb;
      gkyl_cart_modal_hybrid(&pb, cdim, vdim);
      double in[216], out[216], jac[27];
      for (int k = 0; k < pb.num_basis; ++k) in[k] = tn_rnd(-1.0, 1.0);
      for (int j = 0; j < nqv; ++j) jac[j] = tn_rnd(0.5, 1.5);
      tn_jv[op][i](jac, in, out);
      double err = 0.0, nrm = 0.0;
      for (int r = 0; r < 4; ++r) {
        double z[6];
        for (int d = 0; d < cdim; ++d) z[d] = tn_rnd(-1.0, 1.0);
        // Velocity nodes in jacob_vel_gauss order (first direction slowest).
        for (int j = 0; j < nqv; ++j) {
          for (int d = 0; d < vdim; ++d) {
            int s = 1;
            for (int e = d + 1; e < vdim; ++e) s *= 3;
            z[cdim + d] = tn_gl_x[3][(j / s) % 3];
          }
          double a = pb.eval_expand(z, in), b = pb.eval_expand(z, out);
          double want = op == 0 ? a : b, got = op == 0 ? jac[j] * b : jac[j] * a;
          err = fmax(err, fabs(got - want));
          nrm = fmax(nrm, fabs(want));
        }
      }
      (void)pdim;
      TEST_CHECK(err / nrm < tn_tol);
      TEST_MSG("%s %dx%dv: J*out vs in at the velocity nodes, rel diff %.2e", op == 0 ? "divide" :
        "rescale", cdim, vdim, err / nrm);
    }
}

// Moment kernels, [rep][moment: M1i, five, int_five][cv] over tn_mcv.
typedef void (*tn_mom_t)(const double *w, const double *dxv, const int *idx, const double *jacob_vel,
  const double *hamil, const double *f, double *out);
#define TN_MOM(r, m, cv) mom_vlasov_hamil_##r##_##m##_##cv##_tensor_p1
#define TN_MOMS(r, m) { TN_MOM(r, m, 1x1v), TN_MOM(r, m, 1x2v), TN_MOM(r, m, 2x1v), \
  TN_MOM(r, m, 2x2v), TN_MOM(r, m, 1x3v), TN_MOM(r, m, 2x3v), TN_MOM(r, m, 3x3v) }
#define TN_MOMS_NO1V(r, m) { NULL, TN_MOM(r, m, 1x2v), NULL, TN_MOM(r, m, 2x2v), \
  TN_MOM(r, m, 1x3v), TN_MOM(r, m, 2x3v), TN_MOM(r, m, 3x3v) }
static const tn_mom_t tn_mom[3][3][7] = {
  { TN_MOMS_NO1V(vel_sparse, M1i), TN_MOMS_NO1V(vel_sparse, five_moments),
    TN_MOMS_NO1V(vel_sparse, int_five_moments) },
  { TN_MOMS(vel_dense, M1i), TN_MOMS(vel_dense, five_moments), TN_MOMS(vel_dense, int_five_moments) },
  { TN_MOMS(phase, M1i), TN_MOMS(phase, five_moments), TN_MOMS(phase, int_five_moments) },
};

// H at the phase-space point z.
static double
tn_ham_val(const struct tn_case *cs, const double *z)
{
  return cs->rep == TN_PHASE ? cs->pbasis.eval_expand(z, cs->hamil) :
    cs->vbasis.eval_expand(&z[cs->cdim], cs->hamil);
}

// Brute-force moments (see the file comment): mom 0 = M1i, 1 = five moments, 2 = integrated.
static void
tn_ref_mom(const struct tn_case *cs, int mom, double *out)
{
  const int cdim = cs->cdim, vdim = cs->vdim, pdim = cdim + vdim, nc = cs->cbasis.num_basis;
  int nk = mom == 2 ? 1 : nc, nout = mom == 0 ? vdim * nc : mom == 1 ? (vdim + 2) * nc : vdim + 2;
  double volFact = 1.0;
  for (int d = (mom == 2 ? 0 : cdim); d < pdim; ++d) volFact *= 0.5 * cs->dxv[d];
  for (int k = 0; k < nout; ++k) out[k] = 0.0;
  int nquad[6];
  long ntot = 1;
  for (int d = 0; d < pdim; ++d) {
    nquad[d] = d < cdim ? 3 : 5;
    ntot *= nquad[d];
  }
  for (long p = 0; p < ntot; ++p) {
    double z[6], wgt = 1.0, bc[8];
    long rem = p;
    for (int d = pdim - 1; d >= 0; --d) {
      int i = rem % nquad[d];
      rem /= nquad[d];
      z[d] = tn_gl_x[nquad[d]][i];
      wgt *= tn_gl_w[nquad[d]][i];
    }
    cs->cbasis.eval(z, bc);
    double f = cs->pbasis.eval_expand(z, cs->f), h = tn_ham_val(cs, z);
    for (int k = 0; k < nk; ++k) {
      double psi = mom == 2 ? 1.0 : bc[k], c = wgt * volFact * psi;
      int stride = mom == 2 ? 1 : nc, m1 = mom == 0 ? 0 : 1;
      if (mom > 0) out[k] += c * f;
      for (int d = 0; d < vdim; ++d)
        out[(m1 + d) * stride + k] += c * (2.0 / cs->dxv[cdim + d]) *
          tn_ham_grad(cs, cdim + d, z, &z[cdim]) * tn_f_div(cs, z, d, -1);
      if (mom > 0) out[(vdim + 1) * stride + k] += c * h * f;
    }
  }
}

void
test_mom_factored(void)
{
  static const char *tn_mom_name[3] = { "M1i", "five_moments", "int_five_moments" };
  srand(2718);
  for (int r = 0; r < 3; ++r)
    for (int mom = 0; mom < 3; ++mom)
      for (int i = 0; i < 7; ++i) {
        tn_mom_t kern = tn_mom[r][mom][i];
        if (!kern) continue;
        int cdim = tn_mcv[i][0], vdim = tn_mcv[i][1];
        struct tn_case cs;
        tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
        double out[64] = {0}, ref[64];
        int idx[6] = {0};
        kern(cs.w, cs.dxv, idx, cs.jacob_vel, cs.hamil, cs.f, out);
        tn_ref_mom(&cs, mom, ref);
        int nc = cs.cbasis.num_basis, nout = mom == 0 ? vdim * nc : mom == 1 ? (vdim + 2) * nc : vdim + 2;
        double rel = tn_rel_diff(out, ref, nout);
        TEST_CHECK(rel < tn_tol);
        TEST_MSG("%s %s %dx%dv vs definition: rel diff %.2e", tn_rep_name[r], tn_mom_name[mom], cdim,
          vdim, rel);
      }
}

// Velocity-space Hamiltonian streaming and force kernels, [rep: sparse, dense][cv over tn_mcv].
typedef void (*tn_svol_t)(const double *w, const double *dxv, const double *vmap, const double *jacob_pos,
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double *out);
typedef double (*tn_ssurf_t)(const double *w, const double *dxv, const double *jacob_pos_l,
  const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil, const double *fl, const double *fc,
  const double *fr, double *out);
typedef double (*tn_sbnd_t)(const double *w, const double *dxv, const double *jacob_pos_edge,
  const double *jacob_pos_skin, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const int edge, const double *fedge, const double *fskin, double *out);
typedef void (*tn_fvol_t)(const double *w, const double *dxv, const double *jacob_vel, const double *a,
  const double *f, double *out);
typedef void (*tn_phivol_t)(const double *w, const double *dxv, const double *jacob_pos,
  const double *jacob_vel, const double *phi, const double *f, double *out);
#define TN_SK(r, k, cv) vlasov_hamil_##r##_##k##_##cv##_tensor_p1
#define TN_S3(r, k) { NULL, TN_SK(r, k ## x, 1x2v), NULL, TN_SK(r, k ## x, 2x2v), TN_SK(r, k ## x, 1x3v), \
  TN_SK(r, k ## x, 2x3v), TN_SK(r, k ## x, 3x3v) }, { NULL, NULL, NULL, TN_SK(r, k ## y, 2x2v), NULL, \
  TN_SK(r, k ## y, 2x3v), TN_SK(r, k ## y, 3x3v) }, { NULL, NULL, NULL, NULL, NULL, NULL, TN_SK(r, k ## z, 3x3v) }
#define TN_D3(r, k) { TN_SK(r, k ## x, 1x1v), TN_SK(r, k ## x, 1x2v), TN_SK(r, k ## x, 2x1v), \
  TN_SK(r, k ## x, 2x2v), TN_SK(r, k ## x, 1x3v), TN_SK(r, k ## x, 2x3v), TN_SK(r, k ## x, 3x3v) }, \
  { NULL, NULL, NULL, TN_SK(r, k ## y, 2x2v), NULL, TN_SK(r, k ## y, 2x3v), TN_SK(r, k ## y, 3x3v) }, \
  { NULL, NULL, NULL, NULL, NULL, NULL, TN_SK(r, k ## z, 3x3v) }
static const tn_svol_t tn_svol[2][7] = {
  { NULL, TN_SK(vel_sparse, vol, 1x2v), NULL, TN_SK(vel_sparse, vol, 2x2v), TN_SK(vel_sparse, vol, 1x3v),
    TN_SK(vel_sparse, vol, 2x3v), TN_SK(vel_sparse, vol, 3x3v) },
  { TN_SK(vel_dense, vol, 1x1v), TN_SK(vel_dense, vol, 1x2v), TN_SK(vel_dense, vol, 2x1v),
    TN_SK(vel_dense, vol, 2x2v), TN_SK(vel_dense, vol, 1x3v), TN_SK(vel_dense, vol, 2x3v),
    TN_SK(vel_dense, vol, 3x3v) },
};
static const tn_ssurf_t tn_ssurf[2][3][7] = { { TN_S3(vel_sparse, surf) }, { TN_D3(vel_dense, surf) } };
static const tn_sbnd_t tn_sbnd[2][3][7] = { { TN_S3(vel_sparse, boundary_surf) },
  { TN_D3(vel_dense, boundary_surf) } };
#define TN_FV(k) { vlasov_##k##_vol_1x1v_tensor_p1, vlasov_##k##_vol_1x2v_tensor_p1, \
  vlasov_##k##_vol_2x1v_tensor_p1, vlasov_##k##_vol_2x2v_tensor_p1, vlasov_##k##_vol_1x3v_tensor_p1, \
  vlasov_##k##_vol_2x3v_tensor_p1, vlasov_##k##_vol_3x3v_tensor_p1 }
static const tn_fvol_t tn_evol[7] = TN_FV(E), tn_radvol[7] = TN_FV(rad);
static const tn_phivol_t tn_phivol[7] = TN_FV(phi);

// F_d (f divided by J_vd at the 3 Gauss-Legendre points along v_d, re-interpolated) of the
// expansion fa at z.
static double
tn_fdiv_arr(const struct tn_case *cs, const double *fa, const double *z, int d)
{
  const double *xi = tn_gl_x[3];
  int id = cs->cdim + d;
  double zz[6], sum = 0.0;
  memcpy(zz, z, sizeof zz);
  for (int i = 0; i < 3; ++i) {
    double l = 1.0;
    for (int m = 0; m < 3; ++m)
      if (m != i) l *= (z[id] - xi[m]) / (xi[i] - xi[m]);
    zz[id] = xi[i];
    sum += l * cs->pbasis.eval_expand(zz, fa) / cs->jacob_vel[3 * d + i];
  }
  return sum;
}

// Gauss quadrature over the phase-space cell (3 points per configuration-space, 5 per velocity-space
// direction), with direction fix (if >= 0) held at fixv: calls integrand(z, wgt, ctx).
static void
tn_quad(const struct tn_case *cs, int fix, double fixv, void (*integrand)(const double *z, double wgt, void *ctx),
  void *ctx)
{
  int pdim = cs->cdim + cs->vdim, nquad[6];
  long ntot = 1;
  for (int d = 0; d < pdim; ++d) {
    nquad[d] = d == fix ? 1 : d < cs->cdim ? 3 : 5;
    ntot *= nquad[d];
  }
  for (long p = 0; p < ntot; ++p) {
    double z[6], wgt = 1.0;
    long rem = p;
    for (int d = pdim - 1; d >= 0; --d) {
      int i = rem % nquad[d];
      rem /= nquad[d];
      if (d == fix) { z[d] = fixv; continue; }
      z[d] = tn_gl_x[nquad[d]][i];
      wgt *= tn_gl_w[nquad[d]][i];
    }
    integrand(z, wgt, ctx);
  }
}

struct tn_sctx {
  const struct tn_case *cs;
  int kind, d; // kind: 0 streaming vol, 1 E, 2 phi, 3 rad, 4 face
  const double *a, *U, *fu; // force coefficients; face: upwind expansion fu divided by U[0]
  double scale, *out;
};

static void
tn_sint(const double *z, double wgt, void *ctx)
{
  struct tn_sctx *c = ctx;
  const struct tn_case *cs = c->cs;
  int cdim = cs->cdim, d = c->d, np = cs->pbasis.num_basis;
  double b[216], g;
  if (c->kind == 0 || c->kind == 2) { // derivative of psi along x_d (streaming) or v_d (phi)
    int dd = c->kind == 0 ? d : cdim + d;
    double zp[6], zm[6], bp[216], bm[216];
    memcpy(zp, z, sizeof zp); memcpy(zm, z, sizeof zm);
    zp[dd] += 0.5; zm[dd] -= 0.5;
    cs->pbasis.eval(zp, bp); cs->pbasis.eval(zm, bm);
    for (int k = 0; k < np; ++k) b[k] = bp[k] - bm[k]; // exact derivative (at most quadratic)
  } else if (c->kind == 1 || c->kind == 3) {
    double zp[6], zm[6], bp[216], bm[216];
    memcpy(zp, z, sizeof zp); memcpy(zm, z, sizeof zm);
    zp[cdim + d] += 0.5; zm[cdim + d] -= 0.5;
    cs->pbasis.eval(zp, bp); cs->pbasis.eval(zm, bm);
    for (int k = 0; k < np; ++k) b[k] = bp[k] - bm[k];
  } else {
    cs->pbasis.eval(z, b);
  }
  if (c->kind == 0 || c->kind == 4)
    g = tn_ham_grad(cs, cdim + d, z, &z[cdim]) *
      (c->kind == 4 ? tn_fdiv_arr(cs, c->fu, c->U, d) * c->U[6] : tn_fdiv_arr(cs, cs->f, z, d));
  else if (c->kind == 1)
    g = cs->cbasis.eval_expand(z, &c->a[cs->cbasis.num_basis * d]) * tn_fdiv_arr(cs, cs->f, z, d);
  else if (c->kind == 2)
    g = -cs->cbasis.eval_grad_expand(d, z, c->a) * tn_fdiv_arr(cs, cs->f, z, d);
  else
    g = cs->vbasis.eval_expand(&z[cdim], &c->a[cs->vbasis.num_basis * d]) * tn_fdiv_arr(cs, cs->f, z, d);
  for (int k = 0; k < np; ++k) c->out[k] += wgt * c->scale * b[k] * g;
}

// Face contribution: scale int over x_d = side of psi d_vd(H) U, U = F_d of fu at the point with x_d
// replaced by uside, divided by jinv^-1 (the cell's conf-space Jacobian).
static void
tn_face(const struct tn_case *cs, int d, double side, const double *fu, double uside, double jinv,
  double scale, double *out)
{
  // The integrand evaluates F_d of fu at U (z with x_d -> uside); U[6] carries jinv.
  struct tn_sctx c = { .cs = cs, .kind = 4, .d = d, .fu = fu, .scale = scale, .out = out };
  int pdim = cs->cdim + cs->vdim, nquad[6];
  long ntot = 1;
  for (int e = 0; e < pdim; ++e) {
    nquad[e] = e == d ? 1 : e < cs->cdim ? 3 : 5;
    ntot *= nquad[e];
  }
  for (long p = 0; p < ntot; ++p) {
    double z[6], U[7], wgt = 1.0;
    long rem = p;
    for (int e = pdim - 1; e >= 0; --e) {
      int i = rem % nquad[e];
      rem /= nquad[e];
      if (e == d) { z[e] = side; continue; }
      z[e] = tn_gl_x[nquad[e]][i];
      wgt *= tn_gl_w[nquad[e]][i];
    }
    memcpy(U, z, 6 * sizeof(double));
    U[d] = uside; U[6] = jinv;
    c.U = U;
    tn_sint(z, wgt, &c);
  }
}

void
test_stream_factored(void)
{
  srand(1618);
  for (int i = 0; i < 7; ++i) {
    int cdim = tn_mcv[i][0], vdim = tn_mcv[i][1], nd = cdim < vdim ? cdim : vdim;
    for (int r = 0; r < 2; ++r) {
      struct tn_case cs;
      tn_case_init(&cs, (enum tn_rep)r, cdim, vdim);
      int np = cs.pbasis.num_basis;
      double jp[3][9], fa[3][216], out[216], ref[216];
      for (int s = 0; s < 3; ++s) {
        for (int k = 0; k < 9; ++k) jp[s][k] = tn_rnd(0.5, 1.5);
        for (int k = 0; k < np; ++k) fa[s][k] = tn_rnd(-1.0, 1.0);
      }
      memcpy(fa[1], cs.f, np * sizeof(double));
      // Streaming volume.
      if (tn_svol[r][i]) {
        for (int k = 0; k < np; ++k) out[k] = ref[k] = 0.0;
        tn_svol[r][i](cs.w, cs.dxv, cs.vmap, jp[1], cs.jacob_vel, cs.pt_vol, cs.hamil, cs.f, cs.f, out);
        for (int d = 0; d < nd; ++d) {
          struct tn_sctx c = { .cs = &cs, .kind = 0, .d = d, .out = ref,
            .scale = (2.0 / cs.dxv[d]) * (2.0 / cs.dxv[cdim + d]) / jp[1][3 * d] };
          tn_quad(&cs, -1, 0.0, tn_sint, &c);
        }
        double rel = tn_rel_diff(out, ref, np);
        TEST_CHECK(rel < tn_tol);
        TEST_MSG("%s streaming vol %dx%dv vs definition: rel diff %.2e", tn_rep_name[r], cdim, vdim, rel);
      }
      // Interior and boundary faces, both upwind directions.
      for (int d = 0; d < nd; ++d) for (int sgn = -1; sgn <= 1; sgn += 2) {
        double w[6];
        memcpy(w, cs.w, sizeof w);
        w[cdim + d] = sgn * fabs(w[cdim + d]) + 0.25 * sgn;
        double sc = (2.0 / cs.dxv[d]) * (2.0 / cs.dxv[cdim + d]);
        if (tn_ssurf[r][d][i]) {
          for (int k = 0; k < np; ++k) out[k] = ref[k] = 0.0;
          tn_ssurf[r][d][i](w, cs.dxv, jp[0], jp[1], jp[2], cs.jacob_vel, cs.pt_vol, cs.hamil, fa[0], fa[1],
            fa[2], out);
          // Left face (x_d = -1): upwind from the left cell (+1 side) or this cell (-1 side).
          if (sgn > 0) tn_face(&cs, d, -1.0, fa[0], 1.0, 1.0 / jp[0][3 * d], sc, ref);
          else tn_face(&cs, d, -1.0, fa[1], -1.0, 1.0 / jp[1][3 * d], sc, ref);
          if (sgn > 0) tn_face(&cs, d, 1.0, fa[1], 1.0, 1.0 / jp[1][3 * d], -sc, ref);
          else tn_face(&cs, d, 1.0, fa[2], -1.0, 1.0 / jp[2][3 * d], -sc, ref);
          double rel = tn_rel_diff(out, ref, np);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%s surf%c %dx%dv (w_v %+d) vs definition: rel diff %.2e", tn_rep_name[r], "xyz"[d], cdim,
            vdim, sgn, rel);
        }
        if (tn_sbnd[r][d][i]) for (int edge = -1; edge <= 1; edge += 2) {
          // fedge = fa[0], fskin = fa[1]; the face is the skin cell's upper (edge = -1) or lower face.
          for (int k = 0; k < np; ++k) out[k] = ref[k] = 0.0;
          tn_sbnd[r][d][i](w, cs.dxv, jp[0], jp[1], cs.jacob_vel, cs.pt_vol, cs.hamil, edge, fa[0], fa[1], out);
          if (edge == -1) {
            if (sgn > 0) tn_face(&cs, d, 1.0, fa[1], 1.0, 1.0 / jp[1][3 * d], -sc, ref);
            else tn_face(&cs, d, 1.0, fa[0], -1.0, 1.0 / jp[0][3 * d], -sc, ref);
          } else {
            if (sgn > 0) tn_face(&cs, d, -1.0, fa[0], 1.0, 1.0 / jp[0][3 * d], sc, ref);
            else tn_face(&cs, d, -1.0, fa[1], -1.0, 1.0 / jp[1][3 * d], sc, ref);
          }
          double rel = tn_rel_diff(out, ref, np);
          TEST_CHECK(rel < tn_tol);
          TEST_MSG("%s boundary_surf%c %dx%dv (edge %+d, w_v %+d) vs definition: rel diff %.2e", tn_rep_name[r],
            "xyz"[d], cdim, vdim, edge, sgn, rel);
        }
      }
    }
    // Force kernels (independent of the Hamiltonian representation).
    struct tn_case cs;
    tn_case_init(&cs, TN_DENSE, cdim, vdim);
    int np = cs.pbasis.num_basis;
    double a[96], jp[9], out[216], ref[216];
    for (int k = 0; k < 96; ++k) a[k] = tn_rnd(-1.0, 1.0);
    for (int k = 0; k < 9; ++k) jp[k] = tn_rnd(0.5, 1.5);
    for (int kind = 1; kind <= 3; ++kind) {
      if ((kind == 1 && !tn_evol[i]) || (kind == 2 && !tn_phivol[i]) || (kind == 3 && !tn_radvol[i])) continue;
      for (int k = 0; k < np; ++k) out[k] = ref[k] = 0.0;
      if (kind == 1) tn_evol[i](cs.w, cs.dxv, cs.jacob_vel, a, cs.f, out);
      else if (kind == 2) tn_phivol[i](cs.w, cs.dxv, jp, cs.jacob_vel, a, cs.f, out);
      else tn_radvol[i](cs.w, cs.dxv, cs.jacob_vel, a, cs.f, out);
      for (int d = 0; d < (kind == 2 ? nd : vdim); ++d) {
        struct tn_sctx c = { .cs = &cs, .kind = kind, .d = d, .a = a, .out = ref,
          .scale = (2.0 / cs.dxv[cdim + d]) * (kind == 2 ? (2.0 / cs.dxv[d]) / jp[3 * d] : 1.0) };
        tn_quad(&cs, -1, 0.0, tn_sint, &c);
      }
      double rel = tn_rel_diff(out, ref, np);
      TEST_CHECK(rel < tn_tol);
      TEST_MSG("%s vol %dx%dv vs definition: rel diff %.2e", kind == 1 ? "E" : kind == 2 ? "phi" : "rad", cdim,
        vdim, rel);
    }
  }
}

// Canonical-PB phase-space Hamiltonian volume kernels of the tensor p=1 basis.
static const int tn_canon_cv[7][2] = { {1, 1}, {2, 1}, {1, 2}, {2, 2}, {1, 3}, {2, 3}, {3, 3} };
static const tn_vol_t tn_canon_vol[7] = {
  vlasov_hamil_phase_vol_1x1v_tensor_p1, vlasov_hamil_phase_vol_2x1v_tensor_p1,
  vlasov_hamil_phase_vol_1x2v_tensor_p1, vlasov_hamil_phase_vol_2x2v_tensor_p1,
  vlasov_hamil_phase_vol_1x3v_tensor_p1, vlasov_hamil_phase_vol_2x3v_tensor_p1,
  vlasov_hamil_phase_vol_3x3v_tensor_p1 };

void
test_canonical_phase_vol(void)
{
  srand(4321);
  for (int i = 0; i < 7; ++i) {
    if (!tn_canon_vol[i]) continue;
    int cdim = tn_canon_cv[i][0], vdim = tn_canon_cv[i][1];
    for (int n = 0; n < 3; ++n) {
      struct tn_case cs;
      tn_case_init(&cs, TN_PHASE, cdim, vdim);
      cs.canonical = true;
      double out[216], ref[216];
      tn_call_vol(&cs, tn_canon_vol[i], out);
      tn_ref_vol(&cs, ref);
      double rel = tn_rel_diff(out, ref, cs.pbasis.num_basis);
      TEST_CHECK(rel < tn_tol);
      TEST_MSG("canonical phase %dx%dv volume vs definition: rel diff %.2e", cdim, vdim, rel);
    }
  }
}

#ifdef GKYL_HAVE_CUDA

// Monotone nonuniform velocity map.
static void
tn_eval_vmap(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = xn[0] + 0.05 * xn[0] * xn[0] * xn[0];
}

static struct gkyl_array *
tn_mkarr(bool use_gpu, int ncomp, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
                   gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

static void
tn_fill_rnd(struct gkyl_array *arr)
{
  double *d = (double *)arr->data;
  for (size_t i = 0; i < arr->size * arr->ncomp; ++i) d[i] = tn_rnd(-1.0, 1.0);
}

// Device copy of a host array, or the host array itself.
static struct gkyl_array *
tn_on(bool use_gpu, const struct gkyl_array *host)
{
  struct gkyl_array *arr = tn_mkarr(use_gpu, host->ncomp, host->size);
  gkyl_array_copy(arr, host);
  return arr;
}

static double
tn_arr_rel_diff(const struct gkyl_array *dev, const struct gkyl_array *host)
{
  struct gkyl_array *tmp = gkyl_array_new(GKYL_DOUBLE, dev->ncomp, dev->size);
  gkyl_array_copy(tmp, dev);
  double rel = tn_rel_diff((const double *)tmp->data, (const double *)host->data,
    host->size * host->ncomp);
  gkyl_array_release(tmp);
  return rel;
}

// One triad update (conf and vel surface fluxes, then hyper_dg) of random data on a small
// grid, on CPU and on GPU, compared.
static void
tn_updaters_case(enum tn_rep rep, int cdim, int vdim)
{
  const int pdim = cdim + vdim;
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM], vghost[GKYL_MAX_DIM] = {0};
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    cells[d] = pdim == 6 ? 2 : 3;
    ghost[d] = d < cdim ? 1 : 0;
    lower[d] = d < cdim ? 0.0 : -2.0;
    upper[d] = d < cdim ? 1.0 : 2.0;
  }
  struct gkyl_rect_grid cgrid, vgrid, pgrid;
  struct gkyl_range crange, crange_ext, vrange, vrange_ext, prange, prange_ext;
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&cgrid, ghost, &crange_ext, &crange);
  gkyl_rect_grid_init(&vgrid, vdim, &lower[cdim], &upper[cdim], &cells[cdim]);
  gkyl_create_grid_ranges(&vgrid, vghost, &vrange_ext, &vrange);
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&pgrid, ghost, &prange_ext, &prange);

  struct gkyl_basis cbasis, pbasis, vbasis;
  gkyl_cart_modal_tensor(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
  gkyl_cart_modal_tensor(&vbasis, vdim, 2);

  enum gkyl_model_id model_id = rep == TN_PHASE ? GKYL_MODEL_TRIAD_GR : GKYL_MODEL_TRIAD;
  enum gkyl_hamil_id hamil_id = rep == TN_PHASE ? GKYL_HAMIL_PHASE :
    (rep == TN_SPARSE ? GKYL_HAMIL_VEL_SPARSE : GKYL_HAMIL_VEL_DENSE);
  const struct gkyl_range *hrange = rep == TN_PHASE ? &prange : &vrange;

  // Host inputs: nodal Poisson tensor (volume and face storage), Hamiltonian, f.
  int nq = 1 << cdim, ns = 1 << (cdim - 1);
  struct gkyl_array *pt_vol = gkyl_array_new(GKYL_DOUBLE, nq * tn_npt[vdim], crange_ext.volume);
  struct gkyl_array *pt_surf = gkyl_array_new(GKYL_DOUBLE, cdim * 2 * ns * vdim, crange_ext.volume);
  struct gkyl_array *hamil = rep == TN_PHASE ?
    gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, prange_ext.volume) :
    gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vrange.volume);
  struct gkyl_array *fin = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, prange_ext.volume);
  tn_fill_rnd(pt_vol); tn_fill_rnd(pt_surf); tn_fill_rnd(hamil); tn_fill_rnd(fin);

  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  for (int b = 0; b < vdim; ++b) inp_vmap[b].eval_vmap = tn_eval_vmap;

  struct gkyl_array *rhs[2], *cfl[2], *cflux[2], *vflux[2];
  for (int g = 0; g < 2; ++g) {
    bool use_gpu = g == 1;
    struct gkyl_vlasov_velocity_map *vel_map =
      gkyl_vlasov_velocity_map_new(&vgrid, &vrange, &vbasis, inp_vmap, false, use_gpu);
    struct gkyl_vlasov_position_map *pos_map =
      gkyl_vlasov_position_map_new(&cgrid, &crange, &crange_ext, &cbasis, inp_pmap, use_gpu);
    struct gkyl_array *pt_vol_g = tn_on(use_gpu, pt_vol), *pt_surf_g = tn_on(use_gpu, pt_surf);
    struct gkyl_array *hamil_g = tn_on(use_gpu, hamil), *fin_g = tn_on(use_gpu, fin);
    struct gkyl_array *qmem = tn_mkarr(use_gpu, 8 * cbasis.num_basis, crange_ext.volume);
    struct gkyl_array *pot_tot = tn_mkarr(use_gpu, 4 * cbasis.num_basis, crange_ext.volume);
    struct gkyl_array *rad = tn_mkarr(use_gpu, vdim * vbasis.num_basis, vrange.volume);
    gkyl_array_clear(qmem, 0.0); gkyl_array_clear(pot_tot, 0.0); gkyl_array_clear(rad, 0.0);

    struct gkyl_dg_vlasov_conf_flux_surf_inp inp_conf = {
      .phase_grid = &pgrid, .conf_basis = &cbasis, .phase_basis = &pbasis, .vel_range = &vrange,
      .vel_map = vel_map, .pos_map = pos_map, .hamil_range = hrange, .skip_cell_thresh = 0.0,
      .model_id = model_id, .hamil_id = hamil_id, .use_lo = false, .use_gpu = use_gpu,
    };
    struct gkyl_dg_vlasov_conf_flux_surf *calc_conf = gkyl_dg_vlasov_conf_flux_surf_inew(&inp_conf);
    struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel = {
      .phase_grid = &pgrid, .conf_basis = &cbasis, .phase_basis = &pbasis, .vel_map = vel_map,
      .pos_map = pos_map, .hamil_range = hrange, .skip_cell_thresh = 0.0, .model_id = model_id,
      .hamil_id = hamil_id, .has_E = false, .has_phi = false, .has_B = false, .has_rad = false,
      .use_lo = false, .use_gpu = use_gpu,
    };
    struct gkyl_dg_vlasov_vel_flux_surf *calc_vel = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel);
    cflux[g] = tn_mkarr(use_gpu, cdim * gkyl_dg_vlasov_conf_flux_surf_num_surf_basis(calc_conf),
      prange_ext.volume);
    vflux[g] = tn_mkarr(use_gpu, vdim * gkyl_dg_vlasov_vel_flux_surf_num_surf_basis(calc_vel),
      prange_ext.volume);
    rhs[g] = tn_mkarr(use_gpu, pbasis.num_basis, prange_ext.volume);
    cfl[g] = tn_mkarr(use_gpu, 1, prange_ext.volume);
    gkyl_array_clear(cflux[g], 0.0); gkyl_array_clear(vflux[g], 0.0);
    gkyl_array_clear(rhs[g], 0.0); gkyl_array_clear(cfl[g], 0.0);

    struct gkyl_dg_vlasov_inp inp_eqn = {
      .conf_basis = &cbasis, .phase_basis = &pbasis, .conf_range = &crange, .hamil_range = hrange,
      .phase_range = &prange, .vel_map = vel_map, .pos_map = pos_map, .skip_cell_thresh = 0.0,
      .model_id = model_id, .hamil_id = hamil_id, .has_E = false, .has_phi = false,
      .has_B = false, .has_rad = false, .poisson_tensor_conf = pt_vol_g, .hamil = hamil_g,
      .qmem = qmem, .pot_tot = pot_tot, .conf_flux_surf = cflux[g], .vel_flux_surf = vflux[g],
      .f_no_J = fin_g, .rad = rad, .use_lo = false, .use_gpu = use_gpu,
    };
    struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn);
    int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2 * GKYL_MAX_DIM];
    for (int d = 0; d < pdim; ++d) {
      up_dirs[d] = d;
      zero_flux_flags[d] = zero_flux_flags[d + pdim] = d < cdim ? 0 : 1;
    }
    struct gkyl_hyper_dg *slvr =
      gkyl_hyper_dg_new(&pgrid, &pbasis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

    gkyl_dg_vlasov_conf_flux_surf_advance(calc_conf, &crange, &prange, &prange_ext, pt_surf_g,
      hamil_g, fin_g, cfl[g], cflux[g]);
    gkyl_dg_vlasov_vel_flux_surf_advance(calc_vel, &crange, &prange, pt_vol_g, hamil_g, qmem,
      pot_tot, rad, fin_g, cfl[g], vflux[g]);
    gkyl_hyper_dg_advance(slvr, &prange, fin_g, cfl[g], rhs[g]);

    gkyl_hyper_dg_release(slvr);
    gkyl_dg_eqn_release(eqn);
    gkyl_dg_vlasov_vel_flux_surf_release(calc_vel);
    gkyl_dg_vlasov_conf_flux_surf_release(calc_conf);
    gkyl_array_release(qmem); gkyl_array_release(pot_tot); gkyl_array_release(rad);
    gkyl_array_release(pt_vol_g); gkyl_array_release(pt_surf_g);
    gkyl_array_release(hamil_g); gkyl_array_release(fin_g);
    gkyl_vlasov_position_map_release(pos_map);
    gkyl_vlasov_velocity_map_release(vel_map);
  }

  const char *what[4] = { "conf fluxes", "vel fluxes", "rhs", "cflrate" };
  struct gkyl_array *gpu[4] = { cflux[1], vflux[1], rhs[1], cfl[1] };
  struct gkyl_array *cpu[4] = { cflux[0], vflux[0], rhs[0], cfl[0] };
  for (int i = 0; i < 4; ++i) {
    double rel = tn_arr_rel_diff(gpu[i], cpu[i]);
    TEST_CHECK(rel < 1e-12);
    TEST_MSG("%s %dx%dv updaters, %s GPU vs CPU: rel diff %.2e", tn_rep_name[rep], cdim, vdim,
      what[i], rel);
  }
  for (int g = 0; g < 2; ++g) {
    gkyl_array_release(cflux[g]); gkyl_array_release(vflux[g]);
    gkyl_array_release(rhs[g]); gkyl_array_release(cfl[g]);
  }
  gkyl_array_release(pt_vol); gkyl_array_release(pt_surf);
  gkyl_array_release(hamil); gkyl_array_release(fin);
}

void
test_triad_nodal_updaters_dev(void)
{
  srand(1009);
  for (int r = 0; r < 3; ++r)
    for (int i = 0; i < 6; ++i)
      if (tn_tab[r][i].vol) tn_updaters_case((enum tn_rep)r, tn_cv[i][0], tn_cv[i][1]);
}

#endif

TEST_LIST = {
  { "test_triad_nodal_vol", test_triad_nodal_vol },
  { "test_triad_nodal_reps", test_triad_nodal_reps },
  { "test_triad_nodal_conf_force", test_triad_nodal_conf_force },
  { "test_triad_nodal_vel_force", test_triad_nodal_vel_force },
  { "test_canonical_phase_vol", test_canonical_phase_vol },
  { "test_triad_nodal_bvol", test_triad_nodal_bvol },
  { "test_jacob_vel_factored", test_jacob_vel_factored },
  { "test_mom_factored", test_mom_factored },
  { "test_stream_factored", test_stream_factored },
#ifdef GKYL_HAVE_CUDA
  { "test_triad_nodal_updaters_dev", test_triad_nodal_updaters_dev },
#endif
  { NULL, NULL },
};
