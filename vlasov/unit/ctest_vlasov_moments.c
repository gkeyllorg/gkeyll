// Systematic test of the Vlasov moment calculations, scanned across
// dimensionality, polynomial order, basis type, uniform vs non-uniform
// velocity grids, and CPU vs GPU.
//
// The distribution is a multivariate Gaussian with a prescribed full
// pressure tensor,
//   f(x,v) = n(x)/((2 pi)^(vdim/2) sqrt(det Theta)) exp(-(1/2) Thinv_ij c_i c_j),
// with c_i = v_i - u_i and Theta_ij = P_ij/(m n) symmetric positive-definite,
// so the moments have sharp analytic targets and M2ij = n (u_i u_j +
// Theta_ij) directly probes the prescribed pressure tensor, including its
// off-diagonal thermal correlations.
//
// Two reference tiers:
// (1) Discrete exactness (machine precision): the test computes each moment
//     by Gauss-Legendre quadrature of the *discrete modal Jf* against weights
//     built from the *stored velocity map*, treating the geometry nodally
//     exactly as the kernels do (the Jacobian is known at the quadrature
//     points and is time-invariant, so the nodal treatment is exact by
//     definition, for both the C^0 linear and C^1 cubic representations).
// (2) Analytic (tail-truncation limited): the domain-integrated moments are
//     compared against the continuum analytic mode sums.
#include <math.h>

#include <acutest.h>
#include <assert.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_velocity_map.h>

#define VMAX 8.0 // Mapped velocity extent.

static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                                   gkyl_array_new(GKYL_DOUBLE, nc, size);
#else
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
#endif
  return a;
}

// Quadratic mapping (as in the nonuniform regression tests): computational
// domain [-1,1] -> [-VMAX, VMAX].
static void
eval_quad_vmap(double t, const double *vc, double *vp, void *ctx)
{
  vp[0] = vc[0] < 0.0 ? -VMAX * vc[0] * vc[0] : VMAX * vc[0] * vc[0];
}

// Multivariate Gaussian parameters: drift velocity and a symmetric,
// positive-definite Theta_ij = P_ij/(m n) with non-trivial off-diagonal
// correlations; for vdim < 3 the leading vdim x vdim block is used (its
// principal minors are positive). Velocities are O(1) compared to VMAX = 8
// so the truncated tails are below ~1e-9.
static const double udrift[3] = {0.7, -0.4, 0.25};
static const double theta[3][3] = {{0.8, 0.2, -0.15}, {0.2, 1.1, 0.25}, {-0.15, 0.25, 0.65}};

// Smooth conf-space density modulation (periodic on [-pi,pi]^cdim).
static double
den_x(int cdim, const double *x)
{
  double n = 1.0;
  for (int d = 0; d < cdim; ++d) {
    n *= 1.0 + 0.3 * cos(x[d]);
  }
  return n;
}

// Determinant and inverse of the leading vdim x vdim block of Theta.
static double
theta_det(int vdim)
{
  if (vdim == 1) {
    return theta[0][0];
  }
  if (vdim == 2) {
    return theta[0][0] * theta[1][1] - theta[0][1] * theta[1][0];
  }
  return theta[0][0] * (theta[1][1] * theta[2][2] - theta[1][2] * theta[2][1]) -
         theta[0][1] * (theta[1][0] * theta[2][2] - theta[1][2] * theta[2][0]) +
         theta[0][2] * (theta[1][0] * theta[2][1] - theta[1][1] * theta[2][0]);
}

static void
theta_inv(int vdim, double ti[3][3])
{
  double det = theta_det(vdim);
  if (vdim == 1) {
    ti[0][0] = 1.0 / det;
  } else if (vdim == 2) {
    ti[0][0] = theta[1][1] / det;
    ti[1][1] = theta[0][0] / det;
    ti[0][1] = -theta[0][1] / det;
    ti[1][0] = -theta[1][0] / det;
  } else {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        int i1 = (i + 1) % 3, i2 = (i + 2) % 3, j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        ti[j][i] = (theta[i1][j1] * theta[i2][j2] - theta[i1][j2] * theta[i2][j1]) / det;
      }
    }
  }
}

static double
f_eval(int cdim, int vdim, const double *x, const double *v)
{
  double ti[3][3] = {{0.0}};
  theta_inv(vdim, ti);
  double quad = 0.0;
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      quad += ti[i][j] * (v[i] - udrift[i]) * (v[j] - udrift[j]);
    }
  }
  return den_x(cdim, x) * exp(-0.5 * quad) / sqrt(pow(2.0 * GKYL_PI, vdim) * theta_det(vdim));
}

// Evaluate the per-direction velocity map (stored as a degenerate cubic) at
// logical coordinate z in [-1,1].
static inline double
vmap_eval(const struct gkyl_basis *b1, const double *vmap_c, int d, double z)
{
  double zv[1] = {z};
  return b1->eval_expand(zv, &vmap_c[4 * d]);
}

// Build phase-space quadrature data: logical ordinates and total weights at
// nq points per dimension.
struct quad {
  int nq, pdim;
  struct gkyl_range qrange;
  const double *ord, *wgt;
};

static void
quad_init(struct quad *q, int pdim, int nq)
{
  q->nq = nq;
  q->pdim = pdim;
  int shape[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    shape[d] = nq;
  }
  gkyl_range_init_from_shape(&q->qrange, pdim, shape);
  q->ord = gkyl_gauss_ordinates[nq];
  q->wgt = gkyl_gauss_weights[nq];
}

// The moment weight w(v, J) at one mapped-velocity point. mom_id selects:
// 0: M0 (=1), 1..3: M1i, 4: M2, 5: M2ij (vi*vj packed), 6: M3i, 7: M3ijk.
// num_mom returns the component count.
static int
num_mom_comp(int mom_id, int vdim)
{
  switch (mom_id) {
    case 0:
      return 1;
    case 1:
      return vdim;
    case 2:
      return 1;
    case 3:
      return vdim * (vdim + 1) / 2;
    case 4:
      return vdim;
    case 5:
      return vdim * (vdim + 1) * (vdim + 2) / 6;
    default:
      assert(false);
      return 0;
  }
}

static void
mom_weight(int mom_id, int vdim, const double *v, double *w)
{
  int c = 0;
  switch (mom_id) {
    case 0:
      w[0] = 1.0;
      break;
    case 1:
      for (int d = 0; d < vdim; ++d) {
        w[d] = v[d];
      }
      break;
    case 2: {
      double s = 0.0;
      for (int d = 0; d < vdim; ++d) {
        s += v[d] * v[d];
      }
      w[0] = s;
    } break;
    case 3:
      for (int i = 0; i < vdim; ++i) {
        for (int j = i; j < vdim; ++j) {
          w[c++] = v[i] * v[j];
        }
      }
      break;
    case 4: {
      double s = 0.0;
      for (int d = 0; d < vdim; ++d) {
        s += v[d] * v[d];
      }
      for (int i = 0; i < vdim; ++i) {
        w[i] = v[i] * s;
      }
    } break;
    case 5:
      for (int i = 0; i < vdim; ++i) {
        for (int j = i; j < vdim; ++j) {
          for (int k = j; k < vdim; ++k) {
            w[c++] = v[i] * v[j] * v[k];
          }
        }
      }
      break;
    default:
      assert(false);
  }
}

// Continuum analytic domain-integrated moments of the multivariate Gaussian:
// integral over x of n(x) times the velocity-moment.
static double
conf_den_integral(int cdim, double lo, double up)
{
  // integral of n(x) dx over [lo,up]^cdim; product structure per dim:
  // int (1 + 0.3 cos(x)) dx = L + 0.3*(sin(up) - sin(lo)).
  double I = 1.0;
  for (int d = 0; d < cdim; ++d) {
    I *= (up - lo) + 0.3 * (sin(up) - sin(lo));
  }
  return I;
}

static void
analytic_int_moments(int cdim, int vdim, double clo, double cup, double *intmom)
{
  // int_five_moments components: M0, M1i (vdim), M2.
  double N = conf_den_integral(cdim, clo, cup);
  intmom[0] = N;
  double usq = 0.0;
  for (int d = 0; d < vdim; ++d) {
    intmom[1 + d] = N * udrift[d];
    usq += udrift[d] * udrift[d] + theta[d][d];
  }
  // The energy moment is the Hamiltonian moment, H = v^2/2.
  intmom[1 + vdim] = 0.5 * N * usq;
}

// Project Jf: nodal Gauss-Legendre projection at (p+1) points per dimension,
// evaluating f at the mapped velocity and multiplying by the nodal Jacobian
// (the same construction the LTE projection uses).
static void
project_Jf(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *pbasis, const struct gkyl_basis *b1,
  int cdim, int vdim, int poly_order, const struct gkyl_range *phase_local,
  const struct gkyl_range *vel_local, const struct gkyl_vlasov_velocity_map *vvm,
  struct gkyl_array *Jf_ho
)
{
  // Dense quadrature in configuration space (so the conf-space projection of
  // the smooth density modulation does not limit the analytic tier), and
  // (p+1)-point Gauss-Legendre in velocity space where the Jacobian is known
  // nodally (jacob_vel_gauss), as in the LTE projection.
  int nq_c = 8, nq_v = poly_order + 1;
  struct quad qc, qv;
  quad_init(&qc, cdim, nq_c);
  quad_init(&qv, vdim, nq_v);
  int nb = pbasis->num_basis;
  double bvals[160];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_local, iter.idx);
    long vloc = gkyl_range_idx(vel_local, &iter.idx[cdim]);
    const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
    const double *jg_c = gkyl_array_cfetch(vvm->jacob_vel_gauss_host, vloc);
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    double *fc = gkyl_array_fetch(Jf_ho, ploc);
    for (int kk = 0; kk < nb; ++kk) {
      fc[kk] = 0.0;
    }

    struct gkyl_range_iter citer;
    gkyl_range_iter_init(&citer, &qc.qrange);
    while (gkyl_range_iter_next(&citer)) {
      double eta[GKYL_MAX_DIM], x[GKYL_MAX_DIM], wc = 1.0;
      for (int d = 0; d < cdim; ++d) {
        eta[d] = qc.ord[citer.idx[d]];
        wc *= qc.wgt[citer.idx[d]];
        x[d] = xc[d] + 0.5 * grid->dx[d] * eta[d];
      }
      struct gkyl_range_iter viter;
      gkyl_range_iter_init(&viter, &qv.qrange);
      while (gkyl_range_iter_next(&viter)) {
        double v[GKYL_MAX_DIM], wtot = wc;
        for (int d = 0; d < vdim; ++d) {
          eta[cdim + d] = qv.ord[viter.idx[d]];
          wtot *= qv.wgt[viter.idx[d]];
          v[d] = vmap_eval(b1, vmap_c, d, eta[cdim + d]);
        }
        // velocity quadrature index matching jacob_vel_gauss ordering
        long vqidx = 0;
        for (int d = 0; d < vdim; ++d) {
          vqidx = vqidx * nq_v + viter.idx[d];
        }
        double Jfq = jg_c[vqidx] * f_eval(cdim, vdim, x, v);
        pbasis->eval(eta, bvals);
        for (int kk = 0; kk < nb; ++kk) {
          fc[kk] += wtot * bvals[kk] * Jfq;
        }
      }
    }
  }
}

// Discrete reference moment: quadrature of the modal Jf against the moment
// weight. Plain velocity-weighted moments (M0, M2ij, M3ijk) use dense
// quadrature with weights from the stored map (the integrand is polynomial in
// the logical coordinate, so this is exact). The Hamiltonian moments (M1i,
// M2, M3i) mimic the kernels' construction exactly: weights built from the
// MODAL Hamiltonian (dH/dp = (dH/d eta)*(2/dv)/J with the Jacobian at the
// (p+1) Gauss-Legendre nodes, where the geometry is known nodally and the
// treatment is exact by definition), with (p+1)-point quadrature whose nodes
// coincide with the kernels' nodal 1/J placement.
static void
ref_moment(
  int mom_id, const struct gkyl_rect_grid *grid, const struct gkyl_basis *pbasis,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *vbasis, const struct gkyl_basis *b1,
  int cdim, int vdim, int poly_order, const struct gkyl_range *phase_local,
  const struct gkyl_range *conf_local, const struct gkyl_range *vel_local,
  const struct gkyl_vlasov_velocity_map *vvm, const struct gkyl_array *Jf_ho,
  const struct gkyl_array *hamil_ho, struct gkyl_array *ref_ho
)
{
  int nmom = num_mom_comp(mom_id, vdim);
  bool is_hamil_mom = (mom_id == 1) || (mom_id == 2) || (mom_id == 4);
  // The plain-moment integrands are polynomial per logical dimension with
  // degree at most weight (<= 6, the diagonal of M2ij with the cubic map)
  // plus f (p), so (p+4)-point Gauss-Legendre (exact to degree 2p+7) is
  // exact while keeping the 3v reference cost bounded.
  int nq_ref = is_hamil_mom ? poly_order + 1 : poly_order + 4;
  struct quad q;
  quad_init(&q, cdim + vdim, nq_ref);
  double bvals[1024], cvals[32], w[10];

  gkyl_array_clear(ref_ho, 0.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_local);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(phase_local, iter.idx);
    long vloc = gkyl_range_idx(vel_local, &iter.idx[cdim]);
    long cloc = gkyl_range_idx(conf_local, iter.idx); // conf indices lead
    const double *vmap_c = gkyl_array_cfetch(vvm->vmap_host, vloc);
    const double *jacob_vel_c = gkyl_array_cfetch(vvm->jacob_vel_host, vloc);
    const double *hamil_c = gkyl_array_cfetch(hamil_ho, vloc);
    const double *fc = gkyl_array_cfetch(Jf_ho, ploc);
    double *ref_c = gkyl_array_fetch(ref_ho, cloc);

    // velocity-space volume factor: prod dv_d/2 (the Jacobian is inside Jf).
    double velfac = 1.0;
    for (int d = 0; d < vdim; ++d) {
      velfac *= 0.5 * grid->dx[cdim + d];
    }

    struct gkyl_range_iter qiter;
    gkyl_range_iter_init(&qiter, &q.qrange);
    while (gkyl_range_iter_next(&qiter)) {
      double eta[GKYL_MAX_DIM], v[GKYL_MAX_DIM], wtot = 1.0;
      for (int d = 0; d < cdim + vdim; ++d) {
        eta[d] = q.ord[qiter.idx[d]];
        wtot *= q.wgt[qiter.idx[d]];
      }
      for (int d = 0; d < vdim; ++d) {
        v[d] = vmap_eval(b1, vmap_c, d, eta[cdim + d]);
      }

      if (is_hamil_mom) {
        // Weights from the modal Hamiltonian (velocity basis) at this node,
        // with the per-direction Jacobian at the SAME (p+1) GL node the
        // kernels use for the nodal 1/J.
        double Hq = vbasis->eval_expand(&eta[cdim], hamil_c);
        double dH_dp[GKYL_MAX_DIM];
        for (int d = 0; d < vdim; ++d) {
          double Jd = jacob_vel_c[d * (poly_order + 1) + qiter.idx[cdim + d]];
          dH_dp[d] =
            vbasis->eval_grad_expand(d, &eta[cdim], hamil_c) * (2.0 / grid->dx[cdim + d]) / Jd;
        }
        if (mom_id == 1) {
          for (int d = 0; d < vdim; ++d) {
            w[d] = dH_dp[d];
          }
        } else if (mom_id == 2) {
          w[0] = Hq;
        } else { // M3i = int (dH/dp_i) H f
          for (int d = 0; d < vdim; ++d) {
            w[d] = dH_dp[d] * Hq;
          }
        }
      } else {
        mom_weight(mom_id, vdim, v, w);
      }

      double fq = pbasis->eval_expand(eta, fc);
      cbasis->eval(eta, cvals); // conf basis at the conf part of eta
      for (int m = 0; m < nmom; ++m) {
        for (int kk = 0; kk < cbasis->num_basis; ++kk) {
          ref_c[m * cbasis->num_basis + kk] += wtot * cvals[kk] * w[m] * fq * velfac;
        }
      }
    }
  }
}

static const enum gkyl_distribution_moments mom_enum[] = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1,
                                                          GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2IJ,
                                                          GKYL_F_MOMENT_M3, GKYL_F_MOMENT_M3IJK};
static const char *mom_names[] = {"M0", "M1i", "M2", "M2ij", "M3i", "M3ijk"};

static void
test_moments(int cdim, int vdim, int poly_order, bool use_tensor, bool use_nonuniform, bool use_gpu)
{
  int pdim = cdim + vdim;
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = -GKYL_PI;
    upper[d] = GKYL_PI;
    cells[d] = 4;
  }
  for (int d = cdim; d < pdim; ++d) {
    if (use_nonuniform) {
      lower[d] = -1.0;
      upper[d] = 1.0;
    } else {
      lower[d] = -VMAX;
      upper[d] = VMAX;
    }
    cells[d] = 16;
  }

  struct gkyl_rect_grid grid, conf_grid, vel_grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vel_grid, vdim, &lower[cdim], &upper[cdim], &cells[cdim]);

  struct gkyl_basis pbasis, cbasis, vbasis;
  if (use_tensor && poly_order == 1) {
    // Tensor p=1 IS the hybrid: p=1 tensor conf x p=2 tensor velocity basis.
    gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);
    gkyl_cart_modal_tensor(&cbasis, cdim, poly_order);
    gkyl_cart_modal_tensor(&vbasis, vdim, 2);
  } else if (use_tensor) {
    gkyl_cart_modal_tensor(&pbasis, pdim, poly_order);
    gkyl_cart_modal_tensor(&cbasis, cdim, poly_order);
    gkyl_cart_modal_tensor(&vbasis, vdim, poly_order);
  } else {
    gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);
    gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
    gkyl_cart_modal_serendip(&vbasis, vdim, poly_order);
  }
  // Velocity-space polynomial order: drives the quadrature counts and the
  // nodal-Jacobian layout in the reference construction (p=2 for the hybrid).
  int vel_poly = vbasis.poly_order;

  int conf_ghost[GKYL_MAX_DIM] = {1, 1, 1};
  struct gkyl_range conf_local, conf_local_ext;
  gkyl_create_grid_ranges(&conf_grid, conf_ghost, &conf_local_ext, &conf_local);

  int vel_ghost[GKYL_MAX_DIM] = {0, 0, 0};
  struct gkyl_range vel_local, vel_local_ext;
  gkyl_create_grid_ranges(&vel_grid, vel_ghost, &vel_local_ext, &vel_local);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = conf_ghost[d];
  }
  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &phase_local_ext, &phase_local);

  // Velocity map: identity (uniform) or quadratic in every direction.
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = {0};
  if (use_nonuniform) {
    for (int d = 0; d < vdim; ++d) {
      inp_vmap[d].eval_vmap = eval_quad_vmap;
    }
  }
  struct gkyl_vlasov_velocity_map *vvm =
    gkyl_vlasov_velocity_map_new(&vel_grid, &vel_local, &vbasis, inp_vmap, false, use_gpu);

  // Hamiltonian (H = v^2/2) on the mapped grid; device-resident when use_gpu
  // (this host/device mismatch was the pre-existing GPU crash in the old test).
  struct gkyl_array *hamil = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  struct gkyl_array *gamma_inv = mkarr(use_gpu, vbasis.num_basis, vel_local.volume);
  gkyl_dg_vlasov_calc_hamil(
    &vel_grid, &vbasis, &vel_local, GKYL_MODEL_DEFAULT, vvm, hamil, gamma_inv, use_gpu
  );

  // Host copy of the modal Hamiltonian for the reference construction.
  struct gkyl_array *hamil_ho = gkyl_array_new(GKYL_DOUBLE, vbasis.num_basis, vel_local.volume);
  gkyl_array_copy(hamil_ho, hamil);

  // 1D cubic basis for evaluating the stored map.
  struct gkyl_basis b1;
  gkyl_cart_modal_tensor(&b1, 1, 3);

  // Project Jf on the host; copy to device if needed.
  struct gkyl_array *Jf_ho = gkyl_array_new(GKYL_DOUBLE, pbasis.num_basis, phase_local_ext.volume);
  project_Jf(&grid, &pbasis, &b1, cdim, vdim, vel_poly, &phase_local, &vel_local, vvm, Jf_ho);
  struct gkyl_array *Jf = use_gpu ? mkarr(use_gpu, Jf_ho->ncomp, Jf_ho->size) :
                                    gkyl_array_acquire(Jf_ho);
  gkyl_array_copy(Jf, Jf_ho);

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &cbasis,
    .phase_basis = &pbasis,
    .vel_range = &vel_local,
    .hamil_range = &vel_local,
    .hamil = hamil,
    .model_id = GKYL_MODEL_DEFAULT,
    .hamil_id = gkyl_hamil_id_from_model_id(GKYL_MODEL_DEFAULT),
    .mom_type = GKYL_F_MOMENT_M0,
    .use_gpu = use_gpu,
    .vel_map = vvm,
  };

  // Tier 1: every computable moment vs the discrete quadrature reference.
  int n_mom_types = sizeof(mom_enum) / sizeof(mom_enum[0]);
  for (int mt = 0; mt < n_mom_types; ++mt) {
    // No M3i/M3ijk kernels exist for the tensor basis.
    if (use_tensor && (mt == 4 || mt == 5)) {
      continue;
    }
    inp_mom.mom_type = mom_enum[mt];
    struct gkyl_mom_type *mtype = gkyl_mom_vlasov_inew(&inp_mom);
    gkyl_mom_calc *mcalc = gkyl_mom_calc_new(&grid, mtype, use_gpu);

    int nmom = num_mom_comp(mt, vdim);
    struct gkyl_array *mout = mkarr(use_gpu, nmom * cbasis.num_basis, conf_local_ext.volume);
    struct gkyl_array *mout_ho = gkyl_array_new(GKYL_DOUBLE, mout->ncomp, mout->size);
    struct gkyl_array *ref_ho = gkyl_array_new(GKYL_DOUBLE, mout->ncomp, mout->size);

#ifdef GKYL_HAVE_CUDA
    if (use_gpu) {
      gkyl_mom_calc_advance_cu(mcalc, &phase_local, &conf_local, Jf, mout);
    } else {
      gkyl_mom_calc_advance(mcalc, &phase_local, &conf_local, Jf, mout);
    }
#else
    gkyl_mom_calc_advance(mcalc, &phase_local, &conf_local, Jf, mout);
#endif
    gkyl_array_copy(mout_ho, mout);

    ref_moment(
      mt, &grid, &pbasis, &cbasis, &vbasis, &b1, cdim, vdim, vel_poly, &phase_local, &conf_local,
      &vel_local, vvm, Jf_ho, hamil_ho, ref_ho
    );

    struct gkyl_range_iter citer;
    gkyl_range_iter_init(&citer, &conf_local);
    while (gkyl_range_iter_next(&citer)) {
      long cloc = gkyl_range_idx(&conf_local, citer.idx);
      const double *mc = gkyl_array_cfetch(mout_ho, cloc);
      const double *rc = gkyl_array_cfetch(ref_ho, cloc);
      double cmax = 0.0;
      for (int c = 0; c < nmom * cbasis.num_basis; ++c) {
        cmax = fmax(cmax, fabs(rc[c]));
      }
      for (int c = 0; c < nmom * cbasis.num_basis; ++c) {
        TEST_CHECK(fabs(mc[c] - rc[c]) < 1e-11 * cmax);
        TEST_MSG(
          "%s %dx%dv p%d %s %s%s cell=%d comp=%d: kernel %.15e vs ref %.15e", mom_names[mt], cdim,
          vdim, poly_order, use_tensor ? "tensor" : "ser",
          use_nonuniform ? "nonuniform" : "uniform", use_gpu ? " gpu" : "", citer.idx[0], c, mc[c],
          rc[c]
        );
      }
    }

    // For M2ij, anchor the domain integral against the prescribed pressure
    // tensor: int M2ij dx = N (u_i u_j + Theta_ij), off-diagonals included.
    if (mom_enum[mt] == GKYL_F_MOMENT_M2IJ) {
      double N = conf_den_integral(cdim, -GKYL_PI, GKYL_PI);
      double mtol = poly_order == 1 ? 0.15 : 1e-4;
      double cellvol = 1.0;
      for (int d = 0; d < cdim; ++d) {
        cellvol *= conf_grid.dx[d];
      }
      double psi0 = 1.0 / sqrt(pow(2.0, cdim));
      int c = 0;
      for (int i = 0; i < vdim; ++i) {
        for (int j = i; j < vdim; ++j) {
          double total = 0.0;
          struct gkyl_range_iter it2;
          gkyl_range_iter_init(&it2, &conf_local);
          while (gkyl_range_iter_next(&it2)) {
            const double *mc2 = gkyl_array_cfetch(mout_ho, gkyl_range_idx(&conf_local, it2.idx));
            total += mc2[c * cbasis.num_basis] * psi0 * cellvol;
          }
          double exact_ij = N * (udrift[i] * udrift[j] + theta[i][j]);
          TEST_CHECK(gkyl_compare_double(total / exact_ij, 1.0, mtol));
          TEST_MSG(
            "int M2ij(%d,%d) %dx%dv p%d: %.15e vs %.15e (rel %.2e)", i, j, cdim, vdim, poly_order,
            total, exact_ij, total / exact_ij - 1.0
          );
          c += 1;
        }
      }
    }

    gkyl_array_release(mout);
    gkyl_array_release(mout_ho);
    gkyl_array_release(ref_ho);
    gkyl_mom_calc_release(mcalc);
    gkyl_mom_type_release(mtype);
  }

  // Tier 2: domain-integrated five moments vs the continuum analytics
  // (tail truncation at +/-8 max thermal speeds dominates the error).
  inp_mom.mom_type = GKYL_F_MOMENT_M0M1M2;
  struct gkyl_mom_type *int_t = gkyl_int_mom_vlasov_inew(&inp_mom);
  gkyl_mom_calc *intcalc = gkyl_mom_calc_new(&grid, int_t, use_gpu);
  struct gkyl_array *int_mom = mkarr(use_gpu, vdim + 2, conf_local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(intcalc, &phase_local, &conf_local, Jf, int_mom);
  } else {
    gkyl_mom_calc_advance(intcalc, &phase_local, &conf_local, Jf, int_mom);
  }
#else
  gkyl_mom_calc_advance(intcalc, &phase_local, &conf_local, Jf, int_mom);
#endif

  double red_mom[GKYL_MAX_DIM + 2], red_mom_ho[GKYL_MAX_DIM + 2];
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    double *red_cu = gkyl_cu_malloc((vdim + 2) * sizeof(double));
    gkyl_array_reduce_range(red_cu, int_mom, GKYL_SUM, &conf_local);
    gkyl_cu_memcpy(red_mom, red_cu, (vdim + 2) * sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(red_cu);
  } else {
    gkyl_array_reduce_range(red_mom, int_mom, GKYL_SUM, &conf_local);
  }
#else
  gkyl_array_reduce_range(red_mom, int_mom, GKYL_SUM, &conf_local);
#endif

  double exact[GKYL_MAX_DIM + 2];
  analytic_int_moments(cdim, vdim, -GKYL_PI, GKYL_PI, exact);
  // The M1/M2 components are Hamiltonian moments: p=1 truncates H = v^2/2
  // (and dH/dp) per cell, so they carry an O((dv/vth)^2) discretization error
  // by design at this resolution (the discrete-exactness tier above is the
  // correctness check; this tier anchors the physics). p>=2 is limited by the
  // tail truncation and, on mapped grids, the projection of H(vmap).
  // M0 carries the (p+1)-point projection-quadrature error (~1e-6 at this
  // resolution), the same construction the LTE projection uses.
  double tol2 = poly_order == 1 ? 0.15 :
                (use_tensor && use_nonuniform) ?
                                  1e-3 :
                                  1e-4; // H(vmap) is degree 6 for the cubic map, truncated at p=2.
  for (int c = 0; c < vdim + 2; ++c) {
    TEST_CHECK(gkyl_compare_double(red_mom[c] / exact[c], 1.0, c == 0 ? 1e-5 : tol2));
    TEST_MSG(
      "int_mom %dx%dv p%d comp=%d: computed %.15e vs analytic %.15e (rel %.2e)", cdim, vdim,
      poly_order, c, red_mom[c], exact[c], red_mom[c] / exact[c] - 1.0
    );
  }

  gkyl_array_release(int_mom);
  gkyl_mom_calc_release(intcalc);
  gkyl_mom_type_release(int_t);

  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(Jf);
  gkyl_array_release(Jf_ho);
  gkyl_array_release(hamil_ho);
  gkyl_vlasov_velocity_map_release(vvm);
}

// Full scan: {1x1v,1x2v,1x3v,2x2v} x {ser p1, ser p2, tensor p2} x {uniform,
// nonuniform} x {CPU, GPU}, plus ser p3 (1x1v) and a 3x3v ser p1 smoke case.
// Tensor p1 is the HYBRID (p=1 tensor conf x p=2 tensor velocity).
static void
t_1x1v_ser_p1_uni(void)
{
  test_moments(1, 1, 1, false, false, false);
}
static void
t_1x1v_ser_p1_non(void)
{
  test_moments(1, 1, 1, false, true, false);
}
static void
t_1x1v_ser_p2_uni(void)
{
  test_moments(1, 1, 2, false, false, false);
}
static void
t_1x1v_ser_p2_non(void)
{
  test_moments(1, 1, 2, false, true, false);
}
static void
t_1x1v_ten_p2_uni(void)
{
  test_moments(1, 1, 2, true, false, false);
}
static void
t_1x1v_ten_p2_non(void)
{
  test_moments(1, 1, 2, true, true, false);
}
static void
t_1x2v_ser_p1_uni(void)
{
  test_moments(1, 2, 1, false, false, false);
}
static void
t_1x2v_ser_p1_non(void)
{
  test_moments(1, 2, 1, false, true, false);
}
static void
t_1x2v_ser_p2_uni(void)
{
  test_moments(1, 2, 2, false, false, false);
}
static void
t_1x2v_ser_p2_non(void)
{
  test_moments(1, 2, 2, false, true, false);
}
static void
t_1x2v_ten_p2_uni(void)
{
  test_moments(1, 2, 2, true, false, false);
}
static void
t_1x2v_ten_p2_non(void)
{
  test_moments(1, 2, 2, true, true, false);
}
static void
t_1x3v_ser_p1_uni(void)
{
  test_moments(1, 3, 1, false, false, false);
}
static void
t_1x3v_ser_p1_non(void)
{
  test_moments(1, 3, 1, false, true, false);
}
static void
t_1x3v_ser_p2_uni(void)
{
  test_moments(1, 3, 2, false, false, false);
}
static void
t_1x3v_ser_p2_non(void)
{
  test_moments(1, 3, 2, false, true, false);
}
static void
t_2x2v_ser_p1_uni(void)
{
  test_moments(2, 2, 1, false, false, false);
}
static void
t_2x2v_ser_p1_non(void)
{
  test_moments(2, 2, 1, false, true, false);
}
static void
t_2x2v_ser_p2_uni(void)
{
  test_moments(2, 2, 2, false, false, false);
}
static void
t_2x2v_ser_p2_non(void)
{
  test_moments(2, 2, 2, false, true, false);
}
static void
t_2x2v_ten_p2_uni(void)
{
  test_moments(2, 2, 2, true, false, false);
}
static void
t_2x2v_ten_p2_non(void)
{
  test_moments(2, 2, 2, true, true, false);
}
static void
t_1x2v_ten_p1_uni(void)
{
  test_moments(1, 2, 1, true, false, false);
}
static void
t_1x2v_ten_p1_non(void)
{
  test_moments(1, 2, 1, true, true, false);
}
static void
t_2x2v_ten_p1_uni(void)
{
  test_moments(2, 2, 1, true, false, false);
}
static void
t_2x2v_ten_p1_non(void)
{
  test_moments(2, 2, 1, true, true, false);
}
static void
t_2x3v_ten_p1_uni(void)
{
  test_moments(2, 3, 1, true, false, false);
}
static void
t_2x3v_ten_p1_non(void)
{
  test_moments(2, 3, 1, true, true, false);
}
static void
t_1x1v_ser_p3_non(void)
{
  test_moments(1, 1, 3, false, true, false);
}
static void
t_3x3v_ser_p1_uni(void)
{
  test_moments(3, 3, 1, false, false, false);
}
#ifdef GKYL_HAVE_CUDA
static void
t_1x1v_ser_p1_uni_gpu(void)
{
  test_moments(1, 1, 1, false, false, true);
}
static void
t_1x1v_ser_p1_non_gpu(void)
{
  test_moments(1, 1, 1, false, true, true);
}
static void
t_1x1v_ser_p2_uni_gpu(void)
{
  test_moments(1, 1, 2, false, false, true);
}
static void
t_1x1v_ser_p2_non_gpu(void)
{
  test_moments(1, 1, 2, false, true, true);
}
static void
t_1x1v_ten_p2_uni_gpu(void)
{
  test_moments(1, 1, 2, true, false, true);
}
static void
t_1x1v_ten_p2_non_gpu(void)
{
  test_moments(1, 1, 2, true, true, true);
}
static void
t_1x2v_ser_p1_uni_gpu(void)
{
  test_moments(1, 2, 1, false, false, true);
}
static void
t_1x2v_ser_p1_non_gpu(void)
{
  test_moments(1, 2, 1, false, true, true);
}
static void
t_1x2v_ser_p2_uni_gpu(void)
{
  test_moments(1, 2, 2, false, false, true);
}
static void
t_1x2v_ser_p2_non_gpu(void)
{
  test_moments(1, 2, 2, false, true, true);
}
static void
t_1x2v_ten_p2_uni_gpu(void)
{
  test_moments(1, 2, 2, true, false, true);
}
static void
t_1x2v_ten_p2_non_gpu(void)
{
  test_moments(1, 2, 2, true, true, true);
}
static void
t_1x3v_ser_p1_uni_gpu(void)
{
  test_moments(1, 3, 1, false, false, true);
}
static void
t_1x3v_ser_p1_non_gpu(void)
{
  test_moments(1, 3, 1, false, true, true);
}
static void
t_1x3v_ser_p2_uni_gpu(void)
{
  test_moments(1, 3, 2, false, false, true);
}
static void
t_1x3v_ser_p2_non_gpu(void)
{
  test_moments(1, 3, 2, false, true, true);
}
static void
t_2x2v_ser_p1_uni_gpu(void)
{
  test_moments(2, 2, 1, false, false, true);
}
static void
t_2x2v_ser_p1_non_gpu(void)
{
  test_moments(2, 2, 1, false, true, true);
}
static void
t_2x2v_ser_p2_uni_gpu(void)
{
  test_moments(2, 2, 2, false, false, true);
}
static void
t_2x2v_ser_p2_non_gpu(void)
{
  test_moments(2, 2, 2, false, true, true);
}
static void
t_2x2v_ten_p2_uni_gpu(void)
{
  test_moments(2, 2, 2, true, false, true);
}
static void
t_2x2v_ten_p2_non_gpu(void)
{
  test_moments(2, 2, 2, true, true, true);
}
static void
t_1x1v_ser_p3_non_gpu(void)
{
  test_moments(1, 1, 3, false, true, true);
}
static void
t_3x3v_ser_p1_uni_gpu(void)
{
  test_moments(3, 3, 1, false, false, true);
}
#endif

TEST_LIST = {
  {"mom_1x1v_ser_p1_uniform", t_1x1v_ser_p1_uni},
  {"mom_1x1v_ser_p1_nonuniform", t_1x1v_ser_p1_non},
  {"mom_1x1v_ser_p2_uniform", t_1x1v_ser_p2_uni},
  {"mom_1x1v_ser_p2_nonuniform", t_1x1v_ser_p2_non},
  {"mom_1x1v_tensor_p2_uniform", t_1x1v_ten_p2_uni},
  {"mom_1x1v_tensor_p2_nonuniform", t_1x1v_ten_p2_non},
  {"mom_1x2v_ser_p1_uniform", t_1x2v_ser_p1_uni},
  {"mom_1x2v_ser_p1_nonuniform", t_1x2v_ser_p1_non},
  {"mom_1x2v_ser_p2_uniform", t_1x2v_ser_p2_uni},
  {"mom_1x2v_ser_p2_nonuniform", t_1x2v_ser_p2_non},
  {"mom_1x2v_tensor_p2_uniform", t_1x2v_ten_p2_uni},
  {"mom_1x2v_tensor_p2_nonuniform", t_1x2v_ten_p2_non},
  {"mom_1x3v_ser_p1_uniform", t_1x3v_ser_p1_uni},
  {"mom_1x3v_ser_p1_nonuniform", t_1x3v_ser_p1_non},
  {"mom_1x3v_ser_p2_uniform", t_1x3v_ser_p2_uni},
  {"mom_1x3v_ser_p2_nonuniform", t_1x3v_ser_p2_non},
  {"mom_2x2v_ser_p1_uniform", t_2x2v_ser_p1_uni},
  {"mom_2x2v_ser_p1_nonuniform", t_2x2v_ser_p1_non},
  {"mom_2x2v_ser_p2_uniform", t_2x2v_ser_p2_uni},
  {"mom_2x2v_ser_p2_nonuniform", t_2x2v_ser_p2_non},
  {"mom_2x2v_tensor_p2_uniform", t_2x2v_ten_p2_uni},
  {"mom_2x2v_tensor_p2_nonuniform", t_2x2v_ten_p2_non},
  {"mom_1x2v_tensor_p1_uniform", t_1x2v_ten_p1_uni},
  {"mom_1x2v_tensor_p1_nonuniform", t_1x2v_ten_p1_non},
  {"mom_2x2v_tensor_p1_uniform", t_2x2v_ten_p1_uni},
  {"mom_2x2v_tensor_p1_nonuniform", t_2x2v_ten_p1_non},
#ifdef GKYL_BUILD_VLASOV_HYB_2X3V // optional kernel set, see ./configure --help
  {"mom_2x3v_tensor_p1_uniform", t_2x3v_ten_p1_uni},
  {"mom_2x3v_tensor_p1_nonuniform", t_2x3v_ten_p1_non},
#endif
#ifdef GKYL_HAVE_CUDA
  {"mom_1x1v_ser_p1_uniform_gpu", t_1x1v_ser_p1_uni_gpu},
  {"mom_1x1v_ser_p1_nonuniform_gpu", t_1x1v_ser_p1_non_gpu},
  {"mom_1x1v_ser_p2_uniform_gpu", t_1x1v_ser_p2_uni_gpu},
  {"mom_1x1v_ser_p2_nonuniform_gpu", t_1x1v_ser_p2_non_gpu},
  {"mom_1x1v_tensor_p2_uniform_gpu", t_1x1v_ten_p2_uni_gpu},
  {"mom_1x1v_tensor_p2_nonuniform_gpu", t_1x1v_ten_p2_non_gpu},
  {"mom_1x2v_ser_p1_uniform_gpu", t_1x2v_ser_p1_uni_gpu},
  {"mom_1x2v_ser_p1_nonuniform_gpu", t_1x2v_ser_p1_non_gpu},
  {"mom_1x2v_ser_p2_uniform_gpu", t_1x2v_ser_p2_uni_gpu},
  {"mom_1x2v_ser_p2_nonuniform_gpu", t_1x2v_ser_p2_non_gpu},
  {"mom_1x2v_tensor_p2_uniform_gpu", t_1x2v_ten_p2_uni_gpu},
  {"mom_1x2v_tensor_p2_nonuniform_gpu", t_1x2v_ten_p2_non_gpu},
  {"mom_1x3v_ser_p1_uniform_gpu", t_1x3v_ser_p1_uni_gpu},
  {"mom_1x3v_ser_p1_nonuniform_gpu", t_1x3v_ser_p1_non_gpu},
  {"mom_1x3v_ser_p2_uniform_gpu", t_1x3v_ser_p2_uni_gpu},
  {"mom_1x3v_ser_p2_nonuniform_gpu", t_1x3v_ser_p2_non_gpu},
  {"mom_2x2v_ser_p1_uniform_gpu", t_2x2v_ser_p1_uni_gpu},
  {"mom_2x2v_ser_p1_nonuniform_gpu", t_2x2v_ser_p1_non_gpu},
  {"mom_2x2v_ser_p2_uniform_gpu", t_2x2v_ser_p2_uni_gpu},
  {"mom_2x2v_ser_p2_nonuniform_gpu", t_2x2v_ser_p2_non_gpu},
  {"mom_2x2v_tensor_p2_uniform_gpu", t_2x2v_ten_p2_uni_gpu},
  {"mom_2x2v_tensor_p2_nonuniform_gpu", t_2x2v_ten_p2_non_gpu},
#endif
  {NULL, NULL}
};
