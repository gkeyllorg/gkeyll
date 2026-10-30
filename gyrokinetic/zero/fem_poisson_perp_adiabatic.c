#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_perp_adiabatic.h>
#include <gkyl_fem_poisson_perp_adiabatic_priv.h>
#include <gkyl_mat_priv.h>

#include <assert.h>
#include <math.h>

static inline struct gkyl_array *
mkarr(bool use_gpu, long ncomp, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
                   gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

static void
fem_poisson_perp_adiabatic_fem_parproj_op_none(struct gkyl_array *rhs, void *ctx)
{
}

// Solve H phi = rhs, projecting the rhs first.
static void
fem_poisson_perp_adiabatic_helmholtz_solve(
  struct gkyl_fem_poisson_perp_adiabatic *up, struct gkyl_array *rhs, struct gkyl_array *phi
)
{
  up->fem_parproj_op(rhs, up->fem_parproj_ctx);
  gkyl_fem_poisson_perp_set_rhs(up->helmholtz, rhs);
  gkyl_fem_poisson_perp_solve(up->helmholtz, phi);
}

// 1-norm (maximum absolute column sum) of a host matrix.
static double
fem_poisson_perp_adiabatic_mat_norm1(const struct gkyl_mat *mat)
{
  double norm = 0.0;
  for (int j = 0; j < mat->nc; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < mat->nr; i++) {
      col_sum += fabs(gkyl_mat_get(mat, i, j));
    }
    norm = fmax(norm, col_sum);
  }
  return norm;
}

// Invert the zonal matrix A = I - G (a host matrix, overwritten by its LU factors) with LAPACK
// (LU with partial pivoting), and store A^{-1} in up->zonal_inv (on the device if use_gpu) to be
// applied with BLAS/cuBLAS. Returns the 1-norm condition number of A.
static double
fem_poisson_perp_adiabatic_zonal_inv_new(
  struct gkyl_fem_poisson_perp_adiabatic *up, struct gkyl_mat *zonal_mat, bool use_gpu
)
{
  int m = zonal_mat->nr;
  double zonal_mat_norm = fem_poisson_perp_adiabatic_mat_norm1(zonal_mat);

  struct gkyl_mat *zonal_mat_inv = gkyl_mat_diag(gkyl_mat_new(m, m, 0.0), 1.0);
  long *ipiv = gkyl_malloc(sizeof(long[m]));
  bool status = gkyl_mat_linsolve_lu(zonal_mat, zonal_mat_inv, ipiv);
  assert(status);
  double zonal_cond = zonal_mat_norm * fem_poisson_perp_adiabatic_mat_norm1(zonal_mat_inv);

  up->zonal_inv = gkyl_mat_mm_array_mem_new(m, m, 1.0, 0.0, GKYL_NO_TRANS, GKYL_NO_TRANS, use_gpu);
  gkyl_mat_copy(up->zonal_inv->A, zonal_mat_inv);

  gkyl_free(ipiv);
  gkyl_mat_release(zonal_mat_inv);
  return zonal_cond;
}

struct gkyl_fem_poisson_perp_adiabatic *
gkyl_fem_poisson_perp_adiabatic_inew(const struct gkyl_fem_poisson_perp_adiabatic_inp *inp)
{
  struct gkyl_fem_poisson_perp_adiabatic *up = gkyl_malloc(sizeof(*up));

  bool use_gpu = inp->use_gpu;
  up->basis = inp->basis;
  up->local = *inp->local;
  up->local_x = *inp->local_x;
  up->fsa = inp->fsa;
  up->kJ = gkyl_array_acquire(inp->kJ);
  up->lap_fac = inp->lap_fac;
  up->lap_frac = inp->lap_fac / (1.0 + inp->lap_fac);
  up->fem_parproj_op = inp->fem_parproj_op ? inp->fem_parproj_op :
                                             fem_poisson_perp_adiabatic_fem_parproj_op_none;
  up->fem_parproj_ctx = inp->fem_parproj_ctx;

  int nb = up->basis.num_basis, nb_x = inp->basis_x.num_basis;
  long vol = inp->local_ext->volume;
  up->phi_buff = mkarr(use_gpu, nb, vol);
  up->rhs_buff = mkarr(use_gpu, nb, vol);
  up->phi_rhs0 = mkarr(use_gpu, nb, vol);
  up->psi = mkarr(use_gpu, nb_x, inp->local_x_ext->volume);
  up->psi_rhs0 = mkarr(use_gpu, nb_x, inp->local_x_ext->volume);

  // Helmholtz solver for H = -(1+alpha) div(eps grad) + K. The solver subtracts kSq*phi, so
  // kSq = -K; it acquires its weights.
  struct gkyl_array *eps = mkarr(use_gpu, inp->epsilon->ncomp, inp->epsilon->size);
  gkyl_array_set(eps, 1.0 + up->lap_fac, inp->epsilon);
  struct gkyl_array *kSq = mkarr(use_gpu, nb, vol);
  gkyl_array_set(kSq, -1.0, up->kJ);
  up->helmholtz = gkyl_fem_poisson_perp_new(
    &up->local, inp->grid, up->basis, inp->bcs, inp->bias_line_list, eps, kSq, use_gpu
  );
  gkyl_array_release(eps);
  gkyl_array_release(kSq);

  // Solution with zero rhs but Dirichlet BCs, and its flux-surface average.
  gkyl_array_clear(up->rhs_buff, 0.0);
  fem_poisson_perp_adiabatic_helmholtz_solve(up, up->rhs_buff, up->phi_rhs0);
  gkyl_flux_surf_avg_advance(up->fsa, up->phi_rhs0, up->psi_rhs0);

  // Zonal system A = I - G, G = R H^{-1} K E. Column j of G is R H^{-1} K E e_j, so assembling
  // it requires m Helmholtz solves, one per zonal degree of freedom.
  int m = up->local_x.volume * nb_x;
  up->zonal_rhs = mkarr(use_gpu, m, 1);
  up->zonal_sol = mkarr(use_gpu, m, 1);
  struct gkyl_array *zonal_ho = mkarr(false, m, 1);
  double *zonal_vec_ho = zonal_ho->data;
  struct gkyl_mat *zonal_mat = gkyl_mat_new(m, m, 0.0);
  for (int j = 0; j < m; j++) {
    // psi = e_j.
    gkyl_array_clear(zonal_ho, 0.0);
    zonal_vec_ho[j] = 1.0;
    gkyl_array_copy(up->zonal_sol, zonal_ho);
    gkyl_array_copy_from_buffer(up->psi, up->zonal_sol->data, &up->local_x);

    // Solve H phi = K E e_j into phi_buff, removing the rhs = 0 solution.
    gkyl_flux_surf_avg_inflate(up->fsa, up->psi, up->phi_buff);
    gkyl_dg_mul_op_range(&up->basis, 0, up->rhs_buff, 0, up->kJ, 0, up->phi_buff, &up->local);
    fem_poisson_perp_adiabatic_helmholtz_solve(up, up->rhs_buff, up->phi_buff);
    gkyl_array_accumulate_range(up->phi_buff, -1.0, up->phi_rhs0, &up->local);

    // Column j of A = I - G, with G e_j = R H^{-1} K E e_j.
    gkyl_flux_surf_avg_advance(up->fsa, up->phi_buff, up->psi);
    gkyl_array_copy_to_buffer(up->zonal_rhs->data, up->psi, &up->local_x);
    gkyl_array_copy(zonal_ho, up->zonal_rhs);
    for (int i = 0; i < m; i++) {
      gkyl_mat_set(zonal_mat, i, j, (i == j ? 1.0 : 0.0) - zonal_vec_ho[i]);
    }
  }

  // Invert A once; the advance method only applies A^{-1}.
  up->zonal_cond = fem_poisson_perp_adiabatic_zonal_inv_new(up, zonal_mat, use_gpu);

  gkyl_mat_release(zonal_mat);
  gkyl_array_release(zonal_ho);

  return up;
}

void
gkyl_fem_poisson_perp_adiabatic_advance(
  struct gkyl_fem_poisson_perp_adiabatic *up, struct gkyl_array *rho, struct gkyl_array *phi
)
{
  // Solve H phi1 = rho.
  fem_poisson_perp_adiabatic_helmholtz_solve(up, rho, phi);

  // Zonal rhs <phi1> - beta <phi_rhs0> into psi.
  gkyl_flux_surf_avg_advance(up->fsa, phi, up->psi);
  gkyl_array_accumulate_range(up->psi, -up->lap_frac, up->psi_rhs0, &up->local_x);

  // Solve (I - G) psi = rhs, applying the precomputed (I - G)^{-1} to psi packed into an
  // m-vector.
  gkyl_array_copy_to_buffer(up->zonal_rhs->data, up->psi, &up->local_x);
  gkyl_mat_mm_array(up->zonal_inv, up->zonal_rhs, up->zonal_sol);
  gkyl_array_copy_from_buffer(up->psi, up->zonal_sol->data, &up->local_x);

  // E psi into phi_buff, adding alpha E psi to phi.
  gkyl_flux_surf_avg_inflate(up->fsa, up->psi, up->phi_buff);
  gkyl_array_accumulate_range(phi, up->lap_fac, up->phi_buff, &up->local);

  // Solve H phi2 = K E psi into phi_buff, removing (1+beta) times the rhs = 0 solution.
  gkyl_dg_mul_op_range(&up->basis, 0, up->rhs_buff, 0, up->kJ, 0, up->phi_buff, &up->local);
  fem_poisson_perp_adiabatic_helmholtz_solve(up, up->rhs_buff, up->phi_buff);
  gkyl_array_accumulate_range(up->phi_buff, -(1.0 + up->lap_frac), up->phi_rhs0, &up->local);

  // phi = phi1 + alpha E psi + phi2.
  gkyl_array_accumulate_range(phi, 1.0, up->phi_buff, &up->local);
}

double
gkyl_fem_poisson_perp_adiabatic_cond(const struct gkyl_fem_poisson_perp_adiabatic *up)
{
  return up->zonal_cond;
}

void
gkyl_fem_poisson_perp_adiabatic_release(struct gkyl_fem_poisson_perp_adiabatic *up)
{
  gkyl_array_release(up->kJ);
  gkyl_fem_poisson_perp_release(up->helmholtz);
  gkyl_array_release(up->phi_buff);
  gkyl_array_release(up->rhs_buff);
  gkyl_array_release(up->phi_rhs0);
  gkyl_array_release(up->psi);
  gkyl_array_release(up->psi_rhs0);
  gkyl_array_release(up->zonal_rhs);
  gkyl_array_release(up->zonal_sol);
  gkyl_mat_mm_array_mem_release(up->zonal_inv);
  gkyl_free(up);
}
