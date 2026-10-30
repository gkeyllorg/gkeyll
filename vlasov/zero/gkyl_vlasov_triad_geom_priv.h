// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

static inline void
log_to_comp(
  int ndim, const double *eta, const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout
)
{
  for (int d = 0; d < ndim; ++d) {
    xout[d] = 0.5 * dx[d] * eta[d] + xc[d];
  }
}

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, int ndim)
{
  for (int d = 0; d < ndim; d++) {
    xphys[d] = xcomp[d];
  }
}

static void
kernel_metric_1v(const double *cov_tangent_basis, double *GKYL_RESTRICT h_ij)
{
  // Finds the metric
  // cov_tangent_basis: components of the covariant tangent basis
  // h_ij: metric

  h_ij[0] = cov_tangent_basis[0] * cov_tangent_basis[0];
}

static void
kernel_metric_2v(const double *cov_tangent_basis, double *GKYL_RESTRICT h_ij)
{
  // Finds the metric
  // cov_tangent_basis: components of the covariant tangent basis
  // h_ij: metric

  int vdim = 2;
  for (int i = 0; i < vdim; ++i) {
    for (int j = i; j < vdim; ++j) {
      int sym_tensor_index = (i * (2 * vdim - i + 1)) / 2 + (j - i);
      h_ij[sym_tensor_index] = 0.0;
      for (int k = 0; k < vdim; ++k) {
        h_ij[sym_tensor_index] += cov_tangent_basis[i * vdim + k] * cov_tangent_basis[j * vdim + k];
      }
    }
  }
}

static void
kernel_metric_3v(const double *cov_tangent_basis, double *GKYL_RESTRICT h_ij)
{
  // Finds the metric
  // cov_tangent_basis: components of the covariant tangent basis
  // h_ij: metric

  int vdim = 3;
  for (int i = 0; i < vdim; ++i) {
    for (int j = i; j < vdim; ++j) {
      int sym_tensor_index = (i * (2 * vdim - i + 1)) / 2 + (j - i);
      h_ij[sym_tensor_index] = 0.0;
      for (int k = 0; k < vdim; ++k) {
        h_ij[sym_tensor_index] += cov_tangent_basis[i * vdim + k] * cov_tangent_basis[j * vdim + k];
      }
    }
  }
}

static void
kernel_metric_inv_1v(const double *h_ij, double *GKYL_RESTRICT h_ij_inv)
{
  // Finds the inverse of the metric
  // h_ij: metric
  // h_ij_inv: metric inverse

  h_ij_inv[0] = 1.0 / h_ij[0];
}

static void
kernel_metric_inv_2v(const double *h_ij, double *GKYL_RESTRICT h_ij_inv)
{
  // Finds the inverse of the metric
  // h_ij: metric
  // h_ij_inv: metric inverse

  double det_h = h_ij[0] * h_ij[2] - h_ij[1] * h_ij[1];
  h_ij_inv[0] = h_ij[2] / det_h;
  h_ij_inv[1] = -h_ij[1] / det_h;
  h_ij_inv[2] = h_ij[0] / det_h;
}

static void
kernel_metric_inv_3v(const double *h_ij, double *GKYL_RESTRICT h_ij_inv)
{
  // Finds the inverse of the metric
  // h_ij: metric
  // h_ij_inv: metric inverse

  double A = h_ij[3] * h_ij[5] - h_ij[4] * h_ij[4];
  double D = -(h_ij[1] * h_ij[5] - h_ij[2] * h_ij[4]);
  double E = h_ij[0] * h_ij[5] - h_ij[2] * h_ij[2];
  double G = h_ij[1] * h_ij[4] - h_ij[2] * h_ij[3];
  double H = -(h_ij[0] * h_ij[4] - h_ij[2] * h_ij[1]);
  double I = h_ij[0] * h_ij[3] - h_ij[1] * h_ij[1];

  double det_h = (h_ij[0] * A) + (h_ij[1] * D) + (h_ij[2] * G);

  h_ij_inv[0] = A / det_h;
  h_ij_inv[1] = D / det_h;
  h_ij_inv[2] = G / det_h;
  h_ij_inv[3] = E / det_h;
  h_ij_inv[4] = H / det_h;
  h_ij_inv[5] = I / det_h;
}

// NOTE: Only holds for relative vectors/ triads (tetrads need a multiplcation my the lorentzian metric)
static void
kernel_vierbein_inv_1v(const double *vierbein, double *GKYL_RESTRICT vierbein_inv)
{
  // Finds the inverse of the vierbein
  // vierbein: vierbein
  // vierbein_inv: vierbein inverse

  vierbein_inv[0] = 1.0 / vierbein[0];
}

// NOTE: Only holds for relative vectors/ triads (tetrads need a multiplcation my the lorentzian metric)
static void
kernel_vierbein_inv_2v(const double *vierbein, double *GKYL_RESTRICT vierbein_inv)
{
  // Finds the inverse of the vierbein
  // vierbein: vierbein
  // vierbein_inv: vierbein inverse

  const double det_h = vierbein[0] * vierbein[3] - vierbein[1] * vierbein[2];
  vierbein_inv[0] = vierbein[3] / det_h;
  vierbein_inv[1] = -vierbein[1] / det_h;
  vierbein_inv[2] = -vierbein[2] / det_h;
  vierbein_inv[3] = vierbein[0] / det_h;
}

// NOTE: Only holds for relative vectors/ triads (tetrads need a multiplcation my the lorentzian metric)
static void
kernel_vierbein_inv_3v(const double *vierbein, double *GKYL_RESTRICT vierbein_inv)
{
  // Finds the inverse of the vierbein
  // vierbein: vierbein
  // vierbein_inv: vierbein inverse

  const double A = vierbein[0], B = vierbein[1], C = vierbein[2];
  const double D = vierbein[3], E = vierbein[4], F = vierbein[5];
  const double G = vierbein[6], H = vierbein[7], I = vierbein[8];

  // Cofactors
  const double C00 = (E * I - F * H);
  const double C01 = -(D * I - F * G);
  const double C02 = (D * H - E * G);

  const double C10 = -(B * I - C * H);
  const double C11 = (A * I - C * G);
  const double C12 = -(A * H - B * G);

  const double C20 = (B * F - C * E);
  const double C21 = -(A * F - C * D);
  const double C22 = (A * E - B * D);

  // Determinant
  const double det = A * C00 + B * C01 + C * C02;

  // Inverse = transpose(cofactor matrix) / det
  vierbein_inv[0] = C00 / det;
  vierbein_inv[1] = C10 / det;
  vierbein_inv[2] = C20 / det;

  vierbein_inv[3] = C01 / det;
  vierbein_inv[4] = C11 / det;
  vierbein_inv[5] = C21 / det;

  vierbein_inv[6] = C02 / det;
  vierbein_inv[7] = C12 / det;
  vierbein_inv[8] = C22 / det;
}

static void
kernel_metric_det_1v(const double *h_ij, double *GKYL_RESTRICT det_h)
{
  // Finds the determinant of the metric
  // h_ij: metric
  // det_h: (squareroot) metric determinant

  det_h[0] = sqrt(h_ij[0]);
}

static void
kernel_metric_det_2v(const double *h_ij, double *GKYL_RESTRICT det_h)
{
  // Finds the determinant of the metric
  // h_ij: metric
  // det_h: (squareroot) metric determinant

  det_h[0] = sqrt(h_ij[0] * h_ij[2] - h_ij[1] * h_ij[1]);
}

static void
kernel_metric_det_3v(const double *h_ij, double *GKYL_RESTRICT det_h)
{
  // Finds the determinant of the metric
  // h_ij: metric
  // det_h: (squareroot) metric determinant

  double A = h_ij[3] * h_ij[5] - h_ij[4] * h_ij[4];
  double D = -(h_ij[1] * h_ij[5] - h_ij[2] * h_ij[4]);
  double G = h_ij[1] * h_ij[4] - h_ij[2] * h_ij[3];
  det_h[0] = sqrt((h_ij[0] * A) + (h_ij[1] * D) + (h_ij[2] * G));
}

static void
asym_ten_indx(const int vdim, const int p, int *m, int *k)
{
  if (vdim == 1) {
    m[0] = 0;
    k[0] = 0;
  } else if (vdim == 2) {
    // map p = 0,1 to index pairs (m,n) = (0,1)
    static const int pairs[1][2] = {{0, 1}};
    m[0] = pairs[p][0];
    k[0] = pairs[p][1];
  } else if (vdim == 3) {
    // map p = 0,1,2 to index pairs (m,n) = (0,1),(0,2),(1,2)
    static const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    m[0] = pairs[p][0];
    k[0] = pairs[p][1];
  }
}

static void
compute_nu_inv_1v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  double *GKYL_RESTRICT nu_inv
)
{
  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  //        = h^{jn} (\sigma^{i,m} g_{n,m})

  /* Shape of nu_inv.
  +----+
  | 00 |
  +----+
  */

  /* Shape of h_ij_inv
  +----+
  | 00 |
  +----+
  */

  int vdim = 1;

  double h_ij_inv_local[1];
  h_ij_inv_local[0] = h_ij_inv[0];

  // Zero out nu_inv components
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      nu_inv[i * vdim + j] = 0;
    }
  }

  // Set nu_inv components
  // nu^{ij} = h^{jn} (\sigma^{i,m} g_{n,m})
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      for (int n = 0; n < vdim; ++n) {
        for (int m = 0; m < vdim; ++m) {
          nu_inv[i * vdim + j] += h_ij_inv_local[j * vdim + n] *
                                  (triad_basis[i * vdim + m] * cov_tangent_basis[n * vdim + m]);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_1v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  const double *triad_basis_gradient, double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+
  | -- | 00 |
  +----+----+
  | -- | -- |
  +----+----+
  */

  /* Shape of nu_inv.
  +----+
  | 00 |
  +----+
  */

  int vdim = 1;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 0;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  double nu_inv[1] = {0.0};

  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  compute_nu_inv_1v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = nu_inv[0];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(sigma)_j/dx^n \cdot g_i - d(sigma)_j/dx^i \cdot g_n )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k = 0;
          int m = 0;
          asym_ten_indx(vdim, p, &k, &m);

          // Diagonal elements are zero
          double d_sigma_j_dx_n_dot_g_i = 0.0;
          double d_sigma_j_dx_i_dot_g_n = 0.0;
          for (int q = 0; q < vdim; ++q) {
            d_sigma_j_dx_n_dot_g_i += triad_basis_gradient[vdim * vdim * n + j * vdim + q] *
                                      cov_tangent_basis[vdim * i + q];
            d_sigma_j_dx_i_dot_g_n += triad_basis_gradient[vdim * vdim * i + j * vdim + q] *
                                      cov_tangent_basis[vdim * n + q];
          }
          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            nu_inv[k * vdim + n] * nu_inv[m * vdim + i] *
            (d_sigma_j_dx_n_dot_g_i - d_sigma_j_dx_i_dot_g_n);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_vierbein_1v(
  const double *vierbein_inv, const double *vierbein_gradient,
  double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // vierbein_inv: inverse metric
  // vierbein_gradient: vierbein gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+
  | -- | 00 |
  +----+----+
  | -- | -- |
  +----+----+
  */

  /* Shape of vierbein_inv.
  +----+
  | 00 |
  +----+
  */

  int vdim = 1;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 0;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = vierbein_inv[0];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(sigma)_j/dx^n \cdot g_i - d(sigma)_j/dx^i \cdot g_n )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k = 0;
          int m = 0;
          asym_ten_indx(vdim, p, &k, &m);

          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            vierbein_inv[k * vdim + n] * vierbein_inv[m * vdim + i] *
            (vierbein_gradient[vdim * vdim * n + i * vdim + j] -
             vierbein_gradient[vdim * vdim * i + n * vdim + j]);
        }
      }
    }
  }
}

static void
compute_nu_inv_2v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  double *GKYL_RESTRICT nu_inv
)
{
  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  //        = h^{jn} (\sigma^{i,m} g_{n,m})

  /* Shape of nu_inv.
  +----+----+
  | 00 | 01 |
  +----+----+
  | 02 | 03 |
  +----+----+
  */

  /* Shape of h_ij_inv
  +----+----+
  | 00 | 01 |
  +----+----+
  | -- | 02 |
  +----+----+
  */

  int vdim = 2;

  double h_ij_inv_local[4];
  h_ij_inv_local[0] = h_ij_inv[0];
  h_ij_inv_local[1] = h_ij_inv[1];
  h_ij_inv_local[2] = h_ij_inv[1];
  h_ij_inv_local[3] = h_ij_inv[2];

  // Zero out nu_inv components
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      nu_inv[i * vdim + j] = 0;
    }
  }

  // Set nu_inv components
  // nu^{ij} = h^{jn} (\sigma^{i,m} g_{n,m})
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      for (int n = 0; n < vdim; ++n) {
        for (int m = 0; m < vdim; ++m) {
          nu_inv[i * vdim + j] += h_ij_inv_local[j * vdim + n] *
                                  (triad_basis[i * vdim + m] * cov_tangent_basis[n * vdim + m]);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_2v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  const double *triad_basis_gradient, double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+----+----+
  | -- | -- | 00 | 01 |
  +----+----+----+----+
  | -- | -- | 02 | 03 |
 +----+----+----+----+
  | -- | -- | -- |04+j|
  +----+----+----+----+
  | -- | -- | -- | -- |
  +----+----+----+----+
  */

  /* Shape of nu_inv.
  +----+----+
  | 00 | 01 |
  +----+----+
  | 02 | 03 |
  +----+----+
  */

  int vdim = 2;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 1;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  double nu_inv[4] = {0.0, 0.0, 0.0, 0.0};

  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  compute_nu_inv_2v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = nu_inv[0];
  conf_poisson_tensor[1] = nu_inv[2];
  conf_poisson_tensor[2] = nu_inv[1];
  conf_poisson_tensor[3] = nu_inv[3];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(sigma)_j/dx^n \cdot g_i - d(sigma)_j/dx^i \cdot g_n )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k;
          int m;
          asym_ten_indx(vdim, p, &k, &m);

          // Diagonal elements are zero
          double d_sigma_j_dx_n_dot_g_i = 0.0;
          double d_sigma_j_dx_i_dot_g_n = 0.0;
          for (int q = 0; q < vdim; ++q) {
            d_sigma_j_dx_n_dot_g_i += triad_basis_gradient[vdim * vdim * n + j * vdim + q] *
                                      cov_tangent_basis[vdim * i + q];
            d_sigma_j_dx_i_dot_g_n += triad_basis_gradient[vdim * vdim * i + j * vdim + q] *
                                      cov_tangent_basis[vdim * n + q];
          }
          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            nu_inv[k * vdim + n] * nu_inv[m * vdim + i] *
            (d_sigma_j_dx_n_dot_g_i - d_sigma_j_dx_i_dot_g_n);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_vierbein_2v(
  const double *vierbein_inv, const double *vierbein_gradient,
  double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // vierbein_inv: inverse metric
  // vierbein_gradient: vierbein gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+----+----+
  | -- | -- | 00 | 01 |
  +----+----+----+----+
  | -- | -- | 02 | 03 |
 +----+----+----+----+
  | -- | -- | -- |04+j|
  +----+----+----+----+
  | -- | -- | -- | -- |
  +----+----+----+----+
  */

  /* Shape of vierbein_inv.
  +----+----+
  | 00 | 01 |
  +----+----+
  | 02 | 03 |
  +----+----+
  */

  int vdim = 2;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 1;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = vierbein_inv[0];
  conf_poisson_tensor[1] = vierbein_inv[2];
  conf_poisson_tensor[2] = vierbein_inv[1];
  conf_poisson_tensor[3] = vierbein_inv[3];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(e_{ji})/dx^n - d(e_{jn})/dx^i )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k;
          int m;
          asym_ten_indx(vdim, p, &k, &m);

          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            vierbein_inv[k * vdim + n] * vierbein_inv[m * vdim + i] *
            (vierbein_gradient[vdim * vdim * n + i * vdim + j] -
             vierbein_gradient[vdim * vdim * i + n * vdim + j]);
        }
      }
    }
  }
}

static void
compute_nu_inv_3v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  double *GKYL_RESTRICT nu_inv
)
{
  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  //        = h^{jn} (\sigma^{i,m} g_{n,m})

  /* Shape of nu_inv.
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | 03 | 04 | 05 |
  +----+----+----+
  | 06 | 07 | 08 |
  +----+----+----+
  */

  /* Shape of h_ij_inv
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | -- | 03 | 04 |
  +----+----+----+
  | -- | -- | 05 |
  +----+----+----+
  */

  int vdim = 3;

  double h_ij_inv_local[9];
  h_ij_inv_local[0] = h_ij_inv[0];
  h_ij_inv_local[1] = h_ij_inv[1];
  h_ij_inv_local[2] = h_ij_inv[2];
  h_ij_inv_local[3] = h_ij_inv[1];
  h_ij_inv_local[4] = h_ij_inv[3];
  h_ij_inv_local[5] = h_ij_inv[4];
  h_ij_inv_local[6] = h_ij_inv[2];
  h_ij_inv_local[7] = h_ij_inv[4];
  h_ij_inv_local[8] = h_ij_inv[5];

  // Zero out nu_inv components
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      nu_inv[i * vdim + j] = 0;
    }
  }

  // Set nu_inv components
  // nu^{ij} = h^{jn} (\sigma^{i,m} g_{n,m})
  for (int i = 0; i < vdim; ++i) {
    for (int j = 0; j < vdim; ++j) {
      for (int n = 0; n < vdim; ++n) {
        for (int m = 0; m < vdim; ++m) {
          nu_inv[i * vdim + j] += h_ij_inv_local[j * vdim + n] *
                                  (triad_basis[i * vdim + m] * cov_tangent_basis[n * vdim + m]);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_3v(
  const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  const double *triad_basis_gradient, double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+----+----+----+----+
  | -- | -- | -- | 00 | 01 | 02 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 03 | 04 | 05 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 06 | 07 | 08 |
  +----+----+----+----+----+----+
  | -- | -- | -- | -- |09+j|10+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- |11+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- | -- |
  +----+----+----+----+----+----+
  */

  /* Shape of nu_inv.
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | 03 | 04 | 05 |
  +----+----+----+
  | 06 | 07 | 08 |
  +----+----+----+
  */

  int vdim = 3;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 3;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  double nu_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  compute_nu_inv_3v(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = nu_inv[0];
  conf_poisson_tensor[1] = nu_inv[3];
  conf_poisson_tensor[2] = nu_inv[6];
  conf_poisson_tensor[3] = nu_inv[1];
  conf_poisson_tensor[4] = nu_inv[4];
  conf_poisson_tensor[5] = nu_inv[7];
  conf_poisson_tensor[6] = nu_inv[2];
  conf_poisson_tensor[7] = nu_inv[5];
  conf_poisson_tensor[8] = nu_inv[8];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(sigma)_j/dx^n \cdot g_i - d(sigma)_j/dx^i \cdot g_n )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k;
          int m;
          asym_ten_indx(vdim, p, &k, &m);

          // Diagonal elements are zero
          double d_sigma_j_dx_n_dot_g_i = 0.0;
          double d_sigma_j_dx_i_dot_g_n = 0.0;
          for (int q = 0; q < vdim; ++q) {
            d_sigma_j_dx_n_dot_g_i += triad_basis_gradient[vdim * vdim * n + j * vdim + q] *
                                      cov_tangent_basis[vdim * i + q];
            d_sigma_j_dx_i_dot_g_n += triad_basis_gradient[vdim * vdim * i + j * vdim + q] *
                                      cov_tangent_basis[vdim * n + q];
          }
          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            nu_inv[k * vdim + n] * nu_inv[m * vdim + i] *
            (d_sigma_j_dx_n_dot_g_i - d_sigma_j_dx_i_dot_g_n);
        }
      }
    }
  }
}

static void
kernel_conf_poisson_tensor_vierbein_3v(
  const double *vierbein_inv, const double *vierbein_gradient,
  double *GKYL_RESTRICT conf_poisson_tensor
)
{
  // Finds the configuration space components of the Poisson tensor
  // vierbein_inv: inverse metric
  // vierbein_gradient: vierbein gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+----+----+----+----+
  | -- | -- | -- | 00 | 01 | 02 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 03 | 04 | 05 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 06 | 07 | 08 |
  +----+----+----+----+----+----+
  | -- | -- | -- | -- |09+j|10+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- |11+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- | -- |
  +----+----+----+----+----+----+
  */

  /* Shape of vierbein_inv.
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | 03 | 04 | 05 |
  +----+----+----+
  | 06 | 07 | 08 |
  +----+----+----+
  */

  int vdim = 3;
  int pdim = 2 * vdim;
  int nonzero_asym_indices = 3;
  int n_nonzero_indices = vdim * vdim + vdim * nonzero_asym_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i = 0; i < n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  // Top right block: \Pi^{ij}_{xp} = e^{ji} (transpose)
  conf_poisson_tensor[0] = vierbein_inv[0];
  conf_poisson_tensor[1] = vierbein_inv[3];
  conf_poisson_tensor[2] = vierbein_inv[6];
  conf_poisson_tensor[3] = vierbein_inv[1];
  conf_poisson_tensor[4] = vierbein_inv[4];
  conf_poisson_tensor[5] = vierbein_inv[7];
  conf_poisson_tensor[6] = vierbein_inv[2];
  conf_poisson_tensor[7] = vierbein_inv[5];
  conf_poisson_tensor[8] = vierbein_inv[8];

  // Bottom right block
  // Pi^{km}_{j,pp} = \Omega^{km}_j
  // \Omega^{km}_j = \sum_{in} \nu^{kn} \nu^{mi} ( d(e_{ji})/dx^n - d(e_{jn})/dx^i )
  // Sum over i and n, j is and independnant coordinated and k(p), m(p) are also independant coordiantes
  for (int i = 0; i < vdim; ++i) {
    for (int n = 0; n < vdim; ++n) {
      for (int j = 0; j < vdim; ++j) {
        for (int p = 0; p < nonzero_asym_indices; ++p) {
          // index the free component of the antisymmetric tensor
          int k;
          int m;
          asym_ten_indx(vdim, p, &k, &m);

          // Offset by vdim^2 which is the space for the Pi_{xx} block
          // plus additional offsets of nonzero_asym_indices*j for each j
          conf_poisson_tensor[vdim * vdim + nonzero_asym_indices * j + p] +=
            vierbein_inv[k * vdim + n] * vierbein_inv[m * vdim + i] *
            (vierbein_gradient[vdim * vdim * n + i * vdim + j] -
             vierbein_gradient[vdim * vdim * i + n * vdim + j]);
        }
      }
    }
  }
}

typedef void (*metric_t)(const double *pn_cov_tangent_basis, double *GKYL_RESTRICT pn_h_ij);
typedef void (*metric_inv_t)(const double *pn_h_ij, double *GKYL_RESTRICT pn_h_ij_inv);
typedef void (*metric_det_t)(const double *pn_h_ij, double *GKYL_RESTRICT pn_det_h);
typedef void (*vierbein_inv_t)(const double *pn_vierbein, double *GKYL_RESTRICT pn_vierbein_inv);
typedef void (*conf_poisson_tensor_t)(
  const double *pn_h_ij_inv, const double *pn_triad_basis, const double *pn_cov_tangent_basis,
  const double *pn_triad_basis_gradient, double *GKYL_RESTRICT pn_conf_poisson_tensor
);
typedef void (*conf_poisson_tensor_vierbein_t)(
  const double *pn_vierbein_inv, const double *pn_vierbein_gradient,
  double *GKYL_RESTRICT pn_conf_poisson_tensor
);

// for use in kernel tables
typedef struct {
  metric_t kernels[3];
} metric_kern_list;
typedef struct {
  metric_inv_t kernels[3];
} metric_inv_kern_list;
typedef struct {
  metric_det_t kernels[3];
} metric_det_kern_list;
typedef struct {
  conf_poisson_tensor_t kernels[6];
} conf_poisson_tensor_kern_list;

typedef struct {
  vierbein_inv_t kernels[3];
} vierbein_inv_kern_list;
typedef struct {
  conf_poisson_tensor_vierbein_t kernels[6];
} conf_poisson_tensor_vierbein_kern_list;

static const metric_kern_list metric_list[] =
  {{kernel_metric_1v}, {kernel_metric_2v}, {kernel_metric_3v}};

static const metric_inv_kern_list metric_inv_list[] =
  {{kernel_metric_inv_1v}, {kernel_metric_inv_2v}, {kernel_metric_inv_3v}};

static const vierbein_inv_kern_list vierbein_inv_list[] =
  {{kernel_vierbein_inv_1v}, {kernel_vierbein_inv_2v}, {kernel_vierbein_inv_3v}};

static const metric_det_kern_list metric_det_list[] =
  {{kernel_metric_det_1v}, {kernel_metric_det_2v}, {kernel_metric_det_3v}};

static const conf_poisson_tensor_kern_list conf_poisson_tensor_list[] = {
  {kernel_conf_poisson_tensor_1v},
  {kernel_conf_poisson_tensor_2v},
  {kernel_conf_poisson_tensor_3v}
};

static const conf_poisson_tensor_vierbein_kern_list conf_poisson_tensor_vierbein_list[] = {
  {kernel_conf_poisson_tensor_vierbein_1v},
  {kernel_conf_poisson_tensor_vierbein_2v},
  {kernel_conf_poisson_tensor_vierbein_3v}
};

static metric_t
choose_metric_kern(int vdim)
{
  return metric_list[vdim - 1].kernels[0];
}

static metric_inv_t
choose_metric_inv_kern(int vdim)
{
  return metric_inv_list[vdim - 1].kernels[0];
}

static metric_det_t
choose_metric_det_kern(int vdim)
{
  return metric_det_list[vdim - 1].kernels[0];
}

static conf_poisson_tensor_t
choose_conf_poisson_tensor_kern(int vdim)
{
  return conf_poisson_tensor_list[vdim - 1].kernels[0];
}

static vierbein_inv_t
choose_vierbein_inv_kern(int vdim)
{
  return vierbein_inv_list[vdim - 1].kernels[0];
}

static conf_poisson_tensor_vierbein_t
choose_conf_poisson_tensor_vierbein_kern(int vdim)
{
  return conf_poisson_tensor_vierbein_list[vdim - 1].kernels[0];
}