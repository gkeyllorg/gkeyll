#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_pkpm_lte.h>

#define PKPM_LTE_NUM_MODES 3

struct gkyl_pkpm_lte {
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range;
  int num_quad;
  double eps;
  int max_iter;
  double ordinates[gkyl_gauss_max], weights[gkyl_gauss_max];
  double basis_at_ords[gkyl_gauss_max][PKPM_LTE_NUM_MODES];
};

struct gkyl_pkpm_lte *
gkyl_pkpm_lte_inew(const struct gkyl_pkpm_lte_inp *inp)
{
  assert(inp->vel_grid->ndim == 1);
  assert(inp->num_quad >= PKPM_LTE_NUM_MODES && inp->num_quad <= gkyl_gauss_max);

  struct gkyl_pkpm_lte *up = gkyl_malloc(sizeof(*up));
  up->vel_grid = *inp->vel_grid;
  up->vel_range = *inp->vel_range;
  up->num_quad = inp->num_quad;
  up->eps = inp->eps;
  up->max_iter = inp->max_iter;

  struct gkyl_basis vel_basis;
  gkyl_cart_modal_serendip(&vel_basis, 1, 2);
  memcpy(up->ordinates, gkyl_gauss_ordinates[up->num_quad], sizeof(double[up->num_quad]));
  memcpy(up->weights, gkyl_gauss_weights[up->num_quad], sizeof(double[up->num_quad]));
  for (int q = 0; q < up->num_quad; ++q) {
    vel_basis.eval(&up->ordinates[q], up->basis_at_ords[q]);
  }
  return up;
}

void
gkyl_pkpm_lte_proj(const struct gkyl_pkpm_lte *up, const double *moms, struct gkyl_array *f_lte)
{
  double n = moms[0], u = moms[1], vtsq = moms[2];
  double fac = n / sqrt(2.0 * M_PI * vtsq), half_dw = 0.5 * up->vel_grid.dx[0];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->vel_range);
  while (gkyl_range_iter_next(&iter)) {
    double wc;
    gkyl_rect_grid_cell_center(&up->vel_grid, iter.idx, &wc);
    double *f_d = gkyl_array_fetch(f_lte, gkyl_range_idx(&up->vel_range, iter.idx));
    for (int j = 0; j < PKPM_LTE_NUM_MODES; ++j) {
      f_d[j] = 0.0;
    }
    for (int q = 0; q < up->num_quad; ++q) {
      double dw = wc + half_dw * up->ordinates[q] - u;
      double f_q = fac * exp(-0.5 * dw * dw / vtsq);
      for (int j = 0; j < PKPM_LTE_NUM_MODES; ++j) {
        f_d[j] += up->weights[q] * up->basis_at_ords[q][j] * f_q;
      }
    }
  }
}

void
gkyl_pkpm_lte_moments(const struct gkyl_pkpm_lte *up, const struct gkyl_array *f_lte, double *moms)
{
  // Integrals over the reference cell of 1, eta and eta^2 against the orthonormal p=2 basis.
  const double i0_0 = sqrt(2.0), i1_1 = sqrt(2.0 / 3.0), i2_0 = sqrt(2.0) / 3.0;
  const double i2_2 = 2.0 * sqrt(10.0) / 15.0;
  double h = 0.5 * up->vel_grid.dx[0];

  moms[0] = moms[1] = moms[2] = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->vel_range);
  while (gkyl_range_iter_next(&iter)) {
    double wc;
    gkyl_rect_grid_cell_center(&up->vel_grid, iter.idx, &wc);
    const double *f_d = gkyl_array_cfetch(f_lte, gkyl_range_idx(&up->vel_range, iter.idx));
    double m0 = i0_0 * f_d[0], m1 = i1_1 * f_d[1], m2 = i2_0 * f_d[0] + i2_2 * f_d[2];
    moms[0] += h * m0;
    moms[1] += h * (wc * m0 + h * m1);
    moms[2] += h * (wc * wc * m0 + 2.0 * wc * h * m1 + h * h * m2);
  }
}

// Largest relative error of the moments [M0, M1, M2] against the primitive moments [n, u, vtsq].
static double
pkpm_lte_error(const double *target, const double *moms)
{
  double n = target[0], u = target[1], vtsq = target[2], energy = n * (u * u + vtsq);
  double err0 = fabs(moms[0] - n) / n;
  double err1 = fabs(moms[1] - n * u) / (n * sqrt(vtsq));
  double err2 = fabs(moms[2] - energy) / energy;
  return fmax(err0, fmax(err1, err2));
}

struct gkyl_pkpm_lte_status
gkyl_pkpm_lte_correct(const struct gkyl_pkpm_lte *up, const double *moms, struct gkyl_array *f_lte)
{
  struct gkyl_pkpm_lte_status status = {.converged = false, .num_iter = 0, .error = DBL_MAX};
  double iterate[3] = {moms[0], moms[1], moms[2]}, computed[3];

  while (status.num_iter <= up->max_iter) {
    gkyl_pkpm_lte_proj(up, iterate, f_lte);
    gkyl_pkpm_lte_moments(up, f_lte, computed);
    status.error = pkpm_lte_error(moms, computed);
    if (status.error <= up->eps) {
      status.converged = true;
      break;
    }
    if (status.num_iter == up->max_iter) {
      break;
    }

    // Move the projected primitive moments by the mismatch of the discrete ones.
    double n_h = computed[0], u_h = computed[1] / n_h, vtsq_h = computed[2] / n_h - u_h * u_h;
    double next[3] = {
      iterate[0] + moms[0] - n_h, iterate[1] + moms[1] - u_h, iterate[2] + moms[2] - vtsq_h
    };
    if (next[0] <= 0.0 || next[2] <= 0.0) {
      break;
    }
    memcpy(iterate, next, sizeof(next));
    status.num_iter += 1;
  }
  return status;
}

void
gkyl_pkpm_lte_release(struct gkyl_pkpm_lte *up)
{
  gkyl_free(up);
}
