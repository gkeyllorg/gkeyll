#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_dens_mu_cells.h>
#include <gkyl_gk_dens_mu_cells_priv.h>

#include <assert.h>
#include <math.h>

struct gkyl_gk_dens_mu_cells *
gkyl_gk_dens_mu_cells_new(
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, double mass, bool use_gpu
)
{
  struct gkyl_gk_dens_mu_cells *up = gkyl_malloc(sizeof(*up));

  up->cdim = conf_basis->ndim;
  up->pdim = phase_basis->ndim;
  up->num_basis = conf_basis->num_basis;
  up->phase_grid = *phase_grid;
  up->mass = mass;
  up->use_gpu = use_gpu;

  // Only implemented for piecewise linear bases with 2 velocity dimensions.
  int cdim = up->cdim;
  assert(up->pdim == cdim + 2);
  assert(conf_basis->poly_order == 1);

  up->kernel = gk_dens_mu_cells_choose_kernel(cdim);

  // The kernel returns the expansion of the integral of f over vpar in the
  // basis of (configuration space, mu). Integrate it over the mu cell with a
  // 2-point Gauss quadrature in each direction, exact for these integrands.
  struct gkyl_basis cmu_basis;
  gkyl_cart_modal_serendip(&cmu_basis, cdim + 1, conf_basis->poly_order);
  up->num_cmu_basis = cmu_basis.num_basis;
  for (int k = 0; k < up->num_basis; k++) {
    for (int j = 0; j < up->num_cmu_basis; j++) {
      up->mu_weight[k][j] = 0.0;
    }
  }
  double ordinates[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  int num_quad = 1 << (cdim + 1);
  for (int q = 0; q < num_quad; q++) {
    double xc[GKYL_MAX_CDIM + 1];
    for (int d = 0; d < cdim + 1; d++) {
      xc[d] = ordinates[(q >> d) & 1];
    }
    double basis_val[GK_DENS_MU_MAX_BASIS], cmu_basis_val[GK_DENS_MU_MAX_CMU_BASIS];
    conf_basis->eval(xc, basis_val);
    cmu_basis.eval(xc, cmu_basis_val);
    for (int k = 0; k < up->num_basis; k++) {
      for (int j = 0; j < up->num_cmu_basis; j++) {
        up->mu_weight[k][j] += basis_val[k] * cmu_basis_val[j];
      }
    }
  }

  up->on_dev = up;
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gk_dens_mu_cells_dev_new(up);
  }
#endif

  return up;
}

void
gkyl_gk_dens_mu_cells_advance(
  const struct gkyl_gk_dens_mu_cells *up, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *fin, struct gkyl_array *out
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_gk_dens_mu_cells_advance_cu(up, phase_range, conf_range, fin, out);
  }
#endif

  int pdim = up->pdim, num_basis = up->num_basis;
  // Factor from the integral over mu of the vpar-integrated distribution.
  double mu_fac = 2.0 * M_PI / up->mass * 0.5 * up->phase_grid.dx[pdim - 1];

  gkyl_array_clear_range(out, 0.0, conf_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    const double *fin_d = gkyl_array_cfetch(fin, gkyl_range_idx(phase_range, iter.idx));
    double *out_d = gkyl_array_fetch(out, gkyl_range_idx(conf_range, iter.idx));

    double vpar_int[GK_DENS_MU_MAX_CMU_BASIS] = {0.0};
    up->kernel(up->phase_grid.dx, 0, up->mass, 0, fin_d, vpar_int);

    int mu_cell = iter.idx[pdim - 1] - phase_range->lower[pdim - 1];
    double *cell_d = &out_d[mu_cell * num_basis];
    for (int k = 0; k < num_basis; k++) {
      double weighted = 0.0;
      for (int j = 0; j < up->num_cmu_basis; j++) {
        weighted += up->mu_weight[k][j] * vpar_int[j];
      }
      cell_d[k] += mu_fac * weighted;
    }
  }
}

void
gkyl_gk_dens_mu_cells_release(struct gkyl_gk_dens_mu_cells *up)
{
  if (up->use_gpu) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}
