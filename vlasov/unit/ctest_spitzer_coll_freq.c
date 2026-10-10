#include <gkyl_array.h>
#include <gkyl_util.h>
#include <acutest.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include "math.h"

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                                   gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
eval_m0s_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx));
}
void
eval_m0r_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx));
}
void
eval_vtsqs_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2. * M_PI;
  fout[0] = -x * x + 1.5 * pow(M_PI, 2);
}
void
eval_vtsqr_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2. * M_PI;
  fout[0] = -x * x + 3.5 * pow(M_PI, 2);
}
void
eval_nu_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2. * M_PI;
  double norm_nu = 1. / 3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_1x(t, xn, vtsqs, ctx);
  eval_m0r_1x(t, xn, m0r, ctx);
  eval_vtsqr_1x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0] / pow(vtsqs[0] + vtsqr[0], 1.5);
}

void
test_1x(int poly_order, bool use_gpu)
{
  double Lx = 2. * M_PI;
  double lower[] = {-Lx / 2.}, upper[] = {Lx / 2.};
  int cells[] = {32};

  double norm_nu = 1. / 3.;

  int ndim = sizeof(lower) / sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqs = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  m0r = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqr = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_ho, *vtsqs_ho, *m0r_ho, *vtsqr_ho;
  if (use_gpu) { // Create device copies
    m0s_ho = mkarr(false, m0s->ncomp, m0s->size);
    vtsqs_ho = mkarr(false, vtsqs->ncomp, vtsqs->size);
    m0r_ho = mkarr(false, m0r->ncomp, m0r->size);
    vtsqr_ho = mkarr(false, vtsqr->ncomp, vtsqr->size);
  } else {
    m0s_ho = gkyl_array_acquire(m0s);
    vtsqs_ho = gkyl_array_acquire(vtsqs);
    m0r_ho = gkyl_array_acquire(m0r);
    vtsqr_ho = gkyl_array_acquire(vtsqr);
  }

  gkyl_proj_on_basis *proj_m0s =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0s_1x, NULL);
  gkyl_proj_on_basis *proj_vtsqs =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqs_1x, NULL);
  gkyl_proj_on_basis *proj_m0r =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0r_1x, NULL);
  gkyl_proj_on_basis *proj_vtsqr =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqr_1x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s_ho);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs_ho);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r_ho);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr_ho);

  // Copy host array to device.
  gkyl_array_copy(m0s, m0s_ho);
  gkyl_array_copy(vtsqs, vtsqs_ho);
  gkyl_array_copy(m0r, m0r_ho);
  gkyl_array_copy(vtsqr, vtsqr_ho);

  // Package moments into a Maxwellian moments array (u =0);
  struct gkyl_array *moms_s = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_s, 1.0, m0s, 0);
  gkyl_array_set_offset(moms_s, 1.0, vtsqs, 2 * basis.num_basis);
  struct gkyl_array *moms_r = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_r, 1.0, m0r, 0);
  gkyl_array_set_offset(moms_r, 1.0, vtsqr, 2 * basis.num_basis);

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_ho;
  nu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nu_ho = use_gpu ? mkarr(false, nu->ncomp, nu->size) : gkyl_array_acquire(nu);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up =
    gkyl_spitzer_coll_freq_new(&basis, poly_order + 1, nufrac, epsilon_0, hbar, use_gpu);

  gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, moms_s, 0., moms_r, 0., norm_nu, nu);
  gkyl_array_copy(nu_ho, nu);

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA_ho;
  nuA_ho = mkarr(false, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_nu_1x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA_ho);

  for (int k = 0; k < cells[0]; k++) {
    int idx[] = {k + 1};
    long linidx = gkyl_range_idx(&local, idx);
    const double *nu_p = gkyl_array_cfetch(nu_ho, linidx);
    const double *nuA_p = gkyl_array_cfetch(nuA_ho, linidx);
    for (int m = 0; m < basis.num_basis; m++) {
      TEST_CHECK(gkyl_compare(nuA_p[m], nu_p[m], 1e-10));
      TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
      TEST_MSG("Produced: %.13e", nu_p[m]);
    }
  }

  gkyl_array_release(m0s);
  gkyl_array_release(vtsqs);
  gkyl_array_release(m0r);
  gkyl_array_release(vtsqr);
  gkyl_array_release(nu);
  gkyl_array_release(nu_ho);
  gkyl_array_release(nuA_ho);
  gkyl_array_release(m0s_ho);
  gkyl_array_release(vtsqs_ho);
  gkyl_array_release(m0r_ho);
  gkyl_array_release(vtsqr_ho);
  gkyl_array_release(moms_s);
  gkyl_array_release(moms_r);

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void
eval_m0s_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx)) * (y + Ly);
}
void
eval_m0r_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx)) * (y + Ly);
}
void
eval_vtsqs_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  fout[0] = (-x * x + 1.5 * pow(M_PI, 2)) * 0.5 * (1. + cos(0.5 * 2. * M_PI * y / Ly));
}
void
eval_vtsqr_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  fout[0] = (-x * x + 3.5 * pow(M_PI, 2)) * 0.5 * (1. + cos(0.5 * 2. * M_PI * y / Ly));
}
void
eval_nu_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  double norm_nu = 1. / 3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_2x(t, xn, vtsqs, ctx);
  eval_m0r_2x(t, xn, m0r, ctx);
  eval_vtsqr_2x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0] / pow(vtsqs[0] + vtsqr[0], 1.5);
}

void
test_2x(int poly_order, bool use_gpu)
{
  double Lx = 2. * M_PI, Ly = 8. * M_PI;
  double lower[] = {-Lx / 2., -Ly / 2.}, upper[] = {Lx / 2., Ly / 2.};
  int cells[] = {32, 128};

  double norm_nu = 1. / 3.;

  int ndim = sizeof(lower) / sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqs = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  m0r = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqr = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_ho, *vtsqs_ho, *m0r_ho, *vtsqr_ho;
  if (use_gpu) { // Create device copies
    m0s_ho = mkarr(false, m0s->ncomp, m0s->size);
    vtsqs_ho = mkarr(false, vtsqs->ncomp, vtsqs->size);
    m0r_ho = mkarr(false, m0r->ncomp, m0r->size);
    vtsqr_ho = mkarr(false, vtsqr->ncomp, vtsqr->size);
  } else {
    m0s_ho = gkyl_array_acquire(m0s);
    vtsqs_ho = gkyl_array_acquire(vtsqs);
    m0r_ho = gkyl_array_acquire(m0r);
    vtsqr_ho = gkyl_array_acquire(vtsqr);
  }

  gkyl_proj_on_basis *proj_m0s =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0s_2x, NULL);
  gkyl_proj_on_basis *proj_vtsqs =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqs_2x, NULL);
  gkyl_proj_on_basis *proj_m0r =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0r_2x, NULL);
  gkyl_proj_on_basis *proj_vtsqr =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqr_2x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s_ho);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs_ho);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r_ho);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr_ho);

  // Copy host array to device.
  gkyl_array_copy(m0s, m0s_ho);
  gkyl_array_copy(vtsqs, vtsqs_ho);
  gkyl_array_copy(m0r, m0r_ho);
  gkyl_array_copy(vtsqr, vtsqr_ho);

  // Package moments into a Maxwellian moments array (u =0);
  struct gkyl_array *moms_s = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_s, 1.0, m0s, 0);
  gkyl_array_set_offset(moms_s, 1.0, vtsqs, 2 * basis.num_basis);
  struct gkyl_array *moms_r = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_r, 1.0, m0r, 0);
  gkyl_array_set_offset(moms_r, 1.0, vtsqr, 2 * basis.num_basis);

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_ho;
  nu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nu_ho = use_gpu ? mkarr(false, nu->ncomp, nu->size) : gkyl_array_acquire(nu);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up =
    gkyl_spitzer_coll_freq_new(&basis, poly_order + 1, nufrac, epsilon_0, hbar, use_gpu);

  gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, moms_s, 0., moms_r, 0., norm_nu, nu);
  gkyl_array_copy(nu_ho, nu);

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA_ho;
  nuA_ho = mkarr(false, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_nu_2x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA_ho);

  for (int j = 0; j < cells[0]; j++) {
    for (int k = 0; k < cells[1]; k++) {
      int idx[] = {j + 1, k + 1};
      long linidx = gkyl_range_idx(&local, idx);
      const double *nu_p = gkyl_array_cfetch(nu_ho, linidx);
      const double *nuA_p = gkyl_array_cfetch(nuA_ho, linidx);
      for (int m = 0; m < basis.num_basis; m++) {
        TEST_CHECK(gkyl_compare(nuA_p[m], nu_p[m], 1e-6));
        TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
        TEST_MSG("Produced: %.13e", nu_p[m]);
      }
    }
  }

  gkyl_array_release(m0s);
  gkyl_array_release(vtsqs);
  gkyl_array_release(m0r);
  gkyl_array_release(vtsqr);
  gkyl_array_release(nu);
  gkyl_array_release(nu_ho);
  gkyl_array_release(nuA_ho);
  gkyl_array_release(m0s_ho);
  gkyl_array_release(vtsqs_ho);
  gkyl_array_release(m0r_ho);
  gkyl_array_release(vtsqr_ho);
  gkyl_array_release(moms_s);
  gkyl_array_release(moms_r);

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void
eval_m0s_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx)) * (y + Ly) * (Lz / 2. - 0.5 * fabs(z));
}
void
eval_m0r_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  fout[0] = 0.5 * (1. + cos(0.5 * 2. * M_PI * x / Lx)) * (y + Ly) * (Lz / 2. - 0.5 * fabs(z));
  ;
}
void
eval_vtsqs_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  fout[0] =
    (-x * x + 1.5 * pow(M_PI, 2)) * 0.5 * (1. + cos(0.5 * 2. * M_PI * y / Ly)) * (1.5 + tanh(z));
}
void
eval_vtsqr_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  fout[0] =
    (-x * x + 3.5 * pow(M_PI, 2)) * 0.5 * (1. + cos(0.5 * 2. * M_PI * y / Ly)) * (1.5 + tanh(z));
  ;
}
void
eval_nu_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  double norm_nu = 1. / 3.;
  double vtsqs[1], m0r[1], vtsqr[1];
  eval_vtsqs_3x(t, xn, vtsqs, ctx);
  eval_m0r_3x(t, xn, m0r, ctx);
  eval_vtsqr_3x(t, xn, vtsqr, ctx);
  fout[0] = norm_nu * m0r[0] / pow(vtsqs[0] + vtsqr[0], 1.5);
}

void
test_3x(int poly_order, bool use_gpu)
{
  double Lx = 2. * M_PI, Ly = 8. * M_PI, Lz = 4. * M_PI;
  double lower[] = {-Lx / 2., -Ly / 2., -Lz / 2.}, upper[] = {Lx / 2., Ly / 2., Lz / 2.};
  int cells[] = {32, 48, 48};

  double norm_nu = 1. / 3.;

  int ndim = sizeof(lower) / sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext; // Local, local-ext phase-space ranges.
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create density and v_t^2 arrays.
  struct gkyl_array *m0s, *vtsqs, *m0r, *vtsqr;
  m0s = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqs = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  m0r = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsqr = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *m0s_ho, *vtsqs_ho, *m0r_ho, *vtsqr_ho;
  if (use_gpu) { // Create device copies
    m0s_ho = mkarr(false, m0s->ncomp, m0s->size);
    vtsqs_ho = mkarr(false, vtsqs->ncomp, vtsqs->size);
    m0r_ho = mkarr(false, m0r->ncomp, m0r->size);
    vtsqr_ho = mkarr(false, vtsqr->ncomp, vtsqr->size);
  } else {
    m0s_ho = gkyl_array_acquire(m0s);
    vtsqs_ho = gkyl_array_acquire(vtsqs);
    m0r_ho = gkyl_array_acquire(m0r);
    vtsqr_ho = gkyl_array_acquire(vtsqr);
  }

  gkyl_proj_on_basis *proj_m0s =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0s_3x, NULL);
  gkyl_proj_on_basis *proj_vtsqs =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqs_3x, NULL);
  gkyl_proj_on_basis *proj_m0r =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0r_3x, NULL);
  gkyl_proj_on_basis *proj_vtsqr =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsqr_3x, NULL);

  gkyl_proj_on_basis_advance(proj_m0s, 0.0, &local, m0s_ho);
  gkyl_proj_on_basis_advance(proj_vtsqs, 0.0, &local, vtsqs_ho);
  gkyl_proj_on_basis_advance(proj_m0r, 0.0, &local, m0r_ho);
  gkyl_proj_on_basis_advance(proj_vtsqr, 0.0, &local, vtsqr_ho);

  // Copy host array to device.
  gkyl_array_copy(m0s, m0s_ho);
  gkyl_array_copy(vtsqs, vtsqs_ho);
  gkyl_array_copy(m0r, m0r_ho);
  gkyl_array_copy(vtsqr, vtsqr_ho);

  // Package moments into a Maxwellian moments array (u =0);
  struct gkyl_array *moms_s = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_s, 1.0, m0s, 0);
  gkyl_array_set_offset(moms_s, 1.0, vtsqs, 2 * basis.num_basis);
  struct gkyl_array *moms_r = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_r, 1.0, m0r, 0);
  gkyl_array_set_offset(moms_r, 1.0, vtsqr, 2 * basis.num_basis);

  // Create collision frequency array.
  struct gkyl_array *nu, *nu_ho;
  nu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nu_ho = use_gpu ? mkarr(false, nu->ncomp, nu->size) : gkyl_array_acquire(nu);

  // Create Spitzer collision frequency updater.
  double nufrac = 1., epsilon_0 = 1., hbar = 1.;
  gkyl_spitzer_coll_freq *spitz_up =
    gkyl_spitzer_coll_freq_new(&basis, poly_order + 1, nufrac, epsilon_0, hbar, use_gpu);

  gkyl_spitzer_coll_freq_advance_normnu(spitz_up, &local, moms_s, 0., moms_r, 0., norm_nu, nu);
  gkyl_array_copy(nu_ho, nu);

  // Project expected collision frequency and compare.
  struct gkyl_array *nuA_ho;
  nuA_ho = mkarr(false, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_nu =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_nu_3x, NULL);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &local, nuA_ho);

  for (int i = 0; i < cells[0]; i++) {
    for (int j = 0; j < cells[1]; j++) {
      for (int k = 0; k < cells[2]; k++) {
        int idx[] = {i + 1, j + 1, k + 1};
        long linidx = gkyl_range_idx(&local, idx);
        const double *nu_p = gkyl_array_cfetch(nu_ho, linidx);
        const double *nuA_p = gkyl_array_cfetch(nuA_ho, linidx);
        for (int m = 0; m < basis.num_basis; m++) {
          TEST_CHECK(gkyl_compare(nuA_p[m], nu_p[m], 1e-4));
          TEST_MSG("Expected: %.13e in cell (%d)", nuA_p[m], idx[0]);
          TEST_MSG("Produced: %.13e", nu_p[m]);
        }
      }
    }
  }

  gkyl_array_release(m0s);
  gkyl_array_release(vtsqs);
  gkyl_array_release(m0r);
  gkyl_array_release(vtsqr);
  gkyl_array_release(nu);
  gkyl_array_release(nu_ho);
  gkyl_array_release(nuA_ho);
  gkyl_array_release(m0s_ho);
  gkyl_array_release(vtsqs_ho);
  gkyl_array_release(m0r_ho);
  gkyl_array_release(vtsqr_ho);
  gkyl_array_release(moms_s);
  gkyl_array_release(moms_r);

  gkyl_proj_on_basis_release(proj_m0s);
  gkyl_proj_on_basis_release(proj_vtsqs);
  gkyl_proj_on_basis_release(proj_m0r);
  gkyl_proj_on_basis_release(proj_vtsqr);
  gkyl_proj_on_basis_release(proj_nu);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void
eval_bmag_phys_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0; // Tesla.
}
void
eval_m0_phys_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0e19; // m^-3.
}
void
eval_vtsq_elc_phys_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 100.0 * GKYL_ELEMENTARY_CHARGE / GKYL_ELECTRON_MASS; // Te/me, Te=100 eV.
}
void
eval_vtsq_ion_phys_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double md = GKYL_PROTON_MASS * 2.01410177811; // Deuterium mass.
  fout[0] = 100.0 * GKYL_ELEMENTARY_CHARGE / md; // Ti/md, Ti=100 eV.
}

void
test_physical_1x(bool use_gpu)
{
  // Test using physical SI units.
  // Expected collision frequencies computed independently in Python using the
  // same formula as coulomb_log and gkyl_spitzer_coll_freq_advance in
  // spitzer_coll_freq.c, for n=1e19 m^-3, Te=Ti=100 eV, B=1 T, deuterium ions.
  double nu_ei_expected = 4.46772581565451401e+05;
  double nu_ie_expected = 1.20808171444521520e+02;

  double qe = -GKYL_ELEMENTARY_CHARGE, qd = GKYL_ELEMENTARY_CHARGE;
  double me = GKYL_ELECTRON_MASS, md = GKYL_PROTON_MASS * 2.01410177811;
  double eps0 = GKYL_EPSILON0;
  double hbar = GKYL_PLANCKS_CONSTANT_H / (2.0 * M_PI);

  double Lx = 1.0;
  double lower[] = {0.}, upper[] = {Lx};
  int cells[] = {4};
  int poly_order = 1;
  int ndim = 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *bmag, *m0, *vtsq_elc, *vtsq_ion;
  bmag = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  m0 = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsq_elc = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vtsq_ion = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *bmag_ho, *m0_ho, *vtsq_elc_ho, *vtsq_ion_ho;
  if (use_gpu) {
    bmag_ho = mkarr(false, bmag->ncomp, bmag->size);
    m0_ho = mkarr(false, m0->ncomp, m0->size);
    vtsq_elc_ho = mkarr(false, vtsq_elc->ncomp, vtsq_elc->size);
    vtsq_ion_ho = mkarr(false, vtsq_ion->ncomp, vtsq_ion->size);
  } else {
    bmag_ho = gkyl_array_acquire(bmag);
    m0_ho = gkyl_array_acquire(m0);
    vtsq_elc_ho = gkyl_array_acquire(vtsq_elc);
    vtsq_ion_ho = gkyl_array_acquire(vtsq_ion);
  }

  gkyl_proj_on_basis *proj_bmag =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_bmag_phys_1x, NULL);
  gkyl_proj_on_basis *proj_m0 =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_m0_phys_1x, NULL);
  gkyl_proj_on_basis *proj_vtsq_elc =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsq_elc_phys_1x, NULL);
  gkyl_proj_on_basis *proj_vtsq_ion =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, eval_vtsq_ion_phys_1x, NULL);

  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &local, bmag_ho);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &local, m0_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_elc, 0.0, &local, vtsq_elc_ho);
  gkyl_proj_on_basis_advance(proj_vtsq_ion, 0.0, &local, vtsq_ion_ho);

  gkyl_array_copy(bmag, bmag_ho);
  gkyl_array_copy(m0, m0_ho);
  gkyl_array_copy(vtsq_elc, vtsq_elc_ho);
  gkyl_array_copy(vtsq_ion, vtsq_ion_ho);

  // Package moments into Maxwellian moments arrays (n, u=0, vtSq).
  struct gkyl_array *moms_elc = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, vtsq_elc, 2 * basis.num_basis);
  struct gkyl_array *moms_ion = mkarr(use_gpu, 3 * basis.num_basis, local_ext.volume);
  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, vtsq_ion, 2 * basis.num_basis);

  struct gkyl_array *nu_ei, *nu_ie, *nu_ei_ho, *nu_ie_ho;
  nu_ei = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nu_ie = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nu_ei_ho = use_gpu ? mkarr(false, nu_ei->ncomp, nu_ei->size) : gkyl_array_acquire(nu_ei);
  nu_ie_ho = use_gpu ? mkarr(false, nu_ie->ncomp, nu_ie->size) : gkyl_array_acquire(nu_ie);

  gkyl_spitzer_coll_freq *spitz_up =
    gkyl_spitzer_coll_freq_new(&basis, poly_order + 1, 1.0, eps0, hbar, use_gpu);

  gkyl_spitzer_coll_freq_advance(
    spitz_up, &local, bmag, qe, me, moms_elc, 0., qd, md, moms_ion, 0., nu_ei
  );
  gkyl_spitzer_coll_freq_advance(
    spitz_up, &local, bmag, qd, md, moms_ion, 0., qe, me, moms_elc, 0., nu_ie
  );
  gkyl_array_copy(nu_ei_ho, nu_ei);
  gkyl_array_copy(nu_ie_ho, nu_ie);

  // Fields are spatially uniform, so for the constant orthonormal basis
  // function b_0=1/sqrt(2^ndim), the cell-average coefficient equals the
  // physical value times sqrt(2^ndim) in every cell.
  double cellav_fac = 1. / sqrt(pow(2., ndim));
  for (int k = 0; k < cells[0]; k++) {
    int idx[] = {k + 1};
    long linidx = gkyl_range_idx(&local, idx);
    const double *nu_ei_p = gkyl_array_cfetch(nu_ei_ho, linidx);
    const double *nu_ie_p = gkyl_array_cfetch(nu_ie_ho, linidx);

    TEST_CHECK(gkyl_compare(nu_ei_expected, nu_ei_p[0] * cellav_fac, 1e-8));
    TEST_MSG("Expected nu_ei: %.13e in cell (%d)", nu_ei_expected, idx[0]);
    TEST_MSG("Produced nu_ei: %.13e", nu_ei_p[0] * cellav_fac);

    TEST_CHECK(gkyl_compare(nu_ie_expected, nu_ie_p[0] * cellav_fac, 1e-8));
    TEST_MSG("Expected nu_ie: %.13e in cell (%d)", nu_ie_expected, idx[0]);
    TEST_MSG("Produced nu_ie: %.13e", nu_ie_p[0] * cellav_fac);
  }

  gkyl_array_release(bmag);
  gkyl_array_release(m0);
  gkyl_array_release(vtsq_elc);
  gkyl_array_release(vtsq_ion);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(m0_ho);
  gkyl_array_release(vtsq_elc_ho);
  gkyl_array_release(vtsq_ion_ho);
  gkyl_array_release(moms_elc);
  gkyl_array_release(moms_ion);
  gkyl_array_release(nu_ei);
  gkyl_array_release(nu_ie);
  gkyl_array_release(nu_ei_ho);
  gkyl_array_release(nu_ie_ho);

  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_vtsq_elc);
  gkyl_proj_on_basis_release(proj_vtsq_ion);

  gkyl_spitzer_coll_freq_release(spitz_up);
}

void
test_spitzer_coll_freq_1x_p1_ho()
{
  test_1x(1, false);
}
void
test_spitzer_coll_freq_1x_p2_ho()
{
  test_1x(2, false);
}

void
test_spitzer_coll_freq_2x_p1_ho()
{
  test_2x(1, false);
}
void
test_spitzer_coll_freq_2x_p2_ho()
{
  test_2x(2, false);
}

void
test_spitzer_coll_freq_3x_p1_ho()
{
  test_3x(1, false);
}
void
test_spitzer_coll_freq_3x_p2_ho()
{
  test_3x(2, false);
}

void
test_spitzer_coll_freq_physical_1x_ho()
{
  test_physical_1x(false);
}

#ifdef GKYL_HAVE_CUDA
void
test_spitzer_coll_freq_1x_p1_dev()
{
  test_1x(1, true);
}
void
test_spitzer_coll_freq_1x_p2_dev()
{
  test_1x(2, true);
}

void
test_spitzer_coll_freq_2x_p1_dev()
{
  test_2x(1, true);
}
void
test_spitzer_coll_freq_2x_p2_dev()
{
  test_2x(2, true);
}

void
test_spitzer_coll_freq_3x_p1_dev()
{
  test_3x(1, true);
}
void
test_spitzer_coll_freq_3x_p2_dev()
{
  test_3x(2, true);
}

void
test_spitzer_coll_freq_physical_1x_dev()
{
  test_physical_1x(true);
}
#endif

TEST_LIST = {
  {"test_spitzer_coll_freq_1x_p1_ho", test_spitzer_coll_freq_1x_p1_ho},
  {"test_spitzer_coll_freq_1x_p2_ho", test_spitzer_coll_freq_1x_p2_ho},

  {"test_spitzer_coll_freq_2x_p1_ho", test_spitzer_coll_freq_2x_p1_ho},
  {"test_spitzer_coll_freq_2x_p2_ho", test_spitzer_coll_freq_2x_p2_ho},

  {"test_spitzer_coll_freq_3x_p1_ho", test_spitzer_coll_freq_3x_p1_ho},
  {"test_spitzer_coll_freq_3x_p2_ho", test_spitzer_coll_freq_3x_p2_ho},

  {"test_spitzer_coll_freq_physical_1x_ho", test_spitzer_coll_freq_physical_1x_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_spitzer_coll_freq_1x_p1_dev", test_spitzer_coll_freq_1x_p1_dev},
  {"test_spitzer_coll_freq_1x_p2_dev", test_spitzer_coll_freq_1x_p2_dev},

  {"test_spitzer_coll_freq_2x_p1_dev", test_spitzer_coll_freq_2x_p1_dev},
  {"test_spitzer_coll_freq_2x_p2_dev", test_spitzer_coll_freq_2x_p2_dev},

  {"test_spitzer_coll_freq_3x_p1_dev", test_spitzer_coll_freq_3x_p1_dev},
  {"test_spitzer_coll_freq_3x_p2_dev", test_spitzer_coll_freq_3x_p2_dev},

  {"test_spitzer_coll_freq_physical_1x_dev", test_spitzer_coll_freq_physical_1x_dev},
#endif
  {NULL, NULL}
};
