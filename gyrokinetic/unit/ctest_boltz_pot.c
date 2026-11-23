// Tests for the ambipolar and zero-current Boltzmann potential updaters.

#include <acutest.h>
#include <float.h>
#include <gkyl_boltz_pot_ambi.h>
#include <gkyl_boltz_pot_0curr.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_dg_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_boltz_pot_ambi_priv.h>
#include <gkyl_util.h>

// Allocate arrays on the chosen backend.
static struct gkyl_array *
mkarr(bool gpu, int ncomp, long size)
{
  return gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size) :
               gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

struct boltz_test_ctx {
  int cdim, poly_order;
  double nref_fac, shift;
};

// A transverse density represented exactly by both supported bases.
static double
boltz_test_density(const double *x, const struct boltz_test_ctx *ctx)
{
  double density = 3.0;
  if (ctx->cdim > 1) {
    density += 0.2 * x[0];
  }
  if (ctx->cdim > 2) {
    density += 0.1 * x[1] + 0.05 * x[0] * x[1];
  }
  return density;
}

// Prescribe a representable transverse reference potential, including a p=2 mode.
static double
boltz_test_potential(const double *x, const struct boltz_test_ctx *ctx)
{
  double phi = 0.7 + ctx->shift;
  if (ctx->cdim > 1) {
    phi += 0.13 * x[0] + (ctx->poly_order == 2 ? 0.04 * x[0] * x[0] : 0.0);
  }
  if (ctx->cdim > 2) {
    phi -= 0.09 * x[1];
  }
  return phi;
}

// Asymmetric density traces at the lower and upper ends.
static void
boltz_test_m0(double t, const double *x, double *out, void *vctx)
{
  struct boltz_test_ctx *ctx = vctx;
  out[0] = boltz_test_density(x, ctx) * (1.0 - 0.2 * x[ctx->cdim - 1]);
}

// Unequal wall biases with a shared gauge shift.
static void
boltz_test_wall(double t, const double *x, double *out, void *vctx)
{
  struct boltz_test_ctx *ctx = vctx;
  out[0] = 0.15 + 0.3 * x[ctx->cdim - 1] + ctx->shift;
  if (ctx->cdim > 1) {
    out[0] += 0.07 * x[0];
  }
}

// Volume density is below the reference in one cell and above it in the other.
static void
boltz_test_volume_density(double t, const double *x, double *out, void *vctx)
{
  struct boltz_test_ctx *ctx = vctx;
  out[0] = ctx->nref_fac * boltz_test_density(x, ctx) * (x[ctx->cdim - 1] < 0.0 ? 0.5 : 2.0);
}

// Analytic clamped volume potential for the manufactured density.
static void
boltz_test_volume_phi(double t, const double *x, double *out, void *vctx)
{
  struct boltz_test_ctx *ctx = vctx;
  out[0] = boltz_test_potential(x, ctx) + (x[ctx->cdim - 1] < 0.0 ? 0.0 : (3.1 / 2.0) * log(2.0));
}

// Project a manufactured profile on the host for independent analytic comparisons.
static void
boltz_test_project(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis, const struct gkyl_range *range,
  evalf_t eval, struct boltz_test_ctx *ctx, struct gkyl_array *out
)
{
  gkyl_proj_on_basis *proj =
    gkyl_proj_on_basis_new(grid, basis, basis->poly_order + 1, 1, eval, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, range, out);
  gkyl_proj_on_basis_release(proj);
}

// Exercise boundary packing, model normalization, and shared volume response with
// asymmetric ends, transverse variation, both edges, a gauge shift, and the clamp.
static void
check_boltz_model(int cdim, int poly_order, bool gpu, bool ambipolar, double shift)
{
  double lower[3] = {-1, -1, -1}, upper[3] = {1, 1, 1};
  int cells[3] = {2, 2, 2}, ng[3] = {1, 1, 1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  int nb = basis.num_basis, nq = 1;
  for (int d = 0; d < cdim - 1; ++d) {
    nq *= poly_order + 1;
  }

  struct gkyl_range local, ext, skin[2], ghost[2];
  gkyl_create_grid_ranges(&grid, ng, &ext, &local);
  for (int edge = 0; edge < 2; ++edge) {
    gkyl_skin_ghost_ranges(
      &skin[edge], &ghost[edge], cdim - 1, edge ? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE, &ext, ng
    );
  }

  struct boltz_test_ctx ctx = {.cdim = cdim, .poly_order = poly_order, .shift = shift};
  gkyl_boltz_pot_ambi *ambi = gkyl_boltz_pot_ambi_new(&basis, 0.7, -2.0, 3.1, gpu);
  gkyl_boltz_pot_0curr *zero = gkyl_boltz_pot_0curr_new(&basis, 0.7, -2.0, 3.1, gpu);

  struct gkyl_array *m0h = mkarr(false, nb, ext.volume);
  struct gkyl_array *wallh = mkarr(false, nb, ext.volume);
  struct gkyl_array *m0 = mkarr(gpu, nb, ext.volume);
  struct gkyl_array *wall = mkarr(gpu, nb, ext.volume);

  boltz_test_project(&grid, &basis, &local, boltz_test_m0, &ctx, m0h);

  boltz_test_project(&grid, &basis, &local, boltz_test_wall, &ctx, wallh);

  gkyl_array_copy(m0, m0h);
  gkyl_array_copy(wall, wallh);

  struct gkyl_array *end_vals[2], *end_h[2];
  double ordinates[18], weights[9];
  if (cdim > 1) {
    gkyl_ndim_ordinates_weights(cdim - 1, ordinates, weights, poly_order + 1);
  }
  double vte_sided = sqrt(3.1 / (2.0 * M_PI * 0.7)), Te_eV = 3.1 / 2.0;
  for (int edge = 0; edge < 2; ++edge) {
    long size = skin[edge].volume;
    end_vals[edge] = mkarr(gpu, 4 * nq, size);
    end_h[edge] = mkarr(false, 4 * nq, size);
    struct gkyl_array *fh = mkarr(false, nq, size), *kh = mkarr(false, nq, size);
    struct gkyl_array *flux = mkarr(gpu, nq, size), *area = mkarr(gpu, nq, size);

    long cell = 0;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &skin[edge]);
    while (gkyl_range_iter_next(&iter)) {
      double xc[3];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      double *f = gkyl_array_fetch(fh, cell), *k = gkyl_array_fetch(kh, cell);

      for (int q = 0; q < nq; ++q) {
        double x[3] = {0};
        for (int d = 0; d < cdim - 1; ++d) {
          x[d] = xc[d] + 0.5 * grid.dx[d] * ordinates[q * (cdim - 1) + d];
        }
        double n = boltz_test_density(x, &ctx), pref = boltz_test_potential(x, &ctx);
        double k_lo = 1.0 + (cdim > 1 ? 0.07 * x[0] : 0.0), k_up = k_lo + 0.1;
        x[cdim - 1] = -1.0;

        double wlo;
        boltz_test_wall(0, x, &wlo, &ctx);

        x[cdim - 1] = 1.0;
        double wup;
        boltz_test_wall(0, x, &wup, &ctx);

        k[q] = edge ? k_up : k_lo;
        // Choose ion losses that recover the prescribed reference potential.
        if (ambipolar) {
          f[q] =
            vte_sided * k[q] * n * (edge ? 0.8 : 1.2) * exp(((edge ? wup : wlo) - pref) / Te_eV);
        } else {
          // Only the sum of the two losses is constrained in the zero-current model.
          double total =
            vte_sided * n * (k_lo * exp((wlo - pref) / Te_eV) + k_up * exp((wup - pref) / Te_eV));
          f[q] = total * (edge ? 0.63 : 0.37);
        }
      }
      ++cell;
    }
    gkyl_array_copy(flux, fh);
    gkyl_array_copy(area, kh);

    enum gkyl_edge_loc loc = edge ? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE;
    if (ambipolar) {
      gkyl_boltz_pot_ambi_boundary(ambi, loc, &skin[edge], flux, area, m0, wall, end_vals[edge]);
    } else {
      gkyl_boltz_pot_0curr_boundary(zero, loc, &skin[edge], flux, area, m0, wall, end_vals[edge]);
    }

    // Check the packed surface data: flux, projected area, density, and wall potential.
    cell = 0;
    gkyl_array_copy(end_h[edge], end_vals[edge]);
    gkyl_range_iter_init(&iter, &skin[edge]);
    while (gkyl_range_iter_next(&iter)) {
      double xc[3];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      const double *packed = gkyl_array_cfetch(end_h[edge], cell);
      const double *f = gkyl_array_cfetch(fh, cell), *k = gkyl_array_cfetch(kh, cell);

      for (int q = 0; q < nq; ++q) {
        double x[3] = {0};
        for (int d = 0; d < cdim - 1; ++d) {
          x[d] = xc[d] + 0.5 * grid.dx[d] * ordinates[q * (cdim - 1) + d];
        }
        x[cdim - 1] = edge ? 1.0 : -1.0;

        double dens, bias;
        boltz_test_m0(0, x, &dens, &ctx);
        boltz_test_wall(0, x, &bias, &ctx);

        TEST_CHECK(gkyl_compare(packed[q], f[q], 1e-12));
        TEST_CHECK(gkyl_compare(packed[nq + q], k[q], 1e-12));
        TEST_CHECK(gkyl_compare(packed[2 * nq + q], dens, 1e-12));
        TEST_CHECK(gkyl_compare(packed[3 * nq + q], bias, 1e-12));
      }
      ++cell;
    }

    // NULL wall must give exactly zero wall values without changing density.
    if (ambipolar) {
      gkyl_boltz_pot_ambi_boundary(ambi, loc, &skin[edge], flux, area, m0, NULL, end_vals[edge]);
    } else {
      gkyl_boltz_pot_0curr_boundary(zero, loc, &skin[edge], flux, area, m0, NULL, end_vals[edge]);
    }

    struct gkyl_array *no_wall = mkarr(false, 4 * nq, size);
    gkyl_array_copy(no_wall, end_vals[edge]);
    for (long i = 0; i < size; ++i) {
      const double *v = gkyl_array_cfetch(no_wall, i);
      for (int q = 0; q < nq; ++q) {
        TEST_CHECK(v[3 * nq + q] == 0.0);
      }
    }
    gkyl_array_copy(end_vals[edge], end_h[edge]);
    gkyl_array_release(no_wall);
    gkyl_array_release(fh);
    gkyl_array_release(kh);
    gkyl_array_release(flux);
    gkyl_array_release(area);
  }

  struct gkyl_array *ref = mkarr(gpu, 2 * nb, ext.volume), *rh = mkarr(false, 2 * nb, ext.volume);
  struct gkyl_array *phi_a = mkarr(gpu, nb, ext.volume), *phi_z = mkarr(gpu, nb, ext.volume);
  struct gkyl_array *ah = mkarr(false, nb, ext.volume), *zh = mkarr(false, nb, ext.volume);
  struct gkyl_array *expected = mkarr(false, nb, ext.volume);

  // Ambipolar normalization uses either end separately; zero-current uses both ends.
  for (int edge = 0; edge < (ambipolar ? 2 : 1); ++edge) {
    if (ambipolar) {
      gkyl_boltz_pot_ambi_reference(ambi, &ghost[0], end_vals[edge], ref);
    } else {
      gkyl_boltz_pot_0curr_reference(zero, &ghost[0], end_vals[0], end_vals[1], ref);
    }
    gkyl_array_copy(rh, ref);

    ctx.nref_fac = ambipolar ? (edge ? 0.8 : 1.2) : 1.0;

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost[0]);
    while (gkyl_range_iter_next(&iter)) {
      double xc[3];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      const double *v = gkyl_array_cfetch(rh, gkyl_range_idx(&ghost[0], iter.idx));

      for (int q = 0; q < nq; ++q) {
        double x[3] = {0}, z[3] = {0};
        for (int d = 0; d < cdim - 1; ++d) {
          z[d] = ordinates[q * (cdim - 1) + d];
          x[d] = xc[d] + 0.5 * grid.dx[d] * z[d];
        }
        // The reconstructed reference must be constant along the parallel direction.
        for (int par = 0; par < 3; ++par) {
          z[cdim - 1] = par - 1;
          TEST_CHECK(
            gkyl_compare(basis.eval_expand(z, v), ctx.nref_fac * boltz_test_density(x, &ctx), 1e-12)
          );
          TEST_CHECK(
            gkyl_compare(basis.eval_expand(z, v + nb), boltz_test_potential(x, &ctx), 1e-12)
          );
        }
      }
    }

    // Test the shared volume response both below and above the density clamp.
    boltz_test_project(&grid, &basis, &local, boltz_test_volume_density, &ctx, m0h);
    boltz_test_project(&grid, &basis, &local, boltz_test_volume_phi, &ctx, expected);
    gkyl_array_copy(m0, m0h);

    gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &ext, m0, ref, phi_a);
    gkyl_boltz_pot_0curr_phi_calc(zero, &local, &ext, m0, ref, phi_z);
    gkyl_array_copy(ah, phi_a);
    gkyl_array_copy(zh, phi_z);

    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local, iter.idx);
      const double *a = gkyl_array_cfetch(ah, loc), *b = gkyl_array_cfetch(zh, loc),
                   *e = gkyl_array_cfetch(expected, loc);
      for (int k = 0; k < nb; ++k) {
        TEST_CHECK(gkyl_compare(a[k], b[k], 1e-12));
        // Cancellation of a large gauge offset affects small modal coefficients.
        TEST_CHECK(fabs(a[k] - e[k]) <= 64 * DBL_EPSILON * (1.0 + fabs(e[0])));
        TEST_MSG(
          "cdim=%d p=%d shift=%g mode=%d: %.16g != %.16g", cdim, poly_order, shift, k, a[k], e[k]
        );
      }
    }
  }

  for (int edge = 0; edge < 2; ++edge) {
    gkyl_array_release(end_vals[edge]);
    gkyl_array_release(end_h[edge]);
  }
  gkyl_array_release(m0h);
  gkyl_array_release(wallh);
  gkyl_array_release(m0);
  gkyl_array_release(wall);
  gkyl_array_release(ref);
  gkyl_array_release(rh);
  gkyl_array_release(phi_a);
  gkyl_array_release(phi_z);
  gkyl_array_release(ah);
  gkyl_array_release(zh);
  gkyl_array_release(expected);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_boltz_pot_0curr_release(zero);
}

// Pack a constant physical flux and unit projected face area for the profile tests.
// The chosen flux equals the old fixtures' reconstructed surface flux, preserving
// their analytic reference potentials without requiring lifted ghost-cell data.
static void
calc_reference(
  gkyl_boltz_pot_ambi *up, const struct gkyl_basis *basis, double physical_flux,
  enum gkyl_edge_loc edge, const struct gkyl_range *skin, const struct gkyl_range *ghost,
  const struct gkyl_array *m0_elc, struct gkyl_array *nphi_ref
)
{
  // Store one value per transverse quadrature node in each skin cell.
  int nq = 1;
  for (int d = 0; d < basis->ndim - 1; ++d) {
    nq *= basis->poly_order + 1;
  }
  struct gkyl_array *flux = gkyl_array_new(GKYL_DOUBLE, nq, skin->volume);
  struct gkyl_array *area = gkyl_array_new(GKYL_DOUBLE, nq, skin->volume);
  struct gkyl_array *end_vals = gkyl_array_new(GKYL_DOUBLE, 4 * nq, skin->volume);

  gkyl_array_clear(flux, physical_flux);
  gkyl_array_clear(area, 1.0);

  // Trace the density to the boundary with zero wall potential, then form the reference.
  gkyl_boltz_pot_ambi_boundary(up, edge, skin, flux, area, m0_elc, NULL, end_vals);
  gkyl_boltz_pot_ambi_reference(up, ghost, end_vals, nphi_ref);

  gkyl_array_release(flux);
  gkyl_array_release(area);
  gkyl_array_release(end_vals);
}

void
eval_one(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void
eval_hat(double t, const double *xn, double *restrict fout, void *ctx)
{
  fout[0] = 2. - fabs(xn[0]);
}

void
eval_hat_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  fout[0] = 2. - fabs(z);
}

void
eval_ramp_sheath_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  fout[0] = 2. + psi;
}

void
eval_parabola_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double psi = xn[0], z = xn[1];
  fout[0] = 10. - psi * psi - z * z;
}

void
eval_parabola_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double psi = xn[0], z = xn[1], y = xn[2];
  fout[0] = 10. - psi * psi - z * z - y * y;
}

void
test_boltz_pot_ambi_init_1x_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  TEST_CHECK(ambi->common.cdim == 1);
  TEST_CHECK(ambi->common.phi != NULL);
  TEST_CHECK(ambi->common.use_gpu == use_gpu);
  TEST_CHECK(gkyl_compare_double(grid.dx[cdim - 1], 2. / 8., 1e-12));
  TEST_CHECK(gkyl_compare_double(ambi->common.mass_e, mass_e, 1e-12));
  TEST_CHECK(gkyl_compare_double(ambi->common.charge_e, charge_e, 1e-12));
  TEST_CHECK(gkyl_compare_double(ambi->common.temp_e, temp_e, 1e-12));

  gkyl_free(basis);
  gkyl_boltz_pot_ambi_release(ambi);
}

void
test_boltz_pot_ambi_reference_1x_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);

  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[0], &lower_ghost[0], M0,
    sheath_vals[0]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[0], &upper_ghost[0], M0,
    sheath_vals[1]
  );

  // Serendipity 1x basis is [1/sqrt(2), sqrt(3/2)x].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part.
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) F))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double *sheath_lower_c = ((double *)gkyl_array_cfetch(sheath_vals[0], 0));
  double *sheath_upper_c = ((double *)gkyl_array_cfetch(sheath_vals[1], 9));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_lower_c[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_MSG(
    " Got: %.9e | Expected: %.9e\n", sheath_lower_c[2],
    log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2)
  );
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_upper_c[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[3], 0, 1e-12));

  // This operation happens after the sheaths are determined in the app, so we should test this
  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  int idx_par = cdim - 1, off = 2 * idx_par;
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  // Only the lower sheath values are used for computing the field
  double *sheath_lower_c_avg = ((double *)gkyl_array_cfetch(sheath_vals[0], 0));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_lower_c_avg[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[3], 0, 1e-12));

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_boltz_pot_ambi_phi_calc_1x_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);

  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[0], &lower_ghost[0], M0,
    sheath_vals[0]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[0], &upper_ghost[0], M0,
    sheath_vals[1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  int idx_par = cdim - 1, off = 2 * idx_par;
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[0], phi);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    double *phi_c = ((double *)gkyl_array_cfetch(phi, iter.idx[0]));
    // phi should be the same value as the sheath potential
    TEST_CHECK(
      gkyl_compare_double(phi_c[0], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12)
    );
    TEST_CHECK(gkyl_compare_double(phi_c[1], 0.0, 1e-12));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_boltz_pot_ambi_reference_1x_hat_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_hat =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_hat, NULL);

  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);

  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[0], &lower_ghost[0], M0,
    sheath_vals[0]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[0], &upper_ghost[0], M0,
    sheath_vals[1]
  );

  // Serendipity 1x basis is [1/sqrt(2), sqrt(3/2)x].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part.
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) F))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double *sheath_lower_c = ((double *)gkyl_array_cfetch(sheath_vals[0], 0));
  double *sheath_upper_c = ((double *)gkyl_array_cfetch(sheath_vals[1], 9));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_lower_c[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_upper_c[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_CHECK(gkyl_compare_double(sheath_upper_c[3], 0, 1e-12));

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  int idx_par = cdim - 1, off = 2 * idx_par;
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  // Only the lower sheath values are used for calculations
  double *sheath_lower_c_avg = ((double *)gkyl_array_cfetch(sheath_vals[0], 0));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[0], sqrt(2), 1e-12));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[1], 0, 1e-12));
  TEST_CHECK(gkyl_compare_double(
    sheath_lower_c_avg[2], log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1])) * sqrt(2), 1e-12
  ));
  TEST_CHECK(gkyl_compare_double(sheath_lower_c_avg[3], 0, 1e-12));

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

void
test_boltz_pot_ambi_phi_calc_1x_hat_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_hat =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_hat, NULL);

  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);

  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[0], &lower_ghost[0], M0,
    sheath_vals[0]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[0], &upper_ghost[0], M0,
    sheath_vals[1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  int idx_par = cdim - 1, off = 2 * idx_par;
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[0], phi);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  while (gkyl_range_iter_next(&iter)) {
    double *phi_c = ((double *)gkyl_array_cfetch(phi, iter.idx[0]));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, iter.idx[0]));
    double ni = ni_c[0] / sqrt(2);
    TEST_CHECK(gkyl_compare_double(phi_c[0] / sqrt(2), phi_sheath + log(ni), 2e-4));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

void
test_boltz_pot_ambi_init_2x_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {8, 16};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  TEST_CHECK(ambi->common.cdim == 2);
  TEST_CHECK(ambi->common.phi != NULL);
  TEST_CHECK(ambi->common.use_gpu == use_gpu);
  TEST_CHECK(gkyl_compare_double(grid.dx[cdim - 1], 2. / 16., 1e-12)
  ); // Second direction is field line length
  TEST_CHECK(gkyl_compare_double(ambi->common.mass_e, mass_e, 1e-12));
  TEST_CHECK(gkyl_compare_double(ambi->common.charge_e, charge_e, 1e-12));
  TEST_CHECK(gkyl_compare_double(ambi->common.temp_e, temp_e, 1e-12));

  gkyl_free(basis);
  gkyl_boltz_pot_ambi_release(ambi);
}

void
test_boltz_pot_ambi_reference_2x_one_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {8, 16};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);

  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  // Serendipity 2x basis is [1/2,(sqrt(3)*x)/2,(sqrt(3)*y)/2,(3*x*y)/2].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part.
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) F))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  struct gkyl_range_iter iter;
  for (int ix_cdim = 0; ix_cdim < cdim; ix_cdim++) {
    // Only the lower cells are used to calculate the field
    gkyl_range_iter_init(&iter, &lower_ghost[ix_cdim]);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&lower_ghost[ix_cdim], iter.idx);
      double *sheath_lower_c = ((double *)gkyl_array_cfetch(sheath_vals[off], lidx));
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], 2, 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], 0, 1e-12));
      }
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[2], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4], phi_sheath * 2, 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4], 0, 1e-12));
      }
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[5], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[6], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[7], 0, 1e-12));
    }
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_boltz_pot_ambi_reference_2x_hat_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {8, 8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_hat =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_hat_2x, NULL);

  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  // Serendipity 2x basis is [1/2,(sqrt(3)*x)/2,(sqrt(3)*y)/2,(3*x*y)/2].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part.
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) F))
  // phi_s = log(1/(sqrt(2*pi) 0.5 * 1/4))
  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  struct gkyl_range_iter iter;
  for (int ix_cdim = 0; ix_cdim < cdim; ix_cdim++) {
    // Only the lower cells are used to calculate the field
    gkyl_range_iter_init(&iter, &lower_ghost[ix_cdim]);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&lower_ghost[ix_cdim], iter.idx);
      double *sheath_lower_c = ((double *)gkyl_array_cfetch(sheath_vals[off], lidx));
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], 2, 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], 0, 1e-12));
      }
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[2], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4], phi_sheath * 2, 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4], 0, 1e-12));
      }
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[5], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[6], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[7], 0, 1e-12));
    }
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

void
test_boltz_pot_ambi_reference_2x_ramp_sheath_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {8, 8};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_ramp =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_ramp_sheath_2x, NULL);

  gkyl_proj_on_basis_advance(proj_ramp, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  // Serendipity 2x basis is [1/2,(sqrt(3)*x)/2,(sqrt(3)*y)/2,(3*x*y)/2].
  // sheath_vals stores both the ion density and sheath value
  // Ion density is the first 2 components. Should be 1
  // Sheath potential is the second part.
  // phi_s = Te/e * log(ni (Te/me) / (sqrt(2*pi) F))
  struct gkyl_range_iter iter;
  for (int ix_cdim = 0; ix_cdim < cdim; ix_cdim++) {
    // Only the lower cells are used to calculate the field
    gkyl_range_iter_init(&iter, &lower_ghost[ix_cdim]);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&lower_ghost[ix_cdim], iter.idx);
      double *sheath_lower_c = ((double *)gkyl_array_cfetch(sheath_vals[off], lidx));
      double *density_c = ((double *)gkyl_array_cfetch(M0, lidx));
      double phi_sheath = 0.5 * log((density_c[0] - density_c[1]) / 2.0) +
                          0.5 * log((density_c[0] + density_c[1]) / 2.0) -
                          log(sqrt(2 * M_PI) * grid.dx[cdim - 1]);
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], density_c[0], 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[0], 0, 1e-12));
      }
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], density_c[1], 1e-12));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[1], 0, 1e-12));
      }
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[2], 0, 1e-12));
      TEST_CHECK(gkyl_compare_double(sheath_lower_c[3], 0, 1e-12));
      // Not exact because division happens at quadrature nodes
      if (ix_cdim == 1) {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4] / 2, phi_sheath, 1e-3));
      } else {
        TEST_CHECK(gkyl_compare_double(sheath_lower_c[4], 0, 1e-12));
      }
      // Slope of sheath potential is difficult to calculate
    }
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_ramp);
}

void
test_boltz_pot_ambi_phi_calc_2x_one_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {8, 16};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);

  gkyl_proj_on_basis_advance(proj_one, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local, iter.idx);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    TEST_CHECK(gkyl_compare_double(phi_c[0] / 2, phi_sheath, 1e-12));
    TEST_CHECK(gkyl_compare_double(phi_c[1], 0.0, 1e-12));
    TEST_CHECK(gkyl_compare_double(phi_c[2], 0.0, 1e-12));
    TEST_CHECK(gkyl_compare_double(phi_c[3], 0.0, 1e-12));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
}

void
test_boltz_pot_ambi_phi_calc_2x_hat_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {32, 32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_hat =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_hat_2x, NULL);

  gkyl_proj_on_basis_advance(proj_hat, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local, iter.idx);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, lidx));
    double ni = ni_c[0] / 2;
    TEST_CHECK(gkyl_compare_double(phi_c[0] / 2, phi_sheath + log(ni), 1e-3));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_hat);
}

void
test_boltz_pot_ambi_phi_calc_2x_ramp_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {32, 32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_ramp =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_ramp_sheath_2x, NULL);

  gkyl_proj_on_basis_advance(proj_ramp, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  double phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&local, iter.idx);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, lidx));
    double ni = ni_c[0] / 2;
    TEST_CHECK(gkyl_compare_double(phi_c[0] / 2, phi_sheath + log(ni), 1e-3));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_ramp);
}

void
test_boltz_pot_ambi_phi_calc_2x_parabola_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {32, 32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_func =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_parabola_2x, NULL);

  gkyl_proj_on_basis_advance(proj_func, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  struct gkyl_range_iter iter;
  int idx_ghost[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx_ghost);
    idx_ghost[cdim - 1] = local_ext.lower[cdim - 1];
    long lidx = gkyl_range_idx(&local, iter.idx);
    long ghost_lidx = gkyl_range_idx(&local_ext, idx_ghost);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, lidx));
    double *sheath_vals_c = ((double *)gkyl_array_cfetch(sheath_vals[off], ghost_lidx));
    double ni = ni_c[0] / 2;
    double ni_sheath = sheath_vals_c[0] / 2;
    double phi_sheath = sheath_vals_c[4] / 2;
    // The modal mean of log(n/n_ref) is the quadrature mean of the log,
    // rather than the log of the modal means.
    double expected = 0.0;
    for (int sx = -1; sx <= 1; sx += 2) {
      for (int sy = -1; sy <= 1; sy += 2) {
        double density = (ni_c[0] + sx * ni_c[1] + sy * ni_c[2] + sx * sy * ni_c[3]) / 2.0;
        double reference = (sheath_vals_c[0] + sx * sheath_vals_c[1]) / 2.0;
        double potential = (sheath_vals_c[4] + sx * sheath_vals_c[5]) / 2.0;
        expected += 0.25 * (potential + log(fmax(density, reference) / reference));
      }
    }
    TEST_CHECK(gkyl_compare_double(phi_c[0] / 2.0, expected, 1e-12));
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_func);
}

void
test_boltz_pot_ambi_phi_calc_3x_one_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {32, 32, 32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_func =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);

  gkyl_proj_on_basis_advance(proj_func, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  struct gkyl_range_iter iter;
  int idx_ghost[GKYL_MAX_CDIM];
  double known_phi_sheath = log(1 / (sqrt(2 * M_PI) * grid.dx[cdim - 1]));
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx_ghost);
    idx_ghost[cdim - 1] = local_ext.lower[cdim - 1];
    long lidx = gkyl_range_idx(&local, iter.idx);
    long ghost_lidx = gkyl_range_idx(&local_ext, idx_ghost);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, lidx));
    double *sheath_vals_c = ((double *)gkyl_array_cfetch(sheath_vals[off], ghost_lidx));
    double ni = ni_c[0] / pow(2., 3. / 2.);
    double ni_sheath = sheath_vals_c[0] / pow(2., 3. / 2.);
    double phi_sheath = sheath_vals_c[8] / pow(2., 3. / 2.); // 8 coefficients in serendipity 3xP1
    TEST_CHECK(gkyl_compare_double(ni_sheath, 1.0, 1e-12));
    TEST_CHECK(gkyl_compare_double(phi_sheath, known_phi_sheath, 1e-12));
    TEST_CHECK(
      gkyl_compare_double(phi_c[0] / pow(2., 3. / 2.), phi_sheath + log(ni / ni_sheath), 1e-12)
    );
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_func);
}

void
test_boltz_pot_ambi_phi_calc_3x_parabola_ho()
{
  int poly_order = 1;
  double lower[] = {-1.0, -1.0, -1.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {32, 32, 32};
  int cdim = sizeof(lower) / sizeof(lower[0]);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = gkyl_cart_modal_serendip_new(cdim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  double mass_e = 1.0, charge_e = -1.0, temp_e = 1.0;
  bool use_gpu = false;

  struct gkyl_boltz_pot_ambi *ambi =
    gkyl_boltz_pot_ambi_new(basis, mass_e, charge_e, temp_e, use_gpu);

  struct gkyl_array *sheath_vals[2 * cdim];
  for (int j = 0; j < cdim; ++j) {
    sheath_vals[2 * j] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    sheath_vals[2 * j + 1] = gkyl_array_new(GKYL_DOUBLE, 2 * basis->num_basis, local_ext.volume);
    gkyl_array_clear(sheath_vals[2 * j], 0.0);
    gkyl_array_clear(sheath_vals[2 * j + 1], 0.0);
  }

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &local_ext, ghost
    );
  }

  struct gkyl_array *M0 = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_proj_on_basis *proj_one =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_one, NULL);
  gkyl_proj_on_basis *proj_func =
    gkyl_proj_on_basis_new(&grid, basis, poly_order + 1, 1, eval_parabola_3x, NULL);

  gkyl_proj_on_basis_advance(proj_func, 0.0, &local_ext, M0);

  int idx_par = cdim - 1, off = 2 * idx_par;
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_LOWER_EDGE, &lower_skin[idx_par], &lower_ghost[idx_par],
    M0, sheath_vals[off]
  );
  calc_reference(
    ambi, basis, grid.dx[cdim - 1], GKYL_UPPER_EDGE, &upper_skin[idx_par], &upper_ghost[idx_par],
    M0, sheath_vals[off + 1]
  );

  // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
  gkyl_array_copy_range_to_range(
    sheath_vals[off + 1], sheath_vals[off + 1], &lower_ghost[idx_par], &upper_ghost[idx_par]
  );
  gkyl_array_accumulate(sheath_vals[off], 1., sheath_vals[off + 1]);
  gkyl_array_scale(sheath_vals[off], 0.5);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis->num_basis, local_ext.volume);

  gkyl_boltz_pot_ambi_phi_calc(ambi, &local, &local_ext, M0, sheath_vals[off], phi);

  struct gkyl_range_iter iter;
  int idx_ghost[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx_ghost);
    idx_ghost[cdim - 1] = local_ext.lower[cdim - 1];
    long lidx = gkyl_range_idx(&local, iter.idx);
    long ghost_lidx = gkyl_range_idx(&local_ext, idx_ghost);
    double *phi_c = ((double *)gkyl_array_cfetch(phi, lidx));
    double *ni_c = ((double *)gkyl_array_cfetch(M0, lidx));
    double *sheath_vals_c = ((double *)gkyl_array_cfetch(sheath_vals[off], ghost_lidx));
    double ni = ni_c[0] / pow(2., 3. / 2.); // 1st coefficient is 1/2^(3/2)
    double ni_sheath = sheath_vals_c[0] / pow(2., 3. / 2.);
    double phi_sheath = sheath_vals_c[8] / pow(2., 3. / 2.); // 8 coefficients in serendipity 3xP1
    TEST_CHECK(
      gkyl_compare_double(phi_c[0] / pow(2., 3. / 2.), phi_sheath + log(ni / ni_sheath), 1e-5)
    );
  }

  gkyl_free(basis);
  for (int j = 0; j < cdim; ++j) {
    gkyl_array_release(sheath_vals[2 * j]);
    gkyl_array_release(sheath_vals[2 * j + 1]);
  }
  gkyl_array_release(M0);
  gkyl_array_release(phi);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_proj_on_basis_release(proj_one);
  gkyl_proj_on_basis_release(proj_func);
}

// Recover a known surface flux from the actual DG boundary lift, at both ends
// and two cell widths. The former trace-based formula returns twice this flux.
static void
test_boltz_boundary_lift(void)
{
  double w[2] = {0.0}, vmap[2] = {0.0, sqrt(2.0 / 3.0)}, prime[1] = {1.0};
  double one[2] = {sqrt(2.0), 0.0}, zero[2] = {0.0};
  for (int resolution = 0; resolution < 2; ++resolution) {
    double dx[2] = {resolution ? 0.125 : 0.5, 2.0};
    for (int edge = -1; edge <= 1; edge += 2) {
      double surface[3] = {2.0 * edge, 2.0 * edge, 2.0 * edge}, lift[6] = {0.0}, moment[2] = {0.0},
             out[4] = {0.0};
      dg_gyrokinetic_boundary_surfx_1x1v_ser_p1(w, dx, prime, prime, surface, surface, edge, lift);
      mom_gyrokinetic_M0_1x1v_ser_p1(dx, vmap, 1.0, 1.0, one, zero, lift, moment);
      double recovered = dx[0] * moment[0] / sqrt(2.0);
      TEST_CHECK(gkyl_compare(recovered, 4.0, 1e-12));
      double trace = dx[0] / 2.0 * (moment[0] / sqrt(2.0) - edge * sqrt(1.5) * moment[1]);
      TEST_CHECK(gkyl_compare(trace, 2.0 * recovered, 1e-12));
      double end_vals[4] = {recovered, 1.0, 1.0, 0.0};
      boltz_pot_ambi_reference_1x_ser_p1(-1, 1, 1, end_vals, end_vals, out);
      TEST_CHECK(gkyl_compare(out[2] / sqrt(2.0), log(1.0 / (sqrt(2.0 * M_PI) * 4.0)), 1e-12));
    }
  }
}

// Check both orders and all dimensions against manufactured boundary and volume fields.
static void
test_boltz_manufactured_ambi(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      for (int shift = 0; shift < 2; ++shift) {
        check_boltz_model(d, p, false, true, shift ? 2000.0 : 0.0);
      }
    }
  }
}

// Check asymmetric end areas/fluxes, unequal biases, gauge shifts, and the
// invalid-loss fallback against the closed-form combined-current relation.
static void
check_combined_reference(int cdim, bool gpu)
{
  double lo[3] = {-1, -1, -1}, hi[3] = {1, 1, 1};
  int cells[3] = {2, 2, 2}, ng[3] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lo, hi, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, 1);
  struct gkyl_range local, ext, skin, ghost;
  gkyl_create_grid_ranges(&grid, ng, &ext, &local);
  gkyl_skin_ghost_ranges(&skin, &ghost, cdim - 1, GKYL_LOWER_EDGE, &ext, ng);
  struct gkyl_boltz_pot_0curr *up = gkyl_boltz_pot_0curr_new(&basis, 1, -1, 2, gpu);
  int nq = 1 << (cdim - 1), nb = basis.num_basis;
  struct gkyl_array *lh = gkyl_array_new(GKYL_DOUBLE, 4 * nq, skin.volume);
  struct gkyl_array *uh = gkyl_array_new(GKYL_DOUBLE, 4 * nq, skin.volume);
  struct gkyl_array *lower = gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, 4 * nq, skin.volume) :
                                   gkyl_array_acquire(lh);
  struct gkyl_array *upper = gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, 4 * nq, skin.volume) :
                                   gkyl_array_acquire(uh);
  struct gkyl_array *ref = gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, 2 * nb, ext.volume) :
                                 gkyl_array_new(GKYL_DOUBLE, 2 * nb, ext.volume);
  struct gkyl_array *rh = gkyl_array_new(GKYL_DOUBLE, 2 * nb, ext.volume);
  // Build nonfinite inputs from IEEE bits so -ffast-math does not discard them.
  uint64_t invalid_bits[] = {UINT64_C(0x7ff0000000000000), UINT64_C(0x7ff8000000000000)};
  double invalid[2];
  memcpy(invalid, invalid_bits, sizeof invalid);
  for (int mode = 0; mode < 7; ++mode) {
    double shift = mode == 1 ? 3.5 : mode == 2 ? 2000.0 : 0.0;
    for (long cell = 0; cell < skin.volume; ++cell) {
      double *l = gkyl_array_fetch(lh, cell), *u = gkyl_array_fetch(uh, cell);
      for (int q = 0; q < nq; ++q) {
        l[q] = mode == 3 ? 0.0 : 2.0 + q;
        u[q] = mode == 3 ? 0.0 : 1.0 + 0.5 * q;
        if (mode == 4) {
          l[q] = invalid[0];
        }
        if (mode == 5) {
          l[q] = invalid[1];
        }
        l[nq + q] = 1.0 + 0.1 * q;
        u[nq + q] = 0.5 + 0.2 * q;
        l[2 * nq + q] = 3.0 + q;
        u[2 * nq + q] = 5.0 + q;
        if (mode == 6) {
          l[2 * nq + q] = invalid[0];
        }
        l[3 * nq + q] = 0.2 * q + shift;
        u[3 * nq + q] = -0.5 + shift;
      }
    }
    gkyl_array_copy(lower, lh);
    gkyl_array_copy(upper, uh);
    gkyl_boltz_pot_0curr_reference(up, &ghost, lower, upper, ref);
    gkyl_array_copy(rh, ref);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost);
    while (gkyl_range_iter_next(&iter)) {
      const double *val = gkyl_array_cfetch(rh, gkyl_range_idx(&ghost, iter.idx));
      for (int q = 0; q < nq; ++q) {
        double density = 0.0, phi = 0.0, ordinates[18], weights[9], coord[3] = {0}, face_basis[20];
        if (cdim > 1) {
          gkyl_ndim_ordinates_weights(cdim - 1, ordinates, weights, 2);
        }
        for (int d = 0; d < cdim - 1; ++d) {
          coord[d] = ordinates[q * (cdim - 1) + d];
        }
        coord[cdim - 1] = -1.0;
        basis.eval(coord, face_basis);
        for (int k = 0; k < nb; ++k) {
          density += face_basis[k] * val[k];
          phi += face_basis[k] * val[nb + k];
        }
        double expected =
          mode >= 3 ?
            0.0 :
            shift + 2.0 * log(
                            (4.0 + q) * sqrt(1.0 / M_PI) *
                            ((1.0 + 0.1 * q) * exp(0.1 * q) + (0.5 + 0.2 * q) * exp(-0.25)) /
                            (3.0 + 1.5 * q)
                          );
        TEST_CHECK(gkyl_compare(density, mode == 6 ? 0.0 : 4.0 + q, 1e-12));
        TEST_CHECK(gkyl_compare(phi, expected, 1e-12));
        TEST_MSG("cdim=%d mode=%d q=%d: phi %.16g != %.16g", cdim, mode, q, phi, expected);
      }
    }
  }
  gkyl_array_release(lh);
  gkyl_array_release(uh);
  gkyl_array_release(lower);
  gkyl_array_release(upper);
  gkyl_array_release(ref);
  gkyl_array_release(rh);
  gkyl_boltz_pot_0curr_release(up);
}

// Exercise the combined closure in each configuration dimension on the CPU.
static void
test_boltz_combined_reference(void)
{
  for (int d = 1; d <= 3; ++d) {
    check_combined_reference(d, false);
  }
}

// Check both orders and all dimensions against manufactured boundary and volume fields.
static void
test_boltz_manufactured_0curr(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      for (int shift = 0; shift < 2; ++shift) {
        check_boltz_model(d, p, false, false, shift ? 2000.0 : 0.0);
      }
    }
  }
}

// Equal-end data must give the same normalization for both models. Invalid
// loss, density, geometry, or bias must yield the documented finite fallback.
static void
check_boltz_fallback(int cdim, int poly_order, bool gpu)
{
  double lo[3] = {-1, -1, -1}, hi[3] = {1, 1, 1};
  int cells[3] = {1, 1, 1}, ng[3] = {1, 1, 1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lo, hi, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  struct gkyl_range local, ext, skin, ghost;
  gkyl_create_grid_ranges(&grid, ng, &ext, &local);
  gkyl_skin_ghost_ranges(&skin, &ghost, cdim - 1, GKYL_LOWER_EDGE, &ext, ng);

  int nb = basis.num_basis, nq = 1;
  for (int d = 0; d < cdim - 1; ++d) {
    nq *= poly_order + 1;
  }

  gkyl_boltz_pot_ambi *ambi = gkyl_boltz_pot_ambi_new(&basis, 1, -1, 2, gpu);
  gkyl_boltz_pot_0curr *zero = gkyl_boltz_pot_0curr_new(&basis, 1, -1, 2, gpu);

  struct gkyl_array *ends = mkarr(gpu, 4 * nq, skin.volume),
                    *eh = mkarr(false, 4 * nq, skin.volume);
  struct gkyl_array *ar = mkarr(gpu, 2 * nb, ext.volume), *zr = mkarr(gpu, 2 * nb, ext.volume);
  struct gkyl_array *ah = mkarr(false, 2 * nb, ext.volume), *zh = mkarr(false, 2 * nb, ext.volume);

  // Construct quiet NaN and positive infinity directly, including in fast-math builds.
  uint64_t invalid_bits[] = {UINT64_C(0x7ff8000000000000), UINT64_C(0x7ff0000000000000)};
  double invalid[2];
  memcpy(invalid, invalid_bits, sizeof invalid);

  // Mode 0 is valid; the remaining cases invalidate flux, density, area, or wall bias.
  for (int mode = 0; mode < 14; ++mode) {
    double *v = gkyl_array_fetch(eh, 0);
    for (int q = 0; q < nq; ++q) {
      v[q] = 2.0;
      v[nq + q] = 1.0;
      v[2 * nq + q] = 3.0;
      v[3 * nq + q] = -0.5;
      if (mode == 1) {
        v[q] = 0.0;
      }
      if (mode == 2) {
        v[q] = -1.0;
      }
      if (mode == 3 || mode == 4) {
        v[q] = invalid[mode - 3];
      }
      if (mode == 5) {
        v[2 * nq + q] = invalid[0];
      }
      if (mode == 6) {
        v[2 * nq + q] = 0.0;
      }
      if (mode == 7) {
        v[2 * nq + q] = -1.0;
      }
      if (mode == 8) {
        v[2 * nq + q] = invalid[1];
      }
      if (mode == 9) {
        v[nq + q] = 0.0;
      }
      if (mode == 10 || mode == 11) {
        v[nq + q] = invalid[mode - 10];
      }
      if (mode == 12 || mode == 13) {
        v[3 * nq + q] = invalid[mode - 12];
      }
    }
    gkyl_array_copy(ends, eh);

    // Identical ends make the two models equivalent, including their fallbacks.
    gkyl_boltz_pot_ambi_reference(ambi, &ghost, ends, ar);
    gkyl_boltz_pot_0curr_reference(zero, &ghost, ends, ends, zr);

    gkyl_array_copy(ah, ar);
    gkyl_array_copy(zh, zr);
    const double *a = gkyl_array_cfetch(ah, gkyl_range_idx(&ghost, ghost.lower));
    const double *z = gkyl_array_cfetch(zh, gkyl_range_idx(&ghost, ghost.lower));
    // Invalid data gives zero reference potential; constant profiles have only mode 0.
    double expected = mode ? 0.0 : -0.5 + 2.0 * log(3.0 * sqrt(1.0 / M_PI) / 2.0);
    double norm = pow(sqrt(2.0), cdim);
    double dens = (mode == 5 || mode == 6 || mode == 8) ? 0.0 : mode == 7 ? -1.0 : 3.0;
    for (int k = 0; k < nb; ++k) {
      TEST_CHECK(gkyl_is_dfinite(a[k]) && gkyl_is_dfinite(a[nb + k]));
      TEST_CHECK(gkyl_is_dfinite(z[k]) && gkyl_is_dfinite(z[nb + k]));
      TEST_CHECK(gkyl_compare(a[k], k ? 0.0 : norm * dens, 1e-12));
      TEST_CHECK(gkyl_compare(a[nb + k], k ? 0.0 : norm * expected, 1e-12));
      TEST_CHECK(gkyl_compare(a[k], z[k], 1e-12));
      TEST_CHECK(gkyl_compare(a[nb + k], z[nb + k], 1e-12));
      TEST_MSG("fallback cdim=%d p=%d mode=%d k=%d", cdim, poly_order, mode, k);
    }
  }
  gkyl_array_release(ends);
  gkyl_array_release(eh);
  gkyl_array_release(ar);
  gkyl_array_release(zr);
  gkyl_array_release(ah);
  gkyl_array_release(zh);
  gkyl_boltz_pot_ambi_release(ambi);
  gkyl_boltz_pot_0curr_release(zero);
}

// Verify equal-end model agreement and invalid-input fallbacks in every supported basis.
static void
test_boltz_equal_ends_and_fallback(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      check_boltz_fallback(d, p, false);
    }
  }
}

#ifdef GKYL_HAVE_CUDA
// Check both orders and all dimensions against manufactured boundary and volume fields.
static void
test_boltz_manufactured_ambi_dev(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      for (int shift = 0; shift < 2; ++shift) {
        check_boltz_model(d, p, true, true, shift ? 2000.0 : 0.0);
      }
    }
  }
}

// Exercise the same combined closure on device arrays.
static void
test_boltz_combined_reference_dev(void)
{
  for (int d = 1; d <= 3; ++d) {
    check_combined_reference(d, true);
  }
}

// Check both orders and all dimensions against manufactured boundary and volume fields.
static void
test_boltz_manufactured_0curr_dev(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      for (int shift = 0; shift < 2; ++shift) {
        check_boltz_model(d, p, true, false, shift ? 2000.0 : 0.0);
      }
    }
  }
}

// Verify equal-end model agreement and invalid-input fallbacks in every supported basis.
static void
test_boltz_equal_ends_and_fallback_dev(void)
{
  for (int d = 1; d <= 3; ++d) {
    for (int p = 1; p <= 2; ++p) {
      check_boltz_fallback(d, p, true);
    }
  }
}
#endif

TEST_LIST = {
  {"test_boltz_manufactured_ambi", test_boltz_manufactured_ambi},
  {"test_boltz_boundary_lift", test_boltz_boundary_lift},
  {"test_boltz_pot_ambi_init_1x_ho", test_boltz_pot_ambi_init_1x_ho},
  {"test_boltz_pot_ambi_reference_1x_ho", test_boltz_pot_ambi_reference_1x_ho},
  {"test_boltz_pot_ambi_phi_calc_1x_ho", test_boltz_pot_ambi_phi_calc_1x_ho},
  {"test_boltz_pot_ambi_reference_1x_hat_ho", test_boltz_pot_ambi_reference_1x_hat_ho},
  {"test_boltz_pot_ambi_phi_calc_1x_hat_ho", test_boltz_pot_ambi_phi_calc_1x_hat_ho},
  {"test_boltz_pot_ambi_init_2x_ho", test_boltz_pot_ambi_init_2x_ho},
  {"test_boltz_pot_ambi_reference_2x_one_ho", test_boltz_pot_ambi_reference_2x_one_ho},
  {"test_boltz_pot_ambi_reference_2x_hat_ho", test_boltz_pot_ambi_reference_2x_hat_ho},
  {"test_boltz_pot_ambi_reference_2x_ramp_sheath_ho",
   test_boltz_pot_ambi_reference_2x_ramp_sheath_ho},
  {"test_boltz_pot_ambi_phi_calc_2x_one_ho", test_boltz_pot_ambi_phi_calc_2x_one_ho},
  {"test_boltz_pot_ambi_phi_calc_2x_hat_ho", test_boltz_pot_ambi_phi_calc_2x_hat_ho},
  {"test_boltz_pot_ambi_phi_calc_2x_ramp_ho", test_boltz_pot_ambi_phi_calc_2x_ramp_ho},
  {"test_boltz_pot_ambi_phi_calc_2x_parabola_ho", test_boltz_pot_ambi_phi_calc_2x_parabola_ho},
  {"test_boltz_pot_ambi_phi_calc_3x_one_ho", test_boltz_pot_ambi_phi_calc_3x_one_ho},
  {"test_boltz_pot_ambi_phi_calc_3x_parabola_ho", test_boltz_pot_ambi_phi_calc_3x_parabola_ho},
  {"test_boltz_equal_ends_and_fallback", test_boltz_equal_ends_and_fallback},
  {"test_boltz_manufactured_0curr", test_boltz_manufactured_0curr},
  {"test_boltz_combined_reference", test_boltz_combined_reference},
#ifdef GKYL_HAVE_CUDA
  {"test_boltz_manufactured_ambi_dev", test_boltz_manufactured_ambi_dev},
  {"test_boltz_equal_ends_and_fallback_dev", test_boltz_equal_ends_and_fallback_dev},
  {"test_boltz_manufactured_0curr_dev", test_boltz_manufactured_0curr_dev},
  {"test_boltz_combined_reference_dev", test_boltz_combined_reference_dev},
#endif
  {NULL, NULL}
};
