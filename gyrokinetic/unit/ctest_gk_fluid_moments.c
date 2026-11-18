#include <acutest.h>
#include <gkyl_gk_fluid_moments.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>
#include <string.h>

static void
test_moments(bool use_gpu)
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    struct gkyl_basis basis;
    gkyl_cart_modal_serendip(&basis, cdim, 1);
    int nb = basis.num_basis;
    double lower[] = {-1, -1, -1}, upper[] = {1, 1, 1};
    int cells[] = {1, 1, 1}, ghost[] = {1, 1, 1};
    struct gkyl_rect_grid grid;
    gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
    struct gkyl_range ext, local;
    gkyl_create_grid_ranges(&grid, ghost, &ext, &local);
    struct gkyl_array *state = gkyl_array_new(GKYL_DOUBLE, 4 * nb, ext.volume);
    struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, nb, ext.volume);
    double norm = pow(sqrt(2.0), cdim);
    // m=2, Jn=2, u=3, Tpar/m=5, Tperp/m=7, B=3.
    double perp = 28;
    gkyl_array_shiftc(state, 4 * norm, 0);
    gkyl_array_shiftc(state, 12 * norm, nb);
    gkyl_array_shiftc(state, (28 + perp) * norm, 2 * nb);
    gkyl_array_shiftc(state, perp / 3 * norm, 3 * nb);
    gkyl_array_shiftc(phi, 11 * norm, 0);
    struct gkyl_array *state_dev = use_gpu ?
                                     gkyl_array_cu_dev_new(GKYL_DOUBLE, state->ncomp, state->size) :
                                     gkyl_array_acquire(state);
    struct gkyl_array *phi_dev =
      use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
    gkyl_array_copy(state_dev, state);
    gkyl_array_copy(phi_dev, phi);
    struct gkyl_array *bmag = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nb, ext.volume) :
                                        gkyl_array_new(GKYL_DOUBLE, nb, ext.volume);
    gkyl_array_clear(bmag, 0);
    gkyl_array_shiftc(bmag, 3 * norm, 0);
    enum gkyl_distribution_moments ids[] = {
      GKYL_F_MOMENT_M0,
      GKYL_F_MOMENT_M1,
      GKYL_F_MOMENT_M2,
      GKYL_F_MOMENT_M2PAR,
      GKYL_F_MOMENT_M2PERP,
      GKYL_F_MOMENT_M3PAR,
      GKYL_F_MOMENT_M3PERP,
      GKYL_F_MOMENT_MAXWELLIAN,
      GKYL_F_MOMENT_BIMAXWELLIAN,
      GKYL_F_MOMENT_M0M1M2,
      GKYL_F_MOMENT_M0M1M2PARM2PERP,
      GKYL_F_MOMENT_HAMILTONIAN
    };
    double expected[][4] = {
      {2},
      {6},
      {28 + perp},
      {28},
      {perp},
      {144},
      {3 * perp},
      {2, 3, 19.0 / 3},
      {2, 3, 5, 7},
      {2, 6, 28 + perp},
      {2, 6, 28, perp},
      {2, 12, 28 + perp - 22}
    };
    for (int m = 0; m < sizeof(ids) / sizeof(ids[0]); ++m) {
      for (int integrated = 0; integrated <= 1; ++integrated) {
        TEST_CASE_("moment=%d, cdim=%d, integrated=%d", ids[m], cdim, integrated);
        struct gkyl_gk_fluid_moments *up =
          gkyl_gk_fluid_moments_inew(&(struct gkyl_gk_fluid_moments_inp){
            .basis = &basis,
            .grid = &grid,
            .range = &local,
            .range_ext = &ext,
            .bmag = bmag,
            .mass = 2,
            .charge = -1,
            .num_moments = 4,
            .moment = ids[m],
            .integrated = integrated,
            .use_gpu = use_gpu,
          });
        int nm = gkyl_gk_fluid_moments_num_mom(up), nout = integrated ? 1 : nb;
        if (ids[m] == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
          TEST_CHECK(nm == 4);
        }
        struct gkyl_array *out = use_gpu ?
                                   gkyl_array_cu_dev_new(GKYL_DOUBLE, nm * nout, ext.volume) :
                                   gkyl_array_new(GKYL_DOUBLE, nm * nout, ext.volume);
        struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, nm * nout, ext.volume);
        gkyl_gk_fluid_moments_advance(up, &local, state_dev, phi_dev, out);
        gkyl_array_copy(host, out);
        const double *data = gkyl_array_cfetch(host, gkyl_range_idx(&local, local.lower));
        for (int c = 0; c < nm; ++c) {
          TEST_CHECK(
            gkyl_compare(data[c * nout], expected[m][c] * (integrated ? norm * norm : norm), 1e-12)
          );
          for (int k = 1; k < nout; ++k) {
            TEST_CHECK(fabs(data[c * nout + k]) < 1e-11);
          }
        }
        gkyl_array_release(host);
        gkyl_array_release(out);
        gkyl_gk_fluid_moments_release(up);
      }
    }
    gkyl_array_release(state_dev);
    gkyl_array_release(phi_dev);
    gkyl_array_release(state);
    gkyl_array_release(phi);
    gkyl_array_release(bmag);
  }
  TEST_CASE((const char *)NULL);
}

// Independent p1 arithmetic in the power representation a+b*x on [-1,1].
struct linear {
  double a, b;
};

static struct linear
mul(struct linear lhs, struct linear rhs)
{
  return (struct linear){lhs.a * rhs.a + lhs.b * rhs.b / 3, lhs.a * rhs.b + lhs.b * rhs.a};
}

static struct linear
add(struct linear lhs, double factor, struct linear rhs)
{
  return (struct linear){lhs.a + factor * rhs.a, lhs.b + factor * rhs.b};
}

static struct linear
scale(double factor, struct linear inp)
{
  return (struct linear){factor * inp.a, factor * inp.b};
}

static struct linear
divide(struct linear num, struct linear den)
{
  double det = den.a * den.a - den.b * den.b / 3;
  return (struct linear){
    (den.a * num.a - den.b * num.b / 3) / det,
    (den.a * num.b - den.b * num.a) / det,
  };
}

static void
check_linear(const double *actual, struct linear expected)
{
  TEST_CHECK(gkyl_compare(actual[0], sqrt(2.0) * expected.a, 1e-11));
  TEST_CHECK(gkyl_compare(actual[1], sqrt(2.0 / 3) * expected.b, 1e-11));
}

static void
set_linear(struct gkyl_array *array, long loc, int comp, struct linear value)
{
  double *out = gkyl_array_fetch(array, loc);
  out[2 * comp] = sqrt(2.0) * value.a;
  out[2 * comp + 1] = sqrt(2.0 / 3) * value.b;
}

static void
test_spatial_products(bool use_gpu, bool isotropic)
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);
  struct gkyl_rect_grid grid;
  double lower[] = {-1}, upper[] = {1};
  int cells[] = {1}, ghost[] = {1};
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range ext, local;
  gkyl_create_grid_ranges(&grid, ghost, &ext, &local);
  struct gkyl_array *(*alloc)(enum gkyl_elem_type, size_t, size_t) =
    use_gpu ? gkyl_array_cu_dev_new : gkyl_array_new;
  struct gkyl_array *prim = alloc(GKYL_DOUBLE, 8, ext.volume);
  struct gkyl_array *jac = alloc(GKYL_DOUBLE, 2, ext.volume);
  struct gkyl_array *bmag = alloc(GKYL_DOUBLE, 2, ext.volume);
  struct gkyl_array *phi = alloc(GKYL_DOUBLE, 2, ext.volume);
  struct gkyl_array *density = alloc(GKYL_DOUBLE, 2, ext.volume);
  struct gkyl_array *state = alloc(GKYL_DOUBLE, 8, ext.volume);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 8, ext.volume);
  struct gkyl_array *scalar = gkyl_array_new(GKYL_DOUBLE, 2, ext.volume);
  long loc = gkyl_range_idx(&local, local.lower);
  struct linear n = {2, 0.1}, jacob = {1.5, 0.2}, magnetic = {3, 0.4};
  struct linear u = {0.3, 0.1}, tpar = {5, 0.2},
                tperp = isotropic ? tpar : (struct linear){7, -0.3};
  struct linear potential = {1.1, -0.2};
  set_linear(host, loc, 0, n);
  set_linear(host, loc, 1, u);
  set_linear(host, loc, 2, tpar);
  set_linear(host, loc, 3, tperp);
  gkyl_array_copy(prim, host);
  set_linear(scalar, loc, 0, jacob);
  gkyl_array_copy(jac, scalar);
  set_linear(scalar, loc, 0, magnetic);
  gkyl_array_copy(bmag, scalar);
  set_linear(scalar, loc, 0, potential);
  gkyl_array_copy(phi, scalar);
  double mass = 2;
  struct linear jn = mul(jacob, n), rho = scale(mass, jn), mom = mul(rho, u);
  struct linear ppar = mul(rho, tpar), pperp = mul(rho, tperp);
  struct linear energy = add(scale(0.5, add(ppar, 1, mul(u, mom))), 1, pperp);
  struct linear conserved[] = {rho, mom, energy, divide(pperp, magnetic)};
  struct gkyl_gk_fluid_moments_inp inp = {
    .basis = &basis,
    .grid = &grid,
    .range = &local,
    .range_ext = &ext,
    .bmag = bmag,
    .mass = mass,
    .charge = -1,
    .num_moments = 4,
    .moment = GKYL_F_MOMENT_BIMAXWELLIAN,
    .use_gpu = use_gpu,
  };
  struct gkyl_gk_fluid_moments *convert = gkyl_gk_fluid_moments_inew(&inp);
  gkyl_gk_fluid_moments_from_prim(convert, &local, prim, jac, state);
  gkyl_array_copy(host, state);
  const double *actual = gkyl_array_cfetch(host, loc);
  for (int c = 0; c < 4; ++c) {
    check_linear(actual + 2 * c, conserved[c]);
  }

  enum gkyl_distribution_moments ids[] = {
    GKYL_F_MOMENT_M0,
    GKYL_F_MOMENT_M1,
    GKYL_F_MOMENT_M2,
    GKYL_F_MOMENT_M2PAR,
    GKYL_F_MOMENT_M2PERP,
    GKYL_F_MOMENT_M3PAR,
    GKYL_F_MOMENT_M3PERP,
    GKYL_F_MOMENT_MAXWELLIAN,
    GKYL_F_MOMENT_BIMAXWELLIAN,
    GKYL_F_MOMENT_M0M1M2,
    GKYL_F_MOMENT_M0M1M2PARM2PERP,
    GKYL_F_MOMENT_HAMILTONIAN
  };
  struct linear m1 = scale(1 / mass, mom), m2 = scale(2 / mass, energy);
  struct linear m2perp = scale(2 / mass, pperp), m2par = add(m2, -1, m2perp);
  struct linear expected[][4] = {
    {jn},
    {m1},
    {m2},
    {m2par},
    {m2perp},
    {mul(u, add(m2par, 2 / mass, ppar))},
    {mul(u, m2perp)},
    {jn, u, scale(1.0 / 3.0, add(tpar, 2, tperp))},
    {jn, u, tpar, tperp},
    {jn, m1, m2},
    {jn, m1, m2par, m2perp},
    {jn, mom, add(energy, -1, mul(potential, jn))}
  };
  for (int m = 0; m < sizeof(ids) / sizeof(ids[0]); ++m) {
    for (int integrated = 0; integrated <= 1; ++integrated) {
      TEST_CASE_("moment=%d, isotropic=%d, integrated=%d", ids[m], isotropic, integrated);
      inp.moment = ids[m];
      inp.integrated = integrated;
      struct gkyl_gk_fluid_moments *up = gkyl_gk_fluid_moments_inew(&inp);
      int nm = gkyl_gk_fluid_moments_num_mom(up);
      struct gkyl_array *out = alloc(GKYL_DOUBLE, nm * (integrated ? 1 : 2), ext.volume);
      struct gkyl_array *result = gkyl_array_new(GKYL_DOUBLE, out->ncomp, ext.volume);
      gkyl_gk_fluid_moments_advance(up, &local, state, phi, out);
      gkyl_array_copy(result, out);
      actual = gkyl_array_cfetch(result, loc);
      for (int c = 0; c < nm; ++c) {
        if (integrated) {
          TEST_CHECK(gkyl_compare(actual[c], 2 * expected[m][c].a, 1e-11));
        } else {
          check_linear(actual + 2 * c, expected[m][c]);
        }
      }
      gkyl_array_release(out);
      gkyl_array_release(result);
      gkyl_gk_fluid_moments_release(up);
    }
  }
  TEST_CASE((const char *)NULL);
  // A spatially varying density ratio multiplies every conserved field weakly.
  struct linear ratio = {1.1, 0.1}, target = mul(jn, ratio);
  set_linear(scalar, loc, 0, target);
  gkyl_array_copy(density, scalar);
  for (int repeat = 0; repeat < 2; ++repeat) {
    gkyl_gk_fluid_moments_set_density(convert, &local, density, state);
    gkyl_array_copy(host, state);
    actual = gkyl_array_cfetch(host, loc);
    for (int c = 0; c < 4; ++c) {
      check_linear(actual + 2 * c, mul(conserved[c], ratio));
    }
  }
  // Constant rescaling preserves recovered primitives even after spatial scaling.
  struct gkyl_array *before = alloc(GKYL_DOUBLE, 8, ext.volume);
  struct gkyl_array *after = alloc(GKYL_DOUBLE, 8, ext.volume);
  gkyl_gk_fluid_moments_advance(convert, &local, state, 0, before);
  gkyl_array_scale(density, 1.4);
  gkyl_gk_fluid_moments_set_density(convert, &local, density, state);
  gkyl_gk_fluid_moments_advance(convert, &local, state, 0, after);
  gkyl_array_copy(host, before);
  double saved[8];
  memcpy(saved, gkyl_array_cfetch(host, loc), sizeof saved);
  gkyl_array_copy(host, after);
  actual = gkyl_array_cfetch(host, loc);
  for (int k = 0; k < 8; ++k) {
    TEST_CHECK(gkyl_compare(actual[k], saved[k] * (k < 2 ? 1.4 : 1), 1e-11));
  }
  gkyl_array_release(before);
  gkyl_array_release(after);
  gkyl_gk_fluid_moments_release(convert);
  gkyl_array_release(prim);
  gkyl_array_release(jac);
  gkyl_array_release(bmag);
  gkyl_array_release(phi);
  gkyl_array_release(density);
  gkyl_array_release(state);
  gkyl_array_release(host);
  gkyl_array_release(scalar);
}

// Verify scaling of every conserved field, including additional fields and untouched ghosts.
// use_gpu selects device arrays for the updater operations.
static void
test_density_scaling(bool use_gpu)
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);
  struct gkyl_rect_grid grid;
  double lower[] = {-1}, upper[] = {1};
  int cells[] = {1}, ghost[] = {1};
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range ext, local;
  gkyl_create_grid_ranges(&grid, ghost, &ext, &local);
  struct gkyl_array *(*alloc)(enum gkyl_elem_type, size_t, size_t) =
    use_gpu ? gkyl_array_cu_dev_new : gkyl_array_new;
  struct gkyl_array *bmag = alloc(GKYL_DOUBLE, basis.num_basis, ext.volume);
  struct gkyl_array *density = alloc(GKYL_DOUBLE, basis.num_basis, ext.volume);
  struct gkyl_array *scalar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, ext.volume);
  gkyl_array_clear(bmag, 0.0);
  gkyl_array_shiftc(bmag, sqrt(2.0), 0);
  double mass = 2.0;
  long loc = gkyl_range_idx(&local, local.lower);
  struct linear rho = {3.0, 0.2}, ratio = {1.1, 0.1};
  set_linear(scalar, loc, 0, scale(1.0 / mass, mul(rho, ratio)));
  gkyl_array_copy(density, scalar);
  for (int num_moments = 4; num_moments <= 6; num_moments += 2) {
    TEST_CASE_("num_moments=%d", num_moments);
    struct gkyl_gk_fluid_moments *up =
      gkyl_gk_fluid_moments_inew(&(struct gkyl_gk_fluid_moments_inp){
        .basis = &basis,
        .grid = &grid,
        .range = &local,
        .range_ext = &ext,
        .bmag = bmag,
        .mass = mass,
        .num_moments = num_moments,
        .moment = GKYL_F_MOMENT_M0,
        .use_gpu = use_gpu,
      });
    struct gkyl_array *state = alloc(GKYL_DOUBLE, num_moments * basis.num_basis, ext.volume);
    struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, state->ncomp, ext.volume);
    gkyl_array_clear(host, -17.0);
    for (int c = 0; c < num_moments; ++c) {
      set_linear(host, loc, c, scale(c + 1, rho));
    }
    gkyl_array_copy(state, host);
    gkyl_gk_fluid_moments_set_density(up, &local, density, state);
    gkyl_array_copy(host, state);
    for (long i = 0; i < ext.volume; ++i) {
      const double *actual = gkyl_array_cfetch(host, i);
      if (i == loc) {
        for (int c = 0; c < num_moments; ++c) {
          check_linear(actual + c * basis.num_basis, mul(scale(c + 1, rho), ratio));
        }
      } else {
        for (int k = 0; k < host->ncomp; ++k) {
          TEST_CHECK(actual[k] == -17.0);
        }
      }
    }
    gkyl_array_release(state);
    gkyl_array_release(host);
    gkyl_gk_fluid_moments_release(up);
  }
  TEST_CASE((const char *)NULL);
  gkyl_array_release(bmag);
  gkyl_array_release(density);
  gkyl_array_release(scalar);
}

static void
test_density_scaling_cpu(void)
{
  test_density_scaling(false);
}

static void
test_spatial_cpu(void)
{
  for (int iso = 0; iso <= 1; ++iso) {
    test_spatial_products(false, iso);
  }
}

static void
test_cpu(void)
{
  test_moments(false);
}
#ifdef GKYL_HAVE_CUDA
static void
test_density_scaling_gpu(void)
{
  test_density_scaling(true);
}

static void
test_gpu(void)
{
  test_moments(true);
  for (int iso = 0; iso <= 1; ++iso) {
    test_spatial_products(true, iso);
  }
}
#endif
TEST_LIST = {
  {"fluid_moments_cpu", test_cpu},
  {"fluid_moments_spatial_products", test_spatial_cpu},
  {"fluid_moments_density_scaling_cpu", test_density_scaling_cpu},
#ifdef GKYL_HAVE_CUDA
  {"fluid_moments_gpu", test_gpu},
  {"fluid_moments_density_scaling_gpu", test_density_scaling_gpu},
#endif
  {NULL, NULL}
};
