// Tests for the fluid part of the PKPM temperature source at the configuration-space nodes. The
// explicit term must add nu_Q (eps_Q - eps) to the energy alone, return the rate nu_Q eps_Q/eps of
// the normalized distributions, add the largest node rate to the CFL rate, and be homogeneous in
// the position-map Jacobian. The implicit step must be backward Euler of the thermal energy with a
// shape factor S such that S eps_new/eps = 1/(1 + nu_Q dt). Inadmissible nodes must stay finite.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_pkpm_temp_source.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define NUM_FLUID 5
#define MAX_NODES 8

typedef void (*map_t)(const double *in, double *GKYL_RESTRICT out);

static const map_t nodal_to_conf[] = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};
static const map_t conf_to_nodal[] = {
  dg_conf_nodal_conf_to_nodal_1x_ser_p1, dg_conf_nodal_conf_to_nodal_2x_ser_p1,
  dg_conf_nodal_conf_to_nodal_3x_ser_p1
};

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

// Node values of one cell.
struct cell_nodal {
  double rho[MAX_NODES], mom[3][MAX_NODES], eps[MAX_NODES]; // density, momentum, thermal energy
  double rate[MAX_NODES], vtsq[MAX_NODES]; // source rate and temperature over the mass
};

// A random admissible state of one cell, scaled by the Jacobian jac, in modal and nodal form.
static void
set_cell(int cdim, int nc, double jac, struct cell_nodal *cell, double *fluid_d)
{
  double energy[MAX_NODES];
  for (int q = 0; q < nc; ++q) {
    cell->rho[q] = jac * (1.0 + 0.5 * rnd());
    double mom_sq = 0.0;
    for (int i = 0; i < 3; ++i) {
      cell->mom[i][q] = cell->rho[q] * rnd();
      mom_sq += cell->mom[i][q] * cell->mom[i][q];
    }
    cell->eps[q] = jac * (1.0 + 0.9 * rnd());
    energy[q] = cell->eps[q] + 0.5 * mom_sq / cell->rho[q];
    cell->rate[q] = 1.0 + 0.5 * fabs(rnd());
    cell->vtsq[q] = 1.0 + 0.9 * rnd();
  }
  nodal_to_conf[cdim - 1](cell->rho, &fluid_d[0]);
  for (int i = 0; i < 3; ++i) {
    nodal_to_conf[cdim - 1](cell->mom[i], &fluid_d[(1 + i) * nc]);
  }
  nodal_to_conf[cdim - 1](energy, &fluid_d[4 * nc]);
}

// Node values of the thermal energy of a modal fluid state.
static void
thermal_energy(int cdim, int nc, const double *fluid_d, double *eps)
{
  double fluid_n[NUM_FLUID][MAX_NODES];
  for (int c = 0; c < NUM_FLUID; ++c) {
    conf_to_nodal[cdim - 1](&fluid_d[c * nc], fluid_n[c]);
  }
  for (int q = 0; q < nc; ++q) {
    double mom_sq =
      fluid_n[1][q] * fluid_n[1][q] + fluid_n[2][q] * fluid_n[2][q] + fluid_n[3][q] * fluid_n[3][q];
    eps[q] = fluid_n[4][q] - 0.5 * mom_sq / fluid_n[0][q];
  }
}

static void
test_temp_source(int cdim, bool implicit)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nc = cbasis.num_basis;

  int clo[3] = {1, 1, 1}, cup[3] = {3, 2, 2};
  struct gkyl_range crange;
  gkyl_range_init(&crange, cdim, clo, cup);

  struct gkyl_array *rate = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *vtsq = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, NUM_FLUID * nc, crange.volume);
  struct gkyl_array *fluid0 = gkyl_array_new(GKYL_DOUBLE, NUM_FLUID * nc, crange.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, NUM_FLUID * nc, crange.volume);
  struct gkyl_array *cfl = gkyl_array_new(GKYL_DOUBLE, 1, crange.volume);
  struct gkyl_array *shape = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct cell_nodal *cells = gkyl_malloc(sizeof(struct cell_nodal[crange.volume]));
  double cfl0 = 0.3, dt = 0.7;
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cfl, cfl0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    // a different Jacobian in each cell
    double jac = 0.5 + 0.25 * loc;
    set_cell(cdim, nc, jac, &cells[loc], gkyl_array_fetch(fluid, loc));
    double *rate_d = gkyl_array_fetch(rate, loc), *vtsq_d = gkyl_array_fetch(vtsq, loc);
    for (int q = 0; q < nc; ++q) {
      rate_d[q] = cells[loc].rate[q];
      vtsq_d[q] = cells[loc].vtsq[q];
    }
  }
  gkyl_array_copy(fluid0, fluid);

  struct gkyl_pkpm_temp_source *up = gkyl_pkpm_temp_source_new(&cbasis, false);
  if (implicit) {
    gkyl_pkpm_temp_source_implicit_advance(up, &crange, dt, rate, vtsq, fluid, shape);
  } else {
    gkyl_pkpm_temp_source_advance(up, &crange, rate, vtsq, fluid, rhs, cfl, shape);
  }

  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    const struct cell_nodal *cell = &cells[loc];
    const double *fluid_d = gkyl_array_cfetch(fluid, loc),
                 *fluid0_d = gkyl_array_cfetch(fluid0, loc);
    const double *shape_d = gkyl_array_cfetch(shape, loc);

    if (implicit) {
      // density and momentum unchanged
      for (int k = 0; k < 4 * nc; ++k) {
        TEST_CHECK(fluid_d[k] == fluid0_d[k]);
      }
      double eps_new[MAX_NODES];
      thermal_energy(cdim, nc, fluid_d, eps_new);
      for (int q = 0; q < nc; ++q) {
        double nu_dt = cell->rate[q] * dt, eps_target = 1.5 * cell->rho[q] * cell->vtsq[q];
        double expected = (cell->eps[q] + nu_dt * eps_target) / (1.0 + nu_dt);
        TEST_CHECK(gkyl_compare_double(eps_new[q], expected, 1.0e-13));
        // the thermal energy stays between its old value and the target
        TEST_CHECK(eps_new[q] >= fmin(cell->eps[q], eps_target) * (1.0 - 1.0e-13));
        TEST_CHECK(eps_new[q] <= fmax(cell->eps[q], eps_target) * (1.0 + 1.0e-13));
        // every pressure component takes the backward-Euler step of the source
        TEST_CHECK(shape_d[q] > 0.0 && shape_d[q] <= 1.0);
        TEST_CHECK(
          gkyl_compare_double(shape_d[q] * eps_new[q] / cell->eps[q], 1.0 / (1.0 + nu_dt), 1.0e-13)
        );
      }
    } else {
      // the state is not modified, and only the energy has a right-hand side
      for (int k = 0; k < NUM_FLUID * nc; ++k) {
        TEST_CHECK(fluid_d[k] == fluid0_d[k]);
      }
      const double *rhs_d = gkyl_array_cfetch(rhs, loc);
      for (int k = 0; k < 4 * nc; ++k) {
        TEST_CHECK(rhs_d[k] == 0.0);
      }
      double rhs_energy[MAX_NODES], rate_max = 0.0;
      conf_to_nodal[cdim - 1](&rhs_d[4 * nc], rhs_energy);
      for (int q = 0; q < nc; ++q) {
        double eps_target = 1.5 * cell->rho[q] * cell->vtsq[q];
        // absolute tolerance: the two energies may nearly cancel
        double expected = cell->rate[q] * (eps_target - cell->eps[q]);
        TEST_CHECK(
          fabs(rhs_energy[q] - expected) < 1.0e-13 * cell->rate[q] * (eps_target + cell->eps[q])
        );
        // rate of the distributions: the Jacobian of the cell cancels
        TEST_CHECK(
          gkyl_compare_double(shape_d[q], cell->rate[q] * eps_target / cell->eps[q], 1.0e-13)
        );
        rate_max = fmax(rate_max, cell->rate[q]);
      }
      const double *cfl_d = gkyl_array_cfetch(cfl, loc);
      TEST_CHECK(gkyl_compare_double(cfl_d[0], cfl0 + rate_max, 1.0e-14));
    }
  }

  gkyl_pkpm_temp_source_release(up);
  gkyl_free(cells);
  gkyl_array_release(rate);
  gkyl_array_release(vtsq);
  gkyl_array_release(fluid);
  gkyl_array_release(fluid0);
  gkyl_array_release(rhs);
  gkyl_array_release(cfl);
  gkyl_array_release(shape);
}

// One 1x cell with a nonpositive density at node 0 and a nonpositive thermal energy at node 1.
static void
test_inadmissible(bool implicit)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  int nc = cbasis.num_basis, clo[1] = {1}, cup[1] = {1};
  struct gkyl_range crange;
  gkyl_range_init(&crange, 1, clo, cup);

  struct gkyl_array *rate = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *vtsq = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, NUM_FLUID * nc, crange.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, NUM_FLUID * nc, crange.volume);
  struct gkyl_array *cfl = gkyl_array_new(GKYL_DOUBLE, 1, crange.volume);
  struct gkyl_array *shape = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cfl, 0.0);

  double nu = 2.0, vtsq_q = 0.5, dt = 0.4;
  double rho[2] = {-0.1, 1.0}, mom[2] = {0.2, 1.0}, zero[2] = {0.0, 0.0};
  // node 1: eps = E - mom^2/(2 rho) = 0.3 - 0.5 = -0.2
  double energy[2] = {1.0, 0.3};
  double *fluid_d = gkyl_array_fetch(fluid, 0);
  dg_conf_nodal_conf_to_modal_1x_ser_p1(rho, &fluid_d[0]);
  dg_conf_nodal_conf_to_modal_1x_ser_p1(mom, &fluid_d[nc]);
  dg_conf_nodal_conf_to_modal_1x_ser_p1(zero, &fluid_d[2 * nc]);
  dg_conf_nodal_conf_to_modal_1x_ser_p1(zero, &fluid_d[3 * nc]);
  dg_conf_nodal_conf_to_modal_1x_ser_p1(energy, &fluid_d[4 * nc]);
  double *rate_d = gkyl_array_fetch(rate, 0), *vtsq_d = gkyl_array_fetch(vtsq, 0);
  for (int q = 0; q < nc; ++q) {
    rate_d[q] = nu;
    vtsq_d[q] = vtsq_q;
  }

  struct gkyl_pkpm_temp_source *up = gkyl_pkpm_temp_source_new(&cbasis, false);
  double incr[2];
  const double *shape_d = gkyl_array_cfetch(shape, 0);
  double eps_target = 1.5 * rho[1] * vtsq_q, eps = -0.2;
  if (implicit) {
    gkyl_pkpm_temp_source_implicit_advance(up, &crange, dt, rate, vtsq, fluid, shape);
    double energy_new[2];
    dg_conf_nodal_conf_to_nodal_1x_ser_p1(&fluid_d[4 * nc], energy_new);
    incr[0] = energy_new[0] - energy[0];
    incr[1] = energy_new[1] - energy[1];
    // no density: nothing happens to the energy; both nodes take the shape factor of T_Q = T
    TEST_CHECK(fabs(incr[0]) < 1.0e-15);
    TEST_CHECK(gkyl_compare_double(incr[1], nu * dt * (eps_target - eps) / (1.0 + nu * dt), 1.0e-13)
    );
    TEST_CHECK(gkyl_compare_double(shape_d[0], 1.0 / (1.0 + nu * dt), 1.0e-14));
    TEST_CHECK(gkyl_compare_double(shape_d[1], 1.0 / (1.0 + nu * dt), 1.0e-14));
  } else {
    gkyl_pkpm_temp_source_advance(up, &crange, rate, vtsq, fluid, rhs, cfl, shape);
    dg_conf_nodal_conf_to_nodal_1x_ser_p1(
      &((const double *)gkyl_array_cfetch(rhs, 0))[4 * nc], incr
    );
    TEST_CHECK(fabs(incr[0]) < 1.0e-15);
    TEST_CHECK(gkyl_compare_double(incr[1], nu * (eps_target - eps), 1.0e-13));
    TEST_CHECK(gkyl_compare_double(shape_d[0], nu, 1.0e-14));
    TEST_CHECK(gkyl_compare_double(shape_d[1], nu, 1.0e-14));
  }

  gkyl_pkpm_temp_source_release(up);
  gkyl_array_release(rate);
  gkyl_array_release(vtsq);
  gkyl_array_release(fluid);
  gkyl_array_release(rhs);
  gkyl_array_release(cfl);
  gkyl_array_release(shape);
}

void
test_explicit_1x(void)
{
  test_temp_source(1, false);
}
void
test_explicit_2x(void)
{
  test_temp_source(2, false);
}
void
test_explicit_3x(void)
{
  test_temp_source(3, false);
}
void
test_implicit_1x(void)
{
  test_temp_source(1, true);
}
void
test_implicit_2x(void)
{
  test_temp_source(2, true);
}
void
test_implicit_3x(void)
{
  test_temp_source(3, true);
}
void
test_inadmissible_explicit(void)
{
  test_inadmissible(false);
}
void
test_inadmissible_implicit(void)
{
  test_inadmissible(true);
}

TEST_LIST = {
  {"test_explicit_1x", test_explicit_1x},
  {"test_explicit_2x", test_explicit_2x},
  {"test_explicit_3x", test_explicit_3x},
  {"test_implicit_1x", test_implicit_1x},
  {"test_implicit_2x", test_implicit_2x},
  {"test_implicit_3x", test_implicit_3x},
  {"test_inadmissible_explicit", test_inadmissible_explicit},
  {"test_inadmissible_implicit", test_inadmissible_implicit},
  {NULL, NULL}
};
