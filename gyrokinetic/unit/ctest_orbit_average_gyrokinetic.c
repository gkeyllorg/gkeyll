#include <acutest.h>

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <gkyl_alloc.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_orbit_average_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_velocity_map.h>

// Tests integrate the evolved Jacobian-weighted distribution independently
// of the updater and do not inspect its reduced matrix or factorization.
struct orbit_test {
  struct gkyl_rect_grid grid, conf_grid, vel_grid;
  struct gkyl_basis basis, conf_basis;
  struct gkyl_range local, local_ext, conf_local, conf_ext, vel_local, vel_ext;
  struct gkyl_velocity_map *vel_map;
  struct gkyl_array *bmag, *phi, *jacobgeo, *mask;
  double mass, charge;
  bool mapped, double_well;
};

static void
orbit_test_bmag(double tm, const double *xn, double *out, void *ctx)
{
  const struct orbit_test *test = ctx;
  out[0] = test->mapped ? 1.0 + 0.05 * xn[0] : 1.0;
}

static void
orbit_test_phi(double tm, const double *xn, double *out, void *ctx)
{
  const struct orbit_test *test = ctx;
  out[0] = test->double_well ? 4.0 * fabs(fabs(xn[0]) - 1.0) : xn[0] * xn[0];
}

static void
orbit_test_jacobgeo(double tm, const double *xn, double *out, void *ctx)
{
  const struct orbit_test *test = ctx;
  out[0] = test->mapped ? 1.0 + 0.15 * xn[0] : 1.0;
}

static void
orbit_test_velocity_map(double tm, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0] + 0.1 * xn[0] * xn[0] * xn[0];
  out[1] = xn[1] + 0.1 * xn[1] * xn[1];
}

static void
orbit_test_opening_bmag(double tm, const double *xn, double *out, void *ctx)
{
  double distance = fabs(xn[0]);
  out[0] = distance < 0.5 ? 10.0 :
           distance < 1.0 ? 10.0 + 10.0 * (distance - 0.5) :
           distance < 1.5 ? 15.0 - 20.0 * (distance - 1.0) :
                            5.0 + 10.0 * (distance - 1.5);
}

static void
orbit_test_opening_phi(double tm, const double *xn, double *out, void *ctx)
{
  const struct orbit_test *test = ctx;
  orbit_test_opening_bmag(tm, xn, out, ctx);
  // At mu=.2, U=mu*B+q*phi=2 everywhere and all nonzero-v orbits escape.
  out[0] = -0.2 * (out[0] - 10.0) / test->charge;
}

static void
orbit_test_init_mu_lower(struct orbit_test *test, bool mapped, bool double_well, double lower_mu)
{
  *test =
    (struct orbit_test){.mass = 1.3, .charge = 0.7, .mapped = mapped, .double_well = double_well};
  const double lower[] = {-2.0, -1.2, lower_mu}, upper[] = {2.0, 1.2, 0.5};
  const int cells[] = {16, 8, 2};
  gkyl_rect_grid_init(&test->grid, 3, lower, upper, cells);
  gkyl_rect_grid_init(&test->conf_grid, 1, lower, upper, cells);
  gkyl_rect_grid_init(&test->vel_grid, 2, lower + 1, upper + 1, cells + 1);
  gkyl_cart_modal_gkhybrid(&test->basis, 1, 2);
  gkyl_cart_modal_serendip(&test->conf_basis, 1, 1);
  gkyl_create_grid_ranges(&test->grid, (int[]){1, 0, 0}, &test->local_ext, &test->local);
  gkyl_create_grid_ranges(&test->conf_grid, (int[]){1}, &test->conf_ext, &test->conf_local);
  gkyl_create_grid_ranges(&test->vel_grid, (int[]){0, 0}, &test->vel_ext, &test->vel_local);
  struct gkyl_mapc2p_inp map = {.mapping = mapped ? orbit_test_velocity_map : 0};
  test->vel_map = gkyl_velocity_map_new(
    map, test->grid, test->vel_grid, test->local, test->local_ext, test->vel_local, test->vel_ext,
    false
  );
  test->bmag = gkyl_array_new(GKYL_DOUBLE, test->conf_basis.num_basis, test->conf_ext.volume);
  test->phi = gkyl_array_new(GKYL_DOUBLE, test->conf_basis.num_basis, test->conf_ext.volume);
  test->jacobgeo = gkyl_array_new(GKYL_DOUBLE, test->conf_basis.num_basis, test->conf_ext.volume);
  test->mask = gkyl_array_new(GKYL_DOUBLE, 1, test->local_ext.volume);
  evalf_t eval[] = {orbit_test_bmag, orbit_test_phi, orbit_test_jacobgeo};
  struct gkyl_array *fields[] = {test->bmag, test->phi, test->jacobgeo};
  for (int fld = 0; fld < 3; ++fld) {
    gkyl_eval_on_nodes *project =
      gkyl_eval_on_nodes_new(&test->conf_grid, &test->conf_basis, 1, eval[fld], test);
    gkyl_eval_on_nodes_advance(project, 0.0, &test->conf_local, fields[fld]);
    gkyl_eval_on_nodes_release(project);
  }
  struct gkyl_loss_cone_mask_gyrokinetic_inp mask_inp = {
    .conf_basis = &test->conf_basis,
    .vel_map = test->vel_map,
    .mass = test->mass,
    .charge = test->charge,
  };
  gkyl_loss_cone_mask_gyrokinetic *mask_up = gkyl_loss_cone_mask_gyrokinetic_inew(&mask_inp);
  gkyl_loss_cone_mask_gyrokinetic_advance(
    mask_up, &test->local, &test->conf_local, test->bmag, test->phi, 0, 0, test->mask
  );
  gkyl_loss_cone_mask_gyrokinetic_release(mask_up);
}

static void
orbit_test_init(struct orbit_test *test, bool mapped, bool double_well)
{
  orbit_test_init_mu_lower(test, mapped, double_well, 0.1);
}

static void
orbit_test_release(struct orbit_test *test)
{
  gkyl_array_release(test->bmag);
  gkyl_array_release(test->phi);
  gkyl_array_release(test->jacobgeo);
  gkyl_array_release(test->mask);
  gkyl_velocity_map_release(test->vel_map);
}

// Evaluate only the public basis/velocity-map representation. This is shared
// by fixture projection and independent high-order diagnostic quadrature.
static void
orbit_test_coordinates(
  const struct orbit_test *test, const int *idx, const double *eta, double *z, double *vpar,
  double *mu, double *hamiltonian, double *jacobian
)
{
  double center[3];
  gkyl_rect_grid_cell_center(&test->grid, idx, center);
  *z = center[0] + 0.5 * test->grid.dx[0] * eta[0];
  long conf_idx = gkyl_range_idx(&test->conf_local, idx);
  long vel_idx = gkyl_range_idx(&test->vel_local, idx + 1);
  long phase_idx = gkyl_range_idx(&test->local, idx);
  const double *map = gkyl_array_cfetch(test->vel_map->vmap, vel_idx);
  *vpar = test->vel_map->vmap_basis->eval_expand(eta + 1, map);
  *mu = test->vel_map->vmap_basis->eval_expand(eta + 2, map + test->vel_map->vmap_basis->num_basis);
  double bmag = test->conf_basis.eval_expand(eta, gkyl_array_cfetch(test->bmag, conf_idx));
  double phi = test->conf_basis.eval_expand(eta, gkyl_array_cfetch(test->phi, conf_idx));
  *hamiltonian = 0.5 * test->mass * (*vpar) * (*vpar) + (*mu) * bmag + test->charge * phi;
  const double *jacobvel = gkyl_array_cfetch(test->vel_map->jacobvel, phase_idx);
  *jacobian = bmag * jacobvel[0] *
              test->conf_basis.eval_expand(eta, gkyl_array_cfetch(test->jacobgeo, conf_idx));
}

enum orbit_test_distribution {
  ORBIT_TEST_SMOOTH,
  ORBIT_TEST_ORBIT_MEAN,
  ORBIT_TEST_INVARIANT,
  ORBIT_TEST_WELLS
};

static struct gkyl_array *
orbit_test_distribution_new(const struct orbit_test *test, enum orbit_test_distribution type)
{
  struct gkyl_array *distf =
    gkyl_array_new(GKYL_DOUBLE, test->basis.num_basis, test->local_ext.volume);
  // Nonzero ghost values make accidental ghost writes observable.
  gkyl_array_clear(distf, 17.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &test->local);
  while (gkyl_range_iter_next(&iter)) {
    double *coeff = gkyl_array_fetch(distf, gkyl_range_idx(&test->local, iter.idx));
    for (int mode = 0; mode < test->basis.num_basis; ++mode) {
      coeff[mode] = 0.0;
    }
    for (int iz = 0; iz < 5; ++iz) {
      for (int iv = 0; iv < 5; ++iv) {
        for (int im = 0; im < 5; ++im) {
          double eta[] = {
            gkyl_gauss_ordinates[5][iz], gkyl_gauss_ordinates[5][iv], gkyl_gauss_ordinates[5][im]
          };
          double weight =
            gkyl_gauss_weights[5][iz] * gkyl_gauss_weights[5][iv] * gkyl_gauss_weights[5][im];
          double z, vpar, mu, hamiltonian, jacobian;
          orbit_test_coordinates(test, iter.idx, eta, &z, &vpar, &mu, &hamiltonian, &jacobian);
          double value = 1.0 + 0.03 * hamiltonian + 0.08 * z + 0.12 * vpar;
          if (type == ORBIT_TEST_ORBIT_MEAN) {
            value = 1.0 + 0.03 * hamiltonian;
          }
          if (type == ORBIT_TEST_INVARIANT) {
            value = 1.0 + 0.1 * hamiltonian + 0.2 * mu;
          }
          if (type == ORBIT_TEST_WELLS) {
            value = z < 0.0 ? 1.0 : 2.0;
          }
          double basis_values[64];
          test->basis.eval(eta, basis_values);
          for (int mode = 0; mode < test->basis.num_basis; ++mode) {
            coeff[mode] += weight * basis_values[mode] * jacobian * value;
          }
        }
      }
    }
  }
  return distf;
}

static void
orbit_test_moments(const struct orbit_test *test, const struct gkyl_array *distf, double moments[3])
{
  moments[0] = moments[1] = moments[2] = 0.0;
  const double volume = test->grid.dx[0] * test->grid.dx[1] * test->grid.dx[2] / 8.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &test->local);
  while (gkyl_range_iter_next(&iter)) {
    const double *coeff = gkyl_array_cfetch(distf, gkyl_range_idx(&test->local, iter.idx));
    for (int iz = 0; iz < 4; ++iz) {
      for (int iv = 0; iv < 4; ++iv) {
        for (int im = 0; im < 4; ++im) {
          double eta[] = {
            gkyl_gauss_ordinates[4][iz], gkyl_gauss_ordinates[4][iv], gkyl_gauss_ordinates[4][im]
          };
          double weight = volume * gkyl_gauss_weights[4][iz] * gkyl_gauss_weights[4][iv] *
                          gkyl_gauss_weights[4][im];
          double z, vpar, mu, hamiltonian, jacobian;
          orbit_test_coordinates(test, iter.idx, eta, &z, &vpar, &mu, &hamiltonian, &jacobian);
          double value = test->basis.eval_expand(eta, coeff);
          moments[0] += weight * value;
          moments[1] += weight * hamiltonian * value;
          moments[2] += weight * mu * value;
        }
      }
    }
  }
}

static double
orbit_test_difference(const struct gkyl_array *lhs, const struct gkyl_array *rhs)
{
  double error = 0.0;
  for (long cell = 0; cell < lhs->size; ++cell) {
    const double *left = gkyl_array_cfetch(lhs, cell), *right = gkyl_array_cfetch(rhs, cell);
    for (int mode = 0; mode < lhs->ncomp; ++mode) {
      error = fmax(error, fabs(left[mode] - right[mode]));
    }
  }
  return error;
}

static void
orbit_test_unchanged_exterior(
  const struct orbit_test *test, const struct gkyl_array *before, const struct gkyl_array *after
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &test->local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long cell = gkyl_range_idx(&test->local_ext, iter.idx);
    const double *mask = gkyl_array_cfetch(test->mask, cell);
    bool ghost = iter.idx[0] < test->local.lower[0] || iter.idx[0] > test->local.upper[0];
    if (!ghost && mask[0] != 0.0) {
      continue;
    }
    const double *initial = gkyl_array_cfetch(before, cell),
                 *final = gkyl_array_cfetch(after, cell);
    for (int mode = 0; mode < test->basis.num_basis; ++mode) {
      TEST_CHECK(initial[mode] == final[mode]);
    }
  }
}

static double
orbit_test_minimum(const struct orbit_test *test, const struct gkyl_array *distf)
{
  double minimum = DBL_MAX;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &test->local);
  while (gkyl_range_iter_next(&iter)) {
    const double *coeff = gkyl_array_cfetch(distf, gkyl_range_idx(&test->local, iter.idx));
    for (int iz = 0; iz < 4; ++iz) {
      for (int iv = 0; iv < 4; ++iv) {
        for (int im = 0; im < 2; ++im) {
          double eta[] = {
            gkyl_gauss_ordinates[4][iz], gkyl_gauss_ordinates[4][iv], gkyl_gauss_ordinates[2][im]
          };
          minimum = fmin(minimum, test->basis.eval_expand(eta, coeff));
        }
      }
    }
    for (int corner = 0; corner < 8; ++corner) {
      double eta[] = {corner & 1 ? 1.0 : -1.0, corner & 2 ? 1.0 : -1.0, corner & 4 ? 1.0 : -1.0};
      minimum = fmin(minimum, test->basis.eval_expand(eta, coeff));
    }
  }
  return minimum;
}

static void
orbit_test_check_moments(
  const struct orbit_test *test, const struct gkyl_array *before, const struct gkyl_array *after
)
{
  double initial[3], final[3];
  orbit_test_moments(test, before, initial);
  orbit_test_moments(test, after, final);
  for (int moment = 0; moment < 3; ++moment) {
    TEST_CHECK(gkyl_compare(final[moment], initial[moment], 2e-11));
    TEST_MSG("Moment %d: before %.17g, after %.17g", moment, initial[moment], final[moment]);
  }
}

static struct gkyl_orbit_average_gyrokinetic_inp
orbit_test_input(const struct orbit_test *test, bool use_gpu, bool positivity)
{
  return (struct gkyl_orbit_average_gyrokinetic_inp){
    .grid = &test->grid,
    .conf_basis = &test->conf_basis,
    .phase_basis = &test->basis,
    .phase_range = &test->local,
    .conf_range = &test->conf_local,
    .vel_map = test->vel_map,
    .bmag = test->bmag,
    .phi = test->phi,
    .jacobgeo = test->jacobgeo,
    .mask = test->mask,
    .mass = test->mass,
    .charge = test->charge,
    .num_energy_nodes = 6,
    .positivity = positivity,
    .use_gpu = use_gpu,
  };
}

static gkyl_orbit_average_gyrokinetic *
orbit_test_updater(const struct orbit_test *test, bool use_gpu, bool positivity)
{
  struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(test, use_gpu, positivity);
  enum gkyl_orbit_average_gyrokinetic_status status;
  gkyl_orbit_average_gyrokinetic *up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_ASSERT(up != 0);
  TEST_CHECK(status == GKYL_ORBIT_AVERAGE_SUCCESS);
  struct gkyl_orbit_average_gyrokinetic_info info = gkyl_orbit_average_gyrokinetic_get_info(up);
  TEST_CHECK(info.num_selected_cells > 0);
  TEST_MSG(
    "Selected cells: %ld; rejected cells: %ld", info.num_selected_cells, info.num_rejected_cells
  );
  TEST_CHECK(info.num_modes > 1);
  return up;
}

static struct gkyl_orbit_average_gyrokinetic_result
orbit_test_update(
  gkyl_orbit_average_gyrokinetic *up, bool use_gpu, double dt, double frequency,
  const struct gkyl_array *initial, struct gkyl_array *final
)
{
  struct gkyl_array *working = use_gpu ?
                                 gkyl_array_cu_dev_new(GKYL_DOUBLE, initial->ncomp, initial->size) :
                                 gkyl_array_acquire(final);
  gkyl_array_copy(working, initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    gkyl_orbit_average_gyrokinetic_advance(up, dt, frequency, working);
  if (use_gpu) {
    gkyl_array_copy(final, working);
  }
  gkyl_array_release(working);
  return status;
}

// Use a large finite rate to reach the projection limit at double precision.
static void
orbit_test_projection(bool use_gpu, bool mapped)
{
  struct orbit_test test;
  orbit_test_init(&test, mapped, false);
  gkyl_orbit_average_gyrokinetic *up = orbit_test_updater(&test, use_gpu, false);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *projected = gkyl_array_clone(initial);
  struct gkyl_array *twice = gkyl_array_clone(initial);
  struct gkyl_array *finite_dt = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(up, use_gpu, 1.0, 1e20, initial, projected);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_difference(initial, projected) > 1e-4);
  orbit_test_check_moments(&test, initial, projected);
  orbit_test_unchanged_exterior(&test, initial, projected);
  if (!mapped) {
    // Symmetry gives zero exact bounce mean for both z and v_parallel.
    // The affine-H remainder is exactly represented by the DG basis.
    struct gkyl_array *mean = orbit_test_distribution_new(&test, ORBIT_TEST_ORBIT_MEAN);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &test.local);
    while (gkyl_range_iter_next(&iter)) {
      long cell = gkyl_range_idx(&test.local, iter.idx);
      const double *input = gkyl_array_cfetch(initial, cell);
      const double *actual = gkyl_array_cfetch(projected, cell);
      const double *expected = gkyl_array_cfetch(mean, cell);
      double change = 0.0;
      for (int mode = 0; mode < test.basis.num_basis; ++mode) {
        change = fmax(change, fabs(actual[mode] - input[mode]));
      }
      if (change > 1e-10) {
        for (int mode = 0; mode < test.basis.num_basis; ++mode) {
          TEST_CHECK(fabs(actual[mode] - expected[mode]) < 2e-10);
          TEST_MSG(
            "Orbit mean at cell %ld, mode %d: expected %.17g, actual %.17g", cell, mode,
            expected[mode], actual[mode]
          );
        }
      }
    }
    gkyl_array_release(mean);
  }

  status = orbit_test_update(up, use_gpu, 1.0, 1e20, projected, twice);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_difference(projected, twice) < 2e-10);
  TEST_MSG("Idempotence error: %.17g", orbit_test_difference(projected, twice));

  // The finite-step oracle uses the defining backward-Euler resolvent,
  // independently of the updater's implementation of its blend factor.
  const double dt = 0.25, frequency = 3.0, theta = 1.0 / (1.0 + dt * frequency);
  status = orbit_test_update(up, use_gpu, dt, frequency, initial, finite_dt);
  TEST_CHECK(status.success);
  for (long cell = 0; cell < initial->size; ++cell) {
    const double *input = gkyl_array_cfetch(initial, cell),
                 *limit = gkyl_array_cfetch(projected, cell);
    const double *actual = gkyl_array_cfetch(finite_dt, cell);
    for (int mode = 0; mode < test.basis.num_basis; ++mode) {
      TEST_CHECK(fabs(actual[mode] - limit[mode] - theta * (input[mode] - limit[mode])) < 2e-11);
    }
  }
  orbit_test_check_moments(&test, initial, finite_dt);

  status = orbit_test_update(up, use_gpu, 1e300, 1e300, initial, twice);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_difference(projected, twice) < 2e-11);
  status = orbit_test_update(up, use_gpu, 1.0, 0.0, initial, twice);
  TEST_CHECK(status.success && orbit_test_difference(initial, twice) == 0.0);
  status = orbit_test_update(up, use_gpu, 0.0, 1e20, initial, twice);
  TEST_CHECK(status.success && orbit_test_difference(initial, twice) == 0.0);
  status = orbit_test_update(up, use_gpu, -1.0, 1.0, initial, twice);
  TEST_CHECK(!status.success && orbit_test_difference(initial, twice) == 0.0);
  status = orbit_test_update(up, use_gpu, 1.0, -1.0, initial, twice);
  TEST_CHECK(!status.success && orbit_test_difference(initial, twice) == 0.0);

  // The explicit infinite-frequency option overrides the numeric rate.
  struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(&test, use_gpu, false);
  inp.frequency_is_infinite = true;
  enum gkyl_orbit_average_gyrokinetic_status construction;
  gkyl_orbit_average_gyrokinetic *infinite =
    gkyl_orbit_average_gyrokinetic_inew(&inp, &construction);
  TEST_ASSERT(infinite != 0);
  TEST_CHECK(construction == GKYL_ORBIT_AVERAGE_SUCCESS);
  status = orbit_test_update(infinite, use_gpu, 1.0, 0.0, initial, twice);
  TEST_CHECK(status.success && status.min_blend == 1.0);
  TEST_CHECK(orbit_test_difference(projected, twice) < 2e-11);
  status = orbit_test_update(infinite, use_gpu, 1.0, 3.0, initial, twice);
  TEST_CHECK(status.success && status.min_blend == 1.0);
  TEST_CHECK(orbit_test_difference(projected, twice) < 2e-11);
  status = orbit_test_update(infinite, use_gpu, 0.0, 0.0, initial, twice);
  TEST_CHECK(status.success && orbit_test_difference(initial, twice) == 0.0);
  gkyl_orbit_average_gyrokinetic_release(infinite);

  if (use_gpu) {
    gkyl_orbit_average_gyrokinetic *cpu = orbit_test_updater(&test, false, false);
    status = orbit_test_update(cpu, false, 1.0, 1e20, initial, twice);
    TEST_CHECK(status.success);
    TEST_CHECK(orbit_test_difference(projected, twice) < 2e-10);
    gkyl_orbit_average_gyrokinetic_release(cpu);
  }
  gkyl_array_release(initial);
  gkyl_array_release(projected);
  gkyl_array_release(twice);
  gkyl_array_release(finite_dt);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
orbit_test_stationarity(bool use_gpu, bool double_well)
{
  struct orbit_test test;
  orbit_test_init(&test, false, double_well);
  gkyl_orbit_average_gyrokinetic *up = orbit_test_updater(&test, use_gpu, false);
  struct gkyl_array *initial =
    orbit_test_distribution_new(&test, double_well ? ORBIT_TEST_WELLS : ORBIT_TEST_INVARIANT);
  struct gkyl_array *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(up, use_gpu, 1.0, 1e20, initial, final);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_difference(initial, final) < 2e-10);
  TEST_MSG("Stationary-state error: %.17g", orbit_test_difference(initial, final));
  orbit_test_check_moments(&test, initial, final);
  orbit_test_unchanged_exterior(&test, initial, final);
  gkyl_array_release(initial);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
orbit_test_positivity(bool use_gpu)
{
  struct orbit_test test;
  orbit_test_init(&test, false, false);
  gkyl_orbit_average_gyrokinetic *limited = orbit_test_updater(&test, use_gpu, true);
  gkyl_orbit_average_gyrokinetic *unlimited = orbit_test_updater(&test, use_gpu, false);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &test.local);
  while (gkyl_range_iter_next(&iter)) {
    double *coeff = gkyl_array_fetch(initial, gkyl_range_idx(&test.local, iter.idx));
    bool spike = iter.idx[0] == 7 && iter.idx[1] == 4;
    coeff[0] = sqrt(8.0) * (spike ? 1.0 : 0.001);
    for (int mode = 1; mode < test.basis.num_basis; ++mode) {
      coeff[mode] = 0.0;
    }
  }
  struct gkyl_array *projected = gkyl_array_clone(initial), *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(unlimited, use_gpu, 1.0, 1e20, initial, projected);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_minimum(&test, projected) < -1e-6);
  TEST_MSG("Unconstrained minimum: %.17g", orbit_test_minimum(&test, projected));
  status = orbit_test_update(limited, use_gpu, 1.0, 1e20, initial, final);
  TEST_CHECK(status.success);
  TEST_CHECK(status.num_limited_blocks > 0);
  TEST_CHECK(status.min_blend < 1.0 && status.min_blend >= 0.0);
  TEST_CHECK(orbit_test_minimum(&test, final) >= -2e-12);
  orbit_test_check_moments(&test, initial, final);
  orbit_test_unchanged_exterior(&test, initial, final);

  // Failed positive-input precondition must not partially modify another block.
  gkyl_range_iter_init(&iter, &test.local);
  while (gkyl_range_iter_next(&iter)) {
    double *coeff = gkyl_array_fetch(initial, gkyl_range_idx(&test.local, iter.idx));
    if (iter.idx[2] == test.local.upper[2]) {
      coeff[0] = -sqrt(8.0);
    }
  }
  status = orbit_test_update(limited, use_gpu, 1.0, 1e20, initial, final);
  TEST_CHECK(!status.success);
  TEST_CHECK(orbit_test_difference(initial, final) == 0.0);
  gkyl_array_release(initial);
  gkyl_array_release(projected);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(limited);
  gkyl_orbit_average_gyrokinetic_release(unlimited);
  orbit_test_release(&test);
}

static void
orbit_test_empty_and_invalid(void)
{
  struct orbit_test test;
  orbit_test_init(&test, false, false);
  struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(&test, false, false);
  enum gkyl_orbit_average_gyrokinetic_status status;
  inp.mass = -1.0;
  gkyl_orbit_average_gyrokinetic *up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_CHECK(up == 0 && status == GKYL_ORBIT_AVERAGE_INVALID_INPUT);
  inp.mass = test.mass;
  inp.num_energy_nodes = 1;
  up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_CHECK(up == 0 && status == GKYL_ORBIT_AVERAGE_INVALID_INPUT);
  inp.num_energy_nodes = 6;
  inp.separatrix_margin = -1.0;
  up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_CHECK(up == 0 && status == GKYL_ORBIT_AVERAGE_INVALID_INPUT);
  inp.separatrix_margin = 0.0;
  gkyl_array_scale(test.jacobgeo, -1.0);
  up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_CHECK(up == 0 && status == GKYL_ORBIT_AVERAGE_INVALID_INPUT);
  gkyl_array_scale(test.jacobgeo, -1.0);
  gkyl_array_clear(test.mask, 0.0);
  up = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_ASSERT(up != 0);
  TEST_CHECK(status == GKYL_ORBIT_AVERAGE_SUCCESS);
  TEST_CHECK(gkyl_orbit_average_gyrokinetic_get_info(up).num_selected_cells == 0);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result update =
    orbit_test_update(up, false, 1.0, 1e20, initial, final);
  TEST_CHECK(update.success && orbit_test_difference(initial, final) == 0.0);
  gkyl_array_release(initial);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
orbit_test_separatrix_margin(void)
{
  struct orbit_test test;
  orbit_test_init(&test, false, false);
  struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(&test, false, false);
  enum gkyl_orbit_average_gyrokinetic_status status;
  gkyl_orbit_average_gyrokinetic *unbuffered = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_ASSERT(unbuffered != 0);
  long before = gkyl_orbit_average_gyrokinetic_get_info(unbuffered).num_selected_cells;
  // A deliberately wide energy margin distinguishes the criterion from
  // the independent finite-velocity-domain guard in this small fixture.
  inp.separatrix_margin = 2.4;
  gkyl_orbit_average_gyrokinetic *buffered = gkyl_orbit_average_gyrokinetic_inew(&inp, &status);
  TEST_ASSERT(buffered != 0);
  long after = gkyl_orbit_average_gyrokinetic_get_info(buffered).num_selected_cells;
  TEST_CHECK(after > 0 && after < before);
  TEST_MSG("Unbuffered cells %ld, buffered cells %ld", before, after);
  gkyl_orbit_average_gyrokinetic_release(unbuffered);
  gkyl_orbit_average_gyrokinetic_release(buffered);
  orbit_test_release(&test);
}

static void
orbit_test_permissive_mask(bool use_gpu)
{
  struct orbit_test test;
  orbit_test_init(&test, false, false);
  // An overpermissive caller mask must not bypass the closure checks.
  gkyl_array_clear(test.mask, 1.0);
  gkyl_orbit_average_gyrokinetic *up = orbit_test_updater(&test, use_gpu, false);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(up, use_gpu, 1.0, 1e20, initial, final);
  TEST_CHECK(status.success);
  TEST_CHECK(gkyl_orbit_average_gyrokinetic_get_info(up).num_rejected_cells > 0);
  int index[] = {test.local.lower[0], test.local.lower[1], test.local.lower[2]};
  long cell = gkyl_range_idx(&test.local, index);
  const double *before = gkyl_array_cfetch(initial, cell), *after = gkyl_array_cfetch(final, cell);
  for (int mode = 0; mode < test.basis.num_basis; ++mode) {
    TEST_CHECK(before[mode] == after[mode]);
  }
  orbit_test_check_moments(&test, initial, final);
  gkyl_array_release(initial);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
orbit_test_midpoint_opening(bool use_gpu)
{
  struct orbit_test test;
  orbit_test_init(&test, false, false);
  evalf_t eval[] = {orbit_test_opening_bmag, orbit_test_opening_phi};
  struct gkyl_array *fields[] = {test.bmag, test.phi};
  for (int field = 0; field < 2; ++field) {
    gkyl_eval_on_nodes *project =
      gkyl_eval_on_nodes_new(&test.conf_grid, &test.conf_basis, 1, eval[field], &test);
    gkyl_eval_on_nodes_advance(project, 0.0, &test.conf_local, fields[field]);
    gkyl_eval_on_nodes_release(project);
  }
  struct gkyl_loss_cone_mask_gyrokinetic_inp mask_inp = {
    .conf_basis = &test.conf_basis,
    .vel_map = test.vel_map,
    .mass = test.mass,
    .charge = test.charge,
  };
  gkyl_loss_cone_mask_gyrokinetic *mask = gkyl_loss_cone_mask_gyrokinetic_inew(&mask_inp);
  gkyl_loss_cone_mask_gyrokinetic_advance(
    mask, &test.local, &test.conf_local, test.bmag, test.phi, 0, 0, test.mask
  );
  gkyl_loss_cone_mask_gyrokinetic_release(mask);
  // This cell has z in [-.25,0], v in [.3,.6] and mu in [.1,.3].
  // Its maximum kinetic energy is .234, below the corner barriers .5
  // and the Gauss-mu barriers .5/sqrt(3), but above the midpoint barrier0.
  int index[] = {8, 6, 1};
  long cell = gkyl_range_idx(&test.local, index);
  TEST_CHECK(((const double *)gkyl_array_cfetch(test.mask, cell))[0] == 1.0);
  struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(&test, use_gpu, false);
  enum gkyl_orbit_average_gyrokinetic_status construction;
  gkyl_orbit_average_gyrokinetic *up = gkyl_orbit_average_gyrokinetic_inew(&inp, &construction);
  TEST_ASSERT(up != 0 && construction == GKYL_ORBIT_AVERAGE_SUCCESS);
  TEST_CHECK(gkyl_orbit_average_gyrokinetic_get_info(up).num_rejected_cells > 0);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(up, use_gpu, 1.0, 1e20, initial, final);
  TEST_CHECK(status.success);
  const double *before = gkyl_array_cfetch(initial, cell), *after = gkyl_array_cfetch(final, cell);
  for (int mode = 0; mode < test.basis.num_basis; ++mode) {
    TEST_CHECK(before[mode] == after[mode]);
  }
  orbit_test_check_moments(&test, initial, final);
  gkyl_array_release(initial);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
orbit_test_zero_mu(void)
{
  struct orbit_test test;
  orbit_test_init_mu_lower(&test, true, false, 0.0);
  gkyl_orbit_average_gyrokinetic *up = orbit_test_updater(&test, false, false);
  TEST_CHECK(gkyl_orbit_average_gyrokinetic_get_info(up).num_blocks == 2);
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *final = gkyl_array_clone(initial);
  struct gkyl_orbit_average_gyrokinetic_result status =
    orbit_test_update(up, false, 1.0, 1e20, initial, final);
  TEST_CHECK(status.success);
  TEST_CHECK(orbit_test_difference(initial, final) > 1e-4);
  orbit_test_check_moments(&test, initial, final);
  orbit_test_unchanged_exterior(&test, initial, final);
  gkyl_array_release(initial);
  gkyl_array_release(final);
  gkyl_orbit_average_gyrokinetic_release(up);
  orbit_test_release(&test);
}

static void
test_projection(void)
{
  orbit_test_projection(false, false);
}
static void
test_mapped_projection(void)
{
  orbit_test_projection(false, true);
}
static void
test_stationarity(void)
{
  orbit_test_stationarity(false, false);
}
static void
test_disconnected_wells(void)
{
  orbit_test_stationarity(false, true);
}
static void
test_positivity(void)
{
  orbit_test_positivity(false);
}
static void
test_permissive_mask(void)
{
  orbit_test_permissive_mask(false);
}
static void
test_midpoint_opening(void)
{
  orbit_test_midpoint_opening(false);
}

#ifdef GKYL_HAVE_CUDA
static void
test_projection_cu(void)
{
  orbit_test_projection(true, false);
}
static void
test_mapped_projection_cu(void)
{
  orbit_test_projection(true, true);
}
static void
test_stationarity_cu(void)
{
  orbit_test_stationarity(true, false);
}
static void
test_disconnected_wells_cu(void)
{
  orbit_test_stationarity(true, true);
}
static void
test_positivity_cu(void)
{
  orbit_test_positivity(true);
}
static void
test_permissive_mask_cu(void)
{
  orbit_test_permissive_mask(true);
}
static void
test_midpoint_opening_cu(void)
{
  orbit_test_midpoint_opening(true);
}
#endif

// Exercise real concurrent collectives without requiring mpirun for unit tests.
// Unequal partitions include the two boundary slabs with no selected cells.
struct orbit_test_collective {
  pthread_barrier_t barrier;
  double *partial[4], *total;
  bool use_gpu, zero_coeff[4];
};
struct orbit_test_rank {
  struct orbit_test_collective *collective;
  int rank;
  gkyl_orbit_average_gyrokinetic *up;
  struct gkyl_array *f;
  struct gkyl_orbit_average_gyrokinetic_result result;
};

static int
orbit_test_reduce(void *ctx, bool sum, long count, const double *input, double *output)
{
  struct orbit_test_rank *rank = ctx;
  struct orbit_test_collective *comm = rank->collective;
  size_t bytes = count * sizeof(double);
#ifdef GKYL_HAVE_CUDA
  if (comm->use_gpu) {
    gkyl_cu_memcpy(comm->partial[rank->rank], input, bytes, GKYL_CU_MEMCPY_D2H);
  } else
#endif
    memcpy(comm->partial[rank->rank], input, bytes);
  if (sum) {
    comm->zero_coeff[rank->rank] = true;
    for (long i = 0; i < count; ++i) {
      comm->zero_coeff[rank->rank] &= comm->partial[rank->rank][i] == 0.0;
    }
  }
  pthread_barrier_wait(&comm->barrier);
  if (rank->rank == 0) {
    for (long i = 0; i < count; ++i) {
      double value = sum ? 0.0 : DBL_MAX;
      for (int r = 0; r < 4; ++r) {
        value = sum ? value + comm->partial[r][i] : fmin(value, comm->partial[r][i]);
      }
      comm->total[i] = value;
    }
  }
  pthread_barrier_wait(&comm->barrier);
#ifdef GKYL_HAVE_CUDA
  if (comm->use_gpu) {
    gkyl_cu_memcpy(output, comm->total, bytes, GKYL_CU_MEMCPY_H2D);
  } else
#endif
    memcpy(output, comm->total, bytes);
  pthread_barrier_wait(&comm->barrier);
  return 0;
}

static void *
orbit_test_rank_advance(void *ctx)
{
  struct orbit_test_rank *rank = ctx;
  rank->result = gkyl_orbit_average_gyrokinetic_advance(rank->up, 1.0, 1e20, rank->f);
  return 0;
}

static void
orbit_test_distributed(bool use_gpu)
{
  struct orbit_test test;
  orbit_test_init(&test, true, false);
  gkyl_orbit_average_gyrokinetic *reference = orbit_test_updater(&test, use_gpu, true);
  struct gkyl_orbit_average_gyrokinetic_info info =
    gkyl_orbit_average_gyrokinetic_get_info(reference);
  struct orbit_test_collective comm = {.use_gpu = use_gpu};
  pthread_barrier_init(&comm.barrier, 0, 4);
  long count = info.num_modes > info.num_blocks ? info.num_modes : info.num_blocks;
  comm.total = gkyl_malloc(count * sizeof(double));
  struct orbit_test_rank ranks[4];
  struct gkyl_range locals[4], exts[4];
  int starts[] = {1, 2, 9, 16}, ends[] = {1, 8, 15, 16};
  for (int r = 0; r < 4; ++r) {
    comm.partial[r] = gkyl_malloc(count * sizeof(double));
    int lower[] = {starts[r], 1, 1}, upper[] = {ends[r], 8, 2};
    gkyl_range_init(&locals[r], 3, lower, upper);
    gkyl_range_extend(&exts[r], &locals[r], (int[]){1, 0, 0}, (int[]){1, 0, 0});
    gkyl_sub_range_init(&locals[r], &exts[r], lower, upper);
    ranks[r] = (struct orbit_test_rank){.collective = &comm, .rank = r};
    struct gkyl_orbit_average_gyrokinetic_inp inp = orbit_test_input(&test, use_gpu, true);
    inp.local_range = &locals[r];
    inp.reduce = orbit_test_reduce;
    inp.reduce_ctx = &ranks[r];
    ranks[r].up = gkyl_orbit_average_gyrokinetic_inew(&inp, 0);
    TEST_ASSERT(ranks[r].up != 0);
    ranks[r].f = use_gpu ?
                   gkyl_array_cu_dev_new(GKYL_DOUBLE, test.basis.num_basis, exts[r].volume) :
                   gkyl_array_new(GKYL_DOUBLE, test.basis.num_basis, exts[r].volume);
  }
  struct gkyl_array *initial = orbit_test_distribution_new(&test, ORBIT_TEST_SMOOTH);
  struct gkyl_array *expected = gkyl_array_clone(initial), *actual = gkyl_array_clone(initial);
  for (int invalid = 0; invalid < 2; ++invalid) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &test.local);
    while (gkyl_range_iter_next(&iter)) {
      double *modal = gkyl_array_fetch(initial, gkyl_range_idx(&test.local, iter.idx));
      modal[0] = sqrt(8.0) * (iter.idx[0] == 7 && iter.idx[1] == 4 ? 1.0 : 0.001);
      if (invalid && iter.idx[0] == 7) {
        modal[0] = -1.0;
      }
      for (int m = 1; m < test.basis.num_basis; ++m) {
        modal[m] = 0.0;
      }
    }
    gkyl_array_copy(actual, initial);
    struct gkyl_orbit_average_gyrokinetic_result result =
      orbit_test_update(reference, use_gpu, 1.0, 1e20, initial, expected);
    TEST_CHECK(result.success == !invalid);
    if (!invalid) {
      TEST_CHECK(result.num_limited_blocks > 0);
    }
    pthread_t threads[4];
    for (int r = 0; r < 4; ++r) {
      struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, test.basis.num_basis, exts[r].volume);
      gkyl_range_iter_init(&iter, &locals[r]);
      while (gkyl_range_iter_next(&iter)) {
        memcpy(
          gkyl_array_fetch(host, gkyl_range_idx(&locals[r], iter.idx)),
          gkyl_array_cfetch(initial, gkyl_range_idx(&test.local, iter.idx)),
          test.basis.num_basis * sizeof(double)
        );
      }
      gkyl_array_copy(ranks[r].f, host);
      gkyl_array_release(host);
      TEST_ASSERT(pthread_create(&threads[r], 0, orbit_test_rank_advance, &ranks[r]) == 0);
    }
    for (int r = 0; r < 4; ++r) {
      pthread_join(threads[r], 0);
      TEST_CHECK(ranks[r].result.success == !invalid);
      if (!invalid) {
        TEST_CHECK(fabs(ranks[r].result.min_blend - result.min_blend) < 1e-12);
        TEST_CHECK(ranks[r].result.num_limited_blocks == result.num_limited_blocks);
      }
      struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, test.basis.num_basis, exts[r].volume);
      gkyl_array_copy(host, ranks[r].f);
      gkyl_range_iter_init(&iter, &locals[r]);
      while (gkyl_range_iter_next(&iter)) {
        memcpy(
          gkyl_array_fetch(actual, gkyl_range_idx(&test.local, iter.idx)),
          gkyl_array_cfetch(host, gkyl_range_idx(&locals[r], iter.idx)),
          test.basis.num_basis * sizeof(double)
        );
      }
      gkyl_array_release(host);
    }
    TEST_CHECK(comm.zero_coeff[0] && comm.zero_coeff[3]);
    TEST_CHECK(orbit_test_difference(expected, actual) < 2e-12);
    if (invalid) {
      TEST_CHECK(orbit_test_difference(initial, actual) == 0.0);
    } else {
      orbit_test_check_moments(&test, initial, actual);
    }
  }
  for (int r = 0; r < 4; ++r) {
    gkyl_orbit_average_gyrokinetic_release(ranks[r].up);
    gkyl_array_release(ranks[r].f);
    gkyl_free(comm.partial[r]);
  }
  pthread_barrier_destroy(&comm.barrier);
  gkyl_free(comm.total);
  gkyl_array_release(initial);
  gkyl_array_release(expected);
  gkyl_array_release(actual);
  gkyl_orbit_average_gyrokinetic_release(reference);
  orbit_test_release(&test);
}

static void
test_distributed(void)
{
  orbit_test_distributed(false);
}
#ifdef GKYL_HAVE_CUDA
static void
test_distributed_cu(void)
{
  orbit_test_distributed(true);
}
#endif

TEST_LIST = {
  {"distributed", test_distributed},
#ifdef GKYL_HAVE_CUDA
  {"distributed_cu", test_distributed_cu},
#endif
  {"projection", test_projection},
  {"mapped_projection", test_mapped_projection},
  {"stationarity", test_stationarity},
  {"disconnected_wells", test_disconnected_wells},
  {"positivity", test_positivity},
  {"empty_and_invalid", orbit_test_empty_and_invalid},
  {"permissive_mask", test_permissive_mask},
  {"separatrix_margin", orbit_test_separatrix_margin},
  {"midpoint_opening", test_midpoint_opening},
  {"zero_mu", orbit_test_zero_mu},
#ifdef GKYL_HAVE_CUDA
  {"projection_cu", test_projection_cu},
  {"mapped_projection_cu", test_mapped_projection_cu},
  {"stationarity_cu", test_stationarity_cu},
  {"disconnected_wells_cu", test_disconnected_wells_cu},
  {"positivity_cu", test_positivity_cu},
  {"permissive_mask_cu", test_permissive_mask_cu},
  {"midpoint_opening_cu", test_midpoint_opening_cu},
#endif
  {0, 0}
};
