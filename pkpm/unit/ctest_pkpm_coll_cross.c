// The cross-species collisions of the PKPM collision operators at the configuration-space nodes:
// the node values of the fluid state, the exchange coefficient of a pair, the rate a partner
// adds to the kinetic operator (LBO and BGK) and its floor, the backward-Euler step of the
// momentum and energy exchange over two and three species (conservation, the backward-Euler
// equations at the nodes, positive temperatures), and the shape factor of the implicit BGK step.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_pkpm_coll_cross.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define NCELL 3 // configuration cells per direction
#define NFLUID 5 // blocks of the fluid state
#define NPRIM 5 // node fields [n, u (3), vtsq]
#define NSPEC 3 // largest number of species
#define NNODE 8 // largest number of nodes

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

struct cross_data {
  int num_nodes;
  int num_species;
  double mass[NSPEC];
  struct gkyl_range crange;
  struct gkyl_basis cbasis;
  struct gkyl_array *jacob;
  struct gkyl_array *fluid[NSPEC], *prim[NSPEC];
  struct gkyl_array *nu[NSPEC][NSPEC]; // frequency of species s on species r
};

// Fluid states (times the Jacobian) with flows of the order of the thermal speeds and
// frequencies that do not satisfy m_s n_s nu_sr = m_r n_r nu_rs, all varying in the cell.
static struct cross_data
cross_data_new(int cdim, int num_species)
{
  struct cross_data dat = {.mass = {1.0, 25.0, 4.0}, .num_species = num_species};
  gkyl_cart_modal_serendip(&dat.cbasis, cdim, 1);
  int nc = dat.num_nodes = dat.cbasis.num_basis;
  double mean = pow(sqrt(2.0), cdim); // expansion coefficient of the constant 1

  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = 1;
    upper[d] = NCELL;
  }
  gkyl_range_init(&dat.crange, cdim, lower, upper);

  dat.jacob = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
  for (int s = 0; s < num_species; ++s) {
    dat.fluid[s] = gkyl_array_new(GKYL_DOUBLE, NFLUID * nc, dat.crange.volume);
    dat.prim[s] = gkyl_array_new(GKYL_DOUBLE, NPRIM * nc, dat.crange.volume);
    for (int r = 0; r < num_species; ++r) {
      dat.nu[s][r] = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
    }
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&dat.crange, iter.idx);
    double *jacob_d = gkyl_array_fetch(dat.jacob, loc);
    jacob_d[0] = 1.5 + 0.4 * rnd();
    for (int s = 0; s < num_species; ++s) {
      double *fluid_d = gkyl_array_fetch(dat.fluid[s], loc);
      double mass = dat.mass[s];
      for (int i = 0; i < nc; ++i) {
        fluid_d[i] = mass * ((i == 0 ? (1.0 + 0.5 * s) * mean : 0.0) + 0.1 * rnd());
        for (int k = 0; k < 3; ++k) {
          fluid_d[(1 + k) * nc + i] =
            mass * ((i == 0 ? 0.3 * (k + 1 - 2 * s) * mean : 0.0) + 0.1 * rnd());
        }
        fluid_d[4 * nc + i] = (i == 0 ? (4.0 + 3.0 * s) * mass * mean : 0.0) + 0.1 * rnd();
      }
      for (int i = 0; i < NFLUID * nc; ++i) {
        fluid_d[i] *= jacob_d[0];
      }
      for (int r = 0; r < num_species; ++r) {
        double *nu_d = gkyl_array_fetch(dat.nu[s][r], loc);
        for (int q = 0; q < nc; ++q) {
          nu_d[q] = 2.0 / (1.0 + 3.0 * s + r) + 0.1 * rnd();
        }
      }
    }
  }
  return dat;
}

static void
cross_data_release(struct cross_data *dat)
{
  gkyl_array_release(dat->jacob);
  for (int s = 0; s < dat->num_species; ++s) {
    gkyl_array_release(dat->fluid[s]);
    gkyl_array_release(dat->prim[s]);
    for (int r = 0; r < dat->num_species; ++r) {
      gkyl_array_release(dat->nu[s][r]);
    }
  }
}

// Node values, exchange coefficient and kinetic rate of a pair, and the floor.
static void
test_rate(int cdim, double heat_fac)
{
  struct cross_data dat = cross_data_new(cdim, 2);
  int nc = dat.num_nodes;
  const double *mass = dat.mass;

  struct gkyl_array *fluid_nodal[2], *coef[2], *nu_kin[2];
  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&dat.cbasis, 0, false);
  struct gkyl_pkpm_coll_cross *up = gkyl_pkpm_coll_cross_new(&dat.cbasis, heat_fac, false);
  for (int s = 0; s < 2; ++s) {
    fluid_nodal[s] = gkyl_array_new(GKYL_DOUBLE, NFLUID * nc, dat.crange.volume);
    coef[s] = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
    nu_kin[s] = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
    gkyl_array_clear(nu_kin[s], 0.0);
    gkyl_dg_conf_nodal_conf_to_nodal(conf_nodal, &dat.crange, dat.fluid[s], fluid_nodal[s]);
    gkyl_pkpm_coll_cross_prim(up, &dat.crange, mass[s], dat.jacob, dat.fluid[s], dat.prim[s]);
  }
  for (int s = 0; s < 2; ++s) {
    int r = 1 - s;
    gkyl_pkpm_coll_cross_coef(
      up, &dat.crange, mass[s], mass[r], dat.nu[s][r], dat.nu[r][s], dat.prim[s], dat.prim[r],
      coef[s]
    );
    gkyl_pkpm_coll_cross_kinetic_rate(
      up, &dat.crange, mass[s], mass[r], dat.nu[s][r], dat.nu[r][s], dat.prim[s], dat.prim[r],
      nu_kin[s]
    );
  }

  double err_prim = 0.0, err_coef = 0.0, err_kin = 0.0, max_coef = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&dat.crange, iter.idx);
    const double *jacob_d = gkyl_array_cfetch(dat.jacob, loc);
    const double *coef0 = gkyl_array_cfetch(coef[0], loc), *coef1 = gkyl_array_cfetch(coef[1], loc);
    for (int s = 0; s < 2; ++s) {
      int r = 1 - s;
      const double *fs = gkyl_array_cfetch(fluid_nodal[s], loc);
      const double *fr = gkyl_array_cfetch(fluid_nodal[r], loc);
      const double *prim_d = gkyl_array_cfetch(dat.prim[s], loc);
      const double *nus = gkyl_array_cfetch(dat.nu[s][r], loc);
      const double *nur = gkyl_array_cfetch(dat.nu[r][s], loc);
      const double *nu_kin_d = gkyl_array_cfetch(nu_kin[s], loc);
      for (int q = 0; q < nc; ++q) {
        // node values of both species from the node values of their fluid states
        double rho_s = fs[q] / jacob_d[0], rho_r = fr[q] / jacob_d[0];
        double us[3], ur[3], ms_sq = 0.0, mr_sq = 0.0, du_sq = 0.0;
        for (int k = 0; k < 3; ++k) {
          us[k] = fs[(1 + k) * nc + q] / fs[q];
          ur[k] = fr[(1 + k) * nc + q] / fr[q];
          ms_sq += fs[(1 + k) * nc + q] * fs[(1 + k) * nc + q];
          mr_sq += fr[(1 + k) * nc + q] * fr[(1 + k) * nc + q];
          du_sq += (ur[k] - us[k]) * (ur[k] - us[k]);
          err_prim = fmax(err_prim, fabs(prim_d[(1 + k) * nc + q] - us[k]));
        }
        double vtsq_s = (2.0 * fs[4 * nc + q] - ms_sq / fs[q]) / (3.0 * fs[q]);
        double vtsq_r = (2.0 * fr[4 * nc + q] - mr_sq / fr[q]) / (3.0 * fr[q]);
        err_prim = fmax(err_prim, fabs(prim_d[q] - rho_s / mass[s]));
        err_prim = fmax(err_prim, fabs(prim_d[4 * nc + q] - vtsq_s));

        // K = 2 A_s A_r / (A_s + A_r), the same for both members of the pair
        double a_s = rho_s * nus[q], a_r = rho_r * nur[q], fac = 2.0 * a_r / (a_s + a_r);
        double coef_d = s == 0 ? coef0[q] : coef1[q];
        err_coef = fmax(err_coef, fabs(coef_d - a_s * fac));
        max_coef = fmax(max_coef, a_s * fac);
        TEST_CHECK(coef0[q] == coef1[q]);

        // LBO: nu_sr vtsq_sr / vtsq_s. BGK: nu_sr plus the rate of the thermal energy,
        // 3 rho_s nu_sr (vtsq_sr - vtsq_s) / ((3/2) rho_s vtsq_s), with the Greene vtsq_sr.
        double vtsq_cross =
          vtsq_s +
          fac * (mass[r] * vtsq_r - mass[s] * vtsq_s + mass[r] * du_sq / 3.0) / (mass[s] + mass[r]);
        double kin_ref = nus[q] * (1.0 + heat_fac * (vtsq_cross - vtsq_s) / vtsq_s);
        err_kin = fmax(err_kin, fabs(nu_kin_d[q] - kin_ref));
      }
    }
  }
  TEST_CHECK(err_prim < 1.0e-13);
  TEST_MSG("node values of the fluid state differ by %g", err_prim);
  TEST_CHECK(err_coef < 1.0e-13 * max_coef);
  TEST_MSG("exchange coefficient differs by %g (largest %g)", err_coef, max_coef);
  TEST_CHECK(err_kin < 1.0e-13);
  TEST_MSG("kinetic rate differs by %g", err_kin);

  // the floor acts on negative rates only: lower the rate at the first node of every cell
  struct gkyl_array *floored = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
  gkyl_array_copy(floored, nu_kin[0]);
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    double *floored_d = gkyl_array_fetch(floored, gkyl_range_idx(&dat.crange, iter.idx));
    floored_d[0] -= 100.0;
  }
  gkyl_pkpm_coll_cross_floor(up, &dat.crange, floored);
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&dat.crange, iter.idx);
    const double *nu_kin_d = gkyl_array_cfetch(nu_kin[0], loc);
    const double *floored_d = gkyl_array_cfetch(floored, loc);
    TEST_CHECK(nu_kin_d[0] < 100.0 && floored_d[0] == 0.0);
    for (int q = 1; q < nc; ++q) {
      TEST_CHECK(nu_kin_d[q] > 0.0 && floored_d[q] == nu_kin_d[q]);
    }
  }
  gkyl_array_release(floored);

  gkyl_pkpm_coll_cross_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  for (int s = 0; s < 2; ++s) {
    gkyl_array_release(fluid_nodal[s]);
    gkyl_array_release(coef[s]);
    gkyl_array_release(nu_kin[s]);
  }
  cross_data_release(&dat);
}

// The backward-Euler step of the exchange over a time step several times the exchange time,
// every species colliding with every other one.
static void
test_exchange(int cdim, int num_species)
{
  struct cross_data dat = cross_data_new(cdim, num_species);
  int nc = dat.num_nodes;
  const double *mass = dat.mass;
  double dt = 3.0;

  struct gkyl_array *fluid_old[NSPEC], *prim_new[NSPEC], *temp_ratio[NSPEC], *shape, *nu_tot;
  struct gkyl_array *coef[NSPEC][NSPEC];
  struct gkyl_pkpm_coll_cross_species species[NSPEC];
  struct gkyl_pkpm_coll_cross *up = gkyl_pkpm_coll_cross_new(&dat.cbasis, 2.0, false);
  for (int s = 0; s < num_species; ++s) {
    fluid_old[s] = gkyl_array_new(GKYL_DOUBLE, NFLUID * nc, dat.crange.volume);
    prim_new[s] = gkyl_array_new(GKYL_DOUBLE, NPRIM * nc, dat.crange.volume);
    temp_ratio[s] = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
    gkyl_array_copy(fluid_old[s], dat.fluid[s]);
    gkyl_pkpm_coll_cross_prim(up, &dat.crange, mass[s], dat.jacob, dat.fluid[s], dat.prim[s]);
  }
  for (int s = 0; s < num_species; ++s) {
    species[s] = (struct gkyl_pkpm_coll_cross_species){
      .mass = mass[s],
      .prim = dat.prim[s],
      .fluid = dat.fluid[s],
      .temp_ratio = temp_ratio[s],
    };
    for (int r = 0; r < num_species; ++r) {
      if (r == s) {
        continue;
      }
      coef[s][r] = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
      gkyl_pkpm_coll_cross_coef(
        up, &dat.crange, mass[s], mass[r], dat.nu[s][r], dat.nu[r][s], dat.prim[s], dat.prim[r],
        coef[s][r]
      );
      species[s].partner[species[s].num_partners] = r;
      species[s].coef[species[s].num_partners++] = coef[s][r];
    }
  }
  gkyl_pkpm_coll_cross_advance(up, &dat.crange, dt, dat.jacob, num_species, species);
  for (int s = 0; s < num_species; ++s) {
    gkyl_pkpm_coll_cross_prim(up, &dat.crange, mass[s], dat.jacob, dat.fluid[s], prim_new[s]);
  }

  double err_cons = 0.0, err_rho = 0.0, err_mom = 0.0, err_temp = 0.0, err_ratio = 0.0;
  double max_fluid = 0.0, max_coef_dt = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&dat.crange, iter.idx);
    // density unchanged; momentum and energy of the species conserved (expansions, so every node)
    for (int i = 0; i < NFLUID * nc; ++i) {
      double total_old = 0.0, total_new = 0.0;
      for (int s = 0; s < num_species; ++s) {
        const double *old_d = gkyl_array_cfetch(fluid_old[s], loc);
        const double *new_d = gkyl_array_cfetch(dat.fluid[s], loc);
        total_old += old_d[i];
        total_new += new_d[i];
        max_fluid = fmax(max_fluid, fabs(old_d[i]));
        if (i < nc) {
          err_rho = fmax(err_rho, fabs(new_d[i] - old_d[i]));
        }
      }
      err_cons = fmax(err_cons, fabs(total_new - total_old));
    }

    // the backward-Euler equations at the nodes
    for (int s = 0; s < num_species; ++s) {
      const double *ps = gkyl_array_cfetch(dat.prim[s], loc),
                   *ns = gkyl_array_cfetch(prim_new[s], loc);
      const double *ratio_d = gkyl_array_cfetch(temp_ratio[s], loc);
      for (int q = 0; q < nc; ++q) {
        double rho_s = mass[s] * ps[q], step_sq = 0.0, heat = 0.0, temp_flux = 0.0;
        double temp_old = mass[s] * ps[4 * nc + q], temp_new = mass[s] * ns[4 * nc + q];
        double res_mom[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; ++k) {
          double step = ns[(1 + k) * nc + q] - ps[(1 + k) * nc + q];
          res_mom[k] = rho_s * step;
          step_sq += step * step;
        }
        for (int r = 0; r < num_species; ++r) {
          if (r == s) {
            continue;
          }
          const double *nr = gkyl_array_cfetch(prim_new[r], loc);
          const double *coef_d = gkyl_array_cfetch(coef[s][r], loc);
          double du_sq = 0.0;
          for (int k = 0; k < 3; ++k) {
            double du = nr[(1 + k) * nc + q] - ns[(1 + k) * nc + q];
            res_mom[k] -= dt * coef_d[q] * du;
            du_sq += du * du;
          }
          double mass_sum = mass[s] + mass[r];
          temp_flux += dt * coef_d[q] * 3.0 * (mass[r] * nr[4 * nc + q] - temp_new) / mass_sum;
          heat += dt * coef_d[q] * mass[r] * du_sq / mass_sum;
          max_coef_dt = fmax(max_coef_dt, dt * coef_d[q] / rho_s);
        }
        for (int k = 0; k < 3; ++k) {
          err_mom = fmax(err_mom, fabs(res_mom[k]));
        }
        double res_temp =
          1.5 * ps[q] * (temp_new - temp_old) - temp_flux - heat - 0.5 * rho_s * step_sq;
        err_temp = fmax(err_temp, fabs(res_temp));
        err_ratio = fmax(err_ratio, fabs(ratio_d[q] - temp_old / temp_new));
        TEST_CHECK(temp_new > 0.0);
      }
    }
  }
  TEST_CHECK(max_coef_dt > 1.0); // a stiff step
  TEST_CHECK(err_rho == 0.0);
  TEST_CHECK(err_cons < 1.0e-13 * max_fluid);
  TEST_MSG("total momentum or energy changed by %g (largest %g)", err_cons, max_fluid);
  TEST_CHECK(err_mom < 1.0e-12 * max_fluid);
  TEST_MSG("flows depart from backward Euler by %g", err_mom);
  TEST_CHECK(err_temp < 1.0e-12 * max_fluid);
  TEST_MSG("temperatures depart from backward Euler by %g", err_temp);
  TEST_CHECK(err_ratio < 1.0e-13);
  TEST_MSG("temperature ratio differs by %g", err_ratio);

  // shape factor of the implicit BGK step: min(1, temp_ratio / (1 + nu dt)); a negative
  // frequency at the first node of every cell reaches the cap
  shape = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
  nu_tot = gkyl_array_new(GKYL_DOUBLE, nc, dat.crange.volume);
  gkyl_array_copy(nu_tot, dat.nu[0][1]);
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    double *nu_tot_d = gkyl_array_fetch(nu_tot, gkyl_range_idx(&dat.crange, iter.idx));
    nu_tot_d[0] = -0.3;
  }
  gkyl_pkpm_coll_cross_shape(up, &dat.crange, dt, nu_tot, temp_ratio[0], shape);
  gkyl_range_iter_init(&iter, &dat.crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&dat.crange, iter.idx);
    const double *nu_tot_d = gkyl_array_cfetch(nu_tot, loc);
    const double *ratio_d = gkyl_array_cfetch(temp_ratio[0], loc);
    const double *shape_d = gkyl_array_cfetch(shape, loc);
    for (int q = 0; q < nc; ++q) {
      double shape_ref = fmin(1.0, ratio_d[q] / (1.0 + nu_tot_d[q] * dt));
      TEST_CHECK(gkyl_compare_double(shape_d[q], shape_ref, 1.0e-14));
      TEST_CHECK(shape_d[q] > 0.0 && shape_d[q] <= 1.0);
    }
  }
  gkyl_array_release(shape);
  gkyl_array_release(nu_tot);

  gkyl_pkpm_coll_cross_release(up);
  for (int s = 0; s < num_species; ++s) {
    gkyl_array_release(fluid_old[s]);
    gkyl_array_release(prim_new[s]);
    gkyl_array_release(temp_ratio[s]);
    for (int r = 0; r < num_species; ++r) {
      if (r != s) {
        gkyl_array_release(coef[s][r]);
      }
    }
  }
  cross_data_release(&dat);
}

static void
test_lbo_1x(void)
{
  test_rate(1, 1.0);
}

static void
test_lbo_2x(void)
{
  test_rate(2, 1.0);
}

static void
test_lbo_3x(void)
{
  test_rate(3, 1.0);
}

static void
test_bgk_1x(void)
{
  test_rate(1, 2.0);
}

static void
test_bgk_2x(void)
{
  test_rate(2, 2.0);
}

static void
test_bgk_3x(void)
{
  test_rate(3, 2.0);
}

static void
test_exchange_pair_1x(void)
{
  test_exchange(1, 2);
}

static void
test_exchange_pair_2x(void)
{
  test_exchange(2, 2);
}

static void
test_exchange_pair_3x(void)
{
  test_exchange(3, 2);
}

static void
test_exchange_three_1x(void)
{
  test_exchange(1, 3);
}

static void
test_exchange_three_3x(void)
{
  test_exchange(3, 3);
}

TEST_LIST = {
  {"lbo_1x", test_lbo_1x},
  {"lbo_2x", test_lbo_2x},
  {"lbo_3x", test_lbo_3x},
  {"bgk_1x", test_bgk_1x},
  {"bgk_2x", test_bgk_2x},
  {"bgk_3x", test_bgk_3x},
  {"exchange_pair_1x", test_exchange_pair_1x},
  {"exchange_pair_2x", test_exchange_pair_2x},
  {"exchange_pair_3x", test_exchange_pair_3x},
  {"exchange_three_1x", test_exchange_three_1x},
  {"exchange_three_3x", test_exchange_three_3x},
  {NULL, NULL}
};
