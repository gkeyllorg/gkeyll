#include <acutest.h>
#include <gkyl_rho_wall_adjust_priv.h>

#include <float.h>
#include <math.h>
#include <stddef.h>

static void
reject_unchanged(double requested, double other, double axis, double sep,
  int family, int step)
{
  double rho = 123.25, psi = -456.75;
  TEST_CHECK(!gkyl_rho_wall_next(requested, other, axis, sep,
    family, step, &rho, &psi));
  TEST_MSG("expected rejection: requested=%g other=%g axis=%g sep=%g family=%d step=%d",
    requested, other, axis, sep, family, step);
  TEST_CHECK(rho == 123.25 && psi == -456.75);
  TEST_MSG("failure modified outputs: rho=%g psi=%g", rho, psi);
}

static void
test_type_ownership(void)
{
  // The enum's complete declared order is checked, including excluded CORE,
  // IWL and unimplemented USN. A new enum requires an explicit policy decision.
  const int families[] = {
    1, 1, 1, 0, 0,  // standalone DN/LSN, USN, CORE
    1, 1, 1, 1, 1, 1,  // six DN sections
    1, 1, 1,  // three LSN sections
    2, 2, 2, 2,  // upper/lower PF halves
    0, 0, 0,  // split cores and IWL
  };
  TEST_CHECK(sizeof families/sizeof families[0] == GKYL_GEOMETRY_TOKAMAK_IWL+1);
  for (size_t i=0; i<sizeof families/sizeof families[0]; ++i) {
    TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) i) == families[i]);
    TEST_MSG("wrong adjustment family for enum %zu", i);
  }
  TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) -1) == 0);
  TEST_CHECK(gkyl_rho_wall_family((enum gkyl_tok_geo_type) 999) == 0);
  reject_unchanged(1.05, 1.0, 0.0, 1.0,
    gkyl_rho_wall_family(GKYL_GEOMETRY_TOKAMAK_CORE), 1);
}

static void
test_sol_direction_and_flux_sign(void)
{
  double rho, psi;
  TEST_ASSERT(gkyl_rho_wall_next(1.051, 1.0, -2.0, 2.0, 1, 1, &rho, &psi));
  TEST_CHECK(fabs(rho-1.05) < 8*DBL_EPSILON);
  TEST_CHECK(fabs(psi-2.41) < 32*DBL_EPSILON);
  double positive = psi;
  TEST_ASSERT(gkyl_rho_wall_next(1.051, 1.0, 2.0, -2.0, 1, 1, &rho, &psi));
  TEST_CHECK(fabs(psi+positive) < 32*DBL_EPSILON);
  TEST_CHECK(psi < -2.0); // still on the SOL side of native separatrix

  // A nonseparatrix fixed boundary also limits shrinkage.
  TEST_ASSERT(gkyl_rho_wall_next(1.25, 1.125, 0.0, 1.0, 1, 124, &rho, &psi));
  TEST_CHECK(rho > 1.125 && rho < 1.25);
  reject_unchanged(1.25, 1.125, 0.0, 1.0, 1, 125); // exact collapse
  reject_unchanged(1.25, 1.125, 0.0, 1.0, 1, 126); // crossing
  reject_unchanged(1.125, 0.9, 0.0, 1.0, 1, 125); // separatrix stop
  reject_unchanged(1.0, 0.9, 0.0, 1.0, 1, 1);
  reject_unchanged(0.9, 0.8, 0.0, 1.0, 1, 1);
}

static void
test_pf_direction_and_flux_sign(void)
{
  double rho, psi;
  TEST_ASSERT(gkyl_rho_wall_next(0.899, 1.0, -2.0, 2.0, 2, 1, &rho, &psi));
  TEST_CHECK(fabs(rho-0.9) < 8*DBL_EPSILON);
  TEST_CHECK(fabs(psi-1.24) < 32*DBL_EPSILON);
  double positive = psi;
  TEST_ASSERT(gkyl_rho_wall_next(0.899, 1.0, 2.0, -2.0, 2, 1, &rho, &psi));
  TEST_CHECK(fabs(psi+positive) < 32*DBL_EPSILON);
  TEST_CHECK(psi > -2.0); // PF remains below rho=1 for either flux sign
  TEST_ASSERT(gkyl_rho_wall_next(0.0, 1.0, 0.0, 1.0, 2, 1, &rho, &psi));
  TEST_CHECK(rho == 0.001);

  TEST_ASSERT(gkyl_rho_wall_next(0.75, 0.875, 0.0, 1.0, 2, 124, &rho, &psi));
  TEST_CHECK(rho < 0.875 && rho > 0.75);
  reject_unchanged(0.75, 0.875, 0.0, 1.0, 2, 125);
  reject_unchanged(0.75, 0.875, 0.0, 1.0, 2, 126);
  reject_unchanged(0.875, 1.2, 0.0, 1.0, 2, 125); // separatrix stop
  reject_unchanged(-0.1, 1.0, 0.0, 1.0, 2, 1);
  reject_unchanged(1.0, 1.1, 0.0, 1.0, 2, 1);
  reject_unchanged(1.05, 1.1, 0.0, 1.0, 2, 1);
}

static void
test_original_request_steps(void)
{
  // Preserve the non-millistep part of the request and avoid accumulated
  // subtraction error. Query order is deliberately nonmonotone.
  const double requested = 3.14159265;
  const int steps[] = { 2000, 1, 1700, 17, 2000 };
  double rho, psi, first = 0.0;
  for (size_t i=0; i<sizeof steps/sizeof steps[0]; ++i) {
    TEST_ASSERT(gkyl_rho_wall_next(requested, 1.0, 0.0, 1.0,
      1, steps[i], &rho, &psi));
    TEST_CHECK(fabs(rho-(requested-steps[i]/1000.0)) < 4*DBL_EPSILON);
    TEST_CHECK(fabs(sqrt(psi)-rho) < 8*DBL_EPSILON);
    if (i == 0) first = rho;
    if (i == 4) TEST_CHECK(rho == first);
  }
  TEST_ASSERT(gkyl_rho_wall_next(0.10025, 1.0, 0.0, 1.0,
    2, 899, &rho, &psi));
  TEST_CHECK(fabs(rho-0.99925) < 4*DBL_EPSILON);

  TEST_ASSERT(gkyl_rho_wall_next(12.0, 1.0, 0.0, 1.0,
    1, GKYL_RHO_WALL_MAX_STEPS, &rho, &psi));
  TEST_CHECK(rho == 2.0 && psi == 4.0);
  reject_unchanged(12.0, 1.0, 0.0, 1.0, 1, GKYL_RHO_WALL_MAX_STEPS+1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 1, 0);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 1, -1);
}

static void
test_invalid_and_nonfinite_inputs(void)
{
  const double bad[] = { NAN, INFINITY, -INFINITY };
  for (size_t i=0; i<sizeof bad/sizeof bad[0]; ++i) {
    TEST_CHECK(!gkyl_rho_wall_finite(bad[i]));
    reject_unchanged(bad[i], 1.0, 0.0, 1.0, 1, 1);
    reject_unchanged(1.05, bad[i], 0.0, 1.0, 1, 1);
    reject_unchanged(1.05, 1.0, bad[i], 1.0, 1, 1);
    reject_unchanged(1.05, 1.0, 0.0, bad[i], 1, 1);
    reject_unchanged(bad[i], 1.0, 0.0, 1.0, 2, 1);
  }
  TEST_CHECK(gkyl_rho_wall_finite(DBL_MAX));
  TEST_CHECK(gkyl_rho_wall_finite(-DBL_MAX));
  TEST_CHECK(gkyl_rho_wall_finite(0.0));
  reject_unchanged(1.05, 1.0, 2.0, 2.0, 1, 1);
  reject_unchanged(1.05, -0.5, 0.0, 1.0, 1, 1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, 3, 1);
  reject_unchanged(1.05, 1.0, 0.0, 1.0, -1, 1);
  reject_unchanged(2.0, 1.0, 0.0, DBL_MAX, 1, 1); // output flux overflow
  reject_unchanged(1.05, 1.0, DBL_MAX, -DBL_MAX, 1, 1); // span overflow
  reject_unchanged(DBL_MAX, 1.0, 0.0, 1.0, 1, 1); // no representable progress

  double value = 123.25;
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, NULL, &value));
  TEST_CHECK(value == 123.25);
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, &value, NULL));
  TEST_CHECK(value == 123.25);
  TEST_CHECK(!gkyl_rho_wall_next(1.05, 1.0, 0.0, 1.0, 1, 1, &value, &value));
  TEST_CHECK(value == 123.25);
}

TEST_LIST = {
  { "type_ownership", test_type_ownership },
  { "sol_direction_and_flux_sign", test_sol_direction_and_flux_sign },
  { "pf_direction_and_flux_sign", test_pf_direction_and_flux_sign },
  { "original_request_steps", test_original_request_steps },
  { "invalid_and_nonfinite_inputs", test_invalid_and_nonfinite_inputs },
  { NULL, NULL },
};
