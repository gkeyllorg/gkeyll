#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_gpython.h>
#include <gkyl_util.h>

static void
check_quad_roundtrip(gpython_basis *basis, int expected_num_quad)
{
  TEST_ASSERT(basis != NULL);
  int num_basis = gpython_basis_num_basis(basis);
  int num_quad = gpython_basis_num_quad(basis);
  TEST_ASSERT(num_quad == expected_num_quad);

  double *modal = gkyl_calloc(num_basis, sizeof(double));
  double *quad = gkyl_malloc(num_quad * sizeof(double));
  double *restored = gkyl_malloc(num_basis * sizeof(double));

  // Every mode, including the highest velocity modes, must survive projection.
  for (int mode = 0; mode < num_basis; ++mode) {
    modal[mode] = 1.0;
    for (int node = 0; node < num_quad; ++node) {
      quad[node] = NAN;
    }
    for (int index = 0; index < num_basis; ++index) {
      restored[index] = NAN;
    }
    TEST_CHECK(gpython_basis_modal_to_quad(basis, modal, quad) == 0);
    TEST_CHECK(gpython_basis_quad_to_modal(basis, quad, restored) == 0);
    for (int index = 0; index < num_basis; ++index) {
      TEST_CHECK(gkyl_compare_double(restored[index], modal[index], 1e-12));
    }
    modal[mode] = 0.0;
  }

  gkyl_free(restored);
  gkyl_free(quad);
  gkyl_free(modal);
  gpython_basis_release(basis);
}

void
test_gpython_quad_cartesian()
{
  const char *types[] = {"serendipity", "tensor"};
  TEST_CHECK(gpython_api_version() == GPYTHON_API_VERSION);
  for (int type = 0; type < 2; ++type) {
    for (int order = 1; order <= 2; ++order) {
      int num_quad = 1;
      for (int ndim = 1; ndim <= 3; ++ndim) {
        num_quad *= order + 1;
        check_quad_roundtrip(gpython_basis_new(types[type], ndim, order), num_quad);
      }
    }
  }
}

void
test_gpython_quad_hybrid()
{
  for (int cdim = 1; cdim <= 3; ++cdim) {
    int num_quad = 1 << cdim;
    for (int vdim = 1; vdim <= 3; ++vdim) {
      num_quad *= 3;
      check_quad_roundtrip(gpython_basis_new_hybrid("hybrid", cdim, vdim), num_quad);
    }
    for (int vdim = 1; vdim <= 2; ++vdim) {
      check_quad_roundtrip(
        gpython_basis_new_hybrid("gkhybrid", cdim, vdim), 3 * (1 << (cdim + vdim - 1))
      );
    }
  }
}

void
test_gpython_quad_unavailable()
{
  const char *types[] = {"serendipity", "tensor"};
  for (int type = 0; type < 2; ++type) {
    for (int order = 0; order <= 3; order += 3) {
      gpython_basis *basis = gpython_basis_new(types[type], 1, order);
      TEST_ASSERT(basis != NULL);
      double modal[4] = {0.0}, quad[4] = {0.0};
      TEST_CHECK(gpython_basis_num_quad(basis) == 0);
      TEST_CHECK(gpython_basis_modal_to_quad(basis, modal, quad) == -1);
      TEST_CHECK(gpython_basis_quad_to_modal(basis, quad, modal) == -1);
      gpython_basis_release(basis);
    }
  }
}

TEST_LIST = {
  {"quad_cartesian", test_gpython_quad_cartesian},
  {"quad_hybrid", test_gpython_quad_hybrid},
  {"quad_unavailable", test_gpython_quad_unavailable},
  {NULL, NULL}
};
