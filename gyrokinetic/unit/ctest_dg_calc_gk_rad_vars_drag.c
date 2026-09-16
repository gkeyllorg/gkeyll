// Test the gyrokinetic radiation drag-coefficient container allocator
// (gkyl_dg_calc_gk_rad_vars_drag_new / _release). Verifies that the returned
// per-collision structs record the correct number of densities, that each
// per-density drag array has the requested number of components and size, and
// that the data arrays are writable host arrays. Release must not crash.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_calc_gk_rad_vars.h>

void
test_rad_drag_alloc()
{
  int num_collisions = 3;
  int num_densities[] = { 1, 2, 4 };
  int ncomp = 6;
  long sz = 12;

  struct gkyl_gk_rad_drag *drag = gkyl_dg_calc_gk_rad_vars_drag_new(
    num_collisions, num_densities, ncomp, sz, false);

  TEST_CHECK( drag != NULL );

  for (int i=0; i<num_collisions; i++) {
    TEST_CHECK( drag[i].num_dens == num_densities[i] );
    TEST_CHECK( drag[i].data != NULL );
    TEST_CHECK( drag[i].on_dev != NULL );

    for (int n=0; n<num_densities[i]; n++) {
      struct gkyl_array *arr = drag[i].data[n].arr;
      TEST_CHECK( arr != NULL );
      TEST_CHECK( arr->ncomp == (unsigned)ncomp );
      TEST_CHECK( arr->size == (size_t)sz );

      // Array should be a usable host array: write and read back.
      gkyl_array_clear(arr, 0.0);
      double *d = gkyl_array_fetch(arr, 0);
      d[0] = 3.25;
      const double *dc = gkyl_array_cfetch(arr, 0);
      TEST_CHECK( dc[0] == 3.25 );
    }
  }

  gkyl_dg_calc_gk_rad_vars_drag_release(drag, num_collisions, false);
}

void
test_rad_drag_alloc_single()
{
  int num_collisions = 1;
  int num_densities[] = { 5 };
  int ncomp = 3;
  long sz = 8;

  struct gkyl_gk_rad_drag *drag = gkyl_dg_calc_gk_rad_vars_drag_new(
    num_collisions, num_densities, ncomp, sz, false);

  TEST_CHECK( drag != NULL );
  TEST_CHECK( drag[0].num_dens == 5 );
  for (int n=0; n<5; n++) {
    TEST_CHECK( drag[0].data[n].arr->ncomp == (unsigned)ncomp );
    TEST_CHECK( drag[0].data[n].arr->size == (size_t)sz );
  }

  gkyl_dg_calc_gk_rad_vars_drag_release(drag, num_collisions, false);
}

TEST_LIST = {
  { "rad_drag_alloc", test_rad_drag_alloc },
  { "rad_drag_alloc_single", test_rad_drag_alloc_single },
  { NULL, NULL },
};
