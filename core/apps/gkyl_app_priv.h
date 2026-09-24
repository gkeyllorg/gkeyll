// Private header for common functions between species. Do not include
// in any public facing header!
#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm.h>

#include <stdio.h>
#include <stdlib.h>

// ranges for use in BCs
struct app_skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// allocate double array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// allocate integer array (filled with zeros)
static struct gkyl_array*
mk_int_arr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_INT, nc, size);
  else
    a = gkyl_array_new(GKYL_INT, nc, size);
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range *rng)
{
  return gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct app_skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}
