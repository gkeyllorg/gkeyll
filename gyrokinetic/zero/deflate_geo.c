#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_deflate_geo_priv.h>

struct gkyl_deflate_geo *
gkyl_deflate_geo_new(
  const struct gkyl_basis *cbasis, const struct gkyl_basis *deflated_cbasis,
  const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid,
  const int *rem_dirs, bool use_gpu
)
{
  gkyl_deflate_geo *up = gkyl_malloc(sizeof(gkyl_deflate_geo));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;
  up->deflated_basis = deflated_cbasis;
  up->rem_dirs = gkyl_malloc(3 * sizeof(int));
  for (int i = 0; i < 3; i++) {
    up->rem_dirs[i] = rem_dirs[i]; // 1 to remove
  }
  up->kernel = deflate_geo_choose_kernel(
    up->rem_dirs, up->deflated_grid->ndim, cbasis->b_type, cbasis->poly_order
  );

  return up;
}

struct gkyl_deflate_geo_surf *
gkyl_deflate_geo_surf_new(
  const struct gkyl_basis *cbasis, int deflated_num_basis, const struct gkyl_rect_grid *grid,
  const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, int dir, bool use_gpu
)
{
  struct gkyl_deflate_geo_surf *up = gkyl_malloc(sizeof(struct gkyl_deflate_geo_surf));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;
  up->deflated_num_basis = deflated_num_basis;
  up->rem_dirs = gkyl_malloc(3 * sizeof(int));
  up->dir = dir;
  for (int i = 0; i < 3; i++) {
    up->rem_dirs[i] = rem_dirs[i]; // 1 to remove
  }
  up->kernel = deflate_geo_surf_choose_kernel(
    up->dir, up->deflated_grid->ndim, cbasis->b_type, cbasis->poly_order
  );

  return up;
}

void
gkyl_deflate_geo_advance(
  const gkyl_deflate_geo *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
)
{
  // Inflated grid will always be 1 long in other directions
  // So deflate reange to the cell at index 1 in the ignored directions (z or x and y). 1 and not 0 because it is a local range

  int loc_dir[3] = {0.};

  for (int i = 0; i < 3; i++) {
    loc_dir[i] = 1;
    if (up->rem_dirs[i] == 1) {
      loc_dir[i] = 1;
    }
  }

  int do_idx[3];
  for (int i = 0; i < 3; i++) {
    do_idx[i] = loc_dir[i];
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_range);

  while (gkyl_range_iter_next(&iter)) {
    int count = 0;
    for (int i = 0; i < 3; i++) {
      if (up->rem_dirs[i] == 0) {
        do_idx[i] = iter.idx[count];
        count += 1;
      }
    }

    long loc = gkyl_range_idx(range, do_idx);
    const double *fld = gkyl_array_cfetch(field, loc);

    long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
    double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
    for (int c = 0; c < ncomp; c++) {
      up->kernel(&fld[c * up->basis->num_basis], &fld_deflated[c * up->deflated_basis->num_basis]);
    }
  }
}

// Evaluate at the center of the removed quadrature coordinates. With p=1
// this is the average of the two Gauss nodes; p=2 has a node at the center.
static void
deflate_geo_nodal_center(
  const int *rem_dirs, const struct gkyl_range *nrange, const struct gkyl_range *deflated_nrange,
  const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp
)
{
  int samples = 1;
  for (int d = 0; d < 3; ++d) {
    if (rem_dirs[d]) {
      int count = gkyl_range_shape(nrange, d);
      assert(count == 2 || count == 3);
      if (count == 2) {
        samples *= 2;
      }
    }
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_nrange);
  while (gkyl_range_iter_next(&iter)) {
    double *out = gkyl_array_fetch(deflated_field, gkyl_range_idx(deflated_nrange, iter.idx));
    for (int c = 0; c < ncomp; ++c) {
      out[c] = 0.0;
    }
    for (int sample = 0; sample < samples; ++sample) {
      int idx[3], bit = 0, kept = 0;
      for (int d = 0; d < 3; ++d) {
        if (rem_dirs[d]) {
          idx[d] =
            nrange->lower[d] + (gkyl_range_shape(nrange, d) == 3 ? 1 : ((sample >> bit++) & 1));
        } else {
          idx[d] = iter.idx[kept++];
        }
      }
      const double *in = gkyl_array_cfetch(field, gkyl_range_idx(nrange, idx));
      for (int c = 0; c < ncomp; ++c) {
        out[c] += in[c] / samples;
      }
    }
  }
}

void
gkyl_deflate_geo_advance_nodal(
  const gkyl_deflate_geo *up, const struct gkyl_range *nrange,
  const struct gkyl_range *deflated_nrange, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
)
{
  deflate_geo_nodal_center(up->rem_dirs, nrange, deflated_nrange, field, deflated_field, ncomp);
}

void
gkyl_deflate_geo_surf_advance(
  const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *range,
  const struct gkyl_range *deflated_range, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
)
{
  // Inflated grid will always be 1 long in other directions
  // So deflate reange to the cell at index 1 in the ignored directions (z or x and y). 1 and not 0 because it is a local range

  int loc_dir[3] = {0.};

  for (int i = 0; i < 3; i++) {
    loc_dir[i] = 1;
    if (up->rem_dirs[i] == 1) {
      loc_dir[i] = 1;
    }
  }

  int do_idx[3];
  for (int i = 0; i < 3; i++) {
    do_idx[i] = loc_dir[i];
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, deflated_range);

  while (gkyl_range_iter_next(&iter)) {
    int count = 0;
    for (int i = 0; i < 3; i++) {
      if (up->rem_dirs[i] == 0) {
        do_idx[i] = iter.idx[count];
        count += 1;
      }
    }

    long loc = gkyl_range_idx(range, do_idx);
    const double *fld = gkyl_array_cfetch(field, loc);

    long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
    double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
    for (int c = 0; c < ncomp; c++) {
      up->kernel(&fld[c * up->basis->num_basis], &fld_deflated[c * up->deflated_num_basis]);
    }
  }
}

void
gkyl_deflate_geo_surf_advance_nodal(
  const struct gkyl_deflate_geo_surf *up, const struct gkyl_range *nrange,
  const struct gkyl_range *deflated_nrange, const struct gkyl_array *field,
  struct gkyl_array *deflated_field, int ncomp
)
{
  deflate_geo_nodal_center(up->rem_dirs, nrange, deflated_nrange, field, deflated_field, ncomp);
}

void
gkyl_deflate_geo_release(gkyl_deflate_geo *up)
{
  gkyl_free(up->rem_dirs);
  gkyl_free(up);
}

void
gkyl_deflate_geo_surf_release(struct gkyl_deflate_geo_surf *up)
{
  gkyl_free(up->rem_dirs);
  gkyl_free(up);
}
