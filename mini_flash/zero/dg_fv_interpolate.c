#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_fv_interpolate.h>
#include <gkyl_dg_fv_interpolate_priv.h>

#include <assert.h>
#include <math.h>

struct gkyl_dg_fv_interpolate*
gkyl_dg_fv_interpolate_new(struct gkyl_basis *basis,
  struct gkyl_rect_grid *grid_dg, struct gkyl_rect_grid *grid_fv)
{
  // Allocate memory for updater.
  struct gkyl_dg_fv_interpolate *up = gkyl_malloc(sizeof(*up));

  up->num_basis = basis->num_basis;

  // Verify that grid dimensions match.
  assert(grid_dg->ndim == grid_fv->ndim);
  assert(grid_dg->ndim <= GKYL_MAX_CDIM);
  up->ndim = grid_dg->ndim;

  // Verify that grid resolutions are correct.
  for (int i=0; i<up->ndim; ++i) {
    assert(2*grid_dg->cells[i] == grid_fv->cells[i]);
  }

  // Set function pointers.
  up->kernels.modal_to_fv = choose_dg_fv_interp_modal_to_fv_kernel(up->ndim);
  up->kernels.fv_to_modal = choose_dg_fv_interp_fv_to_modal_kernel(up->ndim);

  return up;
}

// Determine linear offsets of the 2^ndim finite volume cells corresponding to a
// single DG cell. Returns the number of cells.
// Ordering is row-major,  e.g. in 2D, the offsets are:
// (0,0), (0,1), (1,0), (1,1).
static int
fv_cell_offsets(const struct gkyl_range *range_fv, int ndim, long offsets[8])
{
  // Number of finite volume cells.
  int num_fv = (int)pow(2, ndim);

  for (int k=0; k<num_fv; ++k) {
    int rel[GKYL_MAX_DIM] = { 0 };

    for (int d=0; d<ndim; ++d)
      rel[d] = (k / (int)pow(2, ndim-1-d)) % 2;

    offsets[k] = gkyl_range_offset(range_fv, rel);
  }
  return num_fv;
}

// Determine linear index of the finite volume cell in the lower corner of
// the DG cell with index `idx_dg`.
static inline long
fv_cell_lower_idx(const struct gkyl_range *range_fv, int ndim, const int *idx_dg)
{
  int idx_fv[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d)
    idx_fv[d] = 2*idx_dg[d] - 1;
  return gkyl_range_idx(range_fv, idx_fv);
}

void
gkyl_dg_fv_interpolate_fv_to_modal_advance(gkyl_dg_fv_interpolate *up,
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out, int ncomp)
{
  long offsets[8];
  int num_fv = fv_cell_offsets(range_fv, up->ndim, offsets);
  const double *fv_data[8];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range_dg);
  while (gkyl_range_iter_next(&iter)) {
    long lower_fv = fv_cell_lower_idx(range_fv, up->ndim, iter.idx);
    double *dg_out_d = gkyl_array_fetch(dg_out, gkyl_range_idx(range_dg, iter.idx));

    for (int c=0; c<ncomp; ++c) {
      for (int k=0; k<num_fv; ++k)
        fv_data[k] = (const double *)gkyl_array_cfetch(fv_in, lower_fv+offsets[k]) + comp_in + c;
      up->kernels.fv_to_modal(fv_data, &dg_out_d[(comp_out+c)*up->num_basis]);
    }
  }
}

void
gkyl_dg_fv_interpolate_modal_to_fv_advance(gkyl_dg_fv_interpolate *up,
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *dg_in, int comp_in,
  struct gkyl_array *fv_out, int comp_out, int ncomp)
{
  long offsets[8];
  int num_fv = fv_cell_offsets(range_fv, up->ndim, offsets);
  double *fv_data[8];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range_dg);
  while (gkyl_range_iter_next(&iter)) {
    long lower_fv = fv_cell_lower_idx(range_fv, up->ndim, iter.idx);
    const double *dg_in_d = gkyl_array_cfetch(dg_in, gkyl_range_idx(range_dg, iter.idx));

    for (int c=0; c<ncomp; ++c) {
      for (int k=0; k<num_fv; ++k)
        fv_data[k] = (double *)gkyl_array_fetch(fv_out, lower_fv+offsets[k]) + comp_out + c;
      up->kernels.modal_to_fv(&dg_in_d[(comp_in+c)*up->num_basis], fv_data);
    }
  }
}

void
gkyl_dg_fv_interpolate_release(gkyl_dg_fv_interpolate *up)
{
  gkyl_free(up);
}
