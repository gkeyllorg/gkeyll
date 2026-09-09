#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>

#include <assert.h>

// allocate array (filled with zeros)
static inline struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
bc_twistshift_refine_enabled(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo)
{
  gkyl_dg_interpolate_advance(up->refine, fdo, up->ffine);
}

static void
bc_twistshift_refine_disabled(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo)
{
  gkyl_array_copy_range_to_range(up->ffine, fdo, &up->ghost_r, &up->coarse_ghost_r);
}

static void
bc_twistshift_coarsen_enabled(struct gkyl_bc_twistshift *up, struct gkyl_array *ftar)
{
  gkyl_dg_interpolate_advance(up->coarsen, up->ffine, ftar);
}

static void
bc_twistshift_coarsen_disabled(struct gkyl_bc_twistshift *up, struct gkyl_array *ftar)
{
  gkyl_array_copy_range_to_range(ftar, up->ffine, &up->coarse_ghost_r, &up->ghost_r);
}

static void
bc_twistshift_advance_ts(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo,
  struct gkyl_array *ftar)
{
  gkyl_twistshift_dg_advance(up->ts, fdo, ftar);
}

static void
bc_twistshift_advance_ts_filtered(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo,
  struct gkyl_array *ftar)
{
  up->refine_func(up, fdo);
  gkyl_twistshift_dg_advance(up->ts, up->ffine, up->ffine);
  gkyl_dg_lowpass_filter_advance(up->filter, up->ffine, up->filt_buff);
  gkyl_array_copy_range(up->ffine, up->filt_buff, &up->ghost_r);
  up->coarsen_func(up, ftar);
}

static void
bc_twistshift_refine_shift(const struct gkyl_bc_twistshift_inp *inp,
  const struct gkyl_rect_grid *ts_grid, const struct gkyl_range *shear_r_fine,
  struct gkyl_array *shift_dg_fine)
{
  // Refine the DG shift to match the supersampled shear grid.
  int shear_dir = inp->shear_dir;
  int shift_poly_order = inp->shift_poly_order? inp->shift_poly_order : inp->basis->poly_order;

  struct gkyl_basis shift_b;
  gkyl_cart_modal_serendip(&shift_b, 1, shift_poly_order);
  assert(inp->shift_dg->ncomp == shift_b.num_basis);
  assert(shift_dg_fine->ncomp == shift_b.num_basis);
  assert(shift_dg_fine->size == shear_r_fine->volume);

  struct gkyl_rect_grid grid_do, grid_tar;
  gkyl_rect_grid_init(&grid_do, 1, &inp->grid->lower[shear_dir],
    &inp->grid->upper[shear_dir], &inp->grid->cells[shear_dir]);
  gkyl_rect_grid_init(&grid_tar, 1, &ts_grid->lower[shear_dir],
    &ts_grid->upper[shear_dir], &ts_grid->cells[shear_dir]);

  struct gkyl_range shear_r_do;
  gkyl_range_init(&shear_r_do, 1, (int[]) {inp->bcdir_ext_update_r->lower[shear_dir]},
                                 (int[]) {inp->bcdir_ext_update_r->upper[shear_dir]});
  assert(inp->shift_dg->size == shear_r_do.volume);

  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(1, &shift_b,
    &grid_do, &grid_tar, &shear_r_do, shear_r_fine, (int[]) {0}, false);
  gkyl_dg_interpolate_advance(interp, inp->shift_dg, shift_dg_fine);
  gkyl_dg_interpolate_release(interp);
}

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_inew(const struct gkyl_bc_twistshift_inp *inp)
{
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = inp->use_gpu;
  up->filter_half_width = inp->filter_half_width;
  up->filter_cutoff_wavelength = inp->filter_cutoff_wavelength;
  up->upsample_factor = inp->upsample_factor > 1 ? inp->upsample_factor : 1;

  up->filter = NULL;
  up->filt_buff = NULL;
  up->ffine = NULL;
  up->shift_dg_fine = NULL;
  up->refine = NULL;
  up->coarsen = NULL;
  up->refine_func = bc_twistshift_refine_disabled;
  up->coarsen_func = bc_twistshift_coarsen_disabled;

  // The half-width counts cells of the original grid.
  up->half_width_fine = up->filter_half_width * up->upsample_factor;

  // A stencil one fine cell wide is the identity kernel.
  assert(up->half_width_fine != 1);
  // Supersampling is only useful if there is filtering.
  if (up->upsample_factor > 1)
    assert(up->filter_half_width > 0 && up->filter_cutoff_wavelength > 0.0);

  if (up->filter_half_width == 0) {
    // Plain twist-shift.
    struct gkyl_twistshift_dg_inp tsinp = {
      .bc_dir = inp->bc_dir,
      .shift_dir = inp->shift_dir,
      .shear_dir = inp->shear_dir,
      .edge = inp->edge,
      .cdim = inp->cdim,
      .bcdir_ext_update_r = inp->bcdir_ext_update_r,
      .num_ghost = inp->num_ghost,
      .basis = inp->basis,
      .grid = inp->grid,
      .shift_func = inp->shift_func,
      .shift_func_ctx = inp->shift_func_ctx,
      .shift_dg = inp->shift_dg,
      .use_gpu = inp->use_gpu,
      .shift_poly_order = inp->shift_poly_order,
    };
    up->ts = gkyl_twistshift_dg_inew(&tsinp);
    up->advance_func = bc_twistshift_advance_ts;
    return up;
  }

  // Upsampling and filtering attributes.
  const int ndim = inp->bcdir_ext_update_r->ndim;
  // Ghost plane this BC fills, on the field's own grid.
  if (inp->edge == GKYL_LOWER_EDGE)
    gkyl_range_shorten_from_above(&up->coarse_ghost_r, inp->bcdir_ext_update_r,
      inp->bc_dir, inp->num_ghost[inp->bc_dir]);
  else
    gkyl_range_shorten_from_below(&up->coarse_ghost_r, inp->bcdir_ext_update_r,
      inp->bc_dir, inp->num_ghost[inp->bc_dir]);

  // Grid supersampled along shear_dir.
  int fine_cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) fine_cells[d] = inp->grid->cells[d];
  fine_cells[inp->shear_dir] *= up->upsample_factor;
  gkyl_rect_grid_init(&up->ts_grid, ndim, inp->grid->lower, inp->grid->upper, fine_cells);
  // Range for supersampled cells.
  int flo[GKYL_MAX_DIM], fup[GKYL_MAX_DIM];
  for (int d=0; d<ndim; d++) {
    flo[d] = up->coarse_ghost_r.lower[d];
    fup[d] = up->coarse_ghost_r.upper[d];
  }
  flo[inp->shear_dir] = (flo[inp->shear_dir]-1)*up->upsample_factor + 1;
  fup[inp->shear_dir] = fup[inp->shear_dir]*up->upsample_factor;
  gkyl_range_init(&up->ts_ext_r, ndim, flo, fup);
  gkyl_sub_range_init(&up->ts_update_r, &up->ts_ext_r, flo, fup);

  up->ffine = mkarr(inp->use_gpu, inp->basis->num_basis, up->ts_ext_r.volume);
  up->filt_buff = mkarr(inp->use_gpu, inp->basis->num_basis, up->ts_ext_r.volume);

  // Ghost plane on the supersampled grid.
  if (inp->edge == GKYL_LOWER_EDGE)
    gkyl_range_shorten_from_above(&up->ghost_r, &up->ts_update_r, inp->bc_dir,
      inp->num_ghost[inp->bc_dir]);
  else
    gkyl_range_shorten_from_below(&up->ghost_r, &up->ts_update_r, inp->bc_dir,
      inp->num_ghost[inp->bc_dir]);

  up->filter = gkyl_dg_lowpass_filter_new(inp->shear_dir, up->half_width_fine,
    up->filter_cutoff_wavelength, inp->basis, &up->ts_grid, &up->ghost_r, inp->use_gpu);

  if (up->upsample_factor > 1) {
    up->refine = gkyl_dg_interpolate_new(inp->cdim, inp->basis, inp->grid, &up->ts_grid,
      &up->coarse_ghost_r, &up->ghost_r, inp->num_ghost, inp->use_gpu);
    up->coarsen = gkyl_dg_interpolate_new(inp->cdim, inp->basis, &up->ts_grid, inp->grid,
      &up->ghost_r, &up->coarse_ghost_r, inp->num_ghost, inp->use_gpu);
    up->refine_func = bc_twistshift_refine_enabled;
    up->coarsen_func = bc_twistshift_coarsen_enabled;
  }

  // Need to upsample the DG shift too.
  struct gkyl_array *shift_dg = inp->shift_dg;
  if (inp->shift_dg && up->upsample_factor > 1) {
    struct gkyl_range shear_r_fine;
    gkyl_range_init(&shear_r_fine, 1, (int[]) {flo[inp->shear_dir]},
                                      (int[]) {fup[inp->shear_dir]});
    up->shift_dg_fine = gkyl_array_new(GKYL_DOUBLE, inp->shift_dg->ncomp, shear_r_fine.volume);
    bc_twistshift_refine_shift(inp, &up->ts_grid, &shear_r_fine, up->shift_dg_fine);
    shift_dg = up->shift_dg_fine;
  }

  struct gkyl_twistshift_dg_inp tsinp = {
    .bc_dir = inp->bc_dir,
    .shift_dir = inp->shift_dir,
    .shear_dir = inp->shear_dir,
    .edge = inp->edge,
    .cdim = inp->cdim,
    .bcdir_ext_update_r = &up->ts_update_r,
    .num_ghost = inp->num_ghost,
    .basis = inp->basis,
    .grid = &up->ts_grid,
    .shift_func = inp->shift_func,
    .shift_func_ctx = inp->shift_func_ctx,
    .shift_dg = shift_dg,
    .use_gpu = inp->use_gpu,
    .shift_poly_order = inp->shift_poly_order,
  };
  up->ts = gkyl_twistshift_dg_inew(&tsinp);
  up->advance_func = bc_twistshift_advance_ts_filtered;

  return up;
}

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(int bc_dir, int shift_dir, int shear_dir,
  enum gkyl_edge_loc edge, int cdim, const struct gkyl_range *bcdir_ext_update_r, const int *num_ghost,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, evalf_t shift_func, void *shift_func_ctx,
  struct gkyl_array *shift_dg, int shift_poly_order, int filter_half_width,
  double filter_cutoff_wavelength, int upsample_factor, bool use_gpu)
{
  struct gkyl_bc_twistshift_inp inp = {
    .bc_dir                   = bc_dir                  ,
    .shift_dir                = shift_dir               ,
    .shear_dir                = shear_dir               ,
    .edge                     = edge                    ,
    .cdim                     = cdim                    ,
    .bcdir_ext_update_r       = bcdir_ext_update_r      ,
    .num_ghost                = num_ghost               ,
    .basis                    = basis                   ,
    .grid                     = grid                    ,
    .shift_func               = shift_func              ,
    .shift_func_ctx           = shift_func_ctx          ,
    .shift_dg                 = shift_dg                ,
    .use_gpu                  = use_gpu                 ,
    .shift_poly_order         = shift_poly_order        ,
    .filter_half_width        = filter_half_width       ,
    .filter_cutoff_wavelength = filter_cutoff_wavelength,
    .upsample_factor          = upsample_factor         ,
  };
  return gkyl_bc_twistshift_inew(&inp);
}

void
gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
{
  up->advance_func(up, fdo, ftar);
}

void
gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up)
{
  gkyl_twistshift_dg_release(up->ts);
  if (up->filter) {
    gkyl_dg_lowpass_filter_release(up->filter);
    gkyl_array_release(up->ffine);
    gkyl_array_release(up->filt_buff);
  }
  if (up->shift_dg_fine)
    gkyl_array_release(up->shift_dg_fine);
  if (up->refine) {
    gkyl_dg_interpolate_release(up->refine);
    gkyl_dg_interpolate_release(up->coarsen);
  }
  gkyl_free(up);
}
