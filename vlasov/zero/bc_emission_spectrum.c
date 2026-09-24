#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

// Increment an int vector by fact*del[d] in each direction d.
static inline void
incr_int_array(int ndim, int fact, const int * GKYL_RESTRICT del,
  const int * GKYL_RESTRICT inp, int *GKYL_RESTRICT out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

// KB - Number of cells in the positive and negative directions
// is assumed to be the same symmetrically, half of the cells in that dimension. Need to
// do this more rigorously.
void
gkyl_bc_emission_flux_ranges(struct gkyl_range *flux_r, int dir,
  const struct gkyl_range *parent, const int *nghost, enum gkyl_edge_loc edge)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};
  
  if (edge == GKYL_LOWER_EDGE) {
    incr_int_array(ndim, 0, nghost, parent->lower, lo);
    incr_int_array(ndim, 0, nghost, parent->upper, up);

    int nneg = (up[dir] - lo[dir] + 1)/2;
  
    up[dir] = lo[dir] + nneg;
    gkyl_sub_range_init(flux_r, parent, lo, up);
  } else {
    incr_int_array(ndim, 0, nghost, parent->lower, lo);
    incr_int_array(ndim, 0, nghost, parent->upper, up);

    int npos = (up[dir] - lo[dir] + 1)/2;
    
    lo[dir] = up[dir] - npos;
    gkyl_sub_range_init(flux_r, parent, lo, up);
  }
}

void
gkyl_bc_emission_spectrum_sey_calc(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *impact_buff_r)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_emission_spectrum_sey_calc_cu(up, yield, grid, impact_buff_r);
    return;
  }
#endif
  double xc[GKYL_MAX_DIM];
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, impact_buff_r);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(impact_buff_r, iter.idx);
    double *out = gkyl_array_fetch(yield, loc);
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);
    up->yield_model->function(out, up->yield_model, xc);
  }
}

struct gkyl_bc_emission_spectrum*
gkyl_bc_emission_spectrum_new(struct gkyl_emission_spectrum_model *spectrum_model,
  struct gkyl_emission_yield_model *yield_model, struct gkyl_array *yield,
  struct gkyl_array *spectrum, int dir, enum gkyl_edge_loc edge, int cdim, int vdim,
  double mass_in, double mass_out, struct gkyl_range *impact_buff_r, struct gkyl_range *emit_buff_r,
  struct gkyl_rect_grid *impact_grid, struct gkyl_rect_grid *emit_grid, int poly_order,
  struct gkyl_basis *basis, struct gkyl_array *proj_buffer, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission_spectrum *up = gkyl_malloc(sizeof(struct gkyl_bc_emission_spectrum));

  up->dir = dir;
  up->cdim = cdim;
  up->vdim = vdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->grid = impact_grid;

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    ghost[cdim+d] = 0;
  }

  up->spectrum_model = gkyl_emission_spectrum_model_acquire(spectrum_model);
  up->spectrum_model->cdim = cdim;
  up->spectrum_model->vdim = vdim;
  up->spectrum_model->mass = mass_out;
  up->yield_model = gkyl_emission_yield_model_acquire(yield_model);
  up->yield_model->cdim = cdim;
  up->yield_model->vdim = vdim;
  up->yield_model->mass = mass_in;

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(emit_grid, basis, poly_order + 1, 1,
      up->spectrum_model->distribution, up->spectrum_model);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_bc_emission_spectrum_set_extern_params_cu(up, cdim, vdim, mass_in, mass_out);
    gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, proj_buffer);
    
    gkyl_array_copy(spectrum, proj_buffer);
  } else {
    gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, spectrum);
  }
#else
  gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, spectrum);
#endif
  gkyl_bc_emission_spectrum_sey_calc(up, yield, impact_grid, impact_buff_r);
  gkyl_proj_on_basis_release(proj);
  
  return up;
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

void
gkyl_bc_emission_spectrum_advance(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_bc_emission_spectrum_advance_cu(up, impact_buff_r, impact_cbuff_r, emit_buff_r,
      bflux, f_emit, yield, spectrum, weight, flux, k);
  }
#endif
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<impact_cbuff_r->ndim; ++d) rem_dir[d] = 1;

  struct gkyl_range vel_buff_r;
  struct gkyl_range_iter conf_iter, vel_iter;

  gkyl_array_clear_range(weight, 0.0, impact_cbuff_r);

  gkyl_range_iter_init(&conf_iter, impact_cbuff_r);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(impact_cbuff_r, conf_iter.idx);

    gkyl_range_deflate(&vel_buff_r, impact_buff_r, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_buff_r);

    double *w = gkyl_array_fetch(weight, midx);
    
    while (gkyl_range_iter_next(&vel_iter)) {
      copy_idx_arrays(impact_cbuff_r->ndim, impact_buff_r->ndim, conf_iter.idx, vel_iter.idx, pidx);
      gkyl_rect_grid_cell_center(up->grid, pidx, xc);
      
      long loc = gkyl_range_idx(&vel_buff_r, vel_iter.idx);
      
      const double *inp = gkyl_array_cfetch(bflux, loc);
      const double *gain = gkyl_array_cfetch(yield, loc);
      
      bc_weighted_delta(inp, up->cdim, up->dir, up->edge, xc, gain, w);
    }
    double effective_delta = w[0]/w[1];
    const double *boundary_flux = gkyl_array_cfetch(flux, midx);
    double *out = gkyl_array_fetch(k, midx);
    
    up->spectrum_model->normalization(out, up->spectrum_model, boundary_flux, effective_delta);
    gkyl_array_accumulate(f_emit, out[0], spectrum);
  }
}

void gkyl_bc_emission_spectrum_release(struct gkyl_bc_emission_spectrum *up)
{
  gkyl_emission_spectrum_model_release(up->spectrum_model);
  gkyl_emission_yield_model_release(up->yield_model);
  // Release updater memory.
  gkyl_free(up);
}
