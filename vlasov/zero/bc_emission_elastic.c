#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <assert.h>
#include <math.h>

struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func(int dir, int cdim, const struct gkyl_basis *basis,
  int ncomp, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_bc_emission_elastic_create_arr_copy_func_cu(dir, cdim, basis, ncomp);
#endif

  struct bc_elastic_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->ncomp = ncomp;

  struct gkyl_array_copy_func *fout = gkyl_malloc(sizeof(*fout));
  fout->func = reflection;
  fout->ctx = ctx;
  fout->ctx_on_dev = fout->ctx;

  fout->flags = 0;
  GKYL_CLEAR_CU_ALLOC(fout->flags);
  fout->on_dev = fout; // CPU function obj points to itself.
  return fout;
}

struct gkyl_bc_emission_elastic*
gkyl_bc_emission_elastic_new(struct gkyl_emission_elastic_model *elastic_model,
  struct gkyl_array *elastic_yield, int dir, enum gkyl_edge_loc edge, int cdim,
  int vdim, double mass, int ncomp, struct gkyl_rect_grid *grid, struct gkyl_range *emit_buff_r,
  int poly_order, const struct gkyl_basis *dev_basis, struct gkyl_basis *basis,
  struct gkyl_array *proj_buffer, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_emission_elastic *up = gkyl_malloc(sizeof(struct gkyl_bc_emission_elastic));

  up->dir = dir;
  up->cdim = cdim;
  up->vdim = vdim;
  up->edge = edge;
  up->use_gpu = use_gpu;

  // Need to pass on_dev basis to create_arr_copy_func, but host copy to proj_on_basis.
  // These are stored separately by the app, so new function takes both dev_basis and basis
  // as arguments.
  up->reflect_func = gkyl_bc_emission_elastic_create_arr_copy_func(dir, cdim, dev_basis, ncomp,
    use_gpu);

  
  up->elastic_model = gkyl_emission_elastic_model_acquire(elastic_model);
  up->elastic_model->cdim = cdim;
  up->elastic_model->vdim = vdim;
  up->elastic_model->mass = mass;

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(grid, basis, poly_order + 1, 1,
      up->elastic_model->function, up->elastic_model);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, proj_buffer);
    
    gkyl_array_copy(elastic_yield, proj_buffer);
  } else {
    gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, elastic_yield);
  }
#else
  gkyl_proj_on_basis_advance(proj, 0.0, emit_buff_r, elastic_yield);
#endif
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
gkyl_bc_emission_elastic_advance(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis)
{
  gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_skin, up->dir+up->cdim, emit_skin_r,
    up->reflect_func->on_dev);
  // Basis is passed directly instead of by pointer for bin op, so advance uses host copy.
  gkyl_dg_mul_op(basis, 0, f_emit, 0, buff_arr, 0, elastic_yield);
}

void gkyl_bc_emission_elastic_release(struct gkyl_bc_emission_elastic *up)
{
  gkyl_emission_elastic_model_release(up->elastic_model);
  gkyl_free(up->reflect_func->ctx);
  gkyl_free(up->reflect_func);
  // Release updater memory.
  gkyl_free(up);
}
