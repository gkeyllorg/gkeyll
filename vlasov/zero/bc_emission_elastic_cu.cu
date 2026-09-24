/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

__global__ static void
gkyl_bc_emission_elastic_set_extern_params_cu_ker(struct gkyl_emission_elastic_model *elastic_model,
  int cdim, int vdim, double mass)
{
  elastic_model->cdim = cdim;
  elastic_model->vdim = vdim;
  elastic_model->mass = mass;
}

__global__ static void
gkyl_bc_emission_elastic_create_set_cu_dev_ptrs(int dir, int cdim, const struct gkyl_basis* basis,
  int ncomp, struct bc_elastic_ctx *ctx, struct gkyl_array_copy_func *fout)
{
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->basis = basis;
  ctx->ncomp = ncomp;

  fout->func = reflection;
  fout->ctx = ctx;
}

void
gkyl_bc_emission_elastic_set_extern_params_cu(const struct gkyl_bc_emission_elastic *up,
  int cdim, int vdim, double mass)
{
  gkyl_bc_emission_elastic_set_extern_params_cu_ker<<<1, 1>>>(up->elastic_model->on_dev,
    cdim, vdim, mass);
}

struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func_cu(int dir, int cdim, const struct gkyl_basis *basis,
  int ncomp)
{
  struct bc_elastic_ctx *ctx = (struct bc_elastic_ctx*) gkyl_malloc(sizeof(struct bc_elastic_ctx));
  struct gkyl_array_copy_func *fout = (struct gkyl_array_copy_func*)
    gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  fout->ctx = ctx;
  
  fout->flags = 0;
  GKYL_SET_CU_ALLOC(fout->flags);
  
  struct bc_elastic_ctx *ctx_cu = (struct bc_elastic_ctx*)
    gkyl_cu_malloc(sizeof(struct bc_elastic_ctx));
  struct gkyl_array_copy_func *fout_cu = (struct gkyl_array_copy_func*)
    gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct bc_elastic_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(fout_cu, fout, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  fout->ctx_on_dev = ctx_cu;

  gkyl_bc_emission_elastic_create_set_cu_dev_ptrs<<<1,1>>>(dir, cdim, basis, ncomp, ctx_cu,
    fout_cu);

  fout->on_dev = fout_cu;
  return fout;
}
