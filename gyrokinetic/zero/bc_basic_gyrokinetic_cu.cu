/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_basic_gyrokinetic.h>
#include <gkyl_bc_basic_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
}

__global__ static void
gkyl_bc_basic_gyrokinetic_create_set_cu_dev_ptrs(int dir, enum gkyl_edge_loc edge,
  int cdim, enum gkyl_gyrokinetic_bc_type bctype,
  const struct gkyl_basis* basis, int ncomp, struct dg_bc_ctx *ctx,
  struct gkyl_array_copy_func *fout)
{
  ctx->dir = dir;
  ctx->edge = edge;
  ctx->cdim = cdim;
  ctx->basis = basis;
  ctx->ncomp = ncomp;

  switch (bctype) {
    case GKYL_BC_GK_SPECIES_COPY:
      fout->func = copy_bc;
      break;

    case GKYL_BC_GK_SPECIES_FIXED_FUNC:
      assert(basis->poly_order == 1); // MF 2025/11/03: Hardcoded for now.
      assert(basis->b_type == GKYL_BASIS_MODAL_GKHYBRID); // MF 2025/11/03: Hardcoded for now.
      fout->func = phase_boundary_value_bc;
      break;

    case GKYL_BC_GK_SPECIES_ABSORB:
      fout->func = species_absorb_bc;
      break;

    case GKYL_BC_GK_SPECIES_REFLECT:
      fout->func = species_reflect_bc;
      break;

    case GKYL_BC_GK_SPECIES_BOUNDARY_VALUE:
      assert(basis->poly_order == 1); // MF 2025/11/03: Hardcoded for now.
      assert(basis->b_type == GKYL_BASIS_MODAL_GKHYBRID); // MF 2025/11/03: Hardcoded for now.
      fout->func = phase_boundary_value_bc;
      break;

    case GKYL_BC_GK_FIELD_BOUNDARY_VALUE:
      assert(basis->poly_order == 1); // MF 2025/11/03: Hardcoded for now.
      fout->func = conf_boundary_value_bc;
      break;

    default:
      assert(false); // bctype not supported in bc_basic_gyrokinetic.
      break;
  }
  fout->ctx = ctx;
}

struct gkyl_array_copy_func*
gkyl_bc_basic_gyrokinetic_create_arr_copy_func_cu(int dir, enum gkyl_edge_loc edge,
  int cdim, enum gkyl_gyrokinetic_bc_type bctype,
  const struct gkyl_basis *basis, int ncomp)
{
  // Create host context and bc func structs.
  struct dg_bc_ctx *ctx = (struct dg_bc_ctx*) gkyl_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *fout = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  fout->ctx = ctx;

  fout->flags = 0;
  GKYL_SET_CU_ALLOC(fout->flags);

  // Create device context and bc func structs.
  struct dg_bc_ctx *ctx_cu = (struct dg_bc_ctx*) gkyl_cu_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *fout_cu = (struct gkyl_array_copy_func*) gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct dg_bc_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(fout_cu, fout, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  fout->ctx_on_dev = ctx_cu;

  gkyl_bc_basic_gyrokinetic_create_set_cu_dev_ptrs<<<1,1>>>(dir, edge, cdim, bctype, basis, ncomp, ctx_cu, fout_cu);

  // Set parent on_dev pointer.
  fout->on_dev = fout_cu;
  return fout;
}
