#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_elastic.h>
#include <gkyl_array_ops.h>
#include <math.h>
#include <assert.h>

typedef void (*emission_elastic_yield_func_t)(double t, const double *xn, double *fout, void *ctx);

// Primary struct for the updater
struct gkyl_bc_emission_elastic {
  int dir, cdim, vdim;
  enum gkyl_edge_loc edge;
  struct gkyl_array_copy_func *reflect_func;
  struct gkyl_emission_elastic_model *elastic_model;
  bool use_gpu;
};

struct bc_elastic_ctx {
  int dir; // direction for BCs.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
reflection(size_t nc, double *out, const double *inp, void *ctx)
{
  struct bc_elastic_ctx *bc_ctx = (struct bc_elastic_ctx *) ctx;
  int dir = bc_ctx->dir, cdim = bc_ctx->cdim;
  
  bc_ctx->basis->flip_odd_sign(dir, inp, out);
  bc_ctx->basis->flip_odd_sign(dir+cdim, out, out);
}

struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func_cu(int dir, int cdim, const struct gkyl_basis *basis,
  int ncomp);

void
gkyl_bc_emission_elastic_set_extern_params_cu(const struct gkyl_bc_emission_elastic *up,
  int cdim, int vdim, double mass);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.
 *
 * @param up BC updater
 * @param emit_skin_r Range over the species skin cells
 * @param buff_arr BC buffer array
 * @param f_skin Skin cell distribution
 * @param f_emit Emitted distribution
 * @param elastic_yield Projection of elastic yield model onto basis
 * @param basis Pointer to basis functions on host
 */
void
gkyl_bc_emission_elastic_advance_cu(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis);
#endif
