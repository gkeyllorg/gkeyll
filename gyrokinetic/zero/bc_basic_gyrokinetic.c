#include <gkyl_bc_basic_gyrokinetic.h>
#include <gkyl_bc_basic_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

// Private function to create a pointer to the function that applies the BC,
// i.e., the array_copy_func applied to expansion coefficients in ghost cell.
struct gkyl_array_copy_func*
gkyl_bc_basic_gyrokinetic_create_arr_copy_func(int dir, enum gkyl_edge_loc edge,
  int cdim, enum gkyl_gyrokinetic_bc_type bctype,
  const struct gkyl_basis *basis, int ncomp, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_bc_basic_gyrokinetic_create_arr_copy_func_cu(dir, edge, cdim, bctype, basis, ncomp);
#endif

  struct dg_bc_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->edge = edge;
  ctx->cdim = cdim;
  ctx->ncomp = ncomp;

  struct gkyl_array_copy_func *fout = gkyl_malloc(sizeof(*fout));
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
      fprintf(stderr, "bctype %d not supported in bc_basic_gyrokinetic. Exiting.\n",bctype);
      assert(false);
      break;
  }
  fout->ctx = ctx;
  fout->ctx_on_dev = fout->ctx;

  fout->flags = 0;
  GKYL_CLEAR_CU_ALLOC(fout->flags);
  fout->on_dev = fout; // CPU function obj points to itself.
  return fout;
}

struct gkyl_bc_basic_gyrokinetic*
gkyl_bc_basic_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, enum gkyl_gyrokinetic_bc_type bctype,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, int num_comp, int cdim, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_bc_basic_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_bc_basic_gyrokinetic));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->bctype = bctype;
  up->use_gpu = use_gpu;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;

  // Create function applied to array contents (DG coefficients) when
  // copying to/from buffer.
  up->array_copy_func = gkyl_bc_basic_gyrokinetic_create_arr_copy_func(dir, edge,
    cdim, up->bctype, basis, num_comp, use_gpu);
  return up;
}

void
gkyl_bc_basic_gyrokinetic_buffer_fixed_func(const struct gkyl_bc_basic_gyrokinetic *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  if (up->bctype == GKYL_BC_GK_SPECIES_FIXED_FUNC)
    gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
      up->ghost_r, up->array_copy_func->on_dev);    
}

void
gkyl_bc_basic_gyrokinetic_advance(const struct gkyl_bc_basic_gyrokinetic *up,
  struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  // Apply BC in two steps:
  // 1) Copy skin to buffer while applying array_copy_func.
  switch (up->bctype) {
    case GKYL_BC_GK_SPECIES_COPY:
    case GKYL_BC_GK_SPECIES_ABSORB:
      gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
                                   up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_GK_SPECIES_REFLECT:
      gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_arr, up->cdim,
                                        up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_GK_SPECIES_BOUNDARY_VALUE:
      gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
                                   up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_GK_FIELD_BOUNDARY_VALUE:
      gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
                                   up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_GK_SPECIES_FIXED_FUNC: // if BC is fixed func, do nothing, buffer already full
      break;

    default:
      fprintf(stderr, "bctype %d not supported in bc_basic_gyrokinetic. Exiting.\n",up->bctype);
      assert(false);
      break;
  }
  // 2) Copy from buffer to ghost.
  gkyl_array_copy_from_buffer(f_arr, buff_arr->data, up->ghost_r);
}

void gkyl_bc_basic_gyrokinetic_release(struct gkyl_bc_basic_gyrokinetic *up)
{
  // Release memory associated with array_copy_func.
  if (up->use_gpu) {
    gkyl_cu_free(up->array_copy_func->ctx_on_dev);
    gkyl_cu_free(up->array_copy_func->on_dev);
  }
  gkyl_free(up->array_copy_func->ctx);
  gkyl_free(up->array_copy_func);
  // Release updater memory.
  gkyl_free(up);
}
