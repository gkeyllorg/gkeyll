#include <gkyl_bc_basic.h>
#include <gkyl_bc_basic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

// Private function to create a pointer to the function that applies the BC,
// i.e., the array_copy_func applied to expansion coefficients in ghost cell.
struct gkyl_array_copy_func*
gkyl_bc_basic_create_arr_copy_func(int dir, enum gkyl_edge_loc edge, int cdim, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis *basis, int ncomp, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_bc_basic_create_arr_copy_func_cu(dir, edge, cdim, bctype, basis, ncomp);
#endif

  struct dg_bc_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->edge = edge;
  ctx->cdim = cdim;
  ctx->ncomp = ncomp;

  struct gkyl_array_copy_func *fout = gkyl_malloc(sizeof(*fout));
  switch (bctype) {
    case GKYL_BC_COPY:
    case GKYL_BC_FIXED_FUNC:
      fout->func = copy_bc;
      break;
      
    case GKYL_BC_ABSORB:
      fout->func = species_absorb_bc;
      break;

    case GKYL_BC_REFLECT:
      fout->func = reflect_bc;
      break;

    case GKYL_BC_DISTF_REFLECT:
      fout->func = species_reflect_bc;
      break;

    case GKYL_BC_CONF_BOUNDARY_VALUE:
      fout->func = conf_boundary_value_bc;
      break;

    // Maxwell's perfect electrical conductor (zero normal B and zero tangent E)
    case GKYL_BC_MAXWELL_PEC:
      fout->func = maxwell_pec_bc;
      break;

    // Maxwell's symmetry BC (zero normal E and zero tangent B)
    case GKYL_BC_MAXWELL_SYM:
      fout->func = maxwell_sym_bc;
      break;

    // Reservoir Maxwell's BCs for heat flux problem
    // Based on Roberg-Clark et al. PRL 2018
    // NOTE: ONLY WORKS WITH X BOUNDARY 
    case GKYL_BC_MAXWELL_RESERVOIR:
      fout->func = maxwell_reservoir_bc;
      break;

    // PKPM Reflecting wall for distribution function
    case GKYL_BC_PKPM_SPECIES_REFLECT:
      fout->func = pkpm_species_reflect_bc;
      break;    

    // PKPM Reflecting wall for momentum
    case GKYL_BC_PKPM_MOM_REFLECT:
      fout->func = pkpm_mom_reflect_bc;
      break;    

    // PKPM No-slip wall for momentum
    case GKYL_BC_PKPM_MOM_NO_SLIP:
      fout->func = pkpm_mom_no_slip_bc;
      break;   

    // Euler Reflecting wall 
    case GKYL_BC_EULER_REFLECT:
      fout->func = euler_reflect_bc;
      break;    

    // Euler No-slip wall 
    case GKYL_BC_EULER_NO_SLIP:
      fout->func = euler_no_slip_bc;
      break;  

    default:
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

struct gkyl_bc_basic*
gkyl_bc_basic_new(int dir, enum gkyl_edge_loc edge, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis *basis, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, int num_comp, int cdim, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_basic *up = gkyl_malloc(sizeof(struct gkyl_bc_basic));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->bctype = bctype;
  up->use_gpu = use_gpu;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;

  // Create function applied to array contents (DG coefficients) when
  // copying to/from buffer.
  up->array_copy_func = gkyl_bc_basic_create_arr_copy_func(dir, edge, cdim, up->bctype, basis, num_comp, use_gpu);
  return up;
}

void
gkyl_bc_basic_buffer_fixed_func(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  if (up->bctype == GKYL_BC_FIXED_FUNC)
    gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
                                 up->skin_r, up->array_copy_func->on_dev);    
}

void
gkyl_bc_basic_advance(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  // Apply BC in two steps:
  // 1) Copy skin to buffer while applying array_copy_func.
  switch (up->bctype) {
    case GKYL_BC_COPY:
    case GKYL_BC_ABSORB:
    case GKYL_BC_REFLECT:
    case GKYL_BC_MAXWELL_PEC:
    case GKYL_BC_MAXWELL_SYM:
    case GKYL_BC_MAXWELL_RESERVOIR:
    case GKYL_BC_PKPM_MOM_REFLECT:
    case GKYL_BC_PKPM_MOM_NO_SLIP:
    case GKYL_BC_EULER_REFLECT:
    case GKYL_BC_EULER_NO_SLIP:
      gkyl_array_copy_to_buffer_fn(buff_arr->data, f_arr,
                                   up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_DISTF_REFLECT:
    case GKYL_BC_PKPM_SPECIES_REFLECT:
      gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_arr, up->dir+up->cdim,
                                        up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_CONF_BOUNDARY_VALUE:
      gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_arr, up->dir,
                                        up->skin_r, up->array_copy_func->on_dev);
      break;

    case GKYL_BC_FIXED_FUNC: // if BC is fixed func, do nothing, buffer already full
      break;
  }
  // 2) Copy from buffer to ghost.
  gkyl_array_copy_from_buffer(f_arr, buff_arr->data, up->ghost_r);
}

void gkyl_bc_basic_release(struct gkyl_bc_basic *up)
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
