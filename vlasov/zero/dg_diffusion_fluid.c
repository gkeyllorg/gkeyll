#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_fluid.h>
#include <gkyl_dg_diffusion_fluid_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_diffusion_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_diffusion_fluid *diffusion =
      container_of(base->on_dev, struct dg_diffusion_fluid, eqn);
    gkyl_cu_free(diffusion);
  }

  struct dg_diffusion_fluid *diffusion = container_of(base, struct dg_diffusion_fluid, eqn);
  gkyl_free(diffusion);
}

void
gkyl_dg_diffusion_fluid_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_fluid_auxfields auxin
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.D)) {
    gkyl_dg_diffusion_fluid_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_diffusion_fluid *diffusion = container_of(eqn, struct dg_diffusion_fluid, eqn);
  diffusion->auxfields.D = auxin.D;
}

struct gkyl_dg_eqn *
gkyl_dg_diffusion_fluid_new(
  const struct gkyl_basis *basis, int num_equations, const bool *diff_in_dir, int diff_order,
  bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_diffusion_fluid_cu_dev_new(basis, num_equations, diff_in_dir, diff_order);
  }
#endif

  struct dg_diffusion_fluid *diffusion = gkyl_malloc(sizeof(struct dg_diffusion_fluid));

  int ndim = basis->ndim;
  int poly_order = basis->poly_order;

  diffusion->num_equations = num_equations;
  diffusion->num_basis = basis->num_basis;
  for (int d = 0; d < ndim; d++) {
    diffusion->diff_in_dir[d] = diff_in_dir[d];
  }

  const gkyl_dg_diffusion_fluid_vol_kern_list *vol_kernels;
  const gkyl_dg_diffusion_fluid_surf_kern_list *surfx_kernels;
  const gkyl_dg_diffusion_fluid_surf_kern_list *surfy_kernels;
  const gkyl_dg_diffusion_fluid_surf_kern_list *surfz_kernels;
  const gkyl_dg_diffusion_fluid_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_dg_diffusion_fluid_boundary_surf_kern_list *boundary_surfy_kernels;
  const gkyl_dg_diffusion_fluid_boundary_surf_kern_list *boundary_surfz_kernels;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels_constcoeff;
      surfx_kernels = ser_surfx_kernels_constcoeff;
      surfy_kernels = ser_surfy_kernels_constcoeff;
      surfz_kernels = ser_surfz_kernels_constcoeff;
      boundary_surfx_kernels = ser_boundary_surfx_kernels_constcoeff;
      boundary_surfy_kernels = ser_boundary_surfy_kernels_constcoeff;
      boundary_surfz_kernels = ser_boundary_surfz_kernels_constcoeff;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = tensor_vol_kernels_constcoeff;
      surfx_kernels = tensor_surfx_kernels_constcoeff;
      surfy_kernels = tensor_surfy_kernels_constcoeff;
      surfz_kernels = tensor_surfz_kernels_constcoeff;
      boundary_surfx_kernels = tensor_boundary_surfx_kernels_constcoeff;
      boundary_surfy_kernels = tensor_boundary_surfy_kernels_constcoeff;
      boundary_surfz_kernels = tensor_boundary_surfz_kernels_constcoeff;
      break;

    default:
      assert(false);
      break;
  }

  int dirs_linidx = diffdirs_linidx(diff_in_dir, ndim);

  diffusion->eqn.num_equations = num_equations;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, ndim, diff_order, poly_order, dirs_linidx);

  diffusion->surf[0] = CKSURF(surfx_kernels, diff_order, ndim, poly_order);
  if (ndim > 1) {
    diffusion->surf[1] = CKSURF(surfy_kernels, diff_order, ndim, poly_order);
  }
  if (ndim > 2) {
    diffusion->surf[2] = CKSURF(surfz_kernels, diff_order, ndim, poly_order);
  }

  diffusion->boundary_surf[0] = CKSURF(boundary_surfx_kernels, diff_order, ndim, poly_order);
  if (ndim > 1) {
    diffusion->boundary_surf[1] = CKSURF(boundary_surfy_kernels, diff_order, ndim, poly_order);
  }
  if (ndim > 2) {
    diffusion->boundary_surf[2] = CKSURF(boundary_surfz_kernels, diff_order, ndim, poly_order);
  }

  // Ensure non-NULL pointers.
  for (int i = 0; i < ndim; ++i) {
    assert(diffusion->surf[i]);
  }

  diffusion->auxfields.D = 0;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_diffusion_fluid_free);
  diffusion->eqn.on_dev = &diffusion->eqn;

  return &diffusion->eqn;
}
