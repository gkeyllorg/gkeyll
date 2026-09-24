#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_vlasov.h>
#include <gkyl_dg_diffusion_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_diffusion_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_diffusion_vlasov *diffusion = container_of(base->on_dev, struct dg_diffusion_vlasov, eqn);
    gkyl_cu_free(diffusion);
  }
  
  struct dg_diffusion_vlasov *diffusion = container_of(base, struct dg_diffusion_vlasov, eqn);
  gkyl_free(diffusion);
}

void
gkyl_dg_diffusion_vlasov_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_diffusion_vlasov_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.D)) {
    gkyl_dg_diffusion_vlasov_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct dg_diffusion_vlasov *diffusion = container_of(eqn, struct dg_diffusion_vlasov, eqn);
  diffusion->auxfields.D = auxin.D;
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_vlasov_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  bool is_diff_const, const bool *diff_in_dir, int diff_order,
  const struct gkyl_range *diff_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_diffusion_vlasov_cu_dev_new(basis, cbasis, is_diff_const, diff_in_dir, diff_order, diff_range);
#endif
  
  struct dg_diffusion_vlasov *diffusion = gkyl_malloc(sizeof(struct dg_diffusion_vlasov));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  diffusion->const_coeff = is_diff_const;
  diffusion->num_basis = basis->num_basis;
  for (int d=0; d<cdim; d++) diffusion->diff_in_dir[d] = diff_in_dir[d];

  const gkyl_dg_diffusion_vlasov_vol_kern_list *vol_kernels;
  const gkyl_dg_diffusion_vlasov_surf_kern_list *surfx_kernels;
  const gkyl_dg_diffusion_vlasov_surf_kern_list *surfy_kernels;
  const gkyl_dg_diffusion_vlasov_surf_kern_list *surfz_kernels; 
  const gkyl_dg_diffusion_vlasov_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_dg_diffusion_vlasov_boundary_surf_kern_list *boundary_surfy_kernels;
  const gkyl_dg_diffusion_vlasov_boundary_surf_kern_list *boundary_surfz_kernels; 

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff                   : ser_vol_kernels_varcoeff                  ;
      surfx_kernels          = diffusion->const_coeff? ser_vlasov_surfx_kernels_constcoeff          : ser_vlasov_surfx_kernels_varcoeff         ;
      surfy_kernels          = diffusion->const_coeff? ser_vlasov_surfy_kernels_constcoeff          : ser_vlasov_surfy_kernels_varcoeff         ;
      surfz_kernels          = diffusion->const_coeff? ser_vlasov_surfz_kernels_constcoeff          : ser_vlasov_surfz_kernels_varcoeff         ;
      boundary_surfx_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfx_kernels_constcoeff : ser_vlasov_boundary_surfx_kernels_varcoeff;
      boundary_surfy_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfy_kernels_constcoeff : ser_vlasov_boundary_surfy_kernels_varcoeff;
      boundary_surfz_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfz_kernels_constcoeff : ser_vlasov_boundary_surfz_kernels_varcoeff;
      break;
  
    default:
      assert(false);
      break;    
  } 

  int dirs_linidx = diffdirs_linidx(diff_in_dir, cdim);

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim, diff_order, poly_order, dirs_linidx);

  diffusion->surf[0] = CKSURF(surfx_kernels, diff_order, cdim, vdim, poly_order);
  if (cdim>1)
    diffusion->surf[1] = CKSURF(surfy_kernels, diff_order, cdim, vdim, poly_order);
  if (cdim>2)
    diffusion->surf[2] = CKSURF(surfz_kernels, diff_order, cdim, vdim, poly_order);

  diffusion->boundary_surf[0] = CKSURF(boundary_surfx_kernels, diff_order, cdim, vdim, poly_order);
  if (cdim>1)
    diffusion->boundary_surf[1] = CKSURF(boundary_surfy_kernels, diff_order, cdim, vdim, poly_order);
  if (cdim>2)
    diffusion->boundary_surf[2] = CKSURF(boundary_surfz_kernels, diff_order, cdim, vdim, poly_order);

  // Ensure non-NULL pointers.
  for (int i=0; i<cdim; ++i) assert(diffusion->surf[i]);

  diffusion->auxfields.D = 0;
  diffusion->diff_range = *diff_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_diffusion_vlasov_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
