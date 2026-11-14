#include <gkyl_dg_gyrokinetic_passive_kernels.h>
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double dg_gyrokinetic_passive_surfx_1x1v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *flux_surf_l, const double *flux_surf_r, 
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // flux_surf_l: Surface expansion of phase space flux on the left.
  // flux_surf_r: Surface expansion of phase space flux on the right.
  // out: output increment in center cell.

  double rdx2 = 2.0/dxv[0];

  double GhatL[3]= {0.0}; 
  double GhatR[3]= {0.0}; 

  const double *fnodal_l = &flux_surf_l[0]; 
  const double *fnodal_r = &flux_surf_r[0]; 
  GhatL[0] = 0.39283710065919303*fnodal_l[2]+0.6285393610547091*fnodal_l[1]+0.39283710065919303*fnodal_l[0]; 
  GhatL[1] = 0.5270462766947298*fnodal_l[2]-0.5270462766947298*fnodal_l[0]; 
  GhatL[2] = 0.35136418446315326*fnodal_l[2]-0.7027283689263066*fnodal_l[1]+0.35136418446315326*fnodal_l[0]; 
  GhatR[0] = 0.39283710065919303*fnodal_r[2]+0.6285393610547091*fnodal_r[1]+0.39283710065919303*fnodal_r[0]; 
  GhatR[1] = 0.5270462766947298*fnodal_r[2]-0.5270462766947298*fnodal_r[0]; 
  GhatR[2] = 0.35136418446315326*fnodal_r[2]-0.7027283689263066*fnodal_r[1]+0.35136418446315326*fnodal_r[0]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += (-(1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (-(1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[4] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[5] += (-(1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 

  return 0.0; 

} 
