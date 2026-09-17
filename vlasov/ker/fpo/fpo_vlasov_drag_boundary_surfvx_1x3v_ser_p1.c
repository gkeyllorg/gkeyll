#include <gkyl_fpo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 


GKYL_CU_DH double fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1(const double* dxv,
  const double *alpha_surf_Edge, const double *alpha_surf_Skin,
  const double *sgn_alpha_surf_Edge, const double *sgn_alpha_surf_Skin,
  const int *const_sgn_alpha_Edge, const int *const_sgn_alpha_Skin,
  const int edge, const double *fEdge, const double *fSkin,
  double* GKYL_RESTRICT out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // alpha_surf_Edge,Skin: Surface expansion of drag coefficient on lower boundary of Skin/Edge cell. 
  // sgn_alpha_Edge,Skin: sign(alpha_surf_l,r) at quadrature points. 
  // const_sgn_alpha_Edge,Skin: Boolean array true if sign(alpha_surf_l,r) is only one sign. 
  // fEdge, fSkin: Distribution function in left, center, and right cells. 
  // out: Incremented output. 


  // Index into drag coefficient surface expansion arrays 
  const double *drag_coeff_surf_Edge = &alpha_surf_Edge[0]; 
  const double *drag_coeff_surf_Skin = &alpha_surf_Skin[0]; 
  const double *sgn_drag_coeff_surf_Edge = &sgn_alpha_surf_Edge[0]; 
  const double *sgn_drag_coeff_surf_Skin = &sgn_alpha_surf_Skin[0]; 
  const int *const_sgn_drag_coeff_Edge = &const_sgn_alpha_Edge[0]; 
  const int *const_sgn_drag_coeff_Skin = &const_sgn_alpha_Skin[0]; 
  double dv_inv = 2.0/dxv[1]; 

  double cflFreq = 0.0; 

  if (edge == -1) { 
  const double *drag_coeff_surf_R = drag_coeff_surf_Edge; 
  const double *sgn_drag_coeff_surf_R = sgn_drag_coeff_surf_Edge; 
  const int *const_sgn_drag_coeff_R = const_sgn_drag_coeff_Edge; 
  double fUp_R[8] = {0.0}; 
  if (const_sgn_drag_coeff_R[0] == 1) { 
    if (sgn_drag_coeff_surf_R[0] == 1.0) { 
  fUp_R[0] = 1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0]; 
  fUp_R[1] = 1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[1]; 
  fUp_R[2] = 1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[3]; 
  fUp_R[3] = 1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[4]; 
  fUp_R[4] = 1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[6]; 
  fUp_R[5] = 1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8]; 
  fUp_R[6] = 1.224744871391589*fSkin[14]+0.7071067811865475*fSkin[10]; 
  fUp_R[7] = 1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[13]; 
    } else { 
  fUp_R[0] = 0.7071067811865475*fEdge[0]-1.224744871391589*fEdge[2]; 
  fUp_R[1] = 0.7071067811865475*fEdge[1]-1.224744871391589*fEdge[5]; 
  fUp_R[2] = 0.7071067811865475*fEdge[3]-1.224744871391589*fEdge[7]; 
  fUp_R[3] = 0.7071067811865475*fEdge[4]-1.224744871391589*fEdge[9]; 
  fUp_R[4] = 0.7071067811865475*fEdge[6]-1.224744871391589*fEdge[11]; 
  fUp_R[5] = 0.7071067811865475*fEdge[8]-1.224744871391589*fEdge[12]; 
  fUp_R[6] = 0.7071067811865475*fEdge[10]-1.224744871391589*fEdge[14]; 
  fUp_R[7] = 0.7071067811865475*fEdge[13]-1.224744871391589*fEdge[15]; 
   } 
  } else { 
  double fC_r[8] = {0.0}; 
  double fR_l[8] = {0.0}; 
  double sgn_drag_coeff_Up_R[8] = {0.0}; 
  ser_4x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_R, sgn_drag_coeff_Up_R); 

  fC_r[0] = 1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0]; 
  fC_r[1] = 1.224744871391589*fSkin[5]+0.7071067811865475*fSkin[1]; 
  fC_r[2] = 1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[3]; 
  fC_r[3] = 1.224744871391589*fSkin[9]+0.7071067811865475*fSkin[4]; 
  fC_r[4] = 1.224744871391589*fSkin[11]+0.7071067811865475*fSkin[6]; 
  fC_r[5] = 1.224744871391589*fSkin[12]+0.7071067811865475*fSkin[8]; 
  fC_r[6] = 1.224744871391589*fSkin[14]+0.7071067811865475*fSkin[10]; 
  fC_r[7] = 1.224744871391589*fSkin[15]+0.7071067811865475*fSkin[13]; 

  fR_l[0] = 0.7071067811865475*fEdge[0]-1.224744871391589*fEdge[2]; 
  fR_l[1] = 0.7071067811865475*fEdge[1]-1.224744871391589*fEdge[5]; 
  fR_l[2] = 0.7071067811865475*fEdge[3]-1.224744871391589*fEdge[7]; 
  fR_l[3] = 0.7071067811865475*fEdge[4]-1.224744871391589*fEdge[9]; 
  fR_l[4] = 0.7071067811865475*fEdge[6]-1.224744871391589*fEdge[11]; 
  fR_l[5] = 0.7071067811865475*fEdge[8]-1.224744871391589*fEdge[12]; 
  fR_l[6] = 0.7071067811865475*fEdge[10]-1.224744871391589*fEdge[14]; 
  fR_l[7] = 0.7071067811865475*fEdge[13]-1.224744871391589*fEdge[15]; 

  fUp_R[0] = -(0.1767766952966368*(fR_l[7]*sgn_drag_coeff_Up_R[7]+fR_l[6]*sgn_drag_coeff_Up_R[6]+fR_l[5]*sgn_drag_coeff_Up_R[5]+fR_l[4]*sgn_drag_coeff_Up_R[4]+fR_l[3]*sgn_drag_coeff_Up_R[3]+fR_l[2]*sgn_drag_coeff_Up_R[2]+fR_l[1]*sgn_drag_coeff_Up_R[1]+fR_l[0]*sgn_drag_coeff_Up_R[0]))+0.1767766952966368*(fC_r[7]*sgn_drag_coeff_Up_R[7]+fC_r[6]*sgn_drag_coeff_Up_R[6]+fC_r[5]*sgn_drag_coeff_Up_R[5]+fC_r[4]*sgn_drag_coeff_Up_R[4]+fC_r[3]*sgn_drag_coeff_Up_R[3]+fC_r[2]*sgn_drag_coeff_Up_R[2]+fC_r[1]*sgn_drag_coeff_Up_R[1]+fC_r[0]*sgn_drag_coeff_Up_R[0])+0.5*(fR_l[0]+fC_r[0]); 
  fUp_R[1] = -(0.1767766952966368*(fR_l[6]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[6]*fR_l[7]+fR_l[3]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[3]*fR_l[5]+fR_l[2]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[2]*fR_l[4]+fR_l[0]*sgn_drag_coeff_Up_R[1]+sgn_drag_coeff_Up_R[0]*fR_l[1]))+0.1767766952966368*(fC_r[6]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[6]*fC_r[7]+fC_r[3]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[3]*fC_r[5]+fC_r[2]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[2]*fC_r[4]+fC_r[0]*sgn_drag_coeff_Up_R[1]+sgn_drag_coeff_Up_R[0]*fC_r[1])+0.5*(fR_l[1]+fC_r[1]); 
  fUp_R[2] = -(0.1767766952966368*(fR_l[5]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[5]*fR_l[7]+fR_l[3]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[3]*fR_l[6]+fR_l[1]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[1]*fR_l[4]+fR_l[0]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[0]*fR_l[2]))+0.1767766952966368*(fC_r[5]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[5]*fC_r[7]+fC_r[3]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[3]*fC_r[6]+fC_r[1]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[1]*fC_r[4]+fC_r[0]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[0]*fC_r[2])+0.5*(fR_l[2]+fC_r[2]); 
  fUp_R[3] = -(0.1767766952966368*(fR_l[4]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[4]*fR_l[7]+fR_l[2]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[2]*fR_l[6]+fR_l[1]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[1]*fR_l[5]+fR_l[0]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[0]*fR_l[3]))+0.1767766952966368*(fC_r[4]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[4]*fC_r[7]+fC_r[2]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[2]*fC_r[6]+fC_r[1]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[1]*fC_r[5]+fC_r[0]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[0]*fC_r[3])+0.5*(fR_l[3]+fC_r[3]); 
  fUp_R[4] = -(0.1767766952966368*(fR_l[3]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[3]*fR_l[7]+fR_l[5]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[5]*fR_l[6]+fR_l[0]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[0]*fR_l[4]+fR_l[1]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[1]*fR_l[2]))+0.1767766952966368*(fC_r[3]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[3]*fC_r[7]+fC_r[5]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[5]*fC_r[6]+fC_r[0]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[0]*fC_r[4]+fC_r[1]*sgn_drag_coeff_Up_R[2]+sgn_drag_coeff_Up_R[1]*fC_r[2])+0.5*(fR_l[4]+fC_r[4]); 
  fUp_R[5] = -(0.1767766952966368*(fR_l[2]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[2]*fR_l[7]+fR_l[4]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[4]*fR_l[6]+fR_l[0]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[0]*fR_l[5]+fR_l[1]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[1]*fR_l[3]))+0.1767766952966368*(fC_r[2]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[2]*fC_r[7]+fC_r[4]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[4]*fC_r[6]+fC_r[0]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[0]*fC_r[5]+fC_r[1]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[1]*fC_r[3])+0.5*(fR_l[5]+fC_r[5]); 
  fUp_R[6] = -(0.1767766952966368*(fR_l[1]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[1]*fR_l[7]+fR_l[0]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[0]*fR_l[6]+fR_l[4]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[4]*fR_l[5]+fR_l[2]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[2]*fR_l[3]))+0.1767766952966368*(fC_r[1]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[1]*fC_r[7]+fC_r[0]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[0]*fC_r[6]+fC_r[4]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[4]*fC_r[5]+fC_r[2]*sgn_drag_coeff_Up_R[3]+sgn_drag_coeff_Up_R[2]*fC_r[3])+0.5*(fR_l[6]+fC_r[6]); 
  fUp_R[7] = -(0.1767766952966368*(fR_l[0]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[0]*fR_l[7]+fR_l[1]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[1]*fR_l[6]+fR_l[2]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[2]*fR_l[5]+fR_l[3]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[3]*fR_l[4]))+0.1767766952966368*(fC_r[0]*sgn_drag_coeff_Up_R[7]+sgn_drag_coeff_Up_R[0]*fC_r[7]+fC_r[1]*sgn_drag_coeff_Up_R[6]+sgn_drag_coeff_Up_R[1]*fC_r[6]+fC_r[2]*sgn_drag_coeff_Up_R[5]+sgn_drag_coeff_Up_R[2]*fC_r[5]+fC_r[3]*sgn_drag_coeff_Up_R[4]+sgn_drag_coeff_Up_R[3]*fC_r[4])+0.5*(fR_l[7]+fC_r[7]); 

  } 
  double GhatR[8] = {0.0}; 
  GhatR[0] = 0.3535533905932737*(drag_coeff_surf_Edge[7]*fUp_R[7]+drag_coeff_surf_Edge[6]*fUp_R[6]+drag_coeff_surf_Edge[5]*fUp_R[5]+drag_coeff_surf_Edge[4]*fUp_R[4]+drag_coeff_surf_Edge[3]*fUp_R[3]+drag_coeff_surf_Edge[2]*fUp_R[2]+drag_coeff_surf_Edge[1]*fUp_R[1]+drag_coeff_surf_Edge[0]*fUp_R[0]); 
  GhatR[1] = 0.3535533905932737*(drag_coeff_surf_Edge[6]*fUp_R[7]+fUp_R[6]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[3]*fUp_R[5]+fUp_R[3]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[2]*fUp_R[4]+fUp_R[2]*drag_coeff_surf_Edge[4]+drag_coeff_surf_Edge[0]*fUp_R[1]+fUp_R[0]*drag_coeff_surf_Edge[1]); 
  GhatR[2] = 0.3535533905932737*(drag_coeff_surf_Edge[5]*fUp_R[7]+fUp_R[5]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[3]*fUp_R[6]+fUp_R[3]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[1]*fUp_R[4]+fUp_R[1]*drag_coeff_surf_Edge[4]+drag_coeff_surf_Edge[0]*fUp_R[2]+fUp_R[0]*drag_coeff_surf_Edge[2]); 
  GhatR[3] = 0.3535533905932737*(drag_coeff_surf_Edge[4]*fUp_R[7]+fUp_R[4]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[2]*fUp_R[6]+fUp_R[2]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[1]*fUp_R[5]+fUp_R[1]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[0]*fUp_R[3]+fUp_R[0]*drag_coeff_surf_Edge[3]); 
  GhatR[4] = 0.3535533905932737*(drag_coeff_surf_Edge[3]*fUp_R[7]+fUp_R[3]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[5]*fUp_R[6]+fUp_R[5]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[0]*fUp_R[4]+fUp_R[0]*drag_coeff_surf_Edge[4]+drag_coeff_surf_Edge[1]*fUp_R[2]+fUp_R[1]*drag_coeff_surf_Edge[2]); 
  GhatR[5] = 0.3535533905932737*(drag_coeff_surf_Edge[2]*fUp_R[7]+fUp_R[2]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[4]*fUp_R[6]+fUp_R[4]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[0]*fUp_R[5]+fUp_R[0]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[1]*fUp_R[3]+fUp_R[1]*drag_coeff_surf_Edge[3]); 
  GhatR[6] = 0.3535533905932737*(drag_coeff_surf_Edge[1]*fUp_R[7]+fUp_R[1]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[0]*fUp_R[6]+fUp_R[0]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[4]*fUp_R[5]+fUp_R[4]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[2]*fUp_R[3]+fUp_R[2]*drag_coeff_surf_Edge[3]); 
  GhatR[7] = 0.3535533905932737*(drag_coeff_surf_Edge[0]*fUp_R[7]+fUp_R[0]*drag_coeff_surf_Edge[7]+drag_coeff_surf_Edge[1]*fUp_R[6]+fUp_R[1]*drag_coeff_surf_Edge[6]+drag_coeff_surf_Edge[2]*fUp_R[5]+fUp_R[2]*drag_coeff_surf_Edge[5]+drag_coeff_surf_Edge[3]*fUp_R[4]+fUp_R[3]*drag_coeff_surf_Edge[4]); 

  out[0] += -(0.35355339059327373*GhatR[0]*dv_inv); 
  out[1] += -(0.35355339059327373*GhatR[1]*dv_inv); 
  out[2] += -(0.6123724356957945*GhatR[0]*dv_inv); 
  out[3] += -(0.35355339059327373*GhatR[2]*dv_inv); 
  out[4] += -(0.35355339059327373*GhatR[3]*dv_inv); 
  out[5] += -(0.6123724356957945*GhatR[1]*dv_inv); 
  out[6] += -(0.35355339059327373*GhatR[4]*dv_inv); 
  out[7] += -(0.6123724356957945*GhatR[2]*dv_inv); 
  out[8] += -(0.35355339059327373*GhatR[5]*dv_inv); 
  out[9] += -(0.6123724356957945*GhatR[3]*dv_inv); 
  out[10] += -(0.35355339059327373*GhatR[6]*dv_inv); 
  out[11] += -(0.6123724356957945*GhatR[4]*dv_inv); 
  out[12] += -(0.6123724356957945*GhatR[5]*dv_inv); 
  out[13] += -(0.35355339059327373*GhatR[7]*dv_inv); 
  out[14] += -(0.6123724356957945*GhatR[6]*dv_inv); 
  out[15] += -(0.6123724356957945*GhatR[7]*dv_inv); 
  cflFreq = fabs(drag_coeff_surf_Edge[0]);

  } else { 
  const double *drag_coeff_surf_L = drag_coeff_surf_Skin; 
  const double *sgn_drag_coeff_surf_L = sgn_drag_coeff_surf_Skin; 
  const int *const_sgn_drag_coeff_L = const_sgn_drag_coeff_Skin; 
  double fUp_L[8] = {0.0}; 
  if (const_sgn_drag_coeff_L[0] == 1) { 
    if (sgn_drag_coeff_surf_L[0] == 1.0) { 
  fUp_L[0] = 1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0]; 
  fUp_L[1] = 1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[1]; 
  fUp_L[2] = 1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[3]; 
  fUp_L[3] = 1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[4]; 
  fUp_L[4] = 1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[6]; 
  fUp_L[5] = 1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8]; 
  fUp_L[6] = 1.224744871391589*fEdge[14]+0.7071067811865475*fEdge[10]; 
  fUp_L[7] = 1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[13]; 
    } else { 
  fUp_L[0] = 0.7071067811865475*fSkin[0]-1.224744871391589*fSkin[2]; 
  fUp_L[1] = 0.7071067811865475*fSkin[1]-1.224744871391589*fSkin[5]; 
  fUp_L[2] = 0.7071067811865475*fSkin[3]-1.224744871391589*fSkin[7]; 
  fUp_L[3] = 0.7071067811865475*fSkin[4]-1.224744871391589*fSkin[9]; 
  fUp_L[4] = 0.7071067811865475*fSkin[6]-1.224744871391589*fSkin[11]; 
  fUp_L[5] = 0.7071067811865475*fSkin[8]-1.224744871391589*fSkin[12]; 
  fUp_L[6] = 0.7071067811865475*fSkin[10]-1.224744871391589*fSkin[14]; 
  fUp_L[7] = 0.7071067811865475*fSkin[13]-1.224744871391589*fSkin[15]; 
   } 
  } else { 
  double fL_r[8] = {0.0}; 
  double fC_l[8] = {0.0}; 
  double sgn_drag_coeff_Up_L[8] = {0.0}; 
  ser_4x_p1_upwind_quad_to_modal(sgn_drag_coeff_surf_L, sgn_drag_coeff_Up_L); 

  fL_r[0] = 1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0]; 
  fL_r[1] = 1.224744871391589*fEdge[5]+0.7071067811865475*fEdge[1]; 
  fL_r[2] = 1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[3]; 
  fL_r[3] = 1.224744871391589*fEdge[9]+0.7071067811865475*fEdge[4]; 
  fL_r[4] = 1.224744871391589*fEdge[11]+0.7071067811865475*fEdge[6]; 
  fL_r[5] = 1.224744871391589*fEdge[12]+0.7071067811865475*fEdge[8]; 
  fL_r[6] = 1.224744871391589*fEdge[14]+0.7071067811865475*fEdge[10]; 
  fL_r[7] = 1.224744871391589*fEdge[15]+0.7071067811865475*fEdge[13]; 

  fC_l[0] = 0.7071067811865475*fSkin[0]-1.224744871391589*fSkin[2]; 
  fC_l[1] = 0.7071067811865475*fSkin[1]-1.224744871391589*fSkin[5]; 
  fC_l[2] = 0.7071067811865475*fSkin[3]-1.224744871391589*fSkin[7]; 
  fC_l[3] = 0.7071067811865475*fSkin[4]-1.224744871391589*fSkin[9]; 
  fC_l[4] = 0.7071067811865475*fSkin[6]-1.224744871391589*fSkin[11]; 
  fC_l[5] = 0.7071067811865475*fSkin[8]-1.224744871391589*fSkin[12]; 
  fC_l[6] = 0.7071067811865475*fSkin[10]-1.224744871391589*fSkin[14]; 
  fC_l[7] = 0.7071067811865475*fSkin[13]-1.224744871391589*fSkin[15]; 

  fUp_L[0] = 0.1767766952966368*(fL_r[7]*sgn_drag_coeff_Up_L[7]+fL_r[6]*sgn_drag_coeff_Up_L[6]+fL_r[5]*sgn_drag_coeff_Up_L[5]+fL_r[4]*sgn_drag_coeff_Up_L[4]+fL_r[3]*sgn_drag_coeff_Up_L[3]+fL_r[2]*sgn_drag_coeff_Up_L[2]+fL_r[1]*sgn_drag_coeff_Up_L[1]+fL_r[0]*sgn_drag_coeff_Up_L[0])-0.1767766952966368*(fC_l[7]*sgn_drag_coeff_Up_L[7]+fC_l[6]*sgn_drag_coeff_Up_L[6]+fC_l[5]*sgn_drag_coeff_Up_L[5]+fC_l[4]*sgn_drag_coeff_Up_L[4]+fC_l[3]*sgn_drag_coeff_Up_L[3]+fC_l[2]*sgn_drag_coeff_Up_L[2]+fC_l[1]*sgn_drag_coeff_Up_L[1]+fC_l[0]*sgn_drag_coeff_Up_L[0])+0.5*(fL_r[0]+fC_l[0]); 
  fUp_L[1] = 0.1767766952966368*(fL_r[6]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[6]*fL_r[7]+fL_r[3]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[3]*fL_r[5]+fL_r[2]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[2]*fL_r[4]+fL_r[0]*sgn_drag_coeff_Up_L[1]+sgn_drag_coeff_Up_L[0]*fL_r[1])-0.1767766952966368*(fC_l[6]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[6]*fC_l[7]+fC_l[3]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[3]*fC_l[5]+fC_l[2]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[2]*fC_l[4]+fC_l[0]*sgn_drag_coeff_Up_L[1]+sgn_drag_coeff_Up_L[0]*fC_l[1])+0.5*(fL_r[1]+fC_l[1]); 
  fUp_L[2] = 0.1767766952966368*(fL_r[5]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[5]*fL_r[7]+fL_r[3]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[3]*fL_r[6]+fL_r[1]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[1]*fL_r[4]+fL_r[0]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[0]*fL_r[2])-0.1767766952966368*(fC_l[5]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[5]*fC_l[7]+fC_l[3]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[3]*fC_l[6]+fC_l[1]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[1]*fC_l[4]+fC_l[0]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[0]*fC_l[2])+0.5*(fL_r[2]+fC_l[2]); 
  fUp_L[3] = 0.1767766952966368*(fL_r[4]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[4]*fL_r[7]+fL_r[2]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[2]*fL_r[6]+fL_r[1]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[1]*fL_r[5]+fL_r[0]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[0]*fL_r[3])-0.1767766952966368*(fC_l[4]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[4]*fC_l[7]+fC_l[2]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[2]*fC_l[6]+fC_l[1]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[1]*fC_l[5]+fC_l[0]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[0]*fC_l[3])+0.5*(fL_r[3]+fC_l[3]); 
  fUp_L[4] = 0.1767766952966368*(fL_r[3]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[3]*fL_r[7]+fL_r[5]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[5]*fL_r[6]+fL_r[0]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[0]*fL_r[4]+fL_r[1]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[1]*fL_r[2])-0.1767766952966368*(fC_l[3]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[3]*fC_l[7]+fC_l[5]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[5]*fC_l[6]+fC_l[0]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[0]*fC_l[4]+fC_l[1]*sgn_drag_coeff_Up_L[2]+sgn_drag_coeff_Up_L[1]*fC_l[2])+0.5*(fL_r[4]+fC_l[4]); 
  fUp_L[5] = 0.1767766952966368*(fL_r[2]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[2]*fL_r[7]+fL_r[4]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[4]*fL_r[6]+fL_r[0]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[0]*fL_r[5]+fL_r[1]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[1]*fL_r[3])-0.1767766952966368*(fC_l[2]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[2]*fC_l[7]+fC_l[4]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[4]*fC_l[6]+fC_l[0]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[0]*fC_l[5]+fC_l[1]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[1]*fC_l[3])+0.5*(fL_r[5]+fC_l[5]); 
  fUp_L[6] = 0.1767766952966368*(fL_r[1]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[1]*fL_r[7]+fL_r[0]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[0]*fL_r[6]+fL_r[4]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[4]*fL_r[5]+fL_r[2]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[2]*fL_r[3])-0.1767766952966368*(fC_l[1]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[1]*fC_l[7]+fC_l[0]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[0]*fC_l[6]+fC_l[4]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[4]*fC_l[5]+fC_l[2]*sgn_drag_coeff_Up_L[3]+sgn_drag_coeff_Up_L[2]*fC_l[3])+0.5*(fL_r[6]+fC_l[6]); 
  fUp_L[7] = 0.1767766952966368*(fL_r[0]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[0]*fL_r[7]+fL_r[1]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[1]*fL_r[6]+fL_r[2]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[2]*fL_r[5]+fL_r[3]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[3]*fL_r[4])-0.1767766952966368*(fC_l[0]*sgn_drag_coeff_Up_L[7]+sgn_drag_coeff_Up_L[0]*fC_l[7]+fC_l[1]*sgn_drag_coeff_Up_L[6]+sgn_drag_coeff_Up_L[1]*fC_l[6]+fC_l[2]*sgn_drag_coeff_Up_L[5]+sgn_drag_coeff_Up_L[2]*fC_l[5]+fC_l[3]*sgn_drag_coeff_Up_L[4]+sgn_drag_coeff_Up_L[3]*fC_l[4])+0.5*(fL_r[7]+fC_l[7]); 

  } 
  double GhatL[8] = {0.0}; 
  GhatL[0] = 0.3535533905932737*(drag_coeff_surf_Skin[7]*fUp_L[7]+drag_coeff_surf_Skin[6]*fUp_L[6]+drag_coeff_surf_Skin[5]*fUp_L[5]+drag_coeff_surf_Skin[4]*fUp_L[4]+drag_coeff_surf_Skin[3]*fUp_L[3]+drag_coeff_surf_Skin[2]*fUp_L[2]+drag_coeff_surf_Skin[1]*fUp_L[1]+drag_coeff_surf_Skin[0]*fUp_L[0]); 
  GhatL[1] = 0.3535533905932737*(drag_coeff_surf_Skin[6]*fUp_L[7]+fUp_L[6]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[3]*fUp_L[5]+fUp_L[3]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[2]*fUp_L[4]+fUp_L[2]*drag_coeff_surf_Skin[4]+drag_coeff_surf_Skin[0]*fUp_L[1]+fUp_L[0]*drag_coeff_surf_Skin[1]); 
  GhatL[2] = 0.3535533905932737*(drag_coeff_surf_Skin[5]*fUp_L[7]+fUp_L[5]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[3]*fUp_L[6]+fUp_L[3]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[1]*fUp_L[4]+fUp_L[1]*drag_coeff_surf_Skin[4]+drag_coeff_surf_Skin[0]*fUp_L[2]+fUp_L[0]*drag_coeff_surf_Skin[2]); 
  GhatL[3] = 0.3535533905932737*(drag_coeff_surf_Skin[4]*fUp_L[7]+fUp_L[4]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[2]*fUp_L[6]+fUp_L[2]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[1]*fUp_L[5]+fUp_L[1]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[0]*fUp_L[3]+fUp_L[0]*drag_coeff_surf_Skin[3]); 
  GhatL[4] = 0.3535533905932737*(drag_coeff_surf_Skin[3]*fUp_L[7]+fUp_L[3]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[5]*fUp_L[6]+fUp_L[5]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[0]*fUp_L[4]+fUp_L[0]*drag_coeff_surf_Skin[4]+drag_coeff_surf_Skin[1]*fUp_L[2]+fUp_L[1]*drag_coeff_surf_Skin[2]); 
  GhatL[5] = 0.3535533905932737*(drag_coeff_surf_Skin[2]*fUp_L[7]+fUp_L[2]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[4]*fUp_L[6]+fUp_L[4]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[0]*fUp_L[5]+fUp_L[0]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[1]*fUp_L[3]+fUp_L[1]*drag_coeff_surf_Skin[3]); 
  GhatL[6] = 0.3535533905932737*(drag_coeff_surf_Skin[1]*fUp_L[7]+fUp_L[1]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[0]*fUp_L[6]+fUp_L[0]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[4]*fUp_L[5]+fUp_L[4]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[2]*fUp_L[3]+fUp_L[2]*drag_coeff_surf_Skin[3]); 
  GhatL[7] = 0.3535533905932737*(drag_coeff_surf_Skin[0]*fUp_L[7]+fUp_L[0]*drag_coeff_surf_Skin[7]+drag_coeff_surf_Skin[1]*fUp_L[6]+fUp_L[1]*drag_coeff_surf_Skin[6]+drag_coeff_surf_Skin[2]*fUp_L[5]+fUp_L[2]*drag_coeff_surf_Skin[5]+drag_coeff_surf_Skin[3]*fUp_L[4]+fUp_L[3]*drag_coeff_surf_Skin[4]); 

  out[0] += 0.35355339059327373*GhatL[0]*dv_inv; 
  out[1] += 0.35355339059327373*GhatL[1]*dv_inv; 
  out[2] += -(0.6123724356957945*GhatL[0]*dv_inv); 
  out[3] += 0.35355339059327373*GhatL[2]*dv_inv; 
  out[4] += 0.35355339059327373*GhatL[3]*dv_inv; 
  out[5] += -(0.6123724356957945*GhatL[1]*dv_inv); 
  out[6] += 0.35355339059327373*GhatL[4]*dv_inv; 
  out[7] += -(0.6123724356957945*GhatL[2]*dv_inv); 
  out[8] += 0.35355339059327373*GhatL[5]*dv_inv; 
  out[9] += -(0.6123724356957945*GhatL[3]*dv_inv); 
  out[10] += 0.35355339059327373*GhatL[6]*dv_inv; 
  out[11] += -(0.6123724356957945*GhatL[4]*dv_inv); 
  out[12] += -(0.6123724356957945*GhatL[5]*dv_inv); 
  out[13] += 0.35355339059327373*GhatL[7]*dv_inv; 
  out[14] += -(0.6123724356957945*GhatL[6]*dv_inv); 
  out[15] += -(0.6123724356957945*GhatL[7]*dv_inv); 
  cflFreq = fabs(drag_coeff_surf_Skin[0]);

  } 

  return 0.5303300858899105*dv_inv*cflFreq; 
}
