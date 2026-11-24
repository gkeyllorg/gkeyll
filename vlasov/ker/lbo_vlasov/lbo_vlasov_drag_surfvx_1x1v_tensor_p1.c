#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // Drag on the velocity faces v_x = -1, +1: f upwinded at the face quadrature nodes (evaluated one 
  // direction at a time) times the drag speed, projected back. 
  const double rdv2 = 2.0/dxv[1]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 
  double alpha_l[2]; 
  alpha_l[0] = -(0.7071067811865475*nuSum[1]*w[1])+0.7071067811865475*nuSum[0]*w[1]+0.7071067811865475*sumNuUx[1]+0.35355339059327373*dxv[1]*nuSum[1]-0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  alpha_l[1] = 0.7071067811865475*nuSum[1]*w[1]+0.7071067811865475*nuSum[0]*w[1]-0.7071067811865475*sumNuUx[1]-0.35355339059327373*dxv[1]*nuSum[1]-0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  double alpha_r[2]; 
  alpha_r[0] = -(0.7071067811865475*nuSum[1]*w[1])+0.7071067811865475*nuSum[0]*w[1]+0.7071067811865475*sumNuUx[1]-0.35355339059327373*dxv[1]*nuSum[1]+0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  alpha_r[1] = 0.7071067811865475*nuSum[1]*w[1]+0.7071067811865475*nuSum[0]*w[1]-0.7071067811865475*sumNuUx[1]+0.35355339059327373*dxv[1]*nuSum[1]+0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  double fuL_l[2]; 
  fuL_l[0] = 1.5811388300841895*fl[4]+1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  fuL_l[1] = 1.5811388300841895*fl[5]+1.224744871391589*fl[3]+0.7071067811865475*fl[1]; 
  double fuR_l[2]; 
  fuR_l[0] = 1.5811388300841895*fc[4]-1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fuR_l[1] = 1.5811388300841895*fc[5]-1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  double fqL_l_1[2]; 
  fqL_l_1[0] = 0.7071067811865475*fuL_l[0]-0.7071067811865475*fuL_l[1]; 
  fqL_l_1[1] = 0.7071067811865475*fuL_l[1]+0.7071067811865475*fuL_l[0]; 
  double fqR_l_1[2]; 
  fqR_l_1[0] = 0.7071067811865475*fuR_l[0]-0.7071067811865475*fuR_l[1]; 
  fqR_l_1[1] = 0.7071067811865475*fuR_l[1]+0.7071067811865475*fuR_l[0]; 
  double fq_l[2]; 
  if (alpha_l[0] < 0) { 
    fq_l[0] = alpha_l[0]*fqL_l_1[0]; 
  } else { 
    fq_l[0] = alpha_l[0]*fqR_l_1[0]; 
  } 
  if (alpha_l[1] < 0) { 
    fq_l[1] = alpha_l[1]*fqL_l_1[1]; 
  } else { 
    fq_l[1] = alpha_l[1]*fqR_l_1[1]; 
  } 
  double Ghat_l[2]; 
  Ghat_l[0] = 0.7071067811865475*fq_l[1]+0.7071067811865475*fq_l[0]; 
  Ghat_l[1] = 0.7071067811865475*fq_l[1]-0.7071067811865475*fq_l[0]; 
  double fuL_r[2]; 
  fuL_r[0] = 1.5811388300841895*fc[4]+1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fuL_r[1] = 1.5811388300841895*fc[5]+1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  double fuR_r[2]; 
  fuR_r[0] = 1.5811388300841895*fr[4]-1.224744871391589*fr[2]+0.7071067811865475*fr[0]; 
  fuR_r[1] = 1.5811388300841895*fr[5]-1.224744871391589*fr[3]+0.7071067811865475*fr[1]; 
  double fqL_r_1[2]; 
  fqL_r_1[0] = 0.7071067811865475*fuL_r[0]-0.7071067811865475*fuL_r[1]; 
  fqL_r_1[1] = 0.7071067811865475*fuL_r[1]+0.7071067811865475*fuL_r[0]; 
  double fqR_r_1[2]; 
  fqR_r_1[0] = 0.7071067811865475*fuR_r[0]-0.7071067811865475*fuR_r[1]; 
  fqR_r_1[1] = 0.7071067811865475*fuR_r[1]+0.7071067811865475*fuR_r[0]; 
  double fq_r[2]; 
  if (alpha_r[0] < 0) { 
    fq_r[0] = alpha_r[0]*fqL_r_1[0]; 
  } else { 
    fq_r[0] = alpha_r[0]*fqR_r_1[0]; 
  } 
  if (alpha_r[1] < 0) { 
    fq_r[1] = alpha_r[1]*fqL_r_1[1]; 
  } else { 
    fq_r[1] = alpha_r[1]*fqR_r_1[1]; 
  } 
  double Ghat_r[2]; 
  Ghat_r[0] = 0.7071067811865475*fq_r[1]+0.7071067811865475*fq_r[0]; 
  Ghat_r[1] = 0.7071067811865475*fq_r[1]-0.7071067811865475*fq_r[0]; 
  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[5] += 1.5811388300841895*Ghat_r[1]*rdv2-1.5811388300841895*Ghat_l[1]*rdv2; 

  return 0.;

} 
