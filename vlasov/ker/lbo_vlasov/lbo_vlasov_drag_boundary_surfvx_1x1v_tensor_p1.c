#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  const double rdv2 = 2.0/dxv[1]; 
  const double *sumNuUx = &nuPrimMomsSum[0]; 
  if (edge == -1) { 
  double alpha_e[2]; 
  alpha_e[0] = -(0.7071067811865475*nuSum[1]*w[1])+0.7071067811865475*nuSum[0]*w[1]+0.7071067811865475*sumNuUx[1]-0.35355339059327373*dxv[1]*nuSum[1]+0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  alpha_e[1] = 0.7071067811865475*nuSum[1]*w[1]+0.7071067811865475*nuSum[0]*w[1]-0.7071067811865475*sumNuUx[1]+0.35355339059327373*dxv[1]*nuSum[1]+0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  double fuL_e[2]; 
  fuL_e[0] = 1.5811388300841895*fSkin[4]+1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0]; 
  fuL_e[1] = 1.5811388300841895*fSkin[5]+1.224744871391589*fSkin[3]+0.7071067811865475*fSkin[1]; 
  double fuR_e[2]; 
  fuR_e[0] = 1.5811388300841895*fEdge[4]-1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0]; 
  fuR_e[1] = 1.5811388300841895*fEdge[5]-1.224744871391589*fEdge[3]+0.7071067811865475*fEdge[1]; 
  double fqL_e_1[2]; 
  fqL_e_1[0] = 0.7071067811865475*fuL_e[0]-0.7071067811865475*fuL_e[1]; 
  fqL_e_1[1] = 0.7071067811865475*fuL_e[1]+0.7071067811865475*fuL_e[0]; 
  double fqR_e_1[2]; 
  fqR_e_1[0] = 0.7071067811865475*fuR_e[0]-0.7071067811865475*fuR_e[1]; 
  fqR_e_1[1] = 0.7071067811865475*fuR_e[1]+0.7071067811865475*fuR_e[0]; 
  double fq_e[2]; 
  if (alpha_e[0] < 0) { 
    fq_e[0] = alpha_e[0]*fqL_e_1[0]; 
  } else { 
    fq_e[0] = alpha_e[0]*fqR_e_1[0]; 
  } 
  if (alpha_e[1] < 0) { 
    fq_e[1] = alpha_e[1]*fqL_e_1[1]; 
  } else { 
    fq_e[1] = alpha_e[1]*fqR_e_1[1]; 
  } 
  double Ghat_e[2]; 
  Ghat_e[0] = 0.7071067811865475*fq_e[1]+0.7071067811865475*fq_e[0]; 
  Ghat_e[1] = 0.7071067811865475*fq_e[1]-0.7071067811865475*fq_e[0]; 
  out[0] += 0.7071067811865475*Ghat_e[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_e[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_e[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat_e[1]*rdv2; 
  out[4] += 1.5811388300841895*Ghat_e[0]*rdv2; 
  out[5] += 1.5811388300841895*Ghat_e[1]*rdv2; 
  } else { 
  double alpha_e[2]; 
  alpha_e[0] = -(0.7071067811865475*nuSum[1]*w[1])+0.7071067811865475*nuSum[0]*w[1]+0.7071067811865475*sumNuUx[1]+0.35355339059327373*dxv[1]*nuSum[1]-0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  alpha_e[1] = 0.7071067811865475*nuSum[1]*w[1]+0.7071067811865475*nuSum[0]*w[1]-0.7071067811865475*sumNuUx[1]-0.35355339059327373*dxv[1]*nuSum[1]-0.35355339059327373*nuSum[0]*dxv[1]-0.7071067811865475*sumNuUx[0]; 
  double fuL_e[2]; 
  fuL_e[0] = 1.5811388300841895*fEdge[4]+1.224744871391589*fEdge[2]+0.7071067811865475*fEdge[0]; 
  fuL_e[1] = 1.5811388300841895*fEdge[5]+1.224744871391589*fEdge[3]+0.7071067811865475*fEdge[1]; 
  double fuR_e[2]; 
  fuR_e[0] = 1.5811388300841895*fSkin[4]-1.224744871391589*fSkin[2]+0.7071067811865475*fSkin[0]; 
  fuR_e[1] = 1.5811388300841895*fSkin[5]-1.224744871391589*fSkin[3]+0.7071067811865475*fSkin[1]; 
  double fqL_e_1[2]; 
  fqL_e_1[0] = 0.7071067811865475*fuL_e[0]-0.7071067811865475*fuL_e[1]; 
  fqL_e_1[1] = 0.7071067811865475*fuL_e[1]+0.7071067811865475*fuL_e[0]; 
  double fqR_e_1[2]; 
  fqR_e_1[0] = 0.7071067811865475*fuR_e[0]-0.7071067811865475*fuR_e[1]; 
  fqR_e_1[1] = 0.7071067811865475*fuR_e[1]+0.7071067811865475*fuR_e[0]; 
  double fq_e[2]; 
  if (alpha_e[0] < 0) { 
    fq_e[0] = alpha_e[0]*fqL_e_1[0]; 
  } else { 
    fq_e[0] = alpha_e[0]*fqR_e_1[0]; 
  } 
  if (alpha_e[1] < 0) { 
    fq_e[1] = alpha_e[1]*fqL_e_1[1]; 
  } else { 
    fq_e[1] = alpha_e[1]*fqR_e_1[1]; 
  } 
  double Ghat_e[2]; 
  Ghat_e[0] = 0.7071067811865475*fq_e[1]+0.7071067811865475*fq_e[0]; 
  Ghat_e[1] = 0.7071067811865475*fq_e[1]-0.7071067811865475*fq_e[0]; 
  out[0] += -(0.7071067811865475*Ghat_e[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat_e[1]*rdv2); 
  out[2] += 1.224744871391589*Ghat_e[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat_e[1]*rdv2; 
  out[4] += -(1.5811388300841895*Ghat_e[0]*rdv2); 
  out[5] += -(1.5811388300841895*Ghat_e[1]*rdv2); 
  } 

  return 0.;

} 
