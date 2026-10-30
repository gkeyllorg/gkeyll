#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_drag_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvpar = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrag[12]; 
  alphaDrag[0] = rdvpar*(1.4142135623730951*sumNuUPar[0]-1.4142135623730951*nuSum[0]*wvpar); 
  alphaDrag[1] = rdvpar*(1.4142135623730951*sumNuUPar[1]-1.4142135623730951*nuSum[1]*wvpar); 
  alphaDrag[2] = rdvpar*(1.4142135623730951*sumNuUPar[2]-1.4142135623730951*nuSum[2]*wvpar); 
  alphaDrag[3] = -(0.408248290463863*nuSum[0]*dvpar*rdvpar); 
  alphaDrag[4] = rdvpar*(1.4142135623730951*sumNuUPar[3]-1.4142135623730951*nuSum[3]*wvpar); 
  alphaDrag[5] = -(0.408248290463863*nuSum[1]*dvpar*rdvpar); 
  alphaDrag[6] = -(0.408248290463863*nuSum[2]*dvpar*rdvpar); 
  alphaDrag[8] = -(0.408248290463863*nuSum[3]*dvpar*rdvpar); 

  out_F_0[3] += 0.6123724356957944*F_0[8]*alphaDrag[8]+0.6123724356957944*F_0[6]*alphaDrag[6]+0.6123724356957944*F_0[5]*alphaDrag[5]+0.6123724356957944*F_0[4]*alphaDrag[4]+0.6123724356957944*F_0[3]*alphaDrag[3]+0.6123724356957944*F_0[2]*alphaDrag[2]+0.6123724356957944*F_0[1]*alphaDrag[1]+0.6123724356957944*F_0[0]*alphaDrag[0]; 
  out_F_0[5] += 0.6123724356957944*F_0[6]*alphaDrag[8]+0.6123724356957944*alphaDrag[6]*F_0[8]+0.6123724356957944*F_0[3]*alphaDrag[5]+0.6123724356957944*alphaDrag[3]*F_0[5]+0.6123724356957944*F_0[2]*alphaDrag[4]+0.6123724356957944*alphaDrag[2]*F_0[4]+0.6123724356957944*F_0[0]*alphaDrag[1]+0.6123724356957944*alphaDrag[0]*F_0[1]; 
  out_F_0[6] += 0.6123724356957944*F_0[5]*alphaDrag[8]+0.6123724356957944*alphaDrag[5]*F_0[8]+0.6123724356957944*F_0[3]*alphaDrag[6]+0.6123724356957944*alphaDrag[3]*F_0[6]+0.6123724356957944*F_0[1]*alphaDrag[4]+0.6123724356957944*alphaDrag[1]*F_0[4]+0.6123724356957944*F_0[0]*alphaDrag[2]+0.6123724356957944*alphaDrag[0]*F_0[2]; 
  out_F_0[7] += 1.224744871391589*alphaDrag[8]*F_0[11]+1.224744871391589*alphaDrag[6]*F_0[10]+1.224744871391589*alphaDrag[5]*F_0[9]+1.369306393762915*F_0[4]*alphaDrag[8]+1.369306393762915*alphaDrag[4]*F_0[8]+1.224744871391589*alphaDrag[3]*F_0[7]+1.369306393762915*F_0[2]*alphaDrag[6]+1.369306393762915*alphaDrag[2]*F_0[6]+1.369306393762915*F_0[1]*alphaDrag[5]+1.369306393762915*alphaDrag[1]*F_0[5]+1.369306393762915*F_0[0]*alphaDrag[3]+1.369306393762915*alphaDrag[0]*F_0[3]; 
  out_F_0[8] += 0.6123724356957944*F_0[3]*alphaDrag[8]+0.6123724356957944*alphaDrag[3]*F_0[8]+0.6123724356957944*F_0[5]*alphaDrag[6]+0.6123724356957944*alphaDrag[5]*F_0[6]+0.6123724356957944*F_0[0]*alphaDrag[4]+0.6123724356957944*alphaDrag[0]*F_0[4]+0.6123724356957944*F_0[1]*alphaDrag[2]+0.6123724356957944*alphaDrag[1]*F_0[2]; 
  out_F_0[9] += 1.224744871391589*alphaDrag[6]*F_0[11]+1.224744871391589*alphaDrag[8]*F_0[10]+1.224744871391589*alphaDrag[3]*F_0[9]+1.369306393762915*F_0[2]*alphaDrag[8]+1.369306393762915*alphaDrag[2]*F_0[8]+1.224744871391589*alphaDrag[5]*F_0[7]+1.369306393762915*F_0[4]*alphaDrag[6]+1.369306393762915*alphaDrag[4]*F_0[6]+1.369306393762915*F_0[0]*alphaDrag[5]+1.369306393762915*alphaDrag[0]*F_0[5]+1.369306393762915*F_0[1]*alphaDrag[3]+1.369306393762915*alphaDrag[1]*F_0[3]; 
  out_F_0[10] += 1.224744871391589*alphaDrag[5]*F_0[11]+1.224744871391589*alphaDrag[3]*F_0[10]+1.224744871391589*alphaDrag[8]*F_0[9]+1.369306393762915*F_0[1]*alphaDrag[8]+1.369306393762915*alphaDrag[1]*F_0[8]+1.224744871391589*alphaDrag[6]*F_0[7]+1.369306393762915*F_0[0]*alphaDrag[6]+1.369306393762915*alphaDrag[0]*F_0[6]+1.369306393762915*F_0[4]*alphaDrag[5]+1.369306393762915*alphaDrag[4]*F_0[5]+1.369306393762915*F_0[2]*alphaDrag[3]+1.369306393762915*alphaDrag[2]*F_0[3]; 
  out_F_0[11] += 1.224744871391589*alphaDrag[3]*F_0[11]+1.224744871391589*alphaDrag[5]*F_0[10]+1.224744871391589*alphaDrag[6]*F_0[9]+1.224744871391589*F_0[7]*alphaDrag[8]+1.369306393762915*F_0[0]*alphaDrag[8]+1.369306393762915*alphaDrag[0]*F_0[8]+1.369306393762915*F_0[1]*alphaDrag[6]+1.369306393762915*alphaDrag[1]*F_0[6]+1.369306393762915*F_0[2]*alphaDrag[5]+1.369306393762915*alphaDrag[2]*F_0[5]+1.369306393762915*F_0[3]*alphaDrag[4]+1.369306393762915*alphaDrag[3]*F_0[4]; 
  out_G_1[3] += 0.6123724356957944*G_1[8]*alphaDrag[8]+0.6123724356957944*G_1[6]*alphaDrag[6]+0.6123724356957944*G_1[5]*alphaDrag[5]+0.6123724356957944*G_1[4]*alphaDrag[4]+0.6123724356957944*G_1[3]*alphaDrag[3]+0.6123724356957944*G_1[2]*alphaDrag[2]+0.6123724356957944*G_1[1]*alphaDrag[1]+0.6123724356957944*G_1[0]*alphaDrag[0]; 
  out_G_1[5] += 0.6123724356957944*G_1[6]*alphaDrag[8]+0.6123724356957944*alphaDrag[6]*G_1[8]+0.6123724356957944*G_1[3]*alphaDrag[5]+0.6123724356957944*alphaDrag[3]*G_1[5]+0.6123724356957944*G_1[2]*alphaDrag[4]+0.6123724356957944*alphaDrag[2]*G_1[4]+0.6123724356957944*G_1[0]*alphaDrag[1]+0.6123724356957944*alphaDrag[0]*G_1[1]; 
  out_G_1[6] += 0.6123724356957944*G_1[5]*alphaDrag[8]+0.6123724356957944*alphaDrag[5]*G_1[8]+0.6123724356957944*G_1[3]*alphaDrag[6]+0.6123724356957944*alphaDrag[3]*G_1[6]+0.6123724356957944*G_1[1]*alphaDrag[4]+0.6123724356957944*alphaDrag[1]*G_1[4]+0.6123724356957944*G_1[0]*alphaDrag[2]+0.6123724356957944*alphaDrag[0]*G_1[2]; 
  out_G_1[7] += 1.224744871391589*alphaDrag[8]*G_1[11]+1.224744871391589*alphaDrag[6]*G_1[10]+1.224744871391589*alphaDrag[5]*G_1[9]+1.369306393762915*G_1[4]*alphaDrag[8]+1.369306393762915*alphaDrag[4]*G_1[8]+1.224744871391589*alphaDrag[3]*G_1[7]+1.369306393762915*G_1[2]*alphaDrag[6]+1.369306393762915*alphaDrag[2]*G_1[6]+1.369306393762915*G_1[1]*alphaDrag[5]+1.369306393762915*alphaDrag[1]*G_1[5]+1.369306393762915*G_1[0]*alphaDrag[3]+1.369306393762915*alphaDrag[0]*G_1[3]; 
  out_G_1[8] += 0.6123724356957944*G_1[3]*alphaDrag[8]+0.6123724356957944*alphaDrag[3]*G_1[8]+0.6123724356957944*G_1[5]*alphaDrag[6]+0.6123724356957944*alphaDrag[5]*G_1[6]+0.6123724356957944*G_1[0]*alphaDrag[4]+0.6123724356957944*alphaDrag[0]*G_1[4]+0.6123724356957944*G_1[1]*alphaDrag[2]+0.6123724356957944*alphaDrag[1]*G_1[2]; 
  out_G_1[9] += 1.224744871391589*alphaDrag[6]*G_1[11]+1.224744871391589*alphaDrag[8]*G_1[10]+1.224744871391589*alphaDrag[3]*G_1[9]+1.369306393762915*G_1[2]*alphaDrag[8]+1.369306393762915*alphaDrag[2]*G_1[8]+1.224744871391589*alphaDrag[5]*G_1[7]+1.369306393762915*G_1[4]*alphaDrag[6]+1.369306393762915*alphaDrag[4]*G_1[6]+1.369306393762915*G_1[0]*alphaDrag[5]+1.369306393762915*alphaDrag[0]*G_1[5]+1.369306393762915*G_1[1]*alphaDrag[3]+1.369306393762915*alphaDrag[1]*G_1[3]; 
  out_G_1[10] += 1.224744871391589*alphaDrag[5]*G_1[11]+1.224744871391589*alphaDrag[3]*G_1[10]+1.224744871391589*alphaDrag[8]*G_1[9]+1.369306393762915*G_1[1]*alphaDrag[8]+1.369306393762915*alphaDrag[1]*G_1[8]+1.224744871391589*alphaDrag[6]*G_1[7]+1.369306393762915*G_1[0]*alphaDrag[6]+1.369306393762915*alphaDrag[0]*G_1[6]+1.369306393762915*G_1[4]*alphaDrag[5]+1.369306393762915*alphaDrag[4]*G_1[5]+1.369306393762915*G_1[2]*alphaDrag[3]+1.369306393762915*alphaDrag[2]*G_1[3]; 
  out_G_1[11] += 1.224744871391589*alphaDrag[3]*G_1[11]+1.224744871391589*alphaDrag[5]*G_1[10]+1.224744871391589*alphaDrag[6]*G_1[9]+1.224744871391589*G_1[7]*alphaDrag[8]+1.369306393762915*G_1[0]*alphaDrag[8]+1.369306393762915*alphaDrag[0]*G_1[8]+1.369306393762915*G_1[1]*alphaDrag[6]+1.369306393762915*alphaDrag[1]*G_1[6]+1.369306393762915*G_1[2]*alphaDrag[5]+1.369306393762915*alphaDrag[2]*G_1[5]+1.369306393762915*G_1[3]*alphaDrag[4]+1.369306393762915*alphaDrag[3]*G_1[4]; 

  return fabs(0.8838834764831842*alphaDrag[0]); 

} 
