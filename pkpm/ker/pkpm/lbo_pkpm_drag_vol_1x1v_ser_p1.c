#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1), M0_perp, N_perp].
  // out:           Incremented output distribution functions. 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  const double *M0_x = &f[12]; 
  const double *M0_y = &f[18]; 
  const double *M0_z = &f[24]; 
  const double *N_x = &f[30]; 
  const double *N_y = &f[36]; 
  const double *N_z = &f[42]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double *out_M0_x = &out[12]; 
  double *out_M0_y = &out[18]; 
  double *out_M0_z = &out[24]; 
  double *out_N_x = &out[30]; 
  double *out_N_y = &out[36]; 
  double *out_N_z = &out[42]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrag[6]; 
  alphaDrag[0] = rdvpar*(1.4142135623730951*sumNuUPar[0]-1.4142135623730951*nuSum[0]*wvpar); 
  alphaDrag[1] = rdvpar*(1.4142135623730951*sumNuUPar[1]-1.4142135623730951*nuSum[1]*wvpar); 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dvpar*rdvpar); 
  alphaDrag[3] = -(0.408248290463863*nuSum[1]*dvpar*rdvpar); 

  const double *sumNuVtSq = &nuPrimMomsSum[2]; 
  out_F_0[2] += 0.8660254037844386*F_0[3]*alphaDrag[3]+0.8660254037844386*F_0[2]*alphaDrag[2]+0.8660254037844386*F_0[1]*alphaDrag[1]+0.8660254037844386*F_0[0]*alphaDrag[0]; 
  out_F_0[3] += 0.8660254037844386*F_0[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*F_0[3]+0.8660254037844386*F_0[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*F_0[1]; 
  out_F_0[4] += 1.7320508075688774*alphaDrag[3]*F_0[5]+1.7320508075688772*alphaDrag[2]*F_0[4]+1.9364916731037085*F_0[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*F_0[3]+1.9364916731037085*F_0[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*F_0[2]; 
  out_F_0[5] += 1.7320508075688772*alphaDrag[2]*F_0[5]+1.7320508075688774*alphaDrag[3]*F_0[4]+1.9364916731037085*F_0[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*F_0[3]+1.9364916731037085*F_0[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*F_0[2]; 
  out_G_1[0] += 1.4142135623730951*F_0[1]*sumNuVtSq[1]-1.4142135623730951*G_1[1]*nuSum[1]+1.4142135623730951*F_0[0]*sumNuVtSq[0]-1.4142135623730951*G_1[0]*nuSum[0]; 
  out_G_1[1] += 1.4142135623730951*F_0[0]*sumNuVtSq[1]-1.4142135623730951*G_1[0]*nuSum[1]-1.4142135623730951*nuSum[0]*G_1[1]+1.4142135623730951*sumNuVtSq[0]*F_0[1]; 
  out_G_1[2] += 0.8660254037844386*G_1[3]*alphaDrag[3]-1.4142135623730951*nuSum[1]*G_1[3]+1.4142135623730951*sumNuVtSq[1]*F_0[3]+0.8660254037844386*G_1[2]*alphaDrag[2]-1.4142135623730951*nuSum[0]*G_1[2]+1.4142135623730951*sumNuVtSq[0]*F_0[2]+0.8660254037844386*G_1[1]*alphaDrag[1]+0.8660254037844386*G_1[0]*alphaDrag[0]; 
  out_G_1[3] += 0.8660254037844386*G_1[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*G_1[3]-1.4142135623730951*nuSum[0]*G_1[3]+1.4142135623730951*sumNuVtSq[0]*F_0[3]-1.4142135623730951*nuSum[1]*G_1[2]+1.4142135623730951*sumNuVtSq[1]*F_0[2]+0.8660254037844386*G_1[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*G_1[1]; 
  out_G_1[4] += 1.7320508075688774*alphaDrag[3]*G_1[5]-1.4142135623730951*nuSum[1]*G_1[5]+1.4142135623730951*sumNuVtSq[1]*F_0[5]+1.7320508075688772*alphaDrag[2]*G_1[4]-1.4142135623730951*nuSum[0]*G_1[4]+1.4142135623730951*sumNuVtSq[0]*F_0[4]+1.9364916731037085*G_1[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*G_1[3]+1.9364916731037085*G_1[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*G_1[2]; 
  out_G_1[5] += 1.7320508075688772*alphaDrag[2]*G_1[5]-1.4142135623730951*nuSum[0]*G_1[5]+1.4142135623730951*sumNuVtSq[0]*F_0[5]+1.7320508075688774*alphaDrag[3]*G_1[4]-1.4142135623730951*nuSum[1]*G_1[4]+1.4142135623730951*sumNuVtSq[1]*F_0[4]+1.9364916731037085*G_1[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*G_1[3]+1.9364916731037085*G_1[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*G_1[2]; 
  out_M0_x[0] += -(0.7071067811865475*M0_x[1]*nuSum[1])-0.7071067811865475*M0_x[0]*nuSum[0]; 
  out_M0_x[1] += -(0.7071067811865475*M0_x[0]*nuSum[1])-0.7071067811865475*nuSum[0]*M0_x[1]; 
  out_M0_x[2] += 0.8660254037844386*M0_x[3]*alphaDrag[3]-0.7071067811865475*nuSum[1]*M0_x[3]+0.8660254037844386*M0_x[2]*alphaDrag[2]-0.7071067811865475*nuSum[0]*M0_x[2]+0.8660254037844386*M0_x[1]*alphaDrag[1]+0.8660254037844386*M0_x[0]*alphaDrag[0]; 
  out_M0_x[3] += 0.8660254037844386*M0_x[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*M0_x[3]-0.7071067811865475*nuSum[0]*M0_x[3]-0.7071067811865475*nuSum[1]*M0_x[2]+0.8660254037844386*M0_x[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*M0_x[1]; 
  out_M0_x[4] += 1.7320508075688774*alphaDrag[3]*M0_x[5]-0.7071067811865475*nuSum[1]*M0_x[5]+1.7320508075688772*alphaDrag[2]*M0_x[4]-0.7071067811865475*nuSum[0]*M0_x[4]+1.9364916731037085*M0_x[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*M0_x[3]+1.9364916731037085*M0_x[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*M0_x[2]; 
  out_M0_x[5] += 1.7320508075688772*alphaDrag[2]*M0_x[5]-0.7071067811865475*nuSum[0]*M0_x[5]+1.7320508075688774*alphaDrag[3]*M0_x[4]-0.7071067811865475*nuSum[1]*M0_x[4]+1.9364916731037085*M0_x[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*M0_x[3]+1.9364916731037085*M0_x[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*M0_x[2]; 
  out_M0_y[0] += -(0.7071067811865475*M0_y[1]*nuSum[1])-0.7071067811865475*M0_y[0]*nuSum[0]; 
  out_M0_y[1] += -(0.7071067811865475*M0_y[0]*nuSum[1])-0.7071067811865475*nuSum[0]*M0_y[1]; 
  out_M0_y[2] += 0.8660254037844386*M0_y[3]*alphaDrag[3]-0.7071067811865475*nuSum[1]*M0_y[3]+0.8660254037844386*M0_y[2]*alphaDrag[2]-0.7071067811865475*nuSum[0]*M0_y[2]+0.8660254037844386*M0_y[1]*alphaDrag[1]+0.8660254037844386*M0_y[0]*alphaDrag[0]; 
  out_M0_y[3] += 0.8660254037844386*M0_y[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*M0_y[3]-0.7071067811865475*nuSum[0]*M0_y[3]-0.7071067811865475*nuSum[1]*M0_y[2]+0.8660254037844386*M0_y[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*M0_y[1]; 
  out_M0_y[4] += 1.7320508075688774*alphaDrag[3]*M0_y[5]-0.7071067811865475*nuSum[1]*M0_y[5]+1.7320508075688772*alphaDrag[2]*M0_y[4]-0.7071067811865475*nuSum[0]*M0_y[4]+1.9364916731037085*M0_y[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*M0_y[3]+1.9364916731037085*M0_y[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*M0_y[2]; 
  out_M0_y[5] += 1.7320508075688772*alphaDrag[2]*M0_y[5]-0.7071067811865475*nuSum[0]*M0_y[5]+1.7320508075688774*alphaDrag[3]*M0_y[4]-0.7071067811865475*nuSum[1]*M0_y[4]+1.9364916731037085*M0_y[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*M0_y[3]+1.9364916731037085*M0_y[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*M0_y[2]; 
  out_M0_z[0] += -(0.7071067811865475*M0_z[1]*nuSum[1])-0.7071067811865475*M0_z[0]*nuSum[0]; 
  out_M0_z[1] += -(0.7071067811865475*M0_z[0]*nuSum[1])-0.7071067811865475*nuSum[0]*M0_z[1]; 
  out_M0_z[2] += 0.8660254037844386*M0_z[3]*alphaDrag[3]-0.7071067811865475*nuSum[1]*M0_z[3]+0.8660254037844386*M0_z[2]*alphaDrag[2]-0.7071067811865475*nuSum[0]*M0_z[2]+0.8660254037844386*M0_z[1]*alphaDrag[1]+0.8660254037844386*M0_z[0]*alphaDrag[0]; 
  out_M0_z[3] += 0.8660254037844386*M0_z[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*M0_z[3]-0.7071067811865475*nuSum[0]*M0_z[3]-0.7071067811865475*nuSum[1]*M0_z[2]+0.8660254037844386*M0_z[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*M0_z[1]; 
  out_M0_z[4] += 1.7320508075688774*alphaDrag[3]*M0_z[5]-0.7071067811865475*nuSum[1]*M0_z[5]+1.7320508075688772*alphaDrag[2]*M0_z[4]-0.7071067811865475*nuSum[0]*M0_z[4]+1.9364916731037085*M0_z[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*M0_z[3]+1.9364916731037085*M0_z[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*M0_z[2]; 
  out_M0_z[5] += 1.7320508075688772*alphaDrag[2]*M0_z[5]-0.7071067811865475*nuSum[0]*M0_z[5]+1.7320508075688774*alphaDrag[3]*M0_z[4]-0.7071067811865475*nuSum[1]*M0_z[4]+1.9364916731037085*M0_z[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*M0_z[3]+1.9364916731037085*M0_z[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*M0_z[2]; 
  out_N_x[0] += 2.8284271247461907*M0_x[1]*sumNuVtSq[1]-2.1213203435596424*N_x[1]*nuSum[1]+2.8284271247461907*M0_x[0]*sumNuVtSq[0]-2.1213203435596424*N_x[0]*nuSum[0]; 
  out_N_x[1] += 2.8284271247461907*M0_x[0]*sumNuVtSq[1]-2.1213203435596424*N_x[0]*nuSum[1]-2.1213203435596424*nuSum[0]*N_x[1]+2.8284271247461907*sumNuVtSq[0]*M0_x[1]; 
  out_N_x[2] += 0.8660254037844386*N_x[3]*alphaDrag[3]-2.1213203435596424*nuSum[1]*N_x[3]+2.8284271247461907*sumNuVtSq[1]*M0_x[3]+0.8660254037844386*N_x[2]*alphaDrag[2]-2.1213203435596424*nuSum[0]*N_x[2]+2.8284271247461907*sumNuVtSq[0]*M0_x[2]+0.8660254037844386*N_x[1]*alphaDrag[1]+0.8660254037844386*N_x[0]*alphaDrag[0]; 
  out_N_x[3] += 0.8660254037844386*N_x[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*N_x[3]-2.1213203435596424*nuSum[0]*N_x[3]+2.8284271247461907*sumNuVtSq[0]*M0_x[3]-2.1213203435596424*nuSum[1]*N_x[2]+2.8284271247461907*sumNuVtSq[1]*M0_x[2]+0.8660254037844386*N_x[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*N_x[1]; 
  out_N_x[4] += 1.7320508075688774*alphaDrag[3]*N_x[5]-2.1213203435596424*nuSum[1]*N_x[5]+2.828427124746191*sumNuVtSq[1]*M0_x[5]+1.7320508075688772*alphaDrag[2]*N_x[4]-2.1213203435596424*nuSum[0]*N_x[4]+2.8284271247461907*sumNuVtSq[0]*M0_x[4]+1.9364916731037085*N_x[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*N_x[3]+1.9364916731037085*N_x[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*N_x[2]; 
  out_N_x[5] += 1.7320508075688772*alphaDrag[2]*N_x[5]-2.1213203435596424*nuSum[0]*N_x[5]+2.8284271247461907*sumNuVtSq[0]*M0_x[5]+1.7320508075688774*alphaDrag[3]*N_x[4]-2.1213203435596424*nuSum[1]*N_x[4]+2.828427124746191*sumNuVtSq[1]*M0_x[4]+1.9364916731037085*N_x[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*N_x[3]+1.9364916731037085*N_x[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*N_x[2]; 
  out_N_y[0] += 2.8284271247461907*M0_y[1]*sumNuVtSq[1]-2.1213203435596424*N_y[1]*nuSum[1]+2.8284271247461907*M0_y[0]*sumNuVtSq[0]-2.1213203435596424*N_y[0]*nuSum[0]; 
  out_N_y[1] += 2.8284271247461907*M0_y[0]*sumNuVtSq[1]-2.1213203435596424*N_y[0]*nuSum[1]-2.1213203435596424*nuSum[0]*N_y[1]+2.8284271247461907*sumNuVtSq[0]*M0_y[1]; 
  out_N_y[2] += 0.8660254037844386*N_y[3]*alphaDrag[3]-2.1213203435596424*nuSum[1]*N_y[3]+2.8284271247461907*sumNuVtSq[1]*M0_y[3]+0.8660254037844386*N_y[2]*alphaDrag[2]-2.1213203435596424*nuSum[0]*N_y[2]+2.8284271247461907*sumNuVtSq[0]*M0_y[2]+0.8660254037844386*N_y[1]*alphaDrag[1]+0.8660254037844386*N_y[0]*alphaDrag[0]; 
  out_N_y[3] += 0.8660254037844386*N_y[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*N_y[3]-2.1213203435596424*nuSum[0]*N_y[3]+2.8284271247461907*sumNuVtSq[0]*M0_y[3]-2.1213203435596424*nuSum[1]*N_y[2]+2.8284271247461907*sumNuVtSq[1]*M0_y[2]+0.8660254037844386*N_y[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*N_y[1]; 
  out_N_y[4] += 1.7320508075688774*alphaDrag[3]*N_y[5]-2.1213203435596424*nuSum[1]*N_y[5]+2.828427124746191*sumNuVtSq[1]*M0_y[5]+1.7320508075688772*alphaDrag[2]*N_y[4]-2.1213203435596424*nuSum[0]*N_y[4]+2.8284271247461907*sumNuVtSq[0]*M0_y[4]+1.9364916731037085*N_y[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*N_y[3]+1.9364916731037085*N_y[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*N_y[2]; 
  out_N_y[5] += 1.7320508075688772*alphaDrag[2]*N_y[5]-2.1213203435596424*nuSum[0]*N_y[5]+2.8284271247461907*sumNuVtSq[0]*M0_y[5]+1.7320508075688774*alphaDrag[3]*N_y[4]-2.1213203435596424*nuSum[1]*N_y[4]+2.828427124746191*sumNuVtSq[1]*M0_y[4]+1.9364916731037085*N_y[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*N_y[3]+1.9364916731037085*N_y[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*N_y[2]; 
  out_N_z[0] += 2.8284271247461907*M0_z[1]*sumNuVtSq[1]-2.1213203435596424*N_z[1]*nuSum[1]+2.8284271247461907*M0_z[0]*sumNuVtSq[0]-2.1213203435596424*N_z[0]*nuSum[0]; 
  out_N_z[1] += 2.8284271247461907*M0_z[0]*sumNuVtSq[1]-2.1213203435596424*N_z[0]*nuSum[1]-2.1213203435596424*nuSum[0]*N_z[1]+2.8284271247461907*sumNuVtSq[0]*M0_z[1]; 
  out_N_z[2] += 0.8660254037844386*N_z[3]*alphaDrag[3]-2.1213203435596424*nuSum[1]*N_z[3]+2.8284271247461907*sumNuVtSq[1]*M0_z[3]+0.8660254037844386*N_z[2]*alphaDrag[2]-2.1213203435596424*nuSum[0]*N_z[2]+2.8284271247461907*sumNuVtSq[0]*M0_z[2]+0.8660254037844386*N_z[1]*alphaDrag[1]+0.8660254037844386*N_z[0]*alphaDrag[0]; 
  out_N_z[3] += 0.8660254037844386*N_z[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*N_z[3]-2.1213203435596424*nuSum[0]*N_z[3]+2.8284271247461907*sumNuVtSq[0]*M0_z[3]-2.1213203435596424*nuSum[1]*N_z[2]+2.8284271247461907*sumNuVtSq[1]*M0_z[2]+0.8660254037844386*N_z[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*N_z[1]; 
  out_N_z[4] += 1.7320508075688774*alphaDrag[3]*N_z[5]-2.1213203435596424*nuSum[1]*N_z[5]+2.828427124746191*sumNuVtSq[1]*M0_z[5]+1.7320508075688772*alphaDrag[2]*N_z[4]-2.1213203435596424*nuSum[0]*N_z[4]+2.8284271247461907*sumNuVtSq[0]*M0_z[4]+1.9364916731037085*N_z[1]*alphaDrag[3]+1.9364916731037085*alphaDrag[1]*N_z[3]+1.9364916731037085*N_z[0]*alphaDrag[2]+1.9364916731037085*alphaDrag[0]*N_z[2]; 
  out_N_z[5] += 1.7320508075688772*alphaDrag[2]*N_z[5]-2.1213203435596424*nuSum[0]*N_z[5]+2.8284271247461907*sumNuVtSq[0]*M0_z[5]+1.7320508075688774*alphaDrag[3]*N_z[4]-2.1213203435596424*nuSum[1]*N_z[4]+2.828427124746191*sumNuVtSq[1]*M0_z[4]+1.9364916731037085*N_z[0]*alphaDrag[3]+1.9364916731037085*alphaDrag[0]*N_z[3]+1.9364916731037085*N_z[1]*alphaDrag[2]+1.9364916731037085*alphaDrag[1]*N_z[2]; 

  return fabs(1.25*alphaDrag[0]) + 3.0*fabs(0.7071067811865475*nuSum[0]); 

} 
