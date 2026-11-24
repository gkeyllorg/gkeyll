#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 

  const double nuq0 = 0.7071067811865475*nuSum[0]-0.7071067811865475*nuSum[1]; 
  const double a0_vx = rdvx2*(0.7071067811865475*nuUSum[0]-0.7071067811865475*nuUSum[1] - nuq0*w[1]); 
  const double b0_vx = -0.5*rdvx2*nuq0*dxv[1]; 
  double f0[3]; 
  f0[0] = 0.7071067811865475*f[0]-0.7071067811865475*f[1]; 
  f0[1] = 0.7071067811865475*f[2]-0.7071067811865475*f[3]; 
  f0[2] = 0.7071067811865475*f[4]-0.7071067811865475*f[5]; 
  double G0[3]; 
  G0[1] = f0[1]*b0_vx+1.7320508075688772*f0[0]*a0_vx; 
  G0[2] = 2.0*f0[2]*b0_vx+2.23606797749979*f0[0]*b0_vx+3.872983346207417*f0[1]*a0_vx; 
  const double nuq1 = 0.7071067811865475*nuSum[1]+0.7071067811865475*nuSum[0]; 
  const double a1_vx = rdvx2*(0.7071067811865475*nuUSum[1]+0.7071067811865475*nuUSum[0] - nuq1*w[1]); 
  const double b1_vx = -0.5*rdvx2*nuq1*dxv[1]; 
  double f1[3]; 
  f1[0] = 0.7071067811865475*f[1]+0.7071067811865475*f[0]; 
  f1[1] = 0.7071067811865475*f[3]+0.7071067811865475*f[2]; 
  f1[2] = 0.7071067811865475*f[5]+0.7071067811865475*f[4]; 
  double G1[3]; 
  G1[1] = f1[1]*b1_vx+1.7320508075688772*f1[0]*a1_vx; 
  G1[2] = 2.0*f1[2]*b1_vx+2.23606797749979*f1[0]*b1_vx+3.872983346207417*f1[1]*a1_vx; 
  out[2] += 0.7071067811865475*G1[1]+0.7071067811865475*G0[1]; 
  out[3] += 0.7071067811865475*G1[1]-0.7071067811865475*G0[1]; 
  out[4] += 0.7071067811865475*G1[2]+0.7071067811865475*G0[2]; 
  out[5] += 0.7071067811865475*G1[2]-0.7071067811865475*G0[2]; 

  return fabs(1.7677669529663687*nuUSum[0]*rdvx2-1.7677669529663687*nuSum[0]*w[1]*rdvx2); 

} 
