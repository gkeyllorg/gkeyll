#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[2]; 

  const double nuq0 = 0.5*nuSum[3]-0.5*nuSum[2]-0.5*nuSum[1]+0.5*nuSum[0]; 
  const double a0_vx = rdvx2*(0.5*nuUSum[3]-0.5*nuUSum[2]-0.5*nuUSum[1]+0.5*nuUSum[0] - nuq0*w[2]); 
  const double b0_vx = -0.5*rdvx2*nuq0*dxv[2]; 
  double f0[3]; 
  f0[0] = 0.5*f[4]-0.5*f[2]-0.5*f[1]+0.5*f[0]; 
  f0[1] = 0.5*f[8]-0.5*f[6]-0.5*f[5]+0.5*f[3]; 
  f0[2] = 0.5*f[11]-0.5*f[10]-0.5*f[9]+0.5*f[7]; 
  double G0[3]; 
  G0[1] = f0[1]*b0_vx+1.7320508075688772*f0[0]*a0_vx; 
  G0[2] = 2.0*f0[2]*b0_vx+2.23606797749979*f0[0]*b0_vx+3.872983346207417*f0[1]*a0_vx; 
  const double nuq1 = -(0.5*nuSum[3])+0.5*nuSum[2]-0.5*nuSum[1]+0.5*nuSum[0]; 
  const double a1_vx = rdvx2*(-(0.5*nuUSum[3])+0.5*nuUSum[2]-0.5*nuUSum[1]+0.5*nuUSum[0] - nuq1*w[2]); 
  const double b1_vx = -0.5*rdvx2*nuq1*dxv[2]; 
  double f1[3]; 
  f1[0] = -(0.5*f[4])+0.5*f[2]-0.5*f[1]+0.5*f[0]; 
  f1[1] = -(0.5*f[8])+0.5*f[6]-0.5*f[5]+0.5*f[3]; 
  f1[2] = -(0.5*f[11])+0.5*f[10]-0.5*f[9]+0.5*f[7]; 
  double G1[3]; 
  G1[1] = f1[1]*b1_vx+1.7320508075688772*f1[0]*a1_vx; 
  G1[2] = 2.0*f1[2]*b1_vx+2.23606797749979*f1[0]*b1_vx+3.872983346207417*f1[1]*a1_vx; 
  const double nuq2 = -(0.5*nuSum[3])-0.5*nuSum[2]+0.5*nuSum[1]+0.5*nuSum[0]; 
  const double a2_vx = rdvx2*(-(0.5*nuUSum[3])-0.5*nuUSum[2]+0.5*nuUSum[1]+0.5*nuUSum[0] - nuq2*w[2]); 
  const double b2_vx = -0.5*rdvx2*nuq2*dxv[2]; 
  double f2[3]; 
  f2[0] = -(0.5*f[4])-0.5*f[2]+0.5*f[1]+0.5*f[0]; 
  f2[1] = -(0.5*f[8])-0.5*f[6]+0.5*f[5]+0.5*f[3]; 
  f2[2] = -(0.5*f[11])-0.5*f[10]+0.5*f[9]+0.5*f[7]; 
  double G2[3]; 
  G2[1] = f2[1]*b2_vx+1.7320508075688772*f2[0]*a2_vx; 
  G2[2] = 2.0*f2[2]*b2_vx+2.23606797749979*f2[0]*b2_vx+3.872983346207417*f2[1]*a2_vx; 
  const double nuq3 = 0.5*nuSum[3]+0.5*nuSum[2]+0.5*nuSum[1]+0.5*nuSum[0]; 
  const double a3_vx = rdvx2*(0.5*nuUSum[3]+0.5*nuUSum[2]+0.5*nuUSum[1]+0.5*nuUSum[0] - nuq3*w[2]); 
  const double b3_vx = -0.5*rdvx2*nuq3*dxv[2]; 
  double f3[3]; 
  f3[0] = 0.5*f[4]+0.5*f[2]+0.5*f[1]+0.5*f[0]; 
  f3[1] = 0.5*f[8]+0.5*f[6]+0.5*f[5]+0.5*f[3]; 
  f3[2] = 0.5*f[11]+0.5*f[10]+0.5*f[9]+0.5*f[7]; 
  double G3[3]; 
  G3[1] = f3[1]*b3_vx+1.7320508075688772*f3[0]*a3_vx; 
  G3[2] = 2.0*f3[2]*b3_vx+2.23606797749979*f3[0]*b3_vx+3.872983346207417*f3[1]*a3_vx; 
  out[3] += 0.5*G3[1]+0.5*G2[1]+0.5*G1[1]+0.5*G0[1]; 
  out[5] += 0.5*G3[1]+0.5*G2[1]-0.5*G1[1]-0.5*G0[1]; 
  out[6] += 0.5*G3[1]-0.5*G2[1]+0.5*G1[1]-0.5*G0[1]; 
  out[7] += 0.5*G3[2]+0.5*G2[2]+0.5*G1[2]+0.5*G0[2]; 
  out[8] += 0.5*G3[1]-0.5*G2[1]-0.5*G1[1]+0.5*G0[1]; 
  out[9] += 0.5*G3[2]+0.5*G2[2]-0.5*G1[2]-0.5*G0[2]; 
  out[10] += 0.5*G3[2]-0.5*G2[2]+0.5*G1[2]-0.5*G0[2]; 
  out[11] += 0.5*G3[2]-0.5*G2[2]-0.5*G1[2]+0.5*G0[2]; 

  return fabs(1.25*nuUSum[0]*rdvx2-1.25*nuSum[0]*w[2]*rdvx2); 

} 
