#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 

  double alphaDrag[40]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.0*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.0*nuUSum[1]-2.0*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -(0.5773502691896258*nuSum[0]*dxv[1]*rdvx2); 
  alphaDrag[4] = -(0.5773502691896258*dxv[1]*nuSum[1]*rdvx2); 
  alphaDrag[7] = (2.0*nuUSum[2]-2.0*w[1]*nuSum[2])*rdvx2; 
  alphaDrag[11] = -(0.5773502691896258*dxv[1]*nuSum[2]*rdvx2); 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[20] = (2.0*nuUSum[3]-2.0*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[21] = (2.0*nuUSum[4]-2.0*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[23] = -(0.5773502691896258*nuSum[0]*dxv[2]*rdvy2); 
  alphaDrag[25] = -(0.5773502691896258*nuSum[1]*dxv[2]*rdvy2); 
  alphaDrag[27] = (2.0*nuUSum[5]-2.0*nuSum[2]*w[2])*rdvy2; 
  alphaDrag[33] = -(0.5773502691896258*dxv[2]*nuSum[2]*rdvy2); 

  out[2] += 0.6123724356957944*(alphaDrag[11]*f[11]+alphaDrag[7]*f[7]+alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[13]*alphaDrag[33]+f[7]*alphaDrag[27]+f[5]*alphaDrag[25]+f[3]*alphaDrag[23]+f[1]*alphaDrag[21]+f[0]*alphaDrag[20]); 
  out[4] += 0.5477225575051661*(alphaDrag[4]*f[11]+f[4]*alphaDrag[11]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*(f[5]*alphaDrag[33]+f[1]*alphaDrag[27]+f[13]*alphaDrag[25])+0.6123724356957944*(f[3]*alphaDrag[25]+f[5]*alphaDrag[23])+0.5477225575051661*f[7]*alphaDrag[21]+0.6123724356957944*(f[0]*alphaDrag[21]+f[1]*alphaDrag[20]); 
  out[6] += 0.6123724356957944*(f[17]*alphaDrag[33]+f[11]*alphaDrag[27]+f[10]*alphaDrag[25]+f[6]*alphaDrag[23]+f[4]*alphaDrag[21]+f[2]*alphaDrag[20]+alphaDrag[11]*f[17]+alphaDrag[7]*f[13]+alphaDrag[4]*f[10]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 1.224744871391589*alphaDrag[4]*f[12]+1.369306393762915*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11])+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.369306393762915*(f[7]*alphaDrag[33]+f[13]*alphaDrag[27])+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[25]+1.224744871391589*f[9]*alphaDrag[23]+1.369306393762915*(f[0]*alphaDrag[23]+f[5]*alphaDrag[21]+f[3]*alphaDrag[20]); 
  out[10] += 0.5477225575051661*(f[10]*alphaDrag[33]+f[4]*alphaDrag[27]+f[17]*alphaDrag[25])+0.6123724356957944*(f[6]*alphaDrag[25]+f[10]*alphaDrag[23])+0.5477225575051661*f[11]*alphaDrag[21]+0.6123724356957944*(f[2]*alphaDrag[21]+f[4]*alphaDrag[20])+0.5477225575051661*(alphaDrag[4]*f[17]+alphaDrag[1]*f[13])+f[10]*(0.5477225575051661*alphaDrag[11]+0.6123724356957944*alphaDrag[2])+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.39123039821797573*alphaDrag[11]*f[11]+0.6123724356957944*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11])+0.39123039821797573*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*(alphaDrag[4]*f[4]+alphaDrag[1]*f[1]); 
  out[12] += 1.0954451150103324*alphaDrag[11]*f[12]+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[1]*alphaDrag[11]+alphaDrag[4]*(f[8]+f[7])+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += (0.39123039821797573*f[13]+0.6123724356957944*f[3])*alphaDrag[33]+(0.39123039821797573*f[7]+0.6123724356957944*f[0])*alphaDrag[27]+0.5477225575051661*f[5]*alphaDrag[25]+0.6123724356957944*f[13]*alphaDrag[23]+0.5477225575051661*f[1]*alphaDrag[21]+0.6123724356957944*f[7]*alphaDrag[20]; 
  out[14] += 0.6123724356957944*(f[18]*alphaDrag[25]+f[14]*alphaDrag[23]+f[12]*alphaDrag[21]+f[8]*alphaDrag[20])+1.224744871391589*alphaDrag[4]*f[18]+1.369306393762915*alphaDrag[7]*f[17]+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[11]*f[13]+alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[4]*f[5]+alphaDrag[2]*f[3]); 
  out[15] += 1.0954451150103324*f[15]*alphaDrag[33]+1.224744871391589*(f[1]*alphaDrag[33]+f[5]*alphaDrag[27])+(1.224744871391589*(f[9]+f[7])+1.369306393762915*f[0])*alphaDrag[25]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[23]+1.224744871391589*f[13]*alphaDrag[21]+1.369306393762915*(f[3]*alphaDrag[21]+f[5]*alphaDrag[20]); 
  out[16] += 1.369306393762915*(f[11]*alphaDrag[33]+f[17]*alphaDrag[27])+(1.224744871391589*f[19]+1.369306393762915*f[4])*alphaDrag[25]+1.224744871391589*f[16]*alphaDrag[23]+1.369306393762915*(f[2]*alphaDrag[23]+f[10]*alphaDrag[21]+f[6]*alphaDrag[20])+0.6123724356957944*(alphaDrag[4]*f[19]+alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[17] += (0.39123039821797573*f[17]+0.6123724356957944*f[6])*alphaDrag[33]+(0.39123039821797573*f[11]+0.6123724356957944*f[2])*alphaDrag[27]+0.5477225575051661*f[10]*alphaDrag[25]+0.6123724356957944*f[17]*alphaDrag[23]+0.5477225575051661*f[4]*alphaDrag[21]+0.6123724356957944*f[11]*alphaDrag[20]+(0.39123039821797573*alphaDrag[11]+0.6123724356957944*alphaDrag[2])*f[17]+0.39123039821797573*alphaDrag[7]*f[13]+0.6123724356957944*(alphaDrag[0]*f[13]+f[6]*alphaDrag[11])+0.5477225575051661*alphaDrag[4]*f[10]+0.6123724356957944*f[3]*alphaDrag[7]+0.5477225575051661*alphaDrag[1]*f[5]; 
  out[18] += 0.5477225575051661*(f[18]*alphaDrag[33]+f[12]*alphaDrag[27])+0.6123724356957944*(f[14]*alphaDrag[25]+f[18]*alphaDrag[23]+f[8]*alphaDrag[21]+f[12]*alphaDrag[20])+1.0954451150103324*alphaDrag[11]*f[18]+1.224744871391589*(alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[4]*(f[14]+f[13])+f[5]*alphaDrag[11]+alphaDrag[7]*f[10])+1.369306393762915*(alphaDrag[0]*f[10]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]+f[3]*alphaDrag[4]); 
  out[19] += 1.0954451150103324*f[19]*alphaDrag[33]+1.224744871391589*(f[4]*alphaDrag[33]+f[10]*alphaDrag[27])+(1.224744871391589*(f[16]+f[11])+1.369306393762915*f[2])*alphaDrag[25]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alphaDrag[23]+1.224744871391589*f[17]*alphaDrag[21]+1.369306393762915*(f[6]*alphaDrag[21]+f[10]*alphaDrag[20])+0.5477225575051661*alphaDrag[11]*f[19]+0.6123724356957944*(alphaDrag[2]*f[19]+alphaDrag[4]*f[16])+0.5477225575051661*alphaDrag[7]*f[15]+0.6123724356957944*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]); 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[7])+fabs(0.8838834764831842*alphaDrag[20]-0.9882117688026182*alphaDrag[27]); 

} 
