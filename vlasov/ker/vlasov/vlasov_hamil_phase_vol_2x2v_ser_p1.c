#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  out[1] += (0.75*f[12]*hamil[15]+0.75*f[9]*hamil[14]+0.75*f[8]*hamil[13]+0.75*f[5]*hamil[11]+0.75*f[4]*hamil[10]+0.75*f[2]*hamil[7]+0.75*f[1]*hamil[6]+0.75*f[0]*hamil[3])*dv10*dx10; 
  out[3] += (-(0.75*f[14]*hamil[15])-0.75*f[10]*hamil[13]-0.75*f[9]*hamil[12]-0.75*f[7]*hamil[11]-0.75*f[4]*hamil[8]-0.75*f[3]*hamil[6]-0.75*f[2]*hamil[5]-0.75*f[0]*hamil[1])*dv10*dx10; 
  out[5] += (0.75*f[8]*hamil[15]+0.75*f[4]*hamil[14]+0.75*f[12]*hamil[13]+0.75*f[1]*hamil[11]+0.75*f[9]*hamil[10]+0.75*f[0]*hamil[7]+0.75*f[5]*hamil[6]+0.75*f[2]*hamil[3])*dv10*dx10; 
  out[6] += (0.75*f[14]*hamil[14]-0.75*f[12]*hamil[12]+0.75*f[10]*hamil[10]-0.75*f[8]*hamil[8]+0.75*f[7]*hamil[7]-0.75*f[5]*hamil[5]+0.75*f[3]*hamil[3]-0.75*f[1]*hamil[1])*dv10*dx10; 
  out[7] += (-(0.75*f[10]*hamil[15])-0.75*hamil[13]*f[14]-0.75*f[4]*hamil[12]-0.75*f[3]*hamil[11]-0.75*hamil[8]*f[9]-0.75*hamil[6]*f[7]-0.75*f[0]*hamil[5]-0.75*hamil[1]*f[2])*dv10*dx10; 
  out[8] += (0.75*f[5]*hamil[15]+0.75*f[2]*hamil[14]+0.75*f[1]*hamil[13]+0.75*hamil[11]*f[12]+0.75*f[0]*hamil[10]+0.75*hamil[7]*f[9]+0.75*hamil[6]*f[8]+0.75*hamil[3]*f[4])*dv10*dx10; 
  out[10] += (-(0.75*f[7]*hamil[15])-0.75*hamil[11]*f[14]-0.75*f[3]*hamil[13]-0.75*f[2]*hamil[12]-0.75*hamil[6]*f[10]-0.75*hamil[5]*f[9]-0.75*f[0]*hamil[8]-0.75*hamil[1]*f[4])*dv10*dx10; 
  out[11] += (0.75*f[10]*hamil[14]+0.75*hamil[10]*f[14]-0.75*f[8]*hamil[12]-0.75*hamil[8]*f[12]+0.75*f[3]*hamil[7]+0.75*hamil[3]*f[7]-0.75*f[1]*hamil[5]-0.75*hamil[1]*f[5])*dv10*dx10; 
  out[12] += (0.75*f[1]*hamil[15]+0.75*f[0]*hamil[14]+0.75*f[5]*hamil[13]+0.75*hamil[6]*f[12]+0.75*f[8]*hamil[11]+0.75*f[2]*hamil[10]+0.75*hamil[3]*f[9]+0.75*f[4]*hamil[7])*dv10*dx10; 
  out[13] += (0.75*f[7]*hamil[14]+0.75*hamil[7]*f[14]-0.75*f[5]*hamil[12]-0.75*hamil[5]*f[12]+0.75*f[3]*hamil[10]+0.75*hamil[3]*f[10]-0.75*f[1]*hamil[8]-0.75*hamil[1]*f[8])*dv10*dx10; 
  out[14] += (-(0.75*f[3]*hamil[15])-0.75*hamil[6]*f[14]-0.75*f[7]*hamil[13]-0.75*f[0]*hamil[12]-0.75*f[10]*hamil[11]-0.75*hamil[1]*f[9]-0.75*f[2]*hamil[8]-0.75*f[4]*hamil[5])*dv10*dx10; 
  out[15] += (0.75*f[3]*hamil[14]+0.75*hamil[3]*f[14]-0.75*f[1]*hamil[12]-0.75*hamil[1]*f[12]+0.75*f[7]*hamil[10]+0.75*hamil[7]*f[10]-0.75*f[5]*hamil[8]-0.75*hamil[5]*f[8])*dv10*dx10; 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  out[2] += (0.75*f[11]*hamil[15]+0.75*f[7]*hamil[14]+0.75*f[6]*hamil[13]+0.75*f[5]*hamil[12]+0.75*f[3]*hamil[10]+0.75*f[2]*hamil[9]+0.75*f[1]*hamil[8]+0.75*f[0]*hamil[4])*dv11*dx11; 
  out[4] += (-(0.75*f[13]*hamil[15])-0.75*f[10]*hamil[14]-0.75*f[8]*hamil[12]-0.75*f[6]*hamil[11]-0.75*f[4]*hamil[9]-0.75*f[3]*hamil[7]-0.75*f[1]*hamil[5]-0.75*f[0]*hamil[2])*dv11*dx11; 
  out[5] += (0.75*f[7]*hamil[15]+0.75*f[11]*hamil[14]+0.75*f[3]*hamil[13]+0.75*f[2]*hamil[12]+0.75*f[6]*hamil[10]+0.75*f[5]*hamil[9]+0.75*f[0]*hamil[8]+0.75*f[1]*hamil[4])*dv11*dx11; 
  out[7] += (0.75*f[5]*hamil[15]+0.75*f[2]*hamil[14]+0.75*f[1]*hamil[13]+0.75*f[11]*hamil[12]+0.75*f[0]*hamil[10]+0.75*f[7]*hamil[9]+0.75*f[6]*hamil[8]+0.75*f[3]*hamil[4])*dv11*dx11; 
  out[8] += (-(0.75*f[10]*hamil[15])-0.75*f[13]*hamil[14]-0.75*f[4]*hamil[12]-0.75*f[3]*hamil[11]-0.75*f[8]*hamil[9]-0.75*f[6]*hamil[7]-0.75*f[0]*hamil[5]-0.75*f[1]*hamil[2])*dv11*dx11; 
  out[9] += (0.75*f[13]*hamil[13]-0.75*f[11]*hamil[11]+0.75*f[10]*hamil[10]+0.75*f[8]*hamil[8]-0.75*f[7]*hamil[7]-0.75*f[5]*hamil[5]+0.75*f[4]*hamil[4]-0.75*f[2]*hamil[2])*dv11*dx11; 
  out[10] += (-(0.75*f[8]*hamil[15])-0.75*f[4]*hamil[14]-0.75*hamil[12]*f[13]-0.75*f[1]*hamil[11]-0.75*hamil[9]*f[10]-0.75*f[0]*hamil[7]-0.75*hamil[5]*f[6]-0.75*hamil[2]*f[3])*dv11*dx11; 
  out[11] += (0.75*f[2]*hamil[15]+0.75*f[5]*hamil[14]+0.75*f[0]*hamil[13]+0.75*f[7]*hamil[12]+0.75*hamil[9]*f[11]+0.75*f[1]*hamil[10]+0.75*f[3]*hamil[8]+0.75*hamil[4]*f[6])*dv11*dx11; 
  out[12] += (0.75*f[10]*hamil[13]+0.75*hamil[10]*f[13]-0.75*f[7]*hamil[11]-0.75*hamil[7]*f[11]+0.75*f[4]*hamil[8]+0.75*hamil[4]*f[8]-0.75*f[2]*hamil[5]-0.75*hamil[2]*f[5])*dv11*dx11; 
  out[13] += (-(0.75*f[4]*hamil[15])-0.75*f[8]*hamil[14]-0.75*hamil[9]*f[13]-0.75*f[10]*hamil[12]-0.75*f[0]*hamil[11]-0.75*f[1]*hamil[7]-0.75*hamil[2]*f[6]-0.75*f[3]*hamil[5])*dv11*dx11; 
  out[14] += (0.75*f[8]*hamil[13]+0.75*hamil[8]*f[13]-0.75*f[5]*hamil[11]-0.75*hamil[5]*f[11]+0.75*f[4]*hamil[10]+0.75*hamil[4]*f[10]-0.75*f[2]*hamil[7]-0.75*hamil[2]*f[7])*dv11*dx11; 
  out[15] += (0.75*f[4]*hamil[13]+0.75*hamil[4]*f[13]-0.75*f[2]*hamil[11]-0.75*hamil[2]*f[11]+0.75*f[8]*hamil[10]+0.75*hamil[8]*f[10]-0.75*f[5]*hamil[7]-0.75*hamil[5]*f[7])*dv11*dx11; 

} 
