#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_vx_inv = 1.0/jacob_vel[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel[2]; 
  const double *vmap_vx = &vmap[0]; 
  const double *vmap_vy = &vmap[4]; 
  const double wx1 = w[1]; 
  const double dv10_sq = dv10*dv10; 
  const double wx2 = w[2]; 
  const double dv11_sq = dv11*dv11; 
  
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  out[1] += (1.060660171779821*f[1]*hamil[1]*poisson_tensor_conf_0[1]+1.060660171779821*f[0]*poisson_tensor_conf_0[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[4] += (1.060660171779821*hamil[1]*poisson_tensor_conf_0[1]*f[4]+1.060660171779821*poisson_tensor_conf_0[0]*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (1.060660171779821*hamil[1]*poisson_tensor_conf_0[1]*f[5]+1.060660171779821*poisson_tensor_conf_0[0]*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[7] += (1.060660171779821*hamil[1]*poisson_tensor_conf_0[1]*f[7]+1.060660171779821*poisson_tensor_conf_0[0]*hamil[1]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  out[1] += (1.060660171779821*f[1]*poisson_tensor_conf_1[1]*hamil[2]+1.060660171779821*f[0]*poisson_tensor_conf_1[0]*hamil[2])*dv11*dx10*jacob_cx_inv*jacob_vy_inv; 
  out[4] += (1.060660171779821*poisson_tensor_conf_1[1]*hamil[2]*f[4]+1.060660171779821*poisson_tensor_conf_1[0]*f[2]*hamil[2])*dv11*dx10*jacob_cx_inv*jacob_vy_inv; 
  out[5] += (1.060660171779821*poisson_tensor_conf_1[1]*hamil[2]*f[5]+1.060660171779821*poisson_tensor_conf_1[0]*hamil[2]*f[3])*dv11*dx10*jacob_cx_inv*jacob_vy_inv; 
  out[7] += (1.060660171779821*poisson_tensor_conf_1[1]*hamil[2]*f[7]+1.060660171779821*poisson_tensor_conf_1[0]*hamil[2]*f[6])*dv11*dx10*jacob_cx_inv*jacob_vy_inv; 

  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[8]; 
  out[2] += (0.75*poisson_tensor_conf_4[1]*vmap_vx[1]*hamil[2]*f[4]+0.75*poisson_tensor_conf_4[0]*vmap_vx[1]*f[2]*hamil[2]+0.75*vmap_vx[0]*f[1]*poisson_tensor_conf_4[1]*hamil[2]+0.75*f[0]*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[3] += (-(0.75*hamil[1]*poisson_tensor_conf_4[1]*vmap_vx[1]*f[4])-0.75*poisson_tensor_conf_4[0]*hamil[1]*vmap_vx[1]*f[2]-0.75*vmap_vx[0]*f[1]*hamil[1]*poisson_tensor_conf_4[1]-0.75*f[0]*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[4] += (0.75*poisson_tensor_conf_4[0]*vmap_vx[1]*hamil[2]*f[4]+0.75*poisson_tensor_conf_4[1]*vmap_vx[1]*f[2]*hamil[2]+0.75*f[0]*vmap_vx[0]*poisson_tensor_conf_4[1]*hamil[2]+0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*f[1]*hamil[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[5] += (-(0.75*poisson_tensor_conf_4[0]*hamil[1]*vmap_vx[1]*f[4])-0.75*hamil[1]*poisson_tensor_conf_4[1]*vmap_vx[1]*f[2]-0.75*f[0]*vmap_vx[0]*hamil[1]*poisson_tensor_conf_4[1]-0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*f[1]*hamil[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[6] += (0.75*poisson_tensor_conf_4[1]*vmap_vx[1]*hamil[2]*f[7]+0.75*poisson_tensor_conf_4[0]*vmap_vx[1]*hamil[2]*f[6]+0.75*vmap_vx[0]*poisson_tensor_conf_4[1]*hamil[2]*f[5]-0.75*vmap_vx[0]*hamil[1]*poisson_tensor_conf_4[1]*f[4]+0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[2]*f[3]-0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[1]*f[2]-0.75*f[1]*hamil[1]*poisson_tensor_conf_4[1]*vmap_vx[1]-0.75*f[0]*poisson_tensor_conf_4[0]*hamil[1]*vmap_vx[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[7] += (0.75*poisson_tensor_conf_4[0]*vmap_vx[1]*hamil[2]*f[7]+0.75*poisson_tensor_conf_4[1]*vmap_vx[1]*hamil[2]*f[6]+0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[2]*f[5]-0.75*poisson_tensor_conf_4[0]*vmap_vx[0]*hamil[1]*f[4]+0.75*vmap_vx[0]*poisson_tensor_conf_4[1]*hamil[2]*f[3]-0.75*vmap_vx[0]*hamil[1]*poisson_tensor_conf_4[1]*f[2]-0.75*f[0]*hamil[1]*poisson_tensor_conf_4[1]*vmap_vx[1]-0.75*poisson_tensor_conf_4[0]*f[1]*hamil[1]*vmap_vx[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 

  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[10]; 
  out[2] += (0.75*poisson_tensor_conf_5[1]*vmap_vy[1]*hamil[2]*f[5]+0.75*poisson_tensor_conf_5[0]*vmap_vy[1]*hamil[2]*f[3]+0.75*vmap_vy[0]*f[1]*poisson_tensor_conf_5[1]*hamil[2]+0.75*f[0]*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[3] += (-(0.75*hamil[1]*poisson_tensor_conf_5[1]*vmap_vy[1]*f[5])-0.75*poisson_tensor_conf_5[0]*hamil[1]*vmap_vy[1]*f[3]-0.75*vmap_vy[0]*f[1]*hamil[1]*poisson_tensor_conf_5[1]-0.75*f[0]*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[4] += (0.75*poisson_tensor_conf_5[0]*vmap_vy[1]*hamil[2]*f[5]+0.75*poisson_tensor_conf_5[1]*vmap_vy[1]*hamil[2]*f[3]+0.75*f[0]*vmap_vy[0]*poisson_tensor_conf_5[1]*hamil[2]+0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*f[1]*hamil[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[5] += (-(0.75*poisson_tensor_conf_5[0]*hamil[1]*vmap_vy[1]*f[5])-0.75*hamil[1]*poisson_tensor_conf_5[1]*vmap_vy[1]*f[3]-0.75*f[0]*vmap_vy[0]*hamil[1]*poisson_tensor_conf_5[1]-0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*f[1]*hamil[1])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[6] += (-(0.75*hamil[1]*poisson_tensor_conf_5[1]*vmap_vy[1]*f[7])-0.75*poisson_tensor_conf_5[0]*hamil[1]*vmap_vy[1]*f[6]+0.75*vmap_vy[0]*poisson_tensor_conf_5[1]*hamil[2]*f[5]-0.75*vmap_vy[0]*hamil[1]*poisson_tensor_conf_5[1]*f[4]+0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[2]*f[3]+0.75*f[1]*poisson_tensor_conf_5[1]*vmap_vy[1]*hamil[2]+0.75*f[0]*poisson_tensor_conf_5[0]*vmap_vy[1]*hamil[2]-0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[1]*f[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 
  out[7] += (-(0.75*poisson_tensor_conf_5[0]*hamil[1]*vmap_vy[1]*f[7])-0.75*hamil[1]*poisson_tensor_conf_5[1]*vmap_vy[1]*f[6]+0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[2]*f[5]-0.75*poisson_tensor_conf_5[0]*vmap_vy[0]*hamil[1]*f[4]+0.75*vmap_vy[0]*poisson_tensor_conf_5[1]*hamil[2]*f[3]+0.75*f[0]*poisson_tensor_conf_5[1]*vmap_vy[1]*hamil[2]+0.75*poisson_tensor_conf_5[0]*f[1]*vmap_vy[1]*hamil[2]-0.75*vmap_vy[0]*hamil[1]*poisson_tensor_conf_5[1]*f[2])*dv10*dv11*jacob_vx_inv*jacob_vy_inv; 

} 
