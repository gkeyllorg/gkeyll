#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[1] += (2.371708245126284*f[2]*hamil[7]+1.060660171779821*f[0]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += (2.1213203435596424*hamil[7]*f[12]+2.371708245126284*f[0]*hamil[7]+1.060660171779821*hamil[1]*f[2])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += (2.371708245126284*f[7]*hamil[7]+1.060660171779821*hamil[1]*f[3])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[8] += (2.371708245126284*hamil[7]*f[9]+1.060660171779821*hamil[1]*f[4])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += (5.303300858899105*f[5]*hamil[7]+2.371708245126284*f[1]*hamil[1])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += (2.1213203435596424*hamil[7]*f[22]+2.371708245126284*f[3]*hamil[7]+1.060660171779821*hamil[1]*f[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[16] += (2.1213203435596424*hamil[7]*f[26]+1.060660171779821*hamil[1]*f[9]+2.371708245126284*f[4]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[17] += (2.371708245126284*hamil[7]*f[18]+1.060660171779821*hamil[1]*f[10])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[19] += (4.743416490252569*hamil[7]*f[20]+5.303300858899106*f[1]*hamil[7]+2.371708245126284*hamil[1]*f[5])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[20] += (1.060660171779821*hamil[1]*f[12]+2.1213203435596424*f[2]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[21] += (5.303300858899106*hamil[7]*f[15]+2.371708245126284*hamil[1]*f[6])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[23] += (2.371708245126284*hamil[7]*f[24]+1.060660171779821*hamil[1]*f[13])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[25] += (5.303300858899106*hamil[7]*f[16]+2.371708245126284*hamil[1]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[28] += (2.371708245126284*hamil[7]*f[29]+1.060660171779821*hamil[1]*f[14])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[31] += (2.1213203435596424*hamil[7]*f[38]+1.060660171779821*hamil[1]*f[18]+2.371708245126284*hamil[7]*f[10])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[32] += (4.743416490252569*hamil[7]*f[33]+2.371708245126284*hamil[1]*f[15]+5.303300858899105*f[6]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[33] += (1.060660171779821*hamil[1]*f[22]+2.1213203435596424*f[7]*hamil[7])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[34] += (1.060660171779821*hamil[1]*f[24]+2.371708245126284*hamil[7]*f[13])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[35] += (4.743416490252569*hamil[7]*f[36]+2.371708245126284*hamil[1]*f[16]+5.303300858899105*hamil[7]*f[8])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[36] += (1.060660171779821*hamil[1]*f[26]+2.1213203435596424*hamil[7]*f[9])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[37] += (5.303300858899105*hamil[7]*f[31]+2.371708245126284*hamil[1]*f[17])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[39] += (2.371708245126284*hamil[7]*f[40]+1.060660171779821*hamil[1]*f[27])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[41] += (1.060660171779821*hamil[1]*f[29]+2.371708245126284*hamil[7]*f[14])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[42] += (2.371708245126284*hamil[7]*f[43]+1.060660171779821*hamil[1]*f[30])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[44] += (4.743416490252569*hamil[7]*f[45]+2.371708245126284*hamil[1]*f[31]+5.303300858899106*hamil[7]*f[17])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[45] += (1.060660171779821*hamil[1]*f[38]+2.1213203435596424*hamil[7]*f[18])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[46] += (1.060660171779821*hamil[1]*f[40]+2.371708245126284*hamil[7]*f[27])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[47] += (1.060660171779821*hamil[1]*f[43]+2.371708245126284*hamil[7]*f[30])*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

} 
