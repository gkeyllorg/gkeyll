#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_phase_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv0 = dxv[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  double f_dJv0[6] = {0.0}; 
  f_dJv0[0] = (0.24845199749997662*f[4]+0.37267799624996495*f[2]+0.2777777777777778*f[0])/jacob_vx[2]+(0.24845199749997662*f[4]-0.37267799624996495*f[2]+0.2777777777777778*f[0])/jacob_vx[0]+(0.4444444444444444*f[0]-0.49690399499995325*f[4])/jacob_vx[1]; 
  f_dJv0[1] = (0.24845199749997673*f[5]+0.37267799624996495*f[3]+0.2777777777777778*f[1])/jacob_vx[2]+(0.24845199749997673*f[5]-0.37267799624996495*f[3]+0.2777777777777778*f[1])/jacob_vx[0]+(0.4444444444444444*f[1]-0.49690399499995347*f[5])/jacob_vx[1]; 
  f_dJv0[2] = (0.3333333333333333*f[4]+0.5*f[2]+0.37267799624996495*f[0])/jacob_vx[2]+(-(0.3333333333333333*f[4])+0.5*f[2]-0.37267799624996495*f[0])/jacob_vx[0]; 
  f_dJv0[3] = (0.3333333333333334*f[5]+0.5*f[3]+0.37267799624996495*f[1])/jacob_vx[2]+(-(0.3333333333333334*f[5])+0.5*f[3]-0.37267799624996495*f[1])/jacob_vx[0]; 
  f_dJv0[4] = (0.5555555555555556*f[4]-0.49690399499995325*f[0])/jacob_vx[1]+(0.2222222222222222*f[4]+0.3333333333333333*f[2]+0.24845199749997662*f[0])/jacob_vx[2]+(0.2222222222222222*f[4]-0.3333333333333333*f[2]+0.24845199749997662*f[0])/jacob_vx[0]; 
  f_dJv0[5] = (0.5555555555555556*f[5]-0.49690399499995347*f[1])/jacob_vx[1]+(0.2222222222222222*f[5]+0.3333333333333334*f[3]+0.24845199749997673*f[1])/jacob_vx[2]+(0.2222222222222222*f[5]-0.3333333333333334*f[3]+0.24845199749997673*f[1])/jacob_vx[0]; 

  const double wx1 = w[1]; 
  const double dv10_sq = dv10*dv10; 
  
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  out[1] += (2.371708245126284*poisson_tensor_conf_0[0]*f_dJv0[3]*hamil[5]+2.371708245126284*poisson_tensor_conf_0[1]*f_dJv0[2]*hamil[5]+2.371708245126284*poisson_tensor_conf_0[1]*f_dJv0[3]*hamil[4]+2.371708245126284*poisson_tensor_conf_0[0]*f_dJv0[2]*hamil[4]+1.060660171779821*f_dJv0[0]*poisson_tensor_conf_0[1]*hamil[3]+1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[1]*hamil[3]+1.060660171779821*f_dJv0[1]*poisson_tensor_conf_0[1]*hamil[2]+1.060660171779821*f_dJv0[0]*poisson_tensor_conf_0[0]*hamil[2])*dv10*dx10*jacob_cx_inv; 
  out[2] += (-(1.060660171779821*poisson_tensor_conf_0[1]*f_dJv0[5]*hamil[5])-1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[4]*hamil[5]-1.060660171779821*poisson_tensor_conf_0[1]*f_dJv0[3]*hamil[3]-1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[2]*hamil[3]-1.060660171779821*f_dJv0[1]*hamil[1]*poisson_tensor_conf_0[1]-1.060660171779821*f_dJv0[0]*poisson_tensor_conf_0[0]*hamil[1])*dv10*dx10*jacob_cx_inv; 
  out[3] += (1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[5]*hamil[5]+1.060660171779821*poisson_tensor_conf_0[1]*f_dJv0[4]*hamil[5]+2.371708245126284*f_dJv0[0]*poisson_tensor_conf_0[1]*hamil[5]+2.371708245126284*poisson_tensor_conf_0[0]*f_dJv0[1]*hamil[5]+2.1213203435596424*poisson_tensor_conf_0[1]*hamil[4]*f_dJv0[5]+2.1213203435596424*poisson_tensor_conf_0[0]*f_dJv0[4]*hamil[4]+2.371708245126284*f_dJv0[1]*poisson_tensor_conf_0[1]*hamil[4]+2.371708245126284*f_dJv0[0]*poisson_tensor_conf_0[0]*hamil[4]+1.060660171779821*poisson_tensor_conf_0[1]*hamil[2]*f_dJv0[3]+1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[2]*hamil[2]-1.060660171779821*f_dJv0[0]*hamil[1]*poisson_tensor_conf_0[1]-1.060660171779821*poisson_tensor_conf_0[0]*f_dJv0[1]*hamil[1])*dv10*dx10*jacob_cx_inv; 
  out[4] += (-(2.1213203435596424*poisson_tensor_conf_0[1]*f_dJv0[3]*hamil[5])-2.1213203435596424*poisson_tensor_conf_0[0]*f_dJv0[2]*hamil[5]-2.1213203435596424*poisson_tensor_conf_0[1]*hamil[3]*f_dJv0[5]-2.1213203435596424*poisson_tensor_conf_0[0]*hamil[3]*f_dJv0[4]-2.371708245126284*f_dJv0[1]*poisson_tensor_conf_0[1]*hamil[3]-2.371708245126284*f_dJv0[0]*poisson_tensor_conf_0[0]*hamil[3]-2.371708245126284*hamil[1]*poisson_tensor_conf_0[1]*f_dJv0[3]-2.371708245126284*poisson_tensor_conf_0[0]*hamil[1]*f_dJv0[2])*dv10*dx10*jacob_cx_inv; 
  out[5] += (-(1.060660171779821*poisson_tensor_conf_0[0]*hamil[3]*f_dJv0[5])+1.060660171779821*poisson_tensor_conf_0[1]*hamil[2]*f_dJv0[5]+2.1213203435596424*poisson_tensor_conf_0[1]*f_dJv0[3]*hamil[4]+2.1213203435596424*poisson_tensor_conf_0[0]*f_dJv0[2]*hamil[4]-1.060660171779821*poisson_tensor_conf_0[1]*hamil[3]*f_dJv0[4]+1.060660171779821*poisson_tensor_conf_0[0]*hamil[2]*f_dJv0[4]-2.371708245126284*f_dJv0[0]*poisson_tensor_conf_0[1]*hamil[3]-2.371708245126284*poisson_tensor_conf_0[0]*f_dJv0[1]*hamil[3]-2.371708245126284*poisson_tensor_conf_0[0]*hamil[1]*f_dJv0[3]-2.371708245126284*hamil[1]*poisson_tensor_conf_0[1]*f_dJv0[2])*dv10*dx10*jacob_cx_inv; 

} 
