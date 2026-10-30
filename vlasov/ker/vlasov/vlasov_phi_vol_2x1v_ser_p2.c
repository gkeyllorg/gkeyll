#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x1v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[3] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[7]*f[8])-3.3541019662496843*f[4]*phi[6]-3.3541019662496847*f[1]*phi[4]-1.5*f[2]*phi[3]-1.5*f[0]*phi[1]); 
  out[5] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[7]*f[12])+phi[6]*(-(3.0*f[11])-3.3541019662496843*f[2])+phi[4]*(-(3.0*f[7])-3.3541019662496847*f[0])-1.5*phi[3]*f[4]-1.5*f[1]*phi[1]); 
  out[6] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[6]*(-(3.0*f[12])-3.3541019662496843*f[1])+phi[3]*(-(1.3416407864998738*f[8])-1.5*f[0])-1.3416407864998738*f[2]*phi[7]-3.3541019662496847*f[4]*phi[4]-1.5*phi[1]*f[2]); 
  out[9] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(3.3541019662496847*phi[7]*f[14])-7.500000000000001*phi[6]*f[10]-3.3541019662496847*phi[3]*f[6]-7.5*phi[4]*f[5]-3.3541019662496847*phi[1]*f[3]); 
  out[10] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[3]*(-(1.3416407864998738*f[12])-1.5*f[1])+phi[4]*(-(3.0*f[11])-3.3541019662496847*f[2])+phi[6]*(-(3.0*f[8])-3.0*f[7]-3.3541019662496843*f[0])-1.3416407864998738*f[4]*phi[7]-1.5*phi[1]*f[4]); 
  out[13] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[3]*f[11])-1.5*phi[1]*f[7]-3.0*f[4]*phi[6]-3.0*f[1]*phi[4]); 
  out[14] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(3.3541019662496847*phi[4]*f[12])-1.5*phi[1]*f[8]+phi[7]*(-(0.9583148474999099*f[8])-1.5*f[0])-3.0*f[4]*phi[6]-1.3416407864998738*f[2]*phi[3]); 
  out[15] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(3.3541019662496847*phi[7]*f[18])+phi[6]*(-(6.708203932499369*f[17])-7.5*f[6])+phi[4]*(-(6.708203932499369*f[13])-7.500000000000001*f[3])-3.3541019662496843*phi[3]*f[10]-3.3541019662496843*phi[1]*f[5]); 
  out[16] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[6]*(-(6.708203932499369*f[18])-7.5*f[5])+phi[3]*(-(3.0*f[14])-3.3541019662496843*f[3])-7.500000000000001*phi[4]*f[10]-3.0*f[6]*phi[7]-3.3541019662496843*phi[1]*f[6]); 
  out[17] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[6]*(-(2.6832815729997477*f[12])-3.0*f[1])-1.3416407864998738*phi[7]*f[11]-1.5*phi[1]*f[11]-1.5*phi[3]*f[7]-3.0*f[4]*phi[4]); 
  out[18] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[1]*f[12])+phi[7]*(-(0.9583148474999099*f[12])-1.5*f[1])+phi[6]*(-(2.6832815729997477*f[11])-3.0*f[2])-3.3541019662496847*phi[4]*f[8]-1.3416407864998738*phi[3]*f[4]); 
  out[19] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[3]*(-(3.0*f[18])-3.3541019662496847*f[5])+phi[4]*(-(6.708203932499369*f[17])-7.5*f[6])+phi[6]*(-(6.708203932499369*f[14])-6.708203932499369*f[13]-7.500000000000001*f[3])-3.0*phi[7]*f[10]-3.3541019662496847*phi[1]*f[10]); 

} 
