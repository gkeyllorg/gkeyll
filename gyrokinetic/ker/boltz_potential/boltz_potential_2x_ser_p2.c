#include <gkyl_boltz_potential_kernels.h>

GKYL_CU_DH void boltz_potential_boundary_lower_2x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[3] = area[0];
  out[6] = -(1.4999999999999996*m0_elc[7])-0.7745966692414833*m0_elc[6]+1.118033988749895*m0_elc[5]+0.4472135954999579*m0_elc[4]+1.1618950038622249*m0_elc[3]-0.8660254037844386*m0_elc[2]-0.6708203932499369*m0_elc[1]+0.5*m0_elc[0];
  out[9] = wall ? -(1.4999999999999996*wall[7])-0.7745966692414833*wall[6]+1.118033988749895*wall[5]+0.4472135954999579*wall[4]+1.1618950038622249*wall[3]-0.8660254037844386*wall[2]-0.6708203932499369*wall[1]+0.5*wall[0] : 0.0;
  out[1] = flux[1];
  out[4] = area[1];
  out[7] = 0.9682458365518543*m0_elc[6]+1.118033988749895*m0_elc[5]-0.5590169943749475*m0_elc[4]-0.8660254037844386*m0_elc[2]+0.5*m0_elc[0];
  out[10] = wall ? 0.9682458365518543*wall[6]+1.118033988749895*wall[5]-0.5590169943749475*wall[4]-0.8660254037844386*wall[2]+0.5*wall[0] : 0.0;
  out[2] = flux[2];
  out[5] = area[2];
  out[8] = 1.4999999999999996*m0_elc[7]-0.7745966692414833*m0_elc[6]+1.118033988749895*m0_elc[5]+0.4472135954999579*m0_elc[4]-1.1618950038622249*m0_elc[3]-0.8660254037844386*m0_elc[2]+0.6708203932499369*m0_elc[1]+0.5*m0_elc[0];
  out[11] = wall ? 1.4999999999999996*wall[7]-0.7745966692414833*wall[6]+1.118033988749895*wall[5]+0.4472135954999579*wall[4]-1.1618950038622249*wall[3]-0.8660254037844386*wall[2]+0.6708203932499369*wall[1]+0.5*wall[0] : 0.0;
}

GKYL_CU_DH void boltz_potential_boundary_upper_2x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[3] = area[0];
  out[6] = -(1.4999999999999996*m0_elc[7])+0.7745966692414833*m0_elc[6]+1.118033988749895*m0_elc[5]+0.4472135954999579*m0_elc[4]-1.1618950038622249*m0_elc[3]+0.8660254037844386*m0_elc[2]-0.6708203932499369*m0_elc[1]+0.5*m0_elc[0];
  out[9] = wall ? -(1.4999999999999996*wall[7])+0.7745966692414833*wall[6]+1.118033988749895*wall[5]+0.4472135954999579*wall[4]-1.1618950038622249*wall[3]+0.8660254037844386*wall[2]-0.6708203932499369*wall[1]+0.5*wall[0] : 0.0;
  out[1] = flux[1];
  out[4] = area[1];
  out[7] = -(0.9682458365518543*m0_elc[6])+1.118033988749895*m0_elc[5]-0.5590169943749475*m0_elc[4]+0.8660254037844386*m0_elc[2]+0.5*m0_elc[0];
  out[10] = wall ? -(0.9682458365518543*wall[6])+1.118033988749895*wall[5]-0.5590169943749475*wall[4]+0.8660254037844386*wall[2]+0.5*wall[0] : 0.0;
  out[2] = flux[2];
  out[5] = area[2];
  out[8] = 1.4999999999999996*m0_elc[7]+0.7745966692414833*m0_elc[6]+1.118033988749895*m0_elc[5]+0.4472135954999579*m0_elc[4]+1.1618950038622249*m0_elc[3]+0.8660254037844386*m0_elc[2]+0.6708203932499369*m0_elc[1]+0.5*m0_elc[0];
  out[11] = wall ? 1.4999999999999996*wall[7]+0.7745966692414833*wall[6]+1.118033988749895*wall[5]+0.4472135954999579*wall[4]+1.1618950038622249*wall[3]+0.8660254037844386*wall[2]+0.6708203932499369*wall[1]+0.5*wall[0] : 0.0;
}

GKYL_CU_DH void boltz_potential_phi_2x_ser_p2(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi) 
{ 
  // q_e: signed electron charge.
  // T_e: electron temperature.
  // m0Elc: electron number density without the configuration Jacobian.
  // nphi_ref: reference density and potential in two volume-modal blocks.
  // phi: electrostatic potential in domain volume.

  double phi_qp[9];
  double m0ElcRef_curr;
  m0ElcRef_curr = -(0.5999999999999994*nphi_ref[7])-0.5999999999999999*nphi_ref[6]+0.4472135954999579*nphi_ref[5]+0.4472135954999579*nphi_ref[4]+0.9*nphi_ref[3]-0.6708203932499369*nphi_ref[2]-0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[0] = -((1.0*T_e*log(fmax(-(0.5999999999999994*m0Elc[7])-0.5999999999999999*m0Elc[6]+0.4472135954999579*m0Elc[5]+0.4472135954999579*m0Elc[4]+0.9*m0Elc[3]-0.6708203932499369*m0Elc[2]-0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.5999999999999994*nphi_ref[15]-0.5999999999999999*nphi_ref[14]+0.4472135954999579*nphi_ref[13]+0.4472135954999579*nphi_ref[12]+0.9*nphi_ref[11]-0.6708203932499369*nphi_ref[10]-0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[0] = 0.0;
  }
  m0ElcRef_curr = 0.75*nphi_ref[7]-0.5590169943749475*nphi_ref[5]+0.4472135954999579*nphi_ref[4]-0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[1] = -((1.0*T_e*log(fmax(0.75*m0Elc[7]-0.5590169943749475*m0Elc[5]+0.4472135954999579*m0Elc[4]-0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.75*nphi_ref[15]-0.5590169943749475*nphi_ref[13]+0.4472135954999579*nphi_ref[12]-0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[1] = 0.0;
  }
  m0ElcRef_curr = -(0.5999999999999994*nphi_ref[7])+0.5999999999999999*nphi_ref[6]+0.4472135954999579*nphi_ref[5]+0.4472135954999579*nphi_ref[4]-0.9*nphi_ref[3]+0.6708203932499369*nphi_ref[2]-0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[2] = -((1.0*T_e*log(fmax(-(0.5999999999999994*m0Elc[7])+0.5999999999999999*m0Elc[6]+0.4472135954999579*m0Elc[5]+0.4472135954999579*m0Elc[4]-0.9*m0Elc[3]+0.6708203932499369*m0Elc[2]-0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.5999999999999994*nphi_ref[15]+0.5999999999999999*nphi_ref[14]+0.4472135954999579*nphi_ref[13]+0.4472135954999579*nphi_ref[12]-0.9*nphi_ref[11]+0.6708203932499369*nphi_ref[10]-0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[2] = 0.0;
  }
  m0ElcRef_curr = 0.75*nphi_ref[6]+0.4472135954999579*nphi_ref[5]-0.5590169943749475*nphi_ref[4]-0.6708203932499369*nphi_ref[2]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[3] = -((1.0*T_e*log(fmax(0.75*m0Elc[6]+0.4472135954999579*m0Elc[5]-0.5590169943749475*m0Elc[4]-0.6708203932499369*m0Elc[2]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.75*nphi_ref[14]+0.4472135954999579*nphi_ref[13]-0.5590169943749475*nphi_ref[12]-0.6708203932499369*nphi_ref[10]+0.5*nphi_ref[8];
  } else {
    phi_qp[3] = 0.0;
  }
  m0ElcRef_curr = -(0.5590169943749475*nphi_ref[5])-0.5590169943749475*nphi_ref[4]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[4] = -((1.0*T_e*log(fmax(-(0.5590169943749475*m0Elc[5])-0.5590169943749475*m0Elc[4]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.5590169943749475*nphi_ref[13]-0.5590169943749475*nphi_ref[12]+0.5*nphi_ref[8];
  } else {
    phi_qp[4] = 0.0;
  }
  m0ElcRef_curr = -(0.75*nphi_ref[6])+0.4472135954999579*nphi_ref[5]-0.5590169943749475*nphi_ref[4]+0.6708203932499369*nphi_ref[2]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[5] = -((1.0*T_e*log(fmax(-(0.75*m0Elc[6])+0.4472135954999579*m0Elc[5]-0.5590169943749475*m0Elc[4]+0.6708203932499369*m0Elc[2]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.75*nphi_ref[14]+0.4472135954999579*nphi_ref[13]-0.5590169943749475*nphi_ref[12]+0.6708203932499369*nphi_ref[10]+0.5*nphi_ref[8];
  } else {
    phi_qp[5] = 0.0;
  }
  m0ElcRef_curr = 0.5999999999999994*nphi_ref[7]-0.5999999999999999*nphi_ref[6]+0.4472135954999579*nphi_ref[5]+0.4472135954999579*nphi_ref[4]-0.9*nphi_ref[3]-0.6708203932499369*nphi_ref[2]+0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[6] = -((1.0*T_e*log(fmax(0.5999999999999994*m0Elc[7]-0.5999999999999999*m0Elc[6]+0.4472135954999579*m0Elc[5]+0.4472135954999579*m0Elc[4]-0.9*m0Elc[3]-0.6708203932499369*m0Elc[2]+0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.5999999999999994*nphi_ref[15]-0.5999999999999999*nphi_ref[14]+0.4472135954999579*nphi_ref[13]+0.4472135954999579*nphi_ref[12]-0.9*nphi_ref[11]-0.6708203932499369*nphi_ref[10]+0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[6] = 0.0;
  }
  m0ElcRef_curr = -(0.75*nphi_ref[7])-0.5590169943749475*nphi_ref[5]+0.4472135954999579*nphi_ref[4]+0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[7] = -((1.0*T_e*log(fmax(-(0.75*m0Elc[7])-0.5590169943749475*m0Elc[5]+0.4472135954999579*m0Elc[4]+0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.75*nphi_ref[15]-0.5590169943749475*nphi_ref[13]+0.4472135954999579*nphi_ref[12]+0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[7] = 0.0;
  }
  m0ElcRef_curr = 0.5999999999999994*nphi_ref[7]+0.5999999999999999*nphi_ref[6]+0.4472135954999579*nphi_ref[5]+0.4472135954999579*nphi_ref[4]+0.9*nphi_ref[3]+0.6708203932499369*nphi_ref[2]+0.6708203932499369*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[8] = -((1.0*T_e*log(fmax(0.5999999999999994*m0Elc[7]+0.5999999999999999*m0Elc[6]+0.4472135954999579*m0Elc[5]+0.4472135954999579*m0Elc[4]+0.9*m0Elc[3]+0.6708203932499369*m0Elc[2]+0.6708203932499369*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.5999999999999994*nphi_ref[15]+0.5999999999999999*nphi_ref[14]+0.4472135954999579*nphi_ref[13]+0.4472135954999579*nphi_ref[12]+0.9*nphi_ref[11]+0.6708203932499369*nphi_ref[10]+0.6708203932499369*nphi_ref[9]+0.5*nphi_ref[8];
  } else {
    phi_qp[8] = 0.0;
  }

  phi[0] = 0.154320987654321*phi_qp[8]+0.24691358024691357*phi_qp[7]+0.154320987654321*phi_qp[6]+0.24691358024691357*phi_qp[5]+0.3950617283950617*phi_qp[4]+0.24691358024691357*phi_qp[3]+0.154320987654321*phi_qp[2]+0.24691358024691357*phi_qp[1]+0.154320987654321*phi_qp[0]; 
  phi[1] = 0.20704333124998056*phi_qp[8]+0.33126932999996883*phi_qp[7]+0.20704333124998056*phi_qp[6]-0.20704333124998056*phi_qp[2]-0.33126932999996883*phi_qp[1]-0.20704333124998056*phi_qp[0]; 
  phi[2] = 0.20704333124998056*phi_qp[8]-0.20704333124998056*phi_qp[6]+0.33126932999996883*phi_qp[5]-0.33126932999996883*phi_qp[3]+0.20704333124998056*phi_qp[2]-0.20704333124998056*phi_qp[0]; 
  phi[3] = 0.27777777777777785*phi_qp[8]-0.27777777777777785*(phi_qp[6]+phi_qp[2])+0.27777777777777785*phi_qp[0]; 
  phi[4] = 0.13802888749998704*phi_qp[8]+0.22084621999997922*phi_qp[7]+0.13802888749998704*phi_qp[6]-0.276057774999974*phi_qp[5]-0.44169243999995844*phi_qp[4]-0.276057774999974*phi_qp[3]+0.13802888749998704*phi_qp[2]+0.22084621999997922*phi_qp[1]+0.13802888749998704*phi_qp[0]; 
  phi[5] = 0.13802888749998704*phi_qp[8]-0.276057774999974*phi_qp[7]+0.13802888749998704*phi_qp[6]+0.22084621999997922*phi_qp[5]-0.44169243999995844*phi_qp[4]+0.22084621999997922*phi_qp[3]+0.13802888749998704*phi_qp[2]-0.276057774999974*phi_qp[1]+0.13802888749998704*phi_qp[0]; 
  phi[6] = 0.18518518518518517*phi_qp[8]-0.18518518518518517*phi_qp[6]-0.37037037037037035*phi_qp[5]+0.37037037037037035*phi_qp[3]+0.18518518518518517*phi_qp[2]-0.18518518518518517*phi_qp[0]; 
  phi[7] = 0.18518518518518506*phi_qp[8]-0.37037037037037035*phi_qp[7]+0.18518518518518506*phi_qp[6]-0.18518518518518506*phi_qp[2]+0.37037037037037035*phi_qp[1]-0.18518518518518506*phi_qp[0]; 
}

