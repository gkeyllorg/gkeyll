#include <gkyl_boltz_potential_kernels.h>

GKYL_CU_DH void boltz_potential_boundary_lower_2x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[2] = area[0];
  out[4] = 0.8660254037844386*m0_elc[3]-0.8660254037844386*m0_elc[2]-0.5*m0_elc[1]+0.5*m0_elc[0];
  out[6] = wall ? 0.8660254037844386*wall[3]-0.8660254037844386*wall[2]-0.5*wall[1]+0.5*wall[0] : 0.0;
  out[1] = flux[1];
  out[3] = area[1];
  out[5] = 0.5*(m0_elc[1]+m0_elc[0])-0.8660254037844386*(m0_elc[3]+m0_elc[2]);
  out[7] = wall ? 0.5*(wall[1]+wall[0])-0.8660254037844386*(wall[3]+wall[2]) : 0.0;
}

GKYL_CU_DH void boltz_potential_boundary_upper_2x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[2] = area[0];
  out[4] = -(0.8660254037844386*m0_elc[3])+0.8660254037844386*m0_elc[2]-0.5*m0_elc[1]+0.5*m0_elc[0];
  out[6] = wall ? -(0.8660254037844386*wall[3])+0.8660254037844386*wall[2]-0.5*wall[1]+0.5*wall[0] : 0.0;
  out[1] = flux[1];
  out[3] = area[1];
  out[5] = 0.8660254037844386*(m0_elc[3]+m0_elc[2])+0.5*(m0_elc[1]+m0_elc[0]);
  out[7] = wall ? 0.8660254037844386*(wall[3]+wall[2])+0.5*(wall[1]+wall[0]) : 0.0;
}

GKYL_CU_DH void boltz_potential_phi_2x_ser_p1(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi) 
{ 
  // q_e: signed electron charge.
  // T_e: electron temperature.
  // m0Elc: electron number density without the configuration Jacobian.
  // nphi_ref: reference density and potential in two volume-modal blocks.
  // phi: electrostatic potential in domain volume.

  double phi_qp[4];
  double m0ElcRef_curr;
  m0ElcRef_curr = 0.5*nphi_ref[3]-0.5*nphi_ref[2]-0.5*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[0] = -((1.0*T_e*log(fmax(0.5*m0Elc[3]-0.5*m0Elc[2]-0.5*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.5*nphi_ref[7]-0.5*nphi_ref[6]-0.5*nphi_ref[5]+0.5*nphi_ref[4];
  } else {
    phi_qp[0] = 0.0;
  }
  m0ElcRef_curr = -(0.5*nphi_ref[3])+0.5*nphi_ref[2]-0.5*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[1] = -((1.0*T_e*log(fmax(-(0.5*m0Elc[3])+0.5*m0Elc[2]-0.5*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.5*nphi_ref[7]+0.5*nphi_ref[6]-0.5*nphi_ref[5]+0.5*nphi_ref[4];
  } else {
    phi_qp[1] = 0.0;
  }
  m0ElcRef_curr = -(0.5*nphi_ref[3])-0.5*nphi_ref[2]+0.5*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[2] = -((1.0*T_e*log(fmax(-(0.5*m0Elc[3])-0.5*m0Elc[2]+0.5*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.5*nphi_ref[7]-0.5*nphi_ref[6]+0.5*nphi_ref[5]+0.5*nphi_ref[4];
  } else {
    phi_qp[2] = 0.0;
  }
  m0ElcRef_curr = 0.5*nphi_ref[3]+0.5*nphi_ref[2]+0.5*nphi_ref[1]+0.5*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[3] = -((1.0*T_e*log(fmax(0.5*m0Elc[3]+0.5*m0Elc[2]+0.5*m0Elc[1]+0.5*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.5*nphi_ref[7]+0.5*nphi_ref[6]+0.5*nphi_ref[5]+0.5*nphi_ref[4];
  } else {
    phi_qp[3] = 0.0;
  }

  phi[0] = 0.5*(phi_qp[3]+phi_qp[2]+phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.5*(phi_qp[3]+phi_qp[2])-0.5*(phi_qp[1]+phi_qp[0]); 
  phi[2] = 0.5*phi_qp[3]-0.5*phi_qp[2]+0.5*phi_qp[1]-0.5*phi_qp[0]; 
  phi[3] = 0.5*phi_qp[3]-0.5*(phi_qp[2]+phi_qp[1])+0.5*phi_qp[0]; 
}

