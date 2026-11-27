#include <gkyl_boltz_potential_kernels.h>

GKYL_CU_DH void boltz_potential_boundary_lower_1x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[1] = area[0];
  out[2] = 1.5811388300841895*m0_elc[2]-1.224744871391589*m0_elc[1]+0.7071067811865475*m0_elc[0];
  out[3] = wall ? 1.5811388300841895*wall[2]-1.224744871391589*wall[1]+0.7071067811865475*wall[0] : 0.0;
}

GKYL_CU_DH void boltz_potential_boundary_upper_1x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out)
{
  out[0] = flux[0];
  out[1] = area[0];
  out[2] = 1.5811388300841895*m0_elc[2]+1.224744871391589*m0_elc[1]+0.7071067811865475*m0_elc[0];
  out[3] = wall ? 1.5811388300841895*wall[2]+1.224744871391589*wall[1]+0.7071067811865475*wall[0] : 0.0;
}

GKYL_CU_DH void boltz_potential_phi_1x_ser_p2(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi) 
{ 
  // q_e: signed electron charge.
  // T_e: electron temperature.
  // m0Elc: electron number density without the configuration Jacobian.
  // nphi_ref: reference density and potential in two volume-modal blocks.
  // phi: electrostatic potential in domain volume.

  double phi_qp[3];
  double m0ElcRef_curr;
  m0ElcRef_curr = 0.6324555320336759*nphi_ref[2]-0.9486832980505137*nphi_ref[1]+0.7071067811865475*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[0] = -((1.0*T_e*log(fmax(0.6324555320336759*m0Elc[2]-0.9486832980505137*m0Elc[1]+0.7071067811865475*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.6324555320336759*nphi_ref[5]-0.9486832980505137*nphi_ref[4]+0.7071067811865475*nphi_ref[3];
  } else {
    phi_qp[0] = 0.0;
  }
  m0ElcRef_curr = 0.7071067811865475*nphi_ref[0]-0.7905694150420947*nphi_ref[2];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[1] = -((1.0*T_e*log(fmax(0.7071067811865475*m0Elc[0]-0.7905694150420947*m0Elc[2],m0ElcRef_curr)/m0ElcRef_curr))/q_e)-0.7905694150420947*nphi_ref[5]+0.7071067811865475*nphi_ref[3];
  } else {
    phi_qp[1] = 0.0;
  }
  m0ElcRef_curr = 0.6324555320336759*nphi_ref[2]+0.9486832980505137*nphi_ref[1]+0.7071067811865475*nphi_ref[0];
  if ((gkyl_is_dfinite(m0ElcRef_curr)) && (m0ElcRef_curr>0.)) {
    phi_qp[2] = -((1.0*T_e*log(fmax(0.6324555320336759*m0Elc[2]+0.9486832980505137*m0Elc[1]+0.7071067811865475*m0Elc[0],m0ElcRef_curr)/m0ElcRef_curr))/q_e)+0.6324555320336759*nphi_ref[5]+0.9486832980505137*nphi_ref[4]+0.7071067811865475*nphi_ref[3];
  } else {
    phi_qp[2] = 0.0;
  }

  phi[0] = 0.39283710065919303*phi_qp[2]+0.6285393610547091*phi_qp[1]+0.39283710065919303*phi_qp[0]; 
  phi[1] = 0.5270462766947298*phi_qp[2]-0.5270462766947298*phi_qp[0]; 
  phi[2] = 0.35136418446315326*phi_qp[2]-0.7027283689263066*phi_qp[1]+0.35136418446315326*phi_qp[0]; 
}

