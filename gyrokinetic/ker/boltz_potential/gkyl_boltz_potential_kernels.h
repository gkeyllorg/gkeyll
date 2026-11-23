#pragma once
#include <gkyl_util.h>
#include <math.h>
EXTERN_C_BEG

GKYL_CU_DH void boltz_potential_boundary_lower_1x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_1x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_1x_ser_p1(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_1x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_1x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_potential_boundary_lower_1x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_1x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_1x_ser_p2(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_1x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_1x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_potential_boundary_lower_2x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_2x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_2x_ser_p1(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_2x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_2x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_potential_boundary_lower_2x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_2x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_2x_ser_p2(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_2x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_2x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_potential_boundary_lower_3x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_3x_ser_p1(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_3x_ser_p1(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_3x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_3x_ser_p1(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_potential_boundary_lower_3x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_boundary_upper_3x_ser_p2(const double *flux, const double *area, const double *m0_elc, const double *wall, double *out);
GKYL_CU_DH void boltz_potential_phi_3x_ser_p2(double q_e, double T_e, const double *m0Elc, const double *nphi_ref, double *phi);
GKYL_CU_DH void boltz_pot_ambi_reference_3x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);
GKYL_CU_DH void boltz_pot_0curr_reference_3x_ser_p2(double q_e, double T_e, double m_e, const double *lower, const double *upper, double *out);

EXTERN_C_END
