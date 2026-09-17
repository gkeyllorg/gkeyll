#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p1(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms) 
{
  // lhs: Matrix to be inverted to solve Ax = rhs. 
  // rhs: Right-hand-side of linear system. 
  // fpo_moms: Volume correction moments for FPO conservation 
  // boundary_corrections: Boundary correction moments for FPO conservation 
  // moms: m0, m1i, and m2 
 
  // Index into moment array 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
  const double *m1z = &moms[6]; 
 
  // Index into volume and boundary correction moments 
  const double *vol_corr_ax = &fpo_moms[0]; 
  const double *vol_corr_ay = &fpo_moms[2]; 
  const double *vol_corr_az = &fpo_moms[4]; 
  const double *vol_corr_energy = &fpo_moms[6]; 
 
  const double *bcorr_ax = &boundary_corrections[0]; 
  const double *bcorr_ay = &boundary_corrections[2]; 
  const double *bcorr_az = &boundary_corrections[4]; 
  const double *bcorr_energy = &boundary_corrections[6]; 
  const double *bcorr_ax_D_ij = &boundary_corrections[8]; 
  const double *bcorr_ay_D_ij = &boundary_corrections[10]; 
  const double *bcorr_az_D_ij = &boundary_corrections[12]; 
  const double *bcorr_energy_D_ij = &boundary_corrections[14]; 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 0, 0, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 0, 1, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 1, 0, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 1, 1, 0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 2, 2, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 2, 3, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 3, 2, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 3, 3, 0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 4, 4, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 4, 5, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 5, 4, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 5, 5, 0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M1x with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 6, 0, 0.7071067811865475*m1x[0]); 
  gkyl_mat_set(lhs, 6, 1, 0.7071067811865475*m1x[1]); 
  gkyl_mat_set(lhs, 7, 0, 0.7071067811865475*m1x[1]); 
  gkyl_mat_set(lhs, 7, 1, 0.7071067811865475*m1x[0]); 
 
  // Block from weak multiply of M1y with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 6, 2, 0.7071067811865475*m1y[0]); 
  gkyl_mat_set(lhs, 6, 3, 0.7071067811865475*m1y[1]); 
  gkyl_mat_set(lhs, 7, 2, 0.7071067811865475*m1y[1]); 
  gkyl_mat_set(lhs, 7, 3, 0.7071067811865475*m1y[0]); 
 
  // Block from weak multiply of M1z with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 6, 4, 0.7071067811865475*m1z[0]); 
  gkyl_mat_set(lhs, 6, 5, 0.7071067811865475*m1z[1]); 
  gkyl_mat_set(lhs, 7, 4, 0.7071067811865475*m1z[1]); 
  gkyl_mat_set(lhs, 7, 5, 0.7071067811865475*m1z[0]); 
 
  // Block from weak multiply of -bcorr_ax with diff_coeff_corr. 
  gkyl_mat_set(lhs, 0, 6, -(0.7071067811865475*bcorr_ax[0])); 
  gkyl_mat_set(lhs, 0, 7, -(0.7071067811865475*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 1, 6, -(0.7071067811865475*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 1, 7, -(0.7071067811865475*bcorr_ax[0])); 
 
  // Block from weak multiply of -bcorr_ay with diff_coeff_corr. 
  gkyl_mat_set(lhs, 2, 6, -(0.7071067811865475*bcorr_ay[0])); 
  gkyl_mat_set(lhs, 2, 7, -(0.7071067811865475*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 3, 6, -(0.7071067811865475*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 3, 7, -(0.7071067811865475*bcorr_ay[0])); 
 
  // Block from weak multiply of -bcorr_az with diff_coeff_corr. 
  gkyl_mat_set(lhs, 4, 6, -(0.7071067811865475*bcorr_az[0])); 
  gkyl_mat_set(lhs, 4, 7, -(0.7071067811865475*bcorr_az[1])); 
  gkyl_mat_set(lhs, 5, 6, -(0.7071067811865475*bcorr_az[1])); 
  gkyl_mat_set(lhs, 5, 7, -(0.7071067811865475*bcorr_az[0])); 
 
  // Block from weak multiply of (vdim*M0-bcorr_energy) with diff_coeff_corr. 
  gkyl_mat_set(lhs, 6, 6, 2.1213203435596424*m0[0]-0.7071067811865475*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 6, 7, 2.1213203435596424*m0[1]-0.7071067811865475*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 7, 6, 2.1213203435596424*m0[1]-0.7071067811865475*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 7, 7, 2.1213203435596424*m0[0]-0.7071067811865475*bcorr_energy[0]); 
 
  // Set rhs vector. 
  gkyl_mat_set(rhs, 0, 0, bcorr_ax_D_ij[0]-vol_corr_ax[0]); 
  gkyl_mat_set(rhs, 1, 0, bcorr_ax_D_ij[1]-vol_corr_ax[1]); 
  gkyl_mat_set(rhs, 2, 0, bcorr_ay_D_ij[0]-vol_corr_ay[0]); 
  gkyl_mat_set(rhs, 3, 0, bcorr_ay_D_ij[1]-vol_corr_ay[1]); 
  gkyl_mat_set(rhs, 4, 0, bcorr_az_D_ij[0]-vol_corr_az[0]); 
  gkyl_mat_set(rhs, 5, 0, bcorr_az_D_ij[1]-vol_corr_az[1]); 
  gkyl_mat_set(rhs, 6, 0, bcorr_energy_D_ij[0]-vol_corr_energy[0]); 
  gkyl_mat_set(rhs, 7, 0, bcorr_energy_D_ij[1]-vol_corr_energy[1]); 
}
