#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_2x3v_ser_p1(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms) 
{
  // lhs: Matrix to be inverted to solve Ax = rhs. 
  // rhs: Right-hand-side of linear system. 
  // fpo_moms: Volume correction moments for FPO conservation 
  // boundary_corrections: Boundary correction moments for FPO conservation 
  // moms: m0, m1i, and m2 
 
  // Index into moment array 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[4]; 
  const double *m1y = &moms[8]; 
  const double *m1z = &moms[12]; 
 
  // Index into volume and boundary correction moments 
  const double *vol_corr_ax = &fpo_moms[0]; 
  const double *vol_corr_ay = &fpo_moms[4]; 
  const double *vol_corr_az = &fpo_moms[8]; 
  const double *vol_corr_energy = &fpo_moms[12]; 
 
  const double *bcorr_ax = &boundary_corrections[0]; 
  const double *bcorr_ay = &boundary_corrections[4]; 
  const double *bcorr_az = &boundary_corrections[8]; 
  const double *bcorr_energy = &boundary_corrections[12]; 
  const double *bcorr_ax_D_ij = &boundary_corrections[16]; 
  const double *bcorr_ay_D_ij = &boundary_corrections[20]; 
  const double *bcorr_az_D_ij = &boundary_corrections[24]; 
  const double *bcorr_energy_D_ij = &boundary_corrections[28]; 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 0, 0, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 0, 1, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 0, 2, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 0, 3, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 1, 0, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 1, 1, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 1, 2, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 1, 3, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 2, 0, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 2, 1, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 2, 2, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 2, 3, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 3, 0, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 3, 1, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 3, 2, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 3, 3, 0.5*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 4, 4, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 4, 5, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 4, 6, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 4, 7, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 5, 4, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 5, 5, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 5, 6, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 5, 7, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 6, 4, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 6, 5, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 6, 6, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 6, 7, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 7, 4, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 7, 5, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 7, 6, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 7, 7, 0.5*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 8, 8, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 8, 9, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 8, 10, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 8, 11, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 9, 8, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 9, 9, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 9, 10, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 9, 11, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 10, 8, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 10, 9, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 10, 10, 0.5*m0[0]); 
  gkyl_mat_set(lhs, 10, 11, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 11, 8, 0.5*m0[3]); 
  gkyl_mat_set(lhs, 11, 9, 0.5*m0[2]); 
  gkyl_mat_set(lhs, 11, 10, 0.5*m0[1]); 
  gkyl_mat_set(lhs, 11, 11, 0.5*m0[0]); 
 
  // Block from weak multiply of M1x with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 12, 0, 0.5*m1x[0]); 
  gkyl_mat_set(lhs, 12, 1, 0.5*m1x[1]); 
  gkyl_mat_set(lhs, 12, 2, 0.5*m1x[2]); 
  gkyl_mat_set(lhs, 12, 3, 0.5*m1x[3]); 
  gkyl_mat_set(lhs, 13, 0, 0.5*m1x[1]); 
  gkyl_mat_set(lhs, 13, 1, 0.5*m1x[0]); 
  gkyl_mat_set(lhs, 13, 2, 0.5*m1x[3]); 
  gkyl_mat_set(lhs, 13, 3, 0.5*m1x[2]); 
  gkyl_mat_set(lhs, 14, 0, 0.5*m1x[2]); 
  gkyl_mat_set(lhs, 14, 1, 0.5*m1x[3]); 
  gkyl_mat_set(lhs, 14, 2, 0.5*m1x[0]); 
  gkyl_mat_set(lhs, 14, 3, 0.5*m1x[1]); 
  gkyl_mat_set(lhs, 15, 0, 0.5*m1x[3]); 
  gkyl_mat_set(lhs, 15, 1, 0.5*m1x[2]); 
  gkyl_mat_set(lhs, 15, 2, 0.5*m1x[1]); 
  gkyl_mat_set(lhs, 15, 3, 0.5*m1x[0]); 
 
  // Block from weak multiply of M1y with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 12, 4, 0.5*m1y[0]); 
  gkyl_mat_set(lhs, 12, 5, 0.5*m1y[1]); 
  gkyl_mat_set(lhs, 12, 6, 0.5*m1y[2]); 
  gkyl_mat_set(lhs, 12, 7, 0.5*m1y[3]); 
  gkyl_mat_set(lhs, 13, 4, 0.5*m1y[1]); 
  gkyl_mat_set(lhs, 13, 5, 0.5*m1y[0]); 
  gkyl_mat_set(lhs, 13, 6, 0.5*m1y[3]); 
  gkyl_mat_set(lhs, 13, 7, 0.5*m1y[2]); 
  gkyl_mat_set(lhs, 14, 4, 0.5*m1y[2]); 
  gkyl_mat_set(lhs, 14, 5, 0.5*m1y[3]); 
  gkyl_mat_set(lhs, 14, 6, 0.5*m1y[0]); 
  gkyl_mat_set(lhs, 14, 7, 0.5*m1y[1]); 
  gkyl_mat_set(lhs, 15, 4, 0.5*m1y[3]); 
  gkyl_mat_set(lhs, 15, 5, 0.5*m1y[2]); 
  gkyl_mat_set(lhs, 15, 6, 0.5*m1y[1]); 
  gkyl_mat_set(lhs, 15, 7, 0.5*m1y[0]); 
 
  // Block from weak multiply of M1z with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 12, 8, 0.5*m1z[0]); 
  gkyl_mat_set(lhs, 12, 9, 0.5*m1z[1]); 
  gkyl_mat_set(lhs, 12, 10, 0.5*m1z[2]); 
  gkyl_mat_set(lhs, 12, 11, 0.5*m1z[3]); 
  gkyl_mat_set(lhs, 13, 8, 0.5*m1z[1]); 
  gkyl_mat_set(lhs, 13, 9, 0.5*m1z[0]); 
  gkyl_mat_set(lhs, 13, 10, 0.5*m1z[3]); 
  gkyl_mat_set(lhs, 13, 11, 0.5*m1z[2]); 
  gkyl_mat_set(lhs, 14, 8, 0.5*m1z[2]); 
  gkyl_mat_set(lhs, 14, 9, 0.5*m1z[3]); 
  gkyl_mat_set(lhs, 14, 10, 0.5*m1z[0]); 
  gkyl_mat_set(lhs, 14, 11, 0.5*m1z[1]); 
  gkyl_mat_set(lhs, 15, 8, 0.5*m1z[3]); 
  gkyl_mat_set(lhs, 15, 9, 0.5*m1z[2]); 
  gkyl_mat_set(lhs, 15, 10, 0.5*m1z[1]); 
  gkyl_mat_set(lhs, 15, 11, 0.5*m1z[0]); 
 
  // Block from weak multiply of -bcorr_ax with diff_coeff_corr. 
  gkyl_mat_set(lhs, 0, 12, -(0.5*bcorr_ax[0])); 
  gkyl_mat_set(lhs, 0, 13, -(0.5*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 0, 14, -(0.5*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 0, 15, -(0.5*bcorr_ax[3])); 
  gkyl_mat_set(lhs, 1, 12, -(0.5*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 1, 13, -(0.5*bcorr_ax[0])); 
  gkyl_mat_set(lhs, 1, 14, -(0.5*bcorr_ax[3])); 
  gkyl_mat_set(lhs, 1, 15, -(0.5*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 2, 12, -(0.5*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 2, 13, -(0.5*bcorr_ax[3])); 
  gkyl_mat_set(lhs, 2, 14, -(0.5*bcorr_ax[0])); 
  gkyl_mat_set(lhs, 2, 15, -(0.5*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 3, 12, -(0.5*bcorr_ax[3])); 
  gkyl_mat_set(lhs, 3, 13, -(0.5*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 3, 14, -(0.5*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 3, 15, -(0.5*bcorr_ax[0])); 
 
  // Block from weak multiply of -bcorr_ay with diff_coeff_corr. 
  gkyl_mat_set(lhs, 4, 12, -(0.5*bcorr_ay[0])); 
  gkyl_mat_set(lhs, 4, 13, -(0.5*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 4, 14, -(0.5*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 4, 15, -(0.5*bcorr_ay[3])); 
  gkyl_mat_set(lhs, 5, 12, -(0.5*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 5, 13, -(0.5*bcorr_ay[0])); 
  gkyl_mat_set(lhs, 5, 14, -(0.5*bcorr_ay[3])); 
  gkyl_mat_set(lhs, 5, 15, -(0.5*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 6, 12, -(0.5*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 6, 13, -(0.5*bcorr_ay[3])); 
  gkyl_mat_set(lhs, 6, 14, -(0.5*bcorr_ay[0])); 
  gkyl_mat_set(lhs, 6, 15, -(0.5*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 7, 12, -(0.5*bcorr_ay[3])); 
  gkyl_mat_set(lhs, 7, 13, -(0.5*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 7, 14, -(0.5*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 7, 15, -(0.5*bcorr_ay[0])); 
 
  // Block from weak multiply of -bcorr_az with diff_coeff_corr. 
  gkyl_mat_set(lhs, 8, 12, -(0.5*bcorr_az[0])); 
  gkyl_mat_set(lhs, 8, 13, -(0.5*bcorr_az[1])); 
  gkyl_mat_set(lhs, 8, 14, -(0.5*bcorr_az[2])); 
  gkyl_mat_set(lhs, 8, 15, -(0.5*bcorr_az[3])); 
  gkyl_mat_set(lhs, 9, 12, -(0.5*bcorr_az[1])); 
  gkyl_mat_set(lhs, 9, 13, -(0.5*bcorr_az[0])); 
  gkyl_mat_set(lhs, 9, 14, -(0.5*bcorr_az[3])); 
  gkyl_mat_set(lhs, 9, 15, -(0.5*bcorr_az[2])); 
  gkyl_mat_set(lhs, 10, 12, -(0.5*bcorr_az[2])); 
  gkyl_mat_set(lhs, 10, 13, -(0.5*bcorr_az[3])); 
  gkyl_mat_set(lhs, 10, 14, -(0.5*bcorr_az[0])); 
  gkyl_mat_set(lhs, 10, 15, -(0.5*bcorr_az[1])); 
  gkyl_mat_set(lhs, 11, 12, -(0.5*bcorr_az[3])); 
  gkyl_mat_set(lhs, 11, 13, -(0.5*bcorr_az[2])); 
  gkyl_mat_set(lhs, 11, 14, -(0.5*bcorr_az[1])); 
  gkyl_mat_set(lhs, 11, 15, -(0.5*bcorr_az[0])); 
 
  // Block from weak multiply of (vdim*M0-bcorr_energy) with diff_coeff_corr. 
  gkyl_mat_set(lhs, 12, 12, 1.5*m0[0]-0.5*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 12, 13, 1.5*m0[1]-0.5*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 12, 14, 1.5*m0[2]-0.5*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 12, 15, 1.5*m0[3]-0.5*bcorr_energy[3]); 
  gkyl_mat_set(lhs, 13, 12, 1.5*m0[1]-0.5*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 13, 13, 1.5*m0[0]-0.5*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 13, 14, 1.5*m0[3]-0.5*bcorr_energy[3]); 
  gkyl_mat_set(lhs, 13, 15, 1.5*m0[2]-0.5*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 14, 12, 1.5*m0[2]-0.5*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 14, 13, 1.5*m0[3]-0.5*bcorr_energy[3]); 
  gkyl_mat_set(lhs, 14, 14, 1.5*m0[0]-0.5*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 14, 15, 1.5*m0[1]-0.5*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 15, 12, 1.5*m0[3]-0.5*bcorr_energy[3]); 
  gkyl_mat_set(lhs, 15, 13, 1.5*m0[2]-0.5*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 15, 14, 1.5*m0[1]-0.5*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 15, 15, 1.5*m0[0]-0.5*bcorr_energy[0]); 
 
  // Set rhs vector. 
  gkyl_mat_set(rhs, 0, 0, bcorr_ax_D_ij[0]-vol_corr_ax[0]); 
  gkyl_mat_set(rhs, 1, 0, bcorr_ax_D_ij[1]-vol_corr_ax[1]); 
  gkyl_mat_set(rhs, 2, 0, bcorr_ax_D_ij[2]-vol_corr_ax[2]); 
  gkyl_mat_set(rhs, 3, 0, bcorr_ax_D_ij[3]-vol_corr_ax[3]); 
  gkyl_mat_set(rhs, 4, 0, bcorr_ay_D_ij[0]-vol_corr_ay[0]); 
  gkyl_mat_set(rhs, 5, 0, bcorr_ay_D_ij[1]-vol_corr_ay[1]); 
  gkyl_mat_set(rhs, 6, 0, bcorr_ay_D_ij[2]-vol_corr_ay[2]); 
  gkyl_mat_set(rhs, 7, 0, bcorr_ay_D_ij[3]-vol_corr_ay[3]); 
  gkyl_mat_set(rhs, 8, 0, bcorr_az_D_ij[0]-vol_corr_az[0]); 
  gkyl_mat_set(rhs, 9, 0, bcorr_az_D_ij[1]-vol_corr_az[1]); 
  gkyl_mat_set(rhs, 10, 0, bcorr_az_D_ij[2]-vol_corr_az[2]); 
  gkyl_mat_set(rhs, 11, 0, bcorr_az_D_ij[3]-vol_corr_az[3]); 
  gkyl_mat_set(rhs, 12, 0, bcorr_energy_D_ij[0]-vol_corr_energy[0]); 
  gkyl_mat_set(rhs, 13, 0, bcorr_energy_D_ij[1]-vol_corr_energy[1]); 
  gkyl_mat_set(rhs, 14, 0, bcorr_energy_D_ij[2]-vol_corr_energy[2]); 
  gkyl_mat_set(rhs, 15, 0, bcorr_energy_D_ij[3]-vol_corr_energy[3]); 
}
