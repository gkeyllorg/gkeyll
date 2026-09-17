#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_mat_1x3v_ser_p2(struct gkyl_mat *lhs, struct gkyl_mat *rhs, const double *fpo_moms, const double *boundary_corrections, const double *moms) 
{
  // lhs: Matrix to be inverted to solve Ax = rhs. 
  // rhs: Right-hand-side of linear system. 
  // fpo_moms: Volume correction moments for FPO conservation 
  // boundary_corrections: Boundary correction moments for FPO conservation 
  // moms: m0, m1i, and m2 
 
  // Index into moment array 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[3]; 
  const double *m1y = &moms[6]; 
  const double *m1z = &moms[9]; 
 
  // Index into volume and boundary correction moments 
  const double *vol_corr_ax = &fpo_moms[0]; 
  const double *vol_corr_ay = &fpo_moms[3]; 
  const double *vol_corr_az = &fpo_moms[6]; 
  const double *vol_corr_energy = &fpo_moms[9]; 
 
  const double *bcorr_ax = &boundary_corrections[0]; 
  const double *bcorr_ay = &boundary_corrections[3]; 
  const double *bcorr_az = &boundary_corrections[6]; 
  const double *bcorr_energy = &boundary_corrections[9]; 
  const double *bcorr_ax_D_ij = &boundary_corrections[12]; 
  const double *bcorr_ay_D_ij = &boundary_corrections[15]; 
  const double *bcorr_az_D_ij = &boundary_corrections[18]; 
  const double *bcorr_energy_D_ij = &boundary_corrections[21]; 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 0, 0, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 0, 1, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 0, 2, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 1, 0, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 1, 1, 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 1, 2, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 2, 0, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 2, 1, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 2, 2, 0.45175395145262565*m0[2]+0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 3, 3, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 3, 4, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 3, 5, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 4, 3, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 4, 4, 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 4, 5, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 5, 3, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 5, 4, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 5, 5, 0.45175395145262565*m0[2]+0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M0 with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 6, 6, 0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 6, 7, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 6, 8, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 7, 6, 0.7071067811865475*m0[1]); 
  gkyl_mat_set(lhs, 7, 7, 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]); 
  gkyl_mat_set(lhs, 7, 8, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 8, 6, 0.7071067811865475*m0[2]); 
  gkyl_mat_set(lhs, 8, 7, 0.6324555320336759*m0[1]); 
  gkyl_mat_set(lhs, 8, 8, 0.45175395145262565*m0[2]+0.7071067811865475*m0[0]); 
 
  // Block from weak multiply of M1x with drag_coeff_corr_vx. 
  gkyl_mat_set(lhs, 9, 0, 0.7071067811865475*m1x[0]); 
  gkyl_mat_set(lhs, 9, 1, 0.7071067811865475*m1x[1]); 
  gkyl_mat_set(lhs, 9, 2, 0.7071067811865475*m1x[2]); 
  gkyl_mat_set(lhs, 10, 0, 0.7071067811865475*m1x[1]); 
  gkyl_mat_set(lhs, 10, 1, 0.6324555320336759*m1x[2]+0.7071067811865475*m1x[0]); 
  gkyl_mat_set(lhs, 10, 2, 0.6324555320336759*m1x[1]); 
  gkyl_mat_set(lhs, 11, 0, 0.7071067811865475*m1x[2]); 
  gkyl_mat_set(lhs, 11, 1, 0.6324555320336759*m1x[1]); 
  gkyl_mat_set(lhs, 11, 2, 0.45175395145262565*m1x[2]+0.7071067811865475*m1x[0]); 
 
  // Block from weak multiply of M1y with drag_coeff_corr_vy. 
  gkyl_mat_set(lhs, 9, 3, 0.7071067811865475*m1y[0]); 
  gkyl_mat_set(lhs, 9, 4, 0.7071067811865475*m1y[1]); 
  gkyl_mat_set(lhs, 9, 5, 0.7071067811865475*m1y[2]); 
  gkyl_mat_set(lhs, 10, 3, 0.7071067811865475*m1y[1]); 
  gkyl_mat_set(lhs, 10, 4, 0.6324555320336759*m1y[2]+0.7071067811865475*m1y[0]); 
  gkyl_mat_set(lhs, 10, 5, 0.6324555320336759*m1y[1]); 
  gkyl_mat_set(lhs, 11, 3, 0.7071067811865475*m1y[2]); 
  gkyl_mat_set(lhs, 11, 4, 0.6324555320336759*m1y[1]); 
  gkyl_mat_set(lhs, 11, 5, 0.45175395145262565*m1y[2]+0.7071067811865475*m1y[0]); 
 
  // Block from weak multiply of M1z with drag_coeff_corr_vz. 
  gkyl_mat_set(lhs, 9, 6, 0.7071067811865475*m1z[0]); 
  gkyl_mat_set(lhs, 9, 7, 0.7071067811865475*m1z[1]); 
  gkyl_mat_set(lhs, 9, 8, 0.7071067811865475*m1z[2]); 
  gkyl_mat_set(lhs, 10, 6, 0.7071067811865475*m1z[1]); 
  gkyl_mat_set(lhs, 10, 7, 0.6324555320336759*m1z[2]+0.7071067811865475*m1z[0]); 
  gkyl_mat_set(lhs, 10, 8, 0.6324555320336759*m1z[1]); 
  gkyl_mat_set(lhs, 11, 6, 0.7071067811865475*m1z[2]); 
  gkyl_mat_set(lhs, 11, 7, 0.6324555320336759*m1z[1]); 
  gkyl_mat_set(lhs, 11, 8, 0.45175395145262565*m1z[2]+0.7071067811865475*m1z[0]); 
 
  // Block from weak multiply of -bcorr_ax with diff_coeff_corr. 
  gkyl_mat_set(lhs, 0, 9, -(0.7071067811865475*bcorr_ax[0])); 
  gkyl_mat_set(lhs, 0, 10, -(0.7071067811865475*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 0, 11, -(0.7071067811865475*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 1, 9, -(0.7071067811865475*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 1, 10, -(0.6324555320336759*bcorr_ax[2])-0.7071067811865475*bcorr_ax[0]); 
  gkyl_mat_set(lhs, 1, 11, -(0.6324555320336759*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 2, 9, -(0.7071067811865475*bcorr_ax[2])); 
  gkyl_mat_set(lhs, 2, 10, -(0.6324555320336759*bcorr_ax[1])); 
  gkyl_mat_set(lhs, 2, 11, -(0.45175395145262565*bcorr_ax[2])-0.7071067811865475*bcorr_ax[0]); 
 
  // Block from weak multiply of -bcorr_ay with diff_coeff_corr. 
  gkyl_mat_set(lhs, 3, 9, -(0.7071067811865475*bcorr_ay[0])); 
  gkyl_mat_set(lhs, 3, 10, -(0.7071067811865475*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 3, 11, -(0.7071067811865475*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 4, 9, -(0.7071067811865475*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 4, 10, -(0.6324555320336759*bcorr_ay[2])-0.7071067811865475*bcorr_ay[0]); 
  gkyl_mat_set(lhs, 4, 11, -(0.6324555320336759*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 5, 9, -(0.7071067811865475*bcorr_ay[2])); 
  gkyl_mat_set(lhs, 5, 10, -(0.6324555320336759*bcorr_ay[1])); 
  gkyl_mat_set(lhs, 5, 11, -(0.45175395145262565*bcorr_ay[2])-0.7071067811865475*bcorr_ay[0]); 
 
  // Block from weak multiply of -bcorr_az with diff_coeff_corr. 
  gkyl_mat_set(lhs, 6, 9, -(0.7071067811865475*bcorr_az[0])); 
  gkyl_mat_set(lhs, 6, 10, -(0.7071067811865475*bcorr_az[1])); 
  gkyl_mat_set(lhs, 6, 11, -(0.7071067811865475*bcorr_az[2])); 
  gkyl_mat_set(lhs, 7, 9, -(0.7071067811865475*bcorr_az[1])); 
  gkyl_mat_set(lhs, 7, 10, -(0.6324555320336759*bcorr_az[2])-0.7071067811865475*bcorr_az[0]); 
  gkyl_mat_set(lhs, 7, 11, -(0.6324555320336759*bcorr_az[1])); 
  gkyl_mat_set(lhs, 8, 9, -(0.7071067811865475*bcorr_az[2])); 
  gkyl_mat_set(lhs, 8, 10, -(0.6324555320336759*bcorr_az[1])); 
  gkyl_mat_set(lhs, 8, 11, -(0.45175395145262565*bcorr_az[2])-0.7071067811865475*bcorr_az[0]); 
 
  // Block from weak multiply of (vdim*M0-bcorr_energy) with diff_coeff_corr. 
  gkyl_mat_set(lhs, 9, 9, 2.1213203435596424*m0[0]-0.7071067811865475*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 9, 10, 2.1213203435596424*m0[1]-0.7071067811865475*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 9, 11, 2.1213203435596424*m0[2]-0.7071067811865475*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 10, 9, 2.1213203435596424*m0[1]-0.7071067811865475*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 10, 10, 1.8973665961010278*m0[2]-0.6324555320336759*bcorr_energy[2]+2.1213203435596424*m0[0]-0.7071067811865475*bcorr_energy[0]); 
  gkyl_mat_set(lhs, 10, 11, 1.8973665961010278*m0[1]-0.6324555320336759*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 11, 9, 2.1213203435596424*m0[2]-0.7071067811865475*bcorr_energy[2]); 
  gkyl_mat_set(lhs, 11, 10, 1.8973665961010278*m0[1]-0.6324555320336759*bcorr_energy[1]); 
  gkyl_mat_set(lhs, 11, 11, 1.355261854357877*m0[2]-0.45175395145262565*bcorr_energy[2]+2.1213203435596424*m0[0]-0.7071067811865475*bcorr_energy[0]); 
 
  // Set rhs vector. 
  gkyl_mat_set(rhs, 0, 0, bcorr_ax_D_ij[0]-vol_corr_ax[0]); 
  gkyl_mat_set(rhs, 1, 0, bcorr_ax_D_ij[1]-vol_corr_ax[1]); 
  gkyl_mat_set(rhs, 2, 0, bcorr_ax_D_ij[2]-vol_corr_ax[2]); 
  gkyl_mat_set(rhs, 3, 0, bcorr_ay_D_ij[0]-vol_corr_ay[0]); 
  gkyl_mat_set(rhs, 4, 0, bcorr_ay_D_ij[1]-vol_corr_ay[1]); 
  gkyl_mat_set(rhs, 5, 0, bcorr_ay_D_ij[2]-vol_corr_ay[2]); 
  gkyl_mat_set(rhs, 6, 0, bcorr_az_D_ij[0]-vol_corr_az[0]); 
  gkyl_mat_set(rhs, 7, 0, bcorr_az_D_ij[1]-vol_corr_az[1]); 
  gkyl_mat_set(rhs, 8, 0, bcorr_az_D_ij[2]-vol_corr_az[2]); 
  gkyl_mat_set(rhs, 9, 0, bcorr_energy_D_ij[0]-vol_corr_energy[0]); 
  gkyl_mat_set(rhs, 10, 0, bcorr_energy_D_ij[1]-vol_corr_energy[1]); 
  gkyl_mat_set(rhs, 11, 0, bcorr_energy_D_ij[2]-vol_corr_energy[2]); 
}
