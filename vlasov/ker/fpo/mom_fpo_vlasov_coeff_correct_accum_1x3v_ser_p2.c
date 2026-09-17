#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_1x3v_ser_p2(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf) 
{
  // drag_diff_coeff_corrs: Corrections to be added to coeffs, function of config space only. 
  // drag_coeff: FPO drag coefficient. 
  // diff_coeff: FPO diffusion coefficient. 
 
  // Index into drag and diffusion coefficients. 
  double *ax = &drag_coeff[0]; 
  double *ay = &drag_coeff[48]; 
  double *az = &drag_coeff[96]; 
  double *Dxx = &diff_coeff[0]; 
  double *Dyy = &diff_coeff[192]; 
  double *Dzz = &diff_coeff[384]; 
 
  // Index into surface expansions. 
  double *ax_surf = &drag_coeff_surf[0]; 
  double *ay_surf = &drag_coeff_surf[20]; 
  double *az_surf = &drag_coeff_surf[40]; 
  double *Dxx_surf = &diff_coeff_surf[0]; 
  double *Dyy_surf = &diff_coeff_surf[160]; 
  double *Dzz_surf = &diff_coeff_surf[320]; 
 
  // Index into correction array. 
  const double* ax_corr = &drag_diff_coeff_corrs[0]; 
  const double* ay_corr = &drag_diff_coeff_corrs[3]; 
  const double* az_corr = &drag_diff_coeff_corrs[6]; 
  const double* D_corr = &drag_diff_coeff_corrs[9]; 
 
  ax[0] += 2.8284271247461907*ax_corr[0]; 
  ax[1] += 2.8284271247461907*ax_corr[1]; 
  ax[11] += 2.8284271247461907*ax_corr[2]; 
 
  ay[0] += 2.8284271247461907*ay_corr[0]; 
  ay[1] += 2.8284271247461907*ay_corr[1]; 
  ay[11] += 2.8284271247461907*ay_corr[2]; 
 
  az[0] += 2.8284271247461907*az_corr[0]; 
  az[1] += 2.8284271247461907*az_corr[1]; 
  az[11] += 2.8284271247461907*az_corr[2]; 
 
  Dxx[0] += 2.8284271247461907*D_corr[0]; 
  Dxx[1] += 2.8284271247461907*D_corr[1]; 
  Dxx[11] += 2.8284271247461907*D_corr[2]; 
 
  Dyy[0] += 2.8284271247461907*D_corr[0]; 
  Dyy[1] += 2.8284271247461907*D_corr[1]; 
  Dyy[11] += 2.8284271247461907*D_corr[2]; 
 
  Dzz[0] += 2.8284271247461907*D_corr[0]; 
  Dzz[1] += 2.8284271247461907*D_corr[1]; 
  Dzz[11] += 2.8284271247461907*D_corr[2]; 
 
  ax_surf[0] += 2.0*ax_corr[0]; 
  ax_surf[1] += 2.0*ax_corr[1]; 
  ax_surf[7] += 2.0*ax_corr[2]; 
 
  ay_surf[0] += 2.0*ay_corr[0]; 
  ay_surf[1] += 2.0*ay_corr[1]; 
  ay_surf[7] += 2.0*ay_corr[2]; 
 
  az_surf[0] += 2.0*az_corr[0]; 
  az_surf[1] += 2.0*az_corr[1]; 
  az_surf[7] += 2.0*az_corr[2]; 
 
  Dxx_surf[0] += 2.0*D_corr[0]; 
  Dxx_surf[1] += 2.0*D_corr[1]; 
  Dxx_surf[7] += 2.0*D_corr[2]; 
 
  Dyy_surf[0] += 2.0*D_corr[0]; 
  Dyy_surf[1] += 2.0*D_corr[1]; 
  Dyy_surf[7] += 2.0*D_corr[2]; 
 
  Dzz_surf[0] += 2.0*D_corr[0]; 
  Dzz_surf[1] += 2.0*D_corr[1]; 
  Dzz_surf[7] += 2.0*D_corr[2]; 
} 

