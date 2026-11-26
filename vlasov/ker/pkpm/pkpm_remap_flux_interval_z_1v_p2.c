#include <gkyl_pkpm_remap_flux_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_remap_flux_interval_z_1v_p2(const double *fc, double w_lo, double w_up,
  const double *side_l, const double *side_r, const double *cell_l, const double *cell_r,
  const double *Y_l, const double *Y_r, double* GKYL_RESTRICT P_l, double* GKYL_RESTRICT P_r,
  double* GKYL_RESTRICT mom, double* GKYL_RESTRICT share) 
{ 
  // fc:          Face coefficients of the flux at the transverse node (PKPM_FC_*, stride 1). 
  // w_lo/w_up:   Ends of the interval in the face-frame label. 
  // side_l/r:    Map of the left/right side to the face frame (PKPM_RF_SIDE_*). 
  // cell_l/r:    Velocity cell of the left/right side on the interval (PKPM_RF_CELL_*). 
  // Y_l/Y_r:     Blocks of the left/right side in the face frame, Legendre coefficients in the 
  //              cell coordinate of the side: Y[m*8 + c] of [1, eta, 3 eta^2 - 1]. 
  // P_l/P_r:     Accumulated projection of the flux onto the Legendre modes of the cell of the 
  //              left/right side, face frame: P[m*8 + c]. 
  // mom:         Accumulated face-frame moments of the flux (PKPM_RF_NUM_MOMS). 
  // share:       Accumulated shares of the flux for the blend (PKPM_RF_SHARE_*). 
  // 3-point Gauss-Legendre rule on the interval. 
  const double ord[3] = {-0.7745966692414834, 0.0, 0.7745966692414834}; 
  const double wgt[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556}; 
  const double mid = 0.5*(w_up+w_lo), half = 0.5*(w_up-w_lo); 
  double X_l[8], X_r[8], phi_l[8], phi_r[8], flux[8], sb[2]; 
  for (int g = 0; g < 3; ++g) { 
    const double wb = mid+ord[g]*half, om = wgt[g]*half; 
    const double eta_l = cell_l[0]+cell_l[1]*wb, leg2_l = 3.0*eta_l*eta_l-1.0; 
    const double eta_r = cell_r[0]+cell_r[1]*wb, leg2_r = 3.0*eta_r*eta_r-1.0; 
    for (int c = 0; c < 8; ++c) { 
      X_l[c] = Y_l[c]+Y_l[8+c]*eta_l+Y_l[16+c]*leg2_l; 
      X_r[c] = Y_r[c]+Y_r[8+c]*eta_r+Y_r[16+c]*leg2_r; 
    } 
    pkpm_harm_conf_flux(2, X_l, X_r, wb, fc, 1, phi_l, phi_r, sb); 
    for (int c = 0; c < 8; ++c) { 
      flux[c] = phi_l[c]+phi_r[c]; 
    } 
    mom[0] += om*flux[0]; 
    mom[1] += om*wb*flux[0]; 
    mom[2] += om*wb*wb*flux[0]; 
    mom[3] += om*flux[1]; 
    mom[4] += om*flux[2]; 
    mom[5] += om*flux[3]; 
    mom[6] += om*flux[4]; 
    const double mass = om*fabs(flux[0]); 
    const double w_l = side_l[0]*wb+side_l[1]; 
    const double cm_l = side_l[5]*flux[2]+side_l[6]*flux[3]+side_l[7]*flux[4]; 
    const double yf_l = side_l[2]*flux[0]; 
    const double yg_l = side_l[3]*(flux[1]-cm_l+0.5*side_l[4]*flux[0]); 
    const double mt_l = mass*fabs(w_l); 
    const double en_l = om*(fabs(yf_l)*w_l*w_l+2.0*fabs(yg_l)); 
    share[0] += cell_l[2]*mt_l; 
    share[2] += mt_l; 
    share[4] += cell_l[2]*en_l; 
    share[6] += en_l; 
    const double w_r = side_r[0]*wb+side_r[1]; 
    const double cm_r = side_r[5]*flux[2]+side_r[6]*flux[3]+side_r[7]*flux[4]; 
    const double yf_r = side_r[2]*flux[0]; 
    const double yg_r = side_r[3]*(flux[1]-cm_r+0.5*side_r[4]*flux[0]); 
    const double mt_r = mass*fabs(w_r); 
    const double en_r = om*(fabs(yf_r)*w_r*w_r+2.0*fabs(yg_r)); 
    share[1] += cell_r[2]*mt_r; 
    share[3] += mt_r; 
    share[5] += cell_r[2]*en_r; 
    share[7] += en_r; 
    const double wt_l = om*cell_l[1]; 
    const double t0_l = 0.7071067811865475*wt_l, t1_l = 1.224744871391589*wt_l*eta_l, t2_l = 0.7905694150420947*wt_l*leg2_l; 
    const double wt_r = om*cell_r[1]; 
    const double t0_r = 0.7071067811865475*wt_r, t1_r = 1.224744871391589*wt_r*eta_r, t2_r = 0.7905694150420947*wt_r*leg2_r; 
    for (int c = 0; c < 8; ++c) { 
      P_l[c] += t0_l*flux[c]; 
      P_l[8+c] += t1_l*flux[c]; 
      P_l[16+c] += t2_l*flux[c]; 
      P_r[c] += t0_r*flux[c]; 
      P_r[8+c] += t1_r*flux[c]; 
      P_r[16+c] += t2_r*flux[c]; 
    } 
  } 
} 
