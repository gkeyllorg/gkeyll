#include <gkyl_pkpm_remap_flux_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_remap_flux_label_x_1v_p2(const double *fc, double w_center, double dw,
  const double *X_l, const double *X_r, double* GKYL_RESTRICT proj, double* GKYL_RESTRICT mom_l,
  double* GKYL_RESTRICT mom_r, double* GKYL_RESTRICT speed) 
{ 
  // fc:          Face coefficients of the flux at the transverse node (PKPM_FC_*, stride 1). 
  // w_center/dw: Center and width of the velocity cell. 
  // X_l/X_r:     Blocks of the left/right side at the node, each in the frame of its side, 
  //              Legendre coefficients in w: X[m*8 + c] of [1, eta, 3 eta^2 - 1]. 
  // proj:        Output projection of the flux onto the Legendre modes of the cell: proj[m*8 + c]. 
  // mom_l/r:     Accumulated moments in w of the contribution of the left/right side 
  //              (PKPM_RF_NUM_MOMS). 
  // speed:       Largest speeds of the flux so far, raised to those of this cell at the node: 
  //              [fastest signal in either direction, largest spread of the two bounds]. 
  // 3-point Gauss-Legendre rule on the cell: the ordinates eta, 3 eta^2 - 1 there, and the 
  // weights times the Legendre modes of the cell. 
  const double eta[3] = {-0.7745966692414834, 0.0, 0.7745966692414834}; 
  const double leg2[3] = {0.8, -1.0, 0.8}; 
  const double wgt[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556}; 
  const double t0[3] = {0.39283710065919303, 0.6285393610547089, 0.39283710065919303}; 
  const double t1[3] = {-0.5270462766947299, 0.0, 0.5270462766947299}; 
  const double t2[3] = {0.35136418446315315, -0.7027283689263063, 0.35136418446315315}; 
  double x_l[8], x_r[8], phi_l[8], phi_r[8], sb[2]; 
  for (int k = 0; k < 24; ++k) { 
    proj[k] = 0.0; 
  } 
  for (int g = 0; g < 3; ++g) { 
    const double wq = w_center+0.5*eta[g]*dw, om = 0.5*wgt[g]*dw; 
    for (int c = 0; c < 8; ++c) { 
      x_l[c] = X_l[c]+eta[g]*X_l[8+c]+leg2[g]*X_l[16+c]; 
      x_r[c] = X_r[c]+eta[g]*X_r[8+c]+leg2[g]*X_r[16+c]; 
    } 
    pkpm_harm_conf_flux(0, x_l, x_r, wq, fc, 1, phi_l, phi_r, sb); 
    speed[0] = fmax(speed[0], fmax(-sb[0], sb[1])); 
    speed[1] = fmax(speed[1], sb[1]-sb[0]); 
    for (int c = 0; c < 8; ++c) { 
      const double flux = phi_l[c]+phi_r[c]; 
      proj[c] += t0[g]*flux; 
      proj[8+c] += t1[g]*flux; 
      proj[16+c] += t2[g]*flux; 
    } 
    mom_l[0] += om*phi_l[0]; 
    mom_l[1] += om*wq*phi_l[0]; 
    mom_l[2] += om*wq*wq*phi_l[0]; 
    mom_l[3] += om*phi_l[1]; 
    mom_l[4] += om*phi_l[2]; 
    mom_l[5] += om*phi_l[3]; 
    mom_l[6] += om*phi_l[4]; 
    mom_r[0] += om*phi_r[0]; 
    mom_r[1] += om*wq*phi_r[0]; 
    mom_r[2] += om*wq*wq*phi_r[0]; 
    mom_r[3] += om*phi_r[1]; 
    mom_r[4] += om*phi_r[2]; 
    mom_r[5] += om*phi_r[3]; 
    mom_r[6] += om*phi_r[4]; 
  } 
} 
