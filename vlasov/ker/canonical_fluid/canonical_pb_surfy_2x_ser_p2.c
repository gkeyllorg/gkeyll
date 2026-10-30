#include <gkyl_canonical_fluid_kernels.h>
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points.
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points.
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1.
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1.
  // fl,fc,fr: input state vector in left, center and right cells.
  // out: output increment in center cell.

  double rdy2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_l[3];
  const double *alphaR = &alpha_surf_r[3];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_l[3];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_r[3];
  const int *const_sgn_alphaL = &const_sgn_alpha_l[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_r[1];

  double fUpL[3] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.5811388300841895*fl[5]+1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  fUpL[1] = 1.5811388300841898*fl[7]+1.224744871391589*fl[3]+0.7071067811865475*fl[1]; 
  fUpL[2] = 1.224744871391589*fl[6]+0.7071067811865475*fl[4]; 
    } else { 
  fUpL[0] = 1.5811388300841895*fc[5]-1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fUpL[1] = 1.5811388300841898*fc[7]-1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  fUpL[2] = 0.7071067811865475*fc[4]-1.224744871391589*fc[6]; 
    } 
  } else { 
  double f_lr[3] = {0.};
  double f_cl[3] = {0.};
  double sgn_alphaUpL[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.5811388300841895*fl[5]+1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  f_lr[1] = 1.5811388300841898*fl[7]+1.224744871391589*fl[3]+0.7071067811865475*fl[1]; 
  f_lr[2] = 1.224744871391589*fl[6]+0.7071067811865475*fl[4]; 

  f_cl[0] = 1.5811388300841895*fc[5]-1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  f_cl[1] = 1.5811388300841898*fc[7]-1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  f_cl[2] = 0.7071067811865475*fc[4]-1.224744871391589*fc[6]; 

  fUpL[0] = (0.3535533905932737*f_lr[2]-0.3535533905932737*f_cl[2])*sgn_alphaUpL[2]+(0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.3162277660168379*f_lr[2]-0.3162277660168379*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.22587697572631277*f_lr[2]-0.22587697572631277*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[2]+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[2]+(0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[1]; 

  } 
  double fUpR[3] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.5811388300841895*fc[5]+1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fUpR[1] = 1.5811388300841898*fc[7]+1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  fUpR[2] = 1.224744871391589*fc[6]+0.7071067811865475*fc[4]; 
    } else { 
  fUpR[0] = 1.5811388300841895*fr[5]-1.224744871391589*fr[2]+0.7071067811865475*fr[0]; 
  fUpR[1] = 1.5811388300841898*fr[7]-1.224744871391589*fr[3]+0.7071067811865475*fr[1]; 
  fUpR[2] = 0.7071067811865475*fr[4]-1.224744871391589*fr[6]; 
    } 
  } else { 
  double f_cr[3] = {0.};
  double f_rl[3] = {0.};
  double sgn_alphaUpR[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.5811388300841895*fc[5]+1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  f_cr[1] = 1.5811388300841898*fc[7]+1.224744871391589*fc[3]+0.7071067811865475*fc[1]; 
  f_cr[2] = 1.224744871391589*fc[6]+0.7071067811865475*fc[4]; 

  f_rl[0] = 1.5811388300841895*fr[5]-1.224744871391589*fr[2]+0.7071067811865475*fr[0]; 
  f_rl[1] = 1.5811388300841898*fr[7]-1.224744871391589*fr[3]+0.7071067811865475*fr[1]; 
  f_rl[2] = 0.7071067811865475*fr[4]-1.224744871391589*fr[6]; 

  fUpR[0] = (0.3535533905932737*f_cr[2]-0.3535533905932737*f_rl[2])*sgn_alphaUpR[2]+(0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(-(0.3162277660168379*f_rl[2])+0.3162277660168379*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (-(0.22587697572631277*f_rl[2])+0.22587697572631277*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])*sgn_alphaUpR[2]+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[2]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[2]+(0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[1]; 

  } 
  double GhatL[3] = {0.};
  double GhatR[3] = {0.};
  GhatL[0] = 0.7071067811865475*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.6324555320336759*alphaL[1]*fUpL[2]+0.7071067811865475*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.7071067811865475*alphaL[0]*fUpL[2]+0.6324555320336759*alphaL[1]*fUpL[1]; 

  GhatR[0] = 0.7071067811865475*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.6324555320336759*alphaR[1]*fUpR[2]+0.7071067811865475*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.7071067811865475*alphaR[0]*fUpR[2]+0.6324555320336759*alphaR[1]*fUpR[1]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdy2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdy2; 
  out[2] += (-(1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdy2; 
  out[3] += (-(1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdy2; 
  out[4] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdy2; 
  out[5] += (1.5811388300841895*GhatL[0]-1.5811388300841895*GhatR[0])*rdy2; 
  out[6] += (-(1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdy2; 
  out[7] += (1.5811388300841898*GhatL[1]-1.5811388300841898*GhatR[1])*rdy2; 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 1.7677669529663687*rdy2*cflFreq; 

} 
