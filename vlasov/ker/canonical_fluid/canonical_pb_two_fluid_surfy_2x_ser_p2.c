#include <gkyl_canonical_fluid_kernels.h>
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_two_fluid_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi, 
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

  const double *f1l = &fl[0]; 
  const double *f2l = &fl[8]; 
  const double *f1c = &fc[0]; 
  const double *f2c = &fc[8]; 
  const double *f1r = &fr[0]; 
  const double *f2r = &fr[8]; 
  double *out1 = &out[0]; 
  double *out2 = &out[8]; 
  double f1UpL[3] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  f1UpL[0] = 1.5811388300841895*f1l[5]+1.224744871391589*f1l[2]+0.7071067811865475*f1l[0]; 
  f1UpL[1] = 1.5811388300841898*f1l[7]+1.224744871391589*f1l[3]+0.7071067811865475*f1l[1]; 
  f1UpL[2] = 1.224744871391589*f1l[6]+0.7071067811865475*f1l[4]; 
    } else { 
  f1UpL[0] = 1.5811388300841895*f1c[5]-1.224744871391589*f1c[2]+0.7071067811865475*f1c[0]; 
  f1UpL[1] = 1.5811388300841898*f1c[7]-1.224744871391589*f1c[3]+0.7071067811865475*f1c[1]; 
  f1UpL[2] = 0.7071067811865475*f1c[4]-1.224744871391589*f1c[6]; 
    } 
  } else { 
  double f_lr[3] = {0.};
  double f_cl[3] = {0.};
  double sgn_alphaUpL[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.5811388300841895*f1l[5]+1.224744871391589*f1l[2]+0.7071067811865475*f1l[0]; 
  f_lr[1] = 1.5811388300841898*f1l[7]+1.224744871391589*f1l[3]+0.7071067811865475*f1l[1]; 
  f_lr[2] = 1.224744871391589*f1l[6]+0.7071067811865475*f1l[4]; 

  f_cl[0] = 1.5811388300841895*f1c[5]-1.224744871391589*f1c[2]+0.7071067811865475*f1c[0]; 
  f_cl[1] = 1.5811388300841898*f1c[7]-1.224744871391589*f1c[3]+0.7071067811865475*f1c[1]; 
  f_cl[2] = 0.7071067811865475*f1c[4]-1.224744871391589*f1c[6]; 

  f1UpL[0] = (0.3535533905932737*f_lr[2]-0.3535533905932737*f_cl[2])*sgn_alphaUpL[2]+(0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  f1UpL[1] = (0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.3162277660168379*f_lr[2]-0.3162277660168379*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 
  f1UpL[2] = (0.22587697572631277*f_lr[2]-0.22587697572631277*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[2]+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[2]+(0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[1]; 

  } 
  double f1UpR[3] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  f1UpR[0] = 1.5811388300841895*f1c[5]+1.224744871391589*f1c[2]+0.7071067811865475*f1c[0]; 
  f1UpR[1] = 1.5811388300841898*f1c[7]+1.224744871391589*f1c[3]+0.7071067811865475*f1c[1]; 
  f1UpR[2] = 1.224744871391589*f1c[6]+0.7071067811865475*f1c[4]; 
    } else { 
  f1UpR[0] = 1.5811388300841895*f1r[5]-1.224744871391589*f1r[2]+0.7071067811865475*f1r[0]; 
  f1UpR[1] = 1.5811388300841898*f1r[7]-1.224744871391589*f1r[3]+0.7071067811865475*f1r[1]; 
  f1UpR[2] = 0.7071067811865475*f1r[4]-1.224744871391589*f1r[6]; 
    } 
  } else { 
  double f_cr[3] = {0.};
  double f_rl[3] = {0.};
  double sgn_alphaUpR[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.5811388300841895*f1c[5]+1.224744871391589*f1c[2]+0.7071067811865475*f1c[0]; 
  f_cr[1] = 1.5811388300841898*f1c[7]+1.224744871391589*f1c[3]+0.7071067811865475*f1c[1]; 
  f_cr[2] = 1.224744871391589*f1c[6]+0.7071067811865475*f1c[4]; 

  f_rl[0] = 1.5811388300841895*f1r[5]-1.224744871391589*f1r[2]+0.7071067811865475*f1r[0]; 
  f_rl[1] = 1.5811388300841898*f1r[7]-1.224744871391589*f1r[3]+0.7071067811865475*f1r[1]; 
  f_rl[2] = 0.7071067811865475*f1r[4]-1.224744871391589*f1r[6]; 

  f1UpR[0] = (0.3535533905932737*f_cr[2]-0.3535533905932737*f_rl[2])*sgn_alphaUpR[2]+(0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  f1UpR[1] = (0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(-(0.3162277660168379*f_rl[2])+0.3162277660168379*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  f1UpR[2] = (-(0.22587697572631277*f_rl[2])+0.22587697572631277*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])*sgn_alphaUpR[2]+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[2]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[2]+(0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[1]; 

  } 
  double Ghat1L[3] = {0.};
  double Ghat1R[3] = {0.};
  Ghat1L[0] = 0.7071067811865475*(alphaL[1]*f1UpL[1]+alphaL[0]*f1UpL[0]); 
  Ghat1L[1] = 0.6324555320336759*alphaL[1]*f1UpL[2]+0.7071067811865475*(alphaL[0]*f1UpL[1]+f1UpL[0]*alphaL[1]); 
  Ghat1L[2] = 0.7071067811865475*alphaL[0]*f1UpL[2]+0.6324555320336759*alphaL[1]*f1UpL[1]; 

  Ghat1R[0] = 0.7071067811865475*(alphaR[1]*f1UpR[1]+alphaR[0]*f1UpR[0]); 
  Ghat1R[1] = 0.6324555320336759*alphaR[1]*f1UpR[2]+0.7071067811865475*(alphaR[0]*f1UpR[1]+f1UpR[0]*alphaR[1]); 
  Ghat1R[2] = 0.7071067811865475*alphaR[0]*f1UpR[2]+0.6324555320336759*alphaR[1]*f1UpR[1]; 

  out1[0] += (0.7071067811865475*Ghat1L[0]-0.7071067811865475*Ghat1R[0])*rdy2; 
  out1[1] += (0.7071067811865475*Ghat1L[1]-0.7071067811865475*Ghat1R[1])*rdy2; 
  out1[2] += (-(1.224744871391589*Ghat1R[0])-1.224744871391589*Ghat1L[0])*rdy2; 
  out1[3] += (-(1.224744871391589*Ghat1R[1])-1.224744871391589*Ghat1L[1])*rdy2; 
  out1[4] += (0.7071067811865475*Ghat1L[2]-0.7071067811865475*Ghat1R[2])*rdy2; 
  out1[5] += (1.5811388300841895*Ghat1L[0]-1.5811388300841895*Ghat1R[0])*rdy2; 
  out1[6] += (-(1.224744871391589*Ghat1R[2])-1.224744871391589*Ghat1L[2])*rdy2; 
  out1[7] += (1.5811388300841898*Ghat1L[1]-1.5811388300841898*Ghat1R[1])*rdy2; 

  double f2UpL[3] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  f2UpL[0] = 1.5811388300841895*f2l[5]+1.224744871391589*f2l[2]+0.7071067811865475*f2l[0]; 
  f2UpL[1] = 1.5811388300841898*f2l[7]+1.224744871391589*f2l[3]+0.7071067811865475*f2l[1]; 
  f2UpL[2] = 1.224744871391589*f2l[6]+0.7071067811865475*f2l[4]; 
    } else { 
  f2UpL[0] = 1.5811388300841895*f2c[5]-1.224744871391589*f2c[2]+0.7071067811865475*f2c[0]; 
  f2UpL[1] = 1.5811388300841898*f2c[7]-1.224744871391589*f2c[3]+0.7071067811865475*f2c[1]; 
  f2UpL[2] = 0.7071067811865475*f2c[4]-1.224744871391589*f2c[6]; 
    } 
  } else { 
  double f_lr[3] = {0.};
  double f_cl[3] = {0.};
  double sgn_alphaUpL[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.5811388300841895*f2l[5]+1.224744871391589*f2l[2]+0.7071067811865475*f2l[0]; 
  f_lr[1] = 1.5811388300841898*f2l[7]+1.224744871391589*f2l[3]+0.7071067811865475*f2l[1]; 
  f_lr[2] = 1.224744871391589*f2l[6]+0.7071067811865475*f2l[4]; 

  f_cl[0] = 1.5811388300841895*f2c[5]-1.224744871391589*f2c[2]+0.7071067811865475*f2c[0]; 
  f_cl[1] = 1.5811388300841898*f2c[7]-1.224744871391589*f2c[3]+0.7071067811865475*f2c[1]; 
  f_cl[2] = 0.7071067811865475*f2c[4]-1.224744871391589*f2c[6]; 

  f2UpL[0] = (0.3535533905932737*f_lr[2]-0.3535533905932737*f_cl[2])*sgn_alphaUpL[2]+(0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  f2UpL[1] = (0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.3162277660168379*f_lr[2]-0.3162277660168379*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 
  f2UpL[2] = (0.22587697572631277*f_lr[2]-0.22587697572631277*f_cl[2]+0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[2]+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[2]+(0.3162277660168379*f_lr[1]-0.3162277660168379*f_cl[1])*sgn_alphaUpL[1]; 

  } 
  double f2UpR[3] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  f2UpR[0] = 1.5811388300841895*f2c[5]+1.224744871391589*f2c[2]+0.7071067811865475*f2c[0]; 
  f2UpR[1] = 1.5811388300841898*f2c[7]+1.224744871391589*f2c[3]+0.7071067811865475*f2c[1]; 
  f2UpR[2] = 1.224744871391589*f2c[6]+0.7071067811865475*f2c[4]; 
    } else { 
  f2UpR[0] = 1.5811388300841895*f2r[5]-1.224744871391589*f2r[2]+0.7071067811865475*f2r[0]; 
  f2UpR[1] = 1.5811388300841898*f2r[7]-1.224744871391589*f2r[3]+0.7071067811865475*f2r[1]; 
  f2UpR[2] = 0.7071067811865475*f2r[4]-1.224744871391589*f2r[6]; 
    } 
  } else { 
  double f_cr[3] = {0.};
  double f_rl[3] = {0.};
  double sgn_alphaUpR[3] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.5811388300841895*f2c[5]+1.224744871391589*f2c[2]+0.7071067811865475*f2c[0]; 
  f_cr[1] = 1.5811388300841898*f2c[7]+1.224744871391589*f2c[3]+0.7071067811865475*f2c[1]; 
  f_cr[2] = 1.224744871391589*f2c[6]+0.7071067811865475*f2c[4]; 

  f_rl[0] = 1.5811388300841895*f2r[5]-1.224744871391589*f2r[2]+0.7071067811865475*f2r[0]; 
  f_rl[1] = 1.5811388300841898*f2r[7]-1.224744871391589*f2r[3]+0.7071067811865475*f2r[1]; 
  f_rl[2] = 0.7071067811865475*f2r[4]-1.224744871391589*f2r[6]; 

  f2UpR[0] = (0.3535533905932737*f_cr[2]-0.3535533905932737*f_rl[2])*sgn_alphaUpR[2]+(0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  f2UpR[1] = (0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(-(0.3162277660168379*f_rl[2])+0.3162277660168379*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  f2UpR[2] = (-(0.22587697572631277*f_rl[2])+0.22587697572631277*f_cr[2]-0.3535533905932737*f_rl[0]+0.3535533905932737*f_cr[0])*sgn_alphaUpR[2]+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[2]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[2]+(0.3162277660168379*f_cr[1]-0.3162277660168379*f_rl[1])*sgn_alphaUpR[1]; 

  } 
  double Ghat2L[3] = {0.};
  double Ghat2R[3] = {0.};
  Ghat2L[0] = 0.7071067811865475*(alphaL[1]*f2UpL[1]+alphaL[0]*f2UpL[0]); 
  Ghat2L[1] = 0.6324555320336759*alphaL[1]*f2UpL[2]+0.7071067811865475*(alphaL[0]*f2UpL[1]+f2UpL[0]*alphaL[1]); 
  Ghat2L[2] = 0.7071067811865475*alphaL[0]*f2UpL[2]+0.6324555320336759*alphaL[1]*f2UpL[1]; 

  Ghat2R[0] = 0.7071067811865475*(alphaR[1]*f2UpR[1]+alphaR[0]*f2UpR[0]); 
  Ghat2R[1] = 0.6324555320336759*alphaR[1]*f2UpR[2]+0.7071067811865475*(alphaR[0]*f2UpR[1]+f2UpR[0]*alphaR[1]); 
  Ghat2R[2] = 0.7071067811865475*alphaR[0]*f2UpR[2]+0.6324555320336759*alphaR[1]*f2UpR[1]; 

  out2[0] += (0.7071067811865475*Ghat2L[0]-0.7071067811865475*Ghat2R[0])*rdy2; 
  out2[1] += (0.7071067811865475*Ghat2L[1]-0.7071067811865475*Ghat2R[1])*rdy2; 
  out2[2] += (-(1.224744871391589*Ghat2R[0])-1.224744871391589*Ghat2L[0])*rdy2; 
  out2[3] += (-(1.224744871391589*Ghat2R[1])-1.224744871391589*Ghat2L[1])*rdy2; 
  out2[4] += (0.7071067811865475*Ghat2L[2]-0.7071067811865475*Ghat2R[2])*rdy2; 
  out2[5] += (1.5811388300841895*Ghat2L[0]-1.5811388300841895*Ghat2R[0])*rdy2; 
  out2[6] += (-(1.224744871391589*Ghat2R[2])-1.224744871391589*Ghat2L[2])*rdy2; 
  out2[7] += (1.5811388300841898*Ghat2L[1]-1.5811388300841898*Ghat2R[1])*rdy2; 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 1.7677669529663687*rdy2*cflFreq; 

} 
