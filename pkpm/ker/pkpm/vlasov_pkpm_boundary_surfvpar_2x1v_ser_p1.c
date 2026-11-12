#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_pkpm_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *div_b, const double *pkpm_accel_vars, 
     const double *g_dist_sourceEdge, const double *g_dist_sourceSkin, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates. 
  // dxv[NDIM]:              Cell spacing. 
  // div_b:                  Input volume expansion of div(b). 
  // pkpm_accel_vars:        Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // g_dist_sourceEdge/Skin: Input [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
  //                         (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]. 
  //                         in skin cell/last edge cell. First input is mirror force source, second input is vperp characteristics source. 
  // edge:                   Determines if the update is for the left edge (-1) or right edge (+1). 
  // fSkin/fEdge:            Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)] in skin cell/last edge cell. 
  // out:                    Incremented output distribution functions in center cell. 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[12]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[12]; 
  const double *F_0_sourceSkin = &fSkin[12]; 
  const double *G_1_sourceSkin = &g_dist_sourceSkin[0]; 
  const double *F_0_sourceEdge = &fEdge[12]; 
  const double *G_1_sourceEdge = &g_dist_sourceEdge[0]; 
  const double *bb_grad_u = &pkpm_accel_vars[4]; 
  const double *p_force = &pkpm_accel_vars[8]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  double div_b_Surf[4] = {0.0}; 
  div_b_Surf[0] = div_b[0]; 
  div_b_Surf[1] = div_b[1]; 
  div_b_Surf[2] = div_b[2]; 
  div_b_Surf[3] = div_b[3]; 

  double alphaSurf[4] = {0.0}; 
  double F_0_UpwindQuad[4] = {0.0};
  double F_0_Upwind[4] = {0.0};;
  double Ghat_F_0[4] = {0.0}; 
  double G_1_UpwindQuad[4] = {0.0};
  double G_1_Upwind[4] = {0.0};;
  double Ghat_G_1[4] = {0.0}; 
  double F_0_div_b_UpwindQuad[4] = {0.0};
  double F_0_div_b_Upwind[4] = {0.0};;
  double Ghat_F_0_div_b[4] = {0.0}; 
  double G_1_div_b_UpwindQuad[4] = {0.0};
  double G_1_div_b_Upwind[4] = {0.0};;
  double Ghat_G_1_div_b[4] = {0.0}; 

  // get stable timestep of alpha_v = 1/rho (div(p_par b) - p_perp div(b)) - v_par bb : grad(u) 
  // from the quadrature point evaluation needed to compute upwinded distribution functions 
  double cflFreq = 0.0;
  double alphaOrd = 0.0;

  if (edge == -1) { 

  alphaSurf[0] = -(1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = -(1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = -(1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = -(1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 

  alphaOrd = 0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Edge); 
  } 
  alphaOrd = -(0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (-(0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Edge); 
  } 
  alphaOrd = 0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]);
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]) > 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Edge); 
  } 
  alphaOrd = 0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]);
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Edge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  if (0.5*div_b_Surf[3]-0.5*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1_sourceEdge); 
  } 
  if (-(0.5*div_b_Surf[3])+0.5*div_b_Surf[2]-0.5*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1_sourceEdge); 
  } 
  if (0.5*(div_b_Surf[1]+div_b_Surf[0])-0.5*(div_b_Surf[3]+div_b_Surf[2]) > 0) { 
    F_0_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1_sourceEdge); 
  } 
  if (0.5*(div_b_Surf[3]+div_b_Surf[2]+div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1_sourceSkin); 
  } else { 
    F_0_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1_sourceEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad, F_0_div_b_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad, G_1_div_b_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaSurf[3]+0.5*F_0_Upwind[2]*alphaSurf[2]+0.5*F_0_Upwind[1]*alphaSurf[1]+0.5*F_0_Upwind[0]*alphaSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaSurf[3]+0.5*G_1_Upwind[2]*alphaSurf[2]+0.5*G_1_Upwind[1]*alphaSurf[1]+0.5*G_1_Upwind[0]*alphaSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*G_1_Upwind[2]; 
  Ghat_F_0_div_b[0] = 0.5*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b[1] = 0.5*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[1]; 
  Ghat_F_0_div_b[2] = 0.5*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[3] = 0.5*F_0_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[2]; 
  Ghat_G_1_div_b[0] = 0.5*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b[1] = 0.5*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[1]; 
  Ghat_G_1_div_b[2] = 0.5*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[3] = 0.5*G_1_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[2]; 

  out_F_0[0] += -(0.7071067811865475*Ghat_F_0_div_b[0]*dv1par)-0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += -(0.7071067811865475*Ghat_F_0_div_b[1]*dv1par)-0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += -(0.7071067811865475*Ghat_F_0_div_b[2]*dv1par)-0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += -(1.224744871391589*Ghat_F_0_div_b[0]*dv1par)-1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += -(0.7071067811865475*Ghat_F_0_div_b[3]*dv1par)-0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += -(1.224744871391589*Ghat_F_0_div_b[1]*dv1par)-1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += -(1.224744871391589*Ghat_F_0_div_b[2]*dv1par)-1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += -(1.5811388300841895*Ghat_F_0_div_b[0]*dv1par)-1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[8] += -(1.224744871391589*Ghat_F_0_div_b[3]*dv1par)-1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += -(1.5811388300841898*Ghat_F_0_div_b[1]*dv1par)-1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += -(1.5811388300841898*Ghat_F_0_div_b[2]*dv1par)-1.5811388300841898*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += -(1.5811388300841895*Ghat_F_0_div_b[3]*dv1par)-1.5811388300841895*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += -(0.7071067811865475*Ghat_G_1_div_b[0]*dv1par)-0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += -(0.7071067811865475*Ghat_G_1_div_b[1]*dv1par)-0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += -(0.7071067811865475*Ghat_G_1_div_b[2]*dv1par)-0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += -(1.224744871391589*Ghat_G_1_div_b[0]*dv1par)-1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += -(0.7071067811865475*Ghat_G_1_div_b[3]*dv1par)-0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += -(1.224744871391589*Ghat_G_1_div_b[1]*dv1par)-1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += -(1.224744871391589*Ghat_G_1_div_b[2]*dv1par)-1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += -(1.5811388300841895*Ghat_G_1_div_b[0]*dv1par)-1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[8] += -(1.224744871391589*Ghat_G_1_div_b[3]*dv1par)-1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += -(1.5811388300841898*Ghat_G_1_div_b[1]*dv1par)-1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += -(1.5811388300841898*Ghat_G_1_div_b[2]*dv1par)-1.5811388300841898*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += -(1.5811388300841895*Ghat_G_1_div_b[3]*dv1par)-1.5811388300841895*Ghat_G_1[3]*dv1par; 

  } else { 

  alphaSurf[0] = -(1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = -(1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = -(1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = -(1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 

  alphaOrd = 0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Skin); 
  } 
  alphaOrd = -(0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0];
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (-(0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Skin); 
  } 
  alphaOrd = 0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]);
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]) > 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Skin); 
  } 
  alphaOrd = 0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]);
  cflFreq = fmax(cflFreq, fabs(alphaOrd));
  if (0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Skin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  if (0.5*div_b_Surf[3]-0.5*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1_sourceSkin); 
  } 
  if (-(0.5*div_b_Surf[3])+0.5*div_b_Surf[2]-0.5*div_b_Surf[1]+0.5*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1_sourceSkin); 
  } 
  if (0.5*(div_b_Surf[1]+div_b_Surf[0])-0.5*(div_b_Surf[3]+div_b_Surf[2]) > 0) { 
    F_0_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1_sourceSkin); 
  } 
  if (0.5*(div_b_Surf[3]+div_b_Surf[2]+div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0_sourceEdge); 
    G_1_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1_sourceEdge); 
  } else { 
    F_0_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0_sourceSkin); 
    G_1_div_b_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1_sourceSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad, F_0_div_b_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad, G_1_div_b_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaSurf[3]+0.5*F_0_Upwind[2]*alphaSurf[2]+0.5*F_0_Upwind[1]*alphaSurf[1]+0.5*F_0_Upwind[0]*alphaSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaSurf[3]+0.5*G_1_Upwind[2]*alphaSurf[2]+0.5*G_1_Upwind[1]*alphaSurf[1]+0.5*G_1_Upwind[0]*alphaSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaSurf[3]+0.5*alphaSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[1]+0.5*alphaSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaSurf[3]+0.5*alphaSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaSurf[2]+0.5*alphaSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaSurf[3]+0.5*alphaSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaSurf[2]+0.5*alphaSurf[1]*G_1_Upwind[2]; 
  Ghat_F_0_div_b[0] = 0.5*F_0_div_b_Upwind[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b[1] = 0.5*F_0_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[1]; 
  Ghat_F_0_div_b[2] = 0.5*F_0_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[2]; 
  Ghat_F_0_div_b[3] = 0.5*F_0_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind[3]+0.5*F_0_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind[2]; 
  Ghat_G_1_div_b[0] = 0.5*G_1_div_b_Upwind[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b[1] = 0.5*G_1_div_b_Upwind[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[1]; 
  Ghat_G_1_div_b[2] = 0.5*G_1_div_b_Upwind[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[2]; 
  Ghat_G_1_div_b[3] = 0.5*G_1_div_b_Upwind[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind[3]+0.5*G_1_div_b_Upwind[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind[2]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0_div_b[0]*dv1par+0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0_div_b[1]*dv1par+0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0_div_b[2]*dv1par+0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += -(1.224744871391589*Ghat_F_0_div_b[0]*dv1par)-1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0_div_b[3]*dv1par+0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += -(1.224744871391589*Ghat_F_0_div_b[1]*dv1par)-1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += -(1.224744871391589*Ghat_F_0_div_b[2]*dv1par)-1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.5811388300841895*Ghat_F_0_div_b[0]*dv1par+1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[8] += -(1.224744871391589*Ghat_F_0_div_b[3]*dv1par)-1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += 1.5811388300841898*Ghat_F_0_div_b[1]*dv1par+1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += 1.5811388300841898*Ghat_F_0_div_b[2]*dv1par+1.5811388300841898*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += 1.5811388300841895*Ghat_F_0_div_b[3]*dv1par+1.5811388300841895*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1_div_b[0]*dv1par+0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1_div_b[1]*dv1par+0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1_div_b[2]*dv1par+0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += -(1.224744871391589*Ghat_G_1_div_b[0]*dv1par)-1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1_div_b[3]*dv1par+0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += -(1.224744871391589*Ghat_G_1_div_b[1]*dv1par)-1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += -(1.224744871391589*Ghat_G_1_div_b[2]*dv1par)-1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.5811388300841895*Ghat_G_1_div_b[0]*dv1par+1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[8] += -(1.224744871391589*Ghat_G_1_div_b[3]*dv1par)-1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += 1.5811388300841898*Ghat_G_1_div_b[1]*dv1par+1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += 1.5811388300841898*Ghat_G_1_div_b[2]*dv1par+1.5811388300841898*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += 1.5811388300841895*Ghat_G_1_div_b[3]*dv1par+1.5811388300841895*Ghat_G_1[3]*dv1par; 

  } 

  return 2.5*dv1par*cflFreq;

} 
