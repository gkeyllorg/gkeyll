#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_3x1v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[24]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[24]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf[8] = {0.0}; 
  double F_0_UpwindQuad[8] = {0.0};
  double F_0_Upwind[8] = {0.0};;
  double Ghat_F_0[8] = {0.0}; 
  double G_1_UpwindQuad[8] = {0.0};
  double G_1_Upwind[8] = {0.0};;
  double Ghat_G_1[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar+0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 
  alphaDrSurf[4] = nuSum[4]*wvpar+0.5*nuSum[4]*dvpar-1.0*sumNuUPar[4]; 
  alphaDrSurf[5] = nuSum[5]*wvpar+0.5*nuSum[5]*dvpar-1.0*sumNuUPar[5]; 
  alphaDrSurf[6] = nuSum[6]*wvpar+0.5*nuSum[6]*dvpar-1.0*sumNuUPar[6]; 
  alphaDrSurf[7] = nuSum[7]*wvpar+0.5*nuSum[7]*dvpar-1.0*sumNuUPar[7]; 

  if (-(0.3535533905932737*alphaDrSurf[7])+0.3535533905932737*(alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[4])-0.3535533905932737*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1Edge); 
  } 
  if (0.3535533905932737*alphaDrSurf[7]-0.3535533905932737*(alphaDrSurf[6]+alphaDrSurf[5])+0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[3])-0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1Edge); 
  } 
  if (0.3535533905932737*alphaDrSurf[7]-0.3535533905932737*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]-0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[3])+0.3535533905932737*alphaDrSurf[2]-0.3535533905932737*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1Edge); 
  } 
  if (-(0.3535533905932737*alphaDrSurf[7])+0.3535533905932737*alphaDrSurf[6]-0.3535533905932737*(alphaDrSurf[5]+alphaDrSurf[4])+0.3535533905932737*(alphaDrSurf[3]+alphaDrSurf[2])-0.3535533905932737*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1Edge); 
  } 
  if (0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6])-0.3535533905932737*(alphaDrSurf[5]+alphaDrSurf[4]+alphaDrSurf[3]+alphaDrSurf[2])+0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0Skin); 
    G_1_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0Edge); 
    G_1_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1Edge); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]))+0.3535533905932737*alphaDrSurf[5]-0.3535533905932737*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]-0.3535533905932737*alphaDrSurf[2]+0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0Skin); 
    G_1_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0Edge); 
    G_1_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1Edge); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]))+0.3535533905932737*alphaDrSurf[4]-0.3535533905932737*alphaDrSurf[3]+0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0Skin); 
    G_1_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0Edge); 
    G_1_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1Edge); 
  } 
  if (0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[4]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0Skin); 
    G_1_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0Edge); 
    G_1_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1Edge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.3535533905932737*F_0_Upwind[7]*alphaDrSurf[7]+0.3535533905932737*F_0_Upwind[6]*alphaDrSurf[6]+0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[5]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[4]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[3]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[2]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[1]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.3535533905932737*F_0_Upwind[6]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[6]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[5]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[3]; 
  Ghat_F_0[4] = 0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_F_0[5] = 0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[3]; 
  Ghat_F_0[6] = 0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[3]; 
  Ghat_F_0[7] = 0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[4]; 
  Ghat_G_1[0] = 0.3535533905932737*G_1_Upwind[7]*alphaDrSurf[7]+0.3535533905932737*G_1_Upwind[6]*alphaDrSurf[6]+0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[5]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[4]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[3]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[2]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[1]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.3535533905932737*G_1_Upwind[6]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[6]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[5]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[3]; 
  Ghat_G_1[4] = 0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[2]; 
  Ghat_G_1[5] = 0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[3]; 
  Ghat_G_1[6] = 0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[3]; 
  Ghat_G_1[7] = 0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[4]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += 0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[4] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[5] += 0.7071067811865475*Ghat_F_0[4]*dv1par; 
  out_F_0[6] += 0.7071067811865475*Ghat_F_0[5]*dv1par; 
  out_F_0[7] += 0.7071067811865475*Ghat_F_0[6]*dv1par; 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[9] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[10] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[11] += 1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[12] += 0.7071067811865475*Ghat_F_0[7]*dv1par; 
  out_F_0[13] += 1.224744871391589*Ghat_F_0[4]*dv1par; 
  out_F_0[14] += 1.224744871391589*Ghat_F_0[5]*dv1par; 
  out_F_0[15] += 1.224744871391589*Ghat_F_0[6]*dv1par; 
  out_F_0[16] += 1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_F_0[17] += 1.5811388300841898*Ghat_F_0[2]*dv1par; 
  out_F_0[18] += 1.5811388300841898*Ghat_F_0[3]*dv1par; 
  out_F_0[19] += 1.224744871391589*Ghat_F_0[7]*dv1par; 
  out_F_0[20] += 1.5811388300841895*Ghat_F_0[4]*dv1par; 
  out_F_0[21] += 1.5811388300841895*Ghat_F_0[5]*dv1par; 
  out_F_0[22] += 1.5811388300841895*Ghat_F_0[6]*dv1par; 
  out_F_0[23] += 1.5811388300841898*Ghat_F_0[7]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += 0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[4] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[5] += 0.7071067811865475*Ghat_G_1[4]*dv1par; 
  out_G_1[6] += 0.7071067811865475*Ghat_G_1[5]*dv1par; 
  out_G_1[7] += 0.7071067811865475*Ghat_G_1[6]*dv1par; 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[9] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[10] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[11] += 1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[12] += 0.7071067811865475*Ghat_G_1[7]*dv1par; 
  out_G_1[13] += 1.224744871391589*Ghat_G_1[4]*dv1par; 
  out_G_1[14] += 1.224744871391589*Ghat_G_1[5]*dv1par; 
  out_G_1[15] += 1.224744871391589*Ghat_G_1[6]*dv1par; 
  out_G_1[16] += 1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_G_1[17] += 1.5811388300841898*Ghat_G_1[2]*dv1par; 
  out_G_1[18] += 1.5811388300841898*Ghat_G_1[3]*dv1par; 
  out_G_1[19] += 1.224744871391589*Ghat_G_1[7]*dv1par; 
  out_G_1[20] += 1.5811388300841895*Ghat_G_1[4]*dv1par; 
  out_G_1[21] += 1.5811388300841895*Ghat_G_1[5]*dv1par; 
  out_G_1[22] += 1.5811388300841895*Ghat_G_1[6]*dv1par; 
  out_G_1[23] += 1.5811388300841898*Ghat_G_1[7]*dv1par; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar-0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 
  alphaDrSurf[4] = nuSum[4]*wvpar-0.5*nuSum[4]*dvpar-1.0*sumNuUPar[4]; 
  alphaDrSurf[5] = nuSum[5]*wvpar-0.5*nuSum[5]*dvpar-1.0*sumNuUPar[5]; 
  alphaDrSurf[6] = nuSum[6]*wvpar-0.5*nuSum[6]*dvpar-1.0*sumNuUPar[6]; 
  alphaDrSurf[7] = nuSum[7]*wvpar-0.5*nuSum[7]*dvpar-1.0*sumNuUPar[7]; 

  if (-(0.3535533905932737*alphaDrSurf[7])+0.3535533905932737*(alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[4])-0.3535533905932737*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1Skin); 
  } 
  if (0.3535533905932737*alphaDrSurf[7]-0.3535533905932737*(alphaDrSurf[6]+alphaDrSurf[5])+0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[3])-0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1Skin); 
  } 
  if (0.3535533905932737*alphaDrSurf[7]-0.3535533905932737*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]-0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[3])+0.3535533905932737*alphaDrSurf[2]-0.3535533905932737*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1Skin); 
  } 
  if (-(0.3535533905932737*alphaDrSurf[7])+0.3535533905932737*alphaDrSurf[6]-0.3535533905932737*(alphaDrSurf[5]+alphaDrSurf[4])+0.3535533905932737*(alphaDrSurf[3]+alphaDrSurf[2])-0.3535533905932737*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1Skin); 
  } 
  if (0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6])-0.3535533905932737*(alphaDrSurf[5]+alphaDrSurf[4]+alphaDrSurf[3]+alphaDrSurf[2])+0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0Edge); 
    G_1_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0Skin); 
    G_1_UpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1Skin); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]))+0.3535533905932737*alphaDrSurf[5]-0.3535533905932737*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]-0.3535533905932737*alphaDrSurf[2]+0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0Edge); 
    G_1_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0Skin); 
    G_1_UpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1Skin); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]))+0.3535533905932737*alphaDrSurf[4]-0.3535533905932737*alphaDrSurf[3]+0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0Edge); 
    G_1_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0Skin); 
    G_1_UpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1Skin); 
  } 
  if (0.3535533905932737*(alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[5]+alphaDrSurf[4]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0Edge); 
    G_1_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0Skin); 
    G_1_UpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1Skin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.3535533905932737*F_0_Upwind[7]*alphaDrSurf[7]+0.3535533905932737*F_0_Upwind[6]*alphaDrSurf[6]+0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[5]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[4]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[3]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[2]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[1]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.3535533905932737*F_0_Upwind[6]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[6]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[5]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[3]; 
  Ghat_F_0[4] = 0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[5]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[4]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_F_0[5] = 0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[3]; 
  Ghat_F_0[6] = 0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[4]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[4]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[3]; 
  Ghat_F_0[7] = 0.3535533905932737*F_0_Upwind[0]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*F_0_Upwind[7]+0.3535533905932737*F_0_Upwind[1]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[1]*F_0_Upwind[6]+0.3535533905932737*F_0_Upwind[2]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[2]*F_0_Upwind[5]+0.3535533905932737*F_0_Upwind[3]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]*F_0_Upwind[4]; 
  Ghat_G_1[0] = 0.3535533905932737*G_1_Upwind[7]*alphaDrSurf[7]+0.3535533905932737*G_1_Upwind[6]*alphaDrSurf[6]+0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[5]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[4]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[3]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[2]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[1]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.3535533905932737*G_1_Upwind[6]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[6]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[5]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[3]; 
  Ghat_G_1[4] = 0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[5]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[5]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[4]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[2]; 
  Ghat_G_1[5] = 0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[3]; 
  Ghat_G_1[6] = 0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[4]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[4]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[3]; 
  Ghat_G_1[7] = 0.3535533905932737*G_1_Upwind[0]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*G_1_Upwind[7]+0.3535533905932737*G_1_Upwind[1]*alphaDrSurf[6]+0.3535533905932737*alphaDrSurf[1]*G_1_Upwind[6]+0.3535533905932737*G_1_Upwind[2]*alphaDrSurf[5]+0.3535533905932737*alphaDrSurf[2]*G_1_Upwind[5]+0.3535533905932737*G_1_Upwind[3]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[3]*G_1_Upwind[4]; 

  out_F_0[0] += -(0.7071067811865475*Ghat_F_0[0]*dv1par); 
  out_F_0[1] += -(0.7071067811865475*Ghat_F_0[1]*dv1par); 
  out_F_0[2] += -(0.7071067811865475*Ghat_F_0[2]*dv1par); 
  out_F_0[3] += -(0.7071067811865475*Ghat_F_0[3]*dv1par); 
  out_F_0[4] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[5] += -(0.7071067811865475*Ghat_F_0[4]*dv1par); 
  out_F_0[6] += -(0.7071067811865475*Ghat_F_0[5]*dv1par); 
  out_F_0[7] += -(0.7071067811865475*Ghat_F_0[6]*dv1par); 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[9] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[10] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[11] += -(1.5811388300841895*Ghat_F_0[0]*dv1par); 
  out_F_0[12] += -(0.7071067811865475*Ghat_F_0[7]*dv1par); 
  out_F_0[13] += 1.224744871391589*Ghat_F_0[4]*dv1par; 
  out_F_0[14] += 1.224744871391589*Ghat_F_0[5]*dv1par; 
  out_F_0[15] += 1.224744871391589*Ghat_F_0[6]*dv1par; 
  out_F_0[16] += -(1.5811388300841898*Ghat_F_0[1]*dv1par); 
  out_F_0[17] += -(1.5811388300841898*Ghat_F_0[2]*dv1par); 
  out_F_0[18] += -(1.5811388300841898*Ghat_F_0[3]*dv1par); 
  out_F_0[19] += 1.224744871391589*Ghat_F_0[7]*dv1par; 
  out_F_0[20] += -(1.5811388300841895*Ghat_F_0[4]*dv1par); 
  out_F_0[21] += -(1.5811388300841895*Ghat_F_0[5]*dv1par); 
  out_F_0[22] += -(1.5811388300841895*Ghat_F_0[6]*dv1par); 
  out_F_0[23] += -(1.5811388300841898*Ghat_F_0[7]*dv1par); 
  out_G_1[0] += -(0.7071067811865475*Ghat_G_1[0]*dv1par); 
  out_G_1[1] += -(0.7071067811865475*Ghat_G_1[1]*dv1par); 
  out_G_1[2] += -(0.7071067811865475*Ghat_G_1[2]*dv1par); 
  out_G_1[3] += -(0.7071067811865475*Ghat_G_1[3]*dv1par); 
  out_G_1[4] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[5] += -(0.7071067811865475*Ghat_G_1[4]*dv1par); 
  out_G_1[6] += -(0.7071067811865475*Ghat_G_1[5]*dv1par); 
  out_G_1[7] += -(0.7071067811865475*Ghat_G_1[6]*dv1par); 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[9] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[10] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[11] += -(1.5811388300841895*Ghat_G_1[0]*dv1par); 
  out_G_1[12] += -(0.7071067811865475*Ghat_G_1[7]*dv1par); 
  out_G_1[13] += 1.224744871391589*Ghat_G_1[4]*dv1par; 
  out_G_1[14] += 1.224744871391589*Ghat_G_1[5]*dv1par; 
  out_G_1[15] += 1.224744871391589*Ghat_G_1[6]*dv1par; 
  out_G_1[16] += -(1.5811388300841898*Ghat_G_1[1]*dv1par); 
  out_G_1[17] += -(1.5811388300841898*Ghat_G_1[2]*dv1par); 
  out_G_1[18] += -(1.5811388300841898*Ghat_G_1[3]*dv1par); 
  out_G_1[19] += 1.224744871391589*Ghat_G_1[7]*dv1par; 
  out_G_1[20] += -(1.5811388300841895*Ghat_G_1[4]*dv1par); 
  out_G_1[21] += -(1.5811388300841895*Ghat_G_1[5]*dv1par); 
  out_G_1[22] += -(1.5811388300841895*Ghat_G_1[6]*dv1par); 
  out_G_1[23] += -(1.5811388300841898*Ghat_G_1[7]*dv1par); 

  }

  return 0.;

} 
