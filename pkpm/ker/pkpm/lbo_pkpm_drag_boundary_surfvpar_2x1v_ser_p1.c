#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[12]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf[4] = {0.0}; 
  double F_0_UpwindQuad[4] = {0.0};
  double F_0_Upwind[4] = {0.0};;
  double Ghat_F_0[4] = {0.0}; 
  double G_1_UpwindQuad[4] = {0.0};
  double G_1_Upwind[4] = {0.0};;
  double Ghat_G_1[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar+0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  if (0.5*alphaDrSurf[3]-0.5*(alphaDrSurf[2]+alphaDrSurf[1])+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Edge); 
  } 
  if (-(0.5*alphaDrSurf[3])+0.5*alphaDrSurf[2]-0.5*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Edge); 
  } 
  if (0.5*(alphaDrSurf[1]+alphaDrSurf[0])-0.5*(alphaDrSurf[3]+alphaDrSurf[2]) < 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Edge); 
  } 
  if (0.5*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Skin); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Edge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaDrSurf[3]+0.5*F_0_Upwind[2]*alphaDrSurf[2]+0.5*F_0_Upwind[1]*alphaDrSurf[1]+0.5*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaDrSurf[3]+0.5*G_1_Upwind[2]*alphaDrSurf[2]+0.5*G_1_Upwind[1]*alphaDrSurf[1]+0.5*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*G_1_Upwind[2]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += 1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += 1.5811388300841898*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += 1.5811388300841895*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += 1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += 1.5811388300841898*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += 1.5811388300841895*Ghat_G_1[3]*dv1par; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar-0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  if (0.5*alphaDrSurf[3]-0.5*(alphaDrSurf[2]+alphaDrSurf[1])+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Skin); 
  } 
  if (-(0.5*alphaDrSurf[3])+0.5*alphaDrSurf[2]-0.5*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Skin); 
  } 
  if (0.5*(alphaDrSurf[1]+alphaDrSurf[0])-0.5*(alphaDrSurf[3]+alphaDrSurf[2]) < 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Skin); 
  } 
  if (0.5*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Edge); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Skin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaDrSurf[3]+0.5*F_0_Upwind[2]*alphaDrSurf[2]+0.5*F_0_Upwind[1]*alphaDrSurf[1]+0.5*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaDrSurf[3]+0.5*G_1_Upwind[2]*alphaDrSurf[2]+0.5*G_1_Upwind[1]*alphaDrSurf[1]+0.5*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*G_1_Upwind[2]; 

  out_F_0[0] += -(0.7071067811865475*Ghat_F_0[0]*dv1par); 
  out_F_0[1] += -(0.7071067811865475*Ghat_F_0[1]*dv1par); 
  out_F_0[2] += -(0.7071067811865475*Ghat_F_0[2]*dv1par); 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += -(0.7071067811865475*Ghat_F_0[3]*dv1par); 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += -(1.5811388300841895*Ghat_F_0[0]*dv1par); 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += -(1.5811388300841898*Ghat_F_0[1]*dv1par); 
  out_F_0[10] += -(1.5811388300841898*Ghat_F_0[2]*dv1par); 
  out_F_0[11] += -(1.5811388300841895*Ghat_F_0[3]*dv1par); 
  out_G_1[0] += -(0.7071067811865475*Ghat_G_1[0]*dv1par); 
  out_G_1[1] += -(0.7071067811865475*Ghat_G_1[1]*dv1par); 
  out_G_1[2] += -(0.7071067811865475*Ghat_G_1[2]*dv1par); 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += -(0.7071067811865475*Ghat_G_1[3]*dv1par); 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += -(1.5811388300841895*Ghat_G_1[0]*dv1par); 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += -(1.5811388300841898*Ghat_G_1[1]*dv1par); 
  out_G_1[10] += -(1.5811388300841898*Ghat_G_1[2]*dv1par); 
  out_G_1[11] += -(1.5811388300841895*Ghat_G_1[3]*dv1par); 

  }

  return 0.;

} 
