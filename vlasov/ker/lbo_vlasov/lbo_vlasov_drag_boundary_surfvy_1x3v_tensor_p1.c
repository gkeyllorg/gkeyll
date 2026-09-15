#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuPrimMomsSum[2]; 

  double alphaDrSurf[18] = {0.0}; 
  double fUpwindQuad[18] = {0.0};
  double fUpwind[18] = {0.0};;
  double Ghat[18] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUy[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_r(fSkin); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_l(fEdge); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+alphaDrSurf[1]*fUpwind[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[0]*fUpwind[6]); 
  Ghat[7] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[10]+15.0*alphaDrSurf[0]*fUpwind[7]); 
  Ghat[8] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[12]+15.0*alphaDrSurf[0]*fUpwind[8]); 
  Ghat[9] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[1]*fUpwind[6]); 
  Ghat[10] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[10]+15.000000000000002*alphaDrSurf[1]*fUpwind[7]); 
  Ghat[11] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[14]+15.0*alphaDrSurf[0]*fUpwind[11]); 
  Ghat[12] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[12]+15.000000000000002*alphaDrSurf[1]*fUpwind[8]); 
  Ghat[13] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[15]+15.0*alphaDrSurf[0]*fUpwind[13]); 
  Ghat[14] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[14]+15.000000000000002*alphaDrSurf[1]*fUpwind[11]); 
  Ghat[15] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[15]+15.000000000000002*alphaDrSurf[1]*fUpwind[13]); 
  Ghat[16] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[17]+alphaDrSurf[0]*fUpwind[16]); 
  Ghat[17] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[17]+alphaDrSurf[1]*fUpwind[16]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[14] += 1.224744871391589*Ghat[4]*rdv2; 
  out[15] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[16] += 1.224744871391589*Ghat[5]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[19] += 1.224744871391589*Ghat[7]*rdv2; 
  out[20] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[21] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[22] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[23] += 1.5811388300841898*Ghat[3]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[25] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += 1.224744871391589*Ghat[10]*rdv2; 
  out[29] += 1.5811388300841895*Ghat[4]*rdv2; 
  out[30] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[31] += 1.224744871391589*Ghat[11]*rdv2; 
  out[32] += 1.5811388300841895*Ghat[5]*rdv2; 
  out[33] += 1.5811388300841895*Ghat[6]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[35] += 1.224744871391589*Ghat[12]*rdv2; 
  out[36] += 1.224744871391589*Ghat[13]*rdv2; 
  out[37] += 1.5811388300841895*Ghat[7]*rdv2; 
  out[38] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[39] += 1.5811388300841895*Ghat[8]*rdv2; 
  out[40] += 1.224744871391589*Ghat[14]*rdv2; 
  out[41] += 1.5811388300841898*Ghat[9]*rdv2; 
  out[42] += 1.224744871391589*Ghat[15]*rdv2; 
  out[43] += 1.5811388300841898*Ghat[10]*rdv2; 
  out[44] += 1.5811388300841898*Ghat[11]*rdv2; 
  out[45] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[46] += 1.224744871391589*Ghat[16]*rdv2; 
  out[47] += 1.5811388300841898*Ghat[12]*rdv2; 
  out[48] += 1.5811388300841898*Ghat[13]*rdv2; 
  out[49] += 1.5811388300841895*Ghat[14]*rdv2; 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += 1.5811388300841895*Ghat[15]*rdv2; 
  out[52] += 1.5811388300841895*Ghat[16]*rdv2; 
  out[53] += 1.5811388300841898*Ghat[17]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUy[1]; 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx3_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx3_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_r(fEdge); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx3_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx3_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx3_eval_quad_node_11_l(fSkin); 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx3_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx3_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx3_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx3_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx3_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx3_eval_quad_node_17_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+alphaDrSurf[1]*fUpwind[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[0]*fUpwind[6]); 
  Ghat[7] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[10]+15.0*alphaDrSurf[0]*fUpwind[7]); 
  Ghat[8] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[12]+15.0*alphaDrSurf[0]*fUpwind[8]); 
  Ghat[9] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[1]*fUpwind[6]); 
  Ghat[10] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[10]+15.000000000000002*alphaDrSurf[1]*fUpwind[7]); 
  Ghat[11] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[14]+15.0*alphaDrSurf[0]*fUpwind[11]); 
  Ghat[12] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[12]+15.000000000000002*alphaDrSurf[1]*fUpwind[8]); 
  Ghat[13] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[15]+15.0*alphaDrSurf[0]*fUpwind[13]); 
  Ghat[14] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[14]+15.000000000000002*alphaDrSurf[1]*fUpwind[11]); 
  Ghat[15] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[15]+15.000000000000002*alphaDrSurf[1]*fUpwind[13]); 
  Ghat[16] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[17]+alphaDrSurf[0]*fUpwind[16]); 
  Ghat[17] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[17]+alphaDrSurf[1]*fUpwind[16]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[9] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[12] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[13] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[14] += 1.224744871391589*Ghat[4]*rdv2; 
  out[15] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[16] += 1.224744871391589*Ghat[5]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[19] += 1.224744871391589*Ghat[7]*rdv2; 
  out[20] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[21] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[22] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[23] += -(1.5811388300841898*Ghat[3]*rdv2); 
  out[24] += -(0.7071067811865475*Ghat[12]*rdv2); 
  out[25] += -(0.7071067811865475*Ghat[13]*rdv2); 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += 1.224744871391589*Ghat[10]*rdv2; 
  out[29] += -(1.5811388300841895*Ghat[4]*rdv2); 
  out[30] += -(0.7071067811865475*Ghat[14]*rdv2); 
  out[31] += 1.224744871391589*Ghat[11]*rdv2; 
  out[32] += -(1.5811388300841895*Ghat[5]*rdv2); 
  out[33] += -(1.5811388300841895*Ghat[6]*rdv2); 
  out[34] += -(0.7071067811865475*Ghat[15]*rdv2); 
  out[35] += 1.224744871391589*Ghat[12]*rdv2; 
  out[36] += 1.224744871391589*Ghat[13]*rdv2; 
  out[37] += -(1.5811388300841895*Ghat[7]*rdv2); 
  out[38] += -(0.7071067811865475*Ghat[16]*rdv2); 
  out[39] += -(1.5811388300841895*Ghat[8]*rdv2); 
  out[40] += 1.224744871391589*Ghat[14]*rdv2; 
  out[41] += -(1.5811388300841898*Ghat[9]*rdv2); 
  out[42] += 1.224744871391589*Ghat[15]*rdv2; 
  out[43] += -(1.5811388300841898*Ghat[10]*rdv2); 
  out[44] += -(1.5811388300841898*Ghat[11]*rdv2); 
  out[45] += -(0.7071067811865475*Ghat[17]*rdv2); 
  out[46] += 1.224744871391589*Ghat[16]*rdv2; 
  out[47] += -(1.5811388300841898*Ghat[12]*rdv2); 
  out[48] += -(1.5811388300841898*Ghat[13]*rdv2); 
  out[49] += -(1.5811388300841895*Ghat[14]*rdv2); 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += -(1.5811388300841895*Ghat[15]*rdv2); 
  out[52] += -(1.5811388300841895*Ghat[16]*rdv2); 
  out[53] += -(1.5811388300841898*Ghat[17]*rdv2); 

  } 

  return 0.;

} 
