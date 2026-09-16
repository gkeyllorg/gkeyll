#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]: Cell-center coordinates. 
  // dxv[5]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[36] = {0.0}; 
  double fUpwindQuad[36] = {0.0};
  double fUpwind[36] = {0.0};;
  double Ghat[36] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[1]; 
  alphaDrSurf[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]+dxv[2]*nuSum[2]; 
  alphaDrSurf[5] = (2.0*w[2]+dxv[2])*nuSum[3]-2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (-alphaDrSurf[5]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fSkin); 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fEdge); 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fEdge); 
  } 
  if (-alphaDrSurf[5]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fSkin); 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fSkin); 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fEdge); 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fEdge); 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fSkin); 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fSkin); 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fSkin); 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fSkin); 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fSkin); 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fSkin); 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fSkin); 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fSkin); 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fEdge); 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fEdge); 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fEdge); 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fEdge); 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fEdge); 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fEdge); 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fEdge); 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fEdge); 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alphaDrSurf[5]*fUpwind[5]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alphaDrSurf[2]*fUpwind[5]+fUpwind[2]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.25*(alphaDrSurf[1]*fUpwind[5]+fUpwind[1]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.25*(alphaDrSurf[5]*fUpwind[13]+alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.25*(alphaDrSurf[5]*fUpwind[14]+alphaDrSurf[2]*fUpwind[9]+alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.25*(alphaDrSurf[0]*fUpwind[5]+fUpwind[0]*alphaDrSurf[5]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[6] = 0.25*(alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[5]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[7] = 0.25*(alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[5]*fUpwind[6]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[8] = 0.25*(alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[5]*fUpwind[9]+alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[9] = 0.25*(alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[5]*fUpwind[8]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[10] = 0.25*(alphaDrSurf[5]*fUpwind[23]+alphaDrSurf[2]*fUpwind[16]+alphaDrSurf[1]*fUpwind[15]+alphaDrSurf[0]*fUpwind[10]); 
  Ghat[11] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[24]+15.000000000000002*(alphaDrSurf[2]*fUpwind[18]+alphaDrSurf[1]*fUpwind[17])+15.0*alphaDrSurf[0]*fUpwind[11]); 
  Ghat[12] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[27]+15.000000000000002*(alphaDrSurf[2]*fUpwind[21]+alphaDrSurf[1]*fUpwind[20])+15.0*alphaDrSurf[0]*fUpwind[12]); 
  Ghat[13] = 0.25*(alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[3]*alphaDrSurf[5]); 
  Ghat[14] = 0.25*(alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[2]*fUpwind[8]+fUpwind[4]*alphaDrSurf[5]); 
  Ghat[15] = 0.25*(alphaDrSurf[2]*fUpwind[23]+alphaDrSurf[5]*fUpwind[16]+alphaDrSurf[0]*fUpwind[15]+alphaDrSurf[1]*fUpwind[10]); 
  Ghat[16] = 0.25*(alphaDrSurf[1]*fUpwind[23]+alphaDrSurf[0]*fUpwind[16]+alphaDrSurf[5]*fUpwind[15]+alphaDrSurf[2]*fUpwind[10]); 
  Ghat[17] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[24]+15.0*(alphaDrSurf[5]*fUpwind[18]+alphaDrSurf[0]*fUpwind[17])+15.000000000000002*alphaDrSurf[1]*fUpwind[11]); 
  Ghat[18] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[24]+15.0*(alphaDrSurf[0]*fUpwind[18]+alphaDrSurf[5]*fUpwind[17])+15.000000000000002*alphaDrSurf[2]*fUpwind[11]); 
  Ghat[19] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[31]+15.000000000000002*(alphaDrSurf[2]*fUpwind[26]+alphaDrSurf[1]*fUpwind[25])+15.0*alphaDrSurf[0]*fUpwind[19]); 
  Ghat[20] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[27]+15.0*(alphaDrSurf[5]*fUpwind[21]+alphaDrSurf[0]*fUpwind[20])+15.000000000000002*alphaDrSurf[1]*fUpwind[12]); 
  Ghat[21] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[27]+15.0*(alphaDrSurf[0]*fUpwind[21]+alphaDrSurf[5]*fUpwind[20])+15.000000000000002*alphaDrSurf[2]*fUpwind[12]); 
  Ghat[22] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[32]+15.000000000000002*(alphaDrSurf[2]*fUpwind[29]+alphaDrSurf[1]*fUpwind[28])+15.0*alphaDrSurf[0]*fUpwind[22]); 
  Ghat[23] = 0.25*(alphaDrSurf[0]*fUpwind[23]+alphaDrSurf[1]*fUpwind[16]+alphaDrSurf[2]*fUpwind[15]+alphaDrSurf[5]*fUpwind[10]); 
  Ghat[24] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[24]+15.000000000000002*(alphaDrSurf[1]*fUpwind[18]+alphaDrSurf[2]*fUpwind[17])+15.0*alphaDrSurf[5]*fUpwind[11]); 
  Ghat[25] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[31]+15.0*(alphaDrSurf[5]*fUpwind[26]+alphaDrSurf[0]*fUpwind[25])+15.000000000000002*alphaDrSurf[1]*fUpwind[19]); 
  Ghat[26] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[31]+15.0*(alphaDrSurf[0]*fUpwind[26]+alphaDrSurf[5]*fUpwind[25])+15.000000000000002*alphaDrSurf[2]*fUpwind[19]); 
  Ghat[27] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[27]+15.000000000000002*(alphaDrSurf[1]*fUpwind[21]+alphaDrSurf[2]*fUpwind[20])+15.0*alphaDrSurf[5]*fUpwind[12]); 
  Ghat[28] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[32]+15.0*(alphaDrSurf[5]*fUpwind[29]+alphaDrSurf[0]*fUpwind[28])+15.000000000000002*alphaDrSurf[1]*fUpwind[22]); 
  Ghat[29] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[32]+15.0*(alphaDrSurf[0]*fUpwind[29]+alphaDrSurf[5]*fUpwind[28])+15.000000000000002*alphaDrSurf[2]*fUpwind[22]); 
  Ghat[30] = 0.25*(alphaDrSurf[5]*fUpwind[35]+alphaDrSurf[2]*fUpwind[34]+alphaDrSurf[1]*fUpwind[33]+alphaDrSurf[0]*fUpwind[30]); 
  Ghat[31] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[31]+15.000000000000002*(alphaDrSurf[1]*fUpwind[26]+alphaDrSurf[2]*fUpwind[25])+15.0*alphaDrSurf[5]*fUpwind[19]); 
  Ghat[32] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[32]+15.000000000000002*(alphaDrSurf[1]*fUpwind[29]+alphaDrSurf[2]*fUpwind[28])+15.0*alphaDrSurf[5]*fUpwind[22]); 
  Ghat[33] = 0.25*(alphaDrSurf[2]*fUpwind[35]+alphaDrSurf[5]*fUpwind[34]+alphaDrSurf[0]*fUpwind[33]+alphaDrSurf[1]*fUpwind[30]); 
  Ghat[34] = 0.25*(alphaDrSurf[1]*fUpwind[35]+alphaDrSurf[0]*fUpwind[34]+alphaDrSurf[5]*fUpwind[33]+alphaDrSurf[2]*fUpwind[30]); 
  Ghat[35] = 0.25*(alphaDrSurf[0]*fUpwind[35]+alphaDrSurf[1]*fUpwind[34]+alphaDrSurf[2]*fUpwind[33]+alphaDrSurf[5]*fUpwind[30]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 1.224744871391589*Ghat[1]*rdv2; 
  out[8] += 1.224744871391589*Ghat[2]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 1.224744871391589*Ghat[4]*rdv2; 
  out[15] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[19] += 1.224744871391589*Ghat[5]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[24] += 1.224744871391589*Ghat[8]*rdv2; 
  out[25] += 1.224744871391589*Ghat[9]*rdv2; 
  out[26] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[27] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[28] += 1.224744871391589*Ghat[10]*rdv2; 
  out[29] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[30] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[31] += 1.5811388300841898*Ghat[3]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[34] += 1.224744871391589*Ghat[11]*rdv2; 
  out[35] += 1.5811388300841898*Ghat[4]*rdv2; 
  out[36] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[37] += 0.7071067811865475*Ghat[20]*rdv2; 
  out[38] += 0.7071067811865475*Ghat[21]*rdv2; 
  out[39] += 1.224744871391589*Ghat[12]*rdv2; 
  out[40] += 0.7071067811865475*Ghat[22]*rdv2; 
  out[41] += 1.224744871391589*Ghat[13]*rdv2; 
  out[42] += 1.224744871391589*Ghat[14]*rdv2; 
  out[43] += 0.7071067811865475*Ghat[23]*rdv2; 
  out[44] += 1.224744871391589*Ghat[15]*rdv2; 
  out[45] += 1.224744871391589*Ghat[16]*rdv2; 
  out[46] += 1.5811388300841895*Ghat[5]*rdv2; 
  out[47] += 1.5811388300841895*Ghat[6]*rdv2; 
  out[48] += 1.5811388300841895*Ghat[7]*rdv2; 
  out[49] += 0.7071067811865475*Ghat[24]*rdv2; 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += 1.224744871391589*Ghat[18]*rdv2; 
  out[52] += 1.5811388300841895*Ghat[8]*rdv2; 
  out[53] += 1.5811388300841895*Ghat[9]*rdv2; 
  out[54] += 1.5811388300841895*Ghat[10]*rdv2; 
  out[55] += 0.7071067811865475*Ghat[25]*rdv2; 
  out[56] += 0.7071067811865475*Ghat[26]*rdv2; 
  out[57] += 1.224744871391589*Ghat[19]*rdv2; 
  out[58] += 0.7071067811865475*Ghat[27]*rdv2; 
  out[59] += 1.224744871391589*Ghat[20]*rdv2; 
  out[60] += 1.224744871391589*Ghat[21]*rdv2; 
  out[61] += 0.7071067811865475*Ghat[28]*rdv2; 
  out[62] += 0.7071067811865475*Ghat[29]*rdv2; 
  out[63] += 1.224744871391589*Ghat[22]*rdv2; 
  out[64] += 1.5811388300841895*Ghat[11]*rdv2; 
  out[65] += 1.5811388300841895*Ghat[12]*rdv2; 
  out[66] += 0.7071067811865475*Ghat[30]*rdv2; 
  out[67] += 1.224744871391589*Ghat[23]*rdv2; 
  out[68] += 1.5811388300841898*Ghat[13]*rdv2; 
  out[69] += 1.224744871391589*Ghat[24]*rdv2; 
  out[70] += 1.5811388300841898*Ghat[14]*rdv2; 
  out[71] += 1.5811388300841898*Ghat[15]*rdv2; 
  out[72] += 1.5811388300841898*Ghat[16]*rdv2; 
  out[73] += 0.7071067811865475*Ghat[31]*rdv2; 
  out[74] += 1.224744871391589*Ghat[25]*rdv2; 
  out[75] += 1.224744871391589*Ghat[26]*rdv2; 
  out[76] += 1.224744871391589*Ghat[27]*rdv2; 
  out[77] += 0.7071067811865475*Ghat[32]*rdv2; 
  out[78] += 1.224744871391589*Ghat[28]*rdv2; 
  out[79] += 1.224744871391589*Ghat[29]*rdv2; 
  out[80] += 1.5811388300841898*Ghat[17]*rdv2; 
  out[81] += 1.5811388300841898*Ghat[18]*rdv2; 
  out[82] += 1.5811388300841898*Ghat[19]*rdv2; 
  out[83] += 1.5811388300841898*Ghat[20]*rdv2; 
  out[84] += 1.5811388300841898*Ghat[21]*rdv2; 
  out[85] += 1.5811388300841898*Ghat[22]*rdv2; 
  out[86] += 0.7071067811865475*Ghat[33]*rdv2; 
  out[87] += 0.7071067811865475*Ghat[34]*rdv2; 
  out[88] += 1.224744871391589*Ghat[30]*rdv2; 
  out[89] += 1.5811388300841895*Ghat[23]*rdv2; 
  out[90] += 1.224744871391589*Ghat[31]*rdv2; 
  out[91] += 1.224744871391589*Ghat[32]*rdv2; 
  out[92] += 1.5811388300841895*Ghat[24]*rdv2; 
  out[93] += 1.5811388300841895*Ghat[25]*rdv2; 
  out[94] += 1.5811388300841895*Ghat[26]*rdv2; 
  out[95] += 1.5811388300841895*Ghat[27]*rdv2; 
  out[96] += 1.5811388300841895*Ghat[28]*rdv2; 
  out[97] += 1.5811388300841895*Ghat[29]*rdv2; 
  out[98] += 0.7071067811865475*Ghat[35]*rdv2; 
  out[99] += 1.224744871391589*Ghat[33]*rdv2; 
  out[100] += 1.224744871391589*Ghat[34]*rdv2; 
  out[101] += 1.5811388300841895*Ghat[30]*rdv2; 
  out[102] += 1.5811388300841898*Ghat[31]*rdv2; 
  out[103] += 1.5811388300841898*Ghat[32]*rdv2; 
  out[104] += 1.224744871391589*Ghat[35]*rdv2; 
  out[105] += 1.5811388300841898*Ghat[33]*rdv2; 
  out[106] += 1.5811388300841898*Ghat[34]*rdv2; 
  out[107] += 1.5811388300841895*Ghat[35]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[1]; 
  alphaDrSurf[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]-1.0*dxv[2]*nuSum[2]; 
  alphaDrSurf[5] = (2.0*w[2]-1.0*dxv[2])*nuSum[3]-2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (-alphaDrSurf[5]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fEdge); 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fSkin); 
    fUpwindQuad[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fSkin); 
  } 
  if (-alphaDrSurf[5]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fEdge); 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fEdge); 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fSkin); 
    fUpwindQuad[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fSkin); 
    fUpwindQuad[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fEdge); 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fEdge); 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fEdge); 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fEdge); 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fEdge); 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fEdge); 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fEdge); 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fEdge); 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fSkin); 
    fUpwindQuad[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fSkin); 
    fUpwindQuad[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fSkin); 
    fUpwindQuad[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fSkin); 
    fUpwindQuad[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fSkin); 
    fUpwindQuad[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fSkin); 
    fUpwindQuad[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fSkin); 
    fUpwindQuad[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fSkin); 
    fUpwindQuad[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alphaDrSurf[5]*fUpwind[5]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alphaDrSurf[2]*fUpwind[5]+fUpwind[2]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.25*(alphaDrSurf[1]*fUpwind[5]+fUpwind[1]*alphaDrSurf[5]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.25*(alphaDrSurf[5]*fUpwind[13]+alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.25*(alphaDrSurf[5]*fUpwind[14]+alphaDrSurf[2]*fUpwind[9]+alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.25*(alphaDrSurf[0]*fUpwind[5]+fUpwind[0]*alphaDrSurf[5]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[6] = 0.25*(alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[5]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[7] = 0.25*(alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[5]*fUpwind[6]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[8] = 0.25*(alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[5]*fUpwind[9]+alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[9] = 0.25*(alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[5]*fUpwind[8]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[10] = 0.25*(alphaDrSurf[5]*fUpwind[23]+alphaDrSurf[2]*fUpwind[16]+alphaDrSurf[1]*fUpwind[15]+alphaDrSurf[0]*fUpwind[10]); 
  Ghat[11] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[24]+15.000000000000002*(alphaDrSurf[2]*fUpwind[18]+alphaDrSurf[1]*fUpwind[17])+15.0*alphaDrSurf[0]*fUpwind[11]); 
  Ghat[12] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[27]+15.000000000000002*(alphaDrSurf[2]*fUpwind[21]+alphaDrSurf[1]*fUpwind[20])+15.0*alphaDrSurf[0]*fUpwind[12]); 
  Ghat[13] = 0.25*(alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[3]*alphaDrSurf[5]); 
  Ghat[14] = 0.25*(alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[2]*fUpwind[8]+fUpwind[4]*alphaDrSurf[5]); 
  Ghat[15] = 0.25*(alphaDrSurf[2]*fUpwind[23]+alphaDrSurf[5]*fUpwind[16]+alphaDrSurf[0]*fUpwind[15]+alphaDrSurf[1]*fUpwind[10]); 
  Ghat[16] = 0.25*(alphaDrSurf[1]*fUpwind[23]+alphaDrSurf[0]*fUpwind[16]+alphaDrSurf[5]*fUpwind[15]+alphaDrSurf[2]*fUpwind[10]); 
  Ghat[17] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[24]+15.0*(alphaDrSurf[5]*fUpwind[18]+alphaDrSurf[0]*fUpwind[17])+15.000000000000002*alphaDrSurf[1]*fUpwind[11]); 
  Ghat[18] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[24]+15.0*(alphaDrSurf[0]*fUpwind[18]+alphaDrSurf[5]*fUpwind[17])+15.000000000000002*alphaDrSurf[2]*fUpwind[11]); 
  Ghat[19] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[31]+15.000000000000002*(alphaDrSurf[2]*fUpwind[26]+alphaDrSurf[1]*fUpwind[25])+15.0*alphaDrSurf[0]*fUpwind[19]); 
  Ghat[20] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[27]+15.0*(alphaDrSurf[5]*fUpwind[21]+alphaDrSurf[0]*fUpwind[20])+15.000000000000002*alphaDrSurf[1]*fUpwind[12]); 
  Ghat[21] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[27]+15.0*(alphaDrSurf[0]*fUpwind[21]+alphaDrSurf[5]*fUpwind[20])+15.000000000000002*alphaDrSurf[2]*fUpwind[12]); 
  Ghat[22] = 0.016666666666666666*(15.0*alphaDrSurf[5]*fUpwind[32]+15.000000000000002*(alphaDrSurf[2]*fUpwind[29]+alphaDrSurf[1]*fUpwind[28])+15.0*alphaDrSurf[0]*fUpwind[22]); 
  Ghat[23] = 0.25*(alphaDrSurf[0]*fUpwind[23]+alphaDrSurf[1]*fUpwind[16]+alphaDrSurf[2]*fUpwind[15]+alphaDrSurf[5]*fUpwind[10]); 
  Ghat[24] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[24]+15.000000000000002*(alphaDrSurf[1]*fUpwind[18]+alphaDrSurf[2]*fUpwind[17])+15.0*alphaDrSurf[5]*fUpwind[11]); 
  Ghat[25] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[31]+15.0*(alphaDrSurf[5]*fUpwind[26]+alphaDrSurf[0]*fUpwind[25])+15.000000000000002*alphaDrSurf[1]*fUpwind[19]); 
  Ghat[26] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[31]+15.0*(alphaDrSurf[0]*fUpwind[26]+alphaDrSurf[5]*fUpwind[25])+15.000000000000002*alphaDrSurf[2]*fUpwind[19]); 
  Ghat[27] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[27]+15.000000000000002*(alphaDrSurf[1]*fUpwind[21]+alphaDrSurf[2]*fUpwind[20])+15.0*alphaDrSurf[5]*fUpwind[12]); 
  Ghat[28] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[2]*fUpwind[32]+15.0*(alphaDrSurf[5]*fUpwind[29]+alphaDrSurf[0]*fUpwind[28])+15.000000000000002*alphaDrSurf[1]*fUpwind[22]); 
  Ghat[29] = 0.016666666666666666*(15.000000000000002*alphaDrSurf[1]*fUpwind[32]+15.0*(alphaDrSurf[0]*fUpwind[29]+alphaDrSurf[5]*fUpwind[28])+15.000000000000002*alphaDrSurf[2]*fUpwind[22]); 
  Ghat[30] = 0.25*(alphaDrSurf[5]*fUpwind[35]+alphaDrSurf[2]*fUpwind[34]+alphaDrSurf[1]*fUpwind[33]+alphaDrSurf[0]*fUpwind[30]); 
  Ghat[31] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[31]+15.000000000000002*(alphaDrSurf[1]*fUpwind[26]+alphaDrSurf[2]*fUpwind[25])+15.0*alphaDrSurf[5]*fUpwind[19]); 
  Ghat[32] = 0.016666666666666666*(15.0*alphaDrSurf[0]*fUpwind[32]+15.000000000000002*(alphaDrSurf[1]*fUpwind[29]+alphaDrSurf[2]*fUpwind[28])+15.0*alphaDrSurf[5]*fUpwind[22]); 
  Ghat[33] = 0.25*(alphaDrSurf[2]*fUpwind[35]+alphaDrSurf[5]*fUpwind[34]+alphaDrSurf[0]*fUpwind[33]+alphaDrSurf[1]*fUpwind[30]); 
  Ghat[34] = 0.25*(alphaDrSurf[1]*fUpwind[35]+alphaDrSurf[0]*fUpwind[34]+alphaDrSurf[5]*fUpwind[33]+alphaDrSurf[2]*fUpwind[30]); 
  Ghat[35] = 0.25*(alphaDrSurf[0]*fUpwind[35]+alphaDrSurf[1]*fUpwind[34]+alphaDrSurf[2]*fUpwind[33]+alphaDrSurf[5]*fUpwind[30]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[6] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[7] += 1.224744871391589*Ghat[1]*rdv2; 
  out[8] += 1.224744871391589*Ghat[2]*rdv2; 
  out[9] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[10] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[13] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[14] += 1.224744871391589*Ghat[4]*rdv2; 
  out[15] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[16] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[17] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[18] += -(0.7071067811865475*Ghat[12]*rdv2); 
  out[19] += 1.224744871391589*Ghat[5]*rdv2; 
  out[20] += -(0.7071067811865475*Ghat[13]*rdv2); 
  out[21] += 1.224744871391589*Ghat[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += -(0.7071067811865475*Ghat[14]*rdv2); 
  out[24] += 1.224744871391589*Ghat[8]*rdv2; 
  out[25] += 1.224744871391589*Ghat[9]*rdv2; 
  out[26] += -(0.7071067811865475*Ghat[15]*rdv2); 
  out[27] += -(0.7071067811865475*Ghat[16]*rdv2); 
  out[28] += 1.224744871391589*Ghat[10]*rdv2; 
  out[29] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[30] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[31] += -(1.5811388300841898*Ghat[3]*rdv2); 
  out[32] += -(0.7071067811865475*Ghat[17]*rdv2); 
  out[33] += -(0.7071067811865475*Ghat[18]*rdv2); 
  out[34] += 1.224744871391589*Ghat[11]*rdv2; 
  out[35] += -(1.5811388300841898*Ghat[4]*rdv2); 
  out[36] += -(0.7071067811865475*Ghat[19]*rdv2); 
  out[37] += -(0.7071067811865475*Ghat[20]*rdv2); 
  out[38] += -(0.7071067811865475*Ghat[21]*rdv2); 
  out[39] += 1.224744871391589*Ghat[12]*rdv2; 
  out[40] += -(0.7071067811865475*Ghat[22]*rdv2); 
  out[41] += 1.224744871391589*Ghat[13]*rdv2; 
  out[42] += 1.224744871391589*Ghat[14]*rdv2; 
  out[43] += -(0.7071067811865475*Ghat[23]*rdv2); 
  out[44] += 1.224744871391589*Ghat[15]*rdv2; 
  out[45] += 1.224744871391589*Ghat[16]*rdv2; 
  out[46] += -(1.5811388300841895*Ghat[5]*rdv2); 
  out[47] += -(1.5811388300841895*Ghat[6]*rdv2); 
  out[48] += -(1.5811388300841895*Ghat[7]*rdv2); 
  out[49] += -(0.7071067811865475*Ghat[24]*rdv2); 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += 1.224744871391589*Ghat[18]*rdv2; 
  out[52] += -(1.5811388300841895*Ghat[8]*rdv2); 
  out[53] += -(1.5811388300841895*Ghat[9]*rdv2); 
  out[54] += -(1.5811388300841895*Ghat[10]*rdv2); 
  out[55] += -(0.7071067811865475*Ghat[25]*rdv2); 
  out[56] += -(0.7071067811865475*Ghat[26]*rdv2); 
  out[57] += 1.224744871391589*Ghat[19]*rdv2; 
  out[58] += -(0.7071067811865475*Ghat[27]*rdv2); 
  out[59] += 1.224744871391589*Ghat[20]*rdv2; 
  out[60] += 1.224744871391589*Ghat[21]*rdv2; 
  out[61] += -(0.7071067811865475*Ghat[28]*rdv2); 
  out[62] += -(0.7071067811865475*Ghat[29]*rdv2); 
  out[63] += 1.224744871391589*Ghat[22]*rdv2; 
  out[64] += -(1.5811388300841895*Ghat[11]*rdv2); 
  out[65] += -(1.5811388300841895*Ghat[12]*rdv2); 
  out[66] += -(0.7071067811865475*Ghat[30]*rdv2); 
  out[67] += 1.224744871391589*Ghat[23]*rdv2; 
  out[68] += -(1.5811388300841898*Ghat[13]*rdv2); 
  out[69] += 1.224744871391589*Ghat[24]*rdv2; 
  out[70] += -(1.5811388300841898*Ghat[14]*rdv2); 
  out[71] += -(1.5811388300841898*Ghat[15]*rdv2); 
  out[72] += -(1.5811388300841898*Ghat[16]*rdv2); 
  out[73] += -(0.7071067811865475*Ghat[31]*rdv2); 
  out[74] += 1.224744871391589*Ghat[25]*rdv2; 
  out[75] += 1.224744871391589*Ghat[26]*rdv2; 
  out[76] += 1.224744871391589*Ghat[27]*rdv2; 
  out[77] += -(0.7071067811865475*Ghat[32]*rdv2); 
  out[78] += 1.224744871391589*Ghat[28]*rdv2; 
  out[79] += 1.224744871391589*Ghat[29]*rdv2; 
  out[80] += -(1.5811388300841898*Ghat[17]*rdv2); 
  out[81] += -(1.5811388300841898*Ghat[18]*rdv2); 
  out[82] += -(1.5811388300841898*Ghat[19]*rdv2); 
  out[83] += -(1.5811388300841898*Ghat[20]*rdv2); 
  out[84] += -(1.5811388300841898*Ghat[21]*rdv2); 
  out[85] += -(1.5811388300841898*Ghat[22]*rdv2); 
  out[86] += -(0.7071067811865475*Ghat[33]*rdv2); 
  out[87] += -(0.7071067811865475*Ghat[34]*rdv2); 
  out[88] += 1.224744871391589*Ghat[30]*rdv2; 
  out[89] += -(1.5811388300841895*Ghat[23]*rdv2); 
  out[90] += 1.224744871391589*Ghat[31]*rdv2; 
  out[91] += 1.224744871391589*Ghat[32]*rdv2; 
  out[92] += -(1.5811388300841895*Ghat[24]*rdv2); 
  out[93] += -(1.5811388300841895*Ghat[25]*rdv2); 
  out[94] += -(1.5811388300841895*Ghat[26]*rdv2); 
  out[95] += -(1.5811388300841895*Ghat[27]*rdv2); 
  out[96] += -(1.5811388300841895*Ghat[28]*rdv2); 
  out[97] += -(1.5811388300841895*Ghat[29]*rdv2); 
  out[98] += -(0.7071067811865475*Ghat[35]*rdv2); 
  out[99] += 1.224744871391589*Ghat[33]*rdv2; 
  out[100] += 1.224744871391589*Ghat[34]*rdv2; 
  out[101] += -(1.5811388300841895*Ghat[30]*rdv2); 
  out[102] += -(1.5811388300841898*Ghat[31]*rdv2); 
  out[103] += -(1.5811388300841898*Ghat[32]*rdv2); 
  out[104] += 1.224744871391589*Ghat[35]*rdv2; 
  out[105] += -(1.5811388300841898*Ghat[33]*rdv2); 
  out[106] += -(1.5811388300841898*Ghat[34]*rdv2); 
  out[107] += -(1.5811388300841895*Ghat[35]*rdv2); 

  } 

  return 0.;

} 
