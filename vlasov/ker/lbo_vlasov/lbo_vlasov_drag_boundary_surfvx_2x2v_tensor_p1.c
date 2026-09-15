#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[12]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[12] = {0.0}; 
  double fUpwindQuad[12] = {0.0};
  double fUpwind[12] = {0.0};;
  double Ghat[12] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]+dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -(0.7071067811865475*(2.0*sumNuUx[3]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[3])); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (-alphaDrSurf[4]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (-alphaDrSurf[4]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[8]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[8]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[7] = 0.023570226039551577*(15.0*alphaDrSurf[4]*fUpwind[11]+15.000000000000002*(alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[1]*fUpwind[9])+15.0*alphaDrSurf[0]*fUpwind[7]); 
  Ghat[8] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 
  Ghat[9] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[2]*fUpwind[11]+15.0*(alphaDrSurf[4]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9])+15.000000000000002*alphaDrSurf[1]*fUpwind[7]); 
  Ghat[10] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[11]+15.0*(alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[4]*fUpwind[9])+15.000000000000002*alphaDrSurf[2]*fUpwind[7]); 
  Ghat[11] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[11]+15.000000000000002*(alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9])+15.0*alphaDrSurf[4]*fUpwind[7]); 

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
  out[11] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[15] += 1.224744871391589*Ghat[5]*rdv2; 
  out[16] += 1.224744871391589*Ghat[6]*rdv2; 
  out[17] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[18] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[19] += 1.5811388300841898*Ghat[3]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[21] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += 1.5811388300841895*Ghat[4]*rdv2; 
  out[25] += 1.5811388300841895*Ghat[5]*rdv2; 
  out[26] += 1.5811388300841895*Ghat[6]*rdv2; 
  out[27] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[28] += 1.224744871391589*Ghat[9]*rdv2; 
  out[29] += 1.224744871391589*Ghat[10]*rdv2; 
  out[30] += 1.5811388300841895*Ghat[7]*rdv2; 
  out[31] += 1.5811388300841898*Ghat[8]*rdv2; 
  out[32] += 1.224744871391589*Ghat[11]*rdv2; 
  out[33] += 1.5811388300841898*Ghat[9]*rdv2; 
  out[34] += 1.5811388300841898*Ghat[10]*rdv2; 
  out[35] += 1.5811388300841895*Ghat[11]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]-1.0*dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -(0.7071067811865475*(2.0*sumNuUx[3]+(dxv[2]-2.0*w[2])*nuSum[3])); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (-alphaDrSurf[4]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (-alphaDrSurf[4]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[8]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[8]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[8]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[7] = 0.023570226039551577*(15.0*alphaDrSurf[4]*fUpwind[11]+15.000000000000002*(alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[1]*fUpwind[9])+15.0*alphaDrSurf[0]*fUpwind[7]); 
  Ghat[8] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[8]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 
  Ghat[9] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[2]*fUpwind[11]+15.0*(alphaDrSurf[4]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9])+15.000000000000002*alphaDrSurf[1]*fUpwind[7]); 
  Ghat[10] = 0.023570226039551577*(15.000000000000002*alphaDrSurf[1]*fUpwind[11]+15.0*(alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[4]*fUpwind[9])+15.000000000000002*alphaDrSurf[2]*fUpwind[7]); 
  Ghat[11] = 0.023570226039551577*(15.0*alphaDrSurf[0]*fUpwind[11]+15.000000000000002*(alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9])+15.0*alphaDrSurf[4]*fUpwind[7]); 

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
  out[11] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[12] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[15] += 1.224744871391589*Ghat[5]*rdv2; 
  out[16] += 1.224744871391589*Ghat[6]*rdv2; 
  out[17] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[18] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[19] += -(1.5811388300841898*Ghat[3]*rdv2); 
  out[20] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[21] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[22] += 1.224744871391589*Ghat[7]*rdv2; 
  out[23] += 1.224744871391589*Ghat[8]*rdv2; 
  out[24] += -(1.5811388300841895*Ghat[4]*rdv2); 
  out[25] += -(1.5811388300841895*Ghat[5]*rdv2); 
  out[26] += -(1.5811388300841895*Ghat[6]*rdv2); 
  out[27] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[28] += 1.224744871391589*Ghat[9]*rdv2; 
  out[29] += 1.224744871391589*Ghat[10]*rdv2; 
  out[30] += -(1.5811388300841895*Ghat[7]*rdv2); 
  out[31] += -(1.5811388300841898*Ghat[8]*rdv2); 
  out[32] += 1.224744871391589*Ghat[11]*rdv2; 
  out[33] += -(1.5811388300841898*Ghat[9]*rdv2); 
  out[34] += -(1.5811388300841898*Ghat[10]*rdv2); 
  out[35] += -(1.5811388300841895*Ghat[11]*rdv2); 

  } 

  return 0.;

} 
