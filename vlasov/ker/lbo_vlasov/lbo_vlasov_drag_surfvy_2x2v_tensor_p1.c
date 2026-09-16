#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[12]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuPrimMomsSum[4]; 

  double alphaDrSurf_l[12] = {0.0}; 
  alphaDrSurf_l[0] = 1.4142135623730951*nuSum[0]*w[3]-0.7071067811865475*nuSum[0]*dxv[3]-1.4142135623730951*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.4142135623730951*nuSum[1]*w[3]-0.7071067811865475*nuSum[1]*dxv[3]-1.4142135623730951*sumNuUy[1]; 
  alphaDrSurf_l[2] = 1.4142135623730951*nuSum[2]*w[3]-0.7071067811865475*nuSum[2]*dxv[3]-1.4142135623730951*sumNuUy[2]; 
  alphaDrSurf_l[4] = 1.4142135623730951*nuSum[3]*w[3]-1.4142135623730951*sumNuUy[3]-0.7071067811865475*dxv[3]*nuSum[3]; 

  double alphaDrSurf_r[12] = {0.0}; 
  alphaDrSurf_r[0] = 1.4142135623730951*nuSum[0]*w[3]+0.7071067811865475*nuSum[0]*dxv[3]-1.4142135623730951*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.4142135623730951*nuSum[1]*w[3]+0.7071067811865475*nuSum[1]*dxv[3]-1.4142135623730951*sumNuUy[1]; 
  alphaDrSurf_r[2] = 1.4142135623730951*nuSum[2]*w[3]+0.7071067811865475*nuSum[2]*dxv[3]-1.4142135623730951*sumNuUy[2]; 
  alphaDrSurf_r[4] = 1.4142135623730951*nuSum[3]*w[3]-1.4142135623730951*sumNuUy[3]+0.7071067811865475*dxv[3]*nuSum[3]; 

  double fUpwindQuad_l[12] = {0.0};
  double fUpwindQuad_r[12] = {0.0};
  double fUpwind_l[12] = {0.0};
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fc); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fr); 
  } 
  if (-alphaDrSurf_l[4]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fc); 
  } 
  if (-alphaDrSurf_r[4]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fr); 
  } 
  if (-alphaDrSurf_l[4]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fc); 
  } 
  if (-alphaDrSurf_r[4]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fr); 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fc); 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[4]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[4]+fUpwind_l[2]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[4]+fUpwind_l[1]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.3535533905932737*(alphaDrSurf_l[4]*fUpwind_l[8]+alphaDrSurf_l[2]*fUpwind_l[6]+alphaDrSurf_l[1]*fUpwind_l[5]+alphaDrSurf_l[0]*fUpwind_l[3]); 
  Ghat_l[4] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[4]+fUpwind_l[0]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[5] = 0.3535533905932737*(alphaDrSurf_l[2]*fUpwind_l[8]+alphaDrSurf_l[4]*fUpwind_l[6]+alphaDrSurf_l[0]*fUpwind_l[5]+alphaDrSurf_l[1]*fUpwind_l[3]); 
  Ghat_l[6] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[8]+alphaDrSurf_l[0]*fUpwind_l[6]+alphaDrSurf_l[4]*fUpwind_l[5]+alphaDrSurf_l[2]*fUpwind_l[3]); 
  Ghat_l[7] = 0.023570226039551577*(15.0*alphaDrSurf_l[4]*fUpwind_l[11]+15.000000000000002*(alphaDrSurf_l[2]*fUpwind_l[10]+alphaDrSurf_l[1]*fUpwind_l[9])+15.0*alphaDrSurf_l[0]*fUpwind_l[7]); 
  Ghat_l[8] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[8]+alphaDrSurf_l[1]*fUpwind_l[6]+alphaDrSurf_l[2]*fUpwind_l[5]+fUpwind_l[3]*alphaDrSurf_l[4]); 
  Ghat_l[9] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[11]+15.0*(alphaDrSurf_l[4]*fUpwind_l[10]+alphaDrSurf_l[0]*fUpwind_l[9])+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[7]); 
  Ghat_l[10] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[11]+15.0*(alphaDrSurf_l[0]*fUpwind_l[10]+alphaDrSurf_l[4]*fUpwind_l[9])+15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[7]); 
  Ghat_l[11] = 0.023570226039551577*(15.0*alphaDrSurf_l[0]*fUpwind_l[11]+15.000000000000002*(alphaDrSurf_l[1]*fUpwind_l[10]+alphaDrSurf_l[2]*fUpwind_l[9])+15.0*alphaDrSurf_l[4]*fUpwind_l[7]); 

  Ghat_r[0] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[4]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[4]+fUpwind_r[2]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[4]+fUpwind_r[1]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.3535533905932737*(alphaDrSurf_r[4]*fUpwind_r[8]+alphaDrSurf_r[2]*fUpwind_r[6]+alphaDrSurf_r[1]*fUpwind_r[5]+alphaDrSurf_r[0]*fUpwind_r[3]); 
  Ghat_r[4] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[4]+fUpwind_r[0]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[5] = 0.3535533905932737*(alphaDrSurf_r[2]*fUpwind_r[8]+alphaDrSurf_r[4]*fUpwind_r[6]+alphaDrSurf_r[0]*fUpwind_r[5]+alphaDrSurf_r[1]*fUpwind_r[3]); 
  Ghat_r[6] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[8]+alphaDrSurf_r[0]*fUpwind_r[6]+alphaDrSurf_r[4]*fUpwind_r[5]+alphaDrSurf_r[2]*fUpwind_r[3]); 
  Ghat_r[7] = 0.023570226039551577*(15.0*alphaDrSurf_r[4]*fUpwind_r[11]+15.000000000000002*(alphaDrSurf_r[2]*fUpwind_r[10]+alphaDrSurf_r[1]*fUpwind_r[9])+15.0*alphaDrSurf_r[0]*fUpwind_r[7]); 
  Ghat_r[8] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[8]+alphaDrSurf_r[1]*fUpwind_r[6]+alphaDrSurf_r[2]*fUpwind_r[5]+fUpwind_r[3]*alphaDrSurf_r[4]); 
  Ghat_r[9] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[11]+15.0*(alphaDrSurf_r[4]*fUpwind_r[10]+alphaDrSurf_r[0]*fUpwind_r[9])+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[7]); 
  Ghat_r[10] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[11]+15.0*(alphaDrSurf_r[0]*fUpwind_r[10]+alphaDrSurf_r[4]*fUpwind_r[9])+15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[7]); 
  Ghat_r[11] = 0.023570226039551577*(15.0*alphaDrSurf_r[0]*fUpwind_r[11]+15.000000000000002*(alphaDrSurf_r[1]*fUpwind_r[10]+alphaDrSurf_r[2]*fUpwind_r[9])+15.0*alphaDrSurf_r[4]*fUpwind_r[7]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[12] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[16] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[19] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[20] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[21] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[22] += 1.5811388300841898*Ghat_r[3]*rdv2-1.5811388300841898*Ghat_l[3]*rdv2; 
  out[23] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[24] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[25] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[26] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[27] += 1.5811388300841895*Ghat_r[4]*rdv2-1.5811388300841895*Ghat_l[4]*rdv2; 
  out[28] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[29] += 1.5811388300841895*Ghat_r[6]*rdv2-1.5811388300841895*Ghat_l[6]*rdv2; 
  out[30] += 1.5811388300841895*Ghat_r[7]*rdv2-1.5811388300841895*Ghat_l[7]*rdv2; 
  out[31] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[32] += 1.5811388300841898*Ghat_r[8]*rdv2-1.5811388300841898*Ghat_l[8]*rdv2; 
  out[33] += 1.5811388300841898*Ghat_r[9]*rdv2-1.5811388300841898*Ghat_l[9]*rdv2; 
  out[34] += 1.5811388300841898*Ghat_r[10]*rdv2-1.5811388300841898*Ghat_l[10]*rdv2; 
  out[35] += 1.5811388300841895*Ghat_r[11]*rdv2-1.5811388300841895*Ghat_l[11]*rdv2; 

  return 0.;

} 
