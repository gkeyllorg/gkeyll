#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[18] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[1]-1.0*nuSum[0]*dxv[1]-2.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]-1.0*dxv[1]*nuSum[1]; 

  double alphaDrSurf_r[18] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[1]+nuSum[0]*dxv[1]-2.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]; 

  double fUpwindQuad_l[18] = {0.0};
  double fUpwindQuad_r[18] = {0.0};
  double fUpwind_l[18] = {0.0};
  double fUpwind_r[18] = {0.0};
  double Ghat_l[18] = {0.0}; 
  double Ghat_r[18] = {0.0}; 

  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] < 0) { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fl); 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fl); 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fc); 
    fUpwindQuad_l[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fc); 
    fUpwindQuad_l[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fc); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] < 0) { 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_r(fc); 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_r(fc); 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_1x3v_p1_surfx2_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_1x3v_p1_surfx2_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_1x3v_p1_surfx2_eval_quad_node_2_l(fr); 
    fUpwindQuad_r[3] = hyb_1x3v_p1_surfx2_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_1x3v_p1_surfx2_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_1x3v_p1_surfx2_eval_quad_node_5_l(fr); 
    fUpwindQuad_r[6] = hyb_1x3v_p1_surfx2_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = hyb_1x3v_p1_surfx2_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = hyb_1x3v_p1_surfx2_eval_quad_node_8_l(fr); 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fl); 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fl); 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fl); 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fl); 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fl); 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fl); 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fc); 
    fUpwindQuad_l[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fc); 
    fUpwindQuad_l[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fc); 
    fUpwindQuad_l[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fc); 
    fUpwindQuad_l[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fc); 
    fUpwindQuad_l[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fc); 
    fUpwindQuad_l[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fc); 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_r(fc); 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_r(fc); 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_r(fc); 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_r(fc); 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_r(fc); 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_r(fc); 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_1x3v_p1_surfx2_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = hyb_1x3v_p1_surfx2_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = hyb_1x3v_p1_surfx2_eval_quad_node_11_l(fr); 
    fUpwindQuad_r[12] = hyb_1x3v_p1_surfx2_eval_quad_node_12_l(fr); 
    fUpwindQuad_r[13] = hyb_1x3v_p1_surfx2_eval_quad_node_13_l(fr); 
    fUpwindQuad_r[14] = hyb_1x3v_p1_surfx2_eval_quad_node_14_l(fr); 
    fUpwindQuad_r[15] = hyb_1x3v_p1_surfx2_eval_quad_node_15_l(fr); 
    fUpwindQuad_r[16] = hyb_1x3v_p1_surfx2_eval_quad_node_16_l(fr); 
    fUpwindQuad_r[17] = hyb_1x3v_p1_surfx2_eval_quad_node_17_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[4]+alphaDrSurf_l[0]*fUpwind_l[2]); 
  Ghat_l[3] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[5]+alphaDrSurf_l[0]*fUpwind_l[3]); 
  Ghat_l[4] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[4]+alphaDrSurf_l[1]*fUpwind_l[2]); 
  Ghat_l[5] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[5]+alphaDrSurf_l[1]*fUpwind_l[3]); 
  Ghat_l[6] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[9]+alphaDrSurf_l[0]*fUpwind_l[6]); 
  Ghat_l[7] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[10]+15.0*alphaDrSurf_l[0]*fUpwind_l[7]); 
  Ghat_l[8] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[12]+15.0*alphaDrSurf_l[0]*fUpwind_l[8]); 
  Ghat_l[9] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[9]+alphaDrSurf_l[1]*fUpwind_l[6]); 
  Ghat_l[10] = 0.023570226039551577*(15.0*alphaDrSurf_l[0]*fUpwind_l[10]+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[7]); 
  Ghat_l[11] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[14]+15.0*alphaDrSurf_l[0]*fUpwind_l[11]); 
  Ghat_l[12] = 0.023570226039551577*(15.0*alphaDrSurf_l[0]*fUpwind_l[12]+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[8]); 
  Ghat_l[13] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[15]+15.0*alphaDrSurf_l[0]*fUpwind_l[13]); 
  Ghat_l[14] = 0.023570226039551577*(15.0*alphaDrSurf_l[0]*fUpwind_l[14]+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[11]); 
  Ghat_l[15] = 0.023570226039551577*(15.0*alphaDrSurf_l[0]*fUpwind_l[15]+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[13]); 
  Ghat_l[16] = 0.3535533905932737*(alphaDrSurf_l[1]*fUpwind_l[17]+alphaDrSurf_l[0]*fUpwind_l[16]); 
  Ghat_l[17] = 0.3535533905932737*(alphaDrSurf_l[0]*fUpwind_l[17]+alphaDrSurf_l[1]*fUpwind_l[16]); 

  Ghat_r[0] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[4]+alphaDrSurf_r[0]*fUpwind_r[2]); 
  Ghat_r[3] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[5]+alphaDrSurf_r[0]*fUpwind_r[3]); 
  Ghat_r[4] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[4]+alphaDrSurf_r[1]*fUpwind_r[2]); 
  Ghat_r[5] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[5]+alphaDrSurf_r[1]*fUpwind_r[3]); 
  Ghat_r[6] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[9]+alphaDrSurf_r[0]*fUpwind_r[6]); 
  Ghat_r[7] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[10]+15.0*alphaDrSurf_r[0]*fUpwind_r[7]); 
  Ghat_r[8] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[12]+15.0*alphaDrSurf_r[0]*fUpwind_r[8]); 
  Ghat_r[9] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[9]+alphaDrSurf_r[1]*fUpwind_r[6]); 
  Ghat_r[10] = 0.023570226039551577*(15.0*alphaDrSurf_r[0]*fUpwind_r[10]+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[7]); 
  Ghat_r[11] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[14]+15.0*alphaDrSurf_r[0]*fUpwind_r[11]); 
  Ghat_r[12] = 0.023570226039551577*(15.0*alphaDrSurf_r[0]*fUpwind_r[12]+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[8]); 
  Ghat_r[13] = 0.023570226039551577*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[15]+15.0*alphaDrSurf_r[0]*fUpwind_r[13]); 
  Ghat_r[14] = 0.023570226039551577*(15.0*alphaDrSurf_r[0]*fUpwind_r[14]+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[11]); 
  Ghat_r[15] = 0.023570226039551577*(15.0*alphaDrSurf_r[0]*fUpwind_r[15]+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[13]); 
  Ghat_r[16] = 0.3535533905932737*(alphaDrSurf_r[1]*fUpwind_r[17]+alphaDrSurf_r[0]*fUpwind_r[16]); 
  Ghat_r[17] = 0.3535533905932737*(alphaDrSurf_r[0]*fUpwind_r[17]+alphaDrSurf_r[1]*fUpwind_r[16]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[9] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[11] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[16] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[17] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[18] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[19] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[22] += 1.5811388300841898*Ghat_r[3]*rdv2-1.5811388300841898*Ghat_l[3]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[24] += 0.7071067811865475*Ghat_r[12]*rdv2-0.7071067811865475*Ghat_l[12]*rdv2; 
  out[25] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[26] += 0.7071067811865475*Ghat_r[13]*rdv2-0.7071067811865475*Ghat_l[13]*rdv2; 
  out[27] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[28] += 1.5811388300841895*Ghat_r[4]*rdv2-1.5811388300841895*Ghat_l[4]*rdv2; 
  out[29] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[30] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[31] += 1.5811388300841895*Ghat_r[6]*rdv2-1.5811388300841895*Ghat_l[6]*rdv2; 
  out[32] += 0.7071067811865475*Ghat_r[14]*rdv2-0.7071067811865475*Ghat_l[14]*rdv2; 
  out[33] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[34] += 1.224744871391589*Ghat_r[12]*rdv2+1.224744871391589*Ghat_l[12]*rdv2; 
  out[35] += 0.7071067811865475*Ghat_r[15]*rdv2-0.7071067811865475*Ghat_l[15]*rdv2; 
  out[36] += 1.224744871391589*Ghat_r[13]*rdv2+1.224744871391589*Ghat_l[13]*rdv2; 
  out[37] += 1.5811388300841895*Ghat_r[7]*rdv2-1.5811388300841895*Ghat_l[7]*rdv2; 
  out[38] += 1.5811388300841895*Ghat_r[8]*rdv2-1.5811388300841895*Ghat_l[8]*rdv2; 
  out[39] += 0.7071067811865475*Ghat_r[16]*rdv2-0.7071067811865475*Ghat_l[16]*rdv2; 
  out[40] += 1.5811388300841898*Ghat_r[9]*rdv2-1.5811388300841898*Ghat_l[9]*rdv2; 
  out[41] += 1.224744871391589*Ghat_r[14]*rdv2+1.224744871391589*Ghat_l[14]*rdv2; 
  out[42] += 1.224744871391589*Ghat_r[15]*rdv2+1.224744871391589*Ghat_l[15]*rdv2; 
  out[43] += 1.5811388300841898*Ghat_r[10]*rdv2-1.5811388300841898*Ghat_l[10]*rdv2; 
  out[44] += 1.5811388300841898*Ghat_r[11]*rdv2-1.5811388300841898*Ghat_l[11]*rdv2; 
  out[45] += 1.5811388300841898*Ghat_r[12]*rdv2-1.5811388300841898*Ghat_l[12]*rdv2; 
  out[46] += 1.5811388300841898*Ghat_r[13]*rdv2-1.5811388300841898*Ghat_l[13]*rdv2; 
  out[47] += 0.7071067811865475*Ghat_r[17]*rdv2-0.7071067811865475*Ghat_l[17]*rdv2; 
  out[48] += 1.224744871391589*Ghat_r[16]*rdv2+1.224744871391589*Ghat_l[16]*rdv2; 
  out[49] += 1.5811388300841895*Ghat_r[14]*rdv2-1.5811388300841895*Ghat_l[14]*rdv2; 
  out[50] += 1.5811388300841895*Ghat_r[15]*rdv2-1.5811388300841895*Ghat_l[15]*rdv2; 
  out[51] += 1.224744871391589*Ghat_r[17]*rdv2+1.224744871391589*Ghat_l[17]*rdv2; 
  out[52] += 1.5811388300841895*Ghat_r[16]*rdv2-1.5811388300841895*Ghat_l[16]*rdv2; 
  out[53] += 1.5811388300841898*Ghat_r[17]*rdv2-1.5811388300841898*Ghat_l[17]*rdv2; 

  return 0.;

} 
