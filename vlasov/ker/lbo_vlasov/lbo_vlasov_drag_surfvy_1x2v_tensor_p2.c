#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[9]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuPrimMomsSum[3]; 

  double alphaDrSurf_l[9] = {0.0}; 
  alphaDrSurf_l[0] = 1.4142135623730951*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.4142135623730951*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.4142135623730951*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.4142135623730951*sumNuUy[1]; 
  alphaDrSurf_l[4] = 1.4142135623730951*nuSum[2]*w[2]-1.4142135623730951*sumNuUy[2]-0.7071067811865475*dxv[2]*nuSum[2]; 

  double alphaDrSurf_r[9] = {0.0}; 
  alphaDrSurf_r[0] = 1.4142135623730951*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.4142135623730951*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.4142135623730951*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.4142135623730951*sumNuUy[1]; 
  alphaDrSurf_r[4] = 1.4142135623730951*nuSum[2]*w[2]-1.4142135623730951*sumNuUy[2]+0.7071067811865475*dxv[2]*nuSum[2]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[9] = {0.0};
  double fUpwind_r[9] = {0.0};
  double Ghat_l[9] = {0.0}; 
  double Ghat_r[9] = {0.0}; 

  if (0.44721359549995715*alphaDrSurf_l[4]-0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if (0.44721359549995715*alphaDrSurf_r[4]-0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[0] = tensor_3x_p2_surfx3_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = tensor_3x_p2_surfx3_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = tensor_3x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749465*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749465*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[3] = tensor_3x_p2_surfx3_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = tensor_3x_p2_surfx3_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = tensor_3x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.44721359549995715*alphaDrSurf_l[4]+0.6708203932499357*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (0.44721359549995715*alphaDrSurf_r[4]+0.6708203932499357*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[6] = tensor_3x_p2_surfx3_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = tensor_3x_p2_surfx3_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = tensor_3x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = 0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[4] = 0.31943828249996997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.5*alphaDrSurf_l[4]*fUpwind_l[8]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Ghat_l[6] = 0.31943828249996997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[7] = 0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[8]+0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Ghat_l[8] = 0.31943828249996997*alphaDrSurf_l[4]*fUpwind_l[8]+0.5*alphaDrSurf_l[0]*fUpwind_l[8]+0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[4]*fUpwind_l[5]; 

  Ghat_r[0] = 0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = 0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[4] = 0.31943828249996997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.5*alphaDrSurf_r[4]*fUpwind_r[8]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Ghat_r[6] = 0.31943828249996997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[7] = 0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[8]+0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Ghat_r[8] = 0.31943828249996997*alphaDrSurf_r[4]*fUpwind_r[8]+0.5*alphaDrSurf_r[0]*fUpwind_r[8]+0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[4]*fUpwind_r[5]; 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[9] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[15] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[16] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[19] += 1.5811388300841895*Ghat_r[3]*rdv2-1.5811388300841895*Ghat_l[3]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[21] += 1.5811388300841895*Ghat_r[4]*rdv2-1.5811388300841895*Ghat_l[4]*rdv2; 
  out[22] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[23] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[24] += 1.5811388300841898*Ghat_r[6]*rdv2-1.5811388300841898*Ghat_l[6]*rdv2; 
  out[25] += 1.5811388300841898*Ghat_r[7]*rdv2-1.5811388300841898*Ghat_l[7]*rdv2; 
  out[26] += 1.5811388300841895*Ghat_r[8]*rdv2-1.5811388300841895*Ghat_l[8]*rdv2; 

  return 0.;

} 
