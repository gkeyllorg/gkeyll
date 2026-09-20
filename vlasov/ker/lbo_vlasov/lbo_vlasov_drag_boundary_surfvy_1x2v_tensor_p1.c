#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuPrimMomsSum[2]; 

  double alphaDrSurf[6] = {0.0}; 
  double fUpwindQuad[6] = {0.0};
  double fUpwind[6] = {0.0};;
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.8284271247461907*w[2]+1.4142135623730951*dxv[2])-2.8284271247461907*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.8284271247461907*w[2]+1.4142135623730951*dxv[2])-2.8284271247461907*sumNuUy[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 
  Ghat[4] = 0.03333333333333333*(15.000000000000002*alphaDrSurf[1]*fUpwind[5]+15.0*alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.03333333333333333*(15.0*alphaDrSurf[0]*fUpwind[5]+15.000000000000002*alphaDrSurf[1]*fUpwind[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[9] += 1.224744871391589*Ghat[3]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[13] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += 1.5811388300841895*Ghat[3]*rdv2; 
  out[16] += 1.5811388300841895*Ghat[4]*rdv2; 
  out[17] += 1.5811388300841898*Ghat[5]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.8284271247461907*w[2]-1.4142135623730951*dxv[2])-2.8284271247461907*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.8284271247461907*w[2]-1.4142135623730951*dxv[2])-2.8284271247461907*sumNuUy[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x2v_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = hyb_1x2v_p1_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = hyb_1x2v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_1x2v_p1_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = hyb_1x2v_p1_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = hyb_1x2v_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 
  Ghat[4] = 0.03333333333333333*(15.000000000000002*alphaDrSurf[1]*fUpwind[5]+15.0*alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.03333333333333333*(15.0*alphaDrSurf[0]*fUpwind[5]+15.000000000000002*alphaDrSurf[1]*fUpwind[4]); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[8] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[9] += 1.224744871391589*Ghat[3]*rdv2; 
  out[10] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[13] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += -(1.5811388300841895*Ghat[3]*rdv2); 
  out[16] += -(1.5811388300841895*Ghat[4]*rdv2); 
  out[17] += -(1.5811388300841898*Ghat[5]*rdv2); 

  } 

  return 0.;

} 
