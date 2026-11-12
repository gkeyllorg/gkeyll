#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[3] = {0.0}; 
  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};;
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -(0.5*(2.0*sumNuUx[2]+(-(2.0*w[1])-1.0*dxv[1])*nuSum[2])); 

  if (0.6324555320336768*alphaDrSurf[2]-0.9486832980505135*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865468*alphaDrSurf[0]-0.7905694150420945*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (0.6324555320336768*alphaDrSurf[2]+0.9486832980505135*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.45175395145262565*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[8] += 1.5811388300841895*Ghat[2]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]-1.0*dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -(0.5*(2.0*sumNuUx[2]+(dxv[1]-2.0*w[1])*nuSum[2])); 

  if (0.6324555320336768*alphaDrSurf[2]-0.9486832980505135*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865468*alphaDrSurf[0]-0.7905694150420945*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (0.6324555320336768*alphaDrSurf[2]+0.9486832980505135*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.45175395145262565*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[5] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[8] += -(1.5811388300841895*Ghat[2]*rdv2); 

  } 

  return 0.;

} 
