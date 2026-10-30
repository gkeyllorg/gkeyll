#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -(0.5*(2.0*sumNuUx[2]+(-(2.0*w[1])-1.0*dxv[1])*nuSum[2])); 
  alphaDrSurf[3] = -(0.5*(2.0*sumNuUx[3]+(-(2.0*w[1])-1.0*dxv[1])*nuSum[3])); 

  if (-(0.5701294036773678*alphaDrSurf[3])+0.9681844646844027*alphaDrSurf[2]-1.0546722811938842*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.770272555658882*alphaDrSurf[3]-0.516430513231777*alphaDrSurf[2]-0.41639003950091236*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (-(0.770272555658882*alphaDrSurf[3])-0.516430513231777*alphaDrSurf[2]+0.41639003950091236*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5701294036773678*alphaDrSurf[3]+0.9681844646844027*alphaDrSurf[2]+1.0546722811938842*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[3]*fUpwind[3]+0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6210590034081186*alphaDrSurf[2]*fUpwind[3]+0.6210590034081186*fUpwind[2]*alphaDrSurf[3]+0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.421637021355784*alphaDrSurf[3]*fUpwind[3]+0.6210590034081186*alphaDrSurf[1]*fUpwind[3]+0.6210590034081186*fUpwind[1]*alphaDrSurf[3]+0.45175395145262565*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[3] = 0.421637021355784*alphaDrSurf[2]*fUpwind[3]+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]+0.421637021355784*fUpwind[2]*alphaDrSurf[3]+0.7071067811865475*fUpwind[0]*alphaDrSurf[3]+0.6210590034081186*alphaDrSurf[1]*fUpwind[2]+0.6210590034081186*fUpwind[1]*alphaDrSurf[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[9] += 1.8708286933869707*Ghat[0]*rdv2; 
  out[10] += 1.2247448713915887*Ghat[3]*rdv2; 
  out[11] += 1.8708286933869707*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]-1.0*dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -(0.5*(2.0*sumNuUx[2]+(dxv[1]-2.0*w[1])*nuSum[2])); 
  alphaDrSurf[3] = -(0.5*(2.0*sumNuUx[3]+(dxv[1]-2.0*w[1])*nuSum[3])); 

  if (-(0.5701294036773678*alphaDrSurf[3])+0.9681844646844027*alphaDrSurf[2]-1.0546722811938842*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.770272555658882*alphaDrSurf[3]-0.516430513231777*alphaDrSurf[2]-0.41639003950091236*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (-(0.770272555658882*alphaDrSurf[3])-0.516430513231777*alphaDrSurf[2]+0.41639003950091236*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5701294036773678*alphaDrSurf[3]+0.9681844646844027*alphaDrSurf[2]+1.0546722811938842*alphaDrSurf[1]+0.7071067811865468*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaDrSurf[3]*fUpwind[3]+0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6210590034081186*alphaDrSurf[2]*fUpwind[3]+0.6210590034081186*fUpwind[2]*alphaDrSurf[3]+0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.421637021355784*alphaDrSurf[3]*fUpwind[3]+0.6210590034081186*alphaDrSurf[1]*fUpwind[3]+0.6210590034081186*fUpwind[1]*alphaDrSurf[3]+0.45175395145262565*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[3] = 0.421637021355784*alphaDrSurf[2]*fUpwind[3]+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]+0.421637021355784*fUpwind[2]*alphaDrSurf[3]+0.7071067811865475*fUpwind[0]*alphaDrSurf[3]+0.6210590034081186*alphaDrSurf[1]*fUpwind[2]+0.6210590034081186*fUpwind[1]*alphaDrSurf[2]; 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[5] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[8] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[9] += 1.8708286933869707*Ghat[0]*rdv2; 
  out[10] += 1.2247448713915887*Ghat[3]*rdv2; 
  out[11] += 1.8708286933869707*Ghat[1]*rdv2; 

  } 

  return 0.;

} 
