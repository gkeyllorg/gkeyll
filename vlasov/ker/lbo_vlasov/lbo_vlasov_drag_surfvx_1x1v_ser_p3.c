#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]: cell-center coordinates. 
  // dxv[2]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*w[1]-0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]-0.5*dxv[1]*nuSum[1]; 
  alphaDrSurf_l[2] = -(1.0*sumNuUx[2])+w[1]*nuSum[2]-0.5*dxv[1]*nuSum[2]; 
  alphaDrSurf_l[3] = -(1.0*sumNuUx[3])+w[1]*nuSum[3]-0.5*dxv[1]*nuSum[3]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*w[1]+0.5*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = nuSum[1]*w[1]-1.0*sumNuUx[1]+0.5*dxv[1]*nuSum[1]; 
  alphaDrSurf_r[2] = -(1.0*sumNuUx[2])+w[1]*nuSum[2]+0.5*dxv[1]*nuSum[2]; 
  alphaDrSurf_r[3] = -(1.0*sumNuUx[3])+w[1]*nuSum[3]+0.5*dxv[1]*nuSum[3]; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (-(0.5701294036773678*alphaDrSurf_l[3])+0.9681844646844027*alphaDrSurf_l[2]-1.0546722811938842*alphaDrSurf_l[1]+0.7071067811865468*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (-(0.5701294036773678*alphaDrSurf_r[3])+0.9681844646844027*alphaDrSurf_r[2]-1.0546722811938842*alphaDrSurf_r[1]+0.7071067811865468*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.770272555658882*alphaDrSurf_l[3]-0.516430513231777*alphaDrSurf_l[2]-0.41639003950091236*alphaDrSurf_l[1]+0.7071067811865468*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.770272555658882*alphaDrSurf_r[3]-0.516430513231777*alphaDrSurf_r[2]-0.41639003950091236*alphaDrSurf_r[1]+0.7071067811865468*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (-(0.770272555658882*alphaDrSurf_l[3])-0.516430513231777*alphaDrSurf_l[2]+0.41639003950091236*alphaDrSurf_l[1]+0.7071067811865468*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fc); 
  } 
  if (-(0.770272555658882*alphaDrSurf_r[3])-0.516430513231777*alphaDrSurf_r[2]+0.41639003950091236*alphaDrSurf_r[1]+0.7071067811865468*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (0.5701294036773678*alphaDrSurf_l[3]+0.9681844646844027*alphaDrSurf_l[2]+1.0546722811938842*alphaDrSurf_l[1]+0.7071067811865468*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fc); 
  } 
  if (0.5701294036773678*alphaDrSurf_r[3]+0.9681844646844027*alphaDrSurf_r[2]+1.0546722811938842*alphaDrSurf_r[1]+0.7071067811865468*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaDrSurf_l[3]*fUpwind_l[3]+0.7071067811865475*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6210590034081186*alphaDrSurf_l[2]*fUpwind_l[3]+0.6210590034081186*fUpwind_l[2]*alphaDrSurf_l[3]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaDrSurf_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.421637021355784*alphaDrSurf_l[3]*fUpwind_l[3]+0.6210590034081186*alphaDrSurf_l[1]*fUpwind_l[3]+0.6210590034081186*fUpwind_l[1]*alphaDrSurf_l[3]+0.45175395145262565*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[2]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[3] = 0.421637021355784*alphaDrSurf_l[2]*fUpwind_l[3]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[3]+0.421637021355784*fUpwind_l[2]*alphaDrSurf_l[3]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[3]+0.6210590034081186*alphaDrSurf_l[1]*fUpwind_l[2]+0.6210590034081186*fUpwind_l[1]*alphaDrSurf_l[2]; 

  Ghat_r[0] = 0.7071067811865475*alphaDrSurf_r[3]*fUpwind_r[3]+0.7071067811865475*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6210590034081186*alphaDrSurf_r[2]*fUpwind_r[3]+0.6210590034081186*fUpwind_r[2]*alphaDrSurf_r[3]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaDrSurf_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.421637021355784*alphaDrSurf_r[3]*fUpwind_r[3]+0.6210590034081186*alphaDrSurf_r[1]*fUpwind_r[3]+0.6210590034081186*fUpwind_r[1]*alphaDrSurf_r[3]+0.45175395145262565*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[2]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[3] = 0.421637021355784*alphaDrSurf_r[2]*fUpwind_r[3]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[3]+0.421637021355784*fUpwind_r[2]*alphaDrSurf_r[3]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[3]+0.6210590034081186*alphaDrSurf_r[1]*fUpwind_r[2]+0.6210590034081186*fUpwind_r[1]*alphaDrSurf_r[2]; 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[5] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[9] += 1.8708286933869707*Ghat_r[0]*rdv2+1.8708286933869707*Ghat_l[0]*rdv2; 
  out[10] += 1.2247448713915887*Ghat_r[3]*rdv2+1.2247448713915887*Ghat_l[3]*rdv2; 
  out[11] += 1.8708286933869707*Ghat_r[1]*rdv2+1.8708286933869707*Ghat_l[1]*rdv2; 

  return 0.;

} 
