#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_tensor_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_tensor_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[27]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuPrimMomsSum[9]; 

  double alphaDrSurf[27] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[27] = {0.0};;
  double Ghat[27] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -(0.7071067811865475*(2.0*sumNuUy[4]+(-(2.0*w[3])-1.0*dxv[3])*nuSum[4])); 
  alphaDrSurf[8] = -(0.7071067811865475*(2.0*sumNuUy[5]+(-(2.0*w[3])-1.0*dxv[3])*nuSum[5])); 
  alphaDrSurf[11] = -(0.7071067811865475*(2.0*sumNuUy[6]+(-(2.0*w[3])-1.0*dxv[3])*nuSum[6])); 
  alphaDrSurf[12] = -(0.7071067811865475*(2.0*sumNuUy[7]+(-(2.0*w[3])-1.0*dxv[3])*nuSum[7])); 
  alphaDrSurf[20] = -(0.7071067811865475*(2.0*sumNuUy[8]+(-(2.0*w[3])-1.0*dxv[3])*nuSum[8])); 

  if (0.28284271247461906*alphaDrSurf[20]-0.42426406871192823*alphaDrSurf[12]-0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])+0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]-0.42426406871192823*alphaDrSurf[12]+0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])+0.5303300858899102*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.44194173824159233*alphaDrSurf[20]-0.3952847075210471*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])-0.5303300858899102*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]+0.42426406871192823*alphaDrSurf[12]-0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])-0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]+0.42426406871192823*alphaDrSurf[12]+0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_eval_quad_node_26_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.31622776601683794*alphaDrSurf[12]*fUpwind[20]+0.31622776601683794*fUpwind[12]*alphaDrSurf[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[4]*fUpwind[11]+0.31622776601683794*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.31622776601683794*alphaDrSurf[11]*fUpwind[20]+0.31622776601683794*fUpwind[11]*alphaDrSurf[20]+0.31622776601683794*alphaDrSurf[4]*fUpwind[12]+0.31622776601683794*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[20]+0.28284271247461906*fUpwind[4]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[11]*fUpwind[12]+0.31622776601683794*alphaDrSurf[2]*fUpwind[12]+0.28284271247461906*fUpwind[11]*alphaDrSurf[12]+0.31622776601683794*fUpwind[2]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[1]*fUpwind[11]+0.31622776601683794*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.31622776601683794*alphaDrSurf[12]*fUpwind[23]+0.3162277660168379*fUpwind[18]*alphaDrSurf[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.31622776601683794*alphaDrSurf[1]*fUpwind[13]+0.31622776601683794*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = 0.31622776601683794*alphaDrSurf[11]*fUpwind[23]+0.3162277660168379*fUpwind[17]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.31622776601683794*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.31622776601683794*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[20]+0.3535533905932737*fUpwind[8]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.22587697572631277*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.22587697572631277*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[7]*fUpwind[20]+0.3535533905932737*fUpwind[7]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.22587697572631277*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[12]*fUpwind[25]+0.3535533905932737*alphaDrSurf[11]*fUpwind[24]+0.3535533905932737*alphaDrSurf[8]*fUpwind[22]+0.3535533905932737*alphaDrSurf[7]*fUpwind[21]+0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  Ghat[10] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[23]+0.28284271247461906*fUpwind[10]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.31622776601683794*alphaDrSurf[4]*fUpwind[14]+0.31622776601683794*alphaDrSurf[4]*fUpwind[13]+0.31622776601683794*fUpwind[6]*alphaDrSurf[12]+0.31622776601683794*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.20203050891044216*alphaDrSurf[11]*fUpwind[20]+0.31622776601683794*alphaDrSurf[2]*fUpwind[20]+0.20203050891044216*fUpwind[11]*alphaDrSurf[20]+0.31622776601683794*fUpwind[2]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[12]+0.28284271247461906*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.22587697572631277*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.22587697572631277*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.31622776601683794*alphaDrSurf[1]*fUpwind[4]+0.31622776601683794*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[12] = 0.20203050891044216*alphaDrSurf[12]*fUpwind[20]+0.31622776601683794*alphaDrSurf[1]*fUpwind[20]+0.20203050891044216*fUpwind[12]*alphaDrSurf[20]+0.31622776601683794*fUpwind[1]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.22587697572631277*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.28284271247461906*alphaDrSurf[4]*fUpwind[11]+0.28284271247461906*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.31622776601683794*alphaDrSurf[2]*fUpwind[4]+0.31622776601683794*fUpwind[2]*alphaDrSurf[4]; 
  Ghat[13] = 0.22587697572631282*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[8]*fUpwind[23]+0.3535533905932737*fUpwind[14]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.22587697572631277*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.22587697572631277*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.31622776601683794*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.31622776601683794*alphaDrSurf[1]*fUpwind[5]; 
  Ghat[14] = 0.22587697572631282*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[7]*fUpwind[23]+0.3535533905932737*fUpwind[13]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.22587697572631277*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.31622776601683794*alphaDrSurf[2]*fUpwind[6]; 
  Ghat[15] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[26]+0.31622776601683794*alphaDrSurf[20]*fUpwind[25]+0.3535533905932737*alphaDrSurf[8]*fUpwind[25]+0.31622776601683794*alphaDrSurf[4]*fUpwind[24]+0.3535533905932737*alphaDrSurf[12]*fUpwind[22]+0.31622776601683794*alphaDrSurf[1]*fUpwind[21]+0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  Ghat[16] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[26]+0.31622776601683794*alphaDrSurf[4]*fUpwind[25]+0.31622776601683794*alphaDrSurf[20]*fUpwind[24]+0.3535533905932737*alphaDrSurf[7]*fUpwind[24]+0.31622776601683794*alphaDrSurf[2]*fUpwind[22]+0.3535533905932737*alphaDrSurf[11]*fUpwind[21]+0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[17] = 0.2020305089104422*alphaDrSurf[11]*fUpwind[23]+0.3162277660168379*alphaDrSurf[2]*fUpwind[23]+0.20203050891044216*fUpwind[17]*alphaDrSurf[20]+0.3162277660168379*fUpwind[6]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.22587697572631277*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.22587697572631277*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  Ghat[18] = 0.2020305089104422*alphaDrSurf[12]*fUpwind[23]+0.3162277660168379*alphaDrSurf[1]*fUpwind[23]+0.20203050891044216*fUpwind[18]*alphaDrSurf[20]+0.3162277660168379*fUpwind[5]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.28284271247461906*alphaDrSurf[4]*fUpwind[17]+0.22587697572631277*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  Ghat[19] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[26]+0.282842712474619*alphaDrSurf[11]*fUpwind[25]+0.3162277660168379*alphaDrSurf[2]*fUpwind[25]+0.282842712474619*alphaDrSurf[12]*fUpwind[24]+0.3162277660168379*alphaDrSurf[1]*fUpwind[24]+0.3162277660168379*alphaDrSurf[4]*fUpwind[22]+0.3162277660168379*alphaDrSurf[4]*fUpwind[21]+0.28284271247461906*fUpwind[19]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 
  Ghat[20] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[20]+0.22587697572631277*alphaDrSurf[7]*fUpwind[20]+0.3535533905932737*alphaDrSurf[0]*fUpwind[20]+0.22587697572631277*fUpwind[8]*alphaDrSurf[20]+0.22587697572631277*fUpwind[7]*alphaDrSurf[20]+0.3535533905932737*fUpwind[0]*alphaDrSurf[20]+0.20203050891044216*alphaDrSurf[12]*fUpwind[12]+0.31622776601683794*alphaDrSurf[1]*fUpwind[12]+0.31622776601683794*fUpwind[1]*alphaDrSurf[12]+0.20203050891044216*alphaDrSurf[11]*fUpwind[11]+0.31622776601683794*alphaDrSurf[2]*fUpwind[11]+0.31622776601683794*fUpwind[2]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[7]*fUpwind[8]+0.3535533905932737*fUpwind[7]*alphaDrSurf[8]+0.28284271247461906*alphaDrSurf[4]*fUpwind[4]; 
  Ghat[21] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[8]*fUpwind[26]+0.31622776601683794*alphaDrSurf[12]*fUpwind[25]+0.22587697572631282*alphaDrSurf[11]*fUpwind[24]+0.3535533905932737*alphaDrSurf[2]*fUpwind[24]+0.3535533905932737*alphaDrSurf[20]*fUpwind[22]+0.22587697572631277*alphaDrSurf[7]*fUpwind[21]+0.3535533905932737*alphaDrSurf[0]*fUpwind[21]+0.3162277660168379*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[11]*fUpwind[16]+0.31622776601683794*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[7]*fUpwind[9]; 
  Ghat[22] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[7]*fUpwind[26]+0.22587697572631282*alphaDrSurf[12]*fUpwind[25]+0.3535533905932737*alphaDrSurf[1]*fUpwind[25]+0.31622776601683794*alphaDrSurf[11]*fUpwind[24]+0.22587697572631277*alphaDrSurf[8]*fUpwind[22]+0.3535533905932737*alphaDrSurf[0]*fUpwind[22]+0.3535533905932737*alphaDrSurf[20]*fUpwind[21]+0.3162277660168379*alphaDrSurf[4]*fUpwind[19]+0.31622776601683794*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[12]*fUpwind[15]+0.3535533905932737*alphaDrSurf[8]*fUpwind[9]; 
  Ghat[23] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[23]+0.22587697572631277*alphaDrSurf[8]*fUpwind[23]+0.22587697572631277*alphaDrSurf[7]*fUpwind[23]+0.3535533905932737*alphaDrSurf[0]*fUpwind[23]+0.22587697572631282*fUpwind[14]*alphaDrSurf[20]+0.22587697572631282*fUpwind[13]*alphaDrSurf[20]+0.3535533905932737*fUpwind[3]*alphaDrSurf[20]+0.2020305089104422*alphaDrSurf[12]*fUpwind[18]+0.3162277660168379*alphaDrSurf[1]*fUpwind[18]+0.2020305089104422*alphaDrSurf[11]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[17]+0.3535533905932737*alphaDrSurf[7]*fUpwind[14]+0.3535533905932737*alphaDrSurf[8]*fUpwind[13]+0.31622776601683794*fUpwind[5]*alphaDrSurf[12]+0.31622776601683794*fUpwind[6]*alphaDrSurf[11]+0.28284271247461906*alphaDrSurf[4]*fUpwind[10]; 
  Ghat[24] = 0.2020305089104422*alphaDrSurf[11]*fUpwind[26]+0.3162277660168379*alphaDrSurf[2]*fUpwind[26]+0.28284271247461906*alphaDrSurf[4]*fUpwind[25]+0.20203050891044216*alphaDrSurf[20]*fUpwind[24]+0.3162277660168379*alphaDrSurf[8]*fUpwind[24]+0.22587697572631277*alphaDrSurf[7]*fUpwind[24]+0.3535533905932737*alphaDrSurf[0]*fUpwind[24]+0.31622776601683794*alphaDrSurf[11]*fUpwind[22]+0.22587697572631282*alphaDrSurf[11]*fUpwind[21]+0.3535533905932737*alphaDrSurf[2]*fUpwind[21]+0.31622776601683794*fUpwind[16]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[12]*fUpwind[19]+0.3162277660168379*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[7]*fUpwind[16]+0.31622776601683794*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaDrSurf[11]; 
  Ghat[25] = 0.2020305089104422*alphaDrSurf[12]*fUpwind[26]+0.3162277660168379*alphaDrSurf[1]*fUpwind[26]+0.20203050891044216*alphaDrSurf[20]*fUpwind[25]+0.22587697572631277*alphaDrSurf[8]*fUpwind[25]+0.3162277660168379*alphaDrSurf[7]*fUpwind[25]+0.3535533905932737*alphaDrSurf[0]*fUpwind[25]+0.28284271247461906*alphaDrSurf[4]*fUpwind[24]+0.22587697572631282*alphaDrSurf[12]*fUpwind[22]+0.3535533905932737*alphaDrSurf[1]*fUpwind[22]+0.31622776601683794*alphaDrSurf[12]*fUpwind[21]+0.31622776601683794*fUpwind[15]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[11]*fUpwind[19]+0.3162277660168379*alphaDrSurf[2]*fUpwind[19]+0.31622776601683794*alphaDrSurf[4]*fUpwind[16]+0.3535533905932737*alphaDrSurf[8]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaDrSurf[12]; 
  Ghat[26] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[26]+0.22587697572631277*alphaDrSurf[8]*fUpwind[26]+0.22587697572631277*alphaDrSurf[7]*fUpwind[26]+0.3535533905932737*alphaDrSurf[0]*fUpwind[26]+0.2020305089104422*alphaDrSurf[12]*fUpwind[25]+0.3162277660168379*alphaDrSurf[1]*fUpwind[25]+0.2020305089104422*alphaDrSurf[11]*fUpwind[24]+0.3162277660168379*alphaDrSurf[2]*fUpwind[24]+0.22587697572631277*alphaDrSurf[20]*fUpwind[22]+0.3535533905932737*alphaDrSurf[7]*fUpwind[22]+0.22587697572631277*alphaDrSurf[20]*fUpwind[21]+0.3535533905932737*alphaDrSurf[8]*fUpwind[21]+0.3535533905932737*fUpwind[9]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[19]+0.3162277660168379*alphaDrSurf[11]*fUpwind[16]+0.3162277660168379*alphaDrSurf[12]*fUpwind[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[15] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[22] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[25] += 1.224744871391589*Ghat[7]*rdv2; 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[29] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[30] += 1.5811388300841898*Ghat[3]*rdv2; 
  out[31] += 1.224744871391589*Ghat[10]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[35] += 1.224744871391589*Ghat[11]*rdv2; 
  out[36] += 1.224744871391589*Ghat[12]*rdv2; 
  out[37] += 1.224744871391589*Ghat[13]*rdv2; 
  out[38] += 1.224744871391589*Ghat[14]*rdv2; 
  out[39] += 1.224744871391589*Ghat[15]*rdv2; 
  out[40] += 1.224744871391589*Ghat[16]*rdv2; 
  out[41] += 1.5811388300841895*Ghat[4]*rdv2; 
  out[42] += 1.5811388300841895*Ghat[5]*rdv2; 
  out[43] += 1.5811388300841895*Ghat[6]*rdv2; 
  out[44] += 0.7071067811865475*Ghat[20]*rdv2; 
  out[45] += 0.7071067811865475*Ghat[21]*rdv2; 
  out[46] += 0.7071067811865475*Ghat[22]*rdv2; 
  out[47] += 1.5811388300841895*Ghat[7]*rdv2; 
  out[48] += 1.5811388300841895*Ghat[8]*rdv2; 
  out[49] += 1.5811388300841895*Ghat[9]*rdv2; 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += 1.224744871391589*Ghat[18]*rdv2; 
  out[52] += 1.224744871391589*Ghat[19]*rdv2; 
  out[53] += 1.5811388300841898*Ghat[10]*rdv2; 
  out[54] += 0.7071067811865475*Ghat[23]*rdv2; 
  out[55] += 0.7071067811865475*Ghat[24]*rdv2; 
  out[56] += 0.7071067811865475*Ghat[25]*rdv2; 
  out[57] += 1.224744871391589*Ghat[20]*rdv2; 
  out[58] += 1.224744871391589*Ghat[21]*rdv2; 
  out[59] += 1.224744871391589*Ghat[22]*rdv2; 
  out[60] += 1.5811388300841898*Ghat[11]*rdv2; 
  out[61] += 1.5811388300841898*Ghat[12]*rdv2; 
  out[62] += 1.5811388300841898*Ghat[13]*rdv2; 
  out[63] += 1.5811388300841898*Ghat[14]*rdv2; 
  out[64] += 1.5811388300841898*Ghat[15]*rdv2; 
  out[65] += 1.5811388300841898*Ghat[16]*rdv2; 
  out[66] += 1.224744871391589*Ghat[23]*rdv2; 
  out[67] += 1.224744871391589*Ghat[24]*rdv2; 
  out[68] += 1.224744871391589*Ghat[25]*rdv2; 
  out[69] += 1.5811388300841895*Ghat[17]*rdv2; 
  out[70] += 1.5811388300841895*Ghat[18]*rdv2; 
  out[71] += 1.5811388300841895*Ghat[19]*rdv2; 
  out[72] += 0.7071067811865475*Ghat[26]*rdv2; 
  out[73] += 1.5811388300841895*Ghat[20]*rdv2; 
  out[74] += 1.5811388300841895*Ghat[21]*rdv2; 
  out[75] += 1.5811388300841895*Ghat[22]*rdv2; 
  out[76] += 1.224744871391589*Ghat[26]*rdv2; 
  out[77] += 1.5811388300841898*Ghat[23]*rdv2; 
  out[78] += 1.5811388300841898*Ghat[24]*rdv2; 
  out[79] += 1.5811388300841898*Ghat[25]*rdv2; 
  out[80] += 1.5811388300841895*Ghat[26]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -(0.7071067811865475*(2.0*sumNuUy[4]+(dxv[3]-2.0*w[3])*nuSum[4])); 
  alphaDrSurf[8] = -(0.7071067811865475*(2.0*sumNuUy[5]+(dxv[3]-2.0*w[3])*nuSum[5])); 
  alphaDrSurf[11] = -(0.7071067811865475*(2.0*sumNuUy[6]+(dxv[3]-2.0*w[3])*nuSum[6])); 
  alphaDrSurf[12] = -(0.7071067811865475*(2.0*sumNuUy[7]+(dxv[3]-2.0*w[3])*nuSum[7])); 
  alphaDrSurf[20] = -(0.7071067811865475*(2.0*sumNuUy[8]+(dxv[3]-2.0*w[3])*nuSum[8])); 

  if (0.28284271247461906*alphaDrSurf[20]-0.42426406871192823*alphaDrSurf[12]-0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = tensor_4x_p2_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = tensor_4x_p2_surfx4_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = tensor_4x_p2_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])+0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = tensor_4x_p2_surfx4_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = tensor_4x_p2_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = tensor_4x_p2_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]-0.42426406871192823*alphaDrSurf[12]+0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = tensor_4x_p2_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = tensor_4x_p2_surfx4_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = tensor_4x_p2_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])+0.5303300858899102*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = tensor_4x_p2_surfx4_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = tensor_4x_p2_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = tensor_4x_p2_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.44194173824159233*alphaDrSurf[20]-0.3952847075210471*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[12] = tensor_4x_p2_surfx4_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = tensor_4x_p2_surfx4_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = tensor_4x_p2_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])-0.5303300858899102*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]-0.3952847075210471*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[15] = tensor_4x_p2_surfx4_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = tensor_4x_p2_surfx4_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = tensor_4x_p2_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]+0.42426406871192823*alphaDrSurf[12]-0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[18] = tensor_4x_p2_surfx4_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = tensor_4x_p2_surfx4_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = tensor_4x_p2_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if (-(0.3535533905932734*alphaDrSurf[20])-0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[21] = tensor_4x_p2_surfx4_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = tensor_4x_p2_surfx4_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = tensor_4x_p2_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if (0.28284271247461906*alphaDrSurf[20]+0.42426406871192823*alphaDrSurf[12]+0.42426406871192857*alphaDrSurf[11]+0.31622776601683783*alphaDrSurf[8]+0.31622776601683783*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[24] = tensor_4x_p2_surfx4_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = tensor_4x_p2_surfx4_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = tensor_4x_p2_surfx4_eval_quad_node_26_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.31622776601683794*alphaDrSurf[12]*fUpwind[20]+0.31622776601683794*fUpwind[12]*alphaDrSurf[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[4]*fUpwind[11]+0.31622776601683794*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.31622776601683794*alphaDrSurf[11]*fUpwind[20]+0.31622776601683794*fUpwind[11]*alphaDrSurf[20]+0.31622776601683794*alphaDrSurf[4]*fUpwind[12]+0.31622776601683794*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[20]+0.28284271247461906*fUpwind[4]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[11]*fUpwind[12]+0.31622776601683794*alphaDrSurf[2]*fUpwind[12]+0.28284271247461906*fUpwind[11]*alphaDrSurf[12]+0.31622776601683794*fUpwind[2]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[1]*fUpwind[11]+0.31622776601683794*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.31622776601683794*alphaDrSurf[12]*fUpwind[23]+0.3162277660168379*fUpwind[18]*alphaDrSurf[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.31622776601683794*alphaDrSurf[1]*fUpwind[13]+0.31622776601683794*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = 0.31622776601683794*alphaDrSurf[11]*fUpwind[23]+0.3162277660168379*fUpwind[17]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.31622776601683794*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.31622776601683794*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[8]*fUpwind[20]+0.3535533905932737*fUpwind[8]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.22587697572631277*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.22587697572631277*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[20]+0.3535533905932737*alphaDrSurf[7]*fUpwind[20]+0.3535533905932737*fUpwind[7]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.22587697572631277*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[12]*fUpwind[25]+0.3535533905932737*alphaDrSurf[11]*fUpwind[24]+0.3535533905932737*alphaDrSurf[8]*fUpwind[22]+0.3535533905932737*alphaDrSurf[7]*fUpwind[21]+0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  Ghat[10] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[23]+0.28284271247461906*fUpwind[10]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.31622776601683794*alphaDrSurf[4]*fUpwind[14]+0.31622776601683794*alphaDrSurf[4]*fUpwind[13]+0.31622776601683794*fUpwind[6]*alphaDrSurf[12]+0.31622776601683794*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.20203050891044216*alphaDrSurf[11]*fUpwind[20]+0.31622776601683794*alphaDrSurf[2]*fUpwind[20]+0.20203050891044216*fUpwind[11]*alphaDrSurf[20]+0.31622776601683794*fUpwind[2]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[12]+0.28284271247461906*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.22587697572631277*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.22587697572631277*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.31622776601683794*alphaDrSurf[1]*fUpwind[4]+0.31622776601683794*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[12] = 0.20203050891044216*alphaDrSurf[12]*fUpwind[20]+0.31622776601683794*alphaDrSurf[1]*fUpwind[20]+0.20203050891044216*fUpwind[12]*alphaDrSurf[20]+0.31622776601683794*fUpwind[1]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.22587697572631277*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.28284271247461906*alphaDrSurf[4]*fUpwind[11]+0.28284271247461906*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.31622776601683794*alphaDrSurf[2]*fUpwind[4]+0.31622776601683794*fUpwind[2]*alphaDrSurf[4]; 
  Ghat[13] = 0.22587697572631282*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[8]*fUpwind[23]+0.3535533905932737*fUpwind[14]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.22587697572631277*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.22587697572631277*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.31622776601683794*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.31622776601683794*alphaDrSurf[1]*fUpwind[5]; 
  Ghat[14] = 0.22587697572631282*alphaDrSurf[20]*fUpwind[23]+0.3535533905932737*alphaDrSurf[7]*fUpwind[23]+0.3535533905932737*fUpwind[13]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.22587697572631277*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.31622776601683794*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.31622776601683794*alphaDrSurf[2]*fUpwind[6]; 
  Ghat[15] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[26]+0.31622776601683794*alphaDrSurf[20]*fUpwind[25]+0.3535533905932737*alphaDrSurf[8]*fUpwind[25]+0.31622776601683794*alphaDrSurf[4]*fUpwind[24]+0.3535533905932737*alphaDrSurf[12]*fUpwind[22]+0.31622776601683794*alphaDrSurf[1]*fUpwind[21]+0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  Ghat[16] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[26]+0.31622776601683794*alphaDrSurf[4]*fUpwind[25]+0.31622776601683794*alphaDrSurf[20]*fUpwind[24]+0.3535533905932737*alphaDrSurf[7]*fUpwind[24]+0.31622776601683794*alphaDrSurf[2]*fUpwind[22]+0.3535533905932737*alphaDrSurf[11]*fUpwind[21]+0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[17] = 0.2020305089104422*alphaDrSurf[11]*fUpwind[23]+0.3162277660168379*alphaDrSurf[2]*fUpwind[23]+0.20203050891044216*fUpwind[17]*alphaDrSurf[20]+0.3162277660168379*fUpwind[6]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.22587697572631277*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.22587697572631277*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  Ghat[18] = 0.2020305089104422*alphaDrSurf[12]*fUpwind[23]+0.3162277660168379*alphaDrSurf[1]*fUpwind[23]+0.20203050891044216*fUpwind[18]*alphaDrSurf[20]+0.3162277660168379*fUpwind[5]*alphaDrSurf[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.28284271247461906*alphaDrSurf[4]*fUpwind[17]+0.22587697572631277*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  Ghat[19] = 0.28284271247461906*alphaDrSurf[4]*fUpwind[26]+0.282842712474619*alphaDrSurf[11]*fUpwind[25]+0.3162277660168379*alphaDrSurf[2]*fUpwind[25]+0.282842712474619*alphaDrSurf[12]*fUpwind[24]+0.3162277660168379*alphaDrSurf[1]*fUpwind[24]+0.3162277660168379*alphaDrSurf[4]*fUpwind[22]+0.3162277660168379*alphaDrSurf[4]*fUpwind[21]+0.28284271247461906*fUpwind[19]*alphaDrSurf[20]+0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 
  Ghat[20] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[20]+0.22587697572631277*alphaDrSurf[8]*fUpwind[20]+0.22587697572631277*alphaDrSurf[7]*fUpwind[20]+0.3535533905932737*alphaDrSurf[0]*fUpwind[20]+0.22587697572631277*fUpwind[8]*alphaDrSurf[20]+0.22587697572631277*fUpwind[7]*alphaDrSurf[20]+0.3535533905932737*fUpwind[0]*alphaDrSurf[20]+0.20203050891044216*alphaDrSurf[12]*fUpwind[12]+0.31622776601683794*alphaDrSurf[1]*fUpwind[12]+0.31622776601683794*fUpwind[1]*alphaDrSurf[12]+0.20203050891044216*alphaDrSurf[11]*fUpwind[11]+0.31622776601683794*alphaDrSurf[2]*fUpwind[11]+0.31622776601683794*fUpwind[2]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[7]*fUpwind[8]+0.3535533905932737*fUpwind[7]*alphaDrSurf[8]+0.28284271247461906*alphaDrSurf[4]*fUpwind[4]; 
  Ghat[21] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[8]*fUpwind[26]+0.31622776601683794*alphaDrSurf[12]*fUpwind[25]+0.22587697572631282*alphaDrSurf[11]*fUpwind[24]+0.3535533905932737*alphaDrSurf[2]*fUpwind[24]+0.3535533905932737*alphaDrSurf[20]*fUpwind[22]+0.22587697572631277*alphaDrSurf[7]*fUpwind[21]+0.3535533905932737*alphaDrSurf[0]*fUpwind[21]+0.3162277660168379*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[11]*fUpwind[16]+0.31622776601683794*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[7]*fUpwind[9]; 
  Ghat[22] = 0.22587697572631277*alphaDrSurf[20]*fUpwind[26]+0.3535533905932737*alphaDrSurf[7]*fUpwind[26]+0.22587697572631282*alphaDrSurf[12]*fUpwind[25]+0.3535533905932737*alphaDrSurf[1]*fUpwind[25]+0.31622776601683794*alphaDrSurf[11]*fUpwind[24]+0.22587697572631277*alphaDrSurf[8]*fUpwind[22]+0.3535533905932737*alphaDrSurf[0]*fUpwind[22]+0.3535533905932737*alphaDrSurf[20]*fUpwind[21]+0.3162277660168379*alphaDrSurf[4]*fUpwind[19]+0.31622776601683794*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[12]*fUpwind[15]+0.3535533905932737*alphaDrSurf[8]*fUpwind[9]; 
  Ghat[23] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[23]+0.22587697572631277*alphaDrSurf[8]*fUpwind[23]+0.22587697572631277*alphaDrSurf[7]*fUpwind[23]+0.3535533905932737*alphaDrSurf[0]*fUpwind[23]+0.22587697572631282*fUpwind[14]*alphaDrSurf[20]+0.22587697572631282*fUpwind[13]*alphaDrSurf[20]+0.3535533905932737*fUpwind[3]*alphaDrSurf[20]+0.2020305089104422*alphaDrSurf[12]*fUpwind[18]+0.3162277660168379*alphaDrSurf[1]*fUpwind[18]+0.2020305089104422*alphaDrSurf[11]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[17]+0.3535533905932737*alphaDrSurf[7]*fUpwind[14]+0.3535533905932737*alphaDrSurf[8]*fUpwind[13]+0.31622776601683794*fUpwind[5]*alphaDrSurf[12]+0.31622776601683794*fUpwind[6]*alphaDrSurf[11]+0.28284271247461906*alphaDrSurf[4]*fUpwind[10]; 
  Ghat[24] = 0.2020305089104422*alphaDrSurf[11]*fUpwind[26]+0.3162277660168379*alphaDrSurf[2]*fUpwind[26]+0.28284271247461906*alphaDrSurf[4]*fUpwind[25]+0.20203050891044216*alphaDrSurf[20]*fUpwind[24]+0.3162277660168379*alphaDrSurf[8]*fUpwind[24]+0.22587697572631277*alphaDrSurf[7]*fUpwind[24]+0.3535533905932737*alphaDrSurf[0]*fUpwind[24]+0.31622776601683794*alphaDrSurf[11]*fUpwind[22]+0.22587697572631282*alphaDrSurf[11]*fUpwind[21]+0.3535533905932737*alphaDrSurf[2]*fUpwind[21]+0.31622776601683794*fUpwind[16]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[12]*fUpwind[19]+0.3162277660168379*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[7]*fUpwind[16]+0.31622776601683794*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaDrSurf[11]; 
  Ghat[25] = 0.2020305089104422*alphaDrSurf[12]*fUpwind[26]+0.3162277660168379*alphaDrSurf[1]*fUpwind[26]+0.20203050891044216*alphaDrSurf[20]*fUpwind[25]+0.22587697572631277*alphaDrSurf[8]*fUpwind[25]+0.3162277660168379*alphaDrSurf[7]*fUpwind[25]+0.3535533905932737*alphaDrSurf[0]*fUpwind[25]+0.28284271247461906*alphaDrSurf[4]*fUpwind[24]+0.22587697572631282*alphaDrSurf[12]*fUpwind[22]+0.3535533905932737*alphaDrSurf[1]*fUpwind[22]+0.31622776601683794*alphaDrSurf[12]*fUpwind[21]+0.31622776601683794*fUpwind[15]*alphaDrSurf[20]+0.282842712474619*alphaDrSurf[11]*fUpwind[19]+0.3162277660168379*alphaDrSurf[2]*fUpwind[19]+0.31622776601683794*alphaDrSurf[4]*fUpwind[16]+0.3535533905932737*alphaDrSurf[8]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaDrSurf[12]; 
  Ghat[26] = 0.14430750636460155*alphaDrSurf[20]*fUpwind[26]+0.22587697572631277*alphaDrSurf[8]*fUpwind[26]+0.22587697572631277*alphaDrSurf[7]*fUpwind[26]+0.3535533905932737*alphaDrSurf[0]*fUpwind[26]+0.2020305089104422*alphaDrSurf[12]*fUpwind[25]+0.3162277660168379*alphaDrSurf[1]*fUpwind[25]+0.2020305089104422*alphaDrSurf[11]*fUpwind[24]+0.3162277660168379*alphaDrSurf[2]*fUpwind[24]+0.22587697572631277*alphaDrSurf[20]*fUpwind[22]+0.3535533905932737*alphaDrSurf[7]*fUpwind[22]+0.22587697572631277*alphaDrSurf[20]*fUpwind[21]+0.3535533905932737*alphaDrSurf[8]*fUpwind[21]+0.3535533905932737*fUpwind[9]*alphaDrSurf[20]+0.28284271247461906*alphaDrSurf[4]*fUpwind[19]+0.3162277660168379*alphaDrSurf[11]*fUpwind[16]+0.3162277660168379*alphaDrSurf[12]*fUpwind[15]; 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[6] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[7] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[12] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[13] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[14] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[15] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[20] += -(0.7071067811865475*Ghat[12]*rdv2); 
  out[21] += -(0.7071067811865475*Ghat[13]*rdv2); 
  out[22] += -(0.7071067811865475*Ghat[14]*rdv2); 
  out[23] += -(0.7071067811865475*Ghat[15]*rdv2); 
  out[24] += -(0.7071067811865475*Ghat[16]*rdv2); 
  out[25] += 1.224744871391589*Ghat[7]*rdv2; 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[29] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[30] += -(1.5811388300841898*Ghat[3]*rdv2); 
  out[31] += 1.224744871391589*Ghat[10]*rdv2; 
  out[32] += -(0.7071067811865475*Ghat[17]*rdv2); 
  out[33] += -(0.7071067811865475*Ghat[18]*rdv2); 
  out[34] += -(0.7071067811865475*Ghat[19]*rdv2); 
  out[35] += 1.224744871391589*Ghat[11]*rdv2; 
  out[36] += 1.224744871391589*Ghat[12]*rdv2; 
  out[37] += 1.224744871391589*Ghat[13]*rdv2; 
  out[38] += 1.224744871391589*Ghat[14]*rdv2; 
  out[39] += 1.224744871391589*Ghat[15]*rdv2; 
  out[40] += 1.224744871391589*Ghat[16]*rdv2; 
  out[41] += -(1.5811388300841895*Ghat[4]*rdv2); 
  out[42] += -(1.5811388300841895*Ghat[5]*rdv2); 
  out[43] += -(1.5811388300841895*Ghat[6]*rdv2); 
  out[44] += -(0.7071067811865475*Ghat[20]*rdv2); 
  out[45] += -(0.7071067811865475*Ghat[21]*rdv2); 
  out[46] += -(0.7071067811865475*Ghat[22]*rdv2); 
  out[47] += -(1.5811388300841895*Ghat[7]*rdv2); 
  out[48] += -(1.5811388300841895*Ghat[8]*rdv2); 
  out[49] += -(1.5811388300841895*Ghat[9]*rdv2); 
  out[50] += 1.224744871391589*Ghat[17]*rdv2; 
  out[51] += 1.224744871391589*Ghat[18]*rdv2; 
  out[52] += 1.224744871391589*Ghat[19]*rdv2; 
  out[53] += -(1.5811388300841898*Ghat[10]*rdv2); 
  out[54] += -(0.7071067811865475*Ghat[23]*rdv2); 
  out[55] += -(0.7071067811865475*Ghat[24]*rdv2); 
  out[56] += -(0.7071067811865475*Ghat[25]*rdv2); 
  out[57] += 1.224744871391589*Ghat[20]*rdv2; 
  out[58] += 1.224744871391589*Ghat[21]*rdv2; 
  out[59] += 1.224744871391589*Ghat[22]*rdv2; 
  out[60] += -(1.5811388300841898*Ghat[11]*rdv2); 
  out[61] += -(1.5811388300841898*Ghat[12]*rdv2); 
  out[62] += -(1.5811388300841898*Ghat[13]*rdv2); 
  out[63] += -(1.5811388300841898*Ghat[14]*rdv2); 
  out[64] += -(1.5811388300841898*Ghat[15]*rdv2); 
  out[65] += -(1.5811388300841898*Ghat[16]*rdv2); 
  out[66] += 1.224744871391589*Ghat[23]*rdv2; 
  out[67] += 1.224744871391589*Ghat[24]*rdv2; 
  out[68] += 1.224744871391589*Ghat[25]*rdv2; 
  out[69] += -(1.5811388300841895*Ghat[17]*rdv2); 
  out[70] += -(1.5811388300841895*Ghat[18]*rdv2); 
  out[71] += -(1.5811388300841895*Ghat[19]*rdv2); 
  out[72] += -(0.7071067811865475*Ghat[26]*rdv2); 
  out[73] += -(1.5811388300841895*Ghat[20]*rdv2); 
  out[74] += -(1.5811388300841895*Ghat[21]*rdv2); 
  out[75] += -(1.5811388300841895*Ghat[22]*rdv2); 
  out[76] += 1.224744871391589*Ghat[26]*rdv2; 
  out[77] += -(1.5811388300841898*Ghat[23]*rdv2); 
  out[78] += -(1.5811388300841898*Ghat[24]*rdv2); 
  out[79] += -(1.5811388300841898*Ghat[25]*rdv2); 
  out[80] += -(1.5811388300841895*Ghat[26]*rdv2); 

  } 

  return 0.;

} 
