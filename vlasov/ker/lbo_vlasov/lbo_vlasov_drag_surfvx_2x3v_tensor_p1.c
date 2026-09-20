#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]: cell-center coordinates. 
  // dxv[5]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[16]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[36] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[2]-1.0*nuSum[0]*dxv[2]-2.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[2]-1.0*nuSum[1]*dxv[2]-2.0*sumNuUx[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]-1.0*dxv[2]*nuSum[2]; 
  alphaDrSurf_l[5] = -(2.0*sumNuUx[3])+2.0*w[2]*nuSum[3]-1.0*dxv[2]*nuSum[3]; 

  double alphaDrSurf_r[36] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[2]+nuSum[0]*dxv[2]-2.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[2]+nuSum[1]*dxv[2]-2.0*sumNuUx[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]+dxv[2]*nuSum[2]; 
  alphaDrSurf_r[5] = -(2.0*sumNuUx[3])+2.0*w[2]*nuSum[3]+dxv[2]*nuSum[3]; 

  double fUpwindQuad_l[36] = {0.0};
  double fUpwindQuad_r[36] = {0.0};
  double fUpwind_l[36] = {0.0};
  double fUpwind_r[36] = {0.0};
  double Ghat_l[36] = {0.0}; 
  double Ghat_r[36] = {0.0}; 

  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fc); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fc); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fr); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fr); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (-alphaDrSurf_l[5]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fl); 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fl); 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fl); 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fl); 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fl); 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fl); 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fc); 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fc); 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fc); 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fc); 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fc); 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fc); 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fc); 
  } 
  if (-alphaDrSurf_r[5]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fc); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fc); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fc); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fc); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fc); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fc); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fr); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fr); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fr); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fr); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fr); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fr); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fr); 
  } 
  if (-alphaDrSurf_l[5]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fl); 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fl); 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fl); 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fl); 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fl); 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fl); 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fl); 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fl); 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fc); 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fc); 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fc); 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fc); 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fc); 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fc); 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fc); 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fc); 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fc); 
  } 
  if (-alphaDrSurf_r[5]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fc); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fc); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fc); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fc); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fc); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fc); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fc); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fc); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fr); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fr); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fr); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fr); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fr); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fr); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fr); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fr); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fr); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fl); 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fl); 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fl); 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fl); 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fl); 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fl); 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fl); 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fl); 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fl); 
  } else { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fc); 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fc); 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fc); 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fc); 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fc); 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fc); 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fc); 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fc); 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fc); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fc); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fc); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fc); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fc); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fc); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fc); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fc); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fc); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fr); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fr); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fr); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fr); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fr); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fr); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fr); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fr); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[5]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[5]+fUpwind_l[2]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[5]+fUpwind_l[1]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[13]+alphaDrSurf_l[2]*fUpwind_l[7]+alphaDrSurf_l[1]*fUpwind_l[6]+alphaDrSurf_l[0]*fUpwind_l[3]); 
  Ghat_l[4] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[14]+alphaDrSurf_l[2]*fUpwind_l[9]+alphaDrSurf_l[1]*fUpwind_l[8]+alphaDrSurf_l[0]*fUpwind_l[4]); 
  Ghat_l[5] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[5]+fUpwind_l[0]*alphaDrSurf_l[5]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[6] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[13]+alphaDrSurf_l[5]*fUpwind_l[7]+alphaDrSurf_l[0]*fUpwind_l[6]+alphaDrSurf_l[1]*fUpwind_l[3]); 
  Ghat_l[7] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[13]+alphaDrSurf_l[0]*fUpwind_l[7]+alphaDrSurf_l[5]*fUpwind_l[6]+alphaDrSurf_l[2]*fUpwind_l[3]); 
  Ghat_l[8] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[14]+alphaDrSurf_l[5]*fUpwind_l[9]+alphaDrSurf_l[0]*fUpwind_l[8]+alphaDrSurf_l[1]*fUpwind_l[4]); 
  Ghat_l[9] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[14]+alphaDrSurf_l[0]*fUpwind_l[9]+alphaDrSurf_l[5]*fUpwind_l[8]+alphaDrSurf_l[2]*fUpwind_l[4]); 
  Ghat_l[10] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[23]+alphaDrSurf_l[2]*fUpwind_l[16]+alphaDrSurf_l[1]*fUpwind_l[15]+alphaDrSurf_l[0]*fUpwind_l[10]); 
  Ghat_l[11] = 0.016666666666666666*(15.0*alphaDrSurf_l[5]*fUpwind_l[24]+15.000000000000002*(alphaDrSurf_l[2]*fUpwind_l[18]+alphaDrSurf_l[1]*fUpwind_l[17])+15.0*alphaDrSurf_l[0]*fUpwind_l[11]); 
  Ghat_l[12] = 0.016666666666666666*(15.0*alphaDrSurf_l[5]*fUpwind_l[27]+15.000000000000002*(alphaDrSurf_l[2]*fUpwind_l[21]+alphaDrSurf_l[1]*fUpwind_l[20])+15.0*alphaDrSurf_l[0]*fUpwind_l[12]); 
  Ghat_l[13] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[13]+alphaDrSurf_l[1]*fUpwind_l[7]+alphaDrSurf_l[2]*fUpwind_l[6]+fUpwind_l[3]*alphaDrSurf_l[5]); 
  Ghat_l[14] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[14]+alphaDrSurf_l[1]*fUpwind_l[9]+alphaDrSurf_l[2]*fUpwind_l[8]+fUpwind_l[4]*alphaDrSurf_l[5]); 
  Ghat_l[15] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[23]+alphaDrSurf_l[5]*fUpwind_l[16]+alphaDrSurf_l[0]*fUpwind_l[15]+alphaDrSurf_l[1]*fUpwind_l[10]); 
  Ghat_l[16] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[23]+alphaDrSurf_l[0]*fUpwind_l[16]+alphaDrSurf_l[5]*fUpwind_l[15]+alphaDrSurf_l[2]*fUpwind_l[10]); 
  Ghat_l[17] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[24]+15.0*(alphaDrSurf_l[5]*fUpwind_l[18]+alphaDrSurf_l[0]*fUpwind_l[17])+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[11]); 
  Ghat_l[18] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[24]+15.0*(alphaDrSurf_l[0]*fUpwind_l[18]+alphaDrSurf_l[5]*fUpwind_l[17])+15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[11]); 
  Ghat_l[19] = 0.016666666666666666*(15.0*alphaDrSurf_l[5]*fUpwind_l[31]+15.000000000000002*(alphaDrSurf_l[2]*fUpwind_l[26]+alphaDrSurf_l[1]*fUpwind_l[25])+15.0*alphaDrSurf_l[0]*fUpwind_l[19]); 
  Ghat_l[20] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[27]+15.0*(alphaDrSurf_l[5]*fUpwind_l[21]+alphaDrSurf_l[0]*fUpwind_l[20])+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[12]); 
  Ghat_l[21] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[27]+15.0*(alphaDrSurf_l[0]*fUpwind_l[21]+alphaDrSurf_l[5]*fUpwind_l[20])+15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[12]); 
  Ghat_l[22] = 0.016666666666666666*(15.0*alphaDrSurf_l[5]*fUpwind_l[32]+15.000000000000002*(alphaDrSurf_l[2]*fUpwind_l[29]+alphaDrSurf_l[1]*fUpwind_l[28])+15.0*alphaDrSurf_l[0]*fUpwind_l[22]); 
  Ghat_l[23] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[23]+alphaDrSurf_l[1]*fUpwind_l[16]+alphaDrSurf_l[2]*fUpwind_l[15]+alphaDrSurf_l[5]*fUpwind_l[10]); 
  Ghat_l[24] = 0.016666666666666666*(15.0*alphaDrSurf_l[0]*fUpwind_l[24]+15.000000000000002*(alphaDrSurf_l[1]*fUpwind_l[18]+alphaDrSurf_l[2]*fUpwind_l[17])+15.0*alphaDrSurf_l[5]*fUpwind_l[11]); 
  Ghat_l[25] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[31]+15.0*(alphaDrSurf_l[5]*fUpwind_l[26]+alphaDrSurf_l[0]*fUpwind_l[25])+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[19]); 
  Ghat_l[26] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[31]+15.0*(alphaDrSurf_l[0]*fUpwind_l[26]+alphaDrSurf_l[5]*fUpwind_l[25])+15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[19]); 
  Ghat_l[27] = 0.016666666666666666*(15.0*alphaDrSurf_l[0]*fUpwind_l[27]+15.000000000000002*(alphaDrSurf_l[1]*fUpwind_l[21]+alphaDrSurf_l[2]*fUpwind_l[20])+15.0*alphaDrSurf_l[5]*fUpwind_l[12]); 
  Ghat_l[28] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[32]+15.0*(alphaDrSurf_l[5]*fUpwind_l[29]+alphaDrSurf_l[0]*fUpwind_l[28])+15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[22]); 
  Ghat_l[29] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_l[1]*fUpwind_l[32]+15.0*(alphaDrSurf_l[0]*fUpwind_l[29]+alphaDrSurf_l[5]*fUpwind_l[28])+15.000000000000002*alphaDrSurf_l[2]*fUpwind_l[22]); 
  Ghat_l[30] = 0.25*(alphaDrSurf_l[5]*fUpwind_l[35]+alphaDrSurf_l[2]*fUpwind_l[34]+alphaDrSurf_l[1]*fUpwind_l[33]+alphaDrSurf_l[0]*fUpwind_l[30]); 
  Ghat_l[31] = 0.016666666666666666*(15.0*alphaDrSurf_l[0]*fUpwind_l[31]+15.000000000000002*(alphaDrSurf_l[1]*fUpwind_l[26]+alphaDrSurf_l[2]*fUpwind_l[25])+15.0*alphaDrSurf_l[5]*fUpwind_l[19]); 
  Ghat_l[32] = 0.016666666666666666*(15.0*alphaDrSurf_l[0]*fUpwind_l[32]+15.000000000000002*(alphaDrSurf_l[1]*fUpwind_l[29]+alphaDrSurf_l[2]*fUpwind_l[28])+15.0*alphaDrSurf_l[5]*fUpwind_l[22]); 
  Ghat_l[33] = 0.25*(alphaDrSurf_l[2]*fUpwind_l[35]+alphaDrSurf_l[5]*fUpwind_l[34]+alphaDrSurf_l[0]*fUpwind_l[33]+alphaDrSurf_l[1]*fUpwind_l[30]); 
  Ghat_l[34] = 0.25*(alphaDrSurf_l[1]*fUpwind_l[35]+alphaDrSurf_l[0]*fUpwind_l[34]+alphaDrSurf_l[5]*fUpwind_l[33]+alphaDrSurf_l[2]*fUpwind_l[30]); 
  Ghat_l[35] = 0.25*(alphaDrSurf_l[0]*fUpwind_l[35]+alphaDrSurf_l[1]*fUpwind_l[34]+alphaDrSurf_l[2]*fUpwind_l[33]+alphaDrSurf_l[5]*fUpwind_l[30]); 

  Ghat_r[0] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[5]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[5]+fUpwind_r[2]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[5]+fUpwind_r[1]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[13]+alphaDrSurf_r[2]*fUpwind_r[7]+alphaDrSurf_r[1]*fUpwind_r[6]+alphaDrSurf_r[0]*fUpwind_r[3]); 
  Ghat_r[4] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[14]+alphaDrSurf_r[2]*fUpwind_r[9]+alphaDrSurf_r[1]*fUpwind_r[8]+alphaDrSurf_r[0]*fUpwind_r[4]); 
  Ghat_r[5] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[5]+fUpwind_r[0]*alphaDrSurf_r[5]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[6] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[13]+alphaDrSurf_r[5]*fUpwind_r[7]+alphaDrSurf_r[0]*fUpwind_r[6]+alphaDrSurf_r[1]*fUpwind_r[3]); 
  Ghat_r[7] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[13]+alphaDrSurf_r[0]*fUpwind_r[7]+alphaDrSurf_r[5]*fUpwind_r[6]+alphaDrSurf_r[2]*fUpwind_r[3]); 
  Ghat_r[8] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[14]+alphaDrSurf_r[5]*fUpwind_r[9]+alphaDrSurf_r[0]*fUpwind_r[8]+alphaDrSurf_r[1]*fUpwind_r[4]); 
  Ghat_r[9] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[14]+alphaDrSurf_r[0]*fUpwind_r[9]+alphaDrSurf_r[5]*fUpwind_r[8]+alphaDrSurf_r[2]*fUpwind_r[4]); 
  Ghat_r[10] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[23]+alphaDrSurf_r[2]*fUpwind_r[16]+alphaDrSurf_r[1]*fUpwind_r[15]+alphaDrSurf_r[0]*fUpwind_r[10]); 
  Ghat_r[11] = 0.016666666666666666*(15.0*alphaDrSurf_r[5]*fUpwind_r[24]+15.000000000000002*(alphaDrSurf_r[2]*fUpwind_r[18]+alphaDrSurf_r[1]*fUpwind_r[17])+15.0*alphaDrSurf_r[0]*fUpwind_r[11]); 
  Ghat_r[12] = 0.016666666666666666*(15.0*alphaDrSurf_r[5]*fUpwind_r[27]+15.000000000000002*(alphaDrSurf_r[2]*fUpwind_r[21]+alphaDrSurf_r[1]*fUpwind_r[20])+15.0*alphaDrSurf_r[0]*fUpwind_r[12]); 
  Ghat_r[13] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[13]+alphaDrSurf_r[1]*fUpwind_r[7]+alphaDrSurf_r[2]*fUpwind_r[6]+fUpwind_r[3]*alphaDrSurf_r[5]); 
  Ghat_r[14] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[14]+alphaDrSurf_r[1]*fUpwind_r[9]+alphaDrSurf_r[2]*fUpwind_r[8]+fUpwind_r[4]*alphaDrSurf_r[5]); 
  Ghat_r[15] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[23]+alphaDrSurf_r[5]*fUpwind_r[16]+alphaDrSurf_r[0]*fUpwind_r[15]+alphaDrSurf_r[1]*fUpwind_r[10]); 
  Ghat_r[16] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[23]+alphaDrSurf_r[0]*fUpwind_r[16]+alphaDrSurf_r[5]*fUpwind_r[15]+alphaDrSurf_r[2]*fUpwind_r[10]); 
  Ghat_r[17] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[24]+15.0*(alphaDrSurf_r[5]*fUpwind_r[18]+alphaDrSurf_r[0]*fUpwind_r[17])+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[11]); 
  Ghat_r[18] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[24]+15.0*(alphaDrSurf_r[0]*fUpwind_r[18]+alphaDrSurf_r[5]*fUpwind_r[17])+15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[11]); 
  Ghat_r[19] = 0.016666666666666666*(15.0*alphaDrSurf_r[5]*fUpwind_r[31]+15.000000000000002*(alphaDrSurf_r[2]*fUpwind_r[26]+alphaDrSurf_r[1]*fUpwind_r[25])+15.0*alphaDrSurf_r[0]*fUpwind_r[19]); 
  Ghat_r[20] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[27]+15.0*(alphaDrSurf_r[5]*fUpwind_r[21]+alphaDrSurf_r[0]*fUpwind_r[20])+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[12]); 
  Ghat_r[21] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[27]+15.0*(alphaDrSurf_r[0]*fUpwind_r[21]+alphaDrSurf_r[5]*fUpwind_r[20])+15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[12]); 
  Ghat_r[22] = 0.016666666666666666*(15.0*alphaDrSurf_r[5]*fUpwind_r[32]+15.000000000000002*(alphaDrSurf_r[2]*fUpwind_r[29]+alphaDrSurf_r[1]*fUpwind_r[28])+15.0*alphaDrSurf_r[0]*fUpwind_r[22]); 
  Ghat_r[23] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[23]+alphaDrSurf_r[1]*fUpwind_r[16]+alphaDrSurf_r[2]*fUpwind_r[15]+alphaDrSurf_r[5]*fUpwind_r[10]); 
  Ghat_r[24] = 0.016666666666666666*(15.0*alphaDrSurf_r[0]*fUpwind_r[24]+15.000000000000002*(alphaDrSurf_r[1]*fUpwind_r[18]+alphaDrSurf_r[2]*fUpwind_r[17])+15.0*alphaDrSurf_r[5]*fUpwind_r[11]); 
  Ghat_r[25] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[31]+15.0*(alphaDrSurf_r[5]*fUpwind_r[26]+alphaDrSurf_r[0]*fUpwind_r[25])+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[19]); 
  Ghat_r[26] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[31]+15.0*(alphaDrSurf_r[0]*fUpwind_r[26]+alphaDrSurf_r[5]*fUpwind_r[25])+15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[19]); 
  Ghat_r[27] = 0.016666666666666666*(15.0*alphaDrSurf_r[0]*fUpwind_r[27]+15.000000000000002*(alphaDrSurf_r[1]*fUpwind_r[21]+alphaDrSurf_r[2]*fUpwind_r[20])+15.0*alphaDrSurf_r[5]*fUpwind_r[12]); 
  Ghat_r[28] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[32]+15.0*(alphaDrSurf_r[5]*fUpwind_r[29]+alphaDrSurf_r[0]*fUpwind_r[28])+15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[22]); 
  Ghat_r[29] = 0.016666666666666666*(15.000000000000002*alphaDrSurf_r[1]*fUpwind_r[32]+15.0*(alphaDrSurf_r[0]*fUpwind_r[29]+alphaDrSurf_r[5]*fUpwind_r[28])+15.000000000000002*alphaDrSurf_r[2]*fUpwind_r[22]); 
  Ghat_r[30] = 0.25*(alphaDrSurf_r[5]*fUpwind_r[35]+alphaDrSurf_r[2]*fUpwind_r[34]+alphaDrSurf_r[1]*fUpwind_r[33]+alphaDrSurf_r[0]*fUpwind_r[30]); 
  Ghat_r[31] = 0.016666666666666666*(15.0*alphaDrSurf_r[0]*fUpwind_r[31]+15.000000000000002*(alphaDrSurf_r[1]*fUpwind_r[26]+alphaDrSurf_r[2]*fUpwind_r[25])+15.0*alphaDrSurf_r[5]*fUpwind_r[19]); 
  Ghat_r[32] = 0.016666666666666666*(15.0*alphaDrSurf_r[0]*fUpwind_r[32]+15.000000000000002*(alphaDrSurf_r[1]*fUpwind_r[29]+alphaDrSurf_r[2]*fUpwind_r[28])+15.0*alphaDrSurf_r[5]*fUpwind_r[22]); 
  Ghat_r[33] = 0.25*(alphaDrSurf_r[2]*fUpwind_r[35]+alphaDrSurf_r[5]*fUpwind_r[34]+alphaDrSurf_r[0]*fUpwind_r[33]+alphaDrSurf_r[1]*fUpwind_r[30]); 
  Ghat_r[34] = 0.25*(alphaDrSurf_r[1]*fUpwind_r[35]+alphaDrSurf_r[0]*fUpwind_r[34]+alphaDrSurf_r[5]*fUpwind_r[33]+alphaDrSurf_r[2]*fUpwind_r[30]); 
  Ghat_r[35] = 0.25*(alphaDrSurf_r[0]*fUpwind_r[35]+alphaDrSurf_r[1]*fUpwind_r[34]+alphaDrSurf_r[2]*fUpwind_r[33]+alphaDrSurf_r[5]*fUpwind_r[30]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[8] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[11] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[15] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[16] += 1.5811388300841895*Ghat_r[0]*rdv2-1.5811388300841895*Ghat_l[0]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_r[12]*rdv2-0.7071067811865475*Ghat_l[12]*rdv2; 
  out[19] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[13]*rdv2-0.7071067811865475*Ghat_l[13]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[22] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_r[14]*rdv2-0.7071067811865475*Ghat_l[14]*rdv2; 
  out[24] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[25] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[26] += 0.7071067811865475*Ghat_r[15]*rdv2-0.7071067811865475*Ghat_l[15]*rdv2; 
  out[27] += 0.7071067811865475*Ghat_r[16]*rdv2-0.7071067811865475*Ghat_l[16]*rdv2; 
  out[28] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[29] += 1.5811388300841898*Ghat_r[1]*rdv2-1.5811388300841898*Ghat_l[1]*rdv2; 
  out[30] += 1.5811388300841898*Ghat_r[2]*rdv2-1.5811388300841898*Ghat_l[2]*rdv2; 
  out[31] += 1.5811388300841898*Ghat_r[3]*rdv2-1.5811388300841898*Ghat_l[3]*rdv2; 
  out[32] += 0.7071067811865475*Ghat_r[17]*rdv2-0.7071067811865475*Ghat_l[17]*rdv2; 
  out[33] += 0.7071067811865475*Ghat_r[18]*rdv2-0.7071067811865475*Ghat_l[18]*rdv2; 
  out[34] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[35] += 1.5811388300841898*Ghat_r[4]*rdv2-1.5811388300841898*Ghat_l[4]*rdv2; 
  out[36] += 0.7071067811865475*Ghat_r[19]*rdv2-0.7071067811865475*Ghat_l[19]*rdv2; 
  out[37] += 0.7071067811865475*Ghat_r[20]*rdv2-0.7071067811865475*Ghat_l[20]*rdv2; 
  out[38] += 0.7071067811865475*Ghat_r[21]*rdv2-0.7071067811865475*Ghat_l[21]*rdv2; 
  out[39] += 1.224744871391589*Ghat_r[12]*rdv2+1.224744871391589*Ghat_l[12]*rdv2; 
  out[40] += 0.7071067811865475*Ghat_r[22]*rdv2-0.7071067811865475*Ghat_l[22]*rdv2; 
  out[41] += 1.224744871391589*Ghat_r[13]*rdv2+1.224744871391589*Ghat_l[13]*rdv2; 
  out[42] += 1.224744871391589*Ghat_r[14]*rdv2+1.224744871391589*Ghat_l[14]*rdv2; 
  out[43] += 0.7071067811865475*Ghat_r[23]*rdv2-0.7071067811865475*Ghat_l[23]*rdv2; 
  out[44] += 1.224744871391589*Ghat_r[15]*rdv2+1.224744871391589*Ghat_l[15]*rdv2; 
  out[45] += 1.224744871391589*Ghat_r[16]*rdv2+1.224744871391589*Ghat_l[16]*rdv2; 
  out[46] += 1.5811388300841895*Ghat_r[5]*rdv2-1.5811388300841895*Ghat_l[5]*rdv2; 
  out[47] += 1.5811388300841895*Ghat_r[6]*rdv2-1.5811388300841895*Ghat_l[6]*rdv2; 
  out[48] += 1.5811388300841895*Ghat_r[7]*rdv2-1.5811388300841895*Ghat_l[7]*rdv2; 
  out[49] += 0.7071067811865475*Ghat_r[24]*rdv2-0.7071067811865475*Ghat_l[24]*rdv2; 
  out[50] += 1.224744871391589*Ghat_r[17]*rdv2+1.224744871391589*Ghat_l[17]*rdv2; 
  out[51] += 1.224744871391589*Ghat_r[18]*rdv2+1.224744871391589*Ghat_l[18]*rdv2; 
  out[52] += 1.5811388300841895*Ghat_r[8]*rdv2-1.5811388300841895*Ghat_l[8]*rdv2; 
  out[53] += 1.5811388300841895*Ghat_r[9]*rdv2-1.5811388300841895*Ghat_l[9]*rdv2; 
  out[54] += 1.5811388300841895*Ghat_r[10]*rdv2-1.5811388300841895*Ghat_l[10]*rdv2; 
  out[55] += 0.7071067811865475*Ghat_r[25]*rdv2-0.7071067811865475*Ghat_l[25]*rdv2; 
  out[56] += 0.7071067811865475*Ghat_r[26]*rdv2-0.7071067811865475*Ghat_l[26]*rdv2; 
  out[57] += 1.224744871391589*Ghat_r[19]*rdv2+1.224744871391589*Ghat_l[19]*rdv2; 
  out[58] += 0.7071067811865475*Ghat_r[27]*rdv2-0.7071067811865475*Ghat_l[27]*rdv2; 
  out[59] += 1.224744871391589*Ghat_r[20]*rdv2+1.224744871391589*Ghat_l[20]*rdv2; 
  out[60] += 1.224744871391589*Ghat_r[21]*rdv2+1.224744871391589*Ghat_l[21]*rdv2; 
  out[61] += 0.7071067811865475*Ghat_r[28]*rdv2-0.7071067811865475*Ghat_l[28]*rdv2; 
  out[62] += 0.7071067811865475*Ghat_r[29]*rdv2-0.7071067811865475*Ghat_l[29]*rdv2; 
  out[63] += 1.224744871391589*Ghat_r[22]*rdv2+1.224744871391589*Ghat_l[22]*rdv2; 
  out[64] += 1.5811388300841895*Ghat_r[11]*rdv2-1.5811388300841895*Ghat_l[11]*rdv2; 
  out[65] += 1.5811388300841895*Ghat_r[12]*rdv2-1.5811388300841895*Ghat_l[12]*rdv2; 
  out[66] += 0.7071067811865475*Ghat_r[30]*rdv2-0.7071067811865475*Ghat_l[30]*rdv2; 
  out[67] += 1.224744871391589*Ghat_r[23]*rdv2+1.224744871391589*Ghat_l[23]*rdv2; 
  out[68] += 1.5811388300841898*Ghat_r[13]*rdv2-1.5811388300841898*Ghat_l[13]*rdv2; 
  out[69] += 1.224744871391589*Ghat_r[24]*rdv2+1.224744871391589*Ghat_l[24]*rdv2; 
  out[70] += 1.5811388300841898*Ghat_r[14]*rdv2-1.5811388300841898*Ghat_l[14]*rdv2; 
  out[71] += 1.5811388300841898*Ghat_r[15]*rdv2-1.5811388300841898*Ghat_l[15]*rdv2; 
  out[72] += 1.5811388300841898*Ghat_r[16]*rdv2-1.5811388300841898*Ghat_l[16]*rdv2; 
  out[73] += 0.7071067811865475*Ghat_r[31]*rdv2-0.7071067811865475*Ghat_l[31]*rdv2; 
  out[74] += 1.224744871391589*Ghat_r[25]*rdv2+1.224744871391589*Ghat_l[25]*rdv2; 
  out[75] += 1.224744871391589*Ghat_r[26]*rdv2+1.224744871391589*Ghat_l[26]*rdv2; 
  out[76] += 1.224744871391589*Ghat_r[27]*rdv2+1.224744871391589*Ghat_l[27]*rdv2; 
  out[77] += 0.7071067811865475*Ghat_r[32]*rdv2-0.7071067811865475*Ghat_l[32]*rdv2; 
  out[78] += 1.224744871391589*Ghat_r[28]*rdv2+1.224744871391589*Ghat_l[28]*rdv2; 
  out[79] += 1.224744871391589*Ghat_r[29]*rdv2+1.224744871391589*Ghat_l[29]*rdv2; 
  out[80] += 1.5811388300841898*Ghat_r[17]*rdv2-1.5811388300841898*Ghat_l[17]*rdv2; 
  out[81] += 1.5811388300841898*Ghat_r[18]*rdv2-1.5811388300841898*Ghat_l[18]*rdv2; 
  out[82] += 1.5811388300841898*Ghat_r[19]*rdv2-1.5811388300841898*Ghat_l[19]*rdv2; 
  out[83] += 1.5811388300841898*Ghat_r[20]*rdv2-1.5811388300841898*Ghat_l[20]*rdv2; 
  out[84] += 1.5811388300841898*Ghat_r[21]*rdv2-1.5811388300841898*Ghat_l[21]*rdv2; 
  out[85] += 1.5811388300841898*Ghat_r[22]*rdv2-1.5811388300841898*Ghat_l[22]*rdv2; 
  out[86] += 0.7071067811865475*Ghat_r[33]*rdv2-0.7071067811865475*Ghat_l[33]*rdv2; 
  out[87] += 0.7071067811865475*Ghat_r[34]*rdv2-0.7071067811865475*Ghat_l[34]*rdv2; 
  out[88] += 1.224744871391589*Ghat_r[30]*rdv2+1.224744871391589*Ghat_l[30]*rdv2; 
  out[89] += 1.5811388300841895*Ghat_r[23]*rdv2-1.5811388300841895*Ghat_l[23]*rdv2; 
  out[90] += 1.224744871391589*Ghat_r[31]*rdv2+1.224744871391589*Ghat_l[31]*rdv2; 
  out[91] += 1.224744871391589*Ghat_r[32]*rdv2+1.224744871391589*Ghat_l[32]*rdv2; 
  out[92] += 1.5811388300841895*Ghat_r[24]*rdv2-1.5811388300841895*Ghat_l[24]*rdv2; 
  out[93] += 1.5811388300841895*Ghat_r[25]*rdv2-1.5811388300841895*Ghat_l[25]*rdv2; 
  out[94] += 1.5811388300841895*Ghat_r[26]*rdv2-1.5811388300841895*Ghat_l[26]*rdv2; 
  out[95] += 1.5811388300841895*Ghat_r[27]*rdv2-1.5811388300841895*Ghat_l[27]*rdv2; 
  out[96] += 1.5811388300841895*Ghat_r[28]*rdv2-1.5811388300841895*Ghat_l[28]*rdv2; 
  out[97] += 1.5811388300841895*Ghat_r[29]*rdv2-1.5811388300841895*Ghat_l[29]*rdv2; 
  out[98] += 0.7071067811865475*Ghat_r[35]*rdv2-0.7071067811865475*Ghat_l[35]*rdv2; 
  out[99] += 1.224744871391589*Ghat_r[33]*rdv2+1.224744871391589*Ghat_l[33]*rdv2; 
  out[100] += 1.224744871391589*Ghat_r[34]*rdv2+1.224744871391589*Ghat_l[34]*rdv2; 
  out[101] += 1.5811388300841895*Ghat_r[30]*rdv2-1.5811388300841895*Ghat_l[30]*rdv2; 
  out[102] += 1.5811388300841898*Ghat_r[31]*rdv2-1.5811388300841898*Ghat_l[31]*rdv2; 
  out[103] += 1.5811388300841898*Ghat_r[32]*rdv2-1.5811388300841898*Ghat_l[32]*rdv2; 
  out[104] += 1.224744871391589*Ghat_r[35]*rdv2+1.224744871391589*Ghat_l[35]*rdv2; 
  out[105] += 1.5811388300841898*Ghat_r[33]*rdv2-1.5811388300841898*Ghat_l[33]*rdv2; 
  out[106] += 1.5811388300841898*Ghat_r[34]*rdv2-1.5811388300841898*Ghat_l[34]*rdv2; 
  out[107] += 1.5811388300841895*Ghat_r[35]*rdv2-1.5811388300841895*Ghat_l[35]*rdv2; 

  return 0.;

} 
