#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p3_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[24]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf[12] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[12] = {0.0};;
  double Ghat[12] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.0*w[2]+dxv[2])-2.0*sumNuUx[1]); 
  alphaDrSurf[2] = 0.5*(2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]+dxv[2]*nuSum[2]); 
  alphaDrSurf[3] = -(0.5*(2.0*sumNuUx[3]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[3])); 
  alphaDrSurf[4] = -(0.5*(2.0*sumNuUx[4]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[4])); 
  alphaDrSurf[5] = -(0.5*(2.0*sumNuUx[5]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[5])); 
  alphaDrSurf[6] = -(0.5*(2.0*sumNuUx[6]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[6])); 
  alphaDrSurf[7] = -(0.5*(2.0*sumNuUx[7]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[7])); 
  alphaDrSurf[8] = -(0.5*(2.0*sumNuUx[8]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[8])); 
  alphaDrSurf[9] = -(0.5*(2.0*sumNuUx[9]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[9])); 
  alphaDrSurf[10] = -(0.5*(2.0*sumNuUx[10]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[10])); 
  alphaDrSurf[11] = -(0.5*(2.0*sumNuUx[11]+(-(2.0*w[2])-1.0*dxv[2])*nuSum[11])); 

  if (0.6012996787521189*alphaDrSurf[11]+0.6012996787521189*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]-1.02111731798518*alphaDrSurf[7]-1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+1.1123336207187142*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p3_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p3_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (-(0.8123851134177964*alphaDrSurf[11])+0.23739620491785088*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]+0.5446649474682886*alphaDrSurf[7]-0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x_p3_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_3x_p3_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.8123851134177964*alphaDrSurf[11]-0.23739620491785088*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]+0.5446649474682886*alphaDrSurf[7]+0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x_p3_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p3_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (-(0.6012996787521189*alphaDrSurf[11])-0.6012996787521189*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]-1.02111731798518*alphaDrSurf[7]+1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-1.1123336207187142*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p3_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p3_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (0.23739620491785063*alphaDrSurf[11]-0.8123851134177964*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]-0.40314236749410914*alphaDrSurf[7]+0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x_p3_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_3x_p3_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if (-(0.3207338198772712*alphaDrSurf[11])-0.3207338198772712*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]+0.2150365218040565*alphaDrSurf[7]+0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.17338066499557184*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x_p3_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_3x_p3_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (0.3207338198772712*alphaDrSurf[11]+0.3207338198772712*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]+0.2150365218040565*alphaDrSurf[7]-0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.17338066499557184*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p3_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_3x_p3_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if (-(0.23739620491785063*alphaDrSurf[11])+0.8123851134177964*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]-0.40314236749410914*alphaDrSurf[7]-0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x_p3_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_3x_p3_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (-(0.23739620491785063*alphaDrSurf[11])+0.8123851134177964*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]+0.40314236749410914*alphaDrSurf[7]+0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x_p3_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_3x_p3_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if (0.3207338198772712*alphaDrSurf[11]+0.3207338198772712*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]-0.2150365218040565*alphaDrSurf[7]+0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.17338066499557184*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_3x_p3_surfx3_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_3x_p3_surfx3_eval_quad_node_9_l(fEdge); 
  } 
  if (-(0.3207338198772712*alphaDrSurf[11])-0.3207338198772712*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]-0.2150365218040565*alphaDrSurf[7]-0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.17338066499557184*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_3x_p3_surfx3_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_3x_p3_surfx3_eval_quad_node_10_l(fEdge); 
  } 
  if (0.23739620491785063*alphaDrSurf[11]-0.8123851134177964*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]+0.40314236749410914*alphaDrSurf[7]-0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_3x_p3_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_3x_p3_surfx3_eval_quad_node_11_l(fEdge); 
  } 
  if (-(0.6012996787521189*alphaDrSurf[11])-0.6012996787521189*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]+1.02111731798518*alphaDrSurf[7]-1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-1.1123336207187142*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_3x_p3_surfx3_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_3x_p3_surfx3_eval_quad_node_12_l(fEdge); 
  } 
  if (0.8123851134177964*alphaDrSurf[11]-0.23739620491785088*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]-0.5446649474682886*alphaDrSurf[7]-0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_3x_p3_surfx3_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_3x_p3_surfx3_eval_quad_node_13_l(fEdge); 
  } 
  if (-(0.8123851134177964*alphaDrSurf[11])+0.23739620491785088*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]-0.5446649474682886*alphaDrSurf[7]+0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_3x_p3_surfx3_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_3x_p3_surfx3_eval_quad_node_14_l(fEdge); 
  } 
  if (0.6012996787521189*alphaDrSurf[11]+0.6012996787521189*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]+1.02111731798518*alphaDrSurf[7]+1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+1.1123336207187142*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_3x_p3_surfx3_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_3x_p3_surfx3_eval_quad_node_15_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[11]*fUpwind[11]+0.5*alphaDrSurf[10]*fUpwind[10]+0.5*alphaDrSurf[9]*fUpwind[9]+0.5*alphaDrSurf[8]*fUpwind[8]+0.5*alphaDrSurf[7]*fUpwind[7]+0.5*alphaDrSurf[6]*fUpwind[6]+0.5*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[3]*fUpwind[3]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alphaDrSurf[9]*fUpwind[11]+0.5*fUpwind[9]*alphaDrSurf[11]+0.4391550328268399*alphaDrSurf[6]*fUpwind[10]+0.4391550328268399*fUpwind[6]*alphaDrSurf[10]+0.4391550328268398*alphaDrSurf[4]*fUpwind[8]+0.4391550328268398*fUpwind[4]*alphaDrSurf[8]+0.5000000000000001*alphaDrSurf[5]*fUpwind[7]+0.5000000000000001*fUpwind[5]*alphaDrSurf[7]+0.44721359549995804*alphaDrSurf[3]*fUpwind[6]+0.44721359549995804*fUpwind[3]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.4391550328268399*alphaDrSurf[7]*fUpwind[11]+0.4391550328268399*fUpwind[7]*alphaDrSurf[11]+0.5*alphaDrSurf[8]*fUpwind[10]+0.5*fUpwind[8]*alphaDrSurf[10]+0.4391550328268398*alphaDrSurf[5]*fUpwind[9]+0.4391550328268398*fUpwind[5]*alphaDrSurf[9]+0.44721359549995804*alphaDrSurf[3]*fUpwind[7]+0.44721359549995804*fUpwind[3]*alphaDrSurf[7]+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[2]*fUpwind[5]+0.4472135954999579*fUpwind[2]*alphaDrSurf[5]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.4391550328268399*alphaDrSurf[5]*fUpwind[11]+0.4391550328268399*fUpwind[5]*alphaDrSurf[11]+0.4391550328268399*alphaDrSurf[4]*fUpwind[10]+0.4391550328268399*fUpwind[4]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[7]*fUpwind[9]+0.4391550328268399*fUpwind[7]*alphaDrSurf[9]+0.4391550328268399*alphaDrSurf[6]*fUpwind[8]+0.4391550328268399*fUpwind[6]*alphaDrSurf[8]+0.4*alphaDrSurf[6]*fUpwind[7]+0.44721359549995804*alphaDrSurf[2]*fUpwind[7]+0.4*fUpwind[6]*alphaDrSurf[7]+0.44721359549995804*fUpwind[2]*alphaDrSurf[7]+0.44721359549995804*alphaDrSurf[1]*fUpwind[6]+0.44721359549995804*fUpwind[1]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[3]*fUpwind[5]+0.4472135954999579*fUpwind[3]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[4] = 0.4472135954999579*alphaDrSurf[11]*fUpwind[11]+0.2981423969999719*alphaDrSurf[10]*fUpwind[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[10]+0.4391550328268399*fUpwind[3]*alphaDrSurf[10]+0.2981423969999719*alphaDrSurf[8]*fUpwind[8]+0.4391550328268398*alphaDrSurf[1]*fUpwind[8]+0.4391550328268398*fUpwind[1]*alphaDrSurf[8]+0.4472135954999579*alphaDrSurf[7]*fUpwind[7]+0.31943828249996997*alphaDrSurf[6]*fUpwind[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[6]+0.31943828249996997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.2981423969999719*alphaDrSurf[11]*fUpwind[11]+0.4391550328268399*alphaDrSurf[3]*fUpwind[11]+0.4391550328268399*fUpwind[3]*alphaDrSurf[11]+0.4472135954999579*alphaDrSurf[10]*fUpwind[10]+0.2981423969999719*alphaDrSurf[9]*fUpwind[9]+0.4391550328268398*alphaDrSurf[2]*fUpwind[9]+0.4391550328268398*fUpwind[2]*alphaDrSurf[9]+0.31943828249996997*alphaDrSurf[7]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5000000000000001*fUpwind[1]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[6]*fUpwind[6]+0.31943828249996997*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*fUpwind[0]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[6] = 0.3927922024247863*alphaDrSurf[7]*fUpwind[11]+0.3927922024247863*fUpwind[7]*alphaDrSurf[11]+0.29814239699997186*alphaDrSurf[8]*fUpwind[10]+0.3927922024247863*alphaDrSurf[7]*fUpwind[10]+0.4391550328268399*alphaDrSurf[1]*fUpwind[10]+0.29814239699997186*fUpwind[8]*alphaDrSurf[10]+0.3927922024247863*fUpwind[7]*alphaDrSurf[10]+0.4391550328268399*fUpwind[1]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[8]+0.4391550328268399*fUpwind[3]*alphaDrSurf[8]+0.4*alphaDrSurf[3]*fUpwind[7]+0.4*fUpwind[3]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[5]*fUpwind[6]+0.31943828249996997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alphaDrSurf[6]+0.31943828249996997*fUpwind[4]*alphaDrSurf[6]+0.5*fUpwind[0]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.44721359549995804*alphaDrSurf[1]*fUpwind[3]+0.44721359549995804*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.29814239699997186*alphaDrSurf[9]*fUpwind[11]+0.3927922024247863*alphaDrSurf[6]*fUpwind[11]+0.4391550328268399*alphaDrSurf[2]*fUpwind[11]+0.29814239699997186*fUpwind[9]*alphaDrSurf[11]+0.3927922024247863*fUpwind[6]*alphaDrSurf[11]+0.4391550328268399*fUpwind[2]*alphaDrSurf[11]+0.3927922024247863*alphaDrSurf[6]*fUpwind[10]+0.3927922024247863*fUpwind[6]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[9]+0.4391550328268399*fUpwind[3]*alphaDrSurf[9]+0.31943828249996997*alphaDrSurf[5]*fUpwind[7]+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.31943828249996997*fUpwind[5]*alphaDrSurf[7]+0.4472135954999579*fUpwind[4]*alphaDrSurf[7]+0.5*fUpwind[0]*alphaDrSurf[7]+0.4*alphaDrSurf[3]*fUpwind[6]+0.4*fUpwind[3]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]+0.5000000000000001*fUpwind[1]*alphaDrSurf[5]+0.44721359549995804*alphaDrSurf[2]*fUpwind[3]+0.44721359549995804*fUpwind[2]*alphaDrSurf[3]; 
  Ghat[8] = 0.29814239699997186*alphaDrSurf[6]*fUpwind[10]+0.5*alphaDrSurf[2]*fUpwind[10]+0.29814239699997186*fUpwind[6]*alphaDrSurf[10]+0.5*fUpwind[2]*alphaDrSurf[10]+0.2981423969999719*alphaDrSurf[4]*fUpwind[8]+0.5*alphaDrSurf[0]*fUpwind[8]+0.2981423969999719*fUpwind[4]*alphaDrSurf[8]+0.5*fUpwind[0]*alphaDrSurf[8]+0.4391550328268399*alphaDrSurf[3]*fUpwind[6]+0.4391550328268399*fUpwind[3]*alphaDrSurf[6]+0.4391550328268398*alphaDrSurf[1]*fUpwind[4]+0.4391550328268398*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[9] = 0.29814239699997186*alphaDrSurf[7]*fUpwind[11]+0.5*alphaDrSurf[1]*fUpwind[11]+0.29814239699997186*fUpwind[7]*alphaDrSurf[11]+0.5*fUpwind[1]*alphaDrSurf[11]+0.2981423969999719*alphaDrSurf[5]*fUpwind[9]+0.5*alphaDrSurf[0]*fUpwind[9]+0.2981423969999719*fUpwind[5]*alphaDrSurf[9]+0.5*fUpwind[0]*alphaDrSurf[9]+0.4391550328268399*alphaDrSurf[3]*fUpwind[7]+0.4391550328268399*fUpwind[3]*alphaDrSurf[7]+0.4391550328268398*alphaDrSurf[2]*fUpwind[5]+0.4391550328268398*fUpwind[2]*alphaDrSurf[5]; 
  Ghat[10] = 0.4472135954999579*alphaDrSurf[5]*fUpwind[10]+0.2981423969999719*alphaDrSurf[4]*fUpwind[10]+0.5*alphaDrSurf[0]*fUpwind[10]+0.4472135954999579*fUpwind[5]*alphaDrSurf[10]+0.2981423969999719*fUpwind[4]*alphaDrSurf[10]+0.5*fUpwind[0]*alphaDrSurf[10]+0.29814239699997186*alphaDrSurf[6]*fUpwind[8]+0.5*alphaDrSurf[2]*fUpwind[8]+0.29814239699997186*fUpwind[6]*alphaDrSurf[8]+0.5*fUpwind[2]*alphaDrSurf[8]+0.3927922024247863*alphaDrSurf[6]*fUpwind[7]+0.3927922024247863*fUpwind[6]*alphaDrSurf[7]+0.4391550328268399*alphaDrSurf[1]*fUpwind[6]+0.4391550328268399*fUpwind[1]*alphaDrSurf[6]+0.4391550328268399*alphaDrSurf[3]*fUpwind[4]+0.4391550328268399*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.2981423969999719*alphaDrSurf[5]*fUpwind[11]+0.4472135954999579*alphaDrSurf[4]*fUpwind[11]+0.5*alphaDrSurf[0]*fUpwind[11]+0.2981423969999719*fUpwind[5]*alphaDrSurf[11]+0.4472135954999579*fUpwind[4]*alphaDrSurf[11]+0.5*fUpwind[0]*alphaDrSurf[11]+0.29814239699997186*alphaDrSurf[7]*fUpwind[9]+0.5*alphaDrSurf[1]*fUpwind[9]+0.29814239699997186*fUpwind[7]*alphaDrSurf[9]+0.5*fUpwind[1]*alphaDrSurf[9]+0.3927922024247863*alphaDrSurf[6]*fUpwind[7]+0.4391550328268399*alphaDrSurf[2]*fUpwind[7]+0.3927922024247863*fUpwind[6]*alphaDrSurf[7]+0.4391550328268399*fUpwind[2]*alphaDrSurf[7]+0.4391550328268399*alphaDrSurf[3]*fUpwind[5]+0.4391550328268399*fUpwind[3]*alphaDrSurf[5]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += 1.5811388300841898*Ghat[1]*rdv2; 
  out[16] += 1.5811388300841898*Ghat[2]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[19] += 1.8708286933869707*Ghat[0]*rdv2; 
  out[20] += 1.224744871391589*Ghat[6]*rdv2; 
  out[21] += 1.224744871391589*Ghat[7]*rdv2; 
  out[22] += 1.5811388300841895*Ghat[3]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[25] += 1.2247448713915887*Ghat[8]*rdv2; 
  out[26] += 1.2247448713915887*Ghat[9]*rdv2; 
  out[27] += 1.8708286933869707*Ghat[1]*rdv2; 
  out[28] += 1.8708286933869707*Ghat[2]*rdv2; 
  out[29] += 1.2247448713915887*Ghat[10]*rdv2; 
  out[30] += 1.2247448713915887*Ghat[11]*rdv2; 
  out[31] += 1.8708286933869707*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*sumNuUx[1]); 
  alphaDrSurf[2] = 0.5*(2.0*nuSum[2]*w[2]-2.0*sumNuUx[2]-1.0*dxv[2]*nuSum[2]); 
  alphaDrSurf[3] = -(0.5*(2.0*sumNuUx[3]+(dxv[2]-2.0*w[2])*nuSum[3])); 
  alphaDrSurf[4] = -(0.5*(2.0*sumNuUx[4]+(dxv[2]-2.0*w[2])*nuSum[4])); 
  alphaDrSurf[5] = -(0.5*(2.0*sumNuUx[5]+(dxv[2]-2.0*w[2])*nuSum[5])); 
  alphaDrSurf[6] = -(0.5*(2.0*sumNuUx[6]+(dxv[2]-2.0*w[2])*nuSum[6])); 
  alphaDrSurf[7] = -(0.5*(2.0*sumNuUx[7]+(dxv[2]-2.0*w[2])*nuSum[7])); 
  alphaDrSurf[8] = -(0.5*(2.0*sumNuUx[8]+(dxv[2]-2.0*w[2])*nuSum[8])); 
  alphaDrSurf[9] = -(0.5*(2.0*sumNuUx[9]+(dxv[2]-2.0*w[2])*nuSum[9])); 
  alphaDrSurf[10] = -(0.5*(2.0*sumNuUx[10]+(dxv[2]-2.0*w[2])*nuSum[10])); 
  alphaDrSurf[11] = -(0.5*(2.0*sumNuUx[11]+(dxv[2]-2.0*w[2])*nuSum[11])); 

  if (0.6012996787521189*alphaDrSurf[11]+0.6012996787521189*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]-1.02111731798518*alphaDrSurf[7]-1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+1.1123336207187142*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p3_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p3_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (-(0.8123851134177964*alphaDrSurf[11])+0.23739620491785088*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]+0.5446649474682886*alphaDrSurf[7]-0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x_p3_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_3x_p3_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.8123851134177964*alphaDrSurf[11]-0.23739620491785088*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]+0.5446649474682886*alphaDrSurf[7]+0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x_p3_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p3_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (-(0.6012996787521189*alphaDrSurf[11])-0.6012996787521189*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]-0.40314236749410914*alphaDrSurf[8]-1.02111731798518*alphaDrSurf[7]+1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-1.1123336207187142*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]-0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p3_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p3_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (0.23739620491785063*alphaDrSurf[11]-0.8123851134177964*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]-0.40314236749410914*alphaDrSurf[7]+0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x_p3_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_3x_p3_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if (-(0.3207338198772712*alphaDrSurf[11])-0.3207338198772712*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]+0.2150365218040565*alphaDrSurf[7]+0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.17338066499557184*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x_p3_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_3x_p3_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (0.3207338198772712*alphaDrSurf[11]+0.3207338198772712*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]+0.2150365218040565*alphaDrSurf[7]-0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.17338066499557184*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p3_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x_p3_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if (-(0.23739620491785063*alphaDrSurf[11])+0.8123851134177964*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]+0.5446649474682886*alphaDrSurf[8]-0.40314236749410914*alphaDrSurf[7]-0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]-0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x_p3_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_3x_p3_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (-(0.23739620491785063*alphaDrSurf[11])+0.8123851134177964*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]+0.40314236749410914*alphaDrSurf[7]+0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x_p3_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_3x_p3_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if (0.3207338198772712*alphaDrSurf[11]+0.3207338198772712*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]-0.2150365218040565*alphaDrSurf[7]+0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]-0.17338066499557184*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_3x_p3_surfx3_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_3x_p3_surfx3_eval_quad_node_9_l(fSkin); 
  } 
  if (-(0.3207338198772712*alphaDrSurf[11])-0.3207338198772712*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]-0.2150365218040565*alphaDrSurf[7]-0.2150365218040565*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.17338066499557184*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_3x_p3_surfx3_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_3x_p3_surfx3_eval_quad_node_10_l(fSkin); 
  } 
  if (0.23739620491785063*alphaDrSurf[11]-0.8123851134177964*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]-0.5446649474682886*alphaDrSurf[8]+0.40314236749410914*alphaDrSurf[7]-0.5446649474682886*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]-0.3651715179178388*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]+0.2944322205496301*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_3x_p3_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_3x_p3_surfx3_eval_quad_node_11_l(fSkin); 
  } 
  if (-(0.6012996787521189*alphaDrSurf[11])-0.6012996787521189*alphaDrSurf[10]-0.40314236749410914*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]+1.02111731798518*alphaDrSurf[7]-1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-1.1123336207187142*alphaDrSurf[3]-0.7457659219616816*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_3x_p3_surfx3_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_3x_p3_surfx3_eval_quad_node_12_l(fSkin); 
  } 
  if (0.8123851134177964*alphaDrSurf[11]-0.23739620491785088*alphaDrSurf[10]+0.5446649474682886*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]-0.5446649474682886*alphaDrSurf[7]-0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]-0.4391550328268395*alphaDrSurf[3]-0.2944322205496301*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_3x_p3_surfx3_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_3x_p3_surfx3_eval_quad_node_13_l(fSkin); 
  } 
  if (-(0.8123851134177964*alphaDrSurf[11])+0.23739620491785088*alphaDrSurf[10]-0.5446649474682886*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]-0.5446649474682886*alphaDrSurf[7]+0.40314236749410914*alphaDrSurf[6]-0.3651715179178388*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+0.4391550328268395*alphaDrSurf[3]+0.2944322205496301*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_3x_p3_surfx3_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_3x_p3_surfx3_eval_quad_node_14_l(fSkin); 
  } 
  if (0.6012996787521189*alphaDrSurf[11]+0.6012996787521189*alphaDrSurf[10]+0.40314236749410914*alphaDrSurf[9]+0.40314236749410914*alphaDrSurf[8]+1.02111731798518*alphaDrSurf[7]+1.02111731798518*alphaDrSurf[6]+0.684609800417808*alphaDrSurf[5]+0.684609800417808*alphaDrSurf[4]+1.1123336207187142*alphaDrSurf[3]+0.7457659219616816*alphaDrSurf[2]+0.7457659219616816*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_3x_p3_surfx3_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_3x_p3_surfx3_eval_quad_node_15_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[11]*fUpwind[11]+0.5*alphaDrSurf[10]*fUpwind[10]+0.5*alphaDrSurf[9]*fUpwind[9]+0.5*alphaDrSurf[8]*fUpwind[8]+0.5*alphaDrSurf[7]*fUpwind[7]+0.5*alphaDrSurf[6]*fUpwind[6]+0.5*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[3]*fUpwind[3]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alphaDrSurf[9]*fUpwind[11]+0.5*fUpwind[9]*alphaDrSurf[11]+0.4391550328268399*alphaDrSurf[6]*fUpwind[10]+0.4391550328268399*fUpwind[6]*alphaDrSurf[10]+0.4391550328268398*alphaDrSurf[4]*fUpwind[8]+0.4391550328268398*fUpwind[4]*alphaDrSurf[8]+0.5000000000000001*alphaDrSurf[5]*fUpwind[7]+0.5000000000000001*fUpwind[5]*alphaDrSurf[7]+0.44721359549995804*alphaDrSurf[3]*fUpwind[6]+0.44721359549995804*fUpwind[3]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.4391550328268399*alphaDrSurf[7]*fUpwind[11]+0.4391550328268399*fUpwind[7]*alphaDrSurf[11]+0.5*alphaDrSurf[8]*fUpwind[10]+0.5*fUpwind[8]*alphaDrSurf[10]+0.4391550328268398*alphaDrSurf[5]*fUpwind[9]+0.4391550328268398*fUpwind[5]*alphaDrSurf[9]+0.44721359549995804*alphaDrSurf[3]*fUpwind[7]+0.44721359549995804*fUpwind[3]*alphaDrSurf[7]+0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[2]*fUpwind[5]+0.4472135954999579*fUpwind[2]*alphaDrSurf[5]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.4391550328268399*alphaDrSurf[5]*fUpwind[11]+0.4391550328268399*fUpwind[5]*alphaDrSurf[11]+0.4391550328268399*alphaDrSurf[4]*fUpwind[10]+0.4391550328268399*fUpwind[4]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[7]*fUpwind[9]+0.4391550328268399*fUpwind[7]*alphaDrSurf[9]+0.4391550328268399*alphaDrSurf[6]*fUpwind[8]+0.4391550328268399*fUpwind[6]*alphaDrSurf[8]+0.4*alphaDrSurf[6]*fUpwind[7]+0.44721359549995804*alphaDrSurf[2]*fUpwind[7]+0.4*fUpwind[6]*alphaDrSurf[7]+0.44721359549995804*fUpwind[2]*alphaDrSurf[7]+0.44721359549995804*alphaDrSurf[1]*fUpwind[6]+0.44721359549995804*fUpwind[1]*alphaDrSurf[6]+0.4472135954999579*alphaDrSurf[3]*fUpwind[5]+0.4472135954999579*fUpwind[3]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[4] = 0.4472135954999579*alphaDrSurf[11]*fUpwind[11]+0.2981423969999719*alphaDrSurf[10]*fUpwind[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[10]+0.4391550328268399*fUpwind[3]*alphaDrSurf[10]+0.2981423969999719*alphaDrSurf[8]*fUpwind[8]+0.4391550328268398*alphaDrSurf[1]*fUpwind[8]+0.4391550328268398*fUpwind[1]*alphaDrSurf[8]+0.4472135954999579*alphaDrSurf[7]*fUpwind[7]+0.31943828249996997*alphaDrSurf[6]*fUpwind[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[6]+0.31943828249996997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.2981423969999719*alphaDrSurf[11]*fUpwind[11]+0.4391550328268399*alphaDrSurf[3]*fUpwind[11]+0.4391550328268399*fUpwind[3]*alphaDrSurf[11]+0.4472135954999579*alphaDrSurf[10]*fUpwind[10]+0.2981423969999719*alphaDrSurf[9]*fUpwind[9]+0.4391550328268398*alphaDrSurf[2]*fUpwind[9]+0.4391550328268398*fUpwind[2]*alphaDrSurf[9]+0.31943828249996997*alphaDrSurf[7]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5000000000000001*fUpwind[1]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[6]*fUpwind[6]+0.31943828249996997*alphaDrSurf[5]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*fUpwind[0]*alphaDrSurf[5]+0.4472135954999579*alphaDrSurf[3]*fUpwind[3]+0.4472135954999579*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[6] = 0.3927922024247863*alphaDrSurf[7]*fUpwind[11]+0.3927922024247863*fUpwind[7]*alphaDrSurf[11]+0.29814239699997186*alphaDrSurf[8]*fUpwind[10]+0.3927922024247863*alphaDrSurf[7]*fUpwind[10]+0.4391550328268399*alphaDrSurf[1]*fUpwind[10]+0.29814239699997186*fUpwind[8]*alphaDrSurf[10]+0.3927922024247863*fUpwind[7]*alphaDrSurf[10]+0.4391550328268399*fUpwind[1]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[8]+0.4391550328268399*fUpwind[3]*alphaDrSurf[8]+0.4*alphaDrSurf[3]*fUpwind[7]+0.4*fUpwind[3]*alphaDrSurf[7]+0.4472135954999579*alphaDrSurf[5]*fUpwind[6]+0.31943828249996997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alphaDrSurf[6]+0.31943828249996997*fUpwind[4]*alphaDrSurf[6]+0.5*fUpwind[0]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.44721359549995804*alphaDrSurf[1]*fUpwind[3]+0.44721359549995804*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.29814239699997186*alphaDrSurf[9]*fUpwind[11]+0.3927922024247863*alphaDrSurf[6]*fUpwind[11]+0.4391550328268399*alphaDrSurf[2]*fUpwind[11]+0.29814239699997186*fUpwind[9]*alphaDrSurf[11]+0.3927922024247863*fUpwind[6]*alphaDrSurf[11]+0.4391550328268399*fUpwind[2]*alphaDrSurf[11]+0.3927922024247863*alphaDrSurf[6]*fUpwind[10]+0.3927922024247863*fUpwind[6]*alphaDrSurf[10]+0.4391550328268399*alphaDrSurf[3]*fUpwind[9]+0.4391550328268399*fUpwind[3]*alphaDrSurf[9]+0.31943828249996997*alphaDrSurf[5]*fUpwind[7]+0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.31943828249996997*fUpwind[5]*alphaDrSurf[7]+0.4472135954999579*fUpwind[4]*alphaDrSurf[7]+0.5*fUpwind[0]*alphaDrSurf[7]+0.4*alphaDrSurf[3]*fUpwind[6]+0.4*fUpwind[3]*alphaDrSurf[6]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]+0.5000000000000001*fUpwind[1]*alphaDrSurf[5]+0.44721359549995804*alphaDrSurf[2]*fUpwind[3]+0.44721359549995804*fUpwind[2]*alphaDrSurf[3]; 
  Ghat[8] = 0.29814239699997186*alphaDrSurf[6]*fUpwind[10]+0.5*alphaDrSurf[2]*fUpwind[10]+0.29814239699997186*fUpwind[6]*alphaDrSurf[10]+0.5*fUpwind[2]*alphaDrSurf[10]+0.2981423969999719*alphaDrSurf[4]*fUpwind[8]+0.5*alphaDrSurf[0]*fUpwind[8]+0.2981423969999719*fUpwind[4]*alphaDrSurf[8]+0.5*fUpwind[0]*alphaDrSurf[8]+0.4391550328268399*alphaDrSurf[3]*fUpwind[6]+0.4391550328268399*fUpwind[3]*alphaDrSurf[6]+0.4391550328268398*alphaDrSurf[1]*fUpwind[4]+0.4391550328268398*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[9] = 0.29814239699997186*alphaDrSurf[7]*fUpwind[11]+0.5*alphaDrSurf[1]*fUpwind[11]+0.29814239699997186*fUpwind[7]*alphaDrSurf[11]+0.5*fUpwind[1]*alphaDrSurf[11]+0.2981423969999719*alphaDrSurf[5]*fUpwind[9]+0.5*alphaDrSurf[0]*fUpwind[9]+0.2981423969999719*fUpwind[5]*alphaDrSurf[9]+0.5*fUpwind[0]*alphaDrSurf[9]+0.4391550328268399*alphaDrSurf[3]*fUpwind[7]+0.4391550328268399*fUpwind[3]*alphaDrSurf[7]+0.4391550328268398*alphaDrSurf[2]*fUpwind[5]+0.4391550328268398*fUpwind[2]*alphaDrSurf[5]; 
  Ghat[10] = 0.4472135954999579*alphaDrSurf[5]*fUpwind[10]+0.2981423969999719*alphaDrSurf[4]*fUpwind[10]+0.5*alphaDrSurf[0]*fUpwind[10]+0.4472135954999579*fUpwind[5]*alphaDrSurf[10]+0.2981423969999719*fUpwind[4]*alphaDrSurf[10]+0.5*fUpwind[0]*alphaDrSurf[10]+0.29814239699997186*alphaDrSurf[6]*fUpwind[8]+0.5*alphaDrSurf[2]*fUpwind[8]+0.29814239699997186*fUpwind[6]*alphaDrSurf[8]+0.5*fUpwind[2]*alphaDrSurf[8]+0.3927922024247863*alphaDrSurf[6]*fUpwind[7]+0.3927922024247863*fUpwind[6]*alphaDrSurf[7]+0.4391550328268399*alphaDrSurf[1]*fUpwind[6]+0.4391550328268399*fUpwind[1]*alphaDrSurf[6]+0.4391550328268399*alphaDrSurf[3]*fUpwind[4]+0.4391550328268399*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.2981423969999719*alphaDrSurf[5]*fUpwind[11]+0.4472135954999579*alphaDrSurf[4]*fUpwind[11]+0.5*alphaDrSurf[0]*fUpwind[11]+0.2981423969999719*fUpwind[5]*alphaDrSurf[11]+0.4472135954999579*fUpwind[4]*alphaDrSurf[11]+0.5*fUpwind[0]*alphaDrSurf[11]+0.29814239699997186*alphaDrSurf[7]*fUpwind[9]+0.5*alphaDrSurf[1]*fUpwind[9]+0.29814239699997186*fUpwind[7]*alphaDrSurf[9]+0.5*fUpwind[1]*alphaDrSurf[9]+0.3927922024247863*alphaDrSurf[6]*fUpwind[7]+0.4391550328268399*alphaDrSurf[2]*fUpwind[7]+0.3927922024247863*fUpwind[6]*alphaDrSurf[7]+0.4391550328268399*fUpwind[2]*alphaDrSurf[7]+0.4391550328268399*alphaDrSurf[3]*fUpwind[5]+0.4391550328268399*fUpwind[3]*alphaDrSurf[5]; 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[8] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[9] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -(0.7071067811865475*Ghat[6]*rdv2); 
  out[12] += -(0.7071067811865475*Ghat[7]*rdv2); 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += -(1.5811388300841898*Ghat[1]*rdv2); 
  out[16] += -(1.5811388300841898*Ghat[2]*rdv2); 
  out[17] += -(0.7071067811865475*Ghat[8]*rdv2); 
  out[18] += -(0.7071067811865475*Ghat[9]*rdv2); 
  out[19] += 1.8708286933869707*Ghat[0]*rdv2; 
  out[20] += 1.224744871391589*Ghat[6]*rdv2; 
  out[21] += 1.224744871391589*Ghat[7]*rdv2; 
  out[22] += -(1.5811388300841895*Ghat[3]*rdv2); 
  out[23] += -(0.7071067811865475*Ghat[10]*rdv2); 
  out[24] += -(0.7071067811865475*Ghat[11]*rdv2); 
  out[25] += 1.2247448713915887*Ghat[8]*rdv2; 
  out[26] += 1.2247448713915887*Ghat[9]*rdv2; 
  out[27] += 1.8708286933869707*Ghat[1]*rdv2; 
  out[28] += 1.8708286933869707*Ghat[2]*rdv2; 
  out[29] += 1.2247448713915887*Ghat[10]*rdv2; 
  out[30] += 1.2247448713915887*Ghat[11]*rdv2; 
  out[31] += 1.8708286933869707*Ghat[3]*rdv2; 

  } 

  return 0.;

} 
