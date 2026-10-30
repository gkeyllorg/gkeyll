#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p3_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvx_2x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[24]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUx = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[12] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*w[2]-0.5*nuSum[0]*dxv[2]-1.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = nuSum[1]*w[2]-0.5*nuSum[1]*dxv[2]-1.0*sumNuUx[1]; 
  alphaDrSurf_l[2] = nuSum[2]*w[2]-1.0*sumNuUx[2]-0.5*dxv[2]*nuSum[2]; 
  alphaDrSurf_l[3] = -(1.0*sumNuUx[3])+w[2]*nuSum[3]-0.5*dxv[2]*nuSum[3]; 
  alphaDrSurf_l[4] = -(1.0*sumNuUx[4])+w[2]*nuSum[4]-0.5*dxv[2]*nuSum[4]; 
  alphaDrSurf_l[5] = -(1.0*sumNuUx[5])+w[2]*nuSum[5]-0.5*dxv[2]*nuSum[5]; 
  alphaDrSurf_l[6] = -(1.0*sumNuUx[6])+w[2]*nuSum[6]-0.5*dxv[2]*nuSum[6]; 
  alphaDrSurf_l[7] = -(1.0*sumNuUx[7])+w[2]*nuSum[7]-0.5*dxv[2]*nuSum[7]; 
  alphaDrSurf_l[8] = -(1.0*sumNuUx[8])+w[2]*nuSum[8]-0.5*dxv[2]*nuSum[8]; 
  alphaDrSurf_l[9] = -(1.0*sumNuUx[9])+w[2]*nuSum[9]-0.5*dxv[2]*nuSum[9]; 
  alphaDrSurf_l[10] = -(1.0*sumNuUx[10])+w[2]*nuSum[10]-0.5*dxv[2]*nuSum[10]; 
  alphaDrSurf_l[11] = -(1.0*sumNuUx[11])+w[2]*nuSum[11]-0.5*dxv[2]*nuSum[11]; 

  double alphaDrSurf_r[12] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*w[2]+0.5*nuSum[0]*dxv[2]-1.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = nuSum[1]*w[2]+0.5*nuSum[1]*dxv[2]-1.0*sumNuUx[1]; 
  alphaDrSurf_r[2] = nuSum[2]*w[2]-1.0*sumNuUx[2]+0.5*dxv[2]*nuSum[2]; 
  alphaDrSurf_r[3] = -(1.0*sumNuUx[3])+w[2]*nuSum[3]+0.5*dxv[2]*nuSum[3]; 
  alphaDrSurf_r[4] = -(1.0*sumNuUx[4])+w[2]*nuSum[4]+0.5*dxv[2]*nuSum[4]; 
  alphaDrSurf_r[5] = -(1.0*sumNuUx[5])+w[2]*nuSum[5]+0.5*dxv[2]*nuSum[5]; 
  alphaDrSurf_r[6] = -(1.0*sumNuUx[6])+w[2]*nuSum[6]+0.5*dxv[2]*nuSum[6]; 
  alphaDrSurf_r[7] = -(1.0*sumNuUx[7])+w[2]*nuSum[7]+0.5*dxv[2]*nuSum[7]; 
  alphaDrSurf_r[8] = -(1.0*sumNuUx[8])+w[2]*nuSum[8]+0.5*dxv[2]*nuSum[8]; 
  alphaDrSurf_r[9] = -(1.0*sumNuUx[9])+w[2]*nuSum[9]+0.5*dxv[2]*nuSum[9]; 
  alphaDrSurf_r[10] = -(1.0*sumNuUx[10])+w[2]*nuSum[10]+0.5*dxv[2]*nuSum[10]; 
  alphaDrSurf_r[11] = -(1.0*sumNuUx[11])+w[2]*nuSum[11]+0.5*dxv[2]*nuSum[11]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[12] = {0.0};
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if (0.6012996787521189*alphaDrSurf_l[11]+0.6012996787521189*alphaDrSurf_l[10]-0.40314236749410914*alphaDrSurf_l[9]-0.40314236749410914*alphaDrSurf_l[8]-1.02111731798518*alphaDrSurf_l[7]-1.02111731798518*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]+1.1123336207187142*alphaDrSurf_l[3]-0.7457659219616816*alphaDrSurf_l[2]-0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_3x_p3_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p3_surfx3_eval_quad_node_0_l(fc); 
  } 
  if (0.6012996787521189*alphaDrSurf_r[11]+0.6012996787521189*alphaDrSurf_r[10]-0.40314236749410914*alphaDrSurf_r[9]-0.40314236749410914*alphaDrSurf_r[8]-1.02111731798518*alphaDrSurf_r[7]-1.02111731798518*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]+1.1123336207187142*alphaDrSurf_r[3]-0.7457659219616816*alphaDrSurf_r[2]-0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_3x_p3_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p3_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (-(0.8123851134177964*alphaDrSurf_l[11])+0.23739620491785088*alphaDrSurf_l[10]+0.5446649474682886*alphaDrSurf_l[9]-0.40314236749410914*alphaDrSurf_l[8]+0.5446649474682886*alphaDrSurf_l[7]-0.40314236749410914*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]+0.4391550328268395*alphaDrSurf_l[3]-0.2944322205496301*alphaDrSurf_l[2]-0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_3x_p3_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p3_surfx3_eval_quad_node_1_l(fc); 
  } 
  if (-(0.8123851134177964*alphaDrSurf_r[11])+0.23739620491785088*alphaDrSurf_r[10]+0.5446649474682886*alphaDrSurf_r[9]-0.40314236749410914*alphaDrSurf_r[8]+0.5446649474682886*alphaDrSurf_r[7]-0.40314236749410914*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]+0.4391550328268395*alphaDrSurf_r[3]-0.2944322205496301*alphaDrSurf_r[2]-0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_3x_p3_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p3_surfx3_eval_quad_node_1_l(fr); 
  } 
  if (0.8123851134177964*alphaDrSurf_l[11]-0.23739620491785088*alphaDrSurf_l[10]-0.5446649474682886*alphaDrSurf_l[9]-0.40314236749410914*alphaDrSurf_l[8]+0.5446649474682886*alphaDrSurf_l[7]+0.40314236749410914*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]-0.4391550328268395*alphaDrSurf_l[3]+0.2944322205496301*alphaDrSurf_l[2]-0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_3x_p3_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p3_surfx3_eval_quad_node_2_l(fc); 
  } 
  if (0.8123851134177964*alphaDrSurf_r[11]-0.23739620491785088*alphaDrSurf_r[10]-0.5446649474682886*alphaDrSurf_r[9]-0.40314236749410914*alphaDrSurf_r[8]+0.5446649474682886*alphaDrSurf_r[7]+0.40314236749410914*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]-0.4391550328268395*alphaDrSurf_r[3]+0.2944322205496301*alphaDrSurf_r[2]-0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_3x_p3_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p3_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (-(0.6012996787521189*alphaDrSurf_l[11])-0.6012996787521189*alphaDrSurf_l[10]+0.40314236749410914*alphaDrSurf_l[9]-0.40314236749410914*alphaDrSurf_l[8]-1.02111731798518*alphaDrSurf_l[7]+1.02111731798518*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]-1.1123336207187142*alphaDrSurf_l[3]+0.7457659219616816*alphaDrSurf_l[2]-0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_3x_p3_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p3_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (-(0.6012996787521189*alphaDrSurf_r[11])-0.6012996787521189*alphaDrSurf_r[10]+0.40314236749410914*alphaDrSurf_r[9]-0.40314236749410914*alphaDrSurf_r[8]-1.02111731798518*alphaDrSurf_r[7]+1.02111731798518*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]-1.1123336207187142*alphaDrSurf_r[3]+0.7457659219616816*alphaDrSurf_r[2]-0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_3x_p3_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x_p3_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.23739620491785063*alphaDrSurf_l[11]-0.8123851134177964*alphaDrSurf_l[10]-0.40314236749410914*alphaDrSurf_l[9]+0.5446649474682886*alphaDrSurf_l[8]-0.40314236749410914*alphaDrSurf_l[7]+0.5446649474682886*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]+0.4391550328268395*alphaDrSurf_l[3]-0.7457659219616816*alphaDrSurf_l[2]-0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_3x_p3_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_3x_p3_surfx3_eval_quad_node_4_l(fc); 
  } 
  if (0.23739620491785063*alphaDrSurf_r[11]-0.8123851134177964*alphaDrSurf_r[10]-0.40314236749410914*alphaDrSurf_r[9]+0.5446649474682886*alphaDrSurf_r[8]-0.40314236749410914*alphaDrSurf_r[7]+0.5446649474682886*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]+0.4391550328268395*alphaDrSurf_r[3]-0.7457659219616816*alphaDrSurf_r[2]-0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_3x_p3_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_3x_p3_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (-(0.3207338198772712*alphaDrSurf_l[11])-0.3207338198772712*alphaDrSurf_l[10]+0.5446649474682886*alphaDrSurf_l[9]+0.5446649474682886*alphaDrSurf_l[8]+0.2150365218040565*alphaDrSurf_l[7]+0.2150365218040565*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]+0.17338066499557184*alphaDrSurf_l[3]-0.2944322205496301*alphaDrSurf_l[2]-0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_3x_p3_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_3x_p3_surfx3_eval_quad_node_5_l(fc); 
  } 
  if (-(0.3207338198772712*alphaDrSurf_r[11])-0.3207338198772712*alphaDrSurf_r[10]+0.5446649474682886*alphaDrSurf_r[9]+0.5446649474682886*alphaDrSurf_r[8]+0.2150365218040565*alphaDrSurf_r[7]+0.2150365218040565*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]+0.17338066499557184*alphaDrSurf_r[3]-0.2944322205496301*alphaDrSurf_r[2]-0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_3x_p3_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_3x_p3_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.3207338198772712*alphaDrSurf_l[11]+0.3207338198772712*alphaDrSurf_l[10]-0.5446649474682886*alphaDrSurf_l[9]+0.5446649474682886*alphaDrSurf_l[8]+0.2150365218040565*alphaDrSurf_l[7]-0.2150365218040565*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]-0.17338066499557184*alphaDrSurf_l[3]+0.2944322205496301*alphaDrSurf_l[2]-0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_3x_p3_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_3x_p3_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (0.3207338198772712*alphaDrSurf_r[11]+0.3207338198772712*alphaDrSurf_r[10]-0.5446649474682886*alphaDrSurf_r[9]+0.5446649474682886*alphaDrSurf_r[8]+0.2150365218040565*alphaDrSurf_r[7]-0.2150365218040565*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]-0.17338066499557184*alphaDrSurf_r[3]+0.2944322205496301*alphaDrSurf_r[2]-0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_3x_p3_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_3x_p3_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (-(0.23739620491785063*alphaDrSurf_l[11])+0.8123851134177964*alphaDrSurf_l[10]+0.40314236749410914*alphaDrSurf_l[9]+0.5446649474682886*alphaDrSurf_l[8]-0.40314236749410914*alphaDrSurf_l[7]-0.5446649474682886*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]-0.4391550328268395*alphaDrSurf_l[3]+0.7457659219616816*alphaDrSurf_l[2]-0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_3x_p3_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_3x_p3_surfx3_eval_quad_node_7_l(fc); 
  } 
  if (-(0.23739620491785063*alphaDrSurf_r[11])+0.8123851134177964*alphaDrSurf_r[10]+0.40314236749410914*alphaDrSurf_r[9]+0.5446649474682886*alphaDrSurf_r[8]-0.40314236749410914*alphaDrSurf_r[7]-0.5446649474682886*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]-0.4391550328268395*alphaDrSurf_r[3]+0.7457659219616816*alphaDrSurf_r[2]-0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_3x_p3_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_3x_p3_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (-(0.23739620491785063*alphaDrSurf_l[11])+0.8123851134177964*alphaDrSurf_l[10]-0.40314236749410914*alphaDrSurf_l[9]-0.5446649474682886*alphaDrSurf_l[8]+0.40314236749410914*alphaDrSurf_l[7]+0.5446649474682886*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]-0.4391550328268395*alphaDrSurf_l[3]-0.7457659219616816*alphaDrSurf_l[2]+0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_3x_p3_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_3x_p3_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (-(0.23739620491785063*alphaDrSurf_r[11])+0.8123851134177964*alphaDrSurf_r[10]-0.40314236749410914*alphaDrSurf_r[9]-0.5446649474682886*alphaDrSurf_r[8]+0.40314236749410914*alphaDrSurf_r[7]+0.5446649474682886*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]-0.4391550328268395*alphaDrSurf_r[3]-0.7457659219616816*alphaDrSurf_r[2]+0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_3x_p3_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_3x_p3_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (0.3207338198772712*alphaDrSurf_l[11]+0.3207338198772712*alphaDrSurf_l[10]+0.5446649474682886*alphaDrSurf_l[9]-0.5446649474682886*alphaDrSurf_l[8]-0.2150365218040565*alphaDrSurf_l[7]+0.2150365218040565*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]-0.17338066499557184*alphaDrSurf_l[3]-0.2944322205496301*alphaDrSurf_l[2]+0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = ser_3x_p3_surfx3_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = ser_3x_p3_surfx3_eval_quad_node_9_l(fc); 
  } 
  if (0.3207338198772712*alphaDrSurf_r[11]+0.3207338198772712*alphaDrSurf_r[10]+0.5446649474682886*alphaDrSurf_r[9]-0.5446649474682886*alphaDrSurf_r[8]-0.2150365218040565*alphaDrSurf_r[7]+0.2150365218040565*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]-0.17338066499557184*alphaDrSurf_r[3]-0.2944322205496301*alphaDrSurf_r[2]+0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = ser_3x_p3_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = ser_3x_p3_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (-(0.3207338198772712*alphaDrSurf_l[11])-0.3207338198772712*alphaDrSurf_l[10]-0.5446649474682886*alphaDrSurf_l[9]-0.5446649474682886*alphaDrSurf_l[8]-0.2150365218040565*alphaDrSurf_l[7]-0.2150365218040565*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]+0.17338066499557184*alphaDrSurf_l[3]+0.2944322205496301*alphaDrSurf_l[2]+0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[10] = ser_3x_p3_surfx3_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = ser_3x_p3_surfx3_eval_quad_node_10_l(fc); 
  } 
  if (-(0.3207338198772712*alphaDrSurf_r[11])-0.3207338198772712*alphaDrSurf_r[10]-0.5446649474682886*alphaDrSurf_r[9]-0.5446649474682886*alphaDrSurf_r[8]-0.2150365218040565*alphaDrSurf_r[7]-0.2150365218040565*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]+0.17338066499557184*alphaDrSurf_r[3]+0.2944322205496301*alphaDrSurf_r[2]+0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[10] = ser_3x_p3_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = ser_3x_p3_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (0.23739620491785063*alphaDrSurf_l[11]-0.8123851134177964*alphaDrSurf_l[10]+0.40314236749410914*alphaDrSurf_l[9]-0.5446649474682886*alphaDrSurf_l[8]+0.40314236749410914*alphaDrSurf_l[7]-0.5446649474682886*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]-0.3651715179178388*alphaDrSurf_l[4]+0.4391550328268395*alphaDrSurf_l[3]+0.7457659219616816*alphaDrSurf_l[2]+0.2944322205496301*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[11] = ser_3x_p3_surfx3_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = ser_3x_p3_surfx3_eval_quad_node_11_l(fc); 
  } 
  if (0.23739620491785063*alphaDrSurf_r[11]-0.8123851134177964*alphaDrSurf_r[10]+0.40314236749410914*alphaDrSurf_r[9]-0.5446649474682886*alphaDrSurf_r[8]+0.40314236749410914*alphaDrSurf_r[7]-0.5446649474682886*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]-0.3651715179178388*alphaDrSurf_r[4]+0.4391550328268395*alphaDrSurf_r[3]+0.7457659219616816*alphaDrSurf_r[2]+0.2944322205496301*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[11] = ser_3x_p3_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = ser_3x_p3_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (-(0.6012996787521189*alphaDrSurf_l[11])-0.6012996787521189*alphaDrSurf_l[10]-0.40314236749410914*alphaDrSurf_l[9]+0.40314236749410914*alphaDrSurf_l[8]+1.02111731798518*alphaDrSurf_l[7]-1.02111731798518*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]-1.1123336207187142*alphaDrSurf_l[3]-0.7457659219616816*alphaDrSurf_l[2]+0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = ser_3x_p3_surfx3_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_3x_p3_surfx3_eval_quad_node_12_l(fc); 
  } 
  if (-(0.6012996787521189*alphaDrSurf_r[11])-0.6012996787521189*alphaDrSurf_r[10]-0.40314236749410914*alphaDrSurf_r[9]+0.40314236749410914*alphaDrSurf_r[8]+1.02111731798518*alphaDrSurf_r[7]-1.02111731798518*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]-1.1123336207187142*alphaDrSurf_r[3]-0.7457659219616816*alphaDrSurf_r[2]+0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = ser_3x_p3_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_3x_p3_surfx3_eval_quad_node_12_l(fr); 
  } 
  if (0.8123851134177964*alphaDrSurf_l[11]-0.23739620491785088*alphaDrSurf_l[10]+0.5446649474682886*alphaDrSurf_l[9]+0.40314236749410914*alphaDrSurf_l[8]-0.5446649474682886*alphaDrSurf_l[7]-0.40314236749410914*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]-0.4391550328268395*alphaDrSurf_l[3]-0.2944322205496301*alphaDrSurf_l[2]+0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[13] = ser_3x_p3_surfx3_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = ser_3x_p3_surfx3_eval_quad_node_13_l(fc); 
  } 
  if (0.8123851134177964*alphaDrSurf_r[11]-0.23739620491785088*alphaDrSurf_r[10]+0.5446649474682886*alphaDrSurf_r[9]+0.40314236749410914*alphaDrSurf_r[8]-0.5446649474682886*alphaDrSurf_r[7]-0.40314236749410914*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]-0.4391550328268395*alphaDrSurf_r[3]-0.2944322205496301*alphaDrSurf_r[2]+0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[13] = ser_3x_p3_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = ser_3x_p3_surfx3_eval_quad_node_13_l(fr); 
  } 
  if (-(0.8123851134177964*alphaDrSurf_l[11])+0.23739620491785088*alphaDrSurf_l[10]-0.5446649474682886*alphaDrSurf_l[9]+0.40314236749410914*alphaDrSurf_l[8]-0.5446649474682886*alphaDrSurf_l[7]+0.40314236749410914*alphaDrSurf_l[6]-0.3651715179178388*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]+0.4391550328268395*alphaDrSurf_l[3]+0.2944322205496301*alphaDrSurf_l[2]+0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[14] = ser_3x_p3_surfx3_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = ser_3x_p3_surfx3_eval_quad_node_14_l(fc); 
  } 
  if (-(0.8123851134177964*alphaDrSurf_r[11])+0.23739620491785088*alphaDrSurf_r[10]-0.5446649474682886*alphaDrSurf_r[9]+0.40314236749410914*alphaDrSurf_r[8]-0.5446649474682886*alphaDrSurf_r[7]+0.40314236749410914*alphaDrSurf_r[6]-0.3651715179178388*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]+0.4391550328268395*alphaDrSurf_r[3]+0.2944322205496301*alphaDrSurf_r[2]+0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[14] = ser_3x_p3_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = ser_3x_p3_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (0.6012996787521189*alphaDrSurf_l[11]+0.6012996787521189*alphaDrSurf_l[10]+0.40314236749410914*alphaDrSurf_l[9]+0.40314236749410914*alphaDrSurf_l[8]+1.02111731798518*alphaDrSurf_l[7]+1.02111731798518*alphaDrSurf_l[6]+0.684609800417808*alphaDrSurf_l[5]+0.684609800417808*alphaDrSurf_l[4]+1.1123336207187142*alphaDrSurf_l[3]+0.7457659219616816*alphaDrSurf_l[2]+0.7457659219616816*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[15] = ser_3x_p3_surfx3_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = ser_3x_p3_surfx3_eval_quad_node_15_l(fc); 
  } 
  if (0.6012996787521189*alphaDrSurf_r[11]+0.6012996787521189*alphaDrSurf_r[10]+0.40314236749410914*alphaDrSurf_r[9]+0.40314236749410914*alphaDrSurf_r[8]+1.02111731798518*alphaDrSurf_r[7]+1.02111731798518*alphaDrSurf_r[6]+0.684609800417808*alphaDrSurf_r[5]+0.684609800417808*alphaDrSurf_r[4]+1.1123336207187142*alphaDrSurf_r[3]+0.7457659219616816*alphaDrSurf_r[2]+0.7457659219616816*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[15] = ser_3x_p3_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = ser_3x_p3_surfx3_eval_quad_node_15_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p3_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p3_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaDrSurf_l[11]*fUpwind_l[11]+0.5*alphaDrSurf_l[10]*fUpwind_l[10]+0.5*alphaDrSurf_l[9]*fUpwind_l[9]+0.5*alphaDrSurf_l[8]*fUpwind_l[8]+0.5*alphaDrSurf_l[7]*fUpwind_l[7]+0.5*alphaDrSurf_l[6]*fUpwind_l[6]+0.5*alphaDrSurf_l[5]*fUpwind_l[5]+0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[3]*fUpwind_l[3]+0.5*alphaDrSurf_l[2]*fUpwind_l[2]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5*alphaDrSurf_l[9]*fUpwind_l[11]+0.5*fUpwind_l[9]*alphaDrSurf_l[11]+0.4391550328268399*alphaDrSurf_l[6]*fUpwind_l[10]+0.4391550328268399*fUpwind_l[6]*alphaDrSurf_l[10]+0.4391550328268398*alphaDrSurf_l[4]*fUpwind_l[8]+0.4391550328268398*fUpwind_l[4]*alphaDrSurf_l[8]+0.5000000000000001*alphaDrSurf_l[5]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[5]*alphaDrSurf_l[7]+0.44721359549995804*alphaDrSurf_l[3]*fUpwind_l[6]+0.44721359549995804*fUpwind_l[3]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.4391550328268399*alphaDrSurf_l[7]*fUpwind_l[11]+0.4391550328268399*fUpwind_l[7]*alphaDrSurf_l[11]+0.5*alphaDrSurf_l[8]*fUpwind_l[10]+0.5*fUpwind_l[8]*alphaDrSurf_l[10]+0.4391550328268398*alphaDrSurf_l[5]*fUpwind_l[9]+0.4391550328268398*fUpwind_l[5]*alphaDrSurf_l[9]+0.44721359549995804*alphaDrSurf_l[3]*fUpwind_l[7]+0.44721359549995804*fUpwind_l[3]*alphaDrSurf_l[7]+0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[2]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[2]*alphaDrSurf_l[5]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = 0.4391550328268399*alphaDrSurf_l[5]*fUpwind_l[11]+0.4391550328268399*fUpwind_l[5]*alphaDrSurf_l[11]+0.4391550328268399*alphaDrSurf_l[4]*fUpwind_l[10]+0.4391550328268399*fUpwind_l[4]*alphaDrSurf_l[10]+0.4391550328268399*alphaDrSurf_l[7]*fUpwind_l[9]+0.4391550328268399*fUpwind_l[7]*alphaDrSurf_l[9]+0.4391550328268399*alphaDrSurf_l[6]*fUpwind_l[8]+0.4391550328268399*fUpwind_l[6]*alphaDrSurf_l[8]+0.4*alphaDrSurf_l[6]*fUpwind_l[7]+0.44721359549995804*alphaDrSurf_l[2]*fUpwind_l[7]+0.4*fUpwind_l[6]*alphaDrSurf_l[7]+0.44721359549995804*fUpwind_l[2]*alphaDrSurf_l[7]+0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[6]+0.44721359549995804*fUpwind_l[1]*alphaDrSurf_l[6]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[5]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alphaDrSurf_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Ghat_l[4] = 0.4472135954999579*alphaDrSurf_l[11]*fUpwind_l[11]+0.2981423969999719*alphaDrSurf_l[10]*fUpwind_l[10]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[10]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[10]+0.2981423969999719*alphaDrSurf_l[8]*fUpwind_l[8]+0.4391550328268398*alphaDrSurf_l[1]*fUpwind_l[8]+0.4391550328268398*fUpwind_l[1]*alphaDrSurf_l[8]+0.4472135954999579*alphaDrSurf_l[7]*fUpwind_l[7]+0.31943828249996997*alphaDrSurf_l[6]*fUpwind_l[6]+0.5000000000000001*alphaDrSurf_l[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[6]+0.31943828249996997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.2981423969999719*alphaDrSurf_l[11]*fUpwind_l[11]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[11]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[11]+0.4472135954999579*alphaDrSurf_l[10]*fUpwind_l[10]+0.2981423969999719*alphaDrSurf_l[9]*fUpwind_l[9]+0.4391550328268398*alphaDrSurf_l[2]*fUpwind_l[9]+0.4391550328268398*fUpwind_l[2]*alphaDrSurf_l[9]+0.31943828249996997*alphaDrSurf_l[7]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[1]*alphaDrSurf_l[7]+0.4472135954999579*alphaDrSurf_l[6]*fUpwind_l[6]+0.31943828249996997*alphaDrSurf_l[5]*fUpwind_l[5]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]+0.5*fUpwind_l[0]*alphaDrSurf_l[5]+0.4472135954999579*alphaDrSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaDrSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.3927922024247863*alphaDrSurf_l[7]*fUpwind_l[11]+0.3927922024247863*fUpwind_l[7]*alphaDrSurf_l[11]+0.29814239699997186*alphaDrSurf_l[8]*fUpwind_l[10]+0.3927922024247863*alphaDrSurf_l[7]*fUpwind_l[10]+0.4391550328268399*alphaDrSurf_l[1]*fUpwind_l[10]+0.29814239699997186*fUpwind_l[8]*alphaDrSurf_l[10]+0.3927922024247863*fUpwind_l[7]*alphaDrSurf_l[10]+0.4391550328268399*fUpwind_l[1]*alphaDrSurf_l[10]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[8]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[8]+0.4*alphaDrSurf_l[3]*fUpwind_l[7]+0.4*fUpwind_l[3]*alphaDrSurf_l[7]+0.4472135954999579*alphaDrSurf_l[5]*fUpwind_l[6]+0.31943828249996997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alphaDrSurf_l[6]+0.31943828249996997*fUpwind_l[4]*alphaDrSurf_l[6]+0.5*fUpwind_l[0]*alphaDrSurf_l[6]+0.5000000000000001*alphaDrSurf_l[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.44721359549995804*alphaDrSurf_l[1]*fUpwind_l[3]+0.44721359549995804*fUpwind_l[1]*alphaDrSurf_l[3]; 
  Ghat_l[7] = 0.29814239699997186*alphaDrSurf_l[9]*fUpwind_l[11]+0.3927922024247863*alphaDrSurf_l[6]*fUpwind_l[11]+0.4391550328268399*alphaDrSurf_l[2]*fUpwind_l[11]+0.29814239699997186*fUpwind_l[9]*alphaDrSurf_l[11]+0.3927922024247863*fUpwind_l[6]*alphaDrSurf_l[11]+0.4391550328268399*fUpwind_l[2]*alphaDrSurf_l[11]+0.3927922024247863*alphaDrSurf_l[6]*fUpwind_l[10]+0.3927922024247863*fUpwind_l[6]*alphaDrSurf_l[10]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[9]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[9]+0.31943828249996997*alphaDrSurf_l[5]*fUpwind_l[7]+0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.31943828249996997*fUpwind_l[5]*alphaDrSurf_l[7]+0.4472135954999579*fUpwind_l[4]*alphaDrSurf_l[7]+0.5*fUpwind_l[0]*alphaDrSurf_l[7]+0.4*alphaDrSurf_l[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alphaDrSurf_l[6]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]+0.5000000000000001*fUpwind_l[1]*alphaDrSurf_l[5]+0.44721359549995804*alphaDrSurf_l[2]*fUpwind_l[3]+0.44721359549995804*fUpwind_l[2]*alphaDrSurf_l[3]; 
  Ghat_l[8] = 0.29814239699997186*alphaDrSurf_l[6]*fUpwind_l[10]+0.5*alphaDrSurf_l[2]*fUpwind_l[10]+0.29814239699997186*fUpwind_l[6]*alphaDrSurf_l[10]+0.5*fUpwind_l[2]*alphaDrSurf_l[10]+0.2981423969999719*alphaDrSurf_l[4]*fUpwind_l[8]+0.5*alphaDrSurf_l[0]*fUpwind_l[8]+0.2981423969999719*fUpwind_l[4]*alphaDrSurf_l[8]+0.5*fUpwind_l[0]*alphaDrSurf_l[8]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[6]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[6]+0.4391550328268398*alphaDrSurf_l[1]*fUpwind_l[4]+0.4391550328268398*fUpwind_l[1]*alphaDrSurf_l[4]; 
  Ghat_l[9] = 0.29814239699997186*alphaDrSurf_l[7]*fUpwind_l[11]+0.5*alphaDrSurf_l[1]*fUpwind_l[11]+0.29814239699997186*fUpwind_l[7]*alphaDrSurf_l[11]+0.5*fUpwind_l[1]*alphaDrSurf_l[11]+0.2981423969999719*alphaDrSurf_l[5]*fUpwind_l[9]+0.5*alphaDrSurf_l[0]*fUpwind_l[9]+0.2981423969999719*fUpwind_l[5]*alphaDrSurf_l[9]+0.5*fUpwind_l[0]*alphaDrSurf_l[9]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[7]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[7]+0.4391550328268398*alphaDrSurf_l[2]*fUpwind_l[5]+0.4391550328268398*fUpwind_l[2]*alphaDrSurf_l[5]; 
  Ghat_l[10] = 0.4472135954999579*alphaDrSurf_l[5]*fUpwind_l[10]+0.2981423969999719*alphaDrSurf_l[4]*fUpwind_l[10]+0.5*alphaDrSurf_l[0]*fUpwind_l[10]+0.4472135954999579*fUpwind_l[5]*alphaDrSurf_l[10]+0.2981423969999719*fUpwind_l[4]*alphaDrSurf_l[10]+0.5*fUpwind_l[0]*alphaDrSurf_l[10]+0.29814239699997186*alphaDrSurf_l[6]*fUpwind_l[8]+0.5*alphaDrSurf_l[2]*fUpwind_l[8]+0.29814239699997186*fUpwind_l[6]*alphaDrSurf_l[8]+0.5*fUpwind_l[2]*alphaDrSurf_l[8]+0.3927922024247863*alphaDrSurf_l[6]*fUpwind_l[7]+0.3927922024247863*fUpwind_l[6]*alphaDrSurf_l[7]+0.4391550328268399*alphaDrSurf_l[1]*fUpwind_l[6]+0.4391550328268399*fUpwind_l[1]*alphaDrSurf_l[6]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[4]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[4]; 
  Ghat_l[11] = 0.2981423969999719*alphaDrSurf_l[5]*fUpwind_l[11]+0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[11]+0.5*alphaDrSurf_l[0]*fUpwind_l[11]+0.2981423969999719*fUpwind_l[5]*alphaDrSurf_l[11]+0.4472135954999579*fUpwind_l[4]*alphaDrSurf_l[11]+0.5*fUpwind_l[0]*alphaDrSurf_l[11]+0.29814239699997186*alphaDrSurf_l[7]*fUpwind_l[9]+0.5*alphaDrSurf_l[1]*fUpwind_l[9]+0.29814239699997186*fUpwind_l[7]*alphaDrSurf_l[9]+0.5*fUpwind_l[1]*alphaDrSurf_l[9]+0.3927922024247863*alphaDrSurf_l[6]*fUpwind_l[7]+0.4391550328268399*alphaDrSurf_l[2]*fUpwind_l[7]+0.3927922024247863*fUpwind_l[6]*alphaDrSurf_l[7]+0.4391550328268399*fUpwind_l[2]*alphaDrSurf_l[7]+0.4391550328268399*alphaDrSurf_l[3]*fUpwind_l[5]+0.4391550328268399*fUpwind_l[3]*alphaDrSurf_l[5]; 

  Ghat_r[0] = 0.5*alphaDrSurf_r[11]*fUpwind_r[11]+0.5*alphaDrSurf_r[10]*fUpwind_r[10]+0.5*alphaDrSurf_r[9]*fUpwind_r[9]+0.5*alphaDrSurf_r[8]*fUpwind_r[8]+0.5*alphaDrSurf_r[7]*fUpwind_r[7]+0.5*alphaDrSurf_r[6]*fUpwind_r[6]+0.5*alphaDrSurf_r[5]*fUpwind_r[5]+0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[3]*fUpwind_r[3]+0.5*alphaDrSurf_r[2]*fUpwind_r[2]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5*alphaDrSurf_r[9]*fUpwind_r[11]+0.5*fUpwind_r[9]*alphaDrSurf_r[11]+0.4391550328268399*alphaDrSurf_r[6]*fUpwind_r[10]+0.4391550328268399*fUpwind_r[6]*alphaDrSurf_r[10]+0.4391550328268398*alphaDrSurf_r[4]*fUpwind_r[8]+0.4391550328268398*fUpwind_r[4]*alphaDrSurf_r[8]+0.5000000000000001*alphaDrSurf_r[5]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[5]*alphaDrSurf_r[7]+0.44721359549995804*alphaDrSurf_r[3]*fUpwind_r[6]+0.44721359549995804*fUpwind_r[3]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.4391550328268399*alphaDrSurf_r[7]*fUpwind_r[11]+0.4391550328268399*fUpwind_r[7]*alphaDrSurf_r[11]+0.5*alphaDrSurf_r[8]*fUpwind_r[10]+0.5*fUpwind_r[8]*alphaDrSurf_r[10]+0.4391550328268398*alphaDrSurf_r[5]*fUpwind_r[9]+0.4391550328268398*fUpwind_r[5]*alphaDrSurf_r[9]+0.44721359549995804*alphaDrSurf_r[3]*fUpwind_r[7]+0.44721359549995804*fUpwind_r[3]*alphaDrSurf_r[7]+0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[2]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[2]*alphaDrSurf_r[5]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = 0.4391550328268399*alphaDrSurf_r[5]*fUpwind_r[11]+0.4391550328268399*fUpwind_r[5]*alphaDrSurf_r[11]+0.4391550328268399*alphaDrSurf_r[4]*fUpwind_r[10]+0.4391550328268399*fUpwind_r[4]*alphaDrSurf_r[10]+0.4391550328268399*alphaDrSurf_r[7]*fUpwind_r[9]+0.4391550328268399*fUpwind_r[7]*alphaDrSurf_r[9]+0.4391550328268399*alphaDrSurf_r[6]*fUpwind_r[8]+0.4391550328268399*fUpwind_r[6]*alphaDrSurf_r[8]+0.4*alphaDrSurf_r[6]*fUpwind_r[7]+0.44721359549995804*alphaDrSurf_r[2]*fUpwind_r[7]+0.4*fUpwind_r[6]*alphaDrSurf_r[7]+0.44721359549995804*fUpwind_r[2]*alphaDrSurf_r[7]+0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[6]+0.44721359549995804*fUpwind_r[1]*alphaDrSurf_r[6]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[5]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alphaDrSurf_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Ghat_r[4] = 0.4472135954999579*alphaDrSurf_r[11]*fUpwind_r[11]+0.2981423969999719*alphaDrSurf_r[10]*fUpwind_r[10]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[10]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[10]+0.2981423969999719*alphaDrSurf_r[8]*fUpwind_r[8]+0.4391550328268398*alphaDrSurf_r[1]*fUpwind_r[8]+0.4391550328268398*fUpwind_r[1]*alphaDrSurf_r[8]+0.4472135954999579*alphaDrSurf_r[7]*fUpwind_r[7]+0.31943828249996997*alphaDrSurf_r[6]*fUpwind_r[6]+0.5000000000000001*alphaDrSurf_r[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[6]+0.31943828249996997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.2981423969999719*alphaDrSurf_r[11]*fUpwind_r[11]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[11]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[11]+0.4472135954999579*alphaDrSurf_r[10]*fUpwind_r[10]+0.2981423969999719*alphaDrSurf_r[9]*fUpwind_r[9]+0.4391550328268398*alphaDrSurf_r[2]*fUpwind_r[9]+0.4391550328268398*fUpwind_r[2]*alphaDrSurf_r[9]+0.31943828249996997*alphaDrSurf_r[7]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[1]*alphaDrSurf_r[7]+0.4472135954999579*alphaDrSurf_r[6]*fUpwind_r[6]+0.31943828249996997*alphaDrSurf_r[5]*fUpwind_r[5]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]+0.5*fUpwind_r[0]*alphaDrSurf_r[5]+0.4472135954999579*alphaDrSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaDrSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.3927922024247863*alphaDrSurf_r[7]*fUpwind_r[11]+0.3927922024247863*fUpwind_r[7]*alphaDrSurf_r[11]+0.29814239699997186*alphaDrSurf_r[8]*fUpwind_r[10]+0.3927922024247863*alphaDrSurf_r[7]*fUpwind_r[10]+0.4391550328268399*alphaDrSurf_r[1]*fUpwind_r[10]+0.29814239699997186*fUpwind_r[8]*alphaDrSurf_r[10]+0.3927922024247863*fUpwind_r[7]*alphaDrSurf_r[10]+0.4391550328268399*fUpwind_r[1]*alphaDrSurf_r[10]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[8]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[8]+0.4*alphaDrSurf_r[3]*fUpwind_r[7]+0.4*fUpwind_r[3]*alphaDrSurf_r[7]+0.4472135954999579*alphaDrSurf_r[5]*fUpwind_r[6]+0.31943828249996997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alphaDrSurf_r[6]+0.31943828249996997*fUpwind_r[4]*alphaDrSurf_r[6]+0.5*fUpwind_r[0]*alphaDrSurf_r[6]+0.5000000000000001*alphaDrSurf_r[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.44721359549995804*alphaDrSurf_r[1]*fUpwind_r[3]+0.44721359549995804*fUpwind_r[1]*alphaDrSurf_r[3]; 
  Ghat_r[7] = 0.29814239699997186*alphaDrSurf_r[9]*fUpwind_r[11]+0.3927922024247863*alphaDrSurf_r[6]*fUpwind_r[11]+0.4391550328268399*alphaDrSurf_r[2]*fUpwind_r[11]+0.29814239699997186*fUpwind_r[9]*alphaDrSurf_r[11]+0.3927922024247863*fUpwind_r[6]*alphaDrSurf_r[11]+0.4391550328268399*fUpwind_r[2]*alphaDrSurf_r[11]+0.3927922024247863*alphaDrSurf_r[6]*fUpwind_r[10]+0.3927922024247863*fUpwind_r[6]*alphaDrSurf_r[10]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[9]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[9]+0.31943828249996997*alphaDrSurf_r[5]*fUpwind_r[7]+0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.31943828249996997*fUpwind_r[5]*alphaDrSurf_r[7]+0.4472135954999579*fUpwind_r[4]*alphaDrSurf_r[7]+0.5*fUpwind_r[0]*alphaDrSurf_r[7]+0.4*alphaDrSurf_r[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alphaDrSurf_r[6]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]+0.5000000000000001*fUpwind_r[1]*alphaDrSurf_r[5]+0.44721359549995804*alphaDrSurf_r[2]*fUpwind_r[3]+0.44721359549995804*fUpwind_r[2]*alphaDrSurf_r[3]; 
  Ghat_r[8] = 0.29814239699997186*alphaDrSurf_r[6]*fUpwind_r[10]+0.5*alphaDrSurf_r[2]*fUpwind_r[10]+0.29814239699997186*fUpwind_r[6]*alphaDrSurf_r[10]+0.5*fUpwind_r[2]*alphaDrSurf_r[10]+0.2981423969999719*alphaDrSurf_r[4]*fUpwind_r[8]+0.5*alphaDrSurf_r[0]*fUpwind_r[8]+0.2981423969999719*fUpwind_r[4]*alphaDrSurf_r[8]+0.5*fUpwind_r[0]*alphaDrSurf_r[8]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[6]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[6]+0.4391550328268398*alphaDrSurf_r[1]*fUpwind_r[4]+0.4391550328268398*fUpwind_r[1]*alphaDrSurf_r[4]; 
  Ghat_r[9] = 0.29814239699997186*alphaDrSurf_r[7]*fUpwind_r[11]+0.5*alphaDrSurf_r[1]*fUpwind_r[11]+0.29814239699997186*fUpwind_r[7]*alphaDrSurf_r[11]+0.5*fUpwind_r[1]*alphaDrSurf_r[11]+0.2981423969999719*alphaDrSurf_r[5]*fUpwind_r[9]+0.5*alphaDrSurf_r[0]*fUpwind_r[9]+0.2981423969999719*fUpwind_r[5]*alphaDrSurf_r[9]+0.5*fUpwind_r[0]*alphaDrSurf_r[9]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[7]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[7]+0.4391550328268398*alphaDrSurf_r[2]*fUpwind_r[5]+0.4391550328268398*fUpwind_r[2]*alphaDrSurf_r[5]; 
  Ghat_r[10] = 0.4472135954999579*alphaDrSurf_r[5]*fUpwind_r[10]+0.2981423969999719*alphaDrSurf_r[4]*fUpwind_r[10]+0.5*alphaDrSurf_r[0]*fUpwind_r[10]+0.4472135954999579*fUpwind_r[5]*alphaDrSurf_r[10]+0.2981423969999719*fUpwind_r[4]*alphaDrSurf_r[10]+0.5*fUpwind_r[0]*alphaDrSurf_r[10]+0.29814239699997186*alphaDrSurf_r[6]*fUpwind_r[8]+0.5*alphaDrSurf_r[2]*fUpwind_r[8]+0.29814239699997186*fUpwind_r[6]*alphaDrSurf_r[8]+0.5*fUpwind_r[2]*alphaDrSurf_r[8]+0.3927922024247863*alphaDrSurf_r[6]*fUpwind_r[7]+0.3927922024247863*fUpwind_r[6]*alphaDrSurf_r[7]+0.4391550328268399*alphaDrSurf_r[1]*fUpwind_r[6]+0.4391550328268399*fUpwind_r[1]*alphaDrSurf_r[6]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[4]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[4]; 
  Ghat_r[11] = 0.2981423969999719*alphaDrSurf_r[5]*fUpwind_r[11]+0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[11]+0.5*alphaDrSurf_r[0]*fUpwind_r[11]+0.2981423969999719*fUpwind_r[5]*alphaDrSurf_r[11]+0.4472135954999579*fUpwind_r[4]*alphaDrSurf_r[11]+0.5*fUpwind_r[0]*alphaDrSurf_r[11]+0.29814239699997186*alphaDrSurf_r[7]*fUpwind_r[9]+0.5*alphaDrSurf_r[1]*fUpwind_r[9]+0.29814239699997186*fUpwind_r[7]*alphaDrSurf_r[9]+0.5*fUpwind_r[1]*alphaDrSurf_r[9]+0.3927922024247863*alphaDrSurf_r[6]*fUpwind_r[7]+0.4391550328268399*alphaDrSurf_r[2]*fUpwind_r[7]+0.3927922024247863*fUpwind_r[6]*alphaDrSurf_r[7]+0.4391550328268399*fUpwind_r[2]*alphaDrSurf_r[7]+0.4391550328268399*alphaDrSurf_r[3]*fUpwind_r[5]+0.4391550328268399*fUpwind_r[3]*alphaDrSurf_r[5]; 

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
  out[17] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[19] += 1.8708286933869707*Ghat_r[0]*rdv2+1.8708286933869707*Ghat_l[0]*rdv2; 
  out[20] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[22] += 1.5811388300841895*Ghat_r[3]*rdv2-1.5811388300841895*Ghat_l[3]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[24] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[25] += 1.2247448713915887*Ghat_r[8]*rdv2+1.2247448713915887*Ghat_l[8]*rdv2; 
  out[26] += 1.2247448713915887*Ghat_r[9]*rdv2+1.2247448713915887*Ghat_l[9]*rdv2; 
  out[27] += 1.8708286933869707*Ghat_r[1]*rdv2+1.8708286933869707*Ghat_l[1]*rdv2; 
  out[28] += 1.8708286933869707*Ghat_r[2]*rdv2+1.8708286933869707*Ghat_l[2]*rdv2; 
  out[29] += 1.2247448713915887*Ghat_r[10]*rdv2+1.2247448713915887*Ghat_l[10]*rdv2; 
  out[30] += 1.2247448713915887*Ghat_r[11]*rdv2+1.2247448713915887*Ghat_l[11]*rdv2; 
  out[31] += 1.8708286933869707*Ghat_r[3]*rdv2+1.8708286933869707*Ghat_l[3]*rdv2; 

  return 0.;

} 
