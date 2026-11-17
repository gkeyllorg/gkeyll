#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_3x1v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)] in left/center/right cells. 
  // out:           Incremented output distribution functions in center cell. 

  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[24]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[24]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[24]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_l[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_l[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf_l[3] = nuSum[3]*wvpar-0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 
  alphaDrSurf_l[4] = nuSum[4]*wvpar-0.5*nuSum[4]*dvpar-1.0*sumNuUPar[4]; 
  alphaDrSurf_l[5] = nuSum[5]*wvpar-0.5*nuSum[5]*dvpar-1.0*sumNuUPar[5]; 
  alphaDrSurf_l[6] = nuSum[6]*wvpar-0.5*nuSum[6]*dvpar-1.0*sumNuUPar[6]; 
  alphaDrSurf_l[7] = nuSum[7]*wvpar-0.5*nuSum[7]*dvpar-1.0*sumNuUPar[7]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_r[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_r[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf_r[3] = nuSum[3]*wvpar+0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 
  alphaDrSurf_r[4] = nuSum[4]*wvpar+0.5*nuSum[4]*dvpar-1.0*sumNuUPar[4]; 
  alphaDrSurf_r[5] = nuSum[5]*wvpar+0.5*nuSum[5]*dvpar-1.0*sumNuUPar[5]; 
  alphaDrSurf_r[6] = nuSum[6]*wvpar+0.5*nuSum[6]*dvpar-1.0*sumNuUPar[6]; 
  alphaDrSurf_r[7] = nuSum[7]*wvpar+0.5*nuSum[7]*dvpar-1.0*sumNuUPar[7]; 

  double F_0_UpwindQuad_l[8] = {0.0};
  double F_0_UpwindQuad_r[8] = {0.0};
  double F_0_Upwind_l[8] = {0.0};
  double F_0_Upwind_r[8] = {0.0};
  double Ghat_F_0_l[8] = {0.0}; 
  double Ghat_F_0_r[8] = {0.0}; 
  double G_1_UpwindQuad_l[8] = {0.0};
  double G_1_UpwindQuad_r[8] = {0.0};
  double G_1_Upwind_l[8] = {0.0};
  double G_1_Upwind_r[8] = {0.0};
  double Ghat_G_1_l[8] = {0.0}; 
  double Ghat_G_1_r[8] = {0.0}; 

  if (-(0.3535533905932737*alphaDrSurf_l[7])+0.3535533905932737*(alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[4])-0.3535533905932737*(alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1c); 
  } 
  if (-(0.3535533905932737*alphaDrSurf_r[7])+0.3535533905932737*(alphaDrSurf_r[6]+alphaDrSurf_r[5]+alphaDrSurf_r[4])-0.3535533905932737*(alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1r); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[7]-0.3535533905932737*(alphaDrSurf_l[6]+alphaDrSurf_l[5])+0.3535533905932737*(alphaDrSurf_l[4]+alphaDrSurf_l[3])-0.3535533905932737*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1c); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[7]-0.3535533905932737*(alphaDrSurf_r[6]+alphaDrSurf_r[5])+0.3535533905932737*(alphaDrSurf_r[4]+alphaDrSurf_r[3])-0.3535533905932737*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1r); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[7]-0.3535533905932737*alphaDrSurf_l[6]+0.3535533905932737*alphaDrSurf_l[5]-0.3535533905932737*(alphaDrSurf_l[4]+alphaDrSurf_l[3])+0.3535533905932737*alphaDrSurf_l[2]-0.3535533905932737*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1c); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[7]-0.3535533905932737*alphaDrSurf_r[6]+0.3535533905932737*alphaDrSurf_r[5]-0.3535533905932737*(alphaDrSurf_r[4]+alphaDrSurf_r[3])+0.3535533905932737*alphaDrSurf_r[2]-0.3535533905932737*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1r); 
  } 
  if (-(0.3535533905932737*alphaDrSurf_l[7])+0.3535533905932737*alphaDrSurf_l[6]-0.3535533905932737*(alphaDrSurf_l[5]+alphaDrSurf_l[4])+0.3535533905932737*(alphaDrSurf_l[3]+alphaDrSurf_l[2])-0.3535533905932737*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1c); 
  } 
  if (-(0.3535533905932737*alphaDrSurf_r[7])+0.3535533905932737*alphaDrSurf_r[6]-0.3535533905932737*(alphaDrSurf_r[5]+alphaDrSurf_r[4])+0.3535533905932737*(alphaDrSurf_r[3]+alphaDrSurf_r[2])-0.3535533905932737*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6])-0.3535533905932737*(alphaDrSurf_l[5]+alphaDrSurf_l[4]+alphaDrSurf_l[3]+alphaDrSurf_l[2])+0.3535533905932737*(alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaDrSurf_r[7]+alphaDrSurf_r[6])-0.3535533905932737*(alphaDrSurf_r[5]+alphaDrSurf_r[4]+alphaDrSurf_r[3]+alphaDrSurf_r[2])+0.3535533905932737*(alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1r); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]))+0.3535533905932737*alphaDrSurf_l[5]-0.3535533905932737*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[3]-0.3535533905932737*alphaDrSurf_l[2]+0.3535533905932737*(alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1c); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf_r[7]+alphaDrSurf_r[6]))+0.3535533905932737*alphaDrSurf_r[5]-0.3535533905932737*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[3]-0.3535533905932737*alphaDrSurf_r[2]+0.3535533905932737*(alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1r); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]))+0.3535533905932737*alphaDrSurf_l[4]-0.3535533905932737*alphaDrSurf_l[3]+0.3535533905932737*(alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0l); 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0c); 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1c); 
  } 
  if (-(0.3535533905932737*(alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[5]))+0.3535533905932737*alphaDrSurf_r[4]-0.3535533905932737*alphaDrSurf_r[3]+0.3535533905932737*(alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0c); 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0r); 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[5]+alphaDrSurf_l[4]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0l); 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0c); 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[5]+alphaDrSurf_r[4]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0c); 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0r); 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.3535533905932737*(F_0_Upwind_l[7]*alphaDrSurf_l[7]+F_0_Upwind_l[6]*alphaDrSurf_l[6]+F_0_Upwind_l[5]*alphaDrSurf_l[5]+F_0_Upwind_l[4]*alphaDrSurf_l[4]+F_0_Upwind_l[3]*alphaDrSurf_l[3]+F_0_Upwind_l[2]*alphaDrSurf_l[2]+F_0_Upwind_l[1]*alphaDrSurf_l[1]+F_0_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.3535533905932737*(F_0_Upwind_l[6]*alphaDrSurf_l[7]+alphaDrSurf_l[6]*F_0_Upwind_l[7]+F_0_Upwind_l[3]*alphaDrSurf_l[5]+alphaDrSurf_l[3]*F_0_Upwind_l[5]+F_0_Upwind_l[2]*alphaDrSurf_l[4]+alphaDrSurf_l[2]*F_0_Upwind_l[4]+F_0_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0_Upwind_l[1]); 
  Ghat_F_0_l[2] = 0.3535533905932737*(F_0_Upwind_l[5]*alphaDrSurf_l[7]+alphaDrSurf_l[5]*F_0_Upwind_l[7]+F_0_Upwind_l[3]*alphaDrSurf_l[6]+alphaDrSurf_l[3]*F_0_Upwind_l[6]+F_0_Upwind_l[1]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*F_0_Upwind_l[4]+F_0_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*F_0_Upwind_l[2]); 
  Ghat_F_0_l[3] = 0.3535533905932737*(F_0_Upwind_l[4]*alphaDrSurf_l[7]+alphaDrSurf_l[4]*F_0_Upwind_l[7]+F_0_Upwind_l[2]*alphaDrSurf_l[6]+alphaDrSurf_l[2]*F_0_Upwind_l[6]+F_0_Upwind_l[1]*alphaDrSurf_l[5]+alphaDrSurf_l[1]*F_0_Upwind_l[5]+F_0_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*F_0_Upwind_l[3]); 
  Ghat_F_0_l[4] = 0.3535533905932737*(F_0_Upwind_l[3]*alphaDrSurf_l[7]+alphaDrSurf_l[3]*F_0_Upwind_l[7]+F_0_Upwind_l[5]*alphaDrSurf_l[6]+alphaDrSurf_l[5]*F_0_Upwind_l[6]+F_0_Upwind_l[0]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*F_0_Upwind_l[4]+F_0_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*F_0_Upwind_l[2]); 
  Ghat_F_0_l[5] = 0.3535533905932737*(F_0_Upwind_l[2]*alphaDrSurf_l[7]+alphaDrSurf_l[2]*F_0_Upwind_l[7]+F_0_Upwind_l[4]*alphaDrSurf_l[6]+alphaDrSurf_l[4]*F_0_Upwind_l[6]+F_0_Upwind_l[0]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*F_0_Upwind_l[5]+F_0_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*F_0_Upwind_l[3]); 
  Ghat_F_0_l[6] = 0.3535533905932737*(F_0_Upwind_l[1]*alphaDrSurf_l[7]+alphaDrSurf_l[1]*F_0_Upwind_l[7]+F_0_Upwind_l[0]*alphaDrSurf_l[6]+alphaDrSurf_l[0]*F_0_Upwind_l[6]+F_0_Upwind_l[4]*alphaDrSurf_l[5]+alphaDrSurf_l[4]*F_0_Upwind_l[5]+F_0_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*F_0_Upwind_l[3]); 
  Ghat_F_0_l[7] = 0.3535533905932737*(F_0_Upwind_l[0]*alphaDrSurf_l[7]+alphaDrSurf_l[0]*F_0_Upwind_l[7]+F_0_Upwind_l[1]*alphaDrSurf_l[6]+alphaDrSurf_l[1]*F_0_Upwind_l[6]+F_0_Upwind_l[2]*alphaDrSurf_l[5]+alphaDrSurf_l[2]*F_0_Upwind_l[5]+F_0_Upwind_l[3]*alphaDrSurf_l[4]+alphaDrSurf_l[3]*F_0_Upwind_l[4]); 
  Ghat_G_1_l[0] = 0.3535533905932737*(G_1_Upwind_l[7]*alphaDrSurf_l[7]+G_1_Upwind_l[6]*alphaDrSurf_l[6]+G_1_Upwind_l[5]*alphaDrSurf_l[5]+G_1_Upwind_l[4]*alphaDrSurf_l[4]+G_1_Upwind_l[3]*alphaDrSurf_l[3]+G_1_Upwind_l[2]*alphaDrSurf_l[2]+G_1_Upwind_l[1]*alphaDrSurf_l[1]+G_1_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.3535533905932737*(G_1_Upwind_l[6]*alphaDrSurf_l[7]+alphaDrSurf_l[6]*G_1_Upwind_l[7]+G_1_Upwind_l[3]*alphaDrSurf_l[5]+alphaDrSurf_l[3]*G_1_Upwind_l[5]+G_1_Upwind_l[2]*alphaDrSurf_l[4]+alphaDrSurf_l[2]*G_1_Upwind_l[4]+G_1_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1_Upwind_l[1]); 
  Ghat_G_1_l[2] = 0.3535533905932737*(G_1_Upwind_l[5]*alphaDrSurf_l[7]+alphaDrSurf_l[5]*G_1_Upwind_l[7]+G_1_Upwind_l[3]*alphaDrSurf_l[6]+alphaDrSurf_l[3]*G_1_Upwind_l[6]+G_1_Upwind_l[1]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*G_1_Upwind_l[4]+G_1_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*G_1_Upwind_l[2]); 
  Ghat_G_1_l[3] = 0.3535533905932737*(G_1_Upwind_l[4]*alphaDrSurf_l[7]+alphaDrSurf_l[4]*G_1_Upwind_l[7]+G_1_Upwind_l[2]*alphaDrSurf_l[6]+alphaDrSurf_l[2]*G_1_Upwind_l[6]+G_1_Upwind_l[1]*alphaDrSurf_l[5]+alphaDrSurf_l[1]*G_1_Upwind_l[5]+G_1_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*G_1_Upwind_l[3]); 
  Ghat_G_1_l[4] = 0.3535533905932737*(G_1_Upwind_l[3]*alphaDrSurf_l[7]+alphaDrSurf_l[3]*G_1_Upwind_l[7]+G_1_Upwind_l[5]*alphaDrSurf_l[6]+alphaDrSurf_l[5]*G_1_Upwind_l[6]+G_1_Upwind_l[0]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*G_1_Upwind_l[4]+G_1_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*G_1_Upwind_l[2]); 
  Ghat_G_1_l[5] = 0.3535533905932737*(G_1_Upwind_l[2]*alphaDrSurf_l[7]+alphaDrSurf_l[2]*G_1_Upwind_l[7]+G_1_Upwind_l[4]*alphaDrSurf_l[6]+alphaDrSurf_l[4]*G_1_Upwind_l[6]+G_1_Upwind_l[0]*alphaDrSurf_l[5]+alphaDrSurf_l[0]*G_1_Upwind_l[5]+G_1_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*G_1_Upwind_l[3]); 
  Ghat_G_1_l[6] = 0.3535533905932737*(G_1_Upwind_l[1]*alphaDrSurf_l[7]+alphaDrSurf_l[1]*G_1_Upwind_l[7]+G_1_Upwind_l[0]*alphaDrSurf_l[6]+alphaDrSurf_l[0]*G_1_Upwind_l[6]+G_1_Upwind_l[4]*alphaDrSurf_l[5]+alphaDrSurf_l[4]*G_1_Upwind_l[5]+G_1_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*G_1_Upwind_l[3]); 
  Ghat_G_1_l[7] = 0.3535533905932737*(G_1_Upwind_l[0]*alphaDrSurf_l[7]+alphaDrSurf_l[0]*G_1_Upwind_l[7]+G_1_Upwind_l[1]*alphaDrSurf_l[6]+alphaDrSurf_l[1]*G_1_Upwind_l[6]+G_1_Upwind_l[2]*alphaDrSurf_l[5]+alphaDrSurf_l[2]*G_1_Upwind_l[5]+G_1_Upwind_l[3]*alphaDrSurf_l[4]+alphaDrSurf_l[3]*G_1_Upwind_l[4]); 

  Ghat_F_0_r[0] = 0.3535533905932737*(F_0_Upwind_r[7]*alphaDrSurf_r[7]+F_0_Upwind_r[6]*alphaDrSurf_r[6]+F_0_Upwind_r[5]*alphaDrSurf_r[5]+F_0_Upwind_r[4]*alphaDrSurf_r[4]+F_0_Upwind_r[3]*alphaDrSurf_r[3]+F_0_Upwind_r[2]*alphaDrSurf_r[2]+F_0_Upwind_r[1]*alphaDrSurf_r[1]+F_0_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.3535533905932737*(F_0_Upwind_r[6]*alphaDrSurf_r[7]+alphaDrSurf_r[6]*F_0_Upwind_r[7]+F_0_Upwind_r[3]*alphaDrSurf_r[5]+alphaDrSurf_r[3]*F_0_Upwind_r[5]+F_0_Upwind_r[2]*alphaDrSurf_r[4]+alphaDrSurf_r[2]*F_0_Upwind_r[4]+F_0_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0_Upwind_r[1]); 
  Ghat_F_0_r[2] = 0.3535533905932737*(F_0_Upwind_r[5]*alphaDrSurf_r[7]+alphaDrSurf_r[5]*F_0_Upwind_r[7]+F_0_Upwind_r[3]*alphaDrSurf_r[6]+alphaDrSurf_r[3]*F_0_Upwind_r[6]+F_0_Upwind_r[1]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*F_0_Upwind_r[4]+F_0_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*F_0_Upwind_r[2]); 
  Ghat_F_0_r[3] = 0.3535533905932737*(F_0_Upwind_r[4]*alphaDrSurf_r[7]+alphaDrSurf_r[4]*F_0_Upwind_r[7]+F_0_Upwind_r[2]*alphaDrSurf_r[6]+alphaDrSurf_r[2]*F_0_Upwind_r[6]+F_0_Upwind_r[1]*alphaDrSurf_r[5]+alphaDrSurf_r[1]*F_0_Upwind_r[5]+F_0_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*F_0_Upwind_r[3]); 
  Ghat_F_0_r[4] = 0.3535533905932737*(F_0_Upwind_r[3]*alphaDrSurf_r[7]+alphaDrSurf_r[3]*F_0_Upwind_r[7]+F_0_Upwind_r[5]*alphaDrSurf_r[6]+alphaDrSurf_r[5]*F_0_Upwind_r[6]+F_0_Upwind_r[0]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*F_0_Upwind_r[4]+F_0_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*F_0_Upwind_r[2]); 
  Ghat_F_0_r[5] = 0.3535533905932737*(F_0_Upwind_r[2]*alphaDrSurf_r[7]+alphaDrSurf_r[2]*F_0_Upwind_r[7]+F_0_Upwind_r[4]*alphaDrSurf_r[6]+alphaDrSurf_r[4]*F_0_Upwind_r[6]+F_0_Upwind_r[0]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*F_0_Upwind_r[5]+F_0_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*F_0_Upwind_r[3]); 
  Ghat_F_0_r[6] = 0.3535533905932737*(F_0_Upwind_r[1]*alphaDrSurf_r[7]+alphaDrSurf_r[1]*F_0_Upwind_r[7]+F_0_Upwind_r[0]*alphaDrSurf_r[6]+alphaDrSurf_r[0]*F_0_Upwind_r[6]+F_0_Upwind_r[4]*alphaDrSurf_r[5]+alphaDrSurf_r[4]*F_0_Upwind_r[5]+F_0_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*F_0_Upwind_r[3]); 
  Ghat_F_0_r[7] = 0.3535533905932737*(F_0_Upwind_r[0]*alphaDrSurf_r[7]+alphaDrSurf_r[0]*F_0_Upwind_r[7]+F_0_Upwind_r[1]*alphaDrSurf_r[6]+alphaDrSurf_r[1]*F_0_Upwind_r[6]+F_0_Upwind_r[2]*alphaDrSurf_r[5]+alphaDrSurf_r[2]*F_0_Upwind_r[5]+F_0_Upwind_r[3]*alphaDrSurf_r[4]+alphaDrSurf_r[3]*F_0_Upwind_r[4]); 
  Ghat_G_1_r[0] = 0.3535533905932737*(G_1_Upwind_r[7]*alphaDrSurf_r[7]+G_1_Upwind_r[6]*alphaDrSurf_r[6]+G_1_Upwind_r[5]*alphaDrSurf_r[5]+G_1_Upwind_r[4]*alphaDrSurf_r[4]+G_1_Upwind_r[3]*alphaDrSurf_r[3]+G_1_Upwind_r[2]*alphaDrSurf_r[2]+G_1_Upwind_r[1]*alphaDrSurf_r[1]+G_1_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.3535533905932737*(G_1_Upwind_r[6]*alphaDrSurf_r[7]+alphaDrSurf_r[6]*G_1_Upwind_r[7]+G_1_Upwind_r[3]*alphaDrSurf_r[5]+alphaDrSurf_r[3]*G_1_Upwind_r[5]+G_1_Upwind_r[2]*alphaDrSurf_r[4]+alphaDrSurf_r[2]*G_1_Upwind_r[4]+G_1_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1_Upwind_r[1]); 
  Ghat_G_1_r[2] = 0.3535533905932737*(G_1_Upwind_r[5]*alphaDrSurf_r[7]+alphaDrSurf_r[5]*G_1_Upwind_r[7]+G_1_Upwind_r[3]*alphaDrSurf_r[6]+alphaDrSurf_r[3]*G_1_Upwind_r[6]+G_1_Upwind_r[1]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*G_1_Upwind_r[4]+G_1_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*G_1_Upwind_r[2]); 
  Ghat_G_1_r[3] = 0.3535533905932737*(G_1_Upwind_r[4]*alphaDrSurf_r[7]+alphaDrSurf_r[4]*G_1_Upwind_r[7]+G_1_Upwind_r[2]*alphaDrSurf_r[6]+alphaDrSurf_r[2]*G_1_Upwind_r[6]+G_1_Upwind_r[1]*alphaDrSurf_r[5]+alphaDrSurf_r[1]*G_1_Upwind_r[5]+G_1_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*G_1_Upwind_r[3]); 
  Ghat_G_1_r[4] = 0.3535533905932737*(G_1_Upwind_r[3]*alphaDrSurf_r[7]+alphaDrSurf_r[3]*G_1_Upwind_r[7]+G_1_Upwind_r[5]*alphaDrSurf_r[6]+alphaDrSurf_r[5]*G_1_Upwind_r[6]+G_1_Upwind_r[0]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*G_1_Upwind_r[4]+G_1_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*G_1_Upwind_r[2]); 
  Ghat_G_1_r[5] = 0.3535533905932737*(G_1_Upwind_r[2]*alphaDrSurf_r[7]+alphaDrSurf_r[2]*G_1_Upwind_r[7]+G_1_Upwind_r[4]*alphaDrSurf_r[6]+alphaDrSurf_r[4]*G_1_Upwind_r[6]+G_1_Upwind_r[0]*alphaDrSurf_r[5]+alphaDrSurf_r[0]*G_1_Upwind_r[5]+G_1_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*G_1_Upwind_r[3]); 
  Ghat_G_1_r[6] = 0.3535533905932737*(G_1_Upwind_r[1]*alphaDrSurf_r[7]+alphaDrSurf_r[1]*G_1_Upwind_r[7]+G_1_Upwind_r[0]*alphaDrSurf_r[6]+alphaDrSurf_r[0]*G_1_Upwind_r[6]+G_1_Upwind_r[4]*alphaDrSurf_r[5]+alphaDrSurf_r[4]*G_1_Upwind_r[5]+G_1_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*G_1_Upwind_r[3]); 
  Ghat_G_1_r[7] = 0.3535533905932737*(G_1_Upwind_r[0]*alphaDrSurf_r[7]+alphaDrSurf_r[0]*G_1_Upwind_r[7]+G_1_Upwind_r[1]*alphaDrSurf_r[6]+alphaDrSurf_r[1]*G_1_Upwind_r[6]+G_1_Upwind_r[2]*alphaDrSurf_r[5]+alphaDrSurf_r[2]*G_1_Upwind_r[5]+G_1_Upwind_r[3]*alphaDrSurf_r[4]+alphaDrSurf_r[3]*G_1_Upwind_r[4]); 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_r[2]-0.7071067811865475*Ghat_F_0_l[2])*dv1par; 
  out_F_0[3] += (0.7071067811865475*Ghat_F_0_r[3]-0.7071067811865475*Ghat_F_0_l[3])*dv1par; 
  out_F_0[4] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[5] += (0.7071067811865475*Ghat_F_0_r[4]-0.7071067811865475*Ghat_F_0_l[4])*dv1par; 
  out_F_0[6] += (0.7071067811865475*Ghat_F_0_r[5]-0.7071067811865475*Ghat_F_0_l[5])*dv1par; 
  out_F_0[7] += (0.7071067811865475*Ghat_F_0_r[6]-0.7071067811865475*Ghat_F_0_l[6])*dv1par; 
  out_F_0[8] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[9] += 1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dv1par; 
  out_F_0[10] += 1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dv1par; 
  out_F_0[11] += (1.5811388300841895*Ghat_F_0_r[0]-1.5811388300841895*Ghat_F_0_l[0])*dv1par; 
  out_F_0[12] += (0.7071067811865475*Ghat_F_0_r[7]-0.7071067811865475*Ghat_F_0_l[7])*dv1par; 
  out_F_0[13] += 1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4])*dv1par; 
  out_F_0[14] += 1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5])*dv1par; 
  out_F_0[15] += 1.224744871391589*(Ghat_F_0_r[6]+Ghat_F_0_l[6])*dv1par; 
  out_F_0[16] += (1.5811388300841898*Ghat_F_0_r[1]-1.5811388300841898*Ghat_F_0_l[1])*dv1par; 
  out_F_0[17] += (1.5811388300841898*Ghat_F_0_r[2]-1.5811388300841898*Ghat_F_0_l[2])*dv1par; 
  out_F_0[18] += (1.5811388300841898*Ghat_F_0_r[3]-1.5811388300841898*Ghat_F_0_l[3])*dv1par; 
  out_F_0[19] += 1.224744871391589*(Ghat_F_0_r[7]+Ghat_F_0_l[7])*dv1par; 
  out_F_0[20] += (1.5811388300841895*Ghat_F_0_r[4]-1.5811388300841895*Ghat_F_0_l[4])*dv1par; 
  out_F_0[21] += (1.5811388300841895*Ghat_F_0_r[5]-1.5811388300841895*Ghat_F_0_l[5])*dv1par; 
  out_F_0[22] += (1.5811388300841895*Ghat_F_0_r[6]-1.5811388300841895*Ghat_F_0_l[6])*dv1par; 
  out_F_0[23] += (1.5811388300841898*Ghat_F_0_r[7]-1.5811388300841898*Ghat_F_0_l[7])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_r[2]-0.7071067811865475*Ghat_G_1_l[2])*dv1par; 
  out_G_1[3] += (0.7071067811865475*Ghat_G_1_r[3]-0.7071067811865475*Ghat_G_1_l[3])*dv1par; 
  out_G_1[4] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_r[4]-0.7071067811865475*Ghat_G_1_l[4])*dv1par; 
  out_G_1[6] += (0.7071067811865475*Ghat_G_1_r[5]-0.7071067811865475*Ghat_G_1_l[5])*dv1par; 
  out_G_1[7] += (0.7071067811865475*Ghat_G_1_r[6]-0.7071067811865475*Ghat_G_1_l[6])*dv1par; 
  out_G_1[8] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[9] += 1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dv1par; 
  out_G_1[10] += 1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dv1par; 
  out_G_1[11] += (1.5811388300841895*Ghat_G_1_r[0]-1.5811388300841895*Ghat_G_1_l[0])*dv1par; 
  out_G_1[12] += (0.7071067811865475*Ghat_G_1_r[7]-0.7071067811865475*Ghat_G_1_l[7])*dv1par; 
  out_G_1[13] += 1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4])*dv1par; 
  out_G_1[14] += 1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5])*dv1par; 
  out_G_1[15] += 1.224744871391589*(Ghat_G_1_r[6]+Ghat_G_1_l[6])*dv1par; 
  out_G_1[16] += (1.5811388300841898*Ghat_G_1_r[1]-1.5811388300841898*Ghat_G_1_l[1])*dv1par; 
  out_G_1[17] += (1.5811388300841898*Ghat_G_1_r[2]-1.5811388300841898*Ghat_G_1_l[2])*dv1par; 
  out_G_1[18] += (1.5811388300841898*Ghat_G_1_r[3]-1.5811388300841898*Ghat_G_1_l[3])*dv1par; 
  out_G_1[19] += 1.224744871391589*(Ghat_G_1_r[7]+Ghat_G_1_l[7])*dv1par; 
  out_G_1[20] += (1.5811388300841895*Ghat_G_1_r[4]-1.5811388300841895*Ghat_G_1_l[4])*dv1par; 
  out_G_1[21] += (1.5811388300841895*Ghat_G_1_r[5]-1.5811388300841895*Ghat_G_1_l[5])*dv1par; 
  out_G_1[22] += (1.5811388300841895*Ghat_G_1_r[6]-1.5811388300841895*Ghat_G_1_l[6])*dv1par; 
  out_G_1[23] += (1.5811388300841898*Ghat_G_1_r[7]-1.5811388300841898*Ghat_G_1_l[7])*dv1par; 

  return 0.;

} 
