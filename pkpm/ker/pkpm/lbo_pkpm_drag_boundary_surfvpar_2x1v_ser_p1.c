#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1), M0_perp, N_perp] in skin cell/last edge cell. 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Skin = &fSkin[12]; 
  const double *G_1Edge = &fEdge[12]; 
  const double *M0_xSkin = &fSkin[24]; 
  const double *M0_xEdge = &fEdge[24]; 
  const double *M0_ySkin = &fSkin[36]; 
  const double *M0_yEdge = &fEdge[36]; 
  const double *M0_zSkin = &fSkin[48]; 
  const double *M0_zEdge = &fEdge[48]; 
  const double *N_xSkin = &fSkin[60]; 
  const double *N_xEdge = &fEdge[60]; 
  const double *N_ySkin = &fSkin[72]; 
  const double *N_yEdge = &fEdge[72]; 
  const double *N_zSkin = &fSkin[84]; 
  const double *N_zEdge = &fEdge[84]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  double *out_M0_x = &out[24]; 
  double *out_M0_y = &out[36]; 
  double *out_M0_z = &out[48]; 
  double *out_N_x = &out[60]; 
  double *out_N_y = &out[72]; 
  double *out_N_z = &out[84]; 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf[4] = {0.0}; 
  double F_0_UpwindQuad[4] = {0.0};
  double F_0_Upwind[4] = {0.0};;
  double Ghat_F_0[4] = {0.0}; 
  double G_1_UpwindQuad[4] = {0.0};
  double G_1_Upwind[4] = {0.0};;
  double Ghat_G_1[4] = {0.0}; 
  double M0_x_UpwindQuad[4] = {0.0};
  double M0_x_Upwind[4] = {0.0};;
  double Ghat_M0_x[4] = {0.0}; 
  double M0_y_UpwindQuad[4] = {0.0};
  double M0_y_Upwind[4] = {0.0};;
  double Ghat_M0_y[4] = {0.0}; 
  double M0_z_UpwindQuad[4] = {0.0};
  double M0_z_Upwind[4] = {0.0};;
  double Ghat_M0_z[4] = {0.0}; 
  double N_x_UpwindQuad[4] = {0.0};
  double N_x_Upwind[4] = {0.0};;
  double Ghat_N_x[4] = {0.0}; 
  double N_y_UpwindQuad[4] = {0.0};
  double N_y_Upwind[4] = {0.0};;
  double Ghat_N_y[4] = {0.0}; 
  double N_z_UpwindQuad[4] = {0.0};
  double N_z_Upwind[4] = {0.0};;
  double Ghat_N_z[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar+0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  if (0.5*alphaDrSurf[3]-0.5*(alphaDrSurf[2]+alphaDrSurf[1])+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Skin); 
    M0_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_xSkin); 
    M0_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_ySkin); 
    M0_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_zSkin); 
    N_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_xSkin); 
    N_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_ySkin); 
    N_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Edge); 
    M0_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_xEdge); 
    M0_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_yEdge); 
    M0_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_zEdge); 
    N_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_xEdge); 
    N_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_yEdge); 
    N_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_zEdge); 
  } 
  if (-(0.5*alphaDrSurf[3])+0.5*alphaDrSurf[2]-0.5*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Skin); 
    M0_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_xSkin); 
    M0_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_ySkin); 
    M0_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_zSkin); 
    N_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_xSkin); 
    N_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_ySkin); 
    N_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Edge); 
    M0_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_xEdge); 
    M0_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_yEdge); 
    M0_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_zEdge); 
    N_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_xEdge); 
    N_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_yEdge); 
    N_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_zEdge); 
  } 
  if (0.5*(alphaDrSurf[1]+alphaDrSurf[0])-0.5*(alphaDrSurf[3]+alphaDrSurf[2]) < 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Skin); 
    M0_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_xSkin); 
    M0_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_ySkin); 
    M0_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_zSkin); 
    N_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_xSkin); 
    N_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_ySkin); 
    N_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Edge); 
    M0_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_xEdge); 
    M0_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_yEdge); 
    M0_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_zEdge); 
    N_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_xEdge); 
    N_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_yEdge); 
    N_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_zEdge); 
  } 
  if (0.5*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Skin); 
    M0_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_xSkin); 
    M0_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_ySkin); 
    M0_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_zSkin); 
    N_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_xSkin); 
    N_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_ySkin); 
    N_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Edge); 
    M0_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_xEdge); 
    M0_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_yEdge); 
    M0_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_zEdge); 
    N_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_xEdge); 
    N_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_yEdge); 
    N_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_zEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad, M0_x_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad, M0_y_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad, M0_z_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad, N_x_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad, N_y_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad, N_z_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaDrSurf[3]+0.5*F_0_Upwind[2]*alphaDrSurf[2]+0.5*F_0_Upwind[1]*alphaDrSurf[1]+0.5*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaDrSurf[3]+0.5*G_1_Upwind[2]*alphaDrSurf[2]+0.5*G_1_Upwind[1]*alphaDrSurf[1]+0.5*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*G_1_Upwind[2]; 
  Ghat_M0_x[0] = 0.5*M0_x_Upwind[3]*alphaDrSurf[3]+0.5*M0_x_Upwind[2]*alphaDrSurf[2]+0.5*M0_x_Upwind[1]*alphaDrSurf[1]+0.5*M0_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_x[1] = 0.5*M0_x_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_x_Upwind[3]+0.5*M0_x_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_x_Upwind[1]; 
  Ghat_M0_x[2] = 0.5*M0_x_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_x_Upwind[3]+0.5*M0_x_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_x_Upwind[2]; 
  Ghat_M0_x[3] = 0.5*M0_x_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_x_Upwind[3]+0.5*M0_x_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_x_Upwind[2]; 
  Ghat_M0_y[0] = 0.5*M0_y_Upwind[3]*alphaDrSurf[3]+0.5*M0_y_Upwind[2]*alphaDrSurf[2]+0.5*M0_y_Upwind[1]*alphaDrSurf[1]+0.5*M0_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_y[1] = 0.5*M0_y_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_y_Upwind[3]+0.5*M0_y_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_y_Upwind[1]; 
  Ghat_M0_y[2] = 0.5*M0_y_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_y_Upwind[3]+0.5*M0_y_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_y_Upwind[2]; 
  Ghat_M0_y[3] = 0.5*M0_y_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_y_Upwind[3]+0.5*M0_y_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_y_Upwind[2]; 
  Ghat_M0_z[0] = 0.5*M0_z_Upwind[3]*alphaDrSurf[3]+0.5*M0_z_Upwind[2]*alphaDrSurf[2]+0.5*M0_z_Upwind[1]*alphaDrSurf[1]+0.5*M0_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_z[1] = 0.5*M0_z_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_z_Upwind[3]+0.5*M0_z_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_z_Upwind[1]; 
  Ghat_M0_z[2] = 0.5*M0_z_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_z_Upwind[3]+0.5*M0_z_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_z_Upwind[2]; 
  Ghat_M0_z[3] = 0.5*M0_z_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_z_Upwind[3]+0.5*M0_z_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_z_Upwind[2]; 
  Ghat_N_x[0] = 0.5*N_x_Upwind[3]*alphaDrSurf[3]+0.5*N_x_Upwind[2]*alphaDrSurf[2]+0.5*N_x_Upwind[1]*alphaDrSurf[1]+0.5*N_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_x[1] = 0.5*N_x_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_x_Upwind[3]+0.5*N_x_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_x_Upwind[1]; 
  Ghat_N_x[2] = 0.5*N_x_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_x_Upwind[3]+0.5*N_x_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_x_Upwind[2]; 
  Ghat_N_x[3] = 0.5*N_x_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_x_Upwind[3]+0.5*N_x_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_x_Upwind[2]; 
  Ghat_N_y[0] = 0.5*N_y_Upwind[3]*alphaDrSurf[3]+0.5*N_y_Upwind[2]*alphaDrSurf[2]+0.5*N_y_Upwind[1]*alphaDrSurf[1]+0.5*N_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_y[1] = 0.5*N_y_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_y_Upwind[3]+0.5*N_y_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_y_Upwind[1]; 
  Ghat_N_y[2] = 0.5*N_y_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_y_Upwind[3]+0.5*N_y_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_y_Upwind[2]; 
  Ghat_N_y[3] = 0.5*N_y_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_y_Upwind[3]+0.5*N_y_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_y_Upwind[2]; 
  Ghat_N_z[0] = 0.5*N_z_Upwind[3]*alphaDrSurf[3]+0.5*N_z_Upwind[2]*alphaDrSurf[2]+0.5*N_z_Upwind[1]*alphaDrSurf[1]+0.5*N_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_z[1] = 0.5*N_z_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_z_Upwind[3]+0.5*N_z_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_z_Upwind[1]; 
  Ghat_N_z[2] = 0.5*N_z_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_z_Upwind[3]+0.5*N_z_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_z_Upwind[2]; 
  Ghat_N_z[3] = 0.5*N_z_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_z_Upwind[3]+0.5*N_z_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_z_Upwind[2]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += 1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += 1.5811388300841898*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += 1.5811388300841895*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += 1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += 1.5811388300841898*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += 1.5811388300841895*Ghat_G_1[3]*dv1par; 
  out_M0_x[0] += 0.7071067811865475*Ghat_M0_x[0]*dv1par; 
  out_M0_x[1] += 0.7071067811865475*Ghat_M0_x[1]*dv1par; 
  out_M0_x[2] += 0.7071067811865475*Ghat_M0_x[2]*dv1par; 
  out_M0_x[3] += 1.224744871391589*Ghat_M0_x[0]*dv1par; 
  out_M0_x[4] += 0.7071067811865475*Ghat_M0_x[3]*dv1par; 
  out_M0_x[5] += 1.224744871391589*Ghat_M0_x[1]*dv1par; 
  out_M0_x[6] += 1.224744871391589*Ghat_M0_x[2]*dv1par; 
  out_M0_x[7] += 1.5811388300841895*Ghat_M0_x[0]*dv1par; 
  out_M0_x[8] += 1.224744871391589*Ghat_M0_x[3]*dv1par; 
  out_M0_x[9] += 1.5811388300841898*Ghat_M0_x[1]*dv1par; 
  out_M0_x[10] += 1.5811388300841898*Ghat_M0_x[2]*dv1par; 
  out_M0_x[11] += 1.5811388300841895*Ghat_M0_x[3]*dv1par; 
  out_M0_y[0] += 0.7071067811865475*Ghat_M0_y[0]*dv1par; 
  out_M0_y[1] += 0.7071067811865475*Ghat_M0_y[1]*dv1par; 
  out_M0_y[2] += 0.7071067811865475*Ghat_M0_y[2]*dv1par; 
  out_M0_y[3] += 1.224744871391589*Ghat_M0_y[0]*dv1par; 
  out_M0_y[4] += 0.7071067811865475*Ghat_M0_y[3]*dv1par; 
  out_M0_y[5] += 1.224744871391589*Ghat_M0_y[1]*dv1par; 
  out_M0_y[6] += 1.224744871391589*Ghat_M0_y[2]*dv1par; 
  out_M0_y[7] += 1.5811388300841895*Ghat_M0_y[0]*dv1par; 
  out_M0_y[8] += 1.224744871391589*Ghat_M0_y[3]*dv1par; 
  out_M0_y[9] += 1.5811388300841898*Ghat_M0_y[1]*dv1par; 
  out_M0_y[10] += 1.5811388300841898*Ghat_M0_y[2]*dv1par; 
  out_M0_y[11] += 1.5811388300841895*Ghat_M0_y[3]*dv1par; 
  out_M0_z[0] += 0.7071067811865475*Ghat_M0_z[0]*dv1par; 
  out_M0_z[1] += 0.7071067811865475*Ghat_M0_z[1]*dv1par; 
  out_M0_z[2] += 0.7071067811865475*Ghat_M0_z[2]*dv1par; 
  out_M0_z[3] += 1.224744871391589*Ghat_M0_z[0]*dv1par; 
  out_M0_z[4] += 0.7071067811865475*Ghat_M0_z[3]*dv1par; 
  out_M0_z[5] += 1.224744871391589*Ghat_M0_z[1]*dv1par; 
  out_M0_z[6] += 1.224744871391589*Ghat_M0_z[2]*dv1par; 
  out_M0_z[7] += 1.5811388300841895*Ghat_M0_z[0]*dv1par; 
  out_M0_z[8] += 1.224744871391589*Ghat_M0_z[3]*dv1par; 
  out_M0_z[9] += 1.5811388300841898*Ghat_M0_z[1]*dv1par; 
  out_M0_z[10] += 1.5811388300841898*Ghat_M0_z[2]*dv1par; 
  out_M0_z[11] += 1.5811388300841895*Ghat_M0_z[3]*dv1par; 
  out_N_x[0] += 0.7071067811865475*Ghat_N_x[0]*dv1par; 
  out_N_x[1] += 0.7071067811865475*Ghat_N_x[1]*dv1par; 
  out_N_x[2] += 0.7071067811865475*Ghat_N_x[2]*dv1par; 
  out_N_x[3] += 1.224744871391589*Ghat_N_x[0]*dv1par; 
  out_N_x[4] += 0.7071067811865475*Ghat_N_x[3]*dv1par; 
  out_N_x[5] += 1.224744871391589*Ghat_N_x[1]*dv1par; 
  out_N_x[6] += 1.224744871391589*Ghat_N_x[2]*dv1par; 
  out_N_x[7] += 1.5811388300841895*Ghat_N_x[0]*dv1par; 
  out_N_x[8] += 1.224744871391589*Ghat_N_x[3]*dv1par; 
  out_N_x[9] += 1.5811388300841898*Ghat_N_x[1]*dv1par; 
  out_N_x[10] += 1.5811388300841898*Ghat_N_x[2]*dv1par; 
  out_N_x[11] += 1.5811388300841895*Ghat_N_x[3]*dv1par; 
  out_N_y[0] += 0.7071067811865475*Ghat_N_y[0]*dv1par; 
  out_N_y[1] += 0.7071067811865475*Ghat_N_y[1]*dv1par; 
  out_N_y[2] += 0.7071067811865475*Ghat_N_y[2]*dv1par; 
  out_N_y[3] += 1.224744871391589*Ghat_N_y[0]*dv1par; 
  out_N_y[4] += 0.7071067811865475*Ghat_N_y[3]*dv1par; 
  out_N_y[5] += 1.224744871391589*Ghat_N_y[1]*dv1par; 
  out_N_y[6] += 1.224744871391589*Ghat_N_y[2]*dv1par; 
  out_N_y[7] += 1.5811388300841895*Ghat_N_y[0]*dv1par; 
  out_N_y[8] += 1.224744871391589*Ghat_N_y[3]*dv1par; 
  out_N_y[9] += 1.5811388300841898*Ghat_N_y[1]*dv1par; 
  out_N_y[10] += 1.5811388300841898*Ghat_N_y[2]*dv1par; 
  out_N_y[11] += 1.5811388300841895*Ghat_N_y[3]*dv1par; 
  out_N_z[0] += 0.7071067811865475*Ghat_N_z[0]*dv1par; 
  out_N_z[1] += 0.7071067811865475*Ghat_N_z[1]*dv1par; 
  out_N_z[2] += 0.7071067811865475*Ghat_N_z[2]*dv1par; 
  out_N_z[3] += 1.224744871391589*Ghat_N_z[0]*dv1par; 
  out_N_z[4] += 0.7071067811865475*Ghat_N_z[3]*dv1par; 
  out_N_z[5] += 1.224744871391589*Ghat_N_z[1]*dv1par; 
  out_N_z[6] += 1.224744871391589*Ghat_N_z[2]*dv1par; 
  out_N_z[7] += 1.5811388300841895*Ghat_N_z[0]*dv1par; 
  out_N_z[8] += 1.224744871391589*Ghat_N_z[3]*dv1par; 
  out_N_z[9] += 1.5811388300841898*Ghat_N_z[1]*dv1par; 
  out_N_z[10] += 1.5811388300841898*Ghat_N_z[2]*dv1par; 
  out_N_z[11] += 1.5811388300841895*Ghat_N_z[3]*dv1par; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf[3] = nuSum[3]*wvpar-0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  if (0.5*alphaDrSurf[3]-0.5*(alphaDrSurf[2]+alphaDrSurf[1])+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1Edge); 
    M0_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_xEdge); 
    M0_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_yEdge); 
    M0_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_zEdge); 
    N_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_xEdge); 
    N_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_yEdge); 
    N_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1Skin); 
    M0_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_xSkin); 
    M0_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_ySkin); 
    M0_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_zSkin); 
    N_x_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_xSkin); 
    N_y_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_ySkin); 
    N_z_UpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_zSkin); 
  } 
  if (-(0.5*alphaDrSurf[3])+0.5*alphaDrSurf[2]-0.5*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1Edge); 
    M0_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_xEdge); 
    M0_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_yEdge); 
    M0_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_zEdge); 
    N_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_xEdge); 
    N_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_yEdge); 
    N_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1Skin); 
    M0_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_xSkin); 
    M0_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_ySkin); 
    M0_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_zSkin); 
    N_x_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_xSkin); 
    N_y_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_ySkin); 
    N_z_UpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_zSkin); 
  } 
  if (0.5*(alphaDrSurf[1]+alphaDrSurf[0])-0.5*(alphaDrSurf[3]+alphaDrSurf[2]) < 0) { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0Edge); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1Edge); 
    M0_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_xEdge); 
    M0_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_yEdge); 
    M0_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_zEdge); 
    N_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_xEdge); 
    N_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_yEdge); 
    N_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0Skin); 
    G_1_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1Skin); 
    M0_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_xSkin); 
    M0_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_ySkin); 
    M0_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_zSkin); 
    N_x_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_xSkin); 
    N_y_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_ySkin); 
    N_z_UpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_zSkin); 
  } 
  if (0.5*(alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0Edge); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1Edge); 
    M0_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_xEdge); 
    M0_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_yEdge); 
    M0_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_zEdge); 
    N_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_xEdge); 
    N_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_yEdge); 
    N_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0Skin); 
    G_1_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1Skin); 
    M0_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_xSkin); 
    M0_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_ySkin); 
    M0_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_zSkin); 
    N_x_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_xSkin); 
    N_y_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_ySkin); 
    N_z_UpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_zSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad, M0_x_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad, M0_y_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad, M0_z_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad, N_x_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad, N_y_Upwind); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad, N_z_Upwind); 

  Ghat_F_0[0] = 0.5*F_0_Upwind[3]*alphaDrSurf[3]+0.5*F_0_Upwind[2]*alphaDrSurf[2]+0.5*F_0_Upwind[1]*alphaDrSurf[1]+0.5*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.5*F_0_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_F_0[2] = 0.5*F_0_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*F_0_Upwind[3]+0.5*F_0_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*F_0_Upwind[2]; 
  Ghat_F_0[3] = 0.5*F_0_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*F_0_Upwind[3]+0.5*F_0_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*F_0_Upwind[2]; 
  Ghat_G_1[0] = 0.5*G_1_Upwind[3]*alphaDrSurf[3]+0.5*G_1_Upwind[2]*alphaDrSurf[2]+0.5*G_1_Upwind[1]*alphaDrSurf[1]+0.5*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.5*G_1_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_G_1[2] = 0.5*G_1_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*G_1_Upwind[3]+0.5*G_1_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*G_1_Upwind[2]; 
  Ghat_G_1[3] = 0.5*G_1_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*G_1_Upwind[3]+0.5*G_1_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*G_1_Upwind[2]; 
  Ghat_M0_x[0] = 0.5*M0_x_Upwind[3]*alphaDrSurf[3]+0.5*M0_x_Upwind[2]*alphaDrSurf[2]+0.5*M0_x_Upwind[1]*alphaDrSurf[1]+0.5*M0_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_x[1] = 0.5*M0_x_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_x_Upwind[3]+0.5*M0_x_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_x_Upwind[1]; 
  Ghat_M0_x[2] = 0.5*M0_x_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_x_Upwind[3]+0.5*M0_x_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_x_Upwind[2]; 
  Ghat_M0_x[3] = 0.5*M0_x_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_x_Upwind[3]+0.5*M0_x_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_x_Upwind[2]; 
  Ghat_M0_y[0] = 0.5*M0_y_Upwind[3]*alphaDrSurf[3]+0.5*M0_y_Upwind[2]*alphaDrSurf[2]+0.5*M0_y_Upwind[1]*alphaDrSurf[1]+0.5*M0_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_y[1] = 0.5*M0_y_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_y_Upwind[3]+0.5*M0_y_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_y_Upwind[1]; 
  Ghat_M0_y[2] = 0.5*M0_y_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_y_Upwind[3]+0.5*M0_y_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_y_Upwind[2]; 
  Ghat_M0_y[3] = 0.5*M0_y_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_y_Upwind[3]+0.5*M0_y_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_y_Upwind[2]; 
  Ghat_M0_z[0] = 0.5*M0_z_Upwind[3]*alphaDrSurf[3]+0.5*M0_z_Upwind[2]*alphaDrSurf[2]+0.5*M0_z_Upwind[1]*alphaDrSurf[1]+0.5*M0_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_z[1] = 0.5*M0_z_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*M0_z_Upwind[3]+0.5*M0_z_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*M0_z_Upwind[1]; 
  Ghat_M0_z[2] = 0.5*M0_z_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*M0_z_Upwind[3]+0.5*M0_z_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*M0_z_Upwind[2]; 
  Ghat_M0_z[3] = 0.5*M0_z_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*M0_z_Upwind[3]+0.5*M0_z_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*M0_z_Upwind[2]; 
  Ghat_N_x[0] = 0.5*N_x_Upwind[3]*alphaDrSurf[3]+0.5*N_x_Upwind[2]*alphaDrSurf[2]+0.5*N_x_Upwind[1]*alphaDrSurf[1]+0.5*N_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_x[1] = 0.5*N_x_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_x_Upwind[3]+0.5*N_x_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_x_Upwind[1]; 
  Ghat_N_x[2] = 0.5*N_x_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_x_Upwind[3]+0.5*N_x_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_x_Upwind[2]; 
  Ghat_N_x[3] = 0.5*N_x_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_x_Upwind[3]+0.5*N_x_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_x_Upwind[2]; 
  Ghat_N_y[0] = 0.5*N_y_Upwind[3]*alphaDrSurf[3]+0.5*N_y_Upwind[2]*alphaDrSurf[2]+0.5*N_y_Upwind[1]*alphaDrSurf[1]+0.5*N_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_y[1] = 0.5*N_y_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_y_Upwind[3]+0.5*N_y_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_y_Upwind[1]; 
  Ghat_N_y[2] = 0.5*N_y_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_y_Upwind[3]+0.5*N_y_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_y_Upwind[2]; 
  Ghat_N_y[3] = 0.5*N_y_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_y_Upwind[3]+0.5*N_y_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_y_Upwind[2]; 
  Ghat_N_z[0] = 0.5*N_z_Upwind[3]*alphaDrSurf[3]+0.5*N_z_Upwind[2]*alphaDrSurf[2]+0.5*N_z_Upwind[1]*alphaDrSurf[1]+0.5*N_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_z[1] = 0.5*N_z_Upwind[2]*alphaDrSurf[3]+0.5*alphaDrSurf[2]*N_z_Upwind[3]+0.5*N_z_Upwind[0]*alphaDrSurf[1]+0.5*alphaDrSurf[0]*N_z_Upwind[1]; 
  Ghat_N_z[2] = 0.5*N_z_Upwind[1]*alphaDrSurf[3]+0.5*alphaDrSurf[1]*N_z_Upwind[3]+0.5*N_z_Upwind[0]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*N_z_Upwind[2]; 
  Ghat_N_z[3] = 0.5*N_z_Upwind[0]*alphaDrSurf[3]+0.5*alphaDrSurf[0]*N_z_Upwind[3]+0.5*N_z_Upwind[1]*alphaDrSurf[2]+0.5*alphaDrSurf[1]*N_z_Upwind[2]; 

  out_F_0[0] += -(0.7071067811865475*Ghat_F_0[0]*dv1par); 
  out_F_0[1] += -(0.7071067811865475*Ghat_F_0[1]*dv1par); 
  out_F_0[2] += -(0.7071067811865475*Ghat_F_0[2]*dv1par); 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += -(0.7071067811865475*Ghat_F_0[3]*dv1par); 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += -(1.5811388300841895*Ghat_F_0[0]*dv1par); 
  out_F_0[8] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[9] += -(1.5811388300841898*Ghat_F_0[1]*dv1par); 
  out_F_0[10] += -(1.5811388300841898*Ghat_F_0[2]*dv1par); 
  out_F_0[11] += -(1.5811388300841895*Ghat_F_0[3]*dv1par); 
  out_G_1[0] += -(0.7071067811865475*Ghat_G_1[0]*dv1par); 
  out_G_1[1] += -(0.7071067811865475*Ghat_G_1[1]*dv1par); 
  out_G_1[2] += -(0.7071067811865475*Ghat_G_1[2]*dv1par); 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += -(0.7071067811865475*Ghat_G_1[3]*dv1par); 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += -(1.5811388300841895*Ghat_G_1[0]*dv1par); 
  out_G_1[8] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[9] += -(1.5811388300841898*Ghat_G_1[1]*dv1par); 
  out_G_1[10] += -(1.5811388300841898*Ghat_G_1[2]*dv1par); 
  out_G_1[11] += -(1.5811388300841895*Ghat_G_1[3]*dv1par); 
  out_M0_x[0] += -(0.7071067811865475*Ghat_M0_x[0]*dv1par); 
  out_M0_x[1] += -(0.7071067811865475*Ghat_M0_x[1]*dv1par); 
  out_M0_x[2] += -(0.7071067811865475*Ghat_M0_x[2]*dv1par); 
  out_M0_x[3] += 1.224744871391589*Ghat_M0_x[0]*dv1par; 
  out_M0_x[4] += -(0.7071067811865475*Ghat_M0_x[3]*dv1par); 
  out_M0_x[5] += 1.224744871391589*Ghat_M0_x[1]*dv1par; 
  out_M0_x[6] += 1.224744871391589*Ghat_M0_x[2]*dv1par; 
  out_M0_x[7] += -(1.5811388300841895*Ghat_M0_x[0]*dv1par); 
  out_M0_x[8] += 1.224744871391589*Ghat_M0_x[3]*dv1par; 
  out_M0_x[9] += -(1.5811388300841898*Ghat_M0_x[1]*dv1par); 
  out_M0_x[10] += -(1.5811388300841898*Ghat_M0_x[2]*dv1par); 
  out_M0_x[11] += -(1.5811388300841895*Ghat_M0_x[3]*dv1par); 
  out_M0_y[0] += -(0.7071067811865475*Ghat_M0_y[0]*dv1par); 
  out_M0_y[1] += -(0.7071067811865475*Ghat_M0_y[1]*dv1par); 
  out_M0_y[2] += -(0.7071067811865475*Ghat_M0_y[2]*dv1par); 
  out_M0_y[3] += 1.224744871391589*Ghat_M0_y[0]*dv1par; 
  out_M0_y[4] += -(0.7071067811865475*Ghat_M0_y[3]*dv1par); 
  out_M0_y[5] += 1.224744871391589*Ghat_M0_y[1]*dv1par; 
  out_M0_y[6] += 1.224744871391589*Ghat_M0_y[2]*dv1par; 
  out_M0_y[7] += -(1.5811388300841895*Ghat_M0_y[0]*dv1par); 
  out_M0_y[8] += 1.224744871391589*Ghat_M0_y[3]*dv1par; 
  out_M0_y[9] += -(1.5811388300841898*Ghat_M0_y[1]*dv1par); 
  out_M0_y[10] += -(1.5811388300841898*Ghat_M0_y[2]*dv1par); 
  out_M0_y[11] += -(1.5811388300841895*Ghat_M0_y[3]*dv1par); 
  out_M0_z[0] += -(0.7071067811865475*Ghat_M0_z[0]*dv1par); 
  out_M0_z[1] += -(0.7071067811865475*Ghat_M0_z[1]*dv1par); 
  out_M0_z[2] += -(0.7071067811865475*Ghat_M0_z[2]*dv1par); 
  out_M0_z[3] += 1.224744871391589*Ghat_M0_z[0]*dv1par; 
  out_M0_z[4] += -(0.7071067811865475*Ghat_M0_z[3]*dv1par); 
  out_M0_z[5] += 1.224744871391589*Ghat_M0_z[1]*dv1par; 
  out_M0_z[6] += 1.224744871391589*Ghat_M0_z[2]*dv1par; 
  out_M0_z[7] += -(1.5811388300841895*Ghat_M0_z[0]*dv1par); 
  out_M0_z[8] += 1.224744871391589*Ghat_M0_z[3]*dv1par; 
  out_M0_z[9] += -(1.5811388300841898*Ghat_M0_z[1]*dv1par); 
  out_M0_z[10] += -(1.5811388300841898*Ghat_M0_z[2]*dv1par); 
  out_M0_z[11] += -(1.5811388300841895*Ghat_M0_z[3]*dv1par); 
  out_N_x[0] += -(0.7071067811865475*Ghat_N_x[0]*dv1par); 
  out_N_x[1] += -(0.7071067811865475*Ghat_N_x[1]*dv1par); 
  out_N_x[2] += -(0.7071067811865475*Ghat_N_x[2]*dv1par); 
  out_N_x[3] += 1.224744871391589*Ghat_N_x[0]*dv1par; 
  out_N_x[4] += -(0.7071067811865475*Ghat_N_x[3]*dv1par); 
  out_N_x[5] += 1.224744871391589*Ghat_N_x[1]*dv1par; 
  out_N_x[6] += 1.224744871391589*Ghat_N_x[2]*dv1par; 
  out_N_x[7] += -(1.5811388300841895*Ghat_N_x[0]*dv1par); 
  out_N_x[8] += 1.224744871391589*Ghat_N_x[3]*dv1par; 
  out_N_x[9] += -(1.5811388300841898*Ghat_N_x[1]*dv1par); 
  out_N_x[10] += -(1.5811388300841898*Ghat_N_x[2]*dv1par); 
  out_N_x[11] += -(1.5811388300841895*Ghat_N_x[3]*dv1par); 
  out_N_y[0] += -(0.7071067811865475*Ghat_N_y[0]*dv1par); 
  out_N_y[1] += -(0.7071067811865475*Ghat_N_y[1]*dv1par); 
  out_N_y[2] += -(0.7071067811865475*Ghat_N_y[2]*dv1par); 
  out_N_y[3] += 1.224744871391589*Ghat_N_y[0]*dv1par; 
  out_N_y[4] += -(0.7071067811865475*Ghat_N_y[3]*dv1par); 
  out_N_y[5] += 1.224744871391589*Ghat_N_y[1]*dv1par; 
  out_N_y[6] += 1.224744871391589*Ghat_N_y[2]*dv1par; 
  out_N_y[7] += -(1.5811388300841895*Ghat_N_y[0]*dv1par); 
  out_N_y[8] += 1.224744871391589*Ghat_N_y[3]*dv1par; 
  out_N_y[9] += -(1.5811388300841898*Ghat_N_y[1]*dv1par); 
  out_N_y[10] += -(1.5811388300841898*Ghat_N_y[2]*dv1par); 
  out_N_y[11] += -(1.5811388300841895*Ghat_N_y[3]*dv1par); 
  out_N_z[0] += -(0.7071067811865475*Ghat_N_z[0]*dv1par); 
  out_N_z[1] += -(0.7071067811865475*Ghat_N_z[1]*dv1par); 
  out_N_z[2] += -(0.7071067811865475*Ghat_N_z[2]*dv1par); 
  out_N_z[3] += 1.224744871391589*Ghat_N_z[0]*dv1par; 
  out_N_z[4] += -(0.7071067811865475*Ghat_N_z[3]*dv1par); 
  out_N_z[5] += 1.224744871391589*Ghat_N_z[1]*dv1par; 
  out_N_z[6] += 1.224744871391589*Ghat_N_z[2]*dv1par; 
  out_N_z[7] += -(1.5811388300841895*Ghat_N_z[0]*dv1par); 
  out_N_z[8] += 1.224744871391589*Ghat_N_z[3]*dv1par; 
  out_N_z[9] += -(1.5811388300841898*Ghat_N_z[1]*dv1par); 
  out_N_z[10] += -(1.5811388300841898*Ghat_N_z[2]*dv1par); 
  out_N_z[11] += -(1.5811388300841895*Ghat_N_z[3]*dv1par); 

  }

  return 0.;

} 
