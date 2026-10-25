#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/fEdge:   Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1), M0_perp, N_perp] in skin cell/last edge cell. 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Skin = &fSkin[6]; 
  const double *G_1Edge = &fEdge[6]; 
  const double *M0_xSkin = &fSkin[12]; 
  const double *M0_xEdge = &fEdge[12]; 
  const double *M0_ySkin = &fSkin[18]; 
  const double *M0_yEdge = &fEdge[18]; 
  const double *M0_zSkin = &fSkin[24]; 
  const double *M0_zEdge = &fEdge[24]; 
  const double *N_xSkin = &fSkin[30]; 
  const double *N_xEdge = &fEdge[30]; 
  const double *N_ySkin = &fSkin[36]; 
  const double *N_yEdge = &fEdge[36]; 
  const double *N_zSkin = &fSkin[42]; 
  const double *N_zEdge = &fEdge[42]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double *out_M0_x = &out[12]; 
  double *out_M0_y = &out[18]; 
  double *out_M0_z = &out[24]; 
  double *out_N_x = &out[30]; 
  double *out_N_y = &out[36]; 
  double *out_N_z = &out[42]; 
  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf[2] = {0.0}; 
  double F_0_UpwindQuad[2] = {0.0};
  double F_0_Upwind[2] = {0.0};;
  double Ghat_F_0[2] = {0.0}; 
  double G_1_UpwindQuad[2] = {0.0};
  double G_1_Upwind[2] = {0.0};;
  double Ghat_G_1[2] = {0.0}; 
  double M0_x_UpwindQuad[2] = {0.0};
  double M0_x_Upwind[2] = {0.0};;
  double Ghat_M0_x[2] = {0.0}; 
  double M0_y_UpwindQuad[2] = {0.0};
  double M0_y_Upwind[2] = {0.0};;
  double Ghat_M0_y[2] = {0.0}; 
  double M0_z_UpwindQuad[2] = {0.0};
  double M0_z_Upwind[2] = {0.0};;
  double Ghat_M0_z[2] = {0.0}; 
  double N_x_UpwindQuad[2] = {0.0};
  double N_x_Upwind[2] = {0.0};;
  double Ghat_N_x[2] = {0.0}; 
  double N_y_UpwindQuad[2] = {0.0};
  double N_y_Upwind[2] = {0.0};;
  double Ghat_N_y[2] = {0.0}; 
  double N_z_UpwindQuad[2] = {0.0};
  double N_z_Upwind[2] = {0.0};;
  double Ghat_N_z[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 

  if (0.7071067811865475*alphaDrSurf[0]-0.7071067811865475*alphaDrSurf[1] < 0) { 
    F_0_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(G_1Skin); 
    M0_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_xSkin); 
    M0_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_ySkin); 
    M0_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_zSkin); 
    N_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_xSkin); 
    N_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_ySkin); 
    N_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(G_1Edge); 
    M0_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_xEdge); 
    M0_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_yEdge); 
    M0_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_zEdge); 
    N_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_xEdge); 
    N_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_yEdge); 
    N_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_zEdge); 
  } 
  if (0.7071067811865475*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(G_1Skin); 
    M0_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_xSkin); 
    M0_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_ySkin); 
    M0_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_zSkin); 
    N_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_xSkin); 
    N_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_ySkin); 
    N_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_zSkin); 
  } else { 
    F_0_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(G_1Edge); 
    M0_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_xEdge); 
    M0_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_yEdge); 
    M0_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_zEdge); 
    N_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_xEdge); 
    N_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_yEdge); 
    N_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_zEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad, M0_x_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad, M0_y_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad, M0_z_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad, N_x_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad, N_y_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad, N_z_Upwind); 

  Ghat_F_0[0] = 0.7071067811865475*F_0_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_G_1[0] = 0.7071067811865475*G_1_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_M0_x[0] = 0.7071067811865475*M0_x_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_x[1] = 0.7071067811865475*M0_x_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_x_Upwind[1]; 
  Ghat_M0_y[0] = 0.7071067811865475*M0_y_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_y[1] = 0.7071067811865475*M0_y_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_y_Upwind[1]; 
  Ghat_M0_z[0] = 0.7071067811865475*M0_z_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_z[1] = 0.7071067811865475*M0_z_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_z_Upwind[1]; 
  Ghat_N_x[0] = 0.7071067811865475*N_x_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_x[1] = 0.7071067811865475*N_x_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_x_Upwind[1]; 
  Ghat_N_y[0] = 0.7071067811865475*N_y_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_y[1] = 0.7071067811865475*N_y_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_y_Upwind[1]; 
  Ghat_N_z[0] = 0.7071067811865475*N_z_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_z[1] = 0.7071067811865475*N_z_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_z_Upwind[1]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += 1.5811388300841895*Ghat_F_0[0]*dv1par; 
  out_F_0[5] += 1.5811388300841898*Ghat_F_0[1]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += 1.5811388300841895*Ghat_G_1[0]*dv1par; 
  out_G_1[5] += 1.5811388300841898*Ghat_G_1[1]*dv1par; 
  out_M0_x[0] += 0.7071067811865475*Ghat_M0_x[0]*dv1par; 
  out_M0_x[1] += 0.7071067811865475*Ghat_M0_x[1]*dv1par; 
  out_M0_x[2] += 1.224744871391589*Ghat_M0_x[0]*dv1par; 
  out_M0_x[3] += 1.224744871391589*Ghat_M0_x[1]*dv1par; 
  out_M0_x[4] += 1.5811388300841895*Ghat_M0_x[0]*dv1par; 
  out_M0_x[5] += 1.5811388300841898*Ghat_M0_x[1]*dv1par; 
  out_M0_y[0] += 0.7071067811865475*Ghat_M0_y[0]*dv1par; 
  out_M0_y[1] += 0.7071067811865475*Ghat_M0_y[1]*dv1par; 
  out_M0_y[2] += 1.224744871391589*Ghat_M0_y[0]*dv1par; 
  out_M0_y[3] += 1.224744871391589*Ghat_M0_y[1]*dv1par; 
  out_M0_y[4] += 1.5811388300841895*Ghat_M0_y[0]*dv1par; 
  out_M0_y[5] += 1.5811388300841898*Ghat_M0_y[1]*dv1par; 
  out_M0_z[0] += 0.7071067811865475*Ghat_M0_z[0]*dv1par; 
  out_M0_z[1] += 0.7071067811865475*Ghat_M0_z[1]*dv1par; 
  out_M0_z[2] += 1.224744871391589*Ghat_M0_z[0]*dv1par; 
  out_M0_z[3] += 1.224744871391589*Ghat_M0_z[1]*dv1par; 
  out_M0_z[4] += 1.5811388300841895*Ghat_M0_z[0]*dv1par; 
  out_M0_z[5] += 1.5811388300841898*Ghat_M0_z[1]*dv1par; 
  out_N_x[0] += 0.7071067811865475*Ghat_N_x[0]*dv1par; 
  out_N_x[1] += 0.7071067811865475*Ghat_N_x[1]*dv1par; 
  out_N_x[2] += 1.224744871391589*Ghat_N_x[0]*dv1par; 
  out_N_x[3] += 1.224744871391589*Ghat_N_x[1]*dv1par; 
  out_N_x[4] += 1.5811388300841895*Ghat_N_x[0]*dv1par; 
  out_N_x[5] += 1.5811388300841898*Ghat_N_x[1]*dv1par; 
  out_N_y[0] += 0.7071067811865475*Ghat_N_y[0]*dv1par; 
  out_N_y[1] += 0.7071067811865475*Ghat_N_y[1]*dv1par; 
  out_N_y[2] += 1.224744871391589*Ghat_N_y[0]*dv1par; 
  out_N_y[3] += 1.224744871391589*Ghat_N_y[1]*dv1par; 
  out_N_y[4] += 1.5811388300841895*Ghat_N_y[0]*dv1par; 
  out_N_y[5] += 1.5811388300841898*Ghat_N_y[1]*dv1par; 
  out_N_z[0] += 0.7071067811865475*Ghat_N_z[0]*dv1par; 
  out_N_z[1] += 0.7071067811865475*Ghat_N_z[1]*dv1par; 
  out_N_z[2] += 1.224744871391589*Ghat_N_z[0]*dv1par; 
  out_N_z[3] += 1.224744871391589*Ghat_N_z[1]*dv1par; 
  out_N_z[4] += 1.5811388300841895*Ghat_N_z[0]*dv1par; 
  out_N_z[5] += 1.5811388300841898*Ghat_N_z[1]*dv1par; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 

  if (0.7071067811865475*alphaDrSurf[0]-0.7071067811865475*alphaDrSurf[1] < 0) { 
    F_0_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(F_0Edge); 
    G_1_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(G_1Edge); 
    M0_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_xEdge); 
    M0_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_yEdge); 
    M0_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_zEdge); 
    N_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_xEdge); 
    N_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_yEdge); 
    N_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(F_0Skin); 
    G_1_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(G_1Skin); 
    M0_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_xSkin); 
    M0_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_ySkin); 
    M0_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_zSkin); 
    N_x_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_xSkin); 
    N_y_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_ySkin); 
    N_z_UpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_zSkin); 
  } 
  if (0.7071067811865475*(alphaDrSurf[1]+alphaDrSurf[0]) < 0) { 
    F_0_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(F_0Edge); 
    G_1_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(G_1Edge); 
    M0_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_xEdge); 
    M0_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_yEdge); 
    M0_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_zEdge); 
    N_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_xEdge); 
    N_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_yEdge); 
    N_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_zEdge); 
  } else { 
    F_0_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(F_0Skin); 
    G_1_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(G_1Skin); 
    M0_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_xSkin); 
    M0_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_ySkin); 
    M0_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_zSkin); 
    N_x_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_xSkin); 
    N_y_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_ySkin); 
    N_z_UpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_zSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad, F_0_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad, G_1_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad, M0_x_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad, M0_y_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad, M0_z_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad, N_x_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad, N_y_Upwind); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad, N_z_Upwind); 

  Ghat_F_0[0] = 0.7071067811865475*F_0_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.7071067811865475*F_0_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*F_0_Upwind[1]; 
  Ghat_G_1[0] = 0.7071067811865475*G_1_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.7071067811865475*G_1_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*G_1_Upwind[1]; 
  Ghat_M0_x[0] = 0.7071067811865475*M0_x_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_x[1] = 0.7071067811865475*M0_x_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_x_Upwind[1]; 
  Ghat_M0_y[0] = 0.7071067811865475*M0_y_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_y[1] = 0.7071067811865475*M0_y_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_y_Upwind[1]; 
  Ghat_M0_z[0] = 0.7071067811865475*M0_z_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*M0_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_M0_z[1] = 0.7071067811865475*M0_z_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*M0_z_Upwind[1]; 
  Ghat_N_x[0] = 0.7071067811865475*N_x_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_x_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_x[1] = 0.7071067811865475*N_x_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_x_Upwind[1]; 
  Ghat_N_y[0] = 0.7071067811865475*N_y_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_y_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_y[1] = 0.7071067811865475*N_y_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_y_Upwind[1]; 
  Ghat_N_z[0] = 0.7071067811865475*N_z_Upwind[1]*alphaDrSurf[1]+0.7071067811865475*N_z_Upwind[0]*alphaDrSurf[0]; 
  Ghat_N_z[1] = 0.7071067811865475*N_z_Upwind[0]*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0]*N_z_Upwind[1]; 

  out_F_0[0] += -(0.7071067811865475*Ghat_F_0[0]*dv1par); 
  out_F_0[1] += -(0.7071067811865475*Ghat_F_0[1]*dv1par); 
  out_F_0[2] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[4] += -(1.5811388300841895*Ghat_F_0[0]*dv1par); 
  out_F_0[5] += -(1.5811388300841898*Ghat_F_0[1]*dv1par); 
  out_G_1[0] += -(0.7071067811865475*Ghat_G_1[0]*dv1par); 
  out_G_1[1] += -(0.7071067811865475*Ghat_G_1[1]*dv1par); 
  out_G_1[2] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[4] += -(1.5811388300841895*Ghat_G_1[0]*dv1par); 
  out_G_1[5] += -(1.5811388300841898*Ghat_G_1[1]*dv1par); 
  out_M0_x[0] += -(0.7071067811865475*Ghat_M0_x[0]*dv1par); 
  out_M0_x[1] += -(0.7071067811865475*Ghat_M0_x[1]*dv1par); 
  out_M0_x[2] += 1.224744871391589*Ghat_M0_x[0]*dv1par; 
  out_M0_x[3] += 1.224744871391589*Ghat_M0_x[1]*dv1par; 
  out_M0_x[4] += -(1.5811388300841895*Ghat_M0_x[0]*dv1par); 
  out_M0_x[5] += -(1.5811388300841898*Ghat_M0_x[1]*dv1par); 
  out_M0_y[0] += -(0.7071067811865475*Ghat_M0_y[0]*dv1par); 
  out_M0_y[1] += -(0.7071067811865475*Ghat_M0_y[1]*dv1par); 
  out_M0_y[2] += 1.224744871391589*Ghat_M0_y[0]*dv1par; 
  out_M0_y[3] += 1.224744871391589*Ghat_M0_y[1]*dv1par; 
  out_M0_y[4] += -(1.5811388300841895*Ghat_M0_y[0]*dv1par); 
  out_M0_y[5] += -(1.5811388300841898*Ghat_M0_y[1]*dv1par); 
  out_M0_z[0] += -(0.7071067811865475*Ghat_M0_z[0]*dv1par); 
  out_M0_z[1] += -(0.7071067811865475*Ghat_M0_z[1]*dv1par); 
  out_M0_z[2] += 1.224744871391589*Ghat_M0_z[0]*dv1par; 
  out_M0_z[3] += 1.224744871391589*Ghat_M0_z[1]*dv1par; 
  out_M0_z[4] += -(1.5811388300841895*Ghat_M0_z[0]*dv1par); 
  out_M0_z[5] += -(1.5811388300841898*Ghat_M0_z[1]*dv1par); 
  out_N_x[0] += -(0.7071067811865475*Ghat_N_x[0]*dv1par); 
  out_N_x[1] += -(0.7071067811865475*Ghat_N_x[1]*dv1par); 
  out_N_x[2] += 1.224744871391589*Ghat_N_x[0]*dv1par; 
  out_N_x[3] += 1.224744871391589*Ghat_N_x[1]*dv1par; 
  out_N_x[4] += -(1.5811388300841895*Ghat_N_x[0]*dv1par); 
  out_N_x[5] += -(1.5811388300841898*Ghat_N_x[1]*dv1par); 
  out_N_y[0] += -(0.7071067811865475*Ghat_N_y[0]*dv1par); 
  out_N_y[1] += -(0.7071067811865475*Ghat_N_y[1]*dv1par); 
  out_N_y[2] += 1.224744871391589*Ghat_N_y[0]*dv1par; 
  out_N_y[3] += 1.224744871391589*Ghat_N_y[1]*dv1par; 
  out_N_y[4] += -(1.5811388300841895*Ghat_N_y[0]*dv1par); 
  out_N_y[5] += -(1.5811388300841898*Ghat_N_y[1]*dv1par); 
  out_N_z[0] += -(0.7071067811865475*Ghat_N_z[0]*dv1par); 
  out_N_z[1] += -(0.7071067811865475*Ghat_N_z[1]*dv1par); 
  out_N_z[2] += 1.224744871391589*Ghat_N_z[0]*dv1par; 
  out_N_z[3] += 1.224744871391589*Ghat_N_z[1]*dv1par; 
  out_N_z[4] += -(1.5811388300841895*Ghat_N_z[0]*dv1par); 
  out_N_z[5] += -(1.5811388300841898*Ghat_N_z[1]*dv1par); 

  }

  return 0.;

} 
