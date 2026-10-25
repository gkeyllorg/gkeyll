#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1), M0_perp, N_perp] in left/center/right cells. 
  // out:           Incremented output distribution functions in center cell. 

  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  const double *G_1l = &fl[6]; 
  const double *G_1c = &fc[6]; 
  const double *G_1r = &fr[6]; 
  const double *M0_xl = &fl[12]; 
  const double *M0_xc = &fc[12]; 
  const double *M0_xr = &fr[12]; 
  const double *M0_yl = &fl[18]; 
  const double *M0_yc = &fc[18]; 
  const double *M0_yr = &fr[18]; 
  const double *M0_zl = &fl[24]; 
  const double *M0_zc = &fc[24]; 
  const double *M0_zr = &fr[24]; 
  const double *N_xl = &fl[30]; 
  const double *N_xc = &fc[30]; 
  const double *N_xr = &fr[30]; 
  const double *N_yl = &fl[36]; 
  const double *N_yc = &fc[36]; 
  const double *N_yr = &fr[36]; 
  const double *N_zl = &fl[42]; 
  const double *N_zc = &fc[42]; 
  const double *N_zr = &fr[42]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double *out_M0_x = &out[12]; 
  double *out_M0_y = &out[18]; 
  double *out_M0_z = &out[24]; 
  double *out_N_x = &out[30]; 
  double *out_N_y = &out[36]; 
  double *out_N_z = &out[42]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_l[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_r[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 

  double F_0_UpwindQuad_l[2] = {0.0};
  double F_0_UpwindQuad_r[2] = {0.0};
  double F_0_Upwind_l[2] = {0.0};
  double F_0_Upwind_r[2] = {0.0};
  double Ghat_F_0_l[2] = {0.0}; 
  double Ghat_F_0_r[2] = {0.0}; 

  double G_1_UpwindQuad_l[2] = {0.0};
  double G_1_UpwindQuad_r[2] = {0.0};
  double G_1_Upwind_l[2] = {0.0};
  double G_1_Upwind_r[2] = {0.0};
  double Ghat_G_1_l[2] = {0.0}; 
  double Ghat_G_1_r[2] = {0.0}; 

  double M0_x_UpwindQuad_l[2] = {0.0};
  double M0_x_UpwindQuad_r[2] = {0.0};
  double M0_x_Upwind_l[2] = {0.0};
  double M0_x_Upwind_r[2] = {0.0};
  double Ghat_M0_x_l[2] = {0.0}; 
  double Ghat_M0_x_r[2] = {0.0}; 

  double M0_y_UpwindQuad_l[2] = {0.0};
  double M0_y_UpwindQuad_r[2] = {0.0};
  double M0_y_Upwind_l[2] = {0.0};
  double M0_y_Upwind_r[2] = {0.0};
  double Ghat_M0_y_l[2] = {0.0}; 
  double Ghat_M0_y_r[2] = {0.0}; 

  double M0_z_UpwindQuad_l[2] = {0.0};
  double M0_z_UpwindQuad_r[2] = {0.0};
  double M0_z_Upwind_l[2] = {0.0};
  double M0_z_Upwind_r[2] = {0.0};
  double Ghat_M0_z_l[2] = {0.0}; 
  double Ghat_M0_z_r[2] = {0.0}; 

  double N_x_UpwindQuad_l[2] = {0.0};
  double N_x_UpwindQuad_r[2] = {0.0};
  double N_x_Upwind_l[2] = {0.0};
  double N_x_Upwind_r[2] = {0.0};
  double Ghat_N_x_l[2] = {0.0}; 
  double Ghat_N_x_r[2] = {0.0}; 

  double N_y_UpwindQuad_l[2] = {0.0};
  double N_y_UpwindQuad_r[2] = {0.0};
  double N_y_Upwind_l[2] = {0.0};
  double N_y_Upwind_r[2] = {0.0};
  double Ghat_N_y_l[2] = {0.0}; 
  double Ghat_N_y_r[2] = {0.0}; 

  double N_z_UpwindQuad_l[2] = {0.0};
  double N_z_UpwindQuad_r[2] = {0.0};
  double N_z_Upwind_l[2] = {0.0};
  double N_z_Upwind_r[2] = {0.0};
  double Ghat_N_z_l[2] = {0.0}; 
  double Ghat_N_z_r[2] = {0.0}; 

  if (0.7071067811865475*alphaDrSurf_l[0]-0.7071067811865475*alphaDrSurf_l[1] < 0) { 
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(G_1l); 
    M0_x_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_xl); 
    M0_y_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_yl); 
    M0_z_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_zl); 
    N_x_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_xl); 
    N_y_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_yl); 
    N_z_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(G_1c); 
    M0_x_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_xc); 
    M0_y_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_yc); 
    M0_z_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_zc); 
    N_x_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_xc); 
    N_y_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_yc); 
    N_z_UpwindQuad_l[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_zc); 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7071067811865475*alphaDrSurf_r[1] < 0) { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(G_1c); 
    M0_x_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_xc); 
    M0_y_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_yc); 
    M0_z_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(M0_zc); 
    N_x_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_xc); 
    N_y_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_yc); 
    N_z_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(G_1r); 
    M0_x_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_xr); 
    M0_y_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_yr); 
    M0_z_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(M0_zr); 
    N_x_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_xr); 
    N_y_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_yr); 
    N_z_UpwindQuad_r[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(N_zr); 
  } 
  if (0.7071067811865475*(alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(G_1l); 
    M0_x_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_xl); 
    M0_y_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_yl); 
    M0_z_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_zl); 
    N_x_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_xl); 
    N_y_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_yl); 
    N_z_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(G_1c); 
    M0_x_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_xc); 
    M0_y_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_yc); 
    M0_z_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_zc); 
    N_x_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_xc); 
    N_y_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_yc); 
    N_z_UpwindQuad_l[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_zc); 
  } 
  if (0.7071067811865475*(alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(G_1c); 
    M0_x_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_xc); 
    M0_y_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_yc); 
    M0_z_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(M0_zc); 
    N_x_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_xc); 
    N_y_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_yc); 
    N_z_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(G_1r); 
    M0_x_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_xr); 
    M0_y_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_yr); 
    M0_z_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(M0_zr); 
    N_x_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_xr); 
    N_y_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_yr); 
    N_z_UpwindQuad_r[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(N_zr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad_l, M0_x_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad_r, M0_x_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad_l, M0_y_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad_r, M0_y_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad_l, M0_z_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad_r, M0_z_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad_l, N_x_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad_r, N_x_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad_l, N_y_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad_r, N_y_Upwind_r); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad_l, N_z_Upwind_l); 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad_r, N_z_Upwind_r); 

  Ghat_F_0_l[0] = 0.7071067811865475*(F_0_Upwind_l[1]*alphaDrSurf_l[1]+F_0_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.7071067811865475*(F_0_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0_Upwind_l[1]); 
  Ghat_F_0_r[0] = 0.7071067811865475*(F_0_Upwind_r[1]*alphaDrSurf_r[1]+F_0_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.7071067811865475*(F_0_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0_Upwind_r[1]); 

  Ghat_G_1_l[0] = 0.7071067811865475*(G_1_Upwind_l[1]*alphaDrSurf_l[1]+G_1_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.7071067811865475*(G_1_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1_Upwind_l[1]); 
  Ghat_G_1_r[0] = 0.7071067811865475*(G_1_Upwind_r[1]*alphaDrSurf_r[1]+G_1_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.7071067811865475*(G_1_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1_Upwind_r[1]); 

  Ghat_M0_x_l[0] = 0.7071067811865475*(M0_x_Upwind_l[1]*alphaDrSurf_l[1]+M0_x_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_x_l[1] = 0.7071067811865475*(M0_x_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_x_Upwind_l[1]); 
  Ghat_M0_x_r[0] = 0.7071067811865475*(M0_x_Upwind_r[1]*alphaDrSurf_r[1]+M0_x_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_x_r[1] = 0.7071067811865475*(M0_x_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_x_Upwind_r[1]); 

  Ghat_M0_y_l[0] = 0.7071067811865475*(M0_y_Upwind_l[1]*alphaDrSurf_l[1]+M0_y_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_y_l[1] = 0.7071067811865475*(M0_y_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_y_Upwind_l[1]); 
  Ghat_M0_y_r[0] = 0.7071067811865475*(M0_y_Upwind_r[1]*alphaDrSurf_r[1]+M0_y_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_y_r[1] = 0.7071067811865475*(M0_y_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_y_Upwind_r[1]); 

  Ghat_M0_z_l[0] = 0.7071067811865475*(M0_z_Upwind_l[1]*alphaDrSurf_l[1]+M0_z_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_z_l[1] = 0.7071067811865475*(M0_z_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_z_Upwind_l[1]); 
  Ghat_M0_z_r[0] = 0.7071067811865475*(M0_z_Upwind_r[1]*alphaDrSurf_r[1]+M0_z_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_z_r[1] = 0.7071067811865475*(M0_z_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_z_Upwind_r[1]); 

  Ghat_N_x_l[0] = 0.7071067811865475*(N_x_Upwind_l[1]*alphaDrSurf_l[1]+N_x_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_x_l[1] = 0.7071067811865475*(N_x_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_x_Upwind_l[1]); 
  Ghat_N_x_r[0] = 0.7071067811865475*(N_x_Upwind_r[1]*alphaDrSurf_r[1]+N_x_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_x_r[1] = 0.7071067811865475*(N_x_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_x_Upwind_r[1]); 

  Ghat_N_y_l[0] = 0.7071067811865475*(N_y_Upwind_l[1]*alphaDrSurf_l[1]+N_y_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_y_l[1] = 0.7071067811865475*(N_y_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_y_Upwind_l[1]); 
  Ghat_N_y_r[0] = 0.7071067811865475*(N_y_Upwind_r[1]*alphaDrSurf_r[1]+N_y_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_y_r[1] = 0.7071067811865475*(N_y_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_y_Upwind_r[1]); 

  Ghat_N_z_l[0] = 0.7071067811865475*(N_z_Upwind_l[1]*alphaDrSurf_l[1]+N_z_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_z_l[1] = 0.7071067811865475*(N_z_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_z_Upwind_l[1]); 
  Ghat_N_z_r[0] = 0.7071067811865475*(N_z_Upwind_r[1]*alphaDrSurf_r[1]+N_z_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_z_r[1] = 0.7071067811865475*(N_z_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_z_Upwind_r[1]); 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[3] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[4] += (1.5811388300841895*Ghat_F_0_r[0]-1.5811388300841895*Ghat_F_0_l[0])*dv1par; 
  out_F_0[5] += (1.5811388300841898*Ghat_F_0_r[1]-1.5811388300841898*Ghat_F_0_l[1])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[3] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[4] += (1.5811388300841895*Ghat_G_1_r[0]-1.5811388300841895*Ghat_G_1_l[0])*dv1par; 
  out_G_1[5] += (1.5811388300841898*Ghat_G_1_r[1]-1.5811388300841898*Ghat_G_1_l[1])*dv1par; 
  out_M0_x[0] += (0.7071067811865475*Ghat_M0_x_r[0]-0.7071067811865475*Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[1] += (0.7071067811865475*Ghat_M0_x_r[1]-0.7071067811865475*Ghat_M0_x_l[1])*dv1par; 
  out_M0_x[2] += 1.224744871391589*(Ghat_M0_x_r[0]+Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[3] += 1.224744871391589*(Ghat_M0_x_r[1]+Ghat_M0_x_l[1])*dv1par; 
  out_M0_x[4] += (1.5811388300841895*Ghat_M0_x_r[0]-1.5811388300841895*Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[5] += (1.5811388300841898*Ghat_M0_x_r[1]-1.5811388300841898*Ghat_M0_x_l[1])*dv1par; 
  out_M0_y[0] += (0.7071067811865475*Ghat_M0_y_r[0]-0.7071067811865475*Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[1] += (0.7071067811865475*Ghat_M0_y_r[1]-0.7071067811865475*Ghat_M0_y_l[1])*dv1par; 
  out_M0_y[2] += 1.224744871391589*(Ghat_M0_y_r[0]+Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[3] += 1.224744871391589*(Ghat_M0_y_r[1]+Ghat_M0_y_l[1])*dv1par; 
  out_M0_y[4] += (1.5811388300841895*Ghat_M0_y_r[0]-1.5811388300841895*Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[5] += (1.5811388300841898*Ghat_M0_y_r[1]-1.5811388300841898*Ghat_M0_y_l[1])*dv1par; 
  out_M0_z[0] += (0.7071067811865475*Ghat_M0_z_r[0]-0.7071067811865475*Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[1] += (0.7071067811865475*Ghat_M0_z_r[1]-0.7071067811865475*Ghat_M0_z_l[1])*dv1par; 
  out_M0_z[2] += 1.224744871391589*(Ghat_M0_z_r[0]+Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[3] += 1.224744871391589*(Ghat_M0_z_r[1]+Ghat_M0_z_l[1])*dv1par; 
  out_M0_z[4] += (1.5811388300841895*Ghat_M0_z_r[0]-1.5811388300841895*Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[5] += (1.5811388300841898*Ghat_M0_z_r[1]-1.5811388300841898*Ghat_M0_z_l[1])*dv1par; 
  out_N_x[0] += (0.7071067811865475*Ghat_N_x_r[0]-0.7071067811865475*Ghat_N_x_l[0])*dv1par; 
  out_N_x[1] += (0.7071067811865475*Ghat_N_x_r[1]-0.7071067811865475*Ghat_N_x_l[1])*dv1par; 
  out_N_x[2] += 1.224744871391589*(Ghat_N_x_r[0]+Ghat_N_x_l[0])*dv1par; 
  out_N_x[3] += 1.224744871391589*(Ghat_N_x_r[1]+Ghat_N_x_l[1])*dv1par; 
  out_N_x[4] += (1.5811388300841895*Ghat_N_x_r[0]-1.5811388300841895*Ghat_N_x_l[0])*dv1par; 
  out_N_x[5] += (1.5811388300841898*Ghat_N_x_r[1]-1.5811388300841898*Ghat_N_x_l[1])*dv1par; 
  out_N_y[0] += (0.7071067811865475*Ghat_N_y_r[0]-0.7071067811865475*Ghat_N_y_l[0])*dv1par; 
  out_N_y[1] += (0.7071067811865475*Ghat_N_y_r[1]-0.7071067811865475*Ghat_N_y_l[1])*dv1par; 
  out_N_y[2] += 1.224744871391589*(Ghat_N_y_r[0]+Ghat_N_y_l[0])*dv1par; 
  out_N_y[3] += 1.224744871391589*(Ghat_N_y_r[1]+Ghat_N_y_l[1])*dv1par; 
  out_N_y[4] += (1.5811388300841895*Ghat_N_y_r[0]-1.5811388300841895*Ghat_N_y_l[0])*dv1par; 
  out_N_y[5] += (1.5811388300841898*Ghat_N_y_r[1]-1.5811388300841898*Ghat_N_y_l[1])*dv1par; 
  out_N_z[0] += (0.7071067811865475*Ghat_N_z_r[0]-0.7071067811865475*Ghat_N_z_l[0])*dv1par; 
  out_N_z[1] += (0.7071067811865475*Ghat_N_z_r[1]-0.7071067811865475*Ghat_N_z_l[1])*dv1par; 
  out_N_z[2] += 1.224744871391589*(Ghat_N_z_r[0]+Ghat_N_z_l[0])*dv1par; 
  out_N_z[3] += 1.224744871391589*(Ghat_N_z_r[1]+Ghat_N_z_l[1])*dv1par; 
  out_N_z[4] += (1.5811388300841895*Ghat_N_z_r[0]-1.5811388300841895*Ghat_N_z_l[0])*dv1par; 
  out_N_z[5] += (1.5811388300841898*Ghat_N_z_r[1]-1.5811388300841898*Ghat_N_z_l[1])*dv1par; 

  return 0.;

} 
