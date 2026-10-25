#include <gkyl_lbo_pkpm_kernels.h>  
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_pkpm_drag_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr:      Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1), M0_perp, N_perp] in left/center/right cells. 
  // out:           Incremented output distribution functions in center cell. 

  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  const double *G_1l = &fl[12]; 
  const double *G_1c = &fc[12]; 
  const double *G_1r = &fr[12]; 
  const double *M0_xl = &fl[24]; 
  const double *M0_xc = &fc[24]; 
  const double *M0_xr = &fr[24]; 
  const double *M0_yl = &fl[36]; 
  const double *M0_yc = &fc[36]; 
  const double *M0_yr = &fr[36]; 
  const double *M0_zl = &fl[48]; 
  const double *M0_zc = &fc[48]; 
  const double *M0_zr = &fr[48]; 
  const double *N_xl = &fl[60]; 
  const double *N_xc = &fc[60]; 
  const double *N_xr = &fr[60]; 
  const double *N_yl = &fl[72]; 
  const double *N_yc = &fc[72]; 
  const double *N_yr = &fr[72]; 
  const double *N_zl = &fl[84]; 
  const double *N_zc = &fc[84]; 
  const double *N_zr = &fr[84]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  double *out_M0_x = &out[24]; 
  double *out_M0_y = &out[36]; 
  double *out_M0_z = &out[48]; 
  double *out_N_x = &out[60]; 
  double *out_N_y = &out[72]; 
  double *out_N_z = &out[84]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = nuSum[0]*wvpar-0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_l[1] = nuSum[1]*wvpar-0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_l[2] = nuSum[2]*wvpar-0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf_l[3] = nuSum[3]*wvpar-0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = nuSum[0]*wvpar+0.5*nuSum[0]*dvpar-1.0*sumNuUPar[0]; 
  alphaDrSurf_r[1] = nuSum[1]*wvpar+0.5*nuSum[1]*dvpar-1.0*sumNuUPar[1]; 
  alphaDrSurf_r[2] = nuSum[2]*wvpar+0.5*nuSum[2]*dvpar-1.0*sumNuUPar[2]; 
  alphaDrSurf_r[3] = nuSum[3]*wvpar+0.5*nuSum[3]*dvpar-1.0*sumNuUPar[3]; 

  double F_0_UpwindQuad_l[4] = {0.0};
  double F_0_UpwindQuad_r[4] = {0.0};
  double F_0_Upwind_l[4] = {0.0};
  double F_0_Upwind_r[4] = {0.0};
  double Ghat_F_0_l[4] = {0.0}; 
  double Ghat_F_0_r[4] = {0.0}; 

  double G_1_UpwindQuad_l[4] = {0.0};
  double G_1_UpwindQuad_r[4] = {0.0};
  double G_1_Upwind_l[4] = {0.0};
  double G_1_Upwind_r[4] = {0.0};
  double Ghat_G_1_l[4] = {0.0}; 
  double Ghat_G_1_r[4] = {0.0}; 

  double M0_x_UpwindQuad_l[4] = {0.0};
  double M0_x_UpwindQuad_r[4] = {0.0};
  double M0_x_Upwind_l[4] = {0.0};
  double M0_x_Upwind_r[4] = {0.0};
  double Ghat_M0_x_l[4] = {0.0}; 
  double Ghat_M0_x_r[4] = {0.0}; 

  double M0_y_UpwindQuad_l[4] = {0.0};
  double M0_y_UpwindQuad_r[4] = {0.0};
  double M0_y_Upwind_l[4] = {0.0};
  double M0_y_Upwind_r[4] = {0.0};
  double Ghat_M0_y_l[4] = {0.0}; 
  double Ghat_M0_y_r[4] = {0.0}; 

  double M0_z_UpwindQuad_l[4] = {0.0};
  double M0_z_UpwindQuad_r[4] = {0.0};
  double M0_z_Upwind_l[4] = {0.0};
  double M0_z_Upwind_r[4] = {0.0};
  double Ghat_M0_z_l[4] = {0.0}; 
  double Ghat_M0_z_r[4] = {0.0}; 

  double N_x_UpwindQuad_l[4] = {0.0};
  double N_x_UpwindQuad_r[4] = {0.0};
  double N_x_Upwind_l[4] = {0.0};
  double N_x_Upwind_r[4] = {0.0};
  double Ghat_N_x_l[4] = {0.0}; 
  double Ghat_N_x_r[4] = {0.0}; 

  double N_y_UpwindQuad_l[4] = {0.0};
  double N_y_UpwindQuad_r[4] = {0.0};
  double N_y_Upwind_l[4] = {0.0};
  double N_y_Upwind_r[4] = {0.0};
  double Ghat_N_y_l[4] = {0.0}; 
  double Ghat_N_y_r[4] = {0.0}; 

  double N_z_UpwindQuad_l[4] = {0.0};
  double N_z_UpwindQuad_r[4] = {0.0};
  double N_z_Upwind_l[4] = {0.0};
  double N_z_Upwind_r[4] = {0.0};
  double Ghat_N_z_l[4] = {0.0}; 
  double Ghat_N_z_r[4] = {0.0}; 

  if (0.5*alphaDrSurf_l[3]-0.5*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.5*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1l); 
    M0_x_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_xl); 
    M0_y_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_yl); 
    M0_z_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_zl); 
    N_x_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_xl); 
    N_y_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_yl); 
    N_z_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1c); 
    M0_x_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_xc); 
    M0_y_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_yc); 
    M0_z_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_zc); 
    N_x_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_xc); 
    N_y_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_yc); 
    N_z_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_zc); 
  } 
  if (0.5*alphaDrSurf_r[3]-0.5*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.5*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1c); 
    M0_x_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_xc); 
    M0_y_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_yc); 
    M0_z_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(M0_zc); 
    N_x_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_xc); 
    N_y_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_yc); 
    N_z_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1r); 
    M0_x_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_xr); 
    M0_y_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_yr); 
    M0_z_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(M0_zr); 
    N_x_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_xr); 
    N_y_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_yr); 
    N_z_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(N_zr); 
  } 
  if (-(0.5*alphaDrSurf_l[3])+0.5*alphaDrSurf_l[2]-0.5*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1l); 
    M0_x_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_xl); 
    M0_y_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_yl); 
    M0_z_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_zl); 
    N_x_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_xl); 
    N_y_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_yl); 
    N_z_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1c); 
    M0_x_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_xc); 
    M0_y_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_yc); 
    M0_z_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_zc); 
    N_x_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_xc); 
    N_y_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_yc); 
    N_z_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_zc); 
  } 
  if (-(0.5*alphaDrSurf_r[3])+0.5*alphaDrSurf_r[2]-0.5*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1c); 
    M0_x_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_xc); 
    M0_y_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_yc); 
    M0_z_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(M0_zc); 
    N_x_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_xc); 
    N_y_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_yc); 
    N_z_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1r); 
    M0_x_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_xr); 
    M0_y_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_yr); 
    M0_z_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(M0_zr); 
    N_x_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_xr); 
    N_y_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_yr); 
    N_z_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(N_zr); 
  } 
  if (0.5*(alphaDrSurf_l[1]+alphaDrSurf_l[0])-0.5*(alphaDrSurf_l[3]+alphaDrSurf_l[2]) < 0) { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1l); 
    M0_x_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_xl); 
    M0_y_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_yl); 
    M0_z_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_zl); 
    N_x_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_xl); 
    N_y_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_yl); 
    N_z_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1c); 
    M0_x_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_xc); 
    M0_y_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_yc); 
    M0_z_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_zc); 
    N_x_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_xc); 
    N_y_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_yc); 
    N_z_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_zc); 
  } 
  if (0.5*(alphaDrSurf_r[1]+alphaDrSurf_r[0])-0.5*(alphaDrSurf_r[3]+alphaDrSurf_r[2]) < 0) { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1c); 
    M0_x_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_xc); 
    M0_y_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_yc); 
    M0_z_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(M0_zc); 
    N_x_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_xc); 
    N_y_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_yc); 
    N_z_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1r); 
    M0_x_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_xr); 
    M0_y_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_yr); 
    M0_z_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(M0_zr); 
    N_x_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_xr); 
    N_y_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_yr); 
    N_z_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(N_zr); 
  } 
  if (0.5*(alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0]) < 0) { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1l); 
    M0_x_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_xl); 
    M0_y_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_yl); 
    M0_z_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_zl); 
    N_x_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_xl); 
    N_y_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_yl); 
    N_z_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_zl); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1c); 
    M0_x_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_xc); 
    M0_y_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_yc); 
    M0_z_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_zc); 
    N_x_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_xc); 
    N_y_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_yc); 
    N_z_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_zc); 
  } 
  if (0.5*(alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0]) < 0) { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1c); 
    M0_x_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_xc); 
    M0_y_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_yc); 
    M0_z_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(M0_zc); 
    N_x_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_xc); 
    N_y_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_yc); 
    N_z_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(N_zc); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1r); 
    M0_x_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_xr); 
    M0_y_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_yr); 
    M0_z_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(M0_zr); 
    N_x_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_xr); 
    N_y_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_yr); 
    N_z_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(N_zr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad_l, M0_x_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_x_UpwindQuad_r, M0_x_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad_l, M0_y_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_y_UpwindQuad_r, M0_y_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad_l, M0_z_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(M0_z_UpwindQuad_r, M0_z_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad_l, N_x_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_x_UpwindQuad_r, N_x_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad_l, N_y_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_y_UpwindQuad_r, N_y_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad_l, N_z_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(N_z_UpwindQuad_r, N_z_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*(F_0_Upwind_l[3]*alphaDrSurf_l[3]+F_0_Upwind_l[2]*alphaDrSurf_l[2]+F_0_Upwind_l[1]*alphaDrSurf_l[1]+F_0_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.5*(F_0_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*F_0_Upwind_l[3]+F_0_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0_Upwind_l[1]); 
  Ghat_F_0_l[2] = 0.5*(F_0_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*F_0_Upwind_l[3]+F_0_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*F_0_Upwind_l[2]); 
  Ghat_F_0_l[3] = 0.5*(F_0_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*F_0_Upwind_l[3]+F_0_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*F_0_Upwind_l[2]); 
  Ghat_F_0_r[0] = 0.5*(F_0_Upwind_r[3]*alphaDrSurf_r[3]+F_0_Upwind_r[2]*alphaDrSurf_r[2]+F_0_Upwind_r[1]*alphaDrSurf_r[1]+F_0_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.5*(F_0_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*F_0_Upwind_r[3]+F_0_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0_Upwind_r[1]); 
  Ghat_F_0_r[2] = 0.5*(F_0_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*F_0_Upwind_r[3]+F_0_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*F_0_Upwind_r[2]); 
  Ghat_F_0_r[3] = 0.5*(F_0_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*F_0_Upwind_r[3]+F_0_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*F_0_Upwind_r[2]); 

  Ghat_G_1_l[0] = 0.5*(G_1_Upwind_l[3]*alphaDrSurf_l[3]+G_1_Upwind_l[2]*alphaDrSurf_l[2]+G_1_Upwind_l[1]*alphaDrSurf_l[1]+G_1_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.5*(G_1_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*G_1_Upwind_l[3]+G_1_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1_Upwind_l[1]); 
  Ghat_G_1_l[2] = 0.5*(G_1_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*G_1_Upwind_l[3]+G_1_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*G_1_Upwind_l[2]); 
  Ghat_G_1_l[3] = 0.5*(G_1_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*G_1_Upwind_l[3]+G_1_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*G_1_Upwind_l[2]); 
  Ghat_G_1_r[0] = 0.5*(G_1_Upwind_r[3]*alphaDrSurf_r[3]+G_1_Upwind_r[2]*alphaDrSurf_r[2]+G_1_Upwind_r[1]*alphaDrSurf_r[1]+G_1_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.5*(G_1_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*G_1_Upwind_r[3]+G_1_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1_Upwind_r[1]); 
  Ghat_G_1_r[2] = 0.5*(G_1_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*G_1_Upwind_r[3]+G_1_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*G_1_Upwind_r[2]); 
  Ghat_G_1_r[3] = 0.5*(G_1_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*G_1_Upwind_r[3]+G_1_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*G_1_Upwind_r[2]); 

  Ghat_M0_x_l[0] = 0.5*(M0_x_Upwind_l[3]*alphaDrSurf_l[3]+M0_x_Upwind_l[2]*alphaDrSurf_l[2]+M0_x_Upwind_l[1]*alphaDrSurf_l[1]+M0_x_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_x_l[1] = 0.5*(M0_x_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*M0_x_Upwind_l[3]+M0_x_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_x_Upwind_l[1]); 
  Ghat_M0_x_l[2] = 0.5*(M0_x_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*M0_x_Upwind_l[3]+M0_x_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*M0_x_Upwind_l[2]); 
  Ghat_M0_x_l[3] = 0.5*(M0_x_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*M0_x_Upwind_l[3]+M0_x_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*M0_x_Upwind_l[2]); 
  Ghat_M0_x_r[0] = 0.5*(M0_x_Upwind_r[3]*alphaDrSurf_r[3]+M0_x_Upwind_r[2]*alphaDrSurf_r[2]+M0_x_Upwind_r[1]*alphaDrSurf_r[1]+M0_x_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_x_r[1] = 0.5*(M0_x_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*M0_x_Upwind_r[3]+M0_x_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_x_Upwind_r[1]); 
  Ghat_M0_x_r[2] = 0.5*(M0_x_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*M0_x_Upwind_r[3]+M0_x_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*M0_x_Upwind_r[2]); 
  Ghat_M0_x_r[3] = 0.5*(M0_x_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*M0_x_Upwind_r[3]+M0_x_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*M0_x_Upwind_r[2]); 

  Ghat_M0_y_l[0] = 0.5*(M0_y_Upwind_l[3]*alphaDrSurf_l[3]+M0_y_Upwind_l[2]*alphaDrSurf_l[2]+M0_y_Upwind_l[1]*alphaDrSurf_l[1]+M0_y_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_y_l[1] = 0.5*(M0_y_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*M0_y_Upwind_l[3]+M0_y_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_y_Upwind_l[1]); 
  Ghat_M0_y_l[2] = 0.5*(M0_y_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*M0_y_Upwind_l[3]+M0_y_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*M0_y_Upwind_l[2]); 
  Ghat_M0_y_l[3] = 0.5*(M0_y_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*M0_y_Upwind_l[3]+M0_y_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*M0_y_Upwind_l[2]); 
  Ghat_M0_y_r[0] = 0.5*(M0_y_Upwind_r[3]*alphaDrSurf_r[3]+M0_y_Upwind_r[2]*alphaDrSurf_r[2]+M0_y_Upwind_r[1]*alphaDrSurf_r[1]+M0_y_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_y_r[1] = 0.5*(M0_y_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*M0_y_Upwind_r[3]+M0_y_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_y_Upwind_r[1]); 
  Ghat_M0_y_r[2] = 0.5*(M0_y_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*M0_y_Upwind_r[3]+M0_y_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*M0_y_Upwind_r[2]); 
  Ghat_M0_y_r[3] = 0.5*(M0_y_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*M0_y_Upwind_r[3]+M0_y_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*M0_y_Upwind_r[2]); 

  Ghat_M0_z_l[0] = 0.5*(M0_z_Upwind_l[3]*alphaDrSurf_l[3]+M0_z_Upwind_l[2]*alphaDrSurf_l[2]+M0_z_Upwind_l[1]*alphaDrSurf_l[1]+M0_z_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_M0_z_l[1] = 0.5*(M0_z_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*M0_z_Upwind_l[3]+M0_z_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*M0_z_Upwind_l[1]); 
  Ghat_M0_z_l[2] = 0.5*(M0_z_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*M0_z_Upwind_l[3]+M0_z_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*M0_z_Upwind_l[2]); 
  Ghat_M0_z_l[3] = 0.5*(M0_z_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*M0_z_Upwind_l[3]+M0_z_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*M0_z_Upwind_l[2]); 
  Ghat_M0_z_r[0] = 0.5*(M0_z_Upwind_r[3]*alphaDrSurf_r[3]+M0_z_Upwind_r[2]*alphaDrSurf_r[2]+M0_z_Upwind_r[1]*alphaDrSurf_r[1]+M0_z_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_M0_z_r[1] = 0.5*(M0_z_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*M0_z_Upwind_r[3]+M0_z_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*M0_z_Upwind_r[1]); 
  Ghat_M0_z_r[2] = 0.5*(M0_z_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*M0_z_Upwind_r[3]+M0_z_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*M0_z_Upwind_r[2]); 
  Ghat_M0_z_r[3] = 0.5*(M0_z_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*M0_z_Upwind_r[3]+M0_z_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*M0_z_Upwind_r[2]); 

  Ghat_N_x_l[0] = 0.5*(N_x_Upwind_l[3]*alphaDrSurf_l[3]+N_x_Upwind_l[2]*alphaDrSurf_l[2]+N_x_Upwind_l[1]*alphaDrSurf_l[1]+N_x_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_x_l[1] = 0.5*(N_x_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*N_x_Upwind_l[3]+N_x_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_x_Upwind_l[1]); 
  Ghat_N_x_l[2] = 0.5*(N_x_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*N_x_Upwind_l[3]+N_x_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*N_x_Upwind_l[2]); 
  Ghat_N_x_l[3] = 0.5*(N_x_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*N_x_Upwind_l[3]+N_x_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*N_x_Upwind_l[2]); 
  Ghat_N_x_r[0] = 0.5*(N_x_Upwind_r[3]*alphaDrSurf_r[3]+N_x_Upwind_r[2]*alphaDrSurf_r[2]+N_x_Upwind_r[1]*alphaDrSurf_r[1]+N_x_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_x_r[1] = 0.5*(N_x_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*N_x_Upwind_r[3]+N_x_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_x_Upwind_r[1]); 
  Ghat_N_x_r[2] = 0.5*(N_x_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*N_x_Upwind_r[3]+N_x_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*N_x_Upwind_r[2]); 
  Ghat_N_x_r[3] = 0.5*(N_x_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*N_x_Upwind_r[3]+N_x_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*N_x_Upwind_r[2]); 

  Ghat_N_y_l[0] = 0.5*(N_y_Upwind_l[3]*alphaDrSurf_l[3]+N_y_Upwind_l[2]*alphaDrSurf_l[2]+N_y_Upwind_l[1]*alphaDrSurf_l[1]+N_y_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_y_l[1] = 0.5*(N_y_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*N_y_Upwind_l[3]+N_y_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_y_Upwind_l[1]); 
  Ghat_N_y_l[2] = 0.5*(N_y_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*N_y_Upwind_l[3]+N_y_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*N_y_Upwind_l[2]); 
  Ghat_N_y_l[3] = 0.5*(N_y_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*N_y_Upwind_l[3]+N_y_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*N_y_Upwind_l[2]); 
  Ghat_N_y_r[0] = 0.5*(N_y_Upwind_r[3]*alphaDrSurf_r[3]+N_y_Upwind_r[2]*alphaDrSurf_r[2]+N_y_Upwind_r[1]*alphaDrSurf_r[1]+N_y_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_y_r[1] = 0.5*(N_y_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*N_y_Upwind_r[3]+N_y_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_y_Upwind_r[1]); 
  Ghat_N_y_r[2] = 0.5*(N_y_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*N_y_Upwind_r[3]+N_y_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*N_y_Upwind_r[2]); 
  Ghat_N_y_r[3] = 0.5*(N_y_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*N_y_Upwind_r[3]+N_y_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*N_y_Upwind_r[2]); 

  Ghat_N_z_l[0] = 0.5*(N_z_Upwind_l[3]*alphaDrSurf_l[3]+N_z_Upwind_l[2]*alphaDrSurf_l[2]+N_z_Upwind_l[1]*alphaDrSurf_l[1]+N_z_Upwind_l[0]*alphaDrSurf_l[0]); 
  Ghat_N_z_l[1] = 0.5*(N_z_Upwind_l[2]*alphaDrSurf_l[3]+alphaDrSurf_l[2]*N_z_Upwind_l[3]+N_z_Upwind_l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*N_z_Upwind_l[1]); 
  Ghat_N_z_l[2] = 0.5*(N_z_Upwind_l[1]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*N_z_Upwind_l[3]+N_z_Upwind_l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*N_z_Upwind_l[2]); 
  Ghat_N_z_l[3] = 0.5*(N_z_Upwind_l[0]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*N_z_Upwind_l[3]+N_z_Upwind_l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*N_z_Upwind_l[2]); 
  Ghat_N_z_r[0] = 0.5*(N_z_Upwind_r[3]*alphaDrSurf_r[3]+N_z_Upwind_r[2]*alphaDrSurf_r[2]+N_z_Upwind_r[1]*alphaDrSurf_r[1]+N_z_Upwind_r[0]*alphaDrSurf_r[0]); 
  Ghat_N_z_r[1] = 0.5*(N_z_Upwind_r[2]*alphaDrSurf_r[3]+alphaDrSurf_r[2]*N_z_Upwind_r[3]+N_z_Upwind_r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*N_z_Upwind_r[1]); 
  Ghat_N_z_r[2] = 0.5*(N_z_Upwind_r[1]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*N_z_Upwind_r[3]+N_z_Upwind_r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*N_z_Upwind_r[2]); 
  Ghat_N_z_r[3] = 0.5*(N_z_Upwind_r[0]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*N_z_Upwind_r[3]+N_z_Upwind_r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*N_z_Upwind_r[2]); 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_r[2]-0.7071067811865475*Ghat_F_0_l[2])*dv1par; 
  out_F_0[3] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[4] += (0.7071067811865475*Ghat_F_0_r[3]-0.7071067811865475*Ghat_F_0_l[3])*dv1par; 
  out_F_0[5] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[6] += 1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dv1par; 
  out_F_0[7] += (1.5811388300841895*Ghat_F_0_r[0]-1.5811388300841895*Ghat_F_0_l[0])*dv1par; 
  out_F_0[8] += 1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dv1par; 
  out_F_0[9] += (1.5811388300841898*Ghat_F_0_r[1]-1.5811388300841898*Ghat_F_0_l[1])*dv1par; 
  out_F_0[10] += (1.5811388300841898*Ghat_F_0_r[2]-1.5811388300841898*Ghat_F_0_l[2])*dv1par; 
  out_F_0[11] += (1.5811388300841895*Ghat_F_0_r[3]-1.5811388300841895*Ghat_F_0_l[3])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_r[2]-0.7071067811865475*Ghat_G_1_l[2])*dv1par; 
  out_G_1[3] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[4] += (0.7071067811865475*Ghat_G_1_r[3]-0.7071067811865475*Ghat_G_1_l[3])*dv1par; 
  out_G_1[5] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[6] += 1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dv1par; 
  out_G_1[7] += (1.5811388300841895*Ghat_G_1_r[0]-1.5811388300841895*Ghat_G_1_l[0])*dv1par; 
  out_G_1[8] += 1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dv1par; 
  out_G_1[9] += (1.5811388300841898*Ghat_G_1_r[1]-1.5811388300841898*Ghat_G_1_l[1])*dv1par; 
  out_G_1[10] += (1.5811388300841898*Ghat_G_1_r[2]-1.5811388300841898*Ghat_G_1_l[2])*dv1par; 
  out_G_1[11] += (1.5811388300841895*Ghat_G_1_r[3]-1.5811388300841895*Ghat_G_1_l[3])*dv1par; 
  out_M0_x[0] += (0.7071067811865475*Ghat_M0_x_r[0]-0.7071067811865475*Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[1] += (0.7071067811865475*Ghat_M0_x_r[1]-0.7071067811865475*Ghat_M0_x_l[1])*dv1par; 
  out_M0_x[2] += (0.7071067811865475*Ghat_M0_x_r[2]-0.7071067811865475*Ghat_M0_x_l[2])*dv1par; 
  out_M0_x[3] += 1.224744871391589*(Ghat_M0_x_r[0]+Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[4] += (0.7071067811865475*Ghat_M0_x_r[3]-0.7071067811865475*Ghat_M0_x_l[3])*dv1par; 
  out_M0_x[5] += 1.224744871391589*(Ghat_M0_x_r[1]+Ghat_M0_x_l[1])*dv1par; 
  out_M0_x[6] += 1.224744871391589*(Ghat_M0_x_r[2]+Ghat_M0_x_l[2])*dv1par; 
  out_M0_x[7] += (1.5811388300841895*Ghat_M0_x_r[0]-1.5811388300841895*Ghat_M0_x_l[0])*dv1par; 
  out_M0_x[8] += 1.224744871391589*(Ghat_M0_x_r[3]+Ghat_M0_x_l[3])*dv1par; 
  out_M0_x[9] += (1.5811388300841898*Ghat_M0_x_r[1]-1.5811388300841898*Ghat_M0_x_l[1])*dv1par; 
  out_M0_x[10] += (1.5811388300841898*Ghat_M0_x_r[2]-1.5811388300841898*Ghat_M0_x_l[2])*dv1par; 
  out_M0_x[11] += (1.5811388300841895*Ghat_M0_x_r[3]-1.5811388300841895*Ghat_M0_x_l[3])*dv1par; 
  out_M0_y[0] += (0.7071067811865475*Ghat_M0_y_r[0]-0.7071067811865475*Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[1] += (0.7071067811865475*Ghat_M0_y_r[1]-0.7071067811865475*Ghat_M0_y_l[1])*dv1par; 
  out_M0_y[2] += (0.7071067811865475*Ghat_M0_y_r[2]-0.7071067811865475*Ghat_M0_y_l[2])*dv1par; 
  out_M0_y[3] += 1.224744871391589*(Ghat_M0_y_r[0]+Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[4] += (0.7071067811865475*Ghat_M0_y_r[3]-0.7071067811865475*Ghat_M0_y_l[3])*dv1par; 
  out_M0_y[5] += 1.224744871391589*(Ghat_M0_y_r[1]+Ghat_M0_y_l[1])*dv1par; 
  out_M0_y[6] += 1.224744871391589*(Ghat_M0_y_r[2]+Ghat_M0_y_l[2])*dv1par; 
  out_M0_y[7] += (1.5811388300841895*Ghat_M0_y_r[0]-1.5811388300841895*Ghat_M0_y_l[0])*dv1par; 
  out_M0_y[8] += 1.224744871391589*(Ghat_M0_y_r[3]+Ghat_M0_y_l[3])*dv1par; 
  out_M0_y[9] += (1.5811388300841898*Ghat_M0_y_r[1]-1.5811388300841898*Ghat_M0_y_l[1])*dv1par; 
  out_M0_y[10] += (1.5811388300841898*Ghat_M0_y_r[2]-1.5811388300841898*Ghat_M0_y_l[2])*dv1par; 
  out_M0_y[11] += (1.5811388300841895*Ghat_M0_y_r[3]-1.5811388300841895*Ghat_M0_y_l[3])*dv1par; 
  out_M0_z[0] += (0.7071067811865475*Ghat_M0_z_r[0]-0.7071067811865475*Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[1] += (0.7071067811865475*Ghat_M0_z_r[1]-0.7071067811865475*Ghat_M0_z_l[1])*dv1par; 
  out_M0_z[2] += (0.7071067811865475*Ghat_M0_z_r[2]-0.7071067811865475*Ghat_M0_z_l[2])*dv1par; 
  out_M0_z[3] += 1.224744871391589*(Ghat_M0_z_r[0]+Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[4] += (0.7071067811865475*Ghat_M0_z_r[3]-0.7071067811865475*Ghat_M0_z_l[3])*dv1par; 
  out_M0_z[5] += 1.224744871391589*(Ghat_M0_z_r[1]+Ghat_M0_z_l[1])*dv1par; 
  out_M0_z[6] += 1.224744871391589*(Ghat_M0_z_r[2]+Ghat_M0_z_l[2])*dv1par; 
  out_M0_z[7] += (1.5811388300841895*Ghat_M0_z_r[0]-1.5811388300841895*Ghat_M0_z_l[0])*dv1par; 
  out_M0_z[8] += 1.224744871391589*(Ghat_M0_z_r[3]+Ghat_M0_z_l[3])*dv1par; 
  out_M0_z[9] += (1.5811388300841898*Ghat_M0_z_r[1]-1.5811388300841898*Ghat_M0_z_l[1])*dv1par; 
  out_M0_z[10] += (1.5811388300841898*Ghat_M0_z_r[2]-1.5811388300841898*Ghat_M0_z_l[2])*dv1par; 
  out_M0_z[11] += (1.5811388300841895*Ghat_M0_z_r[3]-1.5811388300841895*Ghat_M0_z_l[3])*dv1par; 
  out_N_x[0] += (0.7071067811865475*Ghat_N_x_r[0]-0.7071067811865475*Ghat_N_x_l[0])*dv1par; 
  out_N_x[1] += (0.7071067811865475*Ghat_N_x_r[1]-0.7071067811865475*Ghat_N_x_l[1])*dv1par; 
  out_N_x[2] += (0.7071067811865475*Ghat_N_x_r[2]-0.7071067811865475*Ghat_N_x_l[2])*dv1par; 
  out_N_x[3] += 1.224744871391589*(Ghat_N_x_r[0]+Ghat_N_x_l[0])*dv1par; 
  out_N_x[4] += (0.7071067811865475*Ghat_N_x_r[3]-0.7071067811865475*Ghat_N_x_l[3])*dv1par; 
  out_N_x[5] += 1.224744871391589*(Ghat_N_x_r[1]+Ghat_N_x_l[1])*dv1par; 
  out_N_x[6] += 1.224744871391589*(Ghat_N_x_r[2]+Ghat_N_x_l[2])*dv1par; 
  out_N_x[7] += (1.5811388300841895*Ghat_N_x_r[0]-1.5811388300841895*Ghat_N_x_l[0])*dv1par; 
  out_N_x[8] += 1.224744871391589*(Ghat_N_x_r[3]+Ghat_N_x_l[3])*dv1par; 
  out_N_x[9] += (1.5811388300841898*Ghat_N_x_r[1]-1.5811388300841898*Ghat_N_x_l[1])*dv1par; 
  out_N_x[10] += (1.5811388300841898*Ghat_N_x_r[2]-1.5811388300841898*Ghat_N_x_l[2])*dv1par; 
  out_N_x[11] += (1.5811388300841895*Ghat_N_x_r[3]-1.5811388300841895*Ghat_N_x_l[3])*dv1par; 
  out_N_y[0] += (0.7071067811865475*Ghat_N_y_r[0]-0.7071067811865475*Ghat_N_y_l[0])*dv1par; 
  out_N_y[1] += (0.7071067811865475*Ghat_N_y_r[1]-0.7071067811865475*Ghat_N_y_l[1])*dv1par; 
  out_N_y[2] += (0.7071067811865475*Ghat_N_y_r[2]-0.7071067811865475*Ghat_N_y_l[2])*dv1par; 
  out_N_y[3] += 1.224744871391589*(Ghat_N_y_r[0]+Ghat_N_y_l[0])*dv1par; 
  out_N_y[4] += (0.7071067811865475*Ghat_N_y_r[3]-0.7071067811865475*Ghat_N_y_l[3])*dv1par; 
  out_N_y[5] += 1.224744871391589*(Ghat_N_y_r[1]+Ghat_N_y_l[1])*dv1par; 
  out_N_y[6] += 1.224744871391589*(Ghat_N_y_r[2]+Ghat_N_y_l[2])*dv1par; 
  out_N_y[7] += (1.5811388300841895*Ghat_N_y_r[0]-1.5811388300841895*Ghat_N_y_l[0])*dv1par; 
  out_N_y[8] += 1.224744871391589*(Ghat_N_y_r[3]+Ghat_N_y_l[3])*dv1par; 
  out_N_y[9] += (1.5811388300841898*Ghat_N_y_r[1]-1.5811388300841898*Ghat_N_y_l[1])*dv1par; 
  out_N_y[10] += (1.5811388300841898*Ghat_N_y_r[2]-1.5811388300841898*Ghat_N_y_l[2])*dv1par; 
  out_N_y[11] += (1.5811388300841895*Ghat_N_y_r[3]-1.5811388300841895*Ghat_N_y_l[3])*dv1par; 
  out_N_z[0] += (0.7071067811865475*Ghat_N_z_r[0]-0.7071067811865475*Ghat_N_z_l[0])*dv1par; 
  out_N_z[1] += (0.7071067811865475*Ghat_N_z_r[1]-0.7071067811865475*Ghat_N_z_l[1])*dv1par; 
  out_N_z[2] += (0.7071067811865475*Ghat_N_z_r[2]-0.7071067811865475*Ghat_N_z_l[2])*dv1par; 
  out_N_z[3] += 1.224744871391589*(Ghat_N_z_r[0]+Ghat_N_z_l[0])*dv1par; 
  out_N_z[4] += (0.7071067811865475*Ghat_N_z_r[3]-0.7071067811865475*Ghat_N_z_l[3])*dv1par; 
  out_N_z[5] += 1.224744871391589*(Ghat_N_z_r[1]+Ghat_N_z_l[1])*dv1par; 
  out_N_z[6] += 1.224744871391589*(Ghat_N_z_r[2]+Ghat_N_z_l[2])*dv1par; 
  out_N_z[7] += (1.5811388300841895*Ghat_N_z_r[0]-1.5811388300841895*Ghat_N_z_l[0])*dv1par; 
  out_N_z[8] += 1.224744871391589*(Ghat_N_z_r[3]+Ghat_N_z_l[3])*dv1par; 
  out_N_z[9] += (1.5811388300841898*Ghat_N_z_r[1]-1.5811388300841898*Ghat_N_z_l[1])*dv1par; 
  out_N_z[10] += (1.5811388300841898*Ghat_N_z_r[2]-1.5811388300841898*Ghat_N_z_l[2])*dv1par; 
  out_N_z[11] += (1.5811388300841895*Ghat_N_z_r[3]-1.5811388300841895*Ghat_N_z_l[3])*dv1par; 

  return 0.;

} 
