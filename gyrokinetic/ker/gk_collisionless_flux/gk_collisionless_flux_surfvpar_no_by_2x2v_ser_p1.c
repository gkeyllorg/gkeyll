#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_no_by_2x2v_ser_p1(
    const double *w, const double *dxv, 
    const double *vmap_prime_l, const double *vmap_prime_r,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv,
    const double *bmag, const double *yfieldL, const double *yfieldR, const double *JfL, const double *JfR,
    double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_l,vmap_prime_r: velocity space mapping derivative in left and right cells.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // dgv: volume DG geometry.
  // gkdgv: gyrokinetic volume DG geometry.
  // bmag: magnetic field amplitude.
  // yfieldL: Yushmanov field in the left cell.
  // yfieldR: Yushmanov field in the right cell.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //            Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = 1.4142135623730951*vmapSq[0]*m_; 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_; 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_; 

  double *flux_surf_nodal = &flux_surf[24]; 
  double cfl = 0.0; 
  double bmag_quad = 0.0; 
  double B3_quad = 0.0; 
  double Jc_quad = 0.0; 
  double dualcurlbhat_quad[3] = {0.0}; 
  double alphaL_quad = 0.0;
  double alphaR_quad = 0.0;
  double alpha_avg = 0.0;
  double alpha_max = 0.0;
  const double alpha_eps = 1.0e-2;
  double sgn_alpha = 0.0;
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double alphaJfL = 0.0;
  double alphaJfR = 0.0;

  bmag_quad = gkdgv[0].bmag; 
  B3_quad = gkdgv[0].B3; 
  Jc_quad = dgv[0].Jc; 
  dualcurlbhat_quad[0] = gkdgv[0].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[0].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[0].dualcurlbhat.x[2]; 


  alphaL_quad = -(-(0.5590169943749476*yfieldL[47])+0.5590169943749475*yfieldL[46]+0.5590169943749475*yfieldL[45]+0.5590169943749475*yfieldL[44]-0.5590169943749476*yfieldL[43]-0.5590169943749476*yfieldL[42]-0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]-0.4330127018922193*yfieldL[39]+0.4330127018922193*yfieldL[38]+0.4330127018922193*yfieldL[37]-0.25*yfieldL[36]+0.4330127018922193*yfieldL[35]-0.4330127018922193*yfieldL[34]+0.25*yfieldL[33]+0.25*yfieldL[32]-0.4330127018922193*yfieldL[31]-0.4330127018922193*yfieldL[30]+0.25*yfieldL[29]-0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]-0.25*yfieldL[26]-0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(-(0.5590169943749476*yfieldR[47])+0.5590169943749475*yfieldR[46]+0.5590169943749475*yfieldR[45]+0.5590169943749475*yfieldR[44]-0.5590169943749476*yfieldR[43]-0.5590169943749476*yfieldR[42]-0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]+0.4330127018922193*yfieldR[39]-0.4330127018922193*yfieldR[38]-0.4330127018922193*yfieldR[37]-0.25*yfieldR[36]-0.4330127018922193*yfieldR[35]+0.4330127018922193*yfieldR[34]+0.25*yfieldR[33]+0.25*yfieldR[32]+0.4330127018922193*yfieldR[31]+0.4330127018922193*yfieldR[30]+0.25*yfieldR[29]-0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]-0.25*yfieldR[26]-0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[0] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(0.5590169943749476*yfieldL[47]-0.5590169943749475*yfieldL[46]-0.5590169943749475*yfieldL[45]+0.5590169943749475*yfieldL[44]+0.5590169943749476*yfieldL[43]-0.5590169943749476*yfieldL[42]-0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]+0.4330127018922193*yfieldL[39]-0.4330127018922193*yfieldL[38]-0.4330127018922193*yfieldL[37]+0.25*yfieldL[36]+0.4330127018922193*yfieldL[35]+0.4330127018922193*yfieldL[34]-0.25*yfieldL[33]-0.25*yfieldL[32]-0.4330127018922193*yfieldL[31]-0.4330127018922193*yfieldL[30]+0.25*yfieldL[29]+0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]-0.25*yfieldL[26]-0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(0.5590169943749476*yfieldR[47]-0.5590169943749475*yfieldR[46]-0.5590169943749475*yfieldR[45]+0.5590169943749475*yfieldR[44]+0.5590169943749476*yfieldR[43]-0.5590169943749476*yfieldR[42]-0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]-0.4330127018922193*yfieldR[39]+0.4330127018922193*yfieldR[38]+0.4330127018922193*yfieldR[37]+0.25*yfieldR[36]-0.4330127018922193*yfieldR[35]-0.4330127018922193*yfieldR[34]-0.25*yfieldR[33]-0.25*yfieldR[32]+0.4330127018922193*yfieldR[31]+0.4330127018922193*yfieldR[30]+0.25*yfieldR[29]+0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]-0.25*yfieldR[26]-0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[1] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  bmag_quad = gkdgv[1].bmag; 
  B3_quad = gkdgv[1].B3; 
  Jc_quad = dgv[1].Jc; 
  dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 


  alphaL_quad = -(0.5590169943749476*yfieldL[47]-0.5590169943749475*yfieldL[46]+0.5590169943749475*yfieldL[45]-0.5590169943749475*yfieldL[44]-0.5590169943749476*yfieldL[43]+0.5590169943749476*yfieldL[42]-0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]+0.4330127018922193*yfieldL[39]-0.4330127018922193*yfieldL[38]+0.4330127018922193*yfieldL[37]+0.25*yfieldL[36]-0.4330127018922193*yfieldL[35]-0.4330127018922193*yfieldL[34]-0.25*yfieldL[33]+0.25*yfieldL[32]+0.4330127018922193*yfieldL[31]-0.4330127018922193*yfieldL[30]-0.25*yfieldL[29]-0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]+0.25*yfieldL[26]-0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(0.5590169943749476*yfieldR[47]-0.5590169943749475*yfieldR[46]+0.5590169943749475*yfieldR[45]-0.5590169943749475*yfieldR[44]-0.5590169943749476*yfieldR[43]+0.5590169943749476*yfieldR[42]-0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]-0.4330127018922193*yfieldR[39]+0.4330127018922193*yfieldR[38]-0.4330127018922193*yfieldR[37]+0.25*yfieldR[36]+0.4330127018922193*yfieldR[35]+0.4330127018922193*yfieldR[34]-0.25*yfieldR[33]+0.25*yfieldR[32]-0.4330127018922193*yfieldR[31]+0.4330127018922193*yfieldR[30]-0.25*yfieldR[29]-0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]+0.25*yfieldR[26]-0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[2] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(-(0.5590169943749476*yfieldL[47])+0.5590169943749475*yfieldL[46]-0.5590169943749475*yfieldL[45]-0.5590169943749475*yfieldL[44]+0.5590169943749476*yfieldL[43]+0.5590169943749476*yfieldL[42]-0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]-0.4330127018922193*yfieldL[39]+0.4330127018922193*yfieldL[38]-0.4330127018922193*yfieldL[37]-0.25*yfieldL[36]-0.4330127018922193*yfieldL[35]+0.4330127018922193*yfieldL[34]+0.25*yfieldL[33]-0.25*yfieldL[32]+0.4330127018922193*yfieldL[31]-0.4330127018922193*yfieldL[30]-0.25*yfieldL[29]+0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]+0.25*yfieldL[26]-0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(-(0.5590169943749476*yfieldR[47])+0.5590169943749475*yfieldR[46]-0.5590169943749475*yfieldR[45]-0.5590169943749475*yfieldR[44]+0.5590169943749476*yfieldR[43]+0.5590169943749476*yfieldR[42]-0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]+0.4330127018922193*yfieldR[39]-0.4330127018922193*yfieldR[38]+0.4330127018922193*yfieldR[37]-0.25*yfieldR[36]+0.4330127018922193*yfieldR[35]-0.4330127018922193*yfieldR[34]+0.25*yfieldR[33]-0.25*yfieldR[32]-0.4330127018922193*yfieldR[31]+0.4330127018922193*yfieldR[30]-0.25*yfieldR[29]+0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]+0.25*yfieldR[26]-0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[3] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  bmag_quad = gkdgv[2].bmag; 
  B3_quad = gkdgv[2].B3; 
  Jc_quad = dgv[2].Jc; 
  dualcurlbhat_quad[0] = gkdgv[2].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[2].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[2].dualcurlbhat.x[2]; 


  alphaL_quad = -(0.5590169943749476*yfieldL[47]+0.5590169943749475*yfieldL[46]-0.5590169943749475*yfieldL[45]-0.5590169943749475*yfieldL[44]-0.5590169943749476*yfieldL[43]-0.5590169943749476*yfieldL[42]+0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]+0.4330127018922193*yfieldL[39]+0.4330127018922193*yfieldL[38]-0.4330127018922193*yfieldL[37]+0.25*yfieldL[36]-0.4330127018922193*yfieldL[35]-0.4330127018922193*yfieldL[34]+0.25*yfieldL[33]-0.25*yfieldL[32]-0.4330127018922193*yfieldL[31]+0.4330127018922193*yfieldL[30]-0.25*yfieldL[29]-0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]-0.25*yfieldL[26]+0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(0.5590169943749476*yfieldR[47]+0.5590169943749475*yfieldR[46]-0.5590169943749475*yfieldR[45]-0.5590169943749475*yfieldR[44]-0.5590169943749476*yfieldR[43]-0.5590169943749476*yfieldR[42]+0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]-0.4330127018922193*yfieldR[39]-0.4330127018922193*yfieldR[38]+0.4330127018922193*yfieldR[37]+0.25*yfieldR[36]+0.4330127018922193*yfieldR[35]+0.4330127018922193*yfieldR[34]+0.25*yfieldR[33]-0.25*yfieldR[32]+0.4330127018922193*yfieldR[31]-0.4330127018922193*yfieldR[30]-0.25*yfieldR[29]-0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]-0.25*yfieldR[26]+0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[4] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(-(0.5590169943749476*yfieldL[47])-0.5590169943749475*yfieldL[46]+0.5590169943749475*yfieldL[45]-0.5590169943749475*yfieldL[44]+0.5590169943749476*yfieldL[43]-0.5590169943749476*yfieldL[42]+0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]-0.4330127018922193*yfieldL[39]-0.4330127018922193*yfieldL[38]+0.4330127018922193*yfieldL[37]-0.25*yfieldL[36]-0.4330127018922193*yfieldL[35]+0.4330127018922193*yfieldL[34]-0.25*yfieldL[33]+0.25*yfieldL[32]-0.4330127018922193*yfieldL[31]+0.4330127018922193*yfieldL[30]-0.25*yfieldL[29]+0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]-0.25*yfieldL[26]+0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(-(0.5590169943749476*yfieldR[47])-0.5590169943749475*yfieldR[46]+0.5590169943749475*yfieldR[45]-0.5590169943749475*yfieldR[44]+0.5590169943749476*yfieldR[43]-0.5590169943749476*yfieldR[42]+0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]+0.4330127018922193*yfieldR[39]+0.4330127018922193*yfieldR[38]-0.4330127018922193*yfieldR[37]-0.25*yfieldR[36]+0.4330127018922193*yfieldR[35]-0.4330127018922193*yfieldR[34]-0.25*yfieldR[33]+0.25*yfieldR[32]+0.4330127018922193*yfieldR[31]-0.4330127018922193*yfieldR[30]-0.25*yfieldR[29]+0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]-0.25*yfieldR[26]+0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[5] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  bmag_quad = gkdgv[3].bmag; 
  B3_quad = gkdgv[3].B3; 
  Jc_quad = dgv[3].Jc; 
  dualcurlbhat_quad[0] = gkdgv[3].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[3].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[3].dualcurlbhat.x[2]; 


  alphaL_quad = -(-(0.5590169943749476*yfieldL[47])-0.5590169943749475*yfieldL[46]-0.5590169943749475*yfieldL[45]+0.5590169943749475*yfieldL[44]-0.5590169943749476*yfieldL[43]+0.5590169943749476*yfieldL[42]+0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]-0.4330127018922193*yfieldL[39]-0.4330127018922193*yfieldL[38]-0.4330127018922193*yfieldL[37]-0.25*yfieldL[36]+0.4330127018922193*yfieldL[35]-0.4330127018922193*yfieldL[34]-0.25*yfieldL[33]-0.25*yfieldL[32]+0.4330127018922193*yfieldL[31]+0.4330127018922193*yfieldL[30]+0.25*yfieldL[29]-0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]+0.25*yfieldL[26]+0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(-(0.5590169943749476*yfieldR[47])-0.5590169943749475*yfieldR[46]-0.5590169943749475*yfieldR[45]+0.5590169943749475*yfieldR[44]-0.5590169943749476*yfieldR[43]+0.5590169943749476*yfieldR[42]+0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]+0.4330127018922193*yfieldR[39]+0.4330127018922193*yfieldR[38]+0.4330127018922193*yfieldR[37]-0.25*yfieldR[36]-0.4330127018922193*yfieldR[35]+0.4330127018922193*yfieldR[34]-0.25*yfieldR[33]-0.25*yfieldR[32]-0.4330127018922193*yfieldR[31]-0.4330127018922193*yfieldR[30]+0.25*yfieldR[29]-0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]+0.25*yfieldR[26]+0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[6] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(0.5590169943749476*yfieldL[47]+0.5590169943749475*yfieldL[46]+0.5590169943749475*yfieldL[45]+0.5590169943749475*yfieldL[44]+0.5590169943749476*yfieldL[43]+0.5590169943749476*yfieldL[42]+0.5590169943749476*yfieldL[41]+0.5590169943749475*yfieldL[40]+0.4330127018922193*yfieldL[39]+0.4330127018922193*yfieldL[38]+0.4330127018922193*yfieldL[37]+0.25*yfieldL[36]+0.4330127018922193*yfieldL[35]+0.4330127018922193*yfieldL[34]+0.25*yfieldL[33]+0.25*yfieldL[32]+0.4330127018922193*yfieldL[31]+0.4330127018922193*yfieldL[30]+0.25*yfieldL[29]+0.25*yfieldL[28]+0.4330127018922193*yfieldL[27]+0.25*yfieldL[26]+0.25*yfieldL[25]+0.25*yfieldL[24])/m_/bmag_quad * B3_quad ;
  alphaR_quad = -(0.5590169943749476*yfieldR[47]+0.5590169943749475*yfieldR[46]+0.5590169943749475*yfieldR[45]+0.5590169943749475*yfieldR[44]+0.5590169943749476*yfieldR[43]+0.5590169943749476*yfieldR[42]+0.5590169943749476*yfieldR[41]+0.5590169943749475*yfieldR[40]-0.4330127018922193*yfieldR[39]-0.4330127018922193*yfieldR[38]-0.4330127018922193*yfieldR[37]+0.25*yfieldR[36]-0.4330127018922193*yfieldR[35]-0.4330127018922193*yfieldR[34]+0.25*yfieldR[33]+0.25*yfieldR[32]-0.4330127018922193*yfieldR[31]-0.4330127018922193*yfieldR[30]+0.25*yfieldR[29]+0.25*yfieldR[28]-0.4330127018922193*yfieldR[27]+0.25*yfieldR[26]+0.25*yfieldR[25]+0.25*yfieldR[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[7] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
