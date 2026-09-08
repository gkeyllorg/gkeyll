#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_no_by_2x2v_ser_p1(
    const double *w, const double *dxv, 
    const double *vmap_prime_l, const double *vmap_prime_r,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv,
    const double *bmag, const double *yfield, const double *JfL, const double *JfR,
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
  // yfield: Yushmanov field (gradient of Yushmanov potential).
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
  double alpha_quad = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 

  bmag_quad = gkdgv[0].bmag; 
  B3_quad = gkdgv[0].B3; 
  Jc_quad = dgv[0].Jc; 
  dualcurlbhat_quad[0] = gkdgv[0].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[0].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[0].dualcurlbhat.x[2]; 


  alpha_quad = -(-(0.5590169943749476*yfield[47])+0.5590169943749475*yfield[46]+0.5590169943749475*yfield[45]+0.5590169943749475*yfield[44]-0.5590169943749476*yfield[43]-0.5590169943749476*yfield[42]-0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]+0.4330127018922193*yfield[39]-0.4330127018922193*yfield[38]-0.4330127018922193*yfield[37]-0.25*yfield[36]-0.4330127018922193*yfield[35]+0.4330127018922193*yfield[34]+0.25*yfield[33]+0.25*yfield[32]+0.4330127018922193*yfield[31]+0.4330127018922193*yfield[30]+0.25*yfield[29]-0.25*yfield[28]-0.4330127018922193*yfield[27]-0.25*yfield[26]-0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -(0.5590169943749476*yfield[47]-0.5590169943749475*yfield[46]-0.5590169943749475*yfield[45]+0.5590169943749475*yfield[44]+0.5590169943749476*yfield[43]-0.5590169943749476*yfield[42]-0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]-0.4330127018922193*yfield[39]+0.4330127018922193*yfield[38]+0.4330127018922193*yfield[37]+0.25*yfield[36]-0.4330127018922193*yfield[35]-0.4330127018922193*yfield[34]-0.25*yfield[33]-0.25*yfield[32]+0.4330127018922193*yfield[31]+0.4330127018922193*yfield[30]+0.25*yfield[29]+0.25*yfield[28]-0.4330127018922193*yfield[27]-0.25*yfield[26]-0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[1].bmag; 
  B3_quad = gkdgv[1].B3; 
  Jc_quad = dgv[1].Jc; 
  dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 


  alpha_quad = -(0.5590169943749476*yfield[47]-0.5590169943749475*yfield[46]+0.5590169943749475*yfield[45]-0.5590169943749475*yfield[44]-0.5590169943749476*yfield[43]+0.5590169943749476*yfield[42]-0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]-0.4330127018922193*yfield[39]+0.4330127018922193*yfield[38]-0.4330127018922193*yfield[37]+0.25*yfield[36]+0.4330127018922193*yfield[35]+0.4330127018922193*yfield[34]-0.25*yfield[33]+0.25*yfield[32]-0.4330127018922193*yfield[31]+0.4330127018922193*yfield[30]-0.25*yfield[29]-0.25*yfield[28]-0.4330127018922193*yfield[27]+0.25*yfield[26]-0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -(-(0.5590169943749476*yfield[47])+0.5590169943749475*yfield[46]-0.5590169943749475*yfield[45]-0.5590169943749475*yfield[44]+0.5590169943749476*yfield[43]+0.5590169943749476*yfield[42]-0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]+0.4330127018922193*yfield[39]-0.4330127018922193*yfield[38]+0.4330127018922193*yfield[37]-0.25*yfield[36]+0.4330127018922193*yfield[35]-0.4330127018922193*yfield[34]+0.25*yfield[33]-0.25*yfield[32]-0.4330127018922193*yfield[31]+0.4330127018922193*yfield[30]-0.25*yfield[29]+0.25*yfield[28]-0.4330127018922193*yfield[27]+0.25*yfield[26]-0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[2].bmag; 
  B3_quad = gkdgv[2].B3; 
  Jc_quad = dgv[2].Jc; 
  dualcurlbhat_quad[0] = gkdgv[2].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[2].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[2].dualcurlbhat.x[2]; 


  alpha_quad = -(0.5590169943749476*yfield[47]+0.5590169943749475*yfield[46]-0.5590169943749475*yfield[45]-0.5590169943749475*yfield[44]-0.5590169943749476*yfield[43]-0.5590169943749476*yfield[42]+0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]-0.4330127018922193*yfield[39]-0.4330127018922193*yfield[38]+0.4330127018922193*yfield[37]+0.25*yfield[36]+0.4330127018922193*yfield[35]+0.4330127018922193*yfield[34]+0.25*yfield[33]-0.25*yfield[32]+0.4330127018922193*yfield[31]-0.4330127018922193*yfield[30]-0.25*yfield[29]-0.25*yfield[28]-0.4330127018922193*yfield[27]-0.25*yfield[26]+0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -(-(0.5590169943749476*yfield[47])-0.5590169943749475*yfield[46]+0.5590169943749475*yfield[45]-0.5590169943749475*yfield[44]+0.5590169943749476*yfield[43]-0.5590169943749476*yfield[42]+0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]+0.4330127018922193*yfield[39]+0.4330127018922193*yfield[38]-0.4330127018922193*yfield[37]-0.25*yfield[36]+0.4330127018922193*yfield[35]-0.4330127018922193*yfield[34]-0.25*yfield[33]+0.25*yfield[32]+0.4330127018922193*yfield[31]-0.4330127018922193*yfield[30]-0.25*yfield[29]+0.25*yfield[28]-0.4330127018922193*yfield[27]-0.25*yfield[26]+0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[3].bmag; 
  B3_quad = gkdgv[3].B3; 
  Jc_quad = dgv[3].Jc; 
  dualcurlbhat_quad[0] = gkdgv[3].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[3].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[3].dualcurlbhat.x[2]; 


  alpha_quad = -(-(0.5590169943749476*yfield[47])-0.5590169943749475*yfield[46]-0.5590169943749475*yfield[45]+0.5590169943749475*yfield[44]-0.5590169943749476*yfield[43]+0.5590169943749476*yfield[42]+0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]+0.4330127018922193*yfield[39]+0.4330127018922193*yfield[38]+0.4330127018922193*yfield[37]-0.25*yfield[36]-0.4330127018922193*yfield[35]+0.4330127018922193*yfield[34]-0.25*yfield[33]-0.25*yfield[32]-0.4330127018922193*yfield[31]-0.4330127018922193*yfield[30]+0.25*yfield[29]-0.25*yfield[28]-0.4330127018922193*yfield[27]+0.25*yfield[26]+0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -(0.5590169943749476*yfield[47]+0.5590169943749475*yfield[46]+0.5590169943749475*yfield[45]+0.5590169943749475*yfield[44]+0.5590169943749476*yfield[43]+0.5590169943749476*yfield[42]+0.5590169943749476*yfield[41]+0.5590169943749475*yfield[40]-0.4330127018922193*yfield[39]-0.4330127018922193*yfield[38]-0.4330127018922193*yfield[37]+0.25*yfield[36]-0.4330127018922193*yfield[35]-0.4330127018922193*yfield[34]+0.25*yfield[33]+0.25*yfield[32]-0.4330127018922193*yfield[31]-0.4330127018922193*yfield[30]+0.25*yfield[29]+0.25*yfield[28]-0.4330127018922193*yfield[27]+0.25*yfield[26]+0.25*yfield[25]+0.25*yfield[24])/m_/bmag_quad * B3_quad ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
