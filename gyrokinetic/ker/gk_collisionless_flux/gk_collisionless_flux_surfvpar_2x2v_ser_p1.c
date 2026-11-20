#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_2x2v_ser_p1(
  const double *w, const double *dxv, 
  const double *vmap_prime_l, const double *vmap_prime_r,
  const double *vmap, const double *vmapSq, const double q_, const double m_,
  const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv,
  const double *bmag, const double *phi, const double *apar, const double *apardot, const double *JfL, const double *JfR,
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
  // phi: electrostatic potential.
  // apar: parallel component of vector potential.
  // apardot: time derivative of parallel component of vector potential.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //            Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+1.4142135623730951*(vmapSq[0]*m_+bmag[0]*vmap[2]); 
  hamil[1] = 2.0*phi[1]*q_+1.4142135623730951*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_+1.4142135623730951*bmag[2]*vmap[2]; 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_; 
  hamil[4] = 1.4142135623730951*bmag[0]*vmap[3]; 
  hamil[5] = 2.0*phi[3]*q_+1.4142135623730951*vmap[2]*bmag[3]; 
  hamil[8] = 1.4142135623730951*bmag[1]*vmap[3]; 
  hamil[9] = 1.4142135623730951*bmag[2]*vmap[3]; 
  hamil[12] = 1.4142135623730951*bmag[3]*vmap[3]; 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_; 

  double gradHxgradA_1[24] = {0.}; 

  double gradHxgradA_2[24] = {0.}; 
  gradHxgradA_2[0] = (1.5*apar[1]*hamil[2]-1.5*hamil[1]*apar[2])*rdx2*rdz2; 
  gradHxgradA_2[1] = (1.5*apar[1]*hamil[5]-1.5*hamil[1]*apar[3])*rdx2*rdz2; 
  gradHxgradA_2[2] = (1.5*hamil[2]*apar[3]-1.5*apar[2]*hamil[5])*rdx2*rdz2; 
  gradHxgradA_2[4] = (1.5*apar[1]*hamil[9]-1.5*apar[2]*hamil[8])*rdx2*rdz2; 
  gradHxgradA_2[8] = (1.5*apar[1]*hamil[12]-1.5*apar[3]*hamil[8])*rdx2*rdz2; 
  gradHxgradA_2[9] = (1.5*apar[3]*hamil[9]-1.5*apar[2]*hamil[12])*rdx2*rdz2; 

  double gradHxgradA_3[24] = {0.}; 

  double apargradH_1[24] = {0.}; 
  apargradH_1[0] = 0.8660254037844386*(apar[2]*hamil[5]+apar[0]*hamil[1])*rdx2; 
  apargradH_1[1] = 0.8660254037844386*(apar[3]*hamil[5]+apar[1]*hamil[1])*rdx2; 
  apargradH_1[2] = 0.8660254037844386*(apar[0]*hamil[5]+hamil[1]*apar[2])*rdx2; 
  apargradH_1[4] = 0.8660254037844386*(apar[2]*hamil[12]+apar[0]*hamil[8])*rdx2; 
  apargradH_1[5] = 0.8660254037844386*(apar[1]*hamil[5]+hamil[1]*apar[3])*rdx2; 
  apargradH_1[8] = 0.8660254037844386*(apar[3]*hamil[12]+apar[1]*hamil[8])*rdx2; 
  apargradH_1[9] = 0.8660254037844386*(apar[0]*hamil[12]+apar[2]*hamil[8])*rdx2; 
  apargradH_1[12] = 0.8660254037844386*(apar[1]*hamil[12]+apar[3]*hamil[8])*rdx2; 

  double apargradH_2[24] = {0.}; 

  double apargradH_3[24] = {0.}; 
  apargradH_3[0] = 0.8660254037844386*(apar[1]*hamil[5]+apar[0]*hamil[2])*rdz2; 
  apargradH_3[1] = 0.8660254037844386*(apar[0]*hamil[5]+apar[1]*hamil[2])*rdz2; 
  apargradH_3[2] = 0.8660254037844386*(apar[3]*hamil[5]+apar[2]*hamil[2])*rdz2; 
  apargradH_3[4] = 0.8660254037844386*(apar[1]*hamil[12]+apar[0]*hamil[9])*rdz2; 
  apargradH_3[5] = 0.8660254037844386*(apar[2]*hamil[5]+hamil[2]*apar[3])*rdz2; 
  apargradH_3[8] = 0.8660254037844386*(apar[0]*hamil[12]+apar[1]*hamil[9])*rdz2; 
  apargradH_3[9] = 0.8660254037844386*(apar[3]*hamil[12]+apar[2]*hamil[9])*rdz2; 
  apargradH_3[12] = 0.8660254037844386*(apar[2]*hamil[12]+apar[3]*hamil[9])*rdz2; 

  double *flux_surf_nodal = &flux_surf[24]; 
  double cfl = 0.0; 
  double bmag_quad = 0.0; 
  double B3_quad = 0.0; 
  double dualcurlbhat_quad[3] = {0.0}; 
  double alpha_quad = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 
  double bioverJB_quad[3] = {0.0}; 

  bmag_quad = gkdgv[0].bmag; 
  B3_quad = gkdgv[0].B3; 
  dualcurlbhat_quad[0] = gkdgv[0].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[0].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[0].dualcurlbhat.x[2]; 
  bioverJB_quad[0] = gkdgv[0].bioverJB.x[0]; 
  bioverJB_quad[1] = gkdgv[0].bioverJB.x[1]; 
  bioverJB_quad[2] = gkdgv[0].bioverJB.x[2]; 

  alpha_quad = -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[8]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (-(0.25*apargradH_1[12])+0.25*apargradH_1[9]+0.25*apargradH_1[8]+0.25*apargradH_1[5]-0.25*apargradH_1[4]-0.25*apargradH_1[2]-0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(0.25*gradHxgradA_2[9]+0.25*gradHxgradA_2[8]-0.25*gradHxgradA_2[4]-0.25*gradHxgradA_2[2]-0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(0.5*apardot[3]-0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[8]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (0.25*apargradH_1[12]-0.25*apargradH_1[9]-0.25*apargradH_1[8]+0.25*apargradH_1[5]+0.25*apargradH_1[4]-0.25*apargradH_1[2]-0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(-(0.25*gradHxgradA_2[9])-0.25*gradHxgradA_2[8]+0.25*gradHxgradA_2[4]-0.25*gradHxgradA_2[2]-0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(0.5*apardot[3]-0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[1].bmag; 
  B3_quad = gkdgv[1].B3; 
  dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 
  bioverJB_quad[0] = gkdgv[1].bioverJB.x[0]; 
  bioverJB_quad[1] = gkdgv[1].bioverJB.x[1]; 
  bioverJB_quad[2] = gkdgv[1].bioverJB.x[2]; 

  alpha_quad = -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (0.25*apargradH_1[12]-0.25*apargradH_1[9]+0.25*apargradH_1[8]-0.25*apargradH_1[5]-0.25*apargradH_1[4]+0.25*apargradH_1[2]-0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(-(0.25*gradHxgradA_2[9])+0.25*gradHxgradA_2[8]-0.25*gradHxgradA_2[4]+0.25*gradHxgradA_2[2]-0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(-(0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[9]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (-(0.25*apargradH_1[12])+0.25*apargradH_1[9]-0.25*apargradH_1[8]-0.25*apargradH_1[5]+0.25*apargradH_1[4]+0.25*apargradH_1[2]-0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(0.25*gradHxgradA_2[9]-0.25*gradHxgradA_2[8]+0.25*gradHxgradA_2[4]+0.25*gradHxgradA_2[2]-0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(-(0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]-0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]-0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]-0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]-0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]+0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]-0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[2].bmag; 
  B3_quad = gkdgv[2].B3; 
  dualcurlbhat_quad[0] = gkdgv[2].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[2].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[2].dualcurlbhat.x[2]; 
  bioverJB_quad[0] = gkdgv[2].bioverJB.x[0]; 
  bioverJB_quad[1] = gkdgv[2].bioverJB.x[1]; 
  bioverJB_quad[2] = gkdgv[2].bioverJB.x[2]; 

  alpha_quad = -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((0.4330127018922193*hamil[12]-0.4330127018922193*hamil[8]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (0.25*apargradH_1[12]+0.25*apargradH_1[9]-0.25*apargradH_1[8]-0.25*apargradH_1[5]-0.25*apargradH_1[4]-0.25*apargradH_1[2]+0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(0.25*gradHxgradA_2[9]-0.25*gradHxgradA_2[8]-0.25*gradHxgradA_2[4]-0.25*gradHxgradA_2[2]+0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(-(0.5*apardot[3])-0.5*apardot[2]+0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]+0.25*JfL[12]-0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]+0.25*JfL[9]-0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]+0.25*JfR[12]+0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]+0.25*JfR[9]-0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((-(0.4330127018922193*hamil[12])+0.4330127018922193*hamil[8]-0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (-(0.25*apargradH_1[12])-0.25*apargradH_1[9]+0.25*apargradH_1[8]-0.25*apargradH_1[5]+0.25*apargradH_1[4]-0.25*apargradH_1[2]+0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(-(0.25*gradHxgradA_2[9])+0.25*gradHxgradA_2[8]+0.25*gradHxgradA_2[4]-0.25*gradHxgradA_2[2]+0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(-(0.5*apardot[3])-0.5*apardot[2]+0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]-0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]-0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]-0.25*JfL[12]-0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]-0.25*JfL[9]+0.25*JfL[8]-0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]-0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]-0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]-0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]-0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]-0.25*JfR[12]+0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]-0.25*JfR[9]+0.25*JfR[8]+0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]-0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]-0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  bmag_quad = gkdgv[3].bmag; 
  B3_quad = gkdgv[3].B3; 
  dualcurlbhat_quad[0] = gkdgv[3].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[3].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[3].dualcurlbhat.x[2]; 
  bioverJB_quad[0] = gkdgv[3].bioverJB.x[0]; 
  bioverJB_quad[1] = gkdgv[3].bioverJB.x[1]; 
  bioverJB_quad[2] = gkdgv[3].bioverJB.x[2]; 

  alpha_quad = -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[8]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((-(0.4330127018922193*hamil[12])-0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (-(0.25*apargradH_1[12])-0.25*apargradH_1[9]-0.25*apargradH_1[8]+0.25*apargradH_1[5]-0.25*apargradH_1[4]+0.25*apargradH_1[2]+0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(-(0.25*gradHxgradA_2[9])-0.25*gradHxgradA_2[8]-0.25*gradHxgradA_2[4]+0.25*gradHxgradA_2[2]+0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(0.5*apardot[3]+0.5*apardot[2]+0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-(0.5590169943749476*JfL[23])-0.5590169943749475*JfL[22]-0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]-0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]-0.4330127018922193*JfL[15]-0.4330127018922193*JfL[14]-0.4330127018922193*JfL[13]-0.25*JfL[12]+0.4330127018922193*JfL[11]-0.4330127018922193*JfL[10]-0.25*JfL[9]-0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]-0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.5590169943749476*JfR[23])-0.5590169943749475*JfR[22]-0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]-0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]+0.4330127018922193*JfR[15]+0.4330127018922193*JfR[14]+0.4330127018922193*JfR[13]-0.25*JfR[12]-0.4330127018922193*JfR[11]+0.4330127018922193*JfR[10]-0.25*JfR[9]-0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]-0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * B3_quad; 
  alpha_quad += -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[8]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[1])*rdx2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[0]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -((0.4330127018922193*hamil[12]+0.4330127018922193*hamil[9]+0.4330127018922193*hamil[5]+0.4330127018922193*hamil[2])*rdz2)/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.4330127018922193*hamil[3]-1.6770509831248424*hamil[16]))/vmap[1]); 
  alpha_quad += -1.0/m_/bmag_quad * (0.25*apargradH_1[12]+0.25*apargradH_1[9]+0.25*apargradH_1[8]+0.25*apargradH_1[5]+0.25*apargradH_1[4]+0.25*apargradH_1[2]+0.25*apargradH_1[1]+0.25*apargradH_1[0])*dualcurlbhat_quad[0]; 
  alpha_quad += -1.0/m_/bmag_quad * (0.0)*dualcurlbhat_quad[2]; 
  alpha_quad += -bioverJB_quad[2]/m_ *(0.25*gradHxgradA_2[9]+0.25*gradHxgradA_2[8]+0.25*gradHxgradA_2[4]+0.25*gradHxgradA_2[2]+0.25*gradHxgradA_2[1]+0.25*gradHxgradA_2[0]); 
  alpha_quad += -q_/m_*(0.5*apardot[3]+0.5*apardot[2]+0.5*apardot[1]+0.5*apardot[0]); 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (0.5590169943749476*JfL[23]+0.5590169943749475*JfL[22]+0.5590169943749475*JfL[21]+0.5590169943749475*JfL[20]+0.5590169943749476*JfL[19]+0.5590169943749476*JfL[18]+0.5590169943749476*JfL[17]+0.5590169943749475*JfL[16]+0.4330127018922193*JfL[15]+0.4330127018922193*JfL[14]+0.4330127018922193*JfL[13]+0.25*JfL[12]+0.4330127018922193*JfL[11]+0.4330127018922193*JfL[10]+0.25*JfL[9]+0.25*JfL[8]+0.4330127018922193*JfL[7]+0.4330127018922193*JfL[6]+0.25*JfL[5]+0.25*JfL[4]+0.4330127018922193*JfL[3]+0.25*JfL[2]+0.25*JfL[1]+0.25*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.5590169943749476*JfR[23]+0.5590169943749475*JfR[22]+0.5590169943749475*JfR[21]+0.5590169943749475*JfR[20]+0.5590169943749476*JfR[19]+0.5590169943749476*JfR[18]+0.5590169943749476*JfR[17]+0.5590169943749475*JfR[16]-0.4330127018922193*JfR[15]-0.4330127018922193*JfR[14]-0.4330127018922193*JfR[13]+0.25*JfR[12]-0.4330127018922193*JfR[11]-0.4330127018922193*JfR[10]+0.25*JfR[9]+0.25*JfR[8]-0.4330127018922193*JfR[7]-0.4330127018922193*JfR[6]+0.25*JfR[5]+0.25*JfR[4]-0.4330127018922193*JfR[3]+0.25*JfR[2]+0.25*JfR[1]+0.25*JfR[0])/vmap_prime_r[0];
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;


  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
