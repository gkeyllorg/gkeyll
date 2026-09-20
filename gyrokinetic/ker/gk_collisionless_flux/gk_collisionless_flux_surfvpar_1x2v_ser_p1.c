#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x2v_ser_p1(
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
  double rdvpar2 = 2.0/dxv[1];

  double hamil[12] = {0.}; 
  hamil[0] = vmapSq[0]*m_; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[8] = vmapSq[2]*m_; 

  double *flux_surf_nodal = &flux_surf[6]; 
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


  alphaL_quad = -(0.7905694150420947*yfieldL[11]-0.7905694150420948*yfieldL[10]-0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]+0.6123724356957944*yfieldL[7]-0.6123724356957944*yfieldL[6]+0.3535533905932737*yfieldL[5]-0.6123724356957944*yfieldL[4]-0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]-0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * B3_quad -(0.7905694150420947*yfieldL[11]-0.7905694150420948*yfieldL[10]-0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]+0.6123724356957944*yfieldL[7]-0.6123724356957944*yfieldL[6]+0.3535533905932737*yfieldL[5]-0.6123724356957944*yfieldL[4]-0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]-0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);
  alphaR_quad = -(0.7905694150420947*yfieldR[11]-0.7905694150420948*yfieldR[10]-0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]-0.6123724356957944*yfieldR[7]+0.6123724356957944*yfieldR[6]+0.3535533905932737*yfieldR[5]+0.6123724356957944*yfieldR[4]-0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]-0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * B3_quad -(0.7905694150420947*yfieldR[11]-0.7905694150420948*yfieldR[10]-0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]-0.6123724356957944*yfieldR[7]+0.6123724356957944*yfieldR[6]+0.3535533905932737*yfieldR[5]+0.6123724356957944*yfieldR[4]-0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]-0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.7905694150420947*JfL[11]-0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]+0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.7905694150420947*JfR[11]-0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[0] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(-(0.7905694150420947*yfieldL[11])+0.7905694150420948*yfieldL[10]-0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]-0.6123724356957944*yfieldL[7]+0.6123724356957944*yfieldL[6]-0.3535533905932737*yfieldL[5]-0.6123724356957944*yfieldL[4]+0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]-0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * B3_quad -(-(0.7905694150420947*yfieldL[11])+0.7905694150420948*yfieldL[10]-0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]-0.6123724356957944*yfieldL[7]+0.6123724356957944*yfieldL[6]-0.3535533905932737*yfieldL[5]-0.6123724356957944*yfieldL[4]+0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]-0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);
  alphaR_quad = -(-(0.7905694150420947*yfieldR[11])+0.7905694150420948*yfieldR[10]-0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]+0.6123724356957944*yfieldR[7]-0.6123724356957944*yfieldR[6]-0.3535533905932737*yfieldR[5]+0.6123724356957944*yfieldR[4]+0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]-0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * B3_quad -(-(0.7905694150420947*yfieldR[11])+0.7905694150420948*yfieldR[10]-0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]+0.6123724356957944*yfieldR[7]-0.6123724356957944*yfieldR[6]-0.3535533905932737*yfieldR[5]+0.6123724356957944*yfieldR[4]+0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]-0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.7905694150420947*JfL[11])+0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]+0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.7905694150420947*JfR[11])+0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0])/vmap_prime_r[0];
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


  alphaL_quad = -(-(0.7905694150420947*yfieldL[11])-0.7905694150420948*yfieldL[10]+0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]-0.6123724356957944*yfieldL[7]-0.6123724356957944*yfieldL[6]-0.3535533905932737*yfieldL[5]+0.6123724356957944*yfieldL[4]-0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]+0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * B3_quad -(-(0.7905694150420947*yfieldL[11])-0.7905694150420948*yfieldL[10]+0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]-0.6123724356957944*yfieldL[7]-0.6123724356957944*yfieldL[6]-0.3535533905932737*yfieldL[5]+0.6123724356957944*yfieldL[4]-0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]+0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);
  alphaR_quad = -(-(0.7905694150420947*yfieldR[11])-0.7905694150420948*yfieldR[10]+0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]+0.6123724356957944*yfieldR[7]+0.6123724356957944*yfieldR[6]-0.3535533905932737*yfieldR[5]-0.6123724356957944*yfieldR[4]-0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]+0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * B3_quad -(-(0.7905694150420947*yfieldR[11])-0.7905694150420948*yfieldR[10]+0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]+0.6123724356957944*yfieldR[7]+0.6123724356957944*yfieldR[6]-0.3535533905932737*yfieldR[5]-0.6123724356957944*yfieldR[4]-0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]+0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(0.7905694150420947*JfL[11])-0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]+0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(0.7905694150420947*JfR[11])-0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[2] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  alphaL_quad = -(0.7905694150420947*yfieldL[11]+0.7905694150420948*yfieldL[10]+0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]+0.6123724356957944*yfieldL[7]+0.6123724356957944*yfieldL[6]+0.3535533905932737*yfieldL[5]+0.6123724356957944*yfieldL[4]+0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]+0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * B3_quad -(0.7905694150420947*yfieldL[11]+0.7905694150420948*yfieldL[10]+0.7905694150420948*yfieldL[9]+0.7905694150420947*yfieldL[8]+0.6123724356957944*yfieldL[7]+0.6123724356957944*yfieldL[6]+0.3535533905932737*yfieldL[5]+0.6123724356957944*yfieldL[4]+0.3535533905932737*yfieldL[3]+0.6123724356957944*yfieldL[2]+0.3535533905932737*yfieldL[1]+0.3535533905932737*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);
  alphaR_quad = -(0.7905694150420947*yfieldR[11]+0.7905694150420948*yfieldR[10]+0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]-0.6123724356957944*yfieldR[7]-0.6123724356957944*yfieldR[6]+0.3535533905932737*yfieldR[5]-0.6123724356957944*yfieldR[4]+0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]+0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * B3_quad -(0.7905694150420947*yfieldR[11]+0.7905694150420948*yfieldR[10]+0.7905694150420948*yfieldR[9]+0.7905694150420947*yfieldR[8]-0.6123724356957944*yfieldR[7]-0.6123724356957944*yfieldR[6]+0.3535533905932737*yfieldR[5]-0.6123724356957944*yfieldR[4]+0.3535533905932737*yfieldR[3]-0.6123724356957944*yfieldR[2]+0.3535533905932737*yfieldR[1]+0.3535533905932737*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.6123724356957944*hamil[2]-2.371708245126284*hamil[8]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (0.7905694150420947*JfL[11]+0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]+0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (0.7905694150420947*JfR[11]+0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[3] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
