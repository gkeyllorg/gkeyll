#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfvpar_1x1v_ser_p1(
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

  double hamil[6] = {0.}; 
  hamil[0] = 0.7071067811865475*vmapSq[0]*m_; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double *flux_surf_nodal = &flux_surf[3]; 
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


  alphaL_quad = -(-(1.1180339887498951*yfieldL[5])+1.118033988749895*yfieldL[4]-0.8660254037844386*yfieldL[3]+0.8660254037844386*yfieldL[2]-0.5*yfieldL[1]+0.5*yfieldL[0])/m_/bmag_quad * B3_quad -(-(1.1180339887498951*yfieldL[5])+1.118033988749895*yfieldL[4]-0.8660254037844386*yfieldL[3]+0.8660254037844386*yfieldL[2]-0.5*yfieldL[1]+0.5*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.8660254037844386*hamil[2]-3.3541019662496847*hamil[4]))/vmap[1]);
  alphaR_quad = -(-(1.1180339887498951*yfieldR[5])+1.118033988749895*yfieldR[4]+0.8660254037844386*yfieldR[3]-0.8660254037844386*yfieldR[2]-0.5*yfieldR[1]+0.5*yfieldR[0])/m_/bmag_quad * B3_quad -(-(1.1180339887498951*yfieldR[5])+1.118033988749895*yfieldR[4]+0.8660254037844386*yfieldR[3]-0.8660254037844386*yfieldR[2]-0.5*yfieldR[1]+0.5*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.8660254037844386*hamil[2]-3.3541019662496847*hamil[4]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (-(1.1180339887498951*JfL[5])+1.118033988749895*JfL[4]-0.8660254037844386*JfL[3]+0.8660254037844386*JfL[2]-0.5*JfL[1]+0.5*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (-(1.1180339887498951*JfR[5])+1.118033988749895*JfR[4]+0.8660254037844386*JfR[3]-0.8660254037844386*JfR[2]-0.5*JfR[1]+0.5*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[0] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  bmag_quad = gkdgv[1].bmag; 
  B3_quad = gkdgv[1].B3; 
  Jc_quad = dgv[1].Jc; 
  dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
  dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
  dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 


  alphaL_quad = -(1.1180339887498951*yfieldL[5]+1.118033988749895*yfieldL[4]+0.8660254037844386*yfieldL[3]+0.8660254037844386*yfieldL[2]+0.5*yfieldL[1]+0.5*yfieldL[0])/m_/bmag_quad * B3_quad -(1.1180339887498951*yfieldL[5]+1.118033988749895*yfieldL[4]+0.8660254037844386*yfieldL[3]+0.8660254037844386*yfieldL[2]+0.5*yfieldL[1]+0.5*yfieldL[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.8660254037844386*hamil[2]-3.3541019662496847*hamil[4]))/vmap[1]);
  alphaR_quad = -(1.1180339887498951*yfieldR[5]+1.118033988749895*yfieldR[4]-0.8660254037844386*yfieldR[3]-0.8660254037844386*yfieldR[2]+0.5*yfieldR[1]+0.5*yfieldR[0])/m_/bmag_quad * B3_quad -(1.1180339887498951*yfieldR[5]+1.118033988749895*yfieldR[4]-0.8660254037844386*yfieldR[3]-0.8660254037844386*yfieldR[2]+0.5*yfieldR[1]+0.5*yfieldR[0])/m_/bmag_quad * 1/q_*dualcurlbhat_quad[2]*((0.8164965809277261*(0.8660254037844386*hamil[2]-3.3541019662496847*hamil[4]))/vmap[1]);

  cfl = fmax(fmax(fabs(alphaL_quad), fabs(alphaR_quad)), fabs(cfl));
  JfL_quad =  (1.1180339887498951*JfL[5]+1.118033988749895*JfL[4]+0.8660254037844386*JfL[3]+0.8660254037844386*JfL[2]+0.5*JfL[1]+0.5*JfL[0])/vmap_prime_l[0];
  JfR_quad =  (1.1180339887498951*JfR[5]+1.118033988749895*JfR[4]-0.8660254037844386*JfR[3]-0.8660254037844386*JfR[2]+0.5*JfR[1]+0.5*JfR[0])/vmap_prime_r[0];
  alpha_avg = 0.5*(alphaL_quad + alphaR_quad);
  alpha_max = fmax(fabs(alphaL_quad), fabs(alphaR_quad));
  sgn_alpha = fabs(alpha_avg) > alpha_eps*alpha_max ? gsign(alpha_avg) : 0.0;
  alphaJfL = alphaL_quad*JfL_quad;
  alphaJfR = alphaR_quad*JfR_quad;
  flux_surf_nodal[1] = 0.5*(alphaJfL + alphaJfR) - 0.5*sgn_alpha*(alphaJfR - alphaJfL);

  double vmap_prime_min = fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_r[0]));

  return cfl/vmap_prime_min*2.5*rdvpar2; 

} 
