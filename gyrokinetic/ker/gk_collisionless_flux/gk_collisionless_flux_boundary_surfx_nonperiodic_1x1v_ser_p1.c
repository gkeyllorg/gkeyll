#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_nonperiodic_1x1v_ser_p1(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, const double *bmag,
    const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, const double *phiL, const double *phiR,
    const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR,
    double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // dgs: surface DG geometry.
  // gkdgs: gyrokinetic surface DG geometry.
  // bmag: bmag represented on the surface.
  // jacobgeo_rat_surfL: Ratio of surface conf-space Jacobians in left cell.
  // jacobgeo_rat_surfR: Ratio of surface conf-space Jacobians in right cell.
  // phiL: electrostatic potential in left cell.
  // phiR: electrostatic potential in right cell.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // yfieldL: Yushmanov field (gradient of Yushmanov potential) in left cell.
  // yfieldR: Yushmanov field (gradient of Yushmanov potential) in right cell.
  // flux_surf: surface phase space flux in each direction.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[6] = {0.}; 
  hamil[0] = (phiR[0]-1.7320508075688772*phiR[1])*q_+0.5*vmapSq[0]*m_; 
  hamil[1] = 0.5*vmapSq[1]*m_; 
  hamil[2] = 0.5*vmapSq[2]*m_; 

  double *ypotderL = &yfieldL[0]; 
  double *ypotderR = &yfieldR[0]; 
  ypotderR[1] += -(2.4494897427831783*phiR[0]*q_*rdx2); 
  double *flux_surf_nodal = &flux_surf[0]; 
  double cfl = 0.0; 
  double bmag_quad = 0.0; 
  double Jc_quad = 0.0; 
  double B3_quad = 0.0; 
  double normcurlbhat_quad = 0.0; 
  double area_elem_quad = 0.0; 
  double bhat_quad[3] = {0.0}; 
  double alpha_quad = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 
  double ypotL_quad = 0.0; 
  double ypotR_quad = 0.0; 
  double ypot_der_avg_quad = 0.0; 
  double ypot_der_jump_quad = 0.0; 
  double ypot_der_surf_n[3] = {0.0}; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(1.224744871391589*hamil[1]-3.6742346141747664*hamil[2]))/vmap[1]; 
  mvpar_quad[1] = (1.0*hamil[1])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(3.6742346141747664*hamil[2]+1.224744871391589*hamil[1]))/vmap[1]; 
  double mvparsq_quad[3] = {0.0}; 
  mvparsq_quad[0] = mvpar_quad[0]*mvpar_quad[0]/m_; 
  mvparsq_quad[1] = mvpar_quad[1]*mvpar_quad[1]/m_; 
  mvparsq_quad[2] = mvpar_quad[2]*mvpar_quad[2]/m_; 

  bmag_quad = gkdgs[0].bmag; 
  Jc_quad = gkdgs[0].Jc; 
  B3_quad = gkdgs[0].B3; 
  normcurlbhat_quad = gkdgs[0].normcurlbhat; 
  bhat_quad[0] = gkdgs[0].bhat.x[0]; 
  bhat_quad[1] = gkdgs[0].bhat.x[1]; 
  bhat_quad[2] = gkdgs[0].bhat.x[2]; 
  area_elem_quad = dgs[0].area_elem; 


  alpha_quad = (mvpar_quad[0]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]-1.1618950038622249*JfL[3]-0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = -(0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]+1.1618950038622249*JfR[3]-0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypotR_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[0] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvpar_quad[1]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.9682458365518543*JfL[5])-0.5590169943749475*JfL[4]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = 0.9682458365518543*JfR[5]-0.5590169943749475*JfR[4]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypotR_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[1] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvpar_quad[2]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]+1.1618950038622249*JfL[3]+0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = -(0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]-1.1618950038622249*JfR[3]+0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypotR_quad = 0.7071067811865475*(phiR[0]-1.7320508075688772*phiR[1])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[2] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  ypotderR[0] += -(0.2777777777777778*ypot_der_surf_n[2]*rdx2)-0.4444444444444444*ypot_der_surf_n[1]*rdx2-0.2777777777777778*ypot_der_surf_n[0]*rdx2; 
  ypotderR[1] += 0.48112522432468824*ypot_der_surf_n[2]*rdx2+0.7698003589195012*ypot_der_surf_n[1]*rdx2+0.48112522432468824*ypot_der_surf_n[0]*rdx2; 
  ypotderR[2] += 0.37267799624996495*ypot_der_surf_n[0]*rdx2-0.37267799624996495*ypot_der_surf_n[2]*rdx2; 
  ypotderR[3] += 0.6454972243679029*ypot_der_surf_n[2]*rdx2-0.6454972243679029*ypot_der_surf_n[0]*rdx2; 
  ypotderR[4] += -(0.24845199749997662*ypot_der_surf_n[2]*rdx2)+0.49690399499995325*ypot_der_surf_n[1]*rdx2-0.24845199749997662*ypot_der_surf_n[0]*rdx2; 
  ypotderR[5] += 0.4303314829119352*ypot_der_surf_n[2]*rdx2-0.8606629658238704*ypot_der_surf_n[1]*rdx2+0.4303314829119352*ypot_der_surf_n[0]*rdx2; 


  return cfl*1.5*rdx2; 

} 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_nonperiodic_1x1v_ser_p1(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, const double *bmag,
    const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, const double *phiL, const double *phiR,
    const double *JfL, const double *JfR, double* GKYL_RESTRICT yfieldL, double* GKYL_RESTRICT yfieldR,
    double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // dgs: surface DG geometry.
  // gkdgs: gyrokinetic surface DG geometry.
  // bmag: bmag represented on the surface.
  // jacobgeo_rat_surfL: Ratio of surface conf-space Jacobians in left cell.
  // jacobgeo_rat_surfR: Ratio of surface conf-space Jacobians in right cell.
  // phiL: electrostatic potential in left cell.
  // phiR: electrostatic potential in right cell.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // yfieldL: Yushmanov field (gradient of Yushmanov potential) in left cell.
  // yfieldR: Yushmanov field (gradient of Yushmanov potential) in right cell.
  // flux_surf: surface phase space flux in each direction.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[6] = {0.}; 
  hamil[0] = (1.7320508075688772*phiL[1]+phiL[0])*q_+0.5*vmapSq[0]*m_; 
  hamil[1] = 0.5*vmapSq[1]*m_; 
  hamil[2] = 0.5*vmapSq[2]*m_; 

  double *ypotderL = &yfieldL[0]; 
  double *ypotderR = &yfieldR[0]; 
  double *flux_surf_nodal = &flux_surf[0]; 
  double cfl = 0.0; 
  double bmag_quad = 0.0; 
  double Jc_quad = 0.0; 
  double B3_quad = 0.0; 
  double normcurlbhat_quad = 0.0; 
  double area_elem_quad = 0.0; 
  double bhat_quad[3] = {0.0}; 
  double alpha_quad = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 
  double ypotL_quad = 0.0; 
  double ypotR_quad = 0.0; 
  double ypot_der_avg_quad = 0.0; 
  double ypot_der_jump_quad = 0.0; 
  double ypot_der_surf_n[3] = {0.0}; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(1.224744871391589*hamil[1]-3.6742346141747664*hamil[2]))/vmap[1]; 
  mvpar_quad[1] = (1.0*hamil[1])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(3.6742346141747664*hamil[2]+1.224744871391589*hamil[1]))/vmap[1]; 
  double mvparsq_quad[3] = {0.0}; 
  mvparsq_quad[0] = mvpar_quad[0]*mvpar_quad[0]/m_; 
  mvparsq_quad[1] = mvpar_quad[1]*mvpar_quad[1]/m_; 
  mvparsq_quad[2] = mvpar_quad[2]*mvpar_quad[2]/m_; 

  bmag_quad = gkdgs[0].bmag; 
  Jc_quad = gkdgs[0].Jc; 
  B3_quad = gkdgs[0].B3; 
  normcurlbhat_quad = gkdgs[0].normcurlbhat; 
  bhat_quad[0] = gkdgs[0].bhat.x[0]; 
  bhat_quad[1] = gkdgs[0].bhat.x[1]; 
  bhat_quad[2] = gkdgs[0].bhat.x[2]; 
  area_elem_quad = dgs[0].area_elem; 


  alpha_quad = (mvpar_quad[0]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]-1.1618950038622249*JfL[3]-0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = -(0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]+1.1618950038622249*JfR[3]-0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypotR_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[0] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvpar_quad[1]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.9682458365518543*JfL[5])-0.5590169943749475*JfL[4]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = 0.9682458365518543*JfR[5]-0.5590169943749475*JfR[4]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypotR_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[1] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvpar_quad[2]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]+1.1618950038622249*JfL[3]+0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = -(0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]-1.1618950038622249*JfR[3]+0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypotR_quad = 0.7071067811865475*(1.7320508075688772*phiL[1]+phiL[0])*q_; 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[2] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  ypotderL[0] += 0.2777777777777778*ypot_der_surf_n[2]*rdx2+0.4444444444444444*ypot_der_surf_n[1]*rdx2+0.2777777777777778*ypot_der_surf_n[0]*rdx2; 
  ypotderL[1] += 0.48112522432468824*ypot_der_surf_n[2]*rdx2+0.7698003589195012*ypot_der_surf_n[1]*rdx2+0.48112522432468824*ypot_der_surf_n[0]*rdx2; 
  ypotderL[2] += 0.37267799624996495*ypot_der_surf_n[2]*rdx2-0.37267799624996495*ypot_der_surf_n[0]*rdx2; 
  ypotderL[3] += 0.6454972243679029*ypot_der_surf_n[2]*rdx2-0.6454972243679029*ypot_der_surf_n[0]*rdx2; 
  ypotderL[4] += 0.24845199749997662*ypot_der_surf_n[2]*rdx2-0.49690399499995325*ypot_der_surf_n[1]*rdx2+0.24845199749997662*ypot_der_surf_n[0]*rdx2; 
  ypotderL[5] += 0.4303314829119352*ypot_der_surf_n[2]*rdx2-0.8606629658238704*ypot_der_surf_n[1]*rdx2+0.4303314829119352*ypot_der_surf_n[0]*rdx2; 


  return cfl*1.5*rdx2; 

} 
