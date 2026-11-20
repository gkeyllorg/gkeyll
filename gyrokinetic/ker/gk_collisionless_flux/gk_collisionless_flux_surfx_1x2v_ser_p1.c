#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_surfx_1x2v_ser_p1(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, 
    const double *phi, const double *apar, const double *apardot,
    const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf) 
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
  // phi: electrostatic potential.
  // apar: parallel component of vector potential.
  // apardot: time derivative of parallel component of vector potential.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //            Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[12] = {0.}; 
  hamil[0] = (1.4142135623730951*phi[0]-2.4494897427831783*phi[1])*q_+0.7071067811865475*vmapSq[0]*m_+(bmag[0]-1.7320508075688772*bmag[1])*vmap[2]; 
  hamil[1] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[2] = (bmag[0]-1.7320508075688772*bmag[1])*vmap[3]; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double apar_surf[6] = {0.}; 
  apar_surf[0] = 1.4142135623730951*apar[0]-2.4494897427831783*apar[1]; 

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
  double dHdvpar_quad[3] = {0.0}; 
  dHdvpar_quad[0] = (0.8164965809277261*(0.8660254037844386*hamil[1]-2.5980762113533156*hamil[4]))/vmap[1]; 
  dHdvpar_quad[1] = (0.7071067811865476*hamil[1])/vmap[1]; 
  dHdvpar_quad[2] = (0.8164965809277261*(2.5980762113533156*hamil[4]+0.8660254037844386*hamil[1]))/vmap[1]; 
  double mvparsq_quad[3] = {0.0}; 
  mvparsq_quad[0] = dHdvpar_quad[0]*dHdvpar_quad[0]/m_; 
  mvparsq_quad[1] = dHdvpar_quad[1]*dHdvpar_quad[1]/m_; 
  mvparsq_quad[2] = dHdvpar_quad[2]*dHdvpar_quad[2]/m_; 

  bmag_quad = gkdgs[0].bmag; 
  Jc_quad = gkdgs[0].Jc; 
  B3_quad = gkdgs[0].B3; 
  normcurlbhat_quad = gkdgs[0].normcurlbhat; 
  bhat_quad[0] = gkdgs[0].bhat.x[0]; 
  bhat_quad[1] = gkdgs[0].bhat.x[1]; 
  bhat_quad[2] = gkdgs[0].bhat.x[2]; 
  area_elem_quad = dgs[0].area_elem; 


  alpha_quad = (dHdvpar_quad[0]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.5477225575051661*JfL[11])-0.3162277660168379*JfL[10]+0.5477225575051661*JfL[9]+0.3162277660168379*JfL[8]+0.8215838362577489*JfL[7]+0.4743416490252568*JfL[6]-0.6123724356957944*JfL[5]-0.8215838362577489*JfL[4]-0.3535533905932737*JfL[3]-0.4743416490252568*JfL[2]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = 0.5477225575051661*JfR[11]-0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]-0.8215838362577489*JfR[7]+0.4743416490252568*JfR[6]+0.6123724356957944*JfR[5]+0.8215838362577489*JfR[4]-0.3535533905932737*JfR[3]-0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = (dHdvpar_quad[1]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.6846531968814573*JfL[11]+0.39528470752104733*JfL[10]-0.6846531968814574*JfL[9]-0.3952847075210473*JfL[8]-0.6123724356957944*JfL[5]-0.3535533905932737*JfL[3]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = -(0.6846531968814573*JfR[11])+0.39528470752104733*JfR[10]+0.6846531968814574*JfR[9]-0.3952847075210473*JfR[8]+0.6123724356957944*JfR[5]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = (dHdvpar_quad[2]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.5477225575051661*JfL[11])-0.3162277660168379*JfL[10]+0.5477225575051661*JfL[9]+0.3162277660168379*JfL[8]-0.8215838362577489*JfL[7]-0.4743416490252568*JfL[6]-0.6123724356957944*JfL[5]+0.8215838362577489*JfL[4]-0.3535533905932737*JfL[3]+0.4743416490252568*JfL[2]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = 0.5477225575051661*JfR[11]-0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]+0.8215838362577489*JfR[7]-0.4743416490252568*JfR[6]+0.6123724356957944*JfR[5]-0.8215838362577489*JfR[4]-0.3535533905932737*JfR[3]+0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = (dHdvpar_quad[0]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.5477225575051661*JfL[11]+0.3162277660168379*JfL[10]+0.5477225575051661*JfL[9]+0.3162277660168379*JfL[8]-0.8215838362577489*JfL[7]-0.4743416490252568*JfL[6]+0.6123724356957944*JfL[5]-0.8215838362577489*JfL[4]+0.3535533905932737*JfL[3]-0.4743416490252568*JfL[2]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = -(0.5477225575051661*JfR[11])+0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]+0.8215838362577489*JfR[7]-0.4743416490252568*JfR[6]-0.6123724356957944*JfR[5]+0.8215838362577489*JfR[4]+0.3535533905932737*JfR[3]-0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = (dHdvpar_quad[1]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.6846531968814573*JfL[11])-0.39528470752104733*JfL[10]-0.6846531968814574*JfL[9]-0.3952847075210473*JfL[8]+0.6123724356957944*JfL[5]+0.3535533905932737*JfL[3]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = 0.6846531968814573*JfR[11]-0.39528470752104733*JfR[10]+0.6846531968814574*JfR[9]-0.3952847075210473*JfR[8]-0.6123724356957944*JfR[5]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = (dHdvpar_quad[2]*B3_quad/(m_*bmag_quad)); 
  alpha_quad += mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_); 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.5477225575051661*JfL[11]+0.3162277660168379*JfL[10]+0.5477225575051661*JfL[9]+0.3162277660168379*JfL[8]+0.8215838362577489*JfL[7]+0.4743416490252568*JfL[6]+0.6123724356957944*JfL[5]+0.8215838362577489*JfL[4]+0.3535533905932737*JfL[3]+0.4743416490252568*JfL[2]+0.6123724356957944*JfL[1]+0.3535533905932737*JfL[0]; 
  JfR_quad = -(0.5477225575051661*JfR[11])+0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]-0.8215838362577489*JfR[7]+0.4743416490252568*JfR[6]-0.6123724356957944*JfR[5]-0.8215838362577489*JfR[4]+0.3535533905932737*JfR[3]+0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 


  return cfl*1.5*rdx2; 

} 
