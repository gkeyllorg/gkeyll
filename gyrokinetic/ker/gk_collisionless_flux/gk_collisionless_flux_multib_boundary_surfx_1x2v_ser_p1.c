#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH void gk_collisionless_flux_multib_boundary_surfx_1x2v_ser_p1_alpha_surf(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs,
    const double *bmag, const double *phi, double *alpha_surf)
{
  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[12] = {0.}; 
  hamil[0] = (1.4142135623730951*phi[0]-2.4494897427831783*phi[1])*q_+0.7071067811865475*vmapSq[0]*m_+(bmag[0]-1.7320508075688772*bmag[1])*vmap[2]; 
  hamil[1] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[2] = (bmag[0]-1.7320508075688772*bmag[1])*vmap[3]; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double bmag_quad = 0.0; 
  double Jc_quad = 0.0; 
  double B3_quad = 0.0; 
  double normcurlbhat_quad = 0.0; 
  double area_elem_quad = 0.0; 
  double bhat_quad[3] = {0.0}; 
  double alpha_quad = 0.0; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(0.8660254037844386*hamil[1]-2.5980762113533156*hamil[4]))/vmap[1]; 
  mvpar_quad[1] = (0.7071067811865476*hamil[1])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(2.5980762113533156*hamil[4]+0.8660254037844386*hamil[1]))/vmap[1]; 
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

  alpha_surf[0] = alpha_quad;

  alpha_quad = (mvpar_quad[1]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  alpha_surf[1] = alpha_quad;

  alpha_quad = (mvpar_quad[2]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  alpha_surf[2] = alpha_quad;

  alpha_quad = (mvpar_quad[0]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  alpha_surf[3] = alpha_quad;

  alpha_quad = (mvpar_quad[1]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  alpha_surf[4] = alpha_quad;

  alpha_quad = (mvpar_quad[2]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad; 

  alpha_surf[5] = alpha_quad;

}

GKYL_CU_DH double gk_collisionless_flux_multib_boundary_surfx_1x2v_ser_p1(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs,
    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR, const double *phi,
    const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf)
{
  double rdx2 = 2.0/dxv[0];
  double alpha_surf[6] = {0.0};
  gk_collisionless_flux_multib_boundary_surfx_1x2v_ser_p1_alpha_surf(w, dxv, vmap, vmapSq, q_, m_, dgs, gkdgs, bmag, phi, alpha_surf);
  double JRatfL[6] = {0.}; 

  double *flux_surf_nodal = &flux_surf[0]; 
  double cfl = 0.0; 
  double JfL_quad = 0.0; 
  double JfR_quad = 0.0; 
  double Jfavg_quad = 0.0; 
  double Jfjump_quad = 0.0; 
  cfl = fmax(fabs(alpha_surf[0]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = 0.5477225575051661*JfR[11]-0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]-0.8215838362577489*JfR[7]+0.4743416490252568*JfR[6]+0.6123724356957944*JfR[5]+0.8215838362577489*JfR[4]-0.3535533905932737*JfR[3]-0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_surf[0]*Jfavg_quad - fabs(alpha_surf[0])*Jfjump_quad; 
  cfl = fmax(fabs(alpha_surf[1]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = -(0.6846531968814573*JfR[11])+0.39528470752104733*JfR[10]+0.6846531968814574*JfR[9]-0.3952847075210473*JfR[8]+0.6123724356957944*JfR[5]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_surf[1]*Jfavg_quad - fabs(alpha_surf[1])*Jfjump_quad; 
  cfl = fmax(fabs(alpha_surf[2]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = 0.5477225575051661*JfR[11]-0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]+0.8215838362577489*JfR[7]-0.4743416490252568*JfR[6]+0.6123724356957944*JfR[5]-0.8215838362577489*JfR[4]-0.3535533905932737*JfR[3]+0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_surf[2]*Jfavg_quad - fabs(alpha_surf[2])*Jfjump_quad; 
  cfl = fmax(fabs(alpha_surf[3]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = -(0.5477225575051661*JfR[11])+0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]+0.8215838362577489*JfR[7]-0.4743416490252568*JfR[6]-0.6123724356957944*JfR[5]+0.8215838362577489*JfR[4]+0.3535533905932737*JfR[3]-0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[3] = alpha_surf[3]*Jfavg_quad - fabs(alpha_surf[3])*Jfjump_quad; 
  cfl = fmax(fabs(alpha_surf[4]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = 0.6846531968814573*JfR[11]-0.39528470752104733*JfR[10]+0.6846531968814574*JfR[9]-0.3952847075210473*JfR[8]-0.6123724356957944*JfR[5]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[4] = alpha_surf[4]*Jfavg_quad - fabs(alpha_surf[4])*Jfjump_quad; 
  cfl = fmax(fabs(alpha_surf[5]), fabs(cfl)); 
  JfL_quad = 0.0; 
  JfR_quad = -(0.5477225575051661*JfR[11])+0.3162277660168379*JfR[10]-0.5477225575051661*JfR[9]+0.3162277660168379*JfR[8]-0.8215838362577489*JfR[7]+0.4743416490252568*JfR[6]-0.6123724356957944*JfR[5]-0.8215838362577489*JfR[4]+0.3535533905932737*JfR[3]+0.4743416490252568*JfR[2]-0.6123724356957944*JfR[1]+0.3535533905932737*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[5] = alpha_surf[5]*Jfavg_quad - fabs(alpha_surf[5])*Jfjump_quad; 

  return cfl*1.5*rdx2; 

} 
