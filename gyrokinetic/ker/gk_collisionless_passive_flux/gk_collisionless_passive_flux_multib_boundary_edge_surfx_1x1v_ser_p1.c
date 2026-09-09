#include <gkyl_gk_collisionless_passive_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_passive_flux_multib_boundary_edge_surfx_1x1v_ser_p1(
    const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
    const double *bmag, const double *jacobgeo_rat_surfL, const double *jacobgeo_rat_surfR,
    const double *speedsL, const double *speedsR,
    const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // dgs: surface DG geometry.
  // gkdgs: gyrokinetic surface DG geometry.
  // bmag: bmag represented on the surface.
  // jacobgeo_rat_surfL: Ratio of surface conf-space Jacobians in left cell.
  // jacobgeo_rat_surfR: Ratio of surface conf-space Jacobians in right cell.
  // speedsL: characteristic speeds in left cell.
  // speedsR: characteristic speeds in right cell.
  // JfL: distribution times total jacobian in left cell.
  // JfR: distribution times total jacobian in right cell.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //            Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double JRatfR[3] = {0.}; 

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
  bmag_quad = gkdgs[0].bmag; 
  Jc_quad = gkdgs[0].Jc; 
  B3_quad = gkdgs[0].B3; 
  normcurlbhat_quad = gkdgs[0].normcurlbhat; 
  bhat_quad[0] = gkdgs[0].bhat.x[0]; 
  bhat_quad[1] = gkdgs[0].bhat.x[1]; 
  bhat_quad[2] = gkdgs[0].bhat.x[2]; 
  area_elem_quad = dgs[0].area_elem; 


  alpha_quad = 0.0;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]-1.1618950038622249*JfL[3]-0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = 0.0; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.9682458365518543*JfL[5])-0.5590169943749475*JfL[4]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = 0.0; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]+1.1618950038622249*JfL[3]+0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0]; 
  JfR_quad = 0.0; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 


  return cfl*1.5*rdx2; 

} 
