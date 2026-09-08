#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_lower_periodic_2x2v_ser_p1(
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = (1.4142135623730951*phiR[0]-2.4494897427831783*phiR[1])*q_+vmapSq[0]*m_+(bmag[0]-1.7320508075688772*bmag[1])*vmap[2]; 
  hamil[1] = (1.4142135623730951*phiR[2]-2.4494897427831783*phiR[3])*q_+vmap[2]*(bmag[2]-1.7320508075688772*bmag[3]); 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = (bmag[0]-1.7320508075688772*bmag[1])*vmap[3]; 
  hamil[5] = (bmag[2]-1.7320508075688772*bmag[3])*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double *ypotderL = &yfieldL[0]; 
  double *ypotderR = &yfieldR[0]; 
  ypotderR[1] += -(3.4641016151377544*phiR[0]*q_*rdx2)-2.4494897427831783*bmag[0]*vmap[2]*rdx2; 
  ypotderR[5] += -(3.4641016151377544*phiR[2]*q_*rdx2)-2.4494897427831783*bmag[2]*vmap[2]*rdx2; 
  ypotderR[8] += -(2.4494897427831783*bmag[0]*vmap[3]*rdx2); 
  ypotderR[12] += -(2.4494897427831783*bmag[2]*vmap[3]*rdx2); 
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
  double ypot_der_surf_n[12] = {0.0}; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1]; 
  mvpar_quad[1] = (0.49999999999999994*hamil[2])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1]; 
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


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]+0.38729833462074165*(JfR[21]+JfR[20])-0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]-0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[0] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]+0.4841229182759271*(JfL[21]+JfL[20])+0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]-0.4330127018922193*(JfL[8]+JfL[5])-0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]-0.4841229182759271*(JfR[21]+JfR[20])+0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]+0.4330127018922193*(JfR[8]+JfR[5])-0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[1] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]+0.38729833462074165*(JfR[21]+JfR[20])-0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]+0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[2] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]-0.38729833462074165*JfR[21]+0.38729833462074165*JfR[20]+0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[3] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]-0.4841229182759271*JfL[21]+0.4841229182759271*JfL[20]-0.2795084971874738*JfL[19]+0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]+0.4841229182759271*JfR[21]-0.4841229182759271*JfR[20]-0.2795084971874738*JfR[19]+0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[4] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]-0.38729833462074165*JfR[21]+0.38729833462074165*JfR[20]+0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]+0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2])-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[5] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  bmag_quad = gkdgs[1].bmag; 
  Jc_quad = gkdgs[1].Jc; 
  B3_quad = gkdgs[1].B3; 
  normcurlbhat_quad = gkdgs[1].normcurlbhat; 
  bhat_quad[0] = gkdgs[1].bhat.x[0]; 
  bhat_quad[1] = gkdgs[1].bhat.x[1]; 
  bhat_quad[2] = gkdgs[1].bhat.x[2]; 
  area_elem_quad = dgs[1].area_elem; 


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]+0.38729833462074165*JfR[21]-0.38729833462074165*JfR[20]-0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]-0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[6] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]+0.4841229182759271*JfL[21]-0.4841229182759271*JfL[20]+0.2795084971874738*JfL[19]-0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]-0.4841229182759271*JfR[21]+0.4841229182759271*JfR[20]+0.2795084971874738*JfR[19]-0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[7] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]+0.38729833462074165*JfR[21]-0.38729833462074165*JfR[20]-0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[8] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])-0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[8] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]-0.38729833462074165*(JfR[21]+JfR[20])+0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]-0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[9] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[9] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]-0.4841229182759271*(JfL[21]+JfL[20])-0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]+0.4330127018922193*(JfL[8]+JfL[5])+0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]+0.4841229182759271*(JfR[21]+JfR[20])-0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]-0.4330127018922193*(JfR[8]+JfR[5])+0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[10] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[10] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]-0.38729833462074165*(JfR[21]+JfR[20])+0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]+0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[11] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiL[3]*q_)-2.0*phiL[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiL[1]*q_)-2.0*phiL[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(2.0*phiR[3]*q_+1.4142135623730951*vmap[2]*bmag[3])-2.0*phiR[2]*q_-1.4142135623730951*bmag[2]*vmap[2]))-0.25*(1.7320508075688772*(2.0*phiR[1]*q_+1.4142135623730951*bmag[1]*vmap[2])-2.0*phiR[0]*q_-1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(bmag[2]*vmap[3]-1.7320508075688772*bmag[3]*vmap[3])+0.3535533905932737*(bmag[0]*vmap[3]-1.7320508075688772*bmag[1]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[11] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  ypotderR[0] += -(0.1388888888888889*ypot_der_surf_n[11]*rdx2)-0.2222222222222222*ypot_der_surf_n[10]*rdx2-0.1388888888888889*ypot_der_surf_n[9]*rdx2-0.1388888888888889*ypot_der_surf_n[8]*rdx2-0.2222222222222222*ypot_der_surf_n[7]*rdx2-0.1388888888888889*ypot_der_surf_n[6]*rdx2-0.1388888888888889*ypot_der_surf_n[5]*rdx2-0.2222222222222222*ypot_der_surf_n[4]*rdx2-0.1388888888888889*ypot_der_surf_n[3]*rdx2-0.1388888888888889*ypot_der_surf_n[2]*rdx2-0.2222222222222222*ypot_der_surf_n[1]*rdx2-0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderR[1] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2+0.24056261216234412*ypot_der_surf_n[8]*rdx2+0.3849001794597506*ypot_der_surf_n[7]*rdx2+0.24056261216234412*ypot_der_surf_n[6]*rdx2+0.24056261216234412*ypot_der_surf_n[5]*rdx2+0.3849001794597506*ypot_der_surf_n[4]*rdx2+0.24056261216234412*ypot_der_surf_n[3]*rdx2+0.24056261216234412*ypot_der_surf_n[2]*rdx2+0.3849001794597506*ypot_der_surf_n[1]*rdx2+0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderR[2] += -(0.1388888888888889*ypot_der_surf_n[11]*rdx2)-0.2222222222222222*ypot_der_surf_n[10]*rdx2-0.1388888888888889*ypot_der_surf_n[9]*rdx2-0.1388888888888889*ypot_der_surf_n[8]*rdx2-0.2222222222222222*ypot_der_surf_n[7]*rdx2-0.1388888888888889*ypot_der_surf_n[6]*rdx2+0.1388888888888889*ypot_der_surf_n[5]*rdx2+0.2222222222222222*ypot_der_surf_n[4]*rdx2+0.1388888888888889*ypot_der_surf_n[3]*rdx2+0.1388888888888889*ypot_der_surf_n[2]*rdx2+0.2222222222222222*ypot_der_surf_n[1]*rdx2+0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderR[3] += -(0.18633899812498247*ypot_der_surf_n[11]*rdx2)+0.18633899812498247*ypot_der_surf_n[9]*rdx2-0.18633899812498247*ypot_der_surf_n[8]*rdx2+0.18633899812498247*ypot_der_surf_n[6]*rdx2-0.18633899812498247*ypot_der_surf_n[5]*rdx2+0.18633899812498247*ypot_der_surf_n[3]*rdx2-0.18633899812498247*ypot_der_surf_n[2]*rdx2+0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderR[4] += -(0.1388888888888889*ypot_der_surf_n[11]*rdx2)-0.2222222222222222*ypot_der_surf_n[10]*rdx2-0.1388888888888889*ypot_der_surf_n[9]*rdx2+0.1388888888888889*ypot_der_surf_n[8]*rdx2+0.2222222222222222*ypot_der_surf_n[7]*rdx2+0.1388888888888889*ypot_der_surf_n[6]*rdx2-0.1388888888888889*ypot_der_surf_n[5]*rdx2-0.2222222222222222*ypot_der_surf_n[4]*rdx2-0.1388888888888889*ypot_der_surf_n[3]*rdx2+0.1388888888888889*ypot_der_surf_n[2]*rdx2+0.2222222222222222*ypot_der_surf_n[1]*rdx2+0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderR[5] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2+0.24056261216234412*ypot_der_surf_n[8]*rdx2+0.3849001794597506*ypot_der_surf_n[7]*rdx2+0.24056261216234412*ypot_der_surf_n[6]*rdx2-0.24056261216234412*ypot_der_surf_n[5]*rdx2-0.3849001794597506*ypot_der_surf_n[4]*rdx2-0.24056261216234412*ypot_der_surf_n[3]*rdx2-0.24056261216234412*ypot_der_surf_n[2]*rdx2-0.3849001794597506*ypot_der_surf_n[1]*rdx2-0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderR[6] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2+0.32274861218395146*ypot_der_surf_n[8]*rdx2-0.32274861218395146*ypot_der_surf_n[6]*rdx2+0.32274861218395146*ypot_der_surf_n[5]*rdx2-0.32274861218395146*ypot_der_surf_n[3]*rdx2+0.32274861218395146*ypot_der_surf_n[2]*rdx2-0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderR[7] += -(0.18633899812498247*ypot_der_surf_n[11]*rdx2)+0.18633899812498247*ypot_der_surf_n[9]*rdx2-0.18633899812498247*ypot_der_surf_n[8]*rdx2+0.18633899812498247*ypot_der_surf_n[6]*rdx2+0.18633899812498247*ypot_der_surf_n[5]*rdx2-0.18633899812498247*ypot_der_surf_n[3]*rdx2+0.18633899812498247*ypot_der_surf_n[2]*rdx2-0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderR[8] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2-0.24056261216234412*ypot_der_surf_n[8]*rdx2-0.3849001794597506*ypot_der_surf_n[7]*rdx2-0.24056261216234412*ypot_der_surf_n[6]*rdx2+0.24056261216234412*ypot_der_surf_n[5]*rdx2+0.3849001794597506*ypot_der_surf_n[4]*rdx2+0.24056261216234412*ypot_der_surf_n[3]*rdx2-0.24056261216234412*ypot_der_surf_n[2]*rdx2-0.3849001794597506*ypot_der_surf_n[1]*rdx2-0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderR[9] += -(0.1388888888888889*ypot_der_surf_n[11]*rdx2)-0.2222222222222222*ypot_der_surf_n[10]*rdx2-0.1388888888888889*ypot_der_surf_n[9]*rdx2+0.1388888888888889*ypot_der_surf_n[8]*rdx2+0.2222222222222222*ypot_der_surf_n[7]*rdx2+0.1388888888888889*ypot_der_surf_n[6]*rdx2+0.1388888888888889*ypot_der_surf_n[5]*rdx2+0.2222222222222222*ypot_der_surf_n[4]*rdx2+0.1388888888888889*ypot_der_surf_n[3]*rdx2-0.1388888888888889*ypot_der_surf_n[2]*rdx2-0.2222222222222222*ypot_der_surf_n[1]*rdx2-0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderR[10] += -(0.18633899812498247*ypot_der_surf_n[11]*rdx2)+0.18633899812498247*ypot_der_surf_n[9]*rdx2+0.18633899812498247*ypot_der_surf_n[8]*rdx2-0.18633899812498247*ypot_der_surf_n[6]*rdx2-0.18633899812498247*ypot_der_surf_n[5]*rdx2+0.18633899812498247*ypot_der_surf_n[3]*rdx2+0.18633899812498247*ypot_der_surf_n[2]*rdx2-0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderR[11] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2+0.32274861218395146*ypot_der_surf_n[8]*rdx2-0.32274861218395146*ypot_der_surf_n[6]*rdx2-0.32274861218395146*ypot_der_surf_n[5]*rdx2+0.32274861218395146*ypot_der_surf_n[3]*rdx2-0.32274861218395146*ypot_der_surf_n[2]*rdx2+0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderR[12] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2-0.24056261216234412*ypot_der_surf_n[8]*rdx2-0.3849001794597506*ypot_der_surf_n[7]*rdx2-0.24056261216234412*ypot_der_surf_n[6]*rdx2-0.24056261216234412*ypot_der_surf_n[5]*rdx2-0.3849001794597506*ypot_der_surf_n[4]*rdx2-0.24056261216234412*ypot_der_surf_n[3]*rdx2+0.24056261216234412*ypot_der_surf_n[2]*rdx2+0.3849001794597506*ypot_der_surf_n[1]*rdx2+0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderR[13] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2-0.32274861218395146*ypot_der_surf_n[8]*rdx2+0.32274861218395146*ypot_der_surf_n[6]*rdx2+0.32274861218395146*ypot_der_surf_n[5]*rdx2-0.32274861218395146*ypot_der_surf_n[3]*rdx2-0.32274861218395146*ypot_der_surf_n[2]*rdx2+0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderR[14] += -(0.18633899812498247*ypot_der_surf_n[11]*rdx2)+0.18633899812498247*ypot_der_surf_n[9]*rdx2+0.18633899812498247*ypot_der_surf_n[8]*rdx2-0.18633899812498247*ypot_der_surf_n[6]*rdx2+0.18633899812498247*ypot_der_surf_n[5]*rdx2-0.18633899812498247*ypot_der_surf_n[3]*rdx2-0.18633899812498247*ypot_der_surf_n[2]*rdx2+0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderR[15] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2-0.32274861218395146*ypot_der_surf_n[8]*rdx2+0.32274861218395146*ypot_der_surf_n[6]*rdx2-0.32274861218395146*ypot_der_surf_n[5]*rdx2+0.32274861218395146*ypot_der_surf_n[3]*rdx2+0.32274861218395146*ypot_der_surf_n[2]*rdx2-0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderR[16] += -(0.12422599874998831*ypot_der_surf_n[11]*rdx2)+0.24845199749997662*ypot_der_surf_n[10]*rdx2-0.12422599874998831*ypot_der_surf_n[9]*rdx2-0.12422599874998831*ypot_der_surf_n[8]*rdx2+0.24845199749997662*ypot_der_surf_n[7]*rdx2-0.12422599874998831*ypot_der_surf_n[6]*rdx2-0.12422599874998831*ypot_der_surf_n[5]*rdx2+0.24845199749997662*ypot_der_surf_n[4]*rdx2-0.12422599874998831*ypot_der_surf_n[3]*rdx2-0.12422599874998831*ypot_der_surf_n[2]*rdx2+0.24845199749997662*ypot_der_surf_n[1]*rdx2-0.12422599874998831*ypot_der_surf_n[0]*rdx2; 
  ypotderR[17] += 0.2151657414559676*ypot_der_surf_n[11]*rdx2-0.4303314829119352*ypot_der_surf_n[10]*rdx2+0.2151657414559676*ypot_der_surf_n[9]*rdx2+0.2151657414559676*ypot_der_surf_n[8]*rdx2-0.4303314829119352*ypot_der_surf_n[7]*rdx2+0.2151657414559676*ypot_der_surf_n[6]*rdx2+0.2151657414559676*ypot_der_surf_n[5]*rdx2-0.4303314829119352*ypot_der_surf_n[4]*rdx2+0.2151657414559676*ypot_der_surf_n[3]*rdx2+0.2151657414559676*ypot_der_surf_n[2]*rdx2-0.4303314829119352*ypot_der_surf_n[1]*rdx2+0.2151657414559676*ypot_der_surf_n[0]*rdx2; 
  ypotderR[18] += -(0.12422599874998837*ypot_der_surf_n[11]*rdx2)+0.24845199749997673*ypot_der_surf_n[10]*rdx2-0.12422599874998837*ypot_der_surf_n[9]*rdx2-0.12422599874998837*ypot_der_surf_n[8]*rdx2+0.24845199749997673*ypot_der_surf_n[7]*rdx2-0.12422599874998837*ypot_der_surf_n[6]*rdx2+0.12422599874998837*ypot_der_surf_n[5]*rdx2-0.24845199749997673*ypot_der_surf_n[4]*rdx2+0.12422599874998837*ypot_der_surf_n[3]*rdx2+0.12422599874998837*ypot_der_surf_n[2]*rdx2-0.24845199749997673*ypot_der_surf_n[1]*rdx2+0.12422599874998837*ypot_der_surf_n[0]*rdx2; 
  ypotderR[19] += -(0.12422599874998837*ypot_der_surf_n[11]*rdx2)+0.24845199749997673*ypot_der_surf_n[10]*rdx2-0.12422599874998837*ypot_der_surf_n[9]*rdx2+0.12422599874998837*ypot_der_surf_n[8]*rdx2-0.24845199749997673*ypot_der_surf_n[7]*rdx2+0.12422599874998837*ypot_der_surf_n[6]*rdx2-0.12422599874998837*ypot_der_surf_n[5]*rdx2+0.24845199749997673*ypot_der_surf_n[4]*rdx2-0.12422599874998837*ypot_der_surf_n[3]*rdx2+0.12422599874998837*ypot_der_surf_n[2]*rdx2-0.24845199749997673*ypot_der_surf_n[1]*rdx2+0.12422599874998837*ypot_der_surf_n[0]*rdx2; 
  ypotderR[20] += 0.21516574145596765*ypot_der_surf_n[11]*rdx2-0.4303314829119353*ypot_der_surf_n[10]*rdx2+0.21516574145596765*ypot_der_surf_n[9]*rdx2+0.21516574145596765*ypot_der_surf_n[8]*rdx2-0.4303314829119353*ypot_der_surf_n[7]*rdx2+0.21516574145596765*ypot_der_surf_n[6]*rdx2-0.21516574145596765*ypot_der_surf_n[5]*rdx2+0.4303314829119353*ypot_der_surf_n[4]*rdx2-0.21516574145596765*ypot_der_surf_n[3]*rdx2-0.21516574145596765*ypot_der_surf_n[2]*rdx2+0.4303314829119353*ypot_der_surf_n[1]*rdx2-0.21516574145596765*ypot_der_surf_n[0]*rdx2; 
  ypotderR[21] += 0.21516574145596765*ypot_der_surf_n[11]*rdx2-0.4303314829119353*ypot_der_surf_n[10]*rdx2+0.21516574145596765*ypot_der_surf_n[9]*rdx2-0.21516574145596765*ypot_der_surf_n[8]*rdx2+0.4303314829119353*ypot_der_surf_n[7]*rdx2-0.21516574145596765*ypot_der_surf_n[6]*rdx2+0.21516574145596765*ypot_der_surf_n[5]*rdx2-0.4303314829119353*ypot_der_surf_n[4]*rdx2+0.21516574145596765*ypot_der_surf_n[3]*rdx2-0.21516574145596765*ypot_der_surf_n[2]*rdx2+0.4303314829119353*ypot_der_surf_n[1]*rdx2-0.21516574145596765*ypot_der_surf_n[0]*rdx2; 
  ypotderR[22] += -(0.12422599874998831*ypot_der_surf_n[11]*rdx2)+0.24845199749997662*ypot_der_surf_n[10]*rdx2-0.12422599874998831*ypot_der_surf_n[9]*rdx2+0.12422599874998831*ypot_der_surf_n[8]*rdx2-0.24845199749997662*ypot_der_surf_n[7]*rdx2+0.12422599874998831*ypot_der_surf_n[6]*rdx2+0.12422599874998831*ypot_der_surf_n[5]*rdx2-0.24845199749997662*ypot_der_surf_n[4]*rdx2+0.12422599874998831*ypot_der_surf_n[3]*rdx2-0.12422599874998831*ypot_der_surf_n[2]*rdx2+0.24845199749997662*ypot_der_surf_n[1]*rdx2-0.12422599874998831*ypot_der_surf_n[0]*rdx2; 
  ypotderR[23] += 0.2151657414559676*ypot_der_surf_n[11]*rdx2-0.4303314829119352*ypot_der_surf_n[10]*rdx2+0.2151657414559676*ypot_der_surf_n[9]*rdx2-0.2151657414559676*ypot_der_surf_n[8]*rdx2+0.4303314829119352*ypot_der_surf_n[7]*rdx2-0.2151657414559676*ypot_der_surf_n[6]*rdx2-0.2151657414559676*ypot_der_surf_n[5]*rdx2+0.4303314829119352*ypot_der_surf_n[4]*rdx2-0.2151657414559676*ypot_der_surf_n[3]*rdx2+0.2151657414559676*ypot_der_surf_n[2]*rdx2-0.4303314829119352*ypot_der_surf_n[1]*rdx2+0.2151657414559676*ypot_der_surf_n[0]*rdx2; 


  return cfl*1.5*rdx2; 

} 
GKYL_CU_DH double gk_collisionless_flux_boundary_surfx_upper_periodic_2x2v_ser_p1(
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = (2.4494897427831783*phiL[1]+1.4142135623730951*phiL[0])*q_+vmapSq[0]*m_+(1.7320508075688772*bmag[1]+bmag[0])*vmap[2]; 
  hamil[1] = (2.4494897427831783*phiL[3]+1.4142135623730951*phiL[2])*q_+vmap[2]*(1.7320508075688772*bmag[3]+bmag[2]); 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = (1.7320508075688772*bmag[1]+bmag[0])*vmap[3]; 
  hamil[5] = (1.7320508075688772*bmag[3]+bmag[2])*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

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
  double ypot_der_surf_n[12] = {0.0}; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1]; 
  mvpar_quad[1] = (0.49999999999999994*hamil[2])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1]; 
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


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]+0.38729833462074165*(JfR[21]+JfR[20])-0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]-0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[0] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]+0.4841229182759271*(JfL[21]+JfL[20])+0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]-0.4330127018922193*(JfL[8]+JfL[5])-0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]-0.4841229182759271*(JfR[21]+JfR[20])+0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]+0.4330127018922193*(JfR[8]+JfR[5])-0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[1] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]+0.38729833462074165*(JfR[21]+JfR[20])-0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]+0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[2] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]-0.38729833462074165*JfR[21]+0.38729833462074165*JfR[20]+0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[3] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]-0.4841229182759271*JfL[21]+0.4841229182759271*JfL[20]-0.2795084971874738*JfL[19]+0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]+0.4841229182759271*JfR[21]-0.4841229182759271*JfR[20]-0.2795084971874738*JfR[19]+0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[4] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]-0.38729833462074165*JfR[21]+0.38729833462074165*JfR[20]+0.22360679774997894*JfR[19]-0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]+0.33541019662496846*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = -(0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = -(0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2]))+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[5] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  bmag_quad = gkdgs[1].bmag; 
  Jc_quad = gkdgs[1].Jc; 
  B3_quad = gkdgs[1].B3; 
  normcurlbhat_quad = gkdgs[1].normcurlbhat; 
  bhat_quad[0] = gkdgs[1].bhat.x[0]; 
  bhat_quad[1] = gkdgs[1].bhat.x[1]; 
  bhat_quad[2] = gkdgs[1].bhat.x[2]; 
  area_elem_quad = dgs[1].area_elem; 


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]+0.38729833462074165*JfR[21]-0.38729833462074165*JfR[20]-0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]-0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[6] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]+0.4841229182759271*JfL[21]-0.4841229182759271*JfL[20]+0.2795084971874738*JfL[19]-0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]-0.4841229182759271*JfR[21]+0.4841229182759271*JfR[20]+0.2795084971874738*JfR[19]-0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[7] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3872983346207417*JfR[23]-0.22360679774997896*JfR[22]+0.38729833462074165*JfR[21]-0.38729833462074165*JfR[20]-0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[8] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])-0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])-0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[8] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]-0.38729833462074165*(JfR[21]+JfR[20])+0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]+0.5809475019311124*JfR[15]-0.33541019662496846*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.33541019662496846*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]-0.33541019662496846*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]-0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[9] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[9] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]-0.4841229182759271*(JfL[21]+JfL[20])-0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]+0.4330127018922193*(JfL[8]+JfL[5])+0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]+0.4841229182759271*(JfR[21]+JfR[20])-0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]-0.4330127018922193*(JfR[8]+JfR[5])+0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[10] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[10] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3872983346207417*JfR[23])+0.22360679774997896*JfR[22]-0.38729833462074165*(JfR[21]+JfR[20])+0.22360679774997894*JfR[19]+0.22360679774997902*JfR[18]-0.3872983346207417*JfR[17]+0.22360679774997896*JfR[16]-0.5809475019311124*JfR[15]+0.33541019662496846*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.33541019662496846*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]+0.33541019662496846*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]+0.33541019662496846*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[11] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 
  ypotL_quad = 0.25*(1.7320508075688772*(2.0*phiL[3]*q_+1.4142135623730951*vmap[2]*bmag[3])+2.0*phiL[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(2.0*phiL[1]*q_+1.4142135623730951*bmag[1]*vmap[2])+2.0*phiL[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypotR_quad = 0.25*(1.7320508075688772*(1.4142135623730951*vmap[2]*bmag[3]-2.0*phiR[3]*q_)+2.0*phiR[2]*q_+1.4142135623730951*bmag[2]*vmap[2])+0.25*(1.7320508075688772*(1.4142135623730951*bmag[1]*vmap[2]-2.0*phiR[1]*q_)+2.0*phiR[0]*q_+1.4142135623730951*bmag[0]*vmap[2])+0.3535533905932737*(1.7320508075688772*bmag[3]*vmap[3]+bmag[2]*vmap[3])+0.3535533905932737*(1.7320508075688772*bmag[1]*vmap[3]+bmag[0]*vmap[3]); 
  ypot_der_avg_quad = (ypotL_quad + ypotR_quad)/2.0; 
  ypot_der_jump_quad = (ypotR_quad - ypotL_quad)/2.0; 
  ypot_der_surf_n[11] = ypot_der_avg_quad + gsign(alpha_quad)*ypot_der_jump_quad; 

  ypotderL[0] += 0.1388888888888889*ypot_der_surf_n[11]*rdx2+0.2222222222222222*ypot_der_surf_n[10]*rdx2+0.1388888888888889*ypot_der_surf_n[9]*rdx2+0.1388888888888889*ypot_der_surf_n[8]*rdx2+0.2222222222222222*ypot_der_surf_n[7]*rdx2+0.1388888888888889*ypot_der_surf_n[6]*rdx2+0.1388888888888889*ypot_der_surf_n[5]*rdx2+0.2222222222222222*ypot_der_surf_n[4]*rdx2+0.1388888888888889*ypot_der_surf_n[3]*rdx2+0.1388888888888889*ypot_der_surf_n[2]*rdx2+0.2222222222222222*ypot_der_surf_n[1]*rdx2+0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderL[1] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2+0.24056261216234412*ypot_der_surf_n[8]*rdx2+0.3849001794597506*ypot_der_surf_n[7]*rdx2+0.24056261216234412*ypot_der_surf_n[6]*rdx2+0.24056261216234412*ypot_der_surf_n[5]*rdx2+0.3849001794597506*ypot_der_surf_n[4]*rdx2+0.24056261216234412*ypot_der_surf_n[3]*rdx2+0.24056261216234412*ypot_der_surf_n[2]*rdx2+0.3849001794597506*ypot_der_surf_n[1]*rdx2+0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderL[2] += 0.1388888888888889*ypot_der_surf_n[11]*rdx2+0.2222222222222222*ypot_der_surf_n[10]*rdx2+0.1388888888888889*ypot_der_surf_n[9]*rdx2+0.1388888888888889*ypot_der_surf_n[8]*rdx2+0.2222222222222222*ypot_der_surf_n[7]*rdx2+0.1388888888888889*ypot_der_surf_n[6]*rdx2-0.1388888888888889*ypot_der_surf_n[5]*rdx2-0.2222222222222222*ypot_der_surf_n[4]*rdx2-0.1388888888888889*ypot_der_surf_n[3]*rdx2-0.1388888888888889*ypot_der_surf_n[2]*rdx2-0.2222222222222222*ypot_der_surf_n[1]*rdx2-0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderL[3] += 0.18633899812498247*ypot_der_surf_n[11]*rdx2-0.18633899812498247*ypot_der_surf_n[9]*rdx2+0.18633899812498247*ypot_der_surf_n[8]*rdx2-0.18633899812498247*ypot_der_surf_n[6]*rdx2+0.18633899812498247*ypot_der_surf_n[5]*rdx2-0.18633899812498247*ypot_der_surf_n[3]*rdx2+0.18633899812498247*ypot_der_surf_n[2]*rdx2-0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderL[4] += 0.1388888888888889*ypot_der_surf_n[11]*rdx2+0.2222222222222222*ypot_der_surf_n[10]*rdx2+0.1388888888888889*ypot_der_surf_n[9]*rdx2-0.1388888888888889*ypot_der_surf_n[8]*rdx2-0.2222222222222222*ypot_der_surf_n[7]*rdx2-0.1388888888888889*ypot_der_surf_n[6]*rdx2+0.1388888888888889*ypot_der_surf_n[5]*rdx2+0.2222222222222222*ypot_der_surf_n[4]*rdx2+0.1388888888888889*ypot_der_surf_n[3]*rdx2-0.1388888888888889*ypot_der_surf_n[2]*rdx2-0.2222222222222222*ypot_der_surf_n[1]*rdx2-0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderL[5] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2+0.24056261216234412*ypot_der_surf_n[8]*rdx2+0.3849001794597506*ypot_der_surf_n[7]*rdx2+0.24056261216234412*ypot_der_surf_n[6]*rdx2-0.24056261216234412*ypot_der_surf_n[5]*rdx2-0.3849001794597506*ypot_der_surf_n[4]*rdx2-0.24056261216234412*ypot_der_surf_n[3]*rdx2-0.24056261216234412*ypot_der_surf_n[2]*rdx2-0.3849001794597506*ypot_der_surf_n[1]*rdx2-0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderL[6] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2+0.32274861218395146*ypot_der_surf_n[8]*rdx2-0.32274861218395146*ypot_der_surf_n[6]*rdx2+0.32274861218395146*ypot_der_surf_n[5]*rdx2-0.32274861218395146*ypot_der_surf_n[3]*rdx2+0.32274861218395146*ypot_der_surf_n[2]*rdx2-0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderL[7] += 0.18633899812498247*ypot_der_surf_n[11]*rdx2-0.18633899812498247*ypot_der_surf_n[9]*rdx2+0.18633899812498247*ypot_der_surf_n[8]*rdx2-0.18633899812498247*ypot_der_surf_n[6]*rdx2-0.18633899812498247*ypot_der_surf_n[5]*rdx2+0.18633899812498247*ypot_der_surf_n[3]*rdx2-0.18633899812498247*ypot_der_surf_n[2]*rdx2+0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderL[8] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2-0.24056261216234412*ypot_der_surf_n[8]*rdx2-0.3849001794597506*ypot_der_surf_n[7]*rdx2-0.24056261216234412*ypot_der_surf_n[6]*rdx2+0.24056261216234412*ypot_der_surf_n[5]*rdx2+0.3849001794597506*ypot_der_surf_n[4]*rdx2+0.24056261216234412*ypot_der_surf_n[3]*rdx2-0.24056261216234412*ypot_der_surf_n[2]*rdx2-0.3849001794597506*ypot_der_surf_n[1]*rdx2-0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderL[9] += 0.1388888888888889*ypot_der_surf_n[11]*rdx2+0.2222222222222222*ypot_der_surf_n[10]*rdx2+0.1388888888888889*ypot_der_surf_n[9]*rdx2-0.1388888888888889*ypot_der_surf_n[8]*rdx2-0.2222222222222222*ypot_der_surf_n[7]*rdx2-0.1388888888888889*ypot_der_surf_n[6]*rdx2-0.1388888888888889*ypot_der_surf_n[5]*rdx2-0.2222222222222222*ypot_der_surf_n[4]*rdx2-0.1388888888888889*ypot_der_surf_n[3]*rdx2+0.1388888888888889*ypot_der_surf_n[2]*rdx2+0.2222222222222222*ypot_der_surf_n[1]*rdx2+0.1388888888888889*ypot_der_surf_n[0]*rdx2; 
  ypotderL[10] += 0.18633899812498247*ypot_der_surf_n[11]*rdx2-0.18633899812498247*ypot_der_surf_n[9]*rdx2-0.18633899812498247*ypot_der_surf_n[8]*rdx2+0.18633899812498247*ypot_der_surf_n[6]*rdx2+0.18633899812498247*ypot_der_surf_n[5]*rdx2-0.18633899812498247*ypot_der_surf_n[3]*rdx2-0.18633899812498247*ypot_der_surf_n[2]*rdx2+0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderL[11] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2+0.32274861218395146*ypot_der_surf_n[8]*rdx2-0.32274861218395146*ypot_der_surf_n[6]*rdx2-0.32274861218395146*ypot_der_surf_n[5]*rdx2+0.32274861218395146*ypot_der_surf_n[3]*rdx2-0.32274861218395146*ypot_der_surf_n[2]*rdx2+0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderL[12] += 0.24056261216234412*ypot_der_surf_n[11]*rdx2+0.3849001794597506*ypot_der_surf_n[10]*rdx2+0.24056261216234412*ypot_der_surf_n[9]*rdx2-0.24056261216234412*ypot_der_surf_n[8]*rdx2-0.3849001794597506*ypot_der_surf_n[7]*rdx2-0.24056261216234412*ypot_der_surf_n[6]*rdx2-0.24056261216234412*ypot_der_surf_n[5]*rdx2-0.3849001794597506*ypot_der_surf_n[4]*rdx2-0.24056261216234412*ypot_der_surf_n[3]*rdx2+0.24056261216234412*ypot_der_surf_n[2]*rdx2+0.3849001794597506*ypot_der_surf_n[1]*rdx2+0.24056261216234412*ypot_der_surf_n[0]*rdx2; 
  ypotderL[13] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2-0.32274861218395146*ypot_der_surf_n[8]*rdx2+0.32274861218395146*ypot_der_surf_n[6]*rdx2+0.32274861218395146*ypot_der_surf_n[5]*rdx2-0.32274861218395146*ypot_der_surf_n[3]*rdx2-0.32274861218395146*ypot_der_surf_n[2]*rdx2+0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderL[14] += 0.18633899812498247*ypot_der_surf_n[11]*rdx2-0.18633899812498247*ypot_der_surf_n[9]*rdx2-0.18633899812498247*ypot_der_surf_n[8]*rdx2+0.18633899812498247*ypot_der_surf_n[6]*rdx2-0.18633899812498247*ypot_der_surf_n[5]*rdx2+0.18633899812498247*ypot_der_surf_n[3]*rdx2+0.18633899812498247*ypot_der_surf_n[2]*rdx2-0.18633899812498247*ypot_der_surf_n[0]*rdx2; 
  ypotderL[15] += 0.32274861218395146*ypot_der_surf_n[11]*rdx2-0.32274861218395146*ypot_der_surf_n[9]*rdx2-0.32274861218395146*ypot_der_surf_n[8]*rdx2+0.32274861218395146*ypot_der_surf_n[6]*rdx2-0.32274861218395146*ypot_der_surf_n[5]*rdx2+0.32274861218395146*ypot_der_surf_n[3]*rdx2+0.32274861218395146*ypot_der_surf_n[2]*rdx2-0.32274861218395146*ypot_der_surf_n[0]*rdx2; 
  ypotderL[16] += 0.12422599874998831*ypot_der_surf_n[11]*rdx2-0.24845199749997662*ypot_der_surf_n[10]*rdx2+0.12422599874998831*ypot_der_surf_n[9]*rdx2+0.12422599874998831*ypot_der_surf_n[8]*rdx2-0.24845199749997662*ypot_der_surf_n[7]*rdx2+0.12422599874998831*ypot_der_surf_n[6]*rdx2+0.12422599874998831*ypot_der_surf_n[5]*rdx2-0.24845199749997662*ypot_der_surf_n[4]*rdx2+0.12422599874998831*ypot_der_surf_n[3]*rdx2+0.12422599874998831*ypot_der_surf_n[2]*rdx2-0.24845199749997662*ypot_der_surf_n[1]*rdx2+0.12422599874998831*ypot_der_surf_n[0]*rdx2; 
  ypotderL[17] += 0.2151657414559676*ypot_der_surf_n[11]*rdx2-0.4303314829119352*ypot_der_surf_n[10]*rdx2+0.2151657414559676*ypot_der_surf_n[9]*rdx2+0.2151657414559676*ypot_der_surf_n[8]*rdx2-0.4303314829119352*ypot_der_surf_n[7]*rdx2+0.2151657414559676*ypot_der_surf_n[6]*rdx2+0.2151657414559676*ypot_der_surf_n[5]*rdx2-0.4303314829119352*ypot_der_surf_n[4]*rdx2+0.2151657414559676*ypot_der_surf_n[3]*rdx2+0.2151657414559676*ypot_der_surf_n[2]*rdx2-0.4303314829119352*ypot_der_surf_n[1]*rdx2+0.2151657414559676*ypot_der_surf_n[0]*rdx2; 
  ypotderL[18] += 0.12422599874998837*ypot_der_surf_n[11]*rdx2-0.24845199749997673*ypot_der_surf_n[10]*rdx2+0.12422599874998837*ypot_der_surf_n[9]*rdx2+0.12422599874998837*ypot_der_surf_n[8]*rdx2-0.24845199749997673*ypot_der_surf_n[7]*rdx2+0.12422599874998837*ypot_der_surf_n[6]*rdx2-0.12422599874998837*ypot_der_surf_n[5]*rdx2+0.24845199749997673*ypot_der_surf_n[4]*rdx2-0.12422599874998837*ypot_der_surf_n[3]*rdx2-0.12422599874998837*ypot_der_surf_n[2]*rdx2+0.24845199749997673*ypot_der_surf_n[1]*rdx2-0.12422599874998837*ypot_der_surf_n[0]*rdx2; 
  ypotderL[19] += 0.12422599874998837*ypot_der_surf_n[11]*rdx2-0.24845199749997673*ypot_der_surf_n[10]*rdx2+0.12422599874998837*ypot_der_surf_n[9]*rdx2-0.12422599874998837*ypot_der_surf_n[8]*rdx2+0.24845199749997673*ypot_der_surf_n[7]*rdx2-0.12422599874998837*ypot_der_surf_n[6]*rdx2+0.12422599874998837*ypot_der_surf_n[5]*rdx2-0.24845199749997673*ypot_der_surf_n[4]*rdx2+0.12422599874998837*ypot_der_surf_n[3]*rdx2-0.12422599874998837*ypot_der_surf_n[2]*rdx2+0.24845199749997673*ypot_der_surf_n[1]*rdx2-0.12422599874998837*ypot_der_surf_n[0]*rdx2; 
  ypotderL[20] += 0.21516574145596765*ypot_der_surf_n[11]*rdx2-0.4303314829119353*ypot_der_surf_n[10]*rdx2+0.21516574145596765*ypot_der_surf_n[9]*rdx2+0.21516574145596765*ypot_der_surf_n[8]*rdx2-0.4303314829119353*ypot_der_surf_n[7]*rdx2+0.21516574145596765*ypot_der_surf_n[6]*rdx2-0.21516574145596765*ypot_der_surf_n[5]*rdx2+0.4303314829119353*ypot_der_surf_n[4]*rdx2-0.21516574145596765*ypot_der_surf_n[3]*rdx2-0.21516574145596765*ypot_der_surf_n[2]*rdx2+0.4303314829119353*ypot_der_surf_n[1]*rdx2-0.21516574145596765*ypot_der_surf_n[0]*rdx2; 
  ypotderL[21] += 0.21516574145596765*ypot_der_surf_n[11]*rdx2-0.4303314829119353*ypot_der_surf_n[10]*rdx2+0.21516574145596765*ypot_der_surf_n[9]*rdx2-0.21516574145596765*ypot_der_surf_n[8]*rdx2+0.4303314829119353*ypot_der_surf_n[7]*rdx2-0.21516574145596765*ypot_der_surf_n[6]*rdx2+0.21516574145596765*ypot_der_surf_n[5]*rdx2-0.4303314829119353*ypot_der_surf_n[4]*rdx2+0.21516574145596765*ypot_der_surf_n[3]*rdx2-0.21516574145596765*ypot_der_surf_n[2]*rdx2+0.4303314829119353*ypot_der_surf_n[1]*rdx2-0.21516574145596765*ypot_der_surf_n[0]*rdx2; 
  ypotderL[22] += 0.12422599874998831*ypot_der_surf_n[11]*rdx2-0.24845199749997662*ypot_der_surf_n[10]*rdx2+0.12422599874998831*ypot_der_surf_n[9]*rdx2-0.12422599874998831*ypot_der_surf_n[8]*rdx2+0.24845199749997662*ypot_der_surf_n[7]*rdx2-0.12422599874998831*ypot_der_surf_n[6]*rdx2-0.12422599874998831*ypot_der_surf_n[5]*rdx2+0.24845199749997662*ypot_der_surf_n[4]*rdx2-0.12422599874998831*ypot_der_surf_n[3]*rdx2+0.12422599874998831*ypot_der_surf_n[2]*rdx2-0.24845199749997662*ypot_der_surf_n[1]*rdx2+0.12422599874998831*ypot_der_surf_n[0]*rdx2; 
  ypotderL[23] += 0.2151657414559676*ypot_der_surf_n[11]*rdx2-0.4303314829119352*ypot_der_surf_n[10]*rdx2+0.2151657414559676*ypot_der_surf_n[9]*rdx2-0.2151657414559676*ypot_der_surf_n[8]*rdx2+0.4303314829119352*ypot_der_surf_n[7]*rdx2-0.2151657414559676*ypot_der_surf_n[6]*rdx2-0.2151657414559676*ypot_der_surf_n[5]*rdx2+0.4303314829119352*ypot_der_surf_n[4]*rdx2-0.2151657414559676*ypot_der_surf_n[3]*rdx2+0.2151657414559676*ypot_der_surf_n[2]*rdx2-0.4303314829119352*ypot_der_surf_n[1]*rdx2+0.2151657414559676*ypot_der_surf_n[0]*rdx2; 


  return cfl*1.5*rdx2; 

} 
