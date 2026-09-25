#include <gkyl_gk_collisionless_flux_kernels.h> 
GKYL_CU_DH double gk_collisionless_flux_no_by_multib_boundary_edge_surfx_2x2v_ser_p1(
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = (2.4494897427831783*phi[1]+1.4142135623730951*phi[0])*q_+vmapSq[0]*m_+(1.7320508075688772*bmag[1]+bmag[0])*vmap[2]; 
  hamil[1] = (2.4494897427831783*phi[3]+1.4142135623730951*phi[2])*q_+vmap[2]*(1.7320508075688772*bmag[3]+bmag[2]); 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = (1.7320508075688772*bmag[1]+bmag[0])*vmap[3]; 
  hamil[5] = (1.7320508075688772*bmag[3]+bmag[2])*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double JRatfR[12] = {0.}; 
  JRatfR[0] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[1]*JfR[5]+jacobgeo_rat_surfR[0]*JfR[1])-1.0*(jacobgeo_rat_surfR[1]*JfR[2]+JfR[0]*jacobgeo_rat_surfR[0]))); 
  JRatfR[1] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[0]*JfR[5]+JfR[1]*jacobgeo_rat_surfR[1])-1.0*(jacobgeo_rat_surfR[0]*JfR[2]+JfR[0]*jacobgeo_rat_surfR[1]))); 
  JRatfR[2] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[1]*JfR[11]+jacobgeo_rat_surfR[0]*JfR[6])-1.0*(jacobgeo_rat_surfR[1]*JfR[7]+jacobgeo_rat_surfR[0]*JfR[3]))); 
  JRatfR[3] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[1]*JfR[12]+jacobgeo_rat_surfR[0]*JfR[8])-1.0*(jacobgeo_rat_surfR[1]*JfR[9]+jacobgeo_rat_surfR[0]*JfR[4]))); 
  JRatfR[4] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[0]*JfR[11]+jacobgeo_rat_surfR[1]*JfR[6])-1.0*(jacobgeo_rat_surfR[0]*JfR[7]+jacobgeo_rat_surfR[1]*JfR[3]))); 
  JRatfR[5] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[0]*JfR[12]+jacobgeo_rat_surfR[1]*JfR[8])-1.0*(jacobgeo_rat_surfR[0]*JfR[9]+jacobgeo_rat_surfR[1]*JfR[4]))); 
  JRatfR[6] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[1]*JfR[15]+jacobgeo_rat_surfR[0]*JfR[13])-1.0*(jacobgeo_rat_surfR[1]*JfR[14]+jacobgeo_rat_surfR[0]*JfR[10]))); 
  JRatfR[7] = -(0.5*(1.7320508075688772*(jacobgeo_rat_surfR[0]*JfR[15]+jacobgeo_rat_surfR[1]*JfR[13])-1.0*(jacobgeo_rat_surfR[0]*JfR[14]+jacobgeo_rat_surfR[1]*JfR[10]))); 
  JRatfR[8] = 0.03333333333333333*(-(25.980762113533157*jacobgeo_rat_surfR[1]*JfR[20])+8.660254037844387*(1.7320508075688772*jacobgeo_rat_surfR[1]*JfR[18]-3.0*jacobgeo_rat_surfR[0]*JfR[17])+15.0*jacobgeo_rat_surfR[0]*JfR[16]); 
  JRatfR[9] = 0.03333333333333333*(8.660254037844387*(1.7320508075688772*jacobgeo_rat_surfR[1]*JfR[16]-3.0*jacobgeo_rat_surfR[0]*JfR[20])+15.0*jacobgeo_rat_surfR[0]*JfR[18]-25.980762113533157*jacobgeo_rat_surfR[1]*JfR[17]); 
  JRatfR[10] = 0.03333333333333333*(-(25.980762113533157*jacobgeo_rat_surfR[1]*JfR[23])+8.660254037844387*(1.7320508075688772*jacobgeo_rat_surfR[1]*JfR[22]-3.0*jacobgeo_rat_surfR[0]*JfR[21])+15.0*jacobgeo_rat_surfR[0]*JfR[19]); 
  JRatfR[11] = 0.03333333333333333*(8.660254037844387*(1.7320508075688772*jacobgeo_rat_surfR[1]*JfR[19]-3.0*jacobgeo_rat_surfR[0]*JfR[23])+15.0*jacobgeo_rat_surfR[0]*JfR[22]-25.980762113533157*jacobgeo_rat_surfR[1]*JfR[21]); 

  double apar_surf[12] = {0.}; 
  apar_surf[0] = 2.4494897427831783*apar[1]+1.4142135623730951*apar[0]; 
  apar_surf[1] = 2.4494897427831783*apar[3]+1.4142135623730951*apar[2]; 

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
  dHdvpar_quad[0] = (0.8164965809277261*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1]; 
  dHdvpar_quad[1] = (0.49999999999999994*hamil[2])/vmap[1]; 
  dHdvpar_quad[2] = (0.8164965809277261*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1]; 
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


  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3162277660168379*JRatfR[11]-0.31622776601683794*(JRatfR[10]+JRatfR[9])+0.3162277660168379*JRatfR[8]-0.4743416490252568*JRatfR[7]+0.4743416490252568*JRatfR[6]+0.3535533905932737*JRatfR[5]+0.4743416490252568*JRatfR[4]-0.3535533905932737*JRatfR[3]-0.4743416490252568*JRatfR[2]-0.3535533905932737*JRatfR[1]+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]+0.4841229182759271*(JfL[21]+JfL[20])+0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]-0.4330127018922193*(JfL[8]+JfL[5])-0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3952847075210473*JRatfR[11])+0.39528470752104733*(JRatfR[10]+JRatfR[9])-0.3952847075210473*JRatfR[8]+0.3535533905932737*JRatfR[5]-0.3535533905932737*(JRatfR[3]+JRatfR[1])+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]-0.38729833462074165*(JfL[21]+JfL[20])-0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3162277660168379*JRatfR[11]-0.31622776601683794*(JRatfR[10]+JRatfR[9])+0.3162277660168379*JRatfR[8]+0.4743416490252568*JRatfR[7]-0.4743416490252568*JRatfR[6]+0.3535533905932737*JRatfR[5]-0.4743416490252568*JRatfR[4]-0.3535533905932737*JRatfR[3]+0.4743416490252568*JRatfR[2]-0.3535533905932737*JRatfR[1]+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3162277660168379*JRatfR[11])+0.31622776601683794*JRatfR[10]-0.31622776601683794*JRatfR[9]+0.3162277660168379*JRatfR[8]+0.4743416490252568*JRatfR[7]-0.4743416490252568*JRatfR[6]-0.3535533905932737*JRatfR[5]+0.4743416490252568*JRatfR[4]+0.3535533905932737*JRatfR[3]-0.4743416490252568*JRatfR[2]-0.3535533905932737*JRatfR[1]+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]-0.4841229182759271*JfL[21]+0.4841229182759271*JfL[20]-0.2795084971874738*JfL[19]+0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3952847075210473*JRatfR[11]-0.39528470752104733*JRatfR[10]+0.39528470752104733*JRatfR[9]-0.3952847075210473*JRatfR[8]-0.3535533905932737*JRatfR[5]+0.3535533905932737*JRatfR[3]-0.3535533905932737*JRatfR[1]+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]+0.38729833462074165*JfL[21]-0.38729833462074165*JfL[20]+0.22360679774997894*JfL[19]-0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3162277660168379*JRatfR[11])+0.31622776601683794*JRatfR[10]-0.31622776601683794*JRatfR[9]+0.3162277660168379*JRatfR[8]-0.4743416490252568*JRatfR[7]+0.4743416490252568*JRatfR[6]-0.3535533905932737*JRatfR[5]-0.4743416490252568*JRatfR[4]+0.3535533905932737*JRatfR[3]+0.4743416490252568*JRatfR[2]-0.3535533905932737*JRatfR[1]+0.3535533905932737*JRatfR[0]; 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  bmag_quad = gkdgs[1].bmag; 
  Jc_quad = gkdgs[1].Jc; 
  B3_quad = gkdgs[1].B3; 
  normcurlbhat_quad = gkdgs[1].normcurlbhat; 
  bhat_quad[0] = gkdgs[1].bhat.x[0]; 
  bhat_quad[1] = gkdgs[1].bhat.x[1]; 
  bhat_quad[2] = gkdgs[1].bhat.x[2]; 
  area_elem_quad = dgs[1].area_elem; 


  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3162277660168379*JRatfR[11])-0.31622776601683794*JRatfR[10]+0.31622776601683794*JRatfR[9]+0.3162277660168379*JRatfR[8]+0.4743416490252568*(JRatfR[7]+JRatfR[6])-0.3535533905932737*JRatfR[5]-0.4743416490252568*JRatfR[4]-0.3535533905932737*JRatfR[3]-0.4743416490252568*JRatfR[2]+0.3535533905932737*(JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]+0.4841229182759271*JfL[21]-0.4841229182759271*JfL[20]+0.2795084971874738*JfL[19]-0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3952847075210473*JRatfR[11]+0.39528470752104733*JRatfR[10]-0.39528470752104733*JRatfR[9]-0.3952847075210473*JRatfR[8]-0.3535533905932737*(JRatfR[5]+JRatfR[3])+0.3535533905932737*(JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.3872983346207417*JfL[23])-0.22360679774997896*JfL[22]-0.38729833462074165*JfL[21]+0.38729833462074165*JfL[20]-0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3162277660168379*JRatfR[11])-0.31622776601683794*JRatfR[10]+0.31622776601683794*JRatfR[9]+0.3162277660168379*JRatfR[8]-0.4743416490252568*(JRatfR[7]+JRatfR[6])-0.3535533905932737*JRatfR[5]+0.4743416490252568*JRatfR[4]-0.3535533905932737*JRatfR[3]+0.4743416490252568*JRatfR[2]+0.3535533905932737*(JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[8] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]-0.5809475019311124*JfL[15]-0.33541019662496846*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]-0.33541019662496846*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]-0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3162277660168379*JRatfR[11]+0.31622776601683794*(JRatfR[10]+JRatfR[9])+0.3162277660168379*JRatfR[8]-0.4743416490252568*(JRatfR[7]+JRatfR[6])+0.3535533905932737*JRatfR[5]-0.4743416490252568*JRatfR[4]+0.3535533905932737*JRatfR[3]-0.4743416490252568*JRatfR[2]+0.3535533905932737*(JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[9] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = -(0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]-0.4841229182759271*(JfL[21]+JfL[20])-0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]+0.4330127018922193*(JfL[8]+JfL[5])+0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = -(0.3952847075210473*JRatfR[11])-0.39528470752104733*(JRatfR[10]+JRatfR[9])-0.3952847075210473*JRatfR[8]+0.3535533905932737*(JRatfR[5]+JRatfR[3]+JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[10] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 

  alpha_quad = 0.0; 
  alpha_quad = alpha_quad * area_elem_quad/Jc_quad; 

  cfl = fmax(fabs(alpha_quad), fabs(cfl)); 
  JfL_quad = 0.3872983346207417*JfL[23]+0.22360679774997896*JfL[22]+0.38729833462074165*(JfL[21]+JfL[20])+0.22360679774997894*JfL[19]+0.22360679774997902*JfL[18]+0.3872983346207417*JfL[17]+0.22360679774997896*JfL[16]+0.5809475019311124*JfL[15]+0.33541019662496846*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.33541019662496846*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]+0.33541019662496846*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]+0.33541019662496846*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0]; 
  JfR_quad = 0.3162277660168379*JRatfR[11]+0.31622776601683794*(JRatfR[10]+JRatfR[9])+0.3162277660168379*JRatfR[8]+0.4743416490252568*(JRatfR[7]+JRatfR[6])+0.3535533905932737*JRatfR[5]+0.4743416490252568*JRatfR[4]+0.3535533905932737*JRatfR[3]+0.4743416490252568*JRatfR[2]+0.3535533905932737*(JRatfR[1]+JRatfR[0]); 
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0; 
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0; 
  flux_surf_nodal[11] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad; 


  return cfl*1.5*rdx2; 

} 
