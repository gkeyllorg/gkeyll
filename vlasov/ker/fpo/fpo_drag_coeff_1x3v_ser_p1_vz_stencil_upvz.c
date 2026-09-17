#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, 
    const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, 
    double *drag_coeff_surf) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *dHdv_surf_C = &fpo_dhdv_surf[16]; 
  
  double *out = &drag_coeff[32]; 
  double *out_surf = &drag_coeff_surf[16]; 
  
  out[0] = ((-(0.18042195912175807*H_L[4])+2.7063293868263707*H_C[4]-0.1875*H_L[0]+0.1875*H_C[0])*dv1+0.2357022603955158*dHdv_surf_C[0])*gamma_avg; 
  out[1] = ((-(0.18042195912175807*H_L[8])+2.7063293868263707*H_C[8]-0.1875*H_L[1]+0.1875*H_C[1])*dv1+0.2357022603955158*dHdv_surf_C[1])*gamma_avg; 
  out[2] = ((-(0.18042195912175807*H_L[9])+2.7063293868263707*H_C[9]-0.1875*H_L[2]+0.1875*H_C[2])*dv1+0.2357022603955158*dHdv_surf_C[2])*gamma_avg; 
  out[3] = ((-(0.18042195912175807*H_L[10])+2.7063293868263707*H_C[10]-0.1875*H_L[3]+0.1875*H_C[3])*dv1+0.2357022603955158*dHdv_surf_C[3])*gamma_avg; 
  out[4] = ((1.6875*H_L[4]+2.6875*H_C[4]+1.4072912811497127*H_L[0]-1.4072912811497127*H_C[0])*dv1+0.408248290463863*dHdv_surf_C[0])*gamma_avg; 
  out[5] = ((-(0.18042195912175807*H_L[12])+2.7063293868263707*H_C[12]-0.1875*H_L[5]+0.1875*H_C[5])*dv1+0.2357022603955158*dHdv_surf_C[4])*gamma_avg; 
  out[6] = ((-(0.18042195912175807*H_L[13])+2.7063293868263707*H_C[13]-0.1875*H_L[6]+0.1875*H_C[6])*dv1+0.2357022603955158*dHdv_surf_C[5])*gamma_avg; 
  out[7] = ((-(0.18042195912175807*H_L[14])+2.7063293868263707*H_C[14]-0.1875*H_L[7]+0.1875*H_C[7])*dv1+0.2357022603955158*dHdv_surf_C[6])*gamma_avg; 
  out[8] = ((1.6875*H_L[8]+2.6875*H_C[8]+1.4072912811497127*H_L[1]-1.4072912811497127*H_C[1])*dv1+0.408248290463863*dHdv_surf_C[1])*gamma_avg; 
  out[9] = ((1.6875*H_L[9]+2.6875*H_C[9]+1.4072912811497127*H_L[2]-1.4072912811497127*H_C[2])*dv1+0.408248290463863*dHdv_surf_C[2])*gamma_avg; 
  out[10] = ((1.6875*H_L[10]+2.6875*H_C[10]+1.4072912811497127*H_L[3]-1.4072912811497127*H_C[3])*dv1+0.408248290463863*dHdv_surf_C[3])*gamma_avg; 
  out[11] = ((-(0.18042195912175807*H_L[15])+2.7063293868263707*H_C[15]-0.1875*H_L[11]+0.1875*H_C[11])*dv1+0.2357022603955158*dHdv_surf_C[7])*gamma_avg; 
  out[12] = ((1.6875*H_L[12]+2.6875*H_C[12]+1.4072912811497127*H_L[5]-1.4072912811497127*H_C[5])*dv1+0.408248290463863*dHdv_surf_C[4])*gamma_avg; 
  out[13] = ((1.6875*H_L[13]+2.6875*H_C[13]+1.4072912811497127*H_L[6]-1.4072912811497127*H_C[6])*dv1+0.408248290463863*dHdv_surf_C[5])*gamma_avg; 
  out[14] = ((1.6875*H_L[14]+2.6875*H_C[14]+1.4072912811497127*H_L[7]-1.4072912811497127*H_C[7])*dv1+0.408248290463863*dHdv_surf_C[6])*gamma_avg; 
  out[15] = ((1.6875*H_L[15]+2.6875*H_C[15]+1.4072912811497127*H_L[11]-1.4072912811497127*H_C[11])*dv1+0.408248290463863*dHdv_surf_C[7])*gamma_avg; 

  out_surf[0] = -(0.1767766952966368*(8.660254037844386*H_L[4]+8.660254037844386*H_C[4]+9.0*H_L[0]-9.0*H_C[0])*dv1*gamma_avg); 
  out_surf[1] = -(0.1767766952966368*(8.660254037844386*H_L[8]+8.660254037844386*H_C[8]+9.0*H_L[1]-9.0*H_C[1])*dv1*gamma_avg); 
  out_surf[2] = -(0.1767766952966368*(8.660254037844386*H_L[9]+8.660254037844386*H_C[9]+9.0*H_L[2]-9.0*H_C[2])*dv1*gamma_avg); 
  out_surf[3] = -(0.1767766952966368*(8.660254037844386*H_L[10]+8.660254037844386*H_C[10]+9.0*H_L[3]-9.0*H_C[3])*dv1*gamma_avg); 
  out_surf[4] = -(0.1767766952966368*(8.660254037844386*H_L[12]+8.660254037844386*H_C[12]+9.0*H_L[5]-9.0*H_C[5])*dv1*gamma_avg); 
  out_surf[5] = -(0.1767766952966368*(8.660254037844386*H_L[13]+8.660254037844386*H_C[13]+9.0*H_L[6]-9.0*H_C[6])*dv1*gamma_avg); 
  out_surf[6] = -(0.1767766952966368*(8.660254037844386*H_L[14]+8.660254037844386*H_C[14]+9.0*H_L[7]-9.0*H_C[7])*dv1*gamma_avg); 
  out_surf[7] = -(0.1767766952966368*(8.660254037844386*H_L[15]+8.660254037844386*H_C[15]+9.0*H_L[11]-9.0*H_C[11])*dv1*gamma_avg); 
}

