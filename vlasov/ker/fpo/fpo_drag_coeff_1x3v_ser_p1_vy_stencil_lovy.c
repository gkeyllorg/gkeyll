#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p1_lovy(const double *dxv, const double *gamma, 
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
  double dv1 = 2.0/dxv[2]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *dHdv_surf_C = &fpo_dhdv_surf[8]; 
  
  double *out = &drag_coeff[16]; 
  double *out_surf = &drag_coeff_surf[8]; 
  
  out[0] = ((-(0.18042195912175807*H_R[3])+2.7063293868263707*H_C[3]+0.1875*H_R[0]-0.1875*H_C[0])*dv1+0.2357022603955158*dHdv_surf_C[0])*gamma_avg; 
  out[1] = ((-(0.18042195912175807*H_R[6])+2.7063293868263707*H_C[6]+0.1875*H_R[1]-0.1875*H_C[1])*dv1+0.2357022603955158*dHdv_surf_C[1])*gamma_avg; 
  out[2] = ((-(0.18042195912175807*H_R[7])+2.7063293868263707*H_C[7]+0.1875*H_R[2]-0.1875*H_C[2])*dv1+0.2357022603955158*dHdv_surf_C[2])*gamma_avg; 
  out[3] = ((-(1.6875*H_R[3])-2.6875*H_C[3]+1.4072912811497127*H_R[0]-1.4072912811497127*H_C[0])*dv1-0.408248290463863*dHdv_surf_C[0])*gamma_avg; 
  out[4] = ((-(0.18042195912175807*H_R[10])+2.7063293868263707*H_C[10]+0.1875*H_R[4]-0.1875*H_C[4])*dv1+0.2357022603955158*dHdv_surf_C[3])*gamma_avg; 
  out[5] = ((-(0.18042195912175807*H_R[11])+2.7063293868263707*H_C[11]+0.1875*H_R[5]-0.1875*H_C[5])*dv1+0.2357022603955158*dHdv_surf_C[4])*gamma_avg; 
  out[6] = ((-(1.6875*H_R[6])-2.6875*H_C[6]+1.4072912811497127*H_R[1]-1.4072912811497127*H_C[1])*dv1-0.408248290463863*dHdv_surf_C[1])*gamma_avg; 
  out[7] = ((-(1.6875*H_R[7])-2.6875*H_C[7]+1.4072912811497127*H_R[2]-1.4072912811497127*H_C[2])*dv1-0.408248290463863*dHdv_surf_C[2])*gamma_avg; 
  out[8] = ((-(0.18042195912175807*H_R[13])+2.7063293868263707*H_C[13]+0.1875*H_R[8]-0.1875*H_C[8])*dv1+0.2357022603955158*dHdv_surf_C[5])*gamma_avg; 
  out[9] = ((-(0.18042195912175807*H_R[14])+2.7063293868263707*H_C[14]+0.1875*H_R[9]-0.1875*H_C[9])*dv1+0.2357022603955158*dHdv_surf_C[6])*gamma_avg; 
  out[10] = ((-(1.6875*H_R[10])-2.6875*H_C[10]+1.4072912811497127*H_R[4]-1.4072912811497127*H_C[4])*dv1-0.408248290463863*dHdv_surf_C[3])*gamma_avg; 
  out[11] = ((-(1.6875*H_R[11])-2.6875*H_C[11]+1.4072912811497127*H_R[5]-1.4072912811497127*H_C[5])*dv1-0.408248290463863*dHdv_surf_C[4])*gamma_avg; 
  out[12] = ((-(0.18042195912175807*H_R[15])+2.7063293868263707*H_C[15]+0.1875*H_R[12]-0.1875*H_C[12])*dv1+0.2357022603955158*dHdv_surf_C[7])*gamma_avg; 
  out[13] = ((-(1.6875*H_R[13])-2.6875*H_C[13]+1.4072912811497127*H_R[8]-1.4072912811497127*H_C[8])*dv1-0.408248290463863*dHdv_surf_C[5])*gamma_avg; 
  out[14] = ((-(1.6875*H_R[14])-2.6875*H_C[14]+1.4072912811497127*H_R[9]-1.4072912811497127*H_C[9])*dv1-0.408248290463863*dHdv_surf_C[6])*gamma_avg; 
  out[15] = ((-(1.6875*H_R[15])-2.6875*H_C[15]+1.4072912811497127*H_R[12]-1.4072912811497127*H_C[12])*dv1-0.408248290463863*dHdv_surf_C[7])*gamma_avg; 

  out_surf[0] = 0.0; 
  out_surf[1] = 0.0; 
  out_surf[2] = 0.0; 
  out_surf[3] = 0.0; 
  out_surf[4] = 0.0; 
  out_surf[5] = 0.0; 
  out_surf[6] = 0.0; 
  out_surf[7] = 0.0; 
}

