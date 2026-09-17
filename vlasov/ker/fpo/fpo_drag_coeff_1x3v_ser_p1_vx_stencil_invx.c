#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, 
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
  double dv1 = 2.0/dxv[1]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *dHdv_surf_C = &fpo_dhdv_surf[0]; 
  
  double *out = &drag_coeff[0]; 
  double *out_surf = &drag_coeff_surf[0]; 
  
  out[0] = (-(0.5773502691896258*(H_R[2]+H_L[2]))+1.1547005383792517*H_C[2]+0.5*H_R[0]-0.5*H_L[0])*dv1*gamma_avg; 
  out[1] = (-(0.5773502691896258*(H_R[5]+H_L[5]))+1.1547005383792517*H_C[5]+0.5*H_R[1]-0.5*H_L[1])*dv1*gamma_avg; 
  out[2] = (-(1.0*H_R[2])+H_L[2]+0.8660254037844386*(H_R[0]+H_L[0])-1.7320508075688772*H_C[0])*dv1*gamma_avg; 
  out[3] = (-(0.5773502691896258*(H_R[7]+H_L[7]))+1.1547005383792517*H_C[7]+0.5*H_R[3]-0.5*H_L[3])*dv1*gamma_avg; 
  out[4] = (-(0.5773502691896258*(H_R[9]+H_L[9]))+1.1547005383792517*H_C[9]+0.5*H_R[4]-0.5*H_L[4])*dv1*gamma_avg; 
  out[5] = (-(1.0*H_R[5])+H_L[5]+0.8660254037844386*(H_R[1]+H_L[1])-1.7320508075688772*H_C[1])*dv1*gamma_avg; 
  out[6] = (-(0.5773502691896258*(H_R[11]+H_L[11]))+1.1547005383792517*H_C[11]+0.5*H_R[6]-0.5*H_L[6])*dv1*gamma_avg; 
  out[7] = (-(1.0*H_R[7])+H_L[7]+0.8660254037844386*(H_R[3]+H_L[3])-1.7320508075688772*H_C[3])*dv1*gamma_avg; 
  out[8] = (-(0.5773502691896258*(H_R[12]+H_L[12]))+1.1547005383792517*H_C[12]+0.5*H_R[8]-0.5*H_L[8])*dv1*gamma_avg; 
  out[9] = (-(1.0*H_R[9])+H_L[9]+0.8660254037844386*(H_R[4]+H_L[4])-1.7320508075688772*H_C[4])*dv1*gamma_avg; 
  out[10] = (-(0.5773502691896258*(H_R[14]+H_L[14]))+1.1547005383792517*H_C[14]+0.5*H_R[10]-0.5*H_L[10])*dv1*gamma_avg; 
  out[11] = (-(1.0*H_R[11])+H_L[11]+0.8660254037844386*(H_R[6]+H_L[6])-1.7320508075688772*H_C[6])*dv1*gamma_avg; 
  out[12] = (-(1.0*H_R[12])+H_L[12]+0.8660254037844386*(H_R[8]+H_L[8])-1.7320508075688772*H_C[8])*dv1*gamma_avg; 
  out[13] = (-(0.5773502691896258*(H_R[15]+H_L[15]))+1.1547005383792517*H_C[15]+0.5*H_R[13]-0.5*H_L[13])*dv1*gamma_avg; 
  out[14] = (-(1.0*H_R[14])+H_L[14]+0.8660254037844386*(H_R[10]+H_L[10])-1.7320508075688772*H_C[10])*dv1*gamma_avg; 
  out[15] = (-(1.0*H_R[15])+H_L[15]+0.8660254037844386*(H_R[13]+H_L[13])-1.7320508075688772*H_C[13])*dv1*gamma_avg; 

  out_surf[0] = -(0.1767766952966368*(8.660254037844386*H_L[2]+8.660254037844386*H_C[2]+9.0*H_L[0]-9.0*H_C[0])*dv1*gamma_avg); 
  out_surf[1] = -(0.1767766952966368*(8.660254037844386*H_L[5]+8.660254037844386*H_C[5]+9.0*H_L[1]-9.0*H_C[1])*dv1*gamma_avg); 
  out_surf[2] = -(0.1767766952966368*(8.660254037844386*H_L[7]+8.660254037844386*H_C[7]+9.0*H_L[3]-9.0*H_C[3])*dv1*gamma_avg); 
  out_surf[3] = -(0.1767766952966368*(8.660254037844386*H_L[9]+8.660254037844386*H_C[9]+9.0*H_L[4]-9.0*H_C[4])*dv1*gamma_avg); 
  out_surf[4] = -(0.1767766952966368*(8.660254037844386*H_L[11]+8.660254037844386*H_C[11]+9.0*H_L[6]-9.0*H_C[6])*dv1*gamma_avg); 
  out_surf[5] = -(0.1767766952966368*(8.660254037844386*H_L[12]+8.660254037844386*H_C[12]+9.0*H_L[8]-9.0*H_C[8])*dv1*gamma_avg); 
  out_surf[6] = -(0.1767766952966368*(8.660254037844386*H_L[14]+8.660254037844386*H_C[14]+9.0*H_L[10]-9.0*H_C[10])*dv1*gamma_avg); 
  out_surf[7] = -(0.1767766952966368*(8.660254037844386*H_L[15]+8.660254037844386*H_C[15]+9.0*H_L[13]-9.0*H_C[13])*dv1*gamma_avg); 
}

