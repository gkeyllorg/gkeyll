#include <gkyl_dg_fv_interpolate_kernels.h> 

void dg_fv_interpolate_modal_to_fv_2x(const double *dg_data, double *fv_out[8]) { 
  // dg_data: Input DG data on coarse mesh. 
  // fv_out: Finite volume data on 2x refined mesh. Row-major ordering. 
 
  fv_out[0][0] = 0.375*dg_data[3]-0.4330127018922193*dg_data[2]-0.4330127018922193*dg_data[1]+0.5*dg_data[0]; 
  fv_out[1][0] = -(0.375*dg_data[3])+0.4330127018922193*dg_data[2]-0.4330127018922193*dg_data[1]+0.5*dg_data[0]; 
  fv_out[2][0] = -(0.375*dg_data[3])-0.4330127018922193*dg_data[2]+0.4330127018922193*dg_data[1]+0.5*dg_data[0]; 
  fv_out[3][0] = 0.375*dg_data[3]+0.4330127018922193*dg_data[2]+0.4330127018922193*dg_data[1]+0.5*dg_data[0]; 
}
