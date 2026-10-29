#include <gkyl_dg_fv_interpolate_kernels.h> 

void dg_fv_interpolate_modal_to_fv_1x(const double *dg_data, double *fv_out[8]) { 
  // dg_data: Input DG data on coarse mesh. 
  // fv_out: Finite volume data on 2x refined mesh. Row-major ordering. 
 
  fv_out[0][0] = 0.7071067811865475*dg_data[0]-0.6123724356957944*dg_data[1]; 
  fv_out[1][0] = 0.6123724356957944*dg_data[1]+0.7071067811865475*dg_data[0]; 
}
