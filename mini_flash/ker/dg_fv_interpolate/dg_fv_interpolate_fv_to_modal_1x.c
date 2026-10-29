#include <gkyl_dg_fv_interpolate_kernels.h> 

void dg_fv_interpolate_fv_to_modal_1x(const double *fv_data[8], double *dg_out) { 
  // fv_data: List of finite volume data in adjacent cells. Row-major ordering. 
  // out: DG data on 2x coarsened mesh. 
 
  dg_out[0] = 0.7071067811865475*fv_data[1][0]+0.7071067811865475*fv_data[0][0]; 
  dg_out[1] = 0.8164965809277261*fv_data[1][0]-0.8164965809277261*fv_data[0][0]; 
}
