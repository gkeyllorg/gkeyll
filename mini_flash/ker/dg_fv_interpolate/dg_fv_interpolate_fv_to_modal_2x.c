#include <gkyl_dg_fv_interpolate_kernels.h> 

void dg_fv_interpolate_fv_to_modal_2x(const double *fv_data[8], double *dg_out) { 
  // fv_data: List of finite volume data in adjacent cells. Row-major ordering. 
  // out: DG data on 2x coarsened mesh. 
 
  dg_out[0] = 0.5*fv_data[3][0]+0.5*fv_data[2][0]+0.5*fv_data[1][0]+0.5*fv_data[0][0]; 
  dg_out[1] = 0.5773502691896258*fv_data[3][0]+0.5773502691896258*fv_data[2][0]-0.5773502691896258*fv_data[1][0]-0.5773502691896258*fv_data[0][0]; 
  dg_out[2] = 0.5773502691896258*fv_data[3][0]-0.5773502691896258*fv_data[2][0]+0.5773502691896258*fv_data[1][0]-0.5773502691896258*fv_data[0][0]; 
  dg_out[3] = 0.6666666666666666*fv_data[3][0]-0.6666666666666666*fv_data[2][0]-0.6666666666666666*fv_data[1][0]+0.6666666666666666*fv_data[0][0]; 
}
