#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_1x2v_tensor_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:   cell length in each direction. 
  // fIn[18]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[2]; 
    out[0] += (-(2.23606797749979*fIn[7])+1.7320508075688772*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (-(2.2360679774997902*fIn[10])+1.7320508075688772*fIn[4]-1.0*fIn[1])*dS; 
    out[4] += vBoundary[0]*(-(2.23606797749979*fIn[7])+1.7320508075688772*fIn[2]-1.0*fIn[0])*dS; 
    out[5] += vBoundary[0]*(-(2.2360679774997902*fIn[10])+1.7320508075688772*fIn[4]-1.0*fIn[1])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[2]; 
    out[0] += (2.23606797749979*fIn[7]+1.7320508075688772*fIn[2]+fIn[0])*dS; 
    out[1] += (2.2360679774997902*fIn[10]+1.7320508075688772*fIn[4]+fIn[1])*dS; 
    out[4] += vBoundary[2]*(2.23606797749979*fIn[7]+1.7320508075688772*fIn[2]+fIn[0])*dS; 
    out[5] += vBoundary[2]*(2.2360679774997902*fIn[10]+1.7320508075688772*fIn[4]+fIn[1])*dS; 
  }
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[1]; 
    out[2] += (-(2.23606797749979*fIn[8])+1.7320508075688772*fIn[3]-1.0*fIn[0])*dS; 
    out[3] += (-(2.2360679774997902*fIn[12])+1.7320508075688772*fIn[5]-1.0*fIn[1])*dS; 
    out[4] += vBoundary[1]*(-(2.23606797749979*fIn[8])+1.7320508075688772*fIn[3]-1.0*fIn[0])*dS; 
    out[5] += vBoundary[1]*(-(2.2360679774997902*fIn[12])+1.7320508075688772*fIn[5]-1.0*fIn[1])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[1]; 
    out[2] += (2.23606797749979*fIn[8]+1.7320508075688772*fIn[3]+fIn[0])*dS; 
    out[3] += (2.2360679774997902*fIn[12]+1.7320508075688772*fIn[5]+fIn[1])*dS; 
    out[4] += vBoundary[3]*(2.23606797749979*fIn[8]+1.7320508075688772*fIn[3]+fIn[0])*dS; 
    out[5] += vBoundary[3]*(2.2360679774997902*fIn[12]+1.7320508075688772*fIn[5]+fIn[1])*dS; 
  }
 
} 
 
