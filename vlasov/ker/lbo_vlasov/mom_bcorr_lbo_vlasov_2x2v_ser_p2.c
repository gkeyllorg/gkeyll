#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[48]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += (-(2.23606797749979*fIn[13])+1.7320508075688772*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (-(2.2360679774997902*fIn[23])+1.7320508075688772*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (-(2.2360679774997902*fIn[24])+1.7320508075688772*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (-(2.23606797749979*fIn[34])+1.7320508075688772*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += (1.7320508075688774*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += (1.7320508075688774*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += (1.7320508075688774*fIn[32]-1.0*fIn[19])*dS; 
    out[7] += (1.7320508075688774*fIn[33]-1.0*fIn[20])*dS; 
    out[16] += vBoundary[0]*(-(2.23606797749979*fIn[13])+1.7320508075688772*fIn[3]-1.0*fIn[0])*dS; 
    out[17] += vBoundary[0]*(-(2.2360679774997902*fIn[23])+1.7320508075688772*fIn[6]-1.0*fIn[1])*dS; 
    out[18] += vBoundary[0]*(-(2.2360679774997902*fIn[24])+1.7320508075688772*fIn[7]-1.0*fIn[2])*dS; 
    out[19] += vBoundary[0]*(-(2.23606797749979*fIn[34])+1.7320508075688772*fIn[15]-1.0*fIn[5])*dS; 
    out[20] += vBoundary[0]*(1.7320508075688774*fIn[21]-1.0*fIn[11])*dS; 
    out[21] += vBoundary[0]*(1.7320508075688774*fIn[22]-1.0*fIn[12])*dS; 
    out[22] += vBoundary[0]*(1.7320508075688774*fIn[32]-1.0*fIn[19])*dS; 
    out[23] += vBoundary[0]*(1.7320508075688774*fIn[33]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += (2.23606797749979*fIn[13]+1.7320508075688772*fIn[3]+fIn[0])*dS; 
    out[1] += (2.2360679774997902*fIn[23]+1.7320508075688772*fIn[6]+fIn[1])*dS; 
    out[2] += (2.2360679774997902*fIn[24]+1.7320508075688772*fIn[7]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[34]+1.7320508075688772*fIn[15]+fIn[5])*dS; 
    out[4] += (1.7320508075688774*fIn[21]+fIn[11])*dS; 
    out[5] += (1.7320508075688774*fIn[22]+fIn[12])*dS; 
    out[6] += (1.7320508075688774*fIn[32]+fIn[19])*dS; 
    out[7] += (1.7320508075688774*fIn[33]+fIn[20])*dS; 
    out[16] += vBoundary[2]*(2.23606797749979*fIn[13]+1.7320508075688772*fIn[3]+fIn[0])*dS; 
    out[17] += vBoundary[2]*(2.2360679774997902*fIn[23]+1.7320508075688772*fIn[6]+fIn[1])*dS; 
    out[18] += vBoundary[2]*(2.2360679774997902*fIn[24]+1.7320508075688772*fIn[7]+fIn[2])*dS; 
    out[19] += vBoundary[2]*(2.23606797749979*fIn[34]+1.7320508075688772*fIn[15]+fIn[5])*dS; 
    out[20] += vBoundary[2]*(1.7320508075688774*fIn[21]+fIn[11])*dS; 
    out[21] += vBoundary[2]*(1.7320508075688774*fIn[22]+fIn[12])*dS; 
    out[22] += vBoundary[2]*(1.7320508075688774*fIn[32]+fIn[19])*dS; 
    out[23] += vBoundary[2]*(1.7320508075688774*fIn[33]+fIn[20])*dS; 
  }
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[2]; 
    out[8] += (-(2.23606797749979*fIn[14])+1.7320508075688772*fIn[4]-1.0*fIn[0])*dS; 
    out[9] += (-(2.2360679774997902*fIn[28])+1.7320508075688772*fIn[8]-1.0*fIn[1])*dS; 
    out[10] += (-(2.2360679774997902*fIn[29])+1.7320508075688772*fIn[9]-1.0*fIn[2])*dS; 
    out[11] += (-(2.23606797749979*fIn[41])+1.7320508075688772*fIn[16]-1.0*fIn[5])*dS; 
    out[12] += (1.7320508075688774*fIn[25]-1.0*fIn[11])*dS; 
    out[13] += (1.7320508075688774*fIn[26]-1.0*fIn[12])*dS; 
    out[14] += (1.7320508075688774*fIn[35]-1.0*fIn[19])*dS; 
    out[15] += (1.7320508075688774*fIn[36]-1.0*fIn[20])*dS; 
    out[16] += vBoundary[1]*(-(2.23606797749979*fIn[14])+1.7320508075688772*fIn[4]-1.0*fIn[0])*dS; 
    out[17] += vBoundary[1]*(-(2.2360679774997902*fIn[28])+1.7320508075688772*fIn[8]-1.0*fIn[1])*dS; 
    out[18] += vBoundary[1]*(-(2.2360679774997902*fIn[29])+1.7320508075688772*fIn[9]-1.0*fIn[2])*dS; 
    out[19] += vBoundary[1]*(-(2.23606797749979*fIn[41])+1.7320508075688772*fIn[16]-1.0*fIn[5])*dS; 
    out[20] += vBoundary[1]*(1.7320508075688774*fIn[25]-1.0*fIn[11])*dS; 
    out[21] += vBoundary[1]*(1.7320508075688774*fIn[26]-1.0*fIn[12])*dS; 
    out[22] += vBoundary[1]*(1.7320508075688774*fIn[35]-1.0*fIn[19])*dS; 
    out[23] += vBoundary[1]*(1.7320508075688774*fIn[36]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[2]; 
    out[8] += (2.23606797749979*fIn[14]+1.7320508075688772*fIn[4]+fIn[0])*dS; 
    out[9] += (2.2360679774997902*fIn[28]+1.7320508075688772*fIn[8]+fIn[1])*dS; 
    out[10] += (2.2360679774997902*fIn[29]+1.7320508075688772*fIn[9]+fIn[2])*dS; 
    out[11] += (2.23606797749979*fIn[41]+1.7320508075688772*fIn[16]+fIn[5])*dS; 
    out[12] += (1.7320508075688774*fIn[25]+fIn[11])*dS; 
    out[13] += (1.7320508075688774*fIn[26]+fIn[12])*dS; 
    out[14] += (1.7320508075688774*fIn[35]+fIn[19])*dS; 
    out[15] += (1.7320508075688774*fIn[36]+fIn[20])*dS; 
    out[16] += vBoundary[3]*(2.23606797749979*fIn[14]+1.7320508075688772*fIn[4]+fIn[0])*dS; 
    out[17] += vBoundary[3]*(2.2360679774997902*fIn[28]+1.7320508075688772*fIn[8]+fIn[1])*dS; 
    out[18] += vBoundary[3]*(2.2360679774997902*fIn[29]+1.7320508075688772*fIn[9]+fIn[2])*dS; 
    out[19] += vBoundary[3]*(2.23606797749979*fIn[41]+1.7320508075688772*fIn[16]+fIn[5])*dS; 
    out[20] += vBoundary[3]*(1.7320508075688774*fIn[25]+fIn[11])*dS; 
    out[21] += vBoundary[3]*(1.7320508075688774*fIn[26]+fIn[12])*dS; 
    out[22] += vBoundary[3]*(1.7320508075688774*fIn[35]+fIn[19])*dS; 
    out[23] += vBoundary[3]*(1.7320508075688774*fIn[36]+fIn[20])*dS; 
  }
 
} 
 
