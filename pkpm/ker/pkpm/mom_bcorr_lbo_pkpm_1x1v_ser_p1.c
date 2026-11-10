#include <gkyl_mom_bcorr_lbo_pkpm_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_pkpm_1x1v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[6]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f^(vmax)_(vmin) and vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += (-(1.5811388300841895*fIn[4])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*mass; 
    out[1] += (-(1.5811388300841898*fIn[5])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[0]*(-(1.5811388300841895*fIn[4])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*mass; 
    out[3] += vBoundary[0]*(-(1.5811388300841898*fIn[5])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*mass; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += (1.5811388300841895*fIn[4]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*mass; 
    out[1] += (1.5811388300841898*fIn[5]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[1]*(1.5811388300841895*fIn[4]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*mass; 
    out[3] += vBoundary[1]*(1.5811388300841898*fIn[5]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*mass; 
  }
 
} 
 
