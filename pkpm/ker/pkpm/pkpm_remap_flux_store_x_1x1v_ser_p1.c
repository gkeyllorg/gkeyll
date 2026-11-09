#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_store_x_1x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux) 
{ 
  // wpsi:  Weight of the transverse node times the transverse modes at the node. 
  // proj:  Legendre projections in w of the eight blocks at the node: proj[m*8 + c]. 
  // flux:  Accumulated surface modes of one side of the face: flux[c*3 + k]. 
  for (int c = 0; c < 8; ++c) { 
    double *out = &flux[c*3]; 
    out[0] += 1.4142135623730951*wpsi[0]*proj[0 + c]; 
    out[1] += 2.041241452319315*wpsi[0]*proj[8 + c]; 
    out[2] += -2.432521277052601*wpsi[0]*proj[16 + c]; 
  } 
} 
