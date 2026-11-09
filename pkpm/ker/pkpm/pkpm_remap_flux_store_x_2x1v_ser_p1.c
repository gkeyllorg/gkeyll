#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_store_x_2x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux) 
{ 
  // wpsi:  Weight of the transverse node times the transverse modes at the node. 
  // proj:  Legendre projections in w of the eight blocks at the node: proj[m*8 + c]. 
  // flux:  Accumulated surface modes of one side of the face: flux[c*6 + k]. 
  for (int c = 0; c < 8; ++c) { 
    double *out = &flux[c*6]; 
    out[0] += 2.0*wpsi[0]*proj[0 + c]; 
    out[1] += 2.886751345948129*wpsi[1]*proj[0 + c]; 
    out[2] += 2.886751345948129*wpsi[0]*proj[8 + c]; 
    out[3] += 4.166666666666667*wpsi[1]*proj[8 + c]; 
    out[4] += -3.4401045807689083*wpsi[0]*proj[16 + c]; 
    out[5] += -4.9653632643684835*wpsi[1]*proj[16 + c]; 
  } 
} 
