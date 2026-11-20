#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_store_y_3x1v_ser_p1(const double *wpsi, const double *proj, double* GKYL_RESTRICT flux) 
{ 
  // wpsi:  Weight of the transverse node times the transverse modes at the node. 
  // proj:  Legendre projections in w of the eight blocks at the node: proj[m*8 + c]. 
  // flux:  Accumulated surface modes of one side of the face: flux[c*12 + k]. 
  for (int c = 0; c < 8; ++c) { 
    double *out = &flux[c*12]; 
    out[0] += 2.8284271247461907*wpsi[0]*proj[0 + c]; 
    out[1] += 4.082482904638631*wpsi[1]*proj[0 + c]; 
    out[2] += 4.082482904638631*wpsi[2]*proj[0 + c]; 
    out[3] += 4.082482904638631*wpsi[0]*proj[8 + c]; 
    out[4] += 5.892556509887896*wpsi[3]*proj[0 + c]; 
    out[5] += 5.892556509887896*wpsi[1]*proj[8 + c]; 
    out[6] += 5.892556509887896*wpsi[2]*proj[8 + c]; 
    out[7] += -4.865042554105202*wpsi[0]*proj[16 + c]; 
    out[8] += 8.505172717997146*wpsi[3]*proj[8 + c]; 
    out[9] += -7.022084070579055*wpsi[1]*proj[16 + c]; 
    out[10] += -7.022084070579055*wpsi[2]*proj[16 + c]; 
    out[11] += -10.1355053210525*wpsi[3]*proj[16 + c]; 
  } 
} 
