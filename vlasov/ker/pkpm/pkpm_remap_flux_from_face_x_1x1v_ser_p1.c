#include <gkyl_pkpm_remap_flux_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_remap_flux_from_face_x_1x1v_ser_p1(const double *wpsi, const double *bck, const double *P,
  double* GKYL_RESTRICT flux) 
{ 
  // wpsi:  Weight of the transverse node times the transverse modes at the node. 
  // bck:   Factors of the change of frame from the face frame to the side (PKPM_FRAME_*). 
  // P:     Legendre projections in w of the face-frame flux of the eight blocks on the 
  //        velocity cell of the side at the node: P[m*8 + c]. 
  // flux:  Accumulated surface modes of the side of the face: flux[c*3 + k]. 
  double Z[24]; 
  pkpm_harm_frame_change(&P[0], bck, &Z[0]); 
  pkpm_harm_frame_change(&P[8], bck, &Z[8]); 
  pkpm_harm_frame_change(&P[16], bck, &Z[16]); 
  for (int c = 0; c < 8; ++c) { 
    double *out = &flux[c*3]; 
    out[0] += 1.4142135623730951*wpsi[0]*Z[0 + c]; 
    out[1] += 2.041241452319315*wpsi[0]*Z[8 + c]; 
    out[2] += -2.432521277052601*wpsi[0]*Z[16 + c]; 
  } 
} 
