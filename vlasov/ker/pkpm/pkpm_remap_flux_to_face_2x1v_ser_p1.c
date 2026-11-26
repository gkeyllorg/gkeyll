#include <gkyl_pkpm_remap_flux_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_remap_flux_to_face_2x1v_ser_p1(const double *psi, const double *fwd, const double *G,
  double* GKYL_RESTRICT X, double* GKYL_RESTRICT Y) 
{ 
  // psi:  Transverse modes at the node. 
  // fwd:  Factors of the change of frame from the side to the face frame (PKPM_FRAME_*). 
  // G:    Modal traces of the eight blocks of one velocity cell of the side: 
  //       G[c*6 + a*3 + m], a transverse mode, m Legendre mode in w. 
  // X:    Output blocks at the node in the frame of the side, Legendre coefficients in w: 
  //       X[m*8 + c] of [1, eta, 3 eta^2 - 1]. 
  // Y:    Output blocks at the node in the face frame, Y[m*8 + c]. 
  for (int c = 0; c < 8; ++c) { 
    const double *g = &G[c*6]; 
    X[c] = 0.7071067811865475*(psi[0]*g[0]+psi[1]*g[3]); 
    X[8+c] = 1.224744871391589*(psi[0]*g[1]+psi[1]*g[4]); 
    X[16+c] = 0.7905694150420947*(psi[0]*g[2]+psi[1]*g[5]); 
  } 
  pkpm_harm_frame_change(&X[0], fwd, &Y[0]); 
  pkpm_harm_frame_change(&X[8], fwd, &Y[8]); 
  pkpm_harm_frame_change(&X[16], fwd, &Y[16]); 
} 
