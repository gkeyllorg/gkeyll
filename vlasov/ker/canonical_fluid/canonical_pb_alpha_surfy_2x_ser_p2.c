#include <gkyl_canonical_fluid_kernels.h> 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x_ser_p2(const double *w, const double *dxv, const double *phi,
   double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];

  double *alphaL = &alpha_surf[3];
  double *sgn_alpha_surfL = &sgn_alpha_surf[3];
  alphaL[0] = -(2.7386127875258306*phi[7]*rdx2)+2.1213203435596424*phi[3]*rdx2-1.224744871391589*phi[1]*rdx2; 
  alphaL[1] = 4.743416490252569*phi[6]*rdx2-2.7386127875258306*phi[4]*rdx2; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.7071067811865468*alphaL[0]-0.9486832980505135*alphaL[1] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.7071067811865468*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.9486832980505135*alphaL[1]+0.7071067811865468*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
