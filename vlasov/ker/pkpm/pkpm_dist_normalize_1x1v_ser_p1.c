#include <gkyl_pkpm_kinetic_kernels.h> 
GKYL_CU_DH int pkpm_dist_normalize_1x1v_ser_p1(const double *nmoms, const double *f, double* GKYL_RESTRICT fbar) 
{ 
  // nmoms: Input moments [M0, M2, M0_G, M1, ...] of the distributions (mass 1). 
  // f:     Input distribution functions [F, G]. 
  // fbar:  Output [F/M0, G s/M0] with s = (3 - M2/M0)/(2 M0_G/M0) (nodal). 
  // Returns the number of nodes where M0 <= 0 or M0_G <= 0 (both factors are zero there). 
 
  const double *n0 = &nmoms[0]; 
  const double *n2 = &nmoms[2]; 
  const double *nG = &nmoms[4]; 
  double fac_F[2] = {0.0}; 
  double fac_G[2] = {0.0}; 
  int num_null = 0; 
  { 
    const double n0_q = 0.7071067811865475*n0[0]-0.7071067811865475*n0[1]; 
    const double n2_q = 0.7071067811865475*n2[0]-0.7071067811865475*n2[1]; 
    const double nG_q = 0.7071067811865475*nG[0]-0.7071067811865475*nG[1]; 
    if (n0_q > 0.0 && nG_q > 0.0) { 
      fac_F[0] = 1.0/n0_q; 
      fac_G[0] = (3.0 - n2_q/n0_q)/(2.0*nG_q/n0_q)/n0_q; 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = 0.7071067811865475*n0[1]+0.7071067811865475*n0[0]; 
    const double n2_q = 0.7071067811865475*n2[1]+0.7071067811865475*n2[0]; 
    const double nG_q = 0.7071067811865475*nG[1]+0.7071067811865475*nG[0]; 
    if (n0_q > 0.0 && nG_q > 0.0) { 
      fac_F[1] = 1.0/n0_q; 
      fac_G[1] = (3.0 - n2_q/n0_q)/(2.0*nG_q/n0_q)/n0_q; 
    } else { 
      num_null += 1; 
    } 
  } 
 
  const double *F_0 = &f[0]; 
  double *fbar_F_0 = &fbar[0]; 
  fbar_F_0[0] = 0.5*F_0[1]*fac_F[1]+0.5*F_0[0]*fac_F[1]-0.5*fac_F[0]*F_0[1]+0.5*F_0[0]*fac_F[0]; 
  fbar_F_0[1] = 0.5*F_0[1]*fac_F[1]+0.5*F_0[0]*fac_F[1]+0.5*fac_F[0]*F_0[1]-0.5*F_0[0]*fac_F[0]; 
  fbar_F_0[2] = 0.5*fac_F[1]*F_0[3]-0.5*fac_F[0]*F_0[3]+0.5*fac_F[1]*F_0[2]+0.5*fac_F[0]*F_0[2]; 
  fbar_F_0[3] = 0.5*fac_F[1]*F_0[3]+0.5*fac_F[0]*F_0[3]+0.5*fac_F[1]*F_0[2]-0.5*fac_F[0]*F_0[2]; 
  fbar_F_0[4] = 0.5000000000000001*fac_F[1]*F_0[5]-0.5000000000000001*fac_F[0]*F_0[5]+0.5*fac_F[1]*F_0[4]+0.5*fac_F[0]*F_0[4]; 
  fbar_F_0[5] = 0.5*fac_F[1]*F_0[5]+0.5*fac_F[0]*F_0[5]+0.5000000000000001*fac_F[1]*F_0[4]-0.5000000000000001*fac_F[0]*F_0[4]; 
  const double *G_1 = &f[6]; 
  double *fbar_G_1 = &fbar[6]; 
  fbar_G_1[0] = 0.5*G_1[1]*fac_G[1]+0.5*G_1[0]*fac_G[1]-0.5*fac_G[0]*G_1[1]+0.5*G_1[0]*fac_G[0]; 
  fbar_G_1[1] = 0.5*G_1[1]*fac_G[1]+0.5*G_1[0]*fac_G[1]+0.5*fac_G[0]*G_1[1]-0.5*G_1[0]*fac_G[0]; 
  fbar_G_1[2] = 0.5*fac_G[1]*G_1[3]-0.5*fac_G[0]*G_1[3]+0.5*fac_G[1]*G_1[2]+0.5*fac_G[0]*G_1[2]; 
  fbar_G_1[3] = 0.5*fac_G[1]*G_1[3]+0.5*fac_G[0]*G_1[3]+0.5*fac_G[1]*G_1[2]-0.5*fac_G[0]*G_1[2]; 
  fbar_G_1[4] = 0.5000000000000001*fac_G[1]*G_1[5]-0.5000000000000001*fac_G[0]*G_1[5]+0.5*fac_G[1]*G_1[4]+0.5*fac_G[0]*G_1[4]; 
  fbar_G_1[5] = 0.5*fac_G[1]*G_1[5]+0.5*fac_G[0]*G_1[5]+0.5000000000000001*fac_G[1]*G_1[4]-0.5000000000000001*fac_G[0]*G_1[4]; 
 
  return num_null; 
} 
