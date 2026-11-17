#include <gkyl_iso_euler_kernels.h> 
#include <math.h> 
GKYL_CU_DH int fluid_vars_iso_positivity_3x_ser_p1(double eps, double* GKYL_RESTRICT fluid) 
{ 
  // eps:   Relative floor: the vertex density is kept above eps times the cell value. 
  // fluid: In/out fluid state [rho, rho ux, rho uy, rho uz]; the density slopes are scaled toward the cell average. 
  // Returns -1 (inadmissible cell average: nothing done), 1 (density slopes scaled), 0 when untouched. 
  double *rho = &fluid[0]; 
 
  // Cell averages. 
  const double rho_avg = 0.3535533905932737*rho[0]; 
  if (!(rho_avg > 0.0)) return -1; 
 
  // Vertex values. 
  double rho_v[8]; 
  rho_v[0] = (0.3535533905932737)*rho[0] + (-0.6123724356957944)*rho[1] + (-0.6123724356957944)*rho[2] + (-0.6123724356957944)*rho[3] + (1.060660171779821)*rho[4] + (1.060660171779821)*rho[5] + (1.060660171779821)*rho[6] + (-1.837117307087383)*rho[7]; 
  rho_v[1] = (0.3535533905932737)*rho[0] + (0.6123724356957944)*rho[1] + (-0.6123724356957944)*rho[2] + (-0.6123724356957944)*rho[3] + (-1.060660171779821)*rho[4] + (-1.060660171779821)*rho[5] + (1.060660171779821)*rho[6] + (1.837117307087383)*rho[7]; 
  rho_v[2] = (0.3535533905932737)*rho[0] + (-0.6123724356957944)*rho[1] + (0.6123724356957944)*rho[2] + (-0.6123724356957944)*rho[3] + (-1.060660171779821)*rho[4] + (1.060660171779821)*rho[5] + (-1.060660171779821)*rho[6] + (1.837117307087383)*rho[7]; 
  rho_v[3] = (0.3535533905932737)*rho[0] + (0.6123724356957944)*rho[1] + (0.6123724356957944)*rho[2] + (-0.6123724356957944)*rho[3] + (1.060660171779821)*rho[4] + (-1.060660171779821)*rho[5] + (-1.060660171779821)*rho[6] + (-1.837117307087383)*rho[7]; 
  rho_v[4] = (0.3535533905932737)*rho[0] + (-0.6123724356957944)*rho[1] + (-0.6123724356957944)*rho[2] + (0.6123724356957944)*rho[3] + (1.060660171779821)*rho[4] + (-1.060660171779821)*rho[5] + (-1.060660171779821)*rho[6] + (1.837117307087383)*rho[7]; 
  rho_v[5] = (0.3535533905932737)*rho[0] + (0.6123724356957944)*rho[1] + (-0.6123724356957944)*rho[2] + (0.6123724356957944)*rho[3] + (-1.060660171779821)*rho[4] + (1.060660171779821)*rho[5] + (-1.060660171779821)*rho[6] + (-1.837117307087383)*rho[7]; 
  rho_v[6] = (0.3535533905932737)*rho[0] + (-0.6123724356957944)*rho[1] + (0.6123724356957944)*rho[2] + (0.6123724356957944)*rho[3] + (-1.060660171779821)*rho[4] + (-1.060660171779821)*rho[5] + (1.060660171779821)*rho[6] + (-1.837117307087383)*rho[7]; 
  rho_v[7] = (0.3535533905932737)*rho[0] + (0.6123724356957944)*rho[1] + (0.6123724356957944)*rho[2] + (0.6123724356957944)*rho[3] + (1.060660171779821)*rho[4] + (1.060660171779821)*rho[5] + (1.060660171779821)*rho[6] + (1.837117307087383)*rho[7]; 
 
  int flag = 0; 
  // (1) Density: scale the density slopes so that rho_v >= eps rho_avg. 
  double theta = 1.0; 
  const double rho_min = eps*rho_avg; 
  for (int k = 0; k < 8; ++k) { 
    if (rho_v[k] < rho_min) theta = fmin(theta, (rho_avg - rho_min)/(rho_avg - rho_v[k])); 
  } 
  if (theta < 1.0) { 
    for (int i = 1; i < 8; ++i) rho[i] *= theta; 
    flag |= 1; 
  } 
  return flag; 
} 
