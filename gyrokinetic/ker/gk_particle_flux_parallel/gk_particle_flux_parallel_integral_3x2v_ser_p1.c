#include <gkyl_gk_particle_flux_parallel_kernels.h>

GKYL_CU_DH void gk_particle_flux_parallel_integral_3x2v_ser_p1(double factor, const double *flux, double *out)
{
  // factor: signed velocity-cell measure, including 2*pi/m for 2v.
  // flux: nodal numerical flux on the parallel configuration face.
  // out: velocity integral at each configuration-space surface node.
  out[0] = 0.1111111111111111*(5.0*flux[5]+8.0*flux[4]+5.0*(flux[3]+flux[2])+8.0*flux[1]+5.0*flux[0])*factor; 
  out[1] = 0.1111111111111111*(5.0*flux[11]+8.0*flux[10]+5.0*(flux[9]+flux[8])+8.0*flux[7]+5.0*flux[6])*factor; 
  out[2] = 0.1111111111111111*(5.0*flux[17]+8.0*flux[16]+5.0*(flux[15]+flux[14])+8.0*flux[13]+5.0*flux[12])*factor; 
  out[3] = 0.1111111111111111*(5.0*flux[23]+8.0*flux[22]+5.0*(flux[21]+flux[20])+8.0*flux[19]+5.0*flux[18])*factor; 
}
