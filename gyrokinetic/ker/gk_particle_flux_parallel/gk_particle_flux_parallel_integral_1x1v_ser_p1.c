#include <gkyl_gk_particle_flux_parallel_kernels.h>

GKYL_CU_DH void gk_particle_flux_parallel_integral_1x1v_ser_p1(double factor, const double *flux, double *out)
{
  // factor: signed velocity-cell measure, including 2*pi/m for 2v.
  // flux: nodal numerical flux on the parallel configuration face.
  // out: velocity integral at each configuration-space surface node.
  out[0] = 0.1111111111111111*(5.0*flux[2]+8.0*flux[1]+5.0*flux[0])*factor; 
}
