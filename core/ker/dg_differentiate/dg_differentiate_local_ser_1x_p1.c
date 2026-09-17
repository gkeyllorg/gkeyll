#include <gkyl_dg_differentiate_kernels.h>

GKYL_CU_DH void 
dg_differentiate_local_ser_1x_p1_dir0_order1(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 1);
  out[0] = 1.7320508075688772*fin[1]*rdx2fac; 
  out[1] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_ser_1x_p1_dir0_order2(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 2);
  out[0] = 0.0; 
  out[1] = 0.0; 

}

