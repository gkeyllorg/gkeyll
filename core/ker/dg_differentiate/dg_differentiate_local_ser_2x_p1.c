#include <gkyl_dg_differentiate_kernels.h>

GKYL_CU_DH void 
dg_differentiate_local_ser_2x_p1_dir0_order1(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 1);
  out[0] = 1.7320508075688772*fin[1]*rdx2fac; 
  out[1] = 0.0; 
  out[2] = 1.7320508075688772*fin[3]*rdx2fac; 
  out[3] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_ser_2x_p1_dir0_order2(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 2);
  out[0] = 0.0; 
  out[1] = 0.0; 
  out[2] = 0.0; 
  out[3] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_ser_2x_p1_dir1_order1(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 1);
  out[0] = 1.7320508075688772*fin[2]*rdx2fac; 
  out[1] = 1.7320508075688772*fin[3]*rdx2fac; 
  out[2] = 0.0; 
  out[3] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_ser_2x_p1_dir1_order2(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 2);
  out[0] = 0.0; 
  out[1] = 0.0; 
  out[2] = 0.0; 
  out[3] = 0.0; 

}

