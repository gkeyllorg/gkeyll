#include <gkyl_dg_differentiate_kernels.h>

GKYL_CU_DH void 
dg_differentiate_local_tensor_2x_p2_dir0_order1(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 1);
  out[0] = 1.7320508075688772*fin[1]*rdx2fac; 
  out[1] = 3.872983346207417*fin[4]*rdx2fac; 
  out[2] = 1.7320508075688772*fin[3]*rdx2fac; 
  out[3] = 3.872983346207417*fin[6]*rdx2fac; 
  out[4] = 0.0; 
  out[5] = 1.7320508075688774*fin[7]*rdx2fac; 
  out[6] = 0.0; 
  out[7] = 3.872983346207417*fin[8]*rdx2fac; 
  out[8] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_tensor_2x_p2_dir0_order2(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 2);
  out[0] = 6.708203932499369*fin[4]*rdx2fac; 
  out[1] = 0.0; 
  out[2] = 6.7082039324993685*fin[6]*rdx2fac; 
  out[3] = 0.0; 
  out[4] = 0.0; 
  out[5] = 6.708203932499369*fin[8]*rdx2fac; 
  out[6] = 0.0; 
  out[7] = 0.0; 
  out[8] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_tensor_2x_p2_dir1_order1(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 1);
  out[0] = 1.7320508075688772*fin[2]*rdx2fac; 
  out[1] = 1.7320508075688772*fin[3]*rdx2fac; 
  out[2] = 3.872983346207417*fin[5]*rdx2fac; 
  out[3] = 3.872983346207417*fin[7]*rdx2fac; 
  out[4] = 1.7320508075688774*fin[6]*rdx2fac; 
  out[5] = 0.0; 
  out[6] = 3.872983346207417*fin[8]*rdx2fac; 
  out[7] = 0.0; 
  out[8] = 0.0; 

}

GKYL_CU_DH void 
dg_differentiate_local_tensor_2x_p2_dir1_order2(double dx, const double *fin, double *out) 
{ 
  // dx: cell size in each direction. 
  // fin: Input scalar DG field. 
  // out: output derivative of fin. 

  const double rdx2fac = pow(2.0/dx, 2);
  out[0] = 6.708203932499369*fin[5]*rdx2fac; 
  out[1] = 6.7082039324993685*fin[7]*rdx2fac; 
  out[2] = 0.0; 
  out[3] = 0.0; 
  out[4] = 6.708203932499369*fin[8]*rdx2fac; 
  out[5] = 0.0; 
  out[6] = 0.0; 
  out[7] = 0.0; 
  out[8] = 0.0; 

}

