#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[4]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[6]+10.392304845413262*coords[0]*fdo[3]+6.0*fdo[2]); 
  ftar[2] = 0.07071067811865474*((33.54101966249685*coords0R2-11.18033988749895)*fdo[8]+17.320508075688775*coords[0]*fdo[7]+10.0*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.0883883476483184*(((45.0*coords0R2-15.0)*coords1R2-15.0*coords0R2+5.0)*fdo[8]+(23.237900077244504*coords[0]*coords1R2-7.745966692414834*coords[0])*fdo[7]+(23.237900077244504*coords0R2-7.745966692414834)*coords[1]*fdo[6]+(13.416407864998739*coords1R2-4.47213595499958)*fdo[5]+(13.416407864998739*coords0R2-4.47213595499958)*fdo[4]+12.0*coords[0]*coords[1]*fdo[3]+6.928203230275509*coords[1]*fdo[2]+6.928203230275509*coords[0]*fdo[1]+4.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[5]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[7]+10.392304845413262*coords[0]*fdo[3]+6.0*fdo[1]); 
  ftar[2] = 0.07071067811865474*((33.54101966249685*coords0R2-11.18033988749895)*fdo[8]+17.320508075688775*coords[0]*fdo[6]+10.0*fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_tensor_p2_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 3;
}

