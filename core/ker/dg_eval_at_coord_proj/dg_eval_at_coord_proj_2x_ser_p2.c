#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[4]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[6]+10.392304845413262*coords[0]*fdo[3]+6.0*fdo[2]); 
  ftar[2] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[7]+5.0*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.1767766952966368*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[7]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[6]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[5]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[4]+6.0*coords[0]*coords[1]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[5]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[7]+10.392304845413262*coords[0]*fdo[3]+6.0*fdo[1]); 
  ftar[2] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[6]+5.0*fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p2_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

