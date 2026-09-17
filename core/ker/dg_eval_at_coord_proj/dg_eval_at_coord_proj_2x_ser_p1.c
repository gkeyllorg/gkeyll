#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[3]+fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(3.0*coords[0]*coords[1]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[3]+fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x_ser_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

