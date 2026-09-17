#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[4]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[5]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[6]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[5]+fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[7]+3.0*coords[1]*coords[2]*fdo[6]+3.0*coords[0]*coords[2]*fdo[5]+3.0*coords[0]*coords[1]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[4]+fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[4]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[4]+fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[3]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[5]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_ser_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

