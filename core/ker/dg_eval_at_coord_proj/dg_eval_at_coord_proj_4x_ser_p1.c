#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[5]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[11]+fdo[6]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[12]+fdo[8]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[14]+fdo[10]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[13]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[5]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[12]+fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[14]+3.0*coords[1]*coords[2]*fdo[10]+3.0*coords[0]*coords[2]*fdo[9]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[15]+3.0*coords[1]*coords[2]*fdo[13]+3.0*coords[0]*coords[2]*fdo[12]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[8]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[5]+fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[15]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[14]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[13]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[12]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[11]+3.0*coords[2]*coords[3]*fdo[10]+3.0*coords[1]*coords[3]*fdo[9]+3.0*coords[0]*coords[3]*fdo[8]+3.0*coords[1]*coords[2]*fdo[7]+3.0*coords[0]*coords[2]*fdo[6]+3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[11]+3.0*coords[1]*coords[2]*fdo[7]+3.0*coords[0]*coords[2]*fdo[6]+3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[15]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[5]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[11]+fdo[6]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[12]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[15]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[10]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[11]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[13]+fdo[8]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[14]+fdo[9]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[11]+fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[13]+3.0*coords[1]*coords[2]*fdo[10]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[15]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[12]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[9]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[5]+fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[5]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[12]+fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[10]+fdo[3]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[12]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[13]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[14]+fdo[7]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[11]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[5]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[11]+fdo[7]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[5]+fdo[2]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[8]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[11]+fdo[7]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[12]+fdo[9]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[13]+fdo[10]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[14]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

