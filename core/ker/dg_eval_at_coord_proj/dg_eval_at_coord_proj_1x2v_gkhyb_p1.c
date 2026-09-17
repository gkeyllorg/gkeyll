#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[3]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[5]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[4]); 
  ftar[4] = 0.1*(12.247448713915892*coords[0]*fdo[10]+7.0710678118654755*fdo[8]); 
  ftar[5] = 0.1*(12.247448713915892*coords[0]*fdo[11]+7.0710678118654755*fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[10]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[8]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[11]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[9]+18.0*coords[0]*coords[1]*fdo[7]+10.392304845413262*coords[1]*fdo[5]+10.392304845413262*coords[0]*fdo[4]+6.0*fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*((20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[11]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[10]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[9]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[8]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[7]+6.0*coords[1]*coords[2]*fdo[6]+6.0*coords[0]*coords[2]*fdo[5]+6.0*coords[0]*coords[1]*fdo[4]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[4]+fdo[2]); 
  ftar[2] = 0.1*(15.0*coords[0]*coords[1]*fdo[11]+8.660254037844387*coords[1]*fdo[10]+8.660254037844387*coords[0]*fdo[9]+5.0*fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[8]+4.898979485566357*coords[0]*fdo[2]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[9]+14.696938456699069*coords[0]*fdo[4]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[10]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[3]); 
  ftar[3] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[11]+4.898979485566357*coords[0]*fdo[7]+2.8284271247461907*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[9]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[8]+6.0*coords[0]*coords[1]*fdo[4]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[11]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[10]+18.0*coords[0]*coords[1]*fdo[7]+10.392304845413262*coords[1]*fdo[6]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[4]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[5]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[6]); 
  ftar[4] = 0.1*(12.247448713915892*coords[0]*fdo[9]+7.0710678118654755*fdo[8]); 
  ftar[5] = 0.1*(12.247448713915892*coords[0]*fdo[11]+7.0710678118654755*fdo[10]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_1x2v_gkhyb_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_GKHYBRID_VEL;
  *poly_order = 1;
  *num_basis = 6;
}

