#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[7]+4.898979485566357*coords[0]*fdo[1]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[11]+14.696938456699069*coords[0]*fdo[4]+8.485281374238571*fdo[2]); 
  ftar[2] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[13]+14.696938456699069*coords[0]*fdo[5]+8.485281374238571*fdo[3]); 
  ftar[3] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[17]+4.898979485566357*coords[0]*fdo[10]+2.8284271247461907*fdo[6]); 
  ftar[4] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[20]+24.49489742783179*coords[0]*fdo[12]+14.142135623730955*fdo[8]); 
  ftar[5] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[21]+24.49489742783179*coords[0]*fdo[15]+14.142135623730955*fdo[9]); 
  ftar[6] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[23]+73.48469228349536*coords[0]*fdo[18]+42.42640687119286*fdo[14]); 
  ftar[7] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[24]+73.48469228349536*coords[0]*fdo[19]+42.42640687119286*fdo[16]); 
  ftar[8] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[26]+4.898979485566357*coords[0]*fdo[25]+2.8284271247461907*fdo[22]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*(((45.0*coords0R2-15.0)*coords1R2-15.0*coords0R2+5.0)*fdo[20]+(23.237900077244504*coords[0]*coords1R2-7.745966692414834*coords[0])*fdo[12]+(23.237900077244504*coords0R2-7.745966692414834)*coords[1]*fdo[11]+(13.416407864998739*coords1R2-4.47213595499958)*fdo[8]+(13.416407864998739*coords0R2-4.47213595499958)*fdo[7]+12.0*coords[0]*coords[1]*fdo[4]+6.928203230275509*coords[1]*fdo[2]+6.928203230275509*coords[0]*fdo[1]+4.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*(((135.0*coords0R2-45.0)*coords1R2-45.0*coords0R2+15.0)*fdo[23]+(69.71370023173348*coords[0]*coords1R2-23.2379000772445*coords[0])*fdo[18]+(69.71370023173348*coords0R2-23.2379000772445)*coords[1]*fdo[17]+(40.24922359499621*coords1R2-13.416407864998737)*fdo[14]+(40.24922359499621*coords0R2-13.416407864998737)*fdo[13]+36.0*coords[0]*coords[1]*fdo[10]+20.784609690826525*coords[1]*fdo[6]+20.784609690826525*coords[0]*fdo[5]+12.0*fdo[3]); 
  ftar[2] = 0.025*(((225.0*coords0R2-75.0)*coords1R2-75.0*coords0R2+25.0)*fdo[26]+(116.18950038622252*coords[0]*coords1R2-38.729833462074176*coords[0])*fdo[25]+(116.18950038622252*coords0R2-38.729833462074176)*coords[1]*fdo[24]+(67.0820393249937*coords1R2-22.3606797749979)*fdo[22]+(67.0820393249937*coords0R2-22.3606797749979)*fdo[21]+60.0*coords[0]*coords[1]*fdo[19]+34.64101615137755*coords[1]*fdo[16]+34.64101615137755*coords[0]*fdo[15]+20.0*fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.03125*((((301.86917696247167*coords0R2-100.62305898749055)*coords1R2-100.62305898749055*coords0R2+33.54101966249685)*coords2R2+(33.54101966249685-100.62305898749055*coords0R2)*coords1R2+33.54101966249685*coords0R2-11.18033988749895)*fdo[26]+((155.8845726811989*coords[0]*coords1R2-51.96152422706631*coords[0])*coords2R2-51.96152422706631*coords[0]*coords1R2+17.32050807568877*coords[0])*fdo[25]+((155.8845726811989*coords0R2-51.96152422706631)*coords[1]*coords2R2+(17.32050807568877-51.96152422706631*coords0R2)*coords[1])*fdo[24]+((155.8845726811989*coords0R2-51.96152422706631)*coords1R2-51.96152422706631*coords0R2+17.32050807568877)*coords[2]*fdo[23]+((90.0*coords1R2-30.0)*coords2R2-30.0*coords1R2+10.0)*fdo[22]+((90.0*coords0R2-30.0)*coords2R2-30.0*coords0R2+10.0)*fdo[21]+((90.0*coords0R2-30.0)*coords1R2-30.0*coords0R2+10.0)*fdo[20]+(80.49844718999243*coords[0]*coords[1]*coords2R2-26.832815729997478*coords[0]*coords[1])*fdo[19]+(80.49844718999243*coords[0]*coords1R2-26.832815729997478*coords[0])*coords[2]*fdo[18]+(80.49844718999243*coords0R2-26.832815729997478)*coords[1]*coords[2]*fdo[17]+(46.47580015448901*coords[1]*coords2R2-15.491933384829668*coords[1])*fdo[16]+(46.47580015448901*coords[0]*coords2R2-15.491933384829668*coords[0])*fdo[15]+(46.47580015448901*coords1R2-15.491933384829668)*coords[2]*fdo[14]+(46.47580015448901*coords0R2-15.491933384829668)*coords[2]*fdo[13]+(46.47580015448901*coords[0]*coords1R2-15.491933384829668*coords[0])*fdo[12]+(46.47580015448901*coords0R2-15.491933384829668)*coords[1]*fdo[11]+41.56921938165305*coords[0]*coords[1]*coords[2]*fdo[10]+(26.832815729997478*coords2R2-8.94427190999916)*fdo[9]+(26.832815729997478*coords1R2-8.94427190999916)*fdo[8]+(26.832815729997478*coords0R2-8.94427190999916)*fdo[7]+24.0*coords[1]*coords[2]*fdo[6]+24.0*coords[0]*coords[2]*fdo[5]+24.0*coords[0]*coords[1]*fdo[4]+13.856406460551018*coords[2]*fdo[3]+13.856406460551018*coords[1]*fdo[2]+13.856406460551018*coords[0]*fdo[1]+8.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*(((45.0*coords0R2-15.0)*coords1R2-15.0*coords0R2+5.0)*fdo[21]+(23.237900077244504*coords[0]*coords1R2-7.745966692414834*coords[0])*fdo[15]+(23.237900077244504*coords0R2-7.745966692414834)*coords[1]*fdo[13]+(13.416407864998739*coords1R2-4.47213595499958)*fdo[9]+(13.416407864998739*coords0R2-4.47213595499958)*fdo[7]+12.0*coords[0]*coords[1]*fdo[5]+6.928203230275509*coords[1]*fdo[3]+6.928203230275509*coords[0]*fdo[1]+4.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*(((135.0*coords0R2-45.0)*coords1R2-45.0*coords0R2+15.0)*fdo[24]+(69.71370023173348*coords[0]*coords1R2-23.2379000772445*coords[0])*fdo[19]+(69.71370023173348*coords0R2-23.2379000772445)*coords[1]*fdo[17]+(40.24922359499621*coords1R2-13.416407864998737)*fdo[16]+(40.24922359499621*coords0R2-13.416407864998737)*fdo[11]+36.0*coords[0]*coords[1]*fdo[10]+20.784609690826525*coords[1]*fdo[6]+20.784609690826525*coords[0]*fdo[4]+12.0*fdo[2]); 
  ftar[2] = 0.025*(((225.0*coords0R2-75.0)*coords1R2-75.0*coords0R2+25.0)*fdo[26]+(116.18950038622252*coords[0]*coords1R2-38.729833462074176*coords[0])*fdo[25]+(116.18950038622252*coords0R2-38.729833462074176)*coords[1]*fdo[23]+(67.0820393249937*coords1R2-22.3606797749979)*fdo[22]+(67.0820393249937*coords0R2-22.3606797749979)*fdo[20]+60.0*coords[0]*coords[1]*fdo[18]+34.64101615137755*coords[1]*fdo[14]+34.64101615137755*coords[0]*fdo[12]+20.0*fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[8]+4.898979485566357*coords[0]*fdo[2]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[12]+14.696938456699069*coords[0]*fdo[4]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[14]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[3]); 
  ftar[3] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[18]+4.898979485566357*coords[0]*fdo[10]+2.8284271247461907*fdo[5]); 
  ftar[4] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[20]+24.49489742783179*coords[0]*fdo[11]+14.142135623730955*fdo[7]); 
  ftar[5] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[22]+24.49489742783179*coords[0]*fdo[16]+14.142135623730955*fdo[9]); 
  ftar[6] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[23]+73.48469228349536*coords[0]*fdo[17]+42.42640687119286*fdo[13]); 
  ftar[7] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[25]+73.48469228349536*coords[0]*fdo[19]+42.42640687119286*fdo[15]); 
  ftar[8] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[26]+4.898979485566357*coords[0]*fdo[24]+2.8284271247461907*fdo[21]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*(((45.0*coords0R2-15.0)*coords1R2-15.0*coords0R2+5.0)*fdo[22]+(23.237900077244504*coords[0]*coords1R2-7.745966692414834*coords[0])*fdo[16]+(23.237900077244504*coords0R2-7.745966692414834)*coords[1]*fdo[14]+(13.416407864998739*coords1R2-4.47213595499958)*fdo[9]+(13.416407864998739*coords0R2-4.47213595499958)*fdo[8]+12.0*coords[0]*coords[1]*fdo[6]+6.928203230275509*coords[1]*fdo[3]+6.928203230275509*coords[0]*fdo[2]+4.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*(((135.0*coords0R2-45.0)*coords1R2-45.0*coords0R2+15.0)*fdo[25]+(69.71370023173348*coords[0]*coords1R2-23.2379000772445*coords[0])*fdo[19]+(69.71370023173348*coords0R2-23.2379000772445)*coords[1]*fdo[18]+(40.24922359499621*coords1R2-13.416407864998737)*fdo[15]+(40.24922359499621*coords0R2-13.416407864998737)*fdo[12]+36.0*coords[0]*coords[1]*fdo[10]+20.784609690826525*coords[1]*fdo[5]+20.784609690826525*coords[0]*fdo[4]+12.0*fdo[1]); 
  ftar[2] = 0.025*(((225.0*coords0R2-75.0)*coords1R2-75.0*coords0R2+25.0)*fdo[26]+(116.18950038622252*coords[0]*coords1R2-38.729833462074176*coords[0])*fdo[24]+(116.18950038622252*coords0R2-38.729833462074176)*coords[1]*fdo[23]+(67.0820393249937*coords1R2-22.3606797749979)*fdo[21]+(67.0820393249937*coords0R2-22.3606797749979)*fdo[20]+60.0*coords[0]*coords[1]*fdo[17]+34.64101615137755*coords[1]*fdo[13]+34.64101615137755*coords[0]*fdo[11]+20.0*fdo[7]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[9]+4.898979485566357*coords[0]*fdo[3]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[15]+14.696938456699069*coords[0]*fdo[5]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[16]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[2]); 
  ftar[3] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[19]+4.898979485566357*coords[0]*fdo[10]+2.8284271247461907*fdo[4]); 
  ftar[4] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[21]+24.49489742783179*coords[0]*fdo[13]+14.142135623730955*fdo[7]); 
  ftar[5] = 0.05*((47.43416490252571*coords0R2-15.811388300841902)*fdo[22]+24.49489742783179*coords[0]*fdo[14]+14.142135623730955*fdo[8]); 
  ftar[6] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[24]+73.48469228349536*coords[0]*fdo[17]+42.42640687119286*fdo[11]); 
  ftar[7] = 0.016666666666666666*((142.30249470757707*coords0R2-47.434164902525694)*fdo[25]+73.48469228349536*coords[0]*fdo[18]+42.42640687119286*fdo[12]); 
  ftar[8] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[26]+4.898979485566357*coords[0]*fdo[23]+2.8284271247461907*fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 9;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 9;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x_tensor_p2_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_TENSOR;
  *poly_order = 2;
  *num_basis = 9;
}

