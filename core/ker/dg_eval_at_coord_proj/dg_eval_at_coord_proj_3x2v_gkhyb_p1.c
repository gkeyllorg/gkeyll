#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[5]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[13]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[14]+1.4142135623730951*fdo[3]); 
  ftar[4] = 0.5*(2.4494897427831783*coords[0]*fdo[15]+1.4142135623730951*fdo[4]); 
  ftar[5] = 0.5*(2.4494897427831783*coords[0]*fdo[20]+1.4142135623730951*fdo[6]); 
  ftar[6] = 0.5*(2.4494897427831783*coords[0]*fdo[21]+1.4142135623730951*fdo[7]); 
  ftar[7] = 0.5*(2.4494897427831783*coords[0]*fdo[22]+1.4142135623730951*fdo[8]); 
  ftar[8] = 0.5*(2.4494897427831783*coords[0]*fdo[23]+1.4142135623730951*fdo[9]); 
  ftar[9] = 0.5*(2.4494897427831783*coords[0]*fdo[24]+1.4142135623730951*fdo[10]); 
  ftar[10] = 0.5*(2.4494897427831783*coords[0]*fdo[25]+1.4142135623730951*fdo[11]); 
  ftar[11] = 0.5*(2.4494897427831783*coords[0]*fdo[27]+1.4142135623730951*fdo[16]); 
  ftar[12] = 0.5*(2.4494897427831783*coords[0]*fdo[28]+1.4142135623730951*fdo[17]); 
  ftar[13] = 0.5*(2.4494897427831783*coords[0]*fdo[29]+1.4142135623730951*fdo[18]); 
  ftar[14] = 0.5*(2.4494897427831783*coords[0]*fdo[30]+1.4142135623730951*fdo[19]); 
  ftar[15] = 0.5*(2.4494897427831783*coords[0]*fdo[31]+1.4142135623730951*fdo[26]); 
  ftar[16] = 0.1*(12.247448713915892*coords[0]*fdo[36]+7.0710678118654755*fdo[32]); 
  ftar[17] = 0.1*(12.247448713915892*coords[0]*fdo[40]+7.0710678118654755*fdo[33]); 
  ftar[18] = 0.1*(12.247448713915892*coords[0]*fdo[41]+7.0710678118654755*fdo[34]); 
  ftar[19] = 0.1*(12.247448713915892*coords[0]*fdo[42]+7.0710678118654755*fdo[35]); 
  ftar[20] = 0.1*(12.247448713915892*coords[0]*fdo[44]+7.0710678118654755*fdo[37]); 
  ftar[21] = 0.1*(12.247448713915892*coords[0]*fdo[45]+7.0710678118654755*fdo[38]); 
  ftar[22] = 0.1*(12.247448713915892*coords[0]*fdo[46]+7.0710678118654755*fdo[39]); 
  ftar[23] = 0.1*(12.247448713915892*coords[0]*fdo[47]+7.0710678118654755*fdo[43]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[36]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[32]+6.0*coords[0]*coords[1]*fdo[15]+3.4641016151377544*coords[1]*fdo[5]+3.4641016151377544*coords[0]*fdo[4]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[40]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[33]+18.0*coords[0]*coords[1]*fdo[23]+10.392304845413262*coords[1]*fdo[12]+10.392304845413262*coords[0]*fdo[9]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[41]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[34]+18.0*coords[0]*coords[1]*fdo[24]+10.392304845413262*coords[1]*fdo[13]+10.392304845413262*coords[0]*fdo[10]+6.0*fdo[2]); 
  ftar[3] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[42]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[35]+18.0*coords[0]*coords[1]*fdo[25]+10.392304845413262*coords[1]*fdo[14]+10.392304845413262*coords[0]*fdo[11]+6.0*fdo[3]); 
  ftar[4] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[44]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[37]+6.0*coords[0]*coords[1]*fdo[28]+3.4641016151377544*coords[1]*fdo[20]+3.4641016151377544*coords[0]*fdo[17]+2.0*fdo[6]); 
  ftar[5] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[45]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[38]+6.0*coords[0]*coords[1]*fdo[29]+3.4641016151377544*coords[1]*fdo[21]+3.4641016151377544*coords[0]*fdo[18]+2.0*fdo[7]); 
  ftar[6] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[46]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[39]+6.0*coords[0]*coords[1]*fdo[30]+3.4641016151377544*coords[1]*fdo[22]+3.4641016151377544*coords[0]*fdo[19]+2.0*fdo[8]); 
  ftar[7] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[47]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[43]+18.0*coords[0]*coords[1]*fdo[31]+10.392304845413262*coords[1]*fdo[27]+10.392304845413262*coords[0]*fdo[26]+6.0*fdo[16]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[40]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[36]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[33]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[23]+8.485281374238571*coords[1]*coords[2]*fdo[15]+8.485281374238571*coords[0]*coords[2]*fdo[12]+8.485281374238571*coords[0]*coords[1]*fdo[9]+4.898979485566357*coords[2]*fdo[5]+4.898979485566357*coords[1]*fdo[4]+4.898979485566357*coords[0]*fdo[1]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[44]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[41]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[37]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[34]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[28]+25.455844122715718*coords[1]*coords[2]*fdo[24]+25.455844122715718*coords[0]*coords[2]*fdo[20]+25.455844122715718*coords[0]*coords[1]*fdo[17]+14.696938456699069*coords[2]*fdo[13]+14.696938456699069*coords[1]*fdo[10]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[2]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[45]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[42]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[38]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[35]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[29]+25.455844122715718*coords[1]*coords[2]*fdo[25]+25.455844122715718*coords[0]*coords[2]*fdo[21]+25.455844122715718*coords[0]*coords[1]*fdo[18]+14.696938456699069*coords[2]*fdo[14]+14.696938456699069*coords[1]*fdo[11]+14.696938456699069*coords[0]*fdo[7]+8.485281374238571*fdo[3]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[47]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[46]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[43]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[39]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[30]+8.485281374238571*coords[0]*coords[2]*fdo[27]+8.485281374238571*coords[0]*coords[1]*fdo[26]+4.898979485566357*coords[2]*fdo[22]+4.898979485566357*coords[1]*fdo[19]+4.898979485566357*coords[0]*fdo[16]+2.8284271247461907*fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((34.856850115866756*coords[0]*coords[1]*coords2R2-11.618950038622252*coords[0]*coords[1])*coords[3]*fdo[44]+(20.12461179749811*coords[1]*coords2R2-6.708203932499369*coords[1])*coords[3]*fdo[41]+(20.12461179749811*coords[0]*coords2R2-6.708203932499369*coords[0])*coords[3]*fdo[40]+(20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[37]+(11.618950038622252*coords2R2-3.872983346207417)*coords[3]*fdo[36]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[34]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[33]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[32]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[28]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[24]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[23]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[20]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[17]+6.0*coords[2]*coords[3]*fdo[15]+6.0*coords[1]*coords[3]*fdo[13]+6.0*coords[0]*coords[3]*fdo[12]+6.0*coords[1]*coords[2]*fdo[10]+6.0*coords[0]*coords[2]*fdo[9]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[3]*fdo[5]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*((104.57055034760022*coords[0]*coords[1]*coords2R2-34.85685011586674*coords[0]*coords[1])*coords[3]*fdo[47]+(60.37383539249431*coords[1]*coords2R2-20.124611797498105*coords[1])*coords[3]*fdo[46]+(60.37383539249431*coords[0]*coords2R2-20.124611797498105*coords[0])*coords[3]*fdo[45]+(60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[43]+(34.85685011586674*coords2R2-11.61895003862225)*coords[3]*fdo[42]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[39]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[38]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[35]+54.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+31.17691453623978*coords[1]*coords[2]*coords[3]*fdo[30]+31.17691453623978*coords[0]*coords[2]*coords[3]*fdo[29]+31.17691453623978*coords[0]*coords[1]*coords[3]*fdo[27]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[26]+18.0*coords[2]*coords[3]*fdo[25]+18.0*coords[1]*coords[3]*fdo[22]+18.0*coords[0]*coords[3]*fdo[21]+18.0*coords[1]*coords[2]*fdo[19]+18.0*coords[0]*coords[2]*fdo[18]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[3]*fdo[14]+10.392304845413262*coords[2]*fdo[11]+10.392304845413262*coords[1]*fdo[8]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords3R2 = pow(coords[3],2);

  ftar[0] = 0.0625*((60.37383539249433*coords[0]*coords[1]*coords[2]*coords3R2-20.12461179749811*coords[0]*coords[1]*coords[2])*coords[4]*fdo[47]+(34.856850115866756*coords[1]*coords[2]*coords3R2-11.618950038622252*coords[1]*coords[2])*coords[4]*fdo[46]+(34.856850115866756*coords[0]*coords[2]*coords3R2-11.618950038622252*coords[0]*coords[2])*coords[4]*fdo[45]+(34.856850115866756*coords[0]*coords[1]*coords3R2-11.618950038622252*coords[0]*coords[1])*coords[4]*fdo[44]+(34.856850115866756*coords[0]*coords[1]*coords[2]*coords3R2-11.618950038622252*coords[0]*coords[1]*coords[2])*fdo[43]+(20.12461179749811*coords[2]*coords3R2-6.708203932499369*coords[2])*coords[4]*fdo[42]+(20.12461179749811*coords[1]*coords3R2-6.708203932499369*coords[1])*coords[4]*fdo[41]+(20.12461179749811*coords[0]*coords3R2-6.708203932499369*coords[0])*coords[4]*fdo[40]+(20.12461179749811*coords[1]*coords[2]*coords3R2-6.708203932499369*coords[1]*coords[2])*fdo[39]+(20.12461179749811*coords[0]*coords[2]*coords3R2-6.708203932499369*coords[0]*coords[2])*fdo[38]+(20.12461179749811*coords[0]*coords[1]*coords3R2-6.708203932499369*coords[0]*coords[1])*fdo[37]+(11.618950038622252*coords3R2-3.872983346207417)*coords[4]*fdo[36]+(11.618950038622252*coords[2]*coords3R2-3.872983346207417*coords[2])*fdo[35]+(11.618950038622252*coords[1]*coords3R2-3.872983346207417*coords[1])*fdo[34]+(11.618950038622252*coords[0]*coords3R2-3.872983346207417*coords[0])*fdo[33]+(6.708203932499369*coords3R2-2.23606797749979)*fdo[32]+31.17691453623978*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[31]+18.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[30]+18.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[29]+18.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[28]+18.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[27]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[26]+10.392304845413262*coords[2]*coords[3]*coords[4]*fdo[25]+10.392304845413262*coords[1]*coords[3]*coords[4]*fdo[24]+10.392304845413262*coords[0]*coords[3]*coords[4]*fdo[23]+10.392304845413262*coords[1]*coords[2]*coords[4]*fdo[22]+10.392304845413262*coords[0]*coords[2]*coords[4]*fdo[21]+10.392304845413262*coords[0]*coords[1]*coords[4]*fdo[20]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[19]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[18]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[17]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[16]+6.0*coords[3]*coords[4]*fdo[15]+6.0*coords[2]*coords[4]*fdo[14]+6.0*coords[1]*coords[4]*fdo[13]+6.0*coords[0]*coords[4]*fdo[12]+6.0*coords[2]*coords[3]*fdo[11]+6.0*coords[1]*coords[3]*fdo[10]+6.0*coords[0]*coords[3]*fdo[9]+6.0*coords[1]*coords[2]*fdo[8]+6.0*coords[0]*coords[2]*fdo[7]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[4]*fdo[5]+3.4641016151377544*coords[3]*fdo[4]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((34.856850115866756*coords[0]*coords[1]*coords2R2-11.618950038622252*coords[0]*coords[1])*coords[3]*fdo[45]+(20.12461179749811*coords[1]*coords2R2-6.708203932499369*coords[1])*coords[3]*fdo[42]+(20.12461179749811*coords[0]*coords2R2-6.708203932499369*coords[0])*coords[3]*fdo[40]+(20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[38]+(11.618950038622252*coords2R2-3.872983346207417)*coords[3]*fdo[36]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[35]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[33]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[32]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[29]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[25]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[23]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[21]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[18]+6.0*coords[2]*coords[3]*fdo[15]+6.0*coords[1]*coords[3]*fdo[14]+6.0*coords[0]*coords[3]*fdo[12]+6.0*coords[1]*coords[2]*fdo[11]+6.0*coords[0]*coords[2]*fdo[9]+6.0*coords[0]*coords[1]*fdo[7]+3.4641016151377544*coords[3]*fdo[5]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*((104.57055034760022*coords[0]*coords[1]*coords2R2-34.85685011586674*coords[0]*coords[1])*coords[3]*fdo[47]+(60.37383539249431*coords[1]*coords2R2-20.124611797498105*coords[1])*coords[3]*fdo[46]+(60.37383539249431*coords[0]*coords2R2-20.124611797498105*coords[0])*coords[3]*fdo[44]+(60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[43]+(34.85685011586674*coords2R2-11.61895003862225)*coords[3]*fdo[41]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[39]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[37]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[34]+54.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+31.17691453623978*coords[1]*coords[2]*coords[3]*fdo[30]+31.17691453623978*coords[0]*coords[2]*coords[3]*fdo[28]+31.17691453623978*coords[0]*coords[1]*coords[3]*fdo[27]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[26]+18.0*coords[2]*coords[3]*fdo[24]+18.0*coords[1]*coords[3]*fdo[22]+18.0*coords[0]*coords[3]*fdo[20]+18.0*coords[1]*coords[2]*fdo[19]+18.0*coords[0]*coords[2]*fdo[17]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[3]*fdo[13]+10.392304845413262*coords[2]*fdo[10]+10.392304845413262*coords[1]*fdo[8]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[41]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[36]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[34]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[24]+8.485281374238571*coords[1]*coords[2]*fdo[15]+8.485281374238571*coords[0]*coords[2]*fdo[13]+8.485281374238571*coords[0]*coords[1]*fdo[10]+4.898979485566357*coords[2]*fdo[5]+4.898979485566357*coords[1]*fdo[4]+4.898979485566357*coords[0]*fdo[2]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[44]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[40]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[37]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[33]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[28]+25.455844122715718*coords[1]*coords[2]*fdo[23]+25.455844122715718*coords[0]*coords[2]*fdo[20]+25.455844122715718*coords[0]*coords[1]*fdo[17]+14.696938456699069*coords[2]*fdo[12]+14.696938456699069*coords[1]*fdo[9]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[46]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[42]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[39]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[35]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[30]+25.455844122715718*coords[1]*coords[2]*fdo[25]+25.455844122715718*coords[0]*coords[2]*fdo[22]+25.455844122715718*coords[0]*coords[1]*fdo[19]+14.696938456699069*coords[2]*fdo[14]+14.696938456699069*coords[1]*fdo[11]+14.696938456699069*coords[0]*fdo[8]+8.485281374238571*fdo[3]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[47]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[45]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[43]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[38]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[29]+8.485281374238571*coords[0]*coords[2]*fdo[27]+8.485281374238571*coords[0]*coords[1]*fdo[26]+4.898979485566357*coords[2]*fdo[21]+4.898979485566357*coords[1]*fdo[18]+4.898979485566357*coords[0]*fdo[16]+2.8284271247461907*fdo[7]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((34.856850115866756*coords[0]*coords[1]*coords2R2-11.618950038622252*coords[0]*coords[1])*coords[3]*fdo[46]+(20.12461179749811*coords[1]*coords2R2-6.708203932499369*coords[1])*coords[3]*fdo[42]+(20.12461179749811*coords[0]*coords2R2-6.708203932499369*coords[0])*coords[3]*fdo[41]+(20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[39]+(11.618950038622252*coords2R2-3.872983346207417)*coords[3]*fdo[36]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[35]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[34]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[32]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[30]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[25]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[24]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[22]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[19]+6.0*coords[2]*coords[3]*fdo[15]+6.0*coords[1]*coords[3]*fdo[14]+6.0*coords[0]*coords[3]*fdo[13]+6.0*coords[1]*coords[2]*fdo[11]+6.0*coords[0]*coords[2]*fdo[10]+6.0*coords[0]*coords[1]*fdo[8]+3.4641016151377544*coords[3]*fdo[5]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*((104.57055034760022*coords[0]*coords[1]*coords2R2-34.85685011586674*coords[0]*coords[1])*coords[3]*fdo[47]+(60.37383539249431*coords[1]*coords2R2-20.124611797498105*coords[1])*coords[3]*fdo[45]+(60.37383539249431*coords[0]*coords2R2-20.124611797498105*coords[0])*coords[3]*fdo[44]+(60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[43]+(34.85685011586674*coords2R2-11.61895003862225)*coords[3]*fdo[40]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[38]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[37]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[33]+54.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+31.17691453623978*coords[1]*coords[2]*coords[3]*fdo[29]+31.17691453623978*coords[0]*coords[2]*coords[3]*fdo[28]+31.17691453623978*coords[0]*coords[1]*coords[3]*fdo[27]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[26]+18.0*coords[2]*coords[3]*fdo[23]+18.0*coords[1]*coords[3]*fdo[21]+18.0*coords[0]*coords[3]*fdo[20]+18.0*coords[1]*coords[2]*fdo[18]+18.0*coords[0]*coords[2]*fdo[17]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[3]*fdo[12]+10.392304845413262*coords[2]*fdo[9]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[42]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[36]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[35]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[25]+8.485281374238571*coords[1]*coords[2]*fdo[15]+8.485281374238571*coords[0]*coords[2]*fdo[14]+8.485281374238571*coords[0]*coords[1]*fdo[11]+4.898979485566357*coords[2]*fdo[5]+4.898979485566357*coords[1]*fdo[4]+4.898979485566357*coords[0]*fdo[3]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[45]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[40]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[38]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[33]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[29]+25.455844122715718*coords[1]*coords[2]*fdo[23]+25.455844122715718*coords[0]*coords[2]*fdo[21]+25.455844122715718*coords[0]*coords[1]*fdo[18]+14.696938456699069*coords[2]*fdo[12]+14.696938456699069*coords[1]*fdo[9]+14.696938456699069*coords[0]*fdo[7]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords1R2-28.46049894151541*coords[0])*coords[2]*fdo[46]+(49.295030175464944*coords1R2-16.431676725154983)*coords[2]*fdo[41]+(49.295030175464944*coords[0]*coords1R2-16.431676725154983*coords[0])*fdo[39]+(28.46049894151541*coords1R2-9.48683298050514)*fdo[34]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[30]+25.455844122715718*coords[1]*coords[2]*fdo[24]+25.455844122715718*coords[0]*coords[2]*fdo[22]+25.455844122715718*coords[0]*coords[1]*fdo[19]+14.696938456699069*coords[2]*fdo[13]+14.696938456699069*coords[1]*fdo[10]+14.696938456699069*coords[0]*fdo[8]+8.485281374238571*fdo[2]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords1R2-9.48683298050514*coords[0])*coords[2]*fdo[47]+(16.431676725154986*coords1R2-5.477225575051662)*coords[2]*fdo[44]+(16.431676725154986*coords[0]*coords1R2-5.477225575051662*coords[0])*fdo[43]+(9.48683298050514*coords1R2-3.1622776601683795)*fdo[37]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[28]+8.485281374238571*coords[0]*coords[2]*fdo[27]+8.485281374238571*coords[0]*coords[1]*fdo[26]+4.898979485566357*coords[2]*fdo[20]+4.898979485566357*coords[1]*fdo[17]+4.898979485566357*coords[0]*fdo[16]+2.8284271247461907*fdo[6]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[20]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[6]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[21]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[22]+1.7320508075688772*coords[0]*fdo[16]+fdo[8]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[17]+fdo[10]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[18]+fdo[11]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[26]+fdo[19]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[40]+8.660254037844387*coords[1]*fdo[36]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[41]+8.660254037844387*coords[0]*fdo[37]+5.0*fdo[34]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[42]+8.660254037844387*coords[0]*fdo[38]+5.0*fdo[35]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[46]+8.660254037844387*coords[0]*fdo[43]+5.0*fdo[39]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[20]+4.242640687119286*coords[1]*coords[2]*fdo[13]+4.242640687119286*coords[0]*coords[2]*fdo[12]+4.242640687119286*coords[0]*coords[1]*fdo[6]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[2]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[3]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[25]+2.4494897427831783*coords[1]*fdo[19]+2.4494897427831783*coords[0]*fdo[18]+1.4142135623730951*fdo[11]); 
  ftar[4] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[44]+21.213203435596427*coords[1]*coords[2]*fdo[41]+21.213203435596427*coords[0]*coords[2]*fdo[40]+21.213203435596427*coords[0]*coords[1]*fdo[37]+12.247448713915892*coords[2]*fdo[36]+12.247448713915892*coords[1]*fdo[34]+12.247448713915892*coords[0]*fdo[33]+7.0710678118654755*fdo[32]); 
  ftar[5] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[47]+21.213203435596427*coords[1]*coords[2]*fdo[46]+21.213203435596427*coords[0]*coords[2]*fdo[45]+21.213203435596427*coords[0]*coords[1]*fdo[43]+12.247448713915892*coords[2]*fdo[42]+12.247448713915892*coords[1]*fdo[39]+12.247448713915892*coords[0]*fdo[38]+7.0710678118654755*fdo[35]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[27]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[22]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[21]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[20]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[16]+3.0*coords[2]*coords[3]*fdo[14]+3.0*coords[1]*coords[3]*fdo[13]+3.0*coords[0]*coords[3]*fdo[12]+3.0*coords[1]*coords[2]*fdo[8]+3.0*coords[0]*coords[2]*fdo[7]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[25]+3.0*coords[1]*coords[3]*fdo[24]+3.0*coords[0]*coords[3]*fdo[23]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[3]*fdo[15]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[2] = 0.05*(45.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[47]+25.98076211353316*coords[1]*coords[2]*coords[3]*fdo[46]+25.98076211353316*coords[0]*coords[2]*coords[3]*fdo[45]+25.98076211353316*coords[0]*coords[1]*coords[3]*fdo[44]+25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[43]+15.0*coords[2]*coords[3]*fdo[42]+15.0*coords[1]*coords[3]*fdo[41]+15.0*coords[0]*coords[3]*fdo[40]+15.0*coords[1]*coords[2]*fdo[39]+15.0*coords[0]*coords[2]*fdo[38]+15.0*coords[0]*coords[1]*fdo[37]+8.660254037844387*coords[3]*fdo[36]+8.660254037844387*coords[2]*fdo[35]+8.660254037844387*coords[1]*fdo[34]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[32]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[21]+4.242640687119286*coords[1]*coords[2]*fdo[14]+4.242640687119286*coords[0]*coords[2]*fdo[12]+4.242640687119286*coords[0]*coords[1]*fdo[7]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[13]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[18]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[24]+2.4494897427831783*coords[1]*fdo[19]+2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[10]); 
  ftar[4] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[45]+21.213203435596427*coords[1]*coords[2]*fdo[42]+21.213203435596427*coords[0]*coords[2]*fdo[40]+21.213203435596427*coords[0]*coords[1]*fdo[38]+12.247448713915892*coords[2]*fdo[36]+12.247448713915892*coords[1]*fdo[35]+12.247448713915892*coords[0]*fdo[33]+7.0710678118654755*fdo[32]); 
  ftar[5] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[47]+21.213203435596427*coords[1]*coords[2]*fdo[46]+21.213203435596427*coords[0]*coords[2]*fdo[44]+21.213203435596427*coords[0]*coords[1]*fdo[43]+12.247448713915892*coords[2]*fdo[41]+12.247448713915892*coords[1]*fdo[39]+12.247448713915892*coords[0]*fdo[37]+7.0710678118654755*fdo[34]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[20]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[16]+fdo[7]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[17]+fdo[9]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[19]+fdo[11]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[26]+fdo[18]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[41]+8.660254037844387*coords[1]*fdo[36]+8.660254037844387*coords[0]*fdo[34]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[40]+8.660254037844387*coords[0]*fdo[37]+5.0*fdo[33]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[42]+8.660254037844387*coords[0]*fdo[39]+5.0*fdo[35]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[45]+8.660254037844387*coords[0]*fdo[43]+5.0*fdo[38]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[1]*coords[2]*fdo[14]+4.242640687119286*coords[0]*coords[2]*fdo[13]+4.242640687119286*coords[0]*coords[1]*fdo[8]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[12]+2.4494897427831783*coords[1]*fdo[7]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[1]*fdo[19]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[10]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[23]+2.4494897427831783*coords[1]*fdo[18]+2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[9]); 
  ftar[4] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[46]+21.213203435596427*coords[1]*coords[2]*fdo[42]+21.213203435596427*coords[0]*coords[2]*fdo[41]+21.213203435596427*coords[0]*coords[1]*fdo[39]+12.247448713915892*coords[2]*fdo[36]+12.247448713915892*coords[1]*fdo[35]+12.247448713915892*coords[0]*fdo[34]+7.0710678118654755*fdo[32]); 
  ftar[5] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[47]+21.213203435596427*coords[1]*coords[2]*fdo[45]+21.213203435596427*coords[0]*coords[2]*fdo[44]+21.213203435596427*coords[0]*coords[1]*fdo[43]+12.247448713915892*coords[2]*fdo[40]+12.247448713915892*coords[1]*fdo[38]+12.247448713915892*coords[0]*fdo[37]+7.0710678118654755*fdo[33]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[21]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[8]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[16]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[18]+fdo[9]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[19]+fdo[10]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[26]+fdo[17]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[42]+8.660254037844387*coords[1]*fdo[36]+8.660254037844387*coords[0]*fdo[35]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[40]+8.660254037844387*coords[0]*fdo[38]+5.0*fdo[33]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[41]+8.660254037844387*coords[0]*fdo[39]+5.0*fdo[34]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[44]+8.660254037844387*coords[0]*fdo[43]+5.0*fdo[37]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[32]+4.898979485566357*coords[0]*fdo[4]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[33]+14.696938456699069*coords[0]*fdo[9]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[34]+14.696938456699069*coords[0]*fdo[10]+8.485281374238571*fdo[2]); 
  ftar[3] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[35]+14.696938456699069*coords[0]*fdo[11]+8.485281374238571*fdo[3]); 
  ftar[4] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[36]+14.696938456699069*coords[0]*fdo[15]+8.485281374238571*fdo[5]); 
  ftar[5] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[37]+4.898979485566357*coords[0]*fdo[17]+2.8284271247461907*fdo[6]); 
  ftar[6] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[38]+4.898979485566357*coords[0]*fdo[18]+2.8284271247461907*fdo[7]); 
  ftar[7] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[39]+4.898979485566357*coords[0]*fdo[19]+2.8284271247461907*fdo[8]); 
  ftar[8] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[40]+4.898979485566357*coords[0]*fdo[23]+2.8284271247461907*fdo[12]); 
  ftar[9] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[41]+4.898979485566357*coords[0]*fdo[24]+2.8284271247461907*fdo[13]); 
  ftar[10] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[42]+4.898979485566357*coords[0]*fdo[25]+2.8284271247461907*fdo[14]); 
  ftar[11] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[43]+14.696938456699069*coords[0]*fdo[26]+8.485281374238571*fdo[16]); 
  ftar[12] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[44]+14.696938456699069*coords[0]*fdo[28]+8.485281374238571*fdo[20]); 
  ftar[13] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[45]+14.696938456699069*coords[0]*fdo[29]+8.485281374238571*fdo[21]); 
  ftar[14] = 0.08333333333333333*((28.46049894151541*coords0R2-9.48683298050514)*fdo[46]+14.696938456699069*coords[0]*fdo[30]+8.485281374238571*fdo[22]); 
  ftar[15] = 0.25*((9.48683298050514*coords0R2-3.1622776601683795)*fdo[47]+4.898979485566357*coords[0]*fdo[31]+2.8284271247461907*fdo[27]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[33]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[32]+6.0*coords[0]*coords[1]*fdo[9]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[37]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[34]+18.0*coords[0]*coords[1]*fdo[17]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[2]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[38]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[35]+18.0*coords[0]*coords[1]*fdo[18]+10.392304845413262*coords[1]*fdo[11]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[3]); 
  ftar[3] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[40]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[36]+18.0*coords[0]*coords[1]*fdo[23]+10.392304845413262*coords[1]*fdo[15]+10.392304845413262*coords[0]*fdo[12]+6.0*fdo[5]); 
  ftar[4] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[43]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[39]+6.0*coords[0]*coords[1]*fdo[26]+3.4641016151377544*coords[1]*fdo[19]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[8]); 
  ftar[5] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[44]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[41]+6.0*coords[0]*coords[1]*fdo[28]+3.4641016151377544*coords[1]*fdo[24]+3.4641016151377544*coords[0]*fdo[20]+2.0*fdo[13]); 
  ftar[6] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[45]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[42]+6.0*coords[0]*coords[1]*fdo[29]+3.4641016151377544*coords[1]*fdo[25]+3.4641016151377544*coords[0]*fdo[21]+2.0*fdo[14]); 
  ftar[7] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[47]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[46]+18.0*coords[0]*coords[1]*fdo[31]+10.392304845413262*coords[1]*fdo[30]+10.392304845413262*coords[0]*fdo[27]+6.0*fdo[22]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[37]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[34]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[33]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[17]+8.485281374238571*coords[1]*coords[2]*fdo[10]+8.485281374238571*coords[0]*coords[2]*fdo[9]+8.485281374238571*coords[0]*coords[1]*fdo[6]+4.898979485566357*coords[2]*fdo[4]+4.898979485566357*coords[1]*fdo[2]+4.898979485566357*coords[0]*fdo[1]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[43]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[39]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[38]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[35]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[26]+25.455844122715718*coords[1]*coords[2]*fdo[19]+25.455844122715718*coords[0]*coords[2]*fdo[18]+25.455844122715718*coords[0]*coords[1]*fdo[16]+14.696938456699069*coords[2]*fdo[11]+14.696938456699069*coords[1]*fdo[8]+14.696938456699069*coords[0]*fdo[7]+8.485281374238571*fdo[3]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[44]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[41]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[40]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[36]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[28]+25.455844122715718*coords[1]*coords[2]*fdo[24]+25.455844122715718*coords[0]*coords[2]*fdo[23]+25.455844122715718*coords[0]*coords[1]*fdo[20]+14.696938456699069*coords[2]*fdo[15]+14.696938456699069*coords[1]*fdo[13]+14.696938456699069*coords[0]*fdo[12]+8.485281374238571*fdo[5]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[47]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[46]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[45]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[42]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[30]+8.485281374238571*coords[0]*coords[2]*fdo[29]+8.485281374238571*coords[0]*coords[1]*fdo[27]+4.898979485566357*coords[2]*fdo[25]+4.898979485566357*coords[1]*fdo[22]+4.898979485566357*coords[0]*fdo[21]+2.8284271247461907*fdo[14]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords3R2 = pow(coords[3],2);

  ftar[0] = 0.125*((34.856850115866756*coords[0]*coords[1]*coords[2]*coords3R2-11.618950038622252*coords[0]*coords[1]*coords[2])*fdo[43]+(20.12461179749811*coords[1]*coords[2]*coords3R2-6.708203932499369*coords[1]*coords[2])*fdo[39]+(20.12461179749811*coords[0]*coords[2]*coords3R2-6.708203932499369*coords[0]*coords[2])*fdo[38]+(20.12461179749811*coords[0]*coords[1]*coords3R2-6.708203932499369*coords[0]*coords[1])*fdo[37]+(11.618950038622252*coords[2]*coords3R2-3.872983346207417*coords[2])*fdo[35]+(11.618950038622252*coords[1]*coords3R2-3.872983346207417*coords[1])*fdo[34]+(11.618950038622252*coords[0]*coords3R2-3.872983346207417*coords[0])*fdo[33]+(6.708203932499369*coords3R2-2.23606797749979)*fdo[32]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[26]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[19]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[18]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[17]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[16]+6.0*coords[2]*coords[3]*fdo[11]+6.0*coords[1]*coords[3]*fdo[10]+6.0*coords[0]*coords[3]*fdo[9]+6.0*coords[1]*coords[2]*fdo[8]+6.0*coords[0]*coords[2]*fdo[7]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[3]*fdo[4]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.041666666666666664*((104.57055034760022*coords[0]*coords[1]*coords[2]*coords3R2-34.85685011586674*coords[0]*coords[1]*coords[2])*fdo[47]+(60.37383539249431*coords[1]*coords[2]*coords3R2-20.124611797498105*coords[1]*coords[2])*fdo[46]+(60.37383539249431*coords[0]*coords[2]*coords3R2-20.124611797498105*coords[0]*coords[2])*fdo[45]+(60.37383539249431*coords[0]*coords[1]*coords3R2-20.124611797498105*coords[0]*coords[1])*fdo[44]+(34.85685011586674*coords[2]*coords3R2-11.61895003862225*coords[2])*fdo[42]+(34.85685011586674*coords[1]*coords3R2-11.61895003862225*coords[1])*fdo[41]+(34.85685011586674*coords[0]*coords3R2-11.61895003862225*coords[0])*fdo[40]+(20.124611797498105*coords3R2-6.7082039324993685)*fdo[36]+54.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+31.17691453623978*coords[1]*coords[2]*coords[3]*fdo[30]+31.17691453623978*coords[0]*coords[2]*coords[3]*fdo[29]+31.17691453623978*coords[0]*coords[1]*coords[3]*fdo[28]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[27]+18.0*coords[2]*coords[3]*fdo[25]+18.0*coords[1]*coords[3]*fdo[24]+18.0*coords[0]*coords[3]*fdo[23]+18.0*coords[1]*coords[2]*fdo[22]+18.0*coords[0]*coords[2]*fdo[21]+18.0*coords[0]*coords[1]*fdo[20]+10.392304845413262*coords[3]*fdo[15]+10.392304845413262*coords[2]*fdo[14]+10.392304845413262*coords[1]*fdo[13]+10.392304845413262*coords[0]*fdo[12]+6.0*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[38]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[35]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[33]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[18]+8.485281374238571*coords[1]*coords[2]*fdo[11]+8.485281374238571*coords[0]*coords[2]*fdo[9]+8.485281374238571*coords[0]*coords[1]*fdo[7]+4.898979485566357*coords[2]*fdo[4]+4.898979485566357*coords[1]*fdo[3]+4.898979485566357*coords[0]*fdo[1]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[43]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[39]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[37]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[34]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[26]+25.455844122715718*coords[1]*coords[2]*fdo[19]+25.455844122715718*coords[0]*coords[2]*fdo[17]+25.455844122715718*coords[0]*coords[1]*fdo[16]+14.696938456699069*coords[2]*fdo[10]+14.696938456699069*coords[1]*fdo[8]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[2]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[45]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[42]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[40]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[36]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[29]+25.455844122715718*coords[1]*coords[2]*fdo[25]+25.455844122715718*coords[0]*coords[2]*fdo[23]+25.455844122715718*coords[0]*coords[1]*fdo[21]+14.696938456699069*coords[2]*fdo[15]+14.696938456699069*coords[1]*fdo[14]+14.696938456699069*coords[0]*fdo[12]+8.485281374238571*fdo[5]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[47]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[46]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[44]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[41]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[30]+8.485281374238571*coords[0]*coords[2]*fdo[28]+8.485281374238571*coords[0]*coords[1]*fdo[27]+4.898979485566357*coords[2]*fdo[24]+4.898979485566357*coords[1]*fdo[22]+4.898979485566357*coords[0]*fdo[20]+2.8284271247461907*fdo[13]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[34]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[32]+6.0*coords[0]*coords[1]*fdo[10]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[37]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[33]+18.0*coords[0]*coords[1]*fdo[17]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[39]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[35]+18.0*coords[0]*coords[1]*fdo[19]+10.392304845413262*coords[1]*fdo[11]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[3]); 
  ftar[3] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[41]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[36]+18.0*coords[0]*coords[1]*fdo[24]+10.392304845413262*coords[1]*fdo[15]+10.392304845413262*coords[0]*fdo[13]+6.0*fdo[5]); 
  ftar[4] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[43]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[38]+6.0*coords[0]*coords[1]*fdo[26]+3.4641016151377544*coords[1]*fdo[18]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[7]); 
  ftar[5] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[44]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[40]+6.0*coords[0]*coords[1]*fdo[28]+3.4641016151377544*coords[1]*fdo[23]+3.4641016151377544*coords[0]*fdo[20]+2.0*fdo[12]); 
  ftar[6] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[46]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[42]+6.0*coords[0]*coords[1]*fdo[30]+3.4641016151377544*coords[1]*fdo[25]+3.4641016151377544*coords[0]*fdo[22]+2.0*fdo[14]); 
  ftar[7] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[47]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[45]+18.0*coords[0]*coords[1]*fdo[31]+10.392304845413262*coords[1]*fdo[29]+10.392304845413262*coords[0]*fdo[27]+6.0*fdo[21]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[39]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[35]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[34]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[32]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[19]+8.485281374238571*coords[1]*coords[2]*fdo[11]+8.485281374238571*coords[0]*coords[2]*fdo[10]+8.485281374238571*coords[0]*coords[1]*fdo[8]+4.898979485566357*coords[2]*fdo[4]+4.898979485566357*coords[1]*fdo[3]+4.898979485566357*coords[0]*fdo[2]+2.8284271247461907*fdo[0]); 
  ftar[1] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[43]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[38]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[37]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[33]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[26]+25.455844122715718*coords[1]*coords[2]*fdo[18]+25.455844122715718*coords[0]*coords[2]*fdo[17]+25.455844122715718*coords[0]*coords[1]*fdo[16]+14.696938456699069*coords[2]*fdo[9]+14.696938456699069*coords[1]*fdo[7]+14.696938456699069*coords[0]*fdo[6]+8.485281374238571*fdo[1]); 
  ftar[2] = 0.041666666666666664*((85.38149682454623*coords[0]*coords[1]*coords2R2-28.46049894151541*coords[0]*coords[1])*fdo[46]+(49.295030175464944*coords[1]*coords2R2-16.431676725154983*coords[1])*fdo[42]+(49.295030175464944*coords[0]*coords2R2-16.431676725154983*coords[0])*fdo[41]+(28.46049894151541*coords2R2-9.48683298050514)*fdo[36]+44.090815370097204*coords[0]*coords[1]*coords[2]*fdo[30]+25.455844122715718*coords[1]*coords[2]*fdo[25]+25.455844122715718*coords[0]*coords[2]*fdo[24]+25.455844122715718*coords[0]*coords[1]*fdo[22]+14.696938456699069*coords[2]*fdo[15]+14.696938456699069*coords[1]*fdo[14]+14.696938456699069*coords[0]*fdo[13]+8.485281374238571*fdo[5]); 
  ftar[3] = 0.125*((28.460498941515418*coords[0]*coords[1]*coords2R2-9.48683298050514*coords[0]*coords[1])*fdo[47]+(16.431676725154986*coords[1]*coords2R2-5.477225575051662*coords[1])*fdo[45]+(16.431676725154986*coords[0]*coords2R2-5.477225575051662*coords[0])*fdo[44]+(9.48683298050514*coords2R2-3.1622776601683795)*fdo[40]+14.696938456699069*coords[0]*coords[1]*coords[2]*fdo[31]+8.485281374238571*coords[1]*coords[2]*fdo[29]+8.485281374238571*coords[0]*coords[2]*fdo[28]+8.485281374238571*coords[0]*coords[1]*fdo[27]+4.898979485566357*coords[2]*fdo[23]+4.898979485566357*coords[1]*fdo[21]+4.898979485566357*coords[0]*fdo[20]+2.8284271247461907*fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[35]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[32]+6.0*coords[0]*coords[1]*fdo[11]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[3]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[38]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[33]+18.0*coords[0]*coords[1]*fdo[18]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[39]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[34]+18.0*coords[0]*coords[1]*fdo[19]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[2]); 
  ftar[3] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[42]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[36]+18.0*coords[0]*coords[1]*fdo[25]+10.392304845413262*coords[1]*fdo[15]+10.392304845413262*coords[0]*fdo[14]+6.0*fdo[5]); 
  ftar[4] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[43]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[37]+6.0*coords[0]*coords[1]*fdo[26]+3.4641016151377544*coords[1]*fdo[17]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[6]); 
  ftar[5] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[45]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[40]+6.0*coords[0]*coords[1]*fdo[29]+3.4641016151377544*coords[1]*fdo[23]+3.4641016151377544*coords[0]*fdo[21]+2.0*fdo[12]); 
  ftar[6] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[46]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[41]+6.0*coords[0]*coords[1]*fdo[30]+3.4641016151377544*coords[1]*fdo[24]+3.4641016151377544*coords[0]*fdo[22]+2.0*fdo[13]); 
  ftar[7] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[47]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[44]+18.0*coords[0]*coords[1]*fdo[31]+10.392304845413262*coords[1]*fdo[28]+10.392304845413262*coords[0]*fdo[27]+6.0*fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[4] = 0.5*(2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[5]); 
  ftar[5] = 0.5*(2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[8]); 
  ftar[6] = 0.5*(2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[10]); 
  ftar[7] = 0.5*(2.4494897427831783*coords[0]*fdo[18]+1.4142135623730951*fdo[11]); 
  ftar[8] = 0.5*(2.4494897427831783*coords[0]*fdo[20]+1.4142135623730951*fdo[13]); 
  ftar[9] = 0.5*(2.4494897427831783*coords[0]*fdo[21]+1.4142135623730951*fdo[14]); 
  ftar[10] = 0.5*(2.4494897427831783*coords[0]*fdo[23]+1.4142135623730951*fdo[15]); 
  ftar[11] = 0.5*(2.4494897427831783*coords[0]*fdo[26]+1.4142135623730951*fdo[19]); 
  ftar[12] = 0.5*(2.4494897427831783*coords[0]*fdo[27]+1.4142135623730951*fdo[22]); 
  ftar[13] = 0.5*(2.4494897427831783*coords[0]*fdo[28]+1.4142135623730951*fdo[24]); 
  ftar[14] = 0.5*(2.4494897427831783*coords[0]*fdo[29]+1.4142135623730951*fdo[25]); 
  ftar[15] = 0.5*(2.4494897427831783*coords[0]*fdo[31]+1.4142135623730951*fdo[30]); 
  ftar[16] = 0.1*(12.247448713915892*coords[0]*fdo[33]+7.0710678118654755*fdo[32]); 
  ftar[17] = 0.1*(12.247448713915892*coords[0]*fdo[37]+7.0710678118654755*fdo[34]); 
  ftar[18] = 0.1*(12.247448713915892*coords[0]*fdo[38]+7.0710678118654755*fdo[35]); 
  ftar[19] = 0.1*(12.247448713915892*coords[0]*fdo[40]+7.0710678118654755*fdo[36]); 
  ftar[20] = 0.1*(12.247448713915892*coords[0]*fdo[43]+7.0710678118654755*fdo[39]); 
  ftar[21] = 0.1*(12.247448713915892*coords[0]*fdo[44]+7.0710678118654755*fdo[41]); 
  ftar[22] = 0.1*(12.247448713915892*coords[0]*fdo[45]+7.0710678118654755*fdo[42]); 
  ftar[23] = 0.1*(12.247448713915892*coords[0]*fdo[47]+7.0710678118654755*fdo[46]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[20]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[12]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[18]+fdo[11]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[22]+1.7320508075688772*coords[0]*fdo[21]+fdo[14]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[23]+fdo[15]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[29]+fdo[25]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[37]+8.660254037844387*coords[1]*fdo[34]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[43]+8.660254037844387*coords[1]*fdo[39]+8.660254037844387*coords[0]*fdo[38]+5.0*fdo[35]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[41]+8.660254037844387*coords[0]*fdo[40]+5.0*fdo[36]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[46]+8.660254037844387*coords[0]*fdo[45]+5.0*fdo[42]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[16]+4.242640687119286*coords[1]*coords[2]*fdo[8]+4.242640687119286*coords[0]*coords[2]*fdo[7]+4.242640687119286*coords[0]*coords[1]*fdo[6]+2.4494897427831783*coords[2]*fdo[3]+2.4494897427831783*coords[1]*fdo[2]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[26]+4.242640687119286*coords[1]*coords[2]*fdo[19]+4.242640687119286*coords[0]*coords[2]*fdo[18]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[11]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[20]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[13]+2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[5]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[1]*fdo[28]+2.4494897427831783*coords[2]*fdo[25]+2.4494897427831783*coords[1]*fdo[24]+2.4494897427831783*coords[0]*fdo[23]+1.4142135623730951*fdo[15]); 
  ftar[4] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[43]+21.213203435596427*coords[1]*coords[2]*fdo[39]+21.213203435596427*coords[0]*coords[2]*fdo[38]+21.213203435596427*coords[0]*coords[1]*fdo[37]+12.247448713915892*coords[2]*fdo[35]+12.247448713915892*coords[1]*fdo[34]+12.247448713915892*coords[0]*fdo[33]+7.0710678118654755*fdo[32]); 
  ftar[5] = 0.05*(36.74234614174768*coords[0]*coords[1]*coords[2]*fdo[47]+21.213203435596427*coords[1]*coords[2]*fdo[46]+21.213203435596427*coords[0]*coords[2]*fdo[45]+21.213203435596427*coords[0]*coords[1]*fdo[44]+12.247448713915892*coords[2]*fdo[42]+12.247448713915892*coords[1]*fdo[41]+12.247448713915892*coords[0]*fdo[40]+7.0710678118654755*fdo[36]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[6]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[18]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[21]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[12]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[17]+fdo[10]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[22]+1.7320508075688772*coords[0]*fdo[20]+fdo[13]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[23]+fdo[15]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[28]+fdo[24]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[38]+8.660254037844387*coords[1]*fdo[35]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[43]+8.660254037844387*coords[1]*fdo[39]+8.660254037844387*coords[0]*fdo[37]+5.0*fdo[34]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[42]+8.660254037844387*coords[0]*fdo[40]+5.0*fdo[36]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[46]+8.660254037844387*coords[0]*fdo[44]+5.0*fdo[41]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[8]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[10]+1.4142135623730951*fdo[4]); 
  ftar[4] = 0.5*(2.4494897427831783*coords[0]*fdo[13]+1.4142135623730951*fdo[5]); 
  ftar[5] = 0.5*(2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[7]); 
  ftar[6] = 0.5*(2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[9]); 
  ftar[7] = 0.5*(2.4494897427831783*coords[0]*fdo[19]+1.4142135623730951*fdo[11]); 
  ftar[8] = 0.5*(2.4494897427831783*coords[0]*fdo[20]+1.4142135623730951*fdo[12]); 
  ftar[9] = 0.5*(2.4494897427831783*coords[0]*fdo[22]+1.4142135623730951*fdo[14]); 
  ftar[10] = 0.5*(2.4494897427831783*coords[0]*fdo[24]+1.4142135623730951*fdo[15]); 
  ftar[11] = 0.5*(2.4494897427831783*coords[0]*fdo[26]+1.4142135623730951*fdo[18]); 
  ftar[12] = 0.5*(2.4494897427831783*coords[0]*fdo[27]+1.4142135623730951*fdo[21]); 
  ftar[13] = 0.5*(2.4494897427831783*coords[0]*fdo[28]+1.4142135623730951*fdo[23]); 
  ftar[14] = 0.5*(2.4494897427831783*coords[0]*fdo[30]+1.4142135623730951*fdo[25]); 
  ftar[15] = 0.5*(2.4494897427831783*coords[0]*fdo[31]+1.4142135623730951*fdo[29]); 
  ftar[16] = 0.1*(12.247448713915892*coords[0]*fdo[34]+7.0710678118654755*fdo[32]); 
  ftar[17] = 0.1*(12.247448713915892*coords[0]*fdo[37]+7.0710678118654755*fdo[33]); 
  ftar[18] = 0.1*(12.247448713915892*coords[0]*fdo[39]+7.0710678118654755*fdo[35]); 
  ftar[19] = 0.1*(12.247448713915892*coords[0]*fdo[41]+7.0710678118654755*fdo[36]); 
  ftar[20] = 0.1*(12.247448713915892*coords[0]*fdo[43]+7.0710678118654755*fdo[38]); 
  ftar[21] = 0.1*(12.247448713915892*coords[0]*fdo[44]+7.0710678118654755*fdo[40]); 
  ftar[22] = 0.1*(12.247448713915892*coords[0]*fdo[46]+7.0710678118654755*fdo[42]); 
  ftar[23] = 0.1*(12.247448713915892*coords[0]*fdo[47]+7.0710678118654755*fdo[45]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[19]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[9]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[20]+fdo[12]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[24]+fdo[15]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[28]+fdo[23]); 
  ftar[8] = 0.1*(15.0*coords[0]*coords[1]*fdo[39]+8.660254037844387*coords[1]*fdo[35]+8.660254037844387*coords[0]*fdo[34]+5.0*fdo[32]); 
  ftar[9] = 0.1*(15.0*coords[0]*coords[1]*fdo[43]+8.660254037844387*coords[1]*fdo[38]+8.660254037844387*coords[0]*fdo[37]+5.0*fdo[33]); 
  ftar[10] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[42]+8.660254037844387*coords[0]*fdo[41]+5.0*fdo[36]); 
  ftar[11] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[45]+8.660254037844387*coords[0]*fdo[44]+5.0*fdo[40]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[3]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[8]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[11]+1.4142135623730951*fdo[4]); 
  ftar[4] = 0.5*(2.4494897427831783*coords[0]*fdo[14]+1.4142135623730951*fdo[5]); 
  ftar[5] = 0.5*(2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[6]); 
  ftar[6] = 0.5*(2.4494897427831783*coords[0]*fdo[18]+1.4142135623730951*fdo[9]); 
  ftar[7] = 0.5*(2.4494897427831783*coords[0]*fdo[19]+1.4142135623730951*fdo[10]); 
  ftar[8] = 0.5*(2.4494897427831783*coords[0]*fdo[21]+1.4142135623730951*fdo[12]); 
  ftar[9] = 0.5*(2.4494897427831783*coords[0]*fdo[22]+1.4142135623730951*fdo[13]); 
  ftar[10] = 0.5*(2.4494897427831783*coords[0]*fdo[25]+1.4142135623730951*fdo[15]); 
  ftar[11] = 0.5*(2.4494897427831783*coords[0]*fdo[26]+1.4142135623730951*fdo[17]); 
  ftar[12] = 0.5*(2.4494897427831783*coords[0]*fdo[27]+1.4142135623730951*fdo[20]); 
  ftar[13] = 0.5*(2.4494897427831783*coords[0]*fdo[29]+1.4142135623730951*fdo[23]); 
  ftar[14] = 0.5*(2.4494897427831783*coords[0]*fdo[30]+1.4142135623730951*fdo[24]); 
  ftar[15] = 0.5*(2.4494897427831783*coords[0]*fdo[31]+1.4142135623730951*fdo[28]); 
  ftar[16] = 0.1*(12.247448713915892*coords[0]*fdo[35]+7.0710678118654755*fdo[32]); 
  ftar[17] = 0.1*(12.247448713915892*coords[0]*fdo[38]+7.0710678118654755*fdo[33]); 
  ftar[18] = 0.1*(12.247448713915892*coords[0]*fdo[39]+7.0710678118654755*fdo[34]); 
  ftar[19] = 0.1*(12.247448713915892*coords[0]*fdo[42]+7.0710678118654755*fdo[36]); 
  ftar[20] = 0.1*(12.247448713915892*coords[0]*fdo[43]+7.0710678118654755*fdo[37]); 
  ftar[21] = 0.1*(12.247448713915892*coords[0]*fdo[45]+7.0710678118654755*fdo[40]); 
  ftar[22] = 0.1*(12.247448713915892*coords[0]*fdo[46]+7.0710678118654755*fdo[41]); 
  ftar[23] = 0.1*(12.247448713915892*coords[0]*fdo[47]+7.0710678118654755*fdo[44]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_4_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 3;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 24;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_34_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 3;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_034_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_04_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_014_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_024_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_14_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_24_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 3;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 24;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_GKHYBRID_VEL;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 24;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_3x2v_gkhyb_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 24;
}

