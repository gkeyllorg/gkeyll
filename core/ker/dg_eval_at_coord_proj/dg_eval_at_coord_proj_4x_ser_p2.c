#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[12]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[20]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
  ftar[2] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[22]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[3]); 
  ftar[3] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[26]+10.392304845413262*coords[0]*fdo[9]+6.0*fdo[4]); 
  ftar[4] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[33]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[6]); 
  ftar[5] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[36]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[8]); 
  ftar[6] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[38]+3.4641016151377544*coords[0]*fdo[18]+2.0*fdo[10]); 
  ftar[7] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[19]+5.0*fdo[11]); 
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[24]+5.0*fdo[13]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[29]+5.0*fdo[14]); 
  ftar[10] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[45]+10.392304845413262*coords[0]*fdo[31]+6.0*fdo[17]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[32]+5.0*fdo[21]); 
  ftar[12] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[34]+5.0*fdo[23]); 
  ftar[13] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[35]+5.0*fdo[25]); 
  ftar[14] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[40]+5.0*fdo[27]); 
  ftar[15] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[41]+5.0*fdo[28]); 
  ftar[16] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[43]+5.0*fdo[30]); 
  ftar[17] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[44]+5.0*fdo[37]); 
  ftar[18] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[46]+5.0*fdo[39]); 
  ftar[19] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[47]+5.0*fdo[42]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[24]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[22]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[12]+6.0*coords[0]*coords[1]*fdo[7]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[34]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[33]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[23]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[20]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[1]*fdo[6]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[40]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[38]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[27]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[26]+18.0*coords[0]*coords[1]*fdo[18]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[9]+6.0*fdo[4]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[46]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[45]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[39]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[36]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[17]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[8]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[32]+8.660254037844387*coords[1]*fdo[21]+8.660254037844387*coords[0]*fdo[19]+5.0*fdo[11]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[43]+8.660254037844387*coords[1]*fdo[30]+8.660254037844387*coords[0]*fdo[29]+5.0*fdo[14]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[37]+8.660254037844387*coords[0]*fdo[35]+5.0*fdo[25]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[42]+8.660254037844387*coords[0]*fdo[41]+5.0*fdo[28]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[43]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[40]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[2]*fdo[38]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[30]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[29]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[27]+(11.618950038622252*coords0R2-3.872983346207417)*coords[2]*fdo[26]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[24]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[22]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[18]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[12]+6.0*coords[1]*coords[2]*fdo[10]+6.0*coords[0]*coords[2]*fdo[9]+6.0*coords[0]*coords[1]*fdo[7]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[47]+(60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[46]+(60.37383539249431*coords0R2-20.124611797498105)*coords[1]*coords[2]*fdo[45]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[42]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[41]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[39]+(34.85685011586674*coords0R2-11.61895003862225)*coords[2]*fdo[36]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[34]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[33]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[31]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[28]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[23]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[20]+18.0*coords[1]*coords[2]*fdo[17]+18.0*coords[0]*coords[2]*fdo[16]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[2]*fdo[8]+10.392304845413262*coords[1]*fdo[6]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
  ftar[2] = 0.07071067811865474*(25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[44]+15.0*coords[1]*coords[2]*fdo[37]+15.0*coords[0]*coords[2]*fdo[35]+15.0*coords[0]*coords[1]*fdo[32]+8.660254037844387*coords[2]*fdo[25]+8.660254037844387*coords[1]*fdo[21]+8.660254037844387*coords[0]*fdo[19]+5.0*fdo[11]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);
  const double coords3R2 = pow(coords[3],2);

  ftar[0] = 0.0883883476483184*((34.856850115866756*coords[0]*coords[1]*coords[2]*coords3R2-11.618950038622252*coords[0]*coords[1]*coords[2])*fdo[47]+(34.856850115866756*coords[0]*coords[1]*coords2R2-11.618950038622252*coords[0]*coords[1])*coords[3]*fdo[46]+(34.856850115866756*coords[0]*coords1R2-11.618950038622252*coords[0])*coords[2]*coords[3]*fdo[45]+(34.856850115866756*coords0R2-11.618950038622252)*coords[1]*coords[2]*coords[3]*fdo[44]+(20.12461179749811*coords[1]*coords[2]*coords3R2-6.708203932499369*coords[1]*coords[2])*fdo[43]+(20.12461179749811*coords[0]*coords[2]*coords3R2-6.708203932499369*coords[0]*coords[2])*fdo[42]+(20.12461179749811*coords[0]*coords[1]*coords3R2-6.708203932499369*coords[0]*coords[1])*fdo[41]+(20.12461179749811*coords[1]*coords2R2-6.708203932499369*coords[1])*coords[3]*fdo[40]+(20.12461179749811*coords[0]*coords2R2-6.708203932499369*coords[0])*coords[3]*fdo[39]+(20.12461179749811*coords1R2-6.708203932499369)*coords[2]*coords[3]*fdo[38]+(20.12461179749811*coords0R2-6.708203932499369)*coords[2]*coords[3]*fdo[37]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[3]*fdo[36]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[3]*fdo[35]+(20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[34]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[33]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[2]*fdo[32]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+(11.618950038622252*coords[2]*coords3R2-3.872983346207417*coords[2])*fdo[30]+(11.618950038622252*coords[1]*coords3R2-3.872983346207417*coords[1])*fdo[29]+(11.618950038622252*coords[0]*coords3R2-3.872983346207417*coords[0])*fdo[28]+(11.618950038622252*coords2R2-3.872983346207417)*coords[3]*fdo[27]+(11.618950038622252*coords1R2-3.872983346207417)*coords[3]*fdo[26]+(11.618950038622252*coords0R2-3.872983346207417)*coords[3]*fdo[25]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[24]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[23]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[22]+(11.618950038622252*coords0R2-3.872983346207417)*coords[2]*fdo[21]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[20]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[19]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[18]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[17]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[16]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[15]+(6.708203932499369*coords3R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[12]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[2]*coords[3]*fdo[10]+6.0*coords[1]*coords[3]*fdo[9]+6.0*coords[0]*coords[3]*fdo[8]+6.0*coords[1]*coords[2]*fdo[7]+6.0*coords[0]*coords[2]*fdo[6]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[3]*fdo[4]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[34]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[33]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[2]*fdo[32]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[24]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[23]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[22]+(11.618950038622252*coords0R2-3.872983346207417)*coords[2]*fdo[21]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[20]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[19]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[15]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[12]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[1]*coords[2]*fdo[7]+6.0*coords[0]*coords[2]*fdo[6]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[46]+(60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[45]+(60.37383539249431*coords0R2-20.124611797498105)*coords[1]*coords[2]*fdo[44]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[40]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[39]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[38]+(34.85685011586674*coords0R2-11.61895003862225)*coords[2]*fdo[37]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[36]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[35]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[31]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[27]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[26]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[25]+18.0*coords[1]*coords[2]*fdo[18]+18.0*coords[0]*coords[2]*fdo[17]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[2]*fdo[10]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
  ftar[2] = 0.07071067811865474*(25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[47]+15.0*coords[1]*coords[2]*fdo[43]+15.0*coords[0]*coords[2]*fdo[42]+15.0*coords[0]*coords[1]*fdo[41]+8.660254037844387*coords[2]*fdo[30]+8.660254037844387*coords[1]*fdo[29]+8.660254037844387*coords[0]*fdo[28]+5.0*fdo[14]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[29]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[26]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[12]+6.0*coords[0]*coords[1]*fdo[9]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[41]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[36]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[28]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[20]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[1]*fdo[8]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[43]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[38]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[30]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[22]+18.0*coords[0]*coords[1]*fdo[18]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[3]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[47]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[45]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[42]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[33]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[17]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[6]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[35]+8.660254037844387*coords[1]*fdo[25]+8.660254037844387*coords[0]*fdo[19]+5.0*fdo[11]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[40]+8.660254037844387*coords[1]*fdo[27]+8.660254037844387*coords[0]*fdo[24]+5.0*fdo[13]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[37]+8.660254037844387*coords[0]*fdo[32]+5.0*fdo[21]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[39]+8.660254037844387*coords[0]*fdo[34]+5.0*fdo[23]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[41]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[36]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[2]*fdo[35]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[29]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[28]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[26]+(11.618950038622252*coords0R2-3.872983346207417)*coords[2]*fdo[25]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[20]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[19]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[16]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[12]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[1]*coords[2]*fdo[9]+6.0*coords[0]*coords[2]*fdo[8]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[47]+(60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[45]+(60.37383539249431*coords0R2-20.124611797498105)*coords[1]*coords[2]*fdo[44]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[43]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[42]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[38]+(34.85685011586674*coords0R2-11.61895003862225)*coords[2]*fdo[37]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[33]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[32]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[31]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[30]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[22]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[21]+18.0*coords[1]*coords[2]*fdo[18]+18.0*coords[0]*coords[2]*fdo[17]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[2]*fdo[10]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[3]); 
  ftar[2] = 0.07071067811865474*(25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[46]+15.0*coords[1]*coords[2]*fdo[40]+15.0*coords[0]*coords[2]*fdo[39]+15.0*coords[0]*coords[1]*fdo[34]+8.660254037844387*coords[2]*fdo[27]+8.660254037844387*coords[1]*fdo[24]+8.660254037844387*coords[0]*fdo[23]+5.0*fdo[13]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[20]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[19]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[12]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[33]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[32]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[22]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[21]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[3]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[36]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[35]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[26]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[25]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[45]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[44]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[38]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[37]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[18]+3.4641016151377544*coords[0]*fdo[17]+2.0*fdo[10]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[34]+8.660254037844387*coords[1]*fdo[24]+8.660254037844387*coords[0]*fdo[23]+5.0*fdo[13]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[41]+8.660254037844387*coords[1]*fdo[29]+8.660254037844387*coords[0]*fdo[28]+5.0*fdo[14]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[40]+8.660254037844387*coords[0]*fdo[39]+5.0*fdo[27]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[43]+8.660254037844387*coords[0]*fdo[42]+5.0*fdo[30]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[13]+3.4641016151377544*coords[0]*fdo[3]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[23]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
  ftar[2] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[24]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[2]); 
  ftar[3] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[27]+10.392304845413262*coords[0]*fdo[10]+6.0*fdo[4]); 
  ftar[4] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[34]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[5]); 
  ftar[5] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[39]+3.4641016151377544*coords[0]*fdo[17]+2.0*fdo[8]); 
  ftar[6] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[40]+3.4641016151377544*coords[0]*fdo[18]+2.0*fdo[9]); 
  ftar[7] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[21]+5.0*fdo[11]); 
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[22]+5.0*fdo[12]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[30]+5.0*fdo[14]); 
  ftar[10] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[46]+10.392304845413262*coords[0]*fdo[31]+6.0*fdo[16]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[32]+5.0*fdo[19]); 
  ftar[12] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[33]+5.0*fdo[20]); 
  ftar[13] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[37]+5.0*fdo[25]); 
  ftar[14] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[38]+5.0*fdo[26]); 
  ftar[15] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[42]+5.0*fdo[28]); 
  ftar[16] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[43]+5.0*fdo[29]); 
  ftar[17] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[44]+5.0*fdo[35]); 
  ftar[18] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[45]+5.0*fdo[36]); 
  ftar[19] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[47]+5.0*fdo[41]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[30]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[27]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[13]+6.0*coords[0]*coords[1]*fdo[10]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[3]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[42]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[39]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[28]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[23]+18.0*coords[0]*coords[1]*fdo[17]+10.392304845413262*coords[1]*fdo[8]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[43]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[40]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[29]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[24]+18.0*coords[0]*coords[1]*fdo[18]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[2]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[47]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[46]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[41]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[34]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[16]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[5]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[37]+8.660254037844387*coords[1]*fdo[25]+8.660254037844387*coords[0]*fdo[21]+5.0*fdo[11]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[38]+8.660254037844387*coords[1]*fdo[26]+8.660254037844387*coords[0]*fdo[22]+5.0*fdo[12]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[44]+8.660254037844387*coords[1]*fdo[35]+8.660254037844387*coords[0]*fdo[32]+5.0*fdo[19]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[36]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);
  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[42]+(20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[39]+(20.12461179749811*coords0R2-6.708203932499369)*coords[1]*coords[2]*fdo[37]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[30]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[28]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[27]+(11.618950038622252*coords0R2-3.872983346207417)*coords[2]*fdo[25]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[23]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[21]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[17]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[1]*coords[2]*fdo[10]+6.0*coords[0]*coords[2]*fdo[8]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[47]+(60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[46]+(60.37383539249431*coords0R2-20.124611797498105)*coords[1]*coords[2]*fdo[44]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[43]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[41]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[40]+(34.85685011586674*coords0R2-11.61895003862225)*coords[2]*fdo[35]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[34]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[32]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[31]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[29]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[24]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[19]+18.0*coords[1]*coords[2]*fdo[18]+18.0*coords[0]*coords[2]*fdo[16]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[2]*fdo[9]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
  ftar[2] = 0.07071067811865474*(25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[45]+15.0*coords[1]*coords[2]*fdo[38]+15.0*coords[0]*coords[2]*fdo[36]+15.0*coords[0]*coords[1]*fdo[33]+8.660254037844387*coords[2]*fdo[26]+8.660254037844387*coords[1]*fdo[22]+8.660254037844387*coords[0]*fdo[20]+5.0*fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[23]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[21]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[13]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[34]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[32]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[24]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[19]+18.0*coords[0]*coords[1]*fdo[15]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[39]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[37]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[27]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[25]+18.0*coords[0]*coords[1]*fdo[17]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[46]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[44]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[40]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[35]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[18]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[9]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[33]+8.660254037844387*coords[1]*fdo[22]+8.660254037844387*coords[0]*fdo[20]+5.0*fdo[12]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[42]+8.660254037844387*coords[1]*fdo[30]+8.660254037844387*coords[0]*fdo[28]+5.0*fdo[14]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[38]+8.660254037844387*coords[0]*fdo[36]+5.0*fdo[26]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[47]+8.660254037844387*coords[1]*fdo[43]+8.660254037844387*coords[0]*fdo[41]+5.0*fdo[29]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[14]+3.4641016151377544*coords[0]*fdo[4]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[28]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[1]); 
  ftar[2] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[29]+10.392304845413262*coords[0]*fdo[9]+6.0*fdo[2]); 
  ftar[3] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[30]+10.392304845413262*coords[0]*fdo[10]+6.0*fdo[3]); 
  ftar[4] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[41]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[5]); 
  ftar[5] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[42]+3.4641016151377544*coords[0]*fdo[17]+2.0*fdo[6]); 
  ftar[6] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[43]+3.4641016151377544*coords[0]*fdo[18]+2.0*fdo[7]); 
  ftar[7] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[25]+5.0*fdo[11]); 
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[26]+5.0*fdo[12]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[27]+5.0*fdo[13]); 
  ftar[10] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[47]+10.392304845413262*coords[0]*fdo[31]+6.0*fdo[15]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[35]+5.0*fdo[19]); 
  ftar[12] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[36]+5.0*fdo[20]); 
  ftar[13] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[37]+5.0*fdo[21]); 
  ftar[14] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[38]+5.0*fdo[22]); 
  ftar[15] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[39]+5.0*fdo[23]); 
  ftar[16] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[40]+5.0*fdo[24]); 
  ftar[17] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[44]+5.0*fdo[32]); 
  ftar[18] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[45]+5.0*fdo[33]); 
  ftar[19] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[46]+5.0*fdo[34]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);
  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[28]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[25]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[14]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+6.0*coords[0]*coords[1]*fdo[8]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[41]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[35]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[29]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[19]+18.0*coords[0]*coords[1]*fdo[16]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[42]+(34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[37]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[30]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[21]+18.0*coords[0]*coords[1]*fdo[17]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[3]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[47]+(11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[44]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[43]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[32]+6.0*coords[0]*coords[1]*fdo[31]+3.4641016151377544*coords[1]*fdo[18]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[7]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[36]+8.660254037844387*coords[1]*fdo[26]+8.660254037844387*coords[0]*fdo[20]+5.0*fdo[12]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[39]+8.660254037844387*coords[1]*fdo[27]+8.660254037844387*coords[0]*fdo[23]+5.0*fdo[13]); 
  ftar[6] = 0.1*(15.0*coords[0]*coords[1]*fdo[45]+8.660254037844387*coords[1]*fdo[38]+8.660254037844387*coords[0]*fdo[33]+5.0*fdo[22]); 
  ftar[7] = 0.1*(15.0*coords[0]*coords[1]*fdo[46]+8.660254037844387*coords[1]*fdo[40]+8.660254037844387*coords[0]*fdo[34]+5.0*fdo[24]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[11]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[19]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
  ftar[2] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[21]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[3]); 
  ftar[3] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[25]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
  ftar[4] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[32]+3.4641016151377544*coords[0]*fdo[15]+2.0*fdo[7]); 
  ftar[5] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[35]+3.4641016151377544*coords[0]*fdo[16]+2.0*fdo[9]); 
  ftar[6] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[37]+3.4641016151377544*coords[0]*fdo[17]+2.0*fdo[10]); 
  ftar[7] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[20]+5.0*fdo[12]); 
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[23]+5.0*fdo[13]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[28]+5.0*fdo[14]); 
  ftar[10] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[44]+10.392304845413262*coords[0]*fdo[31]+6.0*fdo[18]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[33]+5.0*fdo[22]); 
  ftar[12] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[34]+5.0*fdo[24]); 
  ftar[13] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[36]+5.0*fdo[26]); 
  ftar[14] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[39]+5.0*fdo[27]); 
  ftar[15] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[41]+5.0*fdo[29]); 
  ftar[16] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[42]+5.0*fdo[30]); 
  ftar[17] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[45]+5.0*fdo[38]); 
  ftar[18] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[46]+5.0*fdo[40]); 
  ftar[19] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[47]+5.0*fdo[43]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 20;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 20;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 20;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_4x_ser_p2_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 20;
}

