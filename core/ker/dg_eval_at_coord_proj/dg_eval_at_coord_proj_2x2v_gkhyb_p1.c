#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
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
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[19]+5.0*fdo[16]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[21]+5.0*fdo[17]); 
  ftar[10] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[22]+5.0*fdo[18]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[23]+5.0*fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[19]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[16]+6.0*coords[0]*coords[1]*fdo[10]+3.4641016151377544*coords[1]*fdo[4]+3.4641016151377544*coords[0]*fdo[3]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[21]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[17]+18.0*coords[0]*coords[1]*fdo[13]+10.392304845413262*coords[1]*fdo[8]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords0R2-11.61895003862225)*coords[1]*fdo[22]+(20.124611797498105*coords0R2-6.7082039324993685)*fdo[18]+18.0*coords[0]*coords[1]*fdo[14]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[2]); 
  ftar[3] = 0.25*((11.618950038622252*coords0R2-3.872983346207417)*coords[1]*fdo[23]+(6.708203932499369*coords0R2-2.23606797749979)*fdo[20]+6.0*coords[0]*coords[1]*fdo[15]+3.4641016151377544*coords[1]*fdo[12]+3.4641016151377544*coords[0]*fdo[11]+2.0*fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[21]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[19]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[17]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[16]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[13]+6.0*coords[1]*coords[2]*fdo[10]+6.0*coords[0]*coords[2]*fdo[8]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[23]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[22]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[20]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[18]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[15]+18.0*coords[1]*coords[2]*fdo[14]+18.0*coords[0]*coords[2]*fdo[12]+18.0*coords[0]*coords[1]*fdo[11]+10.392304845413262*coords[2]*fdo[9]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.0883883476483184*((34.856850115866756*coords[0]*coords[1]*coords2R2-11.618950038622252*coords[0]*coords[1])*coords[3]*fdo[23]+(20.12461179749811*coords[1]*coords2R2-6.708203932499369*coords[1])*coords[3]*fdo[22]+(20.12461179749811*coords[0]*coords2R2-6.708203932499369*coords[0])*coords[3]*fdo[21]+(20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[20]+(11.618950038622252*coords2R2-3.872983346207417)*coords[3]*fdo[19]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[18]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[17]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[16]+18.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[15]+10.392304845413262*coords[1]*coords[2]*coords[3]*fdo[14]+10.392304845413262*coords[0]*coords[2]*coords[3]*fdo[13]+10.392304845413262*coords[0]*coords[1]*coords[3]*fdo[12]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[11]+6.0*coords[2]*coords[3]*fdo[10]+6.0*coords[1]*coords[3]*fdo[9]+6.0*coords[0]*coords[3]*fdo[8]+6.0*coords[1]*coords[2]*fdo[7]+6.0*coords[0]*coords[2]*fdo[6]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[3]*fdo[4]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords1R2-6.708203932499369*coords[0])*coords[2]*fdo[22]+(11.618950038622252*coords1R2-3.872983346207417)*coords[2]*fdo[19]+(11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[18]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[16]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[14]+6.0*coords[1]*coords[2]*fdo[10]+6.0*coords[0]*coords[2]*fdo[9]+6.0*coords[0]*coords[1]*fdo[7]+3.4641016151377544*coords[2]*fdo[4]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords1R2-20.124611797498105*coords[0])*coords[2]*fdo[23]+(34.85685011586674*coords1R2-11.61895003862225)*coords[2]*fdo[21]+(34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[20]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[17]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[15]+18.0*coords[1]*coords[2]*fdo[13]+18.0*coords[0]*coords[2]*fdo[12]+18.0*coords[0]*coords[1]*fdo[11]+10.392304845413262*coords[2]*fdo[8]+10.392304845413262*coords[1]*fdo[6]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[5]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[11]+fdo[7]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[21]+8.660254037844387*coords[1]*fdo[19]+8.660254037844387*coords[0]*fdo[17]+5.0*fdo[16]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[23]+8.660254037844387*coords[1]*fdo[22]+8.660254037844387*coords[0]*fdo[20]+5.0*fdo[18]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[12]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[15]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[2] = 0.07071067811865474*(25.98076211353316*coords[0]*coords[1]*coords[2]*fdo[23]+15.0*coords[1]*coords[2]*fdo[22]+15.0*coords[0]*coords[2]*fdo[21]+15.0*coords[0]*coords[1]*fdo[20]+8.660254037844387*coords[2]*fdo[19]+8.660254037844387*coords[1]*fdo[18]+8.660254037844387*coords[0]*fdo[17]+5.0*fdo[16]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[5]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[11]+fdo[6]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[22]+8.660254037844387*coords[1]*fdo[19]+8.660254037844387*coords[0]*fdo[18]+5.0*fdo[16]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[23]+8.660254037844387*coords[1]*fdo[21]+8.660254037844387*coords[0]*fdo[20]+5.0*fdo[17]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords0R2 = pow(coords[0],2);

  ftar[0] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[16]+3.4641016151377544*coords[0]*fdo[3]+2.0*fdo[0]); 
  ftar[1] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[17]+10.392304845413262*coords[0]*fdo[6]+6.0*fdo[1]); 
  ftar[2] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[18]+10.392304845413262*coords[0]*fdo[7]+6.0*fdo[2]); 
  ftar[3] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[19]+10.392304845413262*coords[0]*fdo[10]+6.0*fdo[4]); 
  ftar[4] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[20]+3.4641016151377544*coords[0]*fdo[11]+2.0*fdo[5]); 
  ftar[5] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[21]+3.4641016151377544*coords[0]*fdo[13]+2.0*fdo[8]); 
  ftar[6] = 0.3535533905932737*((6.708203932499369*coords0R2-2.23606797749979)*fdo[22]+3.4641016151377544*coords[0]*fdo[14]+2.0*fdo[9]); 
  ftar[7] = 0.11785113019775789*((20.124611797498105*coords0R2-6.7082039324993685)*fdo[23]+10.392304845413262*coords[0]*fdo[15]+6.0*fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[17]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[16]+6.0*coords[0]*coords[1]*fdo[6]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[20]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[18]+18.0*coords[0]*coords[1]*fdo[11]+10.392304845413262*coords[1]*fdo[7]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[2]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[21]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[19]+18.0*coords[0]*coords[1]*fdo[13]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[23]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[22]+6.0*coords[0]*coords[1]*fdo[15]+3.4641016151377544*coords[1]*fdo[14]+3.4641016151377544*coords[0]*fdo[12]+2.0*fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords2R2 = pow(coords[2],2);

  ftar[0] = 0.1767766952966368*((20.12461179749811*coords[0]*coords[1]*coords2R2-6.708203932499369*coords[0]*coords[1])*fdo[20]+(11.618950038622252*coords[1]*coords2R2-3.872983346207417*coords[1])*fdo[18]+(11.618950038622252*coords[0]*coords2R2-3.872983346207417*coords[0])*fdo[17]+(6.708203932499369*coords2R2-2.23606797749979)*fdo[16]+10.392304845413262*coords[0]*coords[1]*coords[2]*fdo[11]+6.0*coords[1]*coords[2]*fdo[7]+6.0*coords[0]*coords[2]*fdo[6]+6.0*coords[0]*coords[1]*fdo[5]+3.4641016151377544*coords[2]*fdo[3]+3.4641016151377544*coords[1]*fdo[2]+3.4641016151377544*coords[0]*fdo[1]+2.0*fdo[0]); 
  ftar[1] = 0.05892556509887893*((60.37383539249431*coords[0]*coords[1]*coords2R2-20.124611797498105*coords[0]*coords[1])*fdo[23]+(34.85685011586674*coords[1]*coords2R2-11.61895003862225*coords[1])*fdo[22]+(34.85685011586674*coords[0]*coords2R2-11.61895003862225*coords[0])*fdo[21]+(20.124611797498105*coords2R2-6.7082039324993685)*fdo[19]+31.17691453623978*coords[0]*coords[1]*coords[2]*fdo[15]+18.0*coords[1]*coords[2]*fdo[14]+18.0*coords[0]*coords[2]*fdo[13]+18.0*coords[0]*coords[1]*fdo[12]+10.392304845413262*coords[2]*fdo[10]+10.392304845413262*coords[1]*fdo[9]+10.392304845413262*coords[0]*fdo[8]+6.0*fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.

  const double coords1R2 = pow(coords[1],2);

  ftar[0] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[18]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[16]+6.0*coords[0]*coords[1]*fdo[7]+3.4641016151377544*coords[1]*fdo[3]+3.4641016151377544*coords[0]*fdo[2]+2.0*fdo[0]); 
  ftar[1] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[20]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[17]+18.0*coords[0]*coords[1]*fdo[11]+10.392304845413262*coords[1]*fdo[6]+10.392304845413262*coords[0]*fdo[5]+6.0*fdo[1]); 
  ftar[2] = 0.08333333333333333*((34.85685011586674*coords[0]*coords1R2-11.61895003862225*coords[0])*fdo[22]+(20.124611797498105*coords1R2-6.7082039324993685)*fdo[19]+18.0*coords[0]*coords[1]*fdo[14]+10.392304845413262*coords[1]*fdo[10]+10.392304845413262*coords[0]*fdo[9]+6.0*fdo[4]); 
  ftar[3] = 0.25*((11.618950038622252*coords[0]*coords1R2-3.872983346207417*coords[0])*fdo[23]+(6.708203932499369*coords1R2-2.23606797749979)*fdo[21]+6.0*coords[0]*coords[1]*fdo[15]+3.4641016151377544*coords[1]*fdo[13]+3.4641016151377544*coords[0]*fdo[12]+2.0*fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
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
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[17]+5.0*fdo[16]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[20]+5.0*fdo[18]); 
  ftar[10] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[21]+5.0*fdo[19]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[23]+5.0*fdo[22]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[5]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[3]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[10]); 
  ftar[4] = 0.1*(15.0*coords[0]*coords[1]*fdo[20]+8.660254037844387*coords[1]*fdo[18]+8.660254037844387*coords[0]*fdo[17]+5.0*fdo[16]); 
  ftar[5] = 0.1*(15.0*coords[0]*coords[1]*fdo[23]+8.660254037844387*coords[1]*fdo[22]+8.660254037844387*coords[0]*fdo[21]+5.0*fdo[19]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
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
  ftar[8] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[18]+5.0*fdo[16]); 
  ftar[9] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[20]+5.0*fdo[17]); 
  ftar[10] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[22]+5.0*fdo[19]); 
  ftar[11] = 0.1414213562373095*(8.660254037844387*coords[0]*fdo[23]+5.0*fdo[21]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 2;
  *num_basis = 3;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_HYBRID;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 2;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 0;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_GKHYBRID_VEL;
  *poly_order = 1;
  *num_basis = 6;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_2x2v_gkhyb_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_GKHYBRID;
  *poly_order = 1;
  *num_basis = 12;
}

