#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[18]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[19]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[17]+1.7320508075688772*coords[0]*fdo[16]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[21]+fdo[12]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[22]+fdo[13]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[27]+fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[1]*coords[2]*fdo[15]+4.242640687119286*coords[0]*coords[2]*fdo[14]+4.242640687119286*coords[0]*coords[1]*fdo[11]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[4]+2.4494897427831783*coords[0]*fdo[3]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[1]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[18]+2.4494897427831783*coords[2]*fdo[12]+2.4494897427831783*coords[1]*fdo[9]+2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[1]*fdo[19]+2.4494897427831783*coords[2]*fdo[13]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[8]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[2]*fdo[27]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[20]+2.4494897427831783*coords[1]*fdo[17]+2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[6]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[21]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[18]+3.0*coords[2]*coords[3]*fdo[15]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[12]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[9]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[24]+3.0*coords[1]*coords[3]*fdo[22]+3.0*coords[0]*coords[3]*fdo[20]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[3]*fdo[13]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[6]+fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.125*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[31]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[30]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[29]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[28]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[27]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[26]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[25]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[24]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[23]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[22]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[21]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[20]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[19]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[18]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[17]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[16]+3.0*coords[3]*coords[4]*fdo[15]+3.0*coords[2]*coords[4]*fdo[14]+3.0*coords[1]*coords[4]*fdo[13]+3.0*coords[0]*coords[4]*fdo[12]+3.0*coords[2]*coords[3]*fdo[11]+3.0*coords[1]*coords[3]*fdo[10]+3.0*coords[0]*coords[3]*fdo[9]+3.0*coords[1]*coords[2]*fdo[8]+3.0*coords[0]*coords[2]*fdo[7]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[4]*fdo[5]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[22]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[19]+3.0*coords[2]*coords[3]*fdo[15]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[23]+3.0*coords[1]*coords[3]*fdo[21]+3.0*coords[0]*coords[3]*fdo[20]+3.0*coords[1]*coords[2]*fdo[18]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[3]*fdo[12]+1.7320508075688772*coords[2]*fdo[9]+1.7320508075688772*coords[1]*fdo[7]+1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[18]+4.242640687119286*coords[1]*coords[2]*fdo[11]+4.242640687119286*coords[0]*coords[2]*fdo[9]+4.242640687119286*coords[0]*coords[1]*fdo[7]+2.4494897427831783*coords[2]*fdo[4]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[26]+4.242640687119286*coords[1]*coords[2]*fdo[19]+4.242640687119286*coords[0]*coords[2]*fdo[17]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[10]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[21]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[14]+2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[5]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[27]+2.4494897427831783*coords[2]*fdo[24]+2.4494897427831783*coords[1]*fdo[22]+2.4494897427831783*coords[0]*fdo[20]+1.4142135623730951*fdo[13]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[26]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[19]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[18]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[17]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[16]+3.0*coords[2]*coords[3]*fdo[11]+3.0*coords[1]*coords[3]*fdo[10]+3.0*coords[0]*coords[3]*fdo[9]+3.0*coords[1]*coords[2]*fdo[8]+3.0*coords[0]*coords[2]*fdo[7]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[27]+3.0*coords[2]*coords[3]*fdo[25]+3.0*coords[1]*coords[3]*fdo[24]+3.0*coords[0]*coords[3]*fdo[23]+3.0*coords[1]*coords[2]*fdo[22]+3.0*coords[0]*coords[2]*fdo[21]+3.0*coords[0]*coords[1]*fdo[20]+1.7320508075688772*coords[3]*fdo[15]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[12]+fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[19]+4.242640687119286*coords[1]*coords[2]*fdo[11]+4.242640687119286*coords[0]*coords[2]*fdo[10]+4.242640687119286*coords[0]*coords[1]*fdo[8]+2.4494897427831783*coords[2]*fdo[4]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[26]+4.242640687119286*coords[1]*coords[2]*fdo[18]+4.242640687119286*coords[0]*coords[2]*fdo[17]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[9]+2.4494897427831783*coords[1]*fdo[7]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[1]*fdo[22]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[14]+2.4494897427831783*coords[0]*fdo[13]+1.4142135623730951*fdo[5]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[27]+2.4494897427831783*coords[2]*fdo[23]+2.4494897427831783*coords[1]*fdo[21]+2.4494897427831783*coords[0]*fdo[20]+1.4142135623730951*fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[21]+4.242640687119286*coords[1]*coords[2]*fdo[14]+4.242640687119286*coords[0]*coords[2]*fdo[12]+4.242640687119286*coords[0]*coords[1]*fdo[7]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[13]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[18]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[24]+2.4494897427831783*coords[1]*fdo[19]+2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[10]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[27]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[22]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[21]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[20]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[16]+3.0*coords[2]*coords[3]*fdo[14]+3.0*coords[1]*coords[3]*fdo[13]+3.0*coords[0]*coords[3]*fdo[12]+3.0*coords[1]*coords[2]*fdo[8]+3.0*coords[0]*coords[2]*fdo[7]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[25]+3.0*coords[1]*coords[3]*fdo[24]+3.0*coords[0]*coords[3]*fdo[23]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[3]*fdo[15]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[9]+fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[1]*coords[2]*fdo[14]+4.242640687119286*coords[0]*coords[2]*fdo[13]+4.242640687119286*coords[0]*coords[1]*fdo[8]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[3]+2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[12]+2.4494897427831783*coords[1]*fdo[7]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[1]*fdo[19]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[10]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[2]*fdo[28]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[23]+2.4494897427831783*coords[1]*fdo[18]+2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[16]+4.242640687119286*coords[1]*coords[2]*fdo[8]+4.242640687119286*coords[0]*coords[2]*fdo[7]+4.242640687119286*coords[0]*coords[1]*fdo[6]+2.4494897427831783*coords[2]*fdo[3]+2.4494897427831783*coords[1]*fdo[2]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[26]+4.242640687119286*coords[1]*coords[2]*fdo[19]+4.242640687119286*coords[0]*coords[2]*fdo[18]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[11]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[20]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[13]+2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[5]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[1]*fdo[28]+2.4494897427831783*coords[2]*fdo[25]+2.4494897427831783*coords[1]*fdo[24]+2.4494897427831783*coords[0]*fdo[23]+1.4142135623730951*fdo[15]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(2.4494897427831783*coords[0]*fdo[4]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.5*(2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.5*(2.4494897427831783*coords[0]*fdo[10]+1.4142135623730951*fdo[2]); 
  ftar[3] = 0.5*(2.4494897427831783*coords[0]*fdo[11]+1.4142135623730951*fdo[3]); 
  ftar[4] = 0.5*(2.4494897427831783*coords[0]*fdo[15]+1.4142135623730951*fdo[5]); 
  ftar[5] = 0.5*(2.4494897427831783*coords[0]*fdo[17]+1.4142135623730951*fdo[6]); 
  ftar[6] = 0.5*(2.4494897427831783*coords[0]*fdo[18]+1.4142135623730951*fdo[7]); 
  ftar[7] = 0.5*(2.4494897427831783*coords[0]*fdo[19]+1.4142135623730951*fdo[8]); 
  ftar[8] = 0.5*(2.4494897427831783*coords[0]*fdo[23]+1.4142135623730951*fdo[12]); 
  ftar[9] = 0.5*(2.4494897427831783*coords[0]*fdo[24]+1.4142135623730951*fdo[13]); 
  ftar[10] = 0.5*(2.4494897427831783*coords[0]*fdo[25]+1.4142135623730951*fdo[14]); 
  ftar[11] = 0.5*(2.4494897427831783*coords[0]*fdo[26]+1.4142135623730951*fdo[16]); 
  ftar[12] = 0.5*(2.4494897427831783*coords[0]*fdo[28]+1.4142135623730951*fdo[20]); 
  ftar[13] = 0.5*(2.4494897427831783*coords[0]*fdo[29]+1.4142135623730951*fdo[21]); 
  ftar[14] = 0.5*(2.4494897427831783*coords[0]*fdo[30]+1.4142135623730951*fdo[22]); 
  ftar[15] = 0.5*(2.4494897427831783*coords[0]*fdo[31]+1.4142135623730951*fdo[27]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[9]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[10]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[11]+fdo[3]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[18]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[22]+1.7320508075688772*coords[0]*fdo[19]+fdo[8]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[26]+fdo[16]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[23]+4.242640687119286*coords[1]*coords[2]*fdo[15]+4.242640687119286*coords[0]*coords[2]*fdo[12]+4.242640687119286*coords[0]*coords[1]*fdo[9]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[4]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[13]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[2]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[18]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[27]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[22]+2.4494897427831783*coords[1]*fdo[19]+2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[28]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[20]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[17]+3.0*coords[2]*coords[3]*fdo[15]+3.0*coords[1]*coords[3]*fdo[13]+3.0*coords[0]*coords[3]*fdo[12]+3.0*coords[1]*coords[2]*fdo[10]+3.0*coords[0]*coords[2]*fdo[9]+3.0*coords[0]*coords[1]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[25]+3.0*coords[1]*coords[3]*fdo[22]+3.0*coords[0]*coords[3]*fdo[21]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[3]*fdo[14]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[1]*coords[2]*fdo[15]+4.242640687119286*coords[0]*coords[2]*fdo[13]+4.242640687119286*coords[0]*coords[1]*fdo[10]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[4]+2.4494897427831783*coords[0]*fdo[2]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[1]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[2]*fdo[20]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[12]+2.4494897427831783*coords[1]*fdo[9]+2.4494897427831783*coords[0]*fdo[6]+1.4142135623730951*fdo[1]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[1]*coords[2]*fdo[25]+4.242640687119286*coords[0]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[1]*fdo[19]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[11]+2.4494897427831783*coords[0]*fdo[8]+1.4142135623730951*fdo[3]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[2]*fdo[27]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[21]+2.4494897427831783*coords[1]*fdo[18]+2.4494897427831783*coords[0]*fdo[16]+1.4142135623730951*fdo[7]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[6]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[18]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[7]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[12]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[16]+fdo[8]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[20]+fdo[13]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[21]+fdo[14]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[27]+fdo[22]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[17]+4.242640687119286*coords[1]*coords[2]*fdo[10]+4.242640687119286*coords[0]*coords[2]*fdo[9]+4.242640687119286*coords[0]*coords[1]*fdo[6]+2.4494897427831783*coords[2]*fdo[4]+2.4494897427831783*coords[1]*fdo[2]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[26]+4.242640687119286*coords[1]*coords[2]*fdo[19]+4.242640687119286*coords[0]*coords[2]*fdo[18]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[11]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[3]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[20]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[13]+2.4494897427831783*coords[0]*fdo[12]+1.4142135623730951*fdo[5]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[1]*fdo[27]+2.4494897427831783*coords[2]*fdo[25]+2.4494897427831783*coords[1]*fdo[22]+2.4494897427831783*coords[0]*fdo[21]+1.4142135623730951*fdo[14]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[6]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[19]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[16]+fdo[7]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[20]+fdo[12]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[22]+fdo[14]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[27]+fdo[21]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[20]+4.242640687119286*coords[1]*coords[2]*fdo[13]+4.242640687119286*coords[0]*coords[2]*fdo[12]+4.242640687119286*coords[0]*coords[1]*fdo[6]+2.4494897427831783*coords[2]*fdo[5]+2.4494897427831783*coords[1]*fdo[2]+2.4494897427831783*coords[0]*fdo[1]+1.4142135623730951*fdo[0]); 
  ftar[1] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[27]+4.242640687119286*coords[1]*coords[2]*fdo[22]+4.242640687119286*coords[0]*coords[2]*fdo[21]+4.242640687119286*coords[0]*coords[1]*fdo[16]+2.4494897427831783*coords[2]*fdo[14]+2.4494897427831783*coords[1]*fdo[8]+2.4494897427831783*coords[0]*fdo[7]+1.4142135623730951*fdo[3]); 
  ftar[2] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[28]+4.242640687119286*coords[1]*coords[2]*fdo[24]+4.242640687119286*coords[0]*coords[2]*fdo[23]+4.242640687119286*coords[0]*coords[1]*fdo[17]+2.4494897427831783*coords[2]*fdo[15]+2.4494897427831783*coords[1]*fdo[10]+2.4494897427831783*coords[0]*fdo[9]+1.4142135623730951*fdo[4]); 
  ftar[3] = 0.25*(7.348469228349534*coords[0]*coords[1]*coords[2]*fdo[31]+4.242640687119286*coords[1]*coords[2]*fdo[30]+4.242640687119286*coords[0]*coords[2]*fdo[29]+4.242640687119286*coords[0]*coords[1]*fdo[26]+2.4494897427831783*coords[2]*fdo[25]+2.4494897427831783*coords[1]*fdo[19]+2.4494897427831783*coords[0]*fdo[18]+1.4142135623730951*fdo[11]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
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
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_24_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_024_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_34_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_034_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_4_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_04_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_014_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_14_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_5x_ser_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

