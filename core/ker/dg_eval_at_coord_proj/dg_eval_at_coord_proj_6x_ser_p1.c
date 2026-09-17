#include <gkyl_dg_eval_at_coord_proj_kernels.h> 

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[10]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[11]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[12]+fdo[3]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[16]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[20]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[23]+fdo[7]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[24]+fdo[8]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[25]+fdo[9]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[29]+fdo[13]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[30]+fdo[14]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[31]+fdo[15]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[35]+fdo[17]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[36]+fdo[18]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[37]+fdo[19]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[41]+fdo[21]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[42]+fdo[22]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[44]+fdo[26]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[45]+fdo[27]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[46]+fdo[28]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[48]+fdo[32]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[49]+fdo[33]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[50]+fdo[34]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[54]+fdo[38]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[55]+fdo[39]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[56]+fdo[40]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[57]+fdo[43]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[58]+fdo[47]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[60]+fdo[51]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[61]+fdo[52]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[62]+fdo[53]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[59]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[10]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[11]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[12]+fdo[3]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[41]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[20]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[26]+1.7320508075688772*coords[0]*fdo[23]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[24]+fdo[8]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[25]+fdo[9]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[35]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[36]+fdo[18]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[37]+fdo[19]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[43]+1.7320508075688772*coords[0]*fdo[42]+fdo[22]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[51]+1.7320508075688772*coords[0]*fdo[48]+fdo[32]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[49]+fdo[33]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[50]+fdo[34]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[59]+1.7320508075688772*coords[0]*fdo[58]+fdo[47]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[41]+3.0*coords[1]*coords[2]*fdo[21]+3.0*coords[0]*coords[2]*fdo[20]+3.0*coords[0]*coords[1]*fdo[16]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[54]+3.0*coords[1]*coords[2]*fdo[38]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[10]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[55]+3.0*coords[1]*coords[2]*fdo[39]+3.0*coords[0]*coords[2]*fdo[36]+3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[11]+fdo[2]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[56]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[37]+3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[12]+fdo[3]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[51]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[32]+1.7320508075688772*coords[1]*fdo[26]+1.7320508075688772*coords[0]*fdo[23]+fdo[7]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[52]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[33]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[24]+fdo[8]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[50]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[34]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[25]+fdo[9]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[59]+3.0*coords[0]*coords[2]*fdo[58]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[47]+1.7320508075688772*coords[1]*fdo[43]+1.7320508075688772*coords[0]*fdo[42]+fdo[22]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[41]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[38]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[35]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[29]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[20]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[2]*coords[3]*fdo[39]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[3]*fdo[18]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[2]*coords[3]*fdo[40]+3.0*coords[1]*coords[3]*fdo[37]+3.0*coords[0]*coords[3]*fdo[33]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[3]*fdo[19]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[53]+3.0*coords[1]*coords[3]*fdo[50]+3.0*coords[0]*coords[3]*fdo[47]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[34]+1.7320508075688772*coords[2]*fdo[28]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[60]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[55]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[54]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[51]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[48]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[41]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[39]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[38]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[36]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[35]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[32]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[26]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[23]+3.0*coords[3]*coords[4]*fdo[21]+3.0*coords[2]*coords[4]*fdo[20]+3.0*coords[1]*coords[4]*fdo[18]+3.0*coords[0]*coords[4]*fdo[17]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[4]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[59]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[58]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[56]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[53]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[52]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[47]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[3]*coords[4]*fdo[40]+3.0*coords[2]*coords[4]*fdo[37]+3.0*coords[1]*coords[4]*fdo[34]+3.0*coords[0]*coords[4]*fdo[33]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[28]+3.0*coords[0]*coords[3]*fdo[27]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[4]*fdo[19]+1.7320508075688772*coords[3]*fdo[15]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.0883883476483184*(27.0*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*coords[5]*fdo[63]+15.58845726811989*coords[1]*coords[2]*coords[3]*coords[4]*coords[5]*fdo[62]+15.58845726811989*coords[0]*coords[2]*coords[3]*coords[4]*coords[5]*fdo[61]+15.58845726811989*coords[0]*coords[1]*coords[3]*coords[4]*coords[5]*fdo[60]+15.58845726811989*coords[0]*coords[1]*coords[2]*coords[4]*coords[5]*fdo[59]+15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[5]*fdo[58]+15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[57]+9.0*coords[2]*coords[3]*coords[4]*coords[5]*fdo[56]+9.0*coords[1]*coords[3]*coords[4]*coords[5]*fdo[55]+9.0*coords[0]*coords[3]*coords[4]*coords[5]*fdo[54]+9.0*coords[1]*coords[2]*coords[4]*coords[5]*fdo[53]+9.0*coords[0]*coords[2]*coords[4]*coords[5]*fdo[52]+9.0*coords[0]*coords[1]*coords[4]*coords[5]*fdo[51]+9.0*coords[1]*coords[2]*coords[3]*coords[5]*fdo[50]+9.0*coords[0]*coords[2]*coords[3]*coords[5]*fdo[49]+9.0*coords[0]*coords[1]*coords[3]*coords[5]*fdo[48]+9.0*coords[0]*coords[1]*coords[2]*coords[5]*fdo[47]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[46]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[45]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[44]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[43]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[42]+5.196152422706631*coords[3]*coords[4]*coords[5]*fdo[41]+5.196152422706631*coords[2]*coords[4]*coords[5]*fdo[40]+5.196152422706631*coords[1]*coords[4]*coords[5]*fdo[39]+5.196152422706631*coords[0]*coords[4]*coords[5]*fdo[38]+5.196152422706631*coords[2]*coords[3]*coords[5]*fdo[37]+5.196152422706631*coords[1]*coords[3]*coords[5]*fdo[36]+5.196152422706631*coords[0]*coords[3]*coords[5]*fdo[35]+5.196152422706631*coords[1]*coords[2]*coords[5]*fdo[34]+5.196152422706631*coords[0]*coords[2]*coords[5]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[5]*fdo[32]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[31]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[30]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[29]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[28]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[26]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[4]*coords[5]*fdo[21]+3.0*coords[3]*coords[5]*fdo[20]+3.0*coords[2]*coords[5]*fdo[19]+3.0*coords[1]*coords[5]*fdo[18]+3.0*coords[0]*coords[5]*fdo[17]+3.0*coords[3]*coords[4]*fdo[16]+3.0*coords[2]*coords[4]*fdo[15]+3.0*coords[1]*coords[4]*fdo[14]+3.0*coords[0]*coords[4]*fdo[13]+3.0*coords[2]*coords[3]*fdo[12]+3.0*coords[1]*coords[3]*fdo[11]+3.0*coords[0]*coords[3]*fdo[10]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[5]*fdo[6]+1.7320508075688772*coords[4]*fdo[5]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[56]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[54]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[52]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[49]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[41]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[40]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[38]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[37]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[35]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[33]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[24]+3.0*coords[3]*coords[4]*fdo[21]+3.0*coords[2]*coords[4]*fdo[20]+3.0*coords[1]*coords[4]*fdo[19]+3.0*coords[0]*coords[4]*fdo[17]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[15]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[4]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[60]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[59]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[58]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[55]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[53]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[51]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[47]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[3]*coords[4]*fdo[39]+3.0*coords[2]*coords[4]*fdo[36]+3.0*coords[1]*coords[4]*fdo[34]+3.0*coords[0]*coords[4]*fdo[32]+3.0*coords[2]*coords[3]*fdo[30]+3.0*coords[1]*coords[3]*fdo[28]+3.0*coords[0]*coords[3]*fdo[26]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[4]*fdo[18]+1.7320508075688772*coords[3]*fdo[14]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[41]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[39]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[36]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[30]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[20]+3.0*coords[0]*coords[3]*fdo[18]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[14]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[2]*coords[3]*fdo[38]+3.0*coords[1]*coords[3]*fdo[35]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[29]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[3]*fdo[17]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[2]*coords[3]*fdo[40]+3.0*coords[1]*coords[3]*fdo[37]+3.0*coords[0]*coords[3]*fdo[34]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[28]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[3]*fdo[19]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[52]+3.0*coords[1]*coords[3]*fdo[49]+3.0*coords[0]*coords[3]*fdo[47]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[33]+1.7320508075688772*coords[2]*fdo[27]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[56]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[55]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[53]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[50]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[41]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[40]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[39]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[37]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[36]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[34]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[25]+3.0*coords[3]*coords[4]*fdo[21]+3.0*coords[2]*coords[4]*fdo[20]+3.0*coords[1]*coords[4]*fdo[19]+3.0*coords[0]*coords[4]*fdo[18]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[15]+3.0*coords[0]*coords[3]*fdo[14]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[11]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[4]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[60]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[59]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[58]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[54]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[52]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[51]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[49]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[47]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[3]*coords[4]*fdo[38]+3.0*coords[2]*coords[4]*fdo[35]+3.0*coords[1]*coords[4]*fdo[33]+3.0*coords[0]*coords[4]*fdo[32]+3.0*coords[2]*coords[3]*fdo[29]+3.0*coords[1]*coords[3]*fdo[27]+3.0*coords[0]*coords[3]*fdo[26]+3.0*coords[1]*coords[2]*fdo[24]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[4]*fdo[17]+1.7320508075688772*coords[3]*fdo[13]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[41]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[40]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[37]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[31]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[20]+3.0*coords[0]*coords[3]*fdo[19]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[15]+3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[2]*coords[3]*fdo[38]+3.0*coords[1]*coords[3]*fdo[35]+3.0*coords[0]*coords[3]*fdo[33]+3.0*coords[1]*coords[2]*fdo[29]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[3]*fdo[17]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[2]*coords[3]*fdo[39]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[34]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[28]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[3]*fdo[18]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[51]+3.0*coords[1]*coords[3]*fdo[48]+3.0*coords[0]*coords[3]*fdo[47]+3.0*coords[1]*coords[2]*fdo[44]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[32]+1.7320508075688772*coords[2]*fdo[26]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[29]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[54]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[35]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[28]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[53]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[47]+fdo[34]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[26]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[23]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[28]+3.0*coords[0]*coords[3]*fdo[27]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[15]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[48]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[39]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[36]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[53]+3.0*coords[0]*coords[3]*fdo[52]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[3]*fdo[40]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[57]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[46]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[45]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[44]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[43]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[42]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[31]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[30]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[29]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[28]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[26]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[3]*coords[4]*fdo[16]+3.0*coords[2]*coords[4]*fdo[15]+3.0*coords[1]*coords[4]*fdo[14]+3.0*coords[0]*coords[4]*fdo[13]+3.0*coords[2]*coords[3]*fdo[12]+3.0*coords[1]*coords[3]*fdo[11]+3.0*coords[0]*coords[3]*fdo[10]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[4]*fdo[5]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[60]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[59]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[56]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[55]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[54]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[51]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[3]*coords[4]*fdo[41]+3.0*coords[2]*coords[4]*fdo[40]+3.0*coords[1]*coords[4]*fdo[39]+3.0*coords[0]*coords[4]*fdo[38]+3.0*coords[2]*coords[3]*fdo[37]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[35]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[4]*fdo[21]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[29]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[24]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[15]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[30]+3.0*coords[1]*coords[3]*fdo[28]+3.0*coords[0]*coords[3]*fdo[26]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[14]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[49]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[40]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[2]*coords[3]*fdo[55]+3.0*coords[1]*coords[3]*fdo[53]+3.0*coords[0]*coords[3]*fdo[51]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[3]*fdo[39]+1.7320508075688772*coords[2]*fdo[36]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[30]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[14]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[1]*coords[2]*fdo[29]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[28]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[55]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[39]+3.0*coords[0]*coords[1]*fdo[36]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[27]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[53]+3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[34]+fdo[19]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[52]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[47]+fdo[33]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[31]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[30]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[25]+3.0*coords[2]*coords[3]*fdo[16]+3.0*coords[1]*coords[3]*fdo[15]+3.0*coords[0]*coords[3]*fdo[14]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[11]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[43]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[29]+3.0*coords[1]*coords[3]*fdo[27]+3.0*coords[0]*coords[3]*fdo[26]+3.0*coords[1]*coords[2]*fdo[24]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[13]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[50]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[40]+3.0*coords[0]*coords[3]*fdo[39]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[36]+3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[2]*coords[3]*fdo[54]+3.0*coords[1]*coords[3]*fdo[52]+3.0*coords[0]*coords[3]*fdo[51]+3.0*coords[1]*coords[2]*fdo[49]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[3]*fdo[38]+1.7320508075688772*coords[2]*fdo[35]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[31]+3.0*coords[1]*coords[2]*fdo[16]+3.0*coords[0]*coords[2]*fdo[15]+3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[1]*coords[2]*fdo[29]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[28]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[56]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[40]+3.0*coords[0]*coords[1]*fdo[37]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[19]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[44]+3.0*coords[0]*coords[2]*fdo[43]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[26]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[33]+fdo[17]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[53]+3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[34]+fdo[18]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[60]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[51]+1.7320508075688772*coords[1]*fdo[48]+1.7320508075688772*coords[0]*fdo[47]+fdo[32]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[20]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[4]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[35]+1.7320508075688772*coords[1]*fdo[17]+1.7320508075688772*coords[0]*fdo[10]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[36]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[11]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[37]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[12]+fdo[3]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[41]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[16]+fdo[5]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[32]+1.7320508075688772*coords[0]*fdo[23]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[24]+fdo[8]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[25]+fdo[9]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[29]+fdo[13]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[30]+fdo[14]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[31]+fdo[15]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[47]+1.7320508075688772*coords[0]*fdo[42]+fdo[22]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[51]+1.7320508075688772*coords[0]*fdo[44]+fdo[26]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[45]+fdo[27]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[46]+fdo[28]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[59]+1.7320508075688772*coords[0]*fdo[57]+fdo[43]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[35]+3.0*coords[1]*coords[2]*fdo[20]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[48]+3.0*coords[1]*coords[2]*fdo[36]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[49]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[54]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[34]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[53]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[43]+fdo[28]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[48]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[36]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[35]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[32]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[23]+3.0*coords[2]*coords[3]*fdo[20]+3.0*coords[1]*coords[3]*fdo[18]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[37]+3.0*coords[1]*coords[3]*fdo[34]+3.0*coords[0]*coords[3]*fdo[33]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[19]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[39]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[53]+3.0*coords[0]*coords[3]*fdo[52]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[3]*fdo[40]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[58]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[50]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[49]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[48]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[47]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[42]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[37]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[36]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[35]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[34]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[32]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[3]*coords[4]*fdo[20]+3.0*coords[2]*coords[4]*fdo[19]+3.0*coords[1]*coords[4]*fdo[18]+3.0*coords[0]*coords[4]*fdo[17]+3.0*coords[2]*coords[3]*fdo[12]+3.0*coords[1]*coords[3]*fdo[11]+3.0*coords[0]*coords[3]*fdo[10]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[4]*fdo[6]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[60]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[59]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[56]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[55]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[54]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[51]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[3]*coords[4]*fdo[41]+3.0*coords[2]*coords[4]*fdo[40]+3.0*coords[1]*coords[4]*fdo[39]+3.0*coords[0]*coords[4]*fdo[38]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[30]+3.0*coords[0]*coords[3]*fdo[29]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[4]*fdo[21]+1.7320508075688772*coords[3]*fdo[16]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[37]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[35]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[24]+3.0*coords[2]*coords[3]*fdo[20]+3.0*coords[1]*coords[3]*fdo[19]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[36]+3.0*coords[1]*coords[3]*fdo[34]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[18]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[40]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[55]+3.0*coords[1]*coords[3]*fdo[53]+3.0*coords[0]*coords[3]*fdo[51]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[3]*fdo[39]+1.7320508075688772*coords[2]*fdo[30]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[36]+3.0*coords[1]*coords[2]*fdo[20]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[48]+3.0*coords[1]*coords[2]*fdo[35]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[50]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[34]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[55]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[39]+3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[49]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[33]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[53]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[28]+fdo[15]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[52]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[43]+fdo[27]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[37]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[36]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[34]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[25]+3.0*coords[2]*coords[3]*fdo[20]+3.0*coords[1]*coords[3]*fdo[19]+3.0*coords[0]*coords[3]*fdo[18]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[11]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[35]+3.0*coords[1]*coords[3]*fdo[33]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[24]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[17]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[40]+3.0*coords[0]*coords[3]*fdo[39]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[30]+3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[59]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[54]+3.0*coords[1]*coords[3]*fdo[52]+3.0*coords[0]*coords[3]*fdo[51]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[3]*fdo[38]+1.7320508075688772*coords[2]*fdo[29]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[37]+3.0*coords[1]*coords[2]*fdo[20]+3.0*coords[0]*coords[2]*fdo[19]+3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[49]+3.0*coords[1]*coords[2]*fdo[35]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[50]+3.0*coords[1]*coords[2]*fdo[36]+3.0*coords[0]*coords[2]*fdo[34]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[56]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[40]+3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[15]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[48]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[32]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[27]+fdo[13]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[53]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[28]+fdo[14]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[60]+3.0*coords[0]*coords[2]*fdo[59]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[51]+1.7320508075688772*coords[1]*fdo[44]+1.7320508075688772*coords[0]*fdo[43]+fdo[26]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[10]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[35]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[43]+fdo[28]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[47]+fdo[34]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[51]+fdo[39]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[52]+fdo[40]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[62]+1.7320508075688772*coords[0]*fdo[59]+fdo[53]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[23]+3.0*coords[1]*coords[2]*fdo[11]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[48]+3.0*coords[1]*coords[2]*fdo[36]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[61]+3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[2]*fdo[56]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[52]+fdo[40]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[42]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[25]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[24]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[23]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[2]*coords[3]*fdo[12]+3.0*coords[1]*coords[3]*fdo[11]+3.0*coords[0]*coords[3]*fdo[10]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[4]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[30]+3.0*coords[0]*coords[3]*fdo[29]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[3]*fdo[16]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[2]*coords[3]*fdo[37]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[35]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[55]+3.0*coords[0]*coords[3]*fdo[54]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[3]*fdo[41]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[24]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[10]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[11]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[49]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[30]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[36]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[2]*fdo[55]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[51]+fdo[39]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[11]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[36]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[28]+fdo[15]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[34]+fdo[19]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[39]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[43]+fdo[27]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[47]+fdo[33]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[51]+fdo[38]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[53]+fdo[40]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[61]+1.7320508075688772*coords[0]*fdo[59]+fdo[52]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[25]+3.0*coords[1]*coords[2]*fdo[12]+3.0*coords[0]*coords[2]*fdo[11]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[2]*fdo[4]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[1]*coords[2]*fdo[24]+3.0*coords[0]*coords[2]*fdo[23]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[10]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[30]+3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[50]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[36]+3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[29]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[49]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[35]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[55]+3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[39]+fdo[21]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[2]*fdo[54]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[51]+fdo[38]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[12]+1.7320508075688772*coords[1]*fdo[4]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[10]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[15]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[37]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[19]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[23]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[27]+fdo[13]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[28]+fdo[14]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[33]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[34]+fdo[18]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[40]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[44]+1.7320508075688772*coords[0]*fdo[43]+fdo[26]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[48]+1.7320508075688772*coords[0]*fdo[47]+fdo[32]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[52]+fdo[38]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[53]+fdo[39]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[60]+1.7320508075688772*coords[0]*fdo[59]+fdo[51]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[5]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[13]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[14]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[3]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[16]+fdo[4]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[21]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[26]+fdo[7]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[27]+fdo[8]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[28]+fdo[9]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[29]+fdo[10]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[30]+fdo[11]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[31]+fdo[12]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[38]+fdo[17]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[39]+fdo[18]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[40]+fdo[19]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[41]+fdo[20]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[43]+fdo[22]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[44]+fdo[23]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[45]+fdo[24]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[46]+fdo[25]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[51]+fdo[32]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[52]+fdo[33]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[53]+fdo[34]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[54]+fdo[35]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[55]+fdo[36]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[56]+fdo[37]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[57]+fdo[42]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[59]+fdo[47]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[60]+fdo[48]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[61]+fdo[49]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[62]+fdo[50]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[58]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[21]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[5]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[38]+1.7320508075688772*coords[1]*fdo[17]+1.7320508075688772*coords[0]*fdo[13]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[39]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[14]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[40]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[15]+fdo[3]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[41]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[16]+fdo[4]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[32]+1.7320508075688772*coords[0]*fdo[26]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[27]+fdo[8]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[28]+fdo[9]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[29]+fdo[10]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[30]+fdo[11]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[31]+fdo[12]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[47]+1.7320508075688772*coords[0]*fdo[43]+fdo[22]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[48]+1.7320508075688772*coords[0]*fdo[44]+fdo[23]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[45]+fdo[24]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[46]+fdo[25]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[58]+1.7320508075688772*coords[0]*fdo[57]+fdo[42]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[38]+3.0*coords[1]*coords[2]*fdo[21]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[51]+3.0*coords[1]*coords[2]*fdo[39]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[52]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[54]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[34]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[36]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[58]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[50]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[42]+fdo[25]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[51]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[39]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[38]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[32]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[18]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[2]*coords[3]*fdo[40]+3.0*coords[1]*coords[3]*fdo[34]+3.0*coords[0]*coords[3]*fdo[33]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[19]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[35]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[50]+3.0*coords[0]*coords[3]*fdo[49]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[37]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[59]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[53]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[52]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[51]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[47]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[43]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[40]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[39]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[38]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[34]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[32]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[26]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[3]*coords[4]*fdo[21]+3.0*coords[2]*coords[4]*fdo[19]+3.0*coords[1]*coords[4]*fdo[18]+3.0*coords[0]*coords[4]*fdo[17]+3.0*coords[2]*coords[3]*fdo[15]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[4]*fdo[6]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.1767766952966368*(15.58845726811989*coords[0]*coords[1]*coords[2]*coords[3]*coords[4]*fdo[63]+9.0*coords[1]*coords[2]*coords[3]*coords[4]*fdo[62]+9.0*coords[0]*coords[2]*coords[3]*coords[4]*fdo[61]+9.0*coords[0]*coords[1]*coords[3]*coords[4]*fdo[60]+9.0*coords[0]*coords[1]*coords[2]*coords[4]*fdo[58]+9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[2]*coords[3]*coords[4]*fdo[56]+5.196152422706631*coords[1]*coords[3]*coords[4]*fdo[55]+5.196152422706631*coords[0]*coords[3]*coords[4]*fdo[54]+5.196152422706631*coords[1]*coords[2]*coords[4]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[4]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[4]*fdo[48]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[3]*coords[4]*fdo[41]+3.0*coords[2]*coords[4]*fdo[37]+3.0*coords[1]*coords[4]*fdo[36]+3.0*coords[0]*coords[4]*fdo[35]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[30]+3.0*coords[0]*coords[3]*fdo[29]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[4]*fdo[20]+1.7320508075688772*coords[3]*fdo[16]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[40]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[38]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[27]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[19]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[15]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[2]*coords[3]*fdo[39]+3.0*coords[1]*coords[3]*fdo[34]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[18]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[54]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[37]+3.0*coords[0]*coords[3]*fdo[35]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[55]+3.0*coords[1]*coords[3]*fdo[50]+3.0*coords[0]*coords[3]*fdo[48]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[36]+1.7320508075688772*coords[2]*fdo[30]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[39]+3.0*coords[1]*coords[2]*fdo[21]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[51]+3.0*coords[1]*coords[2]*fdo[38]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[53]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[34]+3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[55]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[36]+3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[52]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[33]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[35]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[50]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[25]+fdo[12]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[58]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[49]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[42]+fdo[24]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[40]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[39]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[34]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[28]+3.0*coords[2]*coords[3]*fdo[21]+3.0*coords[1]*coords[3]*fdo[19]+3.0*coords[0]*coords[3]*fdo[18]+3.0*coords[1]*coords[2]*fdo[15]+3.0*coords[0]*coords[2]*fdo[14]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[47]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[2]*coords[3]*fdo[38]+3.0*coords[1]*coords[3]*fdo[33]+3.0*coords[0]*coords[3]*fdo[32]+3.0*coords[1]*coords[2]*fdo[27]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[3]*fdo[17]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[56]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[55]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[2]*coords[3]*fdo[41]+3.0*coords[1]*coords[3]*fdo[37]+3.0*coords[0]*coords[3]*fdo[36]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[30]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[58]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[54]+3.0*coords[1]*coords[3]*fdo[49]+3.0*coords[0]*coords[3]*fdo[48]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[3]*fdo[35]+1.7320508075688772*coords[2]*fdo[29]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[40]+3.0*coords[1]*coords[2]*fdo[21]+3.0*coords[0]*coords[2]*fdo[19]+3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[52]+3.0*coords[1]*coords[2]*fdo[38]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[53]+3.0*coords[1]*coords[2]*fdo[39]+3.0*coords[0]*coords[2]*fdo[34]+3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[56]+3.0*coords[1]*coords[2]*fdo[41]+3.0*coords[0]*coords[2]*fdo[37]+3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[12]+fdo[4]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[51]+3.0*coords[0]*coords[2]*fdo[47]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[32]+1.7320508075688772*coords[1]*fdo[26]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[54]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[35]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[24]+fdo[10]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[50]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[36]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[25]+fdo[11]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[60]+3.0*coords[0]*coords[2]*fdo[58]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[48]+1.7320508075688772*coords[1]*fdo[44]+1.7320508075688772*coords[0]*fdo[42]+fdo[23]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[13]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[29]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[38]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[42]+fdo[25]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[47]+fdo[34]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[48]+fdo[36]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[49]+fdo[37]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[62]+1.7320508075688772*coords[0]*fdo[58]+fdo[50]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[26]+3.0*coords[1]*coords[2]*fdo[14]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[44]+3.0*coords[1]*coords[2]*fdo[30]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[51]+3.0*coords[1]*coords[2]*fdo[39]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[61]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[56]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[49]+fdo[37]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[43]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[28]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[27]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[26]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[2]*coords[3]*fdo[15]+3.0*coords[1]*coords[3]*fdo[14]+3.0*coords[0]*coords[3]*fdo[13]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[5]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[57]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[46]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[45]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[44]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[31]+3.0*coords[1]*coords[3]*fdo[30]+3.0*coords[0]*coords[3]*fdo[29]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[3]*fdo[16]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[2]*coords[3]*fdo[40]+3.0*coords[1]*coords[3]*fdo[39]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[55]+3.0*coords[0]*coords[3]*fdo[54]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[3]*fdo[41]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[27]+3.0*coords[1]*coords[2]*fdo[15]+3.0*coords[0]*coords[2]*fdo[13]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[14]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[45]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[29]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[52]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[30]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[55]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[48]+fdo[36]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[14]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[30]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[39]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[25]+fdo[12]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[34]+fdo[19]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[36]+fdo[20]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[42]+fdo[24]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[47]+fdo[33]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[48]+fdo[35]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[50]+fdo[37]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[61]+1.7320508075688772*coords[0]*fdo[58]+fdo[49]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[28]+3.0*coords[1]*coords[2]*fdo[15]+3.0*coords[0]*coords[2]*fdo[14]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[2]*fdo[5]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[1]*coords[2]*fdo[27]+3.0*coords[0]*coords[2]*fdo[26]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[13]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[46]+3.0*coords[1]*coords[2]*fdo[31]+3.0*coords[0]*coords[2]*fdo[30]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[16]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[53]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[39]+3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[45]+3.0*coords[0]*coords[2]*fdo[44]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[29]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[52]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[55]+3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[36]+fdo[20]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[2]*fdo[54]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[48]+fdo[35]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[15]+1.7320508075688772*coords[1]*fdo[5]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[13]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[31]+1.7320508075688772*coords[1]*fdo[16]+1.7320508075688772*coords[0]*fdo[12]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[40]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[19]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[26]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[29]+1.7320508075688772*coords[0]*fdo[24]+fdo[10]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[25]+fdo[11]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[33]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[34]+fdo[18]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[37]+fdo[20]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[44]+1.7320508075688772*coords[0]*fdo[42]+fdo[23]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[51]+1.7320508075688772*coords[0]*fdo[47]+fdo[32]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[49]+fdo[35]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[50]+fdo[36]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[60]+1.7320508075688772*coords[0]*fdo[58]+fdo[48]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[6]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[17]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[18]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[19]+fdo[3]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[20]+fdo[4]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[21]+fdo[5]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[32]+fdo[7]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[33]+fdo[8]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[34]+fdo[9]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[35]+fdo[10]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[36]+fdo[11]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[37]+fdo[12]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[38]+fdo[13]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[39]+fdo[14]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[40]+fdo[15]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[41]+fdo[16]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[47]+fdo[22]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[48]+fdo[23]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[49]+fdo[24]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[50]+fdo[25]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[51]+fdo[26]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[52]+fdo[27]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[53]+fdo[28]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[54]+fdo[29]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[55]+fdo[30]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[56]+fdo[31]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[58]+fdo[42]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[59]+fdo[43]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[60]+fdo[44]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[61]+fdo[45]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[62]+fdo[46]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[57]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[17]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[35]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[38]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[54]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[42]+fdo[25]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[43]+fdo[28]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[44]+fdo[30]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[45]+fdo[31]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[62]+1.7320508075688772*coords[0]*fdo[57]+fdo[46]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[32]+3.0*coords[1]*coords[2]*fdo[18]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[48]+3.0*coords[1]*coords[2]*fdo[36]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[51]+3.0*coords[1]*coords[2]*fdo[39]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[60]+3.0*coords[1]*coords[2]*fdo[55]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[61]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[56]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[45]+fdo[31]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[47]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[34]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[33]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[32]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[2]*coords[3]*fdo[19]+3.0*coords[1]*coords[3]*fdo[18]+3.0*coords[0]*coords[3]*fdo[17]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[3]*fdo[6]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[58]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[50]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[49]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[48]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[2]*coords[3]*fdo[37]+3.0*coords[1]*coords[3]*fdo[36]+3.0*coords[0]*coords[3]*fdo[35]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[3]*fdo[20]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[2] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[59]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[53]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[52]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[51]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[2]*coords[3]*fdo[40]+3.0*coords[1]*coords[3]*fdo[39]+3.0*coords[0]*coords[3]*fdo[38]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[3]*fdo[21]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.25*(9.0*coords[0]*coords[1]*coords[2]*coords[3]*fdo[63]+5.196152422706631*coords[1]*coords[2]*coords[3]*fdo[62]+5.196152422706631*coords[0]*coords[2]*coords[3]*fdo[61]+5.196152422706631*coords[0]*coords[1]*coords[3]*fdo[60]+5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[2]*coords[3]*fdo[56]+3.0*coords[1]*coords[3]*fdo[55]+3.0*coords[0]*coords[3]*fdo[54]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[3]*fdo[41]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[33]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[17]+3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[18]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[49]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[35]+3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[52]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[38]+3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[36]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[39]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[61]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[54]+3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[55]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[44]+fdo[30]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[18]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[1]*fdo[17]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[36]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[39]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[25]+fdo[12]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[28]+fdo[15]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[55]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[30]+fdo[16]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[42]+fdo[24]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[43]+fdo[27]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[44]+fdo[29]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[46]+fdo[31]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[61]+1.7320508075688772*coords[0]*fdo[57]+fdo[45]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[34]+3.0*coords[1]*coords[2]*fdo[19]+3.0*coords[0]*coords[2]*fdo[18]+3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[2]*fdo[6]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[1]*coords[2]*fdo[33]+3.0*coords[0]*coords[2]*fdo[32]+3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[2]*fdo[17]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[50]+3.0*coords[1]*coords[2]*fdo[37]+3.0*coords[0]*coords[2]*fdo[36]+3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[2]*fdo[20]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[53]+3.0*coords[1]*coords[2]*fdo[40]+3.0*coords[0]*coords[2]*fdo[39]+3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[2]*fdo[21]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[49]+3.0*coords[0]*coords[2]*fdo[48]+3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[2]*fdo[35]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[52]+3.0*coords[0]*coords[2]*fdo[51]+3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[2]*fdo[38]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[62]+3.0*coords[1]*coords[2]*fdo[56]+3.0*coords[0]*coords[2]*fdo[55]+3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[2]*fdo[41]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[30]+fdo[16]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[61]+3.0*coords[0]*coords[2]*fdo[60]+3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[2]*fdo[54]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[44]+fdo[29]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[19]+1.7320508075688772*coords[1]*fdo[6]+1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[1]*fdo[17]+1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[37]+1.7320508075688772*coords[1]*fdo[20]+1.7320508075688772*coords[0]*fdo[12]+fdo[4]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[40]+1.7320508075688772*coords[1]*fdo[21]+1.7320508075688772*coords[0]*fdo[15]+fdo[5]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[32]+1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[35]+1.7320508075688772*coords[0]*fdo[24]+fdo[10]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[25]+fdo[11]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[38]+1.7320508075688772*coords[0]*fdo[27]+fdo[13]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[28]+fdo[14]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[56]+1.7320508075688772*coords[1]*fdo[41]+1.7320508075688772*coords[0]*fdo[31]+fdo[16]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[48]+1.7320508075688772*coords[0]*fdo[42]+fdo[23]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[51]+1.7320508075688772*coords[0]*fdo[43]+fdo[26]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[54]+1.7320508075688772*coords[0]*fdo[45]+fdo[29]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[46]+fdo[30]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[60]+1.7320508075688772*coords[0]*fdo[57]+fdo[44]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[22]+fdo[9]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[42]+fdo[25]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[43]+fdo[28]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[44]+fdo[30]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[45]+fdo[31]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[47]+fdo[34]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[48]+fdo[36]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[49]+fdo[37]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[51]+fdo[39]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[52]+fdo[40]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[54]+fdo[41]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[57]+fdo[46]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[58]+fdo[50]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[59]+fdo[53]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[60]+fdo[55]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[61]+fdo[56]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[62]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[8]+fdo[3]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[24]+fdo[12]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[27]+fdo[15]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[33]+fdo[19]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[45]+fdo[31]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[49]+fdo[37]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[52]+fdo[40]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[54]+fdo[41]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[62]+1.7320508075688772*coords[0]*fdo[61]+fdo[56]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[22]+3.0*coords[1]*coords[2]*fdo[9]+3.0*coords[0]*coords[2]*fdo[8]+3.0*coords[0]*coords[1]*fdo[7]+1.7320508075688772*coords[2]*fdo[3]+1.7320508075688772*coords[1]*fdo[2]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[42]+3.0*coords[1]*coords[2]*fdo[25]+3.0*coords[0]*coords[2]*fdo[24]+3.0*coords[0]*coords[1]*fdo[23]+1.7320508075688772*coords[2]*fdo[12]+1.7320508075688772*coords[1]*fdo[11]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[2] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[43]+3.0*coords[1]*coords[2]*fdo[28]+3.0*coords[0]*coords[2]*fdo[27]+3.0*coords[0]*coords[1]*fdo[26]+1.7320508075688772*coords[2]*fdo[15]+1.7320508075688772*coords[1]*fdo[14]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[3] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[47]+3.0*coords[1]*coords[2]*fdo[34]+3.0*coords[0]*coords[2]*fdo[33]+3.0*coords[0]*coords[1]*fdo[32]+1.7320508075688772*coords[2]*fdo[19]+1.7320508075688772*coords[1]*fdo[18]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[4] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[57]+3.0*coords[1]*coords[2]*fdo[46]+3.0*coords[0]*coords[2]*fdo[45]+3.0*coords[0]*coords[1]*fdo[44]+1.7320508075688772*coords[2]*fdo[31]+1.7320508075688772*coords[1]*fdo[30]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[5] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[58]+3.0*coords[1]*coords[2]*fdo[50]+3.0*coords[0]*coords[2]*fdo[49]+3.0*coords[0]*coords[1]*fdo[48]+1.7320508075688772*coords[2]*fdo[37]+1.7320508075688772*coords[1]*fdo[36]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[6] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[59]+3.0*coords[1]*coords[2]*fdo[53]+3.0*coords[0]*coords[2]*fdo[52]+3.0*coords[0]*coords[1]*fdo[51]+1.7320508075688772*coords[2]*fdo[40]+1.7320508075688772*coords[1]*fdo[39]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[7] = 0.3535533905932737*(5.196152422706631*coords[0]*coords[1]*coords[2]*fdo[63]+3.0*coords[1]*coords[2]*fdo[62]+3.0*coords[0]*coords[2]*fdo[61]+3.0*coords[0]*coords[1]*fdo[60]+1.7320508075688772*coords[2]*fdo[56]+1.7320508075688772*coords[1]*fdo[55]+1.7320508075688772*coords[0]*fdo[54]+fdo[41]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[8]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[1]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[9]+1.7320508075688772*coords[0]*fdo[7]+fdo[2]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[24]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[10]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[27]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[13]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[33]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[17]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[25]+1.7320508075688772*coords[0]*fdo[23]+fdo[11]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[28]+1.7320508075688772*coords[0]*fdo[26]+fdo[14]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[45]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[29]+fdo[16]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[34]+1.7320508075688772*coords[0]*fdo[32]+fdo[18]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[49]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[35]+fdo[20]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[52]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[38]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[46]+1.7320508075688772*coords[0]*fdo[44]+fdo[30]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[50]+1.7320508075688772*coords[0]*fdo[48]+fdo[36]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[53]+1.7320508075688772*coords[0]*fdo[51]+fdo[39]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[61]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[54]+fdo[41]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[62]+1.7320508075688772*coords[0]*fdo[60]+fdo[55]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[9]+fdo[3]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[22]+fdo[8]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[25]+fdo[12]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[28]+fdo[15]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[30]+fdo[16]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[34]+fdo[19]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[36]+fdo[20]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[39]+fdo[21]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[42]+fdo[24]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[43]+fdo[27]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[44]+fdo[29]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[46]+fdo[31]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[47]+fdo[33]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[48]+fdo[35]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[50]+fdo[37]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[51]+fdo[38]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[53]+fdo[40]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[55]+fdo[41]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[57]+fdo[45]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[58]+fdo[49]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[59]+fdo[52]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[60]+fdo[54]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[62]+fdo[56]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[61]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.5*(3.0*coords[0]*coords[1]*fdo[9]+1.7320508075688772*coords[1]*fdo[3]+1.7320508075688772*coords[0]*fdo[2]+fdo[0]); 
  ftar[1] = 0.5*(3.0*coords[0]*coords[1]*fdo[22]+1.7320508075688772*coords[1]*fdo[8]+1.7320508075688772*coords[0]*fdo[7]+fdo[1]); 
  ftar[2] = 0.5*(3.0*coords[0]*coords[1]*fdo[25]+1.7320508075688772*coords[1]*fdo[12]+1.7320508075688772*coords[0]*fdo[11]+fdo[4]); 
  ftar[3] = 0.5*(3.0*coords[0]*coords[1]*fdo[28]+1.7320508075688772*coords[1]*fdo[15]+1.7320508075688772*coords[0]*fdo[14]+fdo[5]); 
  ftar[4] = 0.5*(3.0*coords[0]*coords[1]*fdo[34]+1.7320508075688772*coords[1]*fdo[19]+1.7320508075688772*coords[0]*fdo[18]+fdo[6]); 
  ftar[5] = 0.5*(3.0*coords[0]*coords[1]*fdo[42]+1.7320508075688772*coords[1]*fdo[24]+1.7320508075688772*coords[0]*fdo[23]+fdo[10]); 
  ftar[6] = 0.5*(3.0*coords[0]*coords[1]*fdo[43]+1.7320508075688772*coords[1]*fdo[27]+1.7320508075688772*coords[0]*fdo[26]+fdo[13]); 
  ftar[7] = 0.5*(3.0*coords[0]*coords[1]*fdo[46]+1.7320508075688772*coords[1]*fdo[31]+1.7320508075688772*coords[0]*fdo[30]+fdo[16]); 
  ftar[8] = 0.5*(3.0*coords[0]*coords[1]*fdo[47]+1.7320508075688772*coords[1]*fdo[33]+1.7320508075688772*coords[0]*fdo[32]+fdo[17]); 
  ftar[9] = 0.5*(3.0*coords[0]*coords[1]*fdo[50]+1.7320508075688772*coords[1]*fdo[37]+1.7320508075688772*coords[0]*fdo[36]+fdo[20]); 
  ftar[10] = 0.5*(3.0*coords[0]*coords[1]*fdo[53]+1.7320508075688772*coords[1]*fdo[40]+1.7320508075688772*coords[0]*fdo[39]+fdo[21]); 
  ftar[11] = 0.5*(3.0*coords[0]*coords[1]*fdo[57]+1.7320508075688772*coords[1]*fdo[45]+1.7320508075688772*coords[0]*fdo[44]+fdo[29]); 
  ftar[12] = 0.5*(3.0*coords[0]*coords[1]*fdo[58]+1.7320508075688772*coords[1]*fdo[49]+1.7320508075688772*coords[0]*fdo[48]+fdo[35]); 
  ftar[13] = 0.5*(3.0*coords[0]*coords[1]*fdo[59]+1.7320508075688772*coords[1]*fdo[52]+1.7320508075688772*coords[0]*fdo[51]+fdo[38]); 
  ftar[14] = 0.5*(3.0*coords[0]*coords[1]*fdo[62]+1.7320508075688772*coords[1]*fdo[56]+1.7320508075688772*coords[0]*fdo[55]+fdo[41]); 
  ftar[15] = 0.5*(3.0*coords[0]*coords[1]*fdo[63]+1.7320508075688772*coords[1]*fdo[61]+1.7320508075688772*coords[0]*fdo[60]+fdo[54]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2(const double *coords, const double *fdo, double *ftar) 
{ 
  // coords: logical coordinates to evaluate the field at.
  // fdo: donor field to get DG coefficients from.
  // ftar: target field whose DG coefficients to populate.


  ftar[0] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[3]+fdo[0]); 
  ftar[1] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[8]+fdo[1]); 
  ftar[2] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[9]+fdo[2]); 
  ftar[3] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[12]+fdo[4]); 
  ftar[4] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[15]+fdo[5]); 
  ftar[5] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[19]+fdo[6]); 
  ftar[6] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[22]+fdo[7]); 
  ftar[7] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[24]+fdo[10]); 
  ftar[8] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[25]+fdo[11]); 
  ftar[9] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[27]+fdo[13]); 
  ftar[10] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[28]+fdo[14]); 
  ftar[11] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[31]+fdo[16]); 
  ftar[12] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[33]+fdo[17]); 
  ftar[13] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[34]+fdo[18]); 
  ftar[14] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[37]+fdo[20]); 
  ftar[15] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[40]+fdo[21]); 
  ftar[16] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[42]+fdo[23]); 
  ftar[17] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[43]+fdo[26]); 
  ftar[18] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[45]+fdo[29]); 
  ftar[19] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[46]+fdo[30]); 
  ftar[20] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[47]+fdo[32]); 
  ftar[21] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[49]+fdo[35]); 
  ftar[22] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[50]+fdo[36]); 
  ftar[23] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[52]+fdo[38]); 
  ftar[24] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[53]+fdo[39]); 
  ftar[25] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[56]+fdo[41]); 
  ftar[26] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[57]+fdo[44]); 
  ftar[27] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[58]+fdo[48]); 
  ftar[28] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[59]+fdo[51]); 
  ftar[29] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[61]+fdo[54]); 
  ftar[30] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[62]+fdo[55]); 
  ftar[31] = 0.7071067811865475*(1.7320508075688772*coords[0]*fdo[63]+fdo[60]); 
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_3_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_34_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 0;
  *num_basis = 1;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2345_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_034_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_134_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_234_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_35_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_035_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0135_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01235_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0235_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_135_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1235_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_235_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_03_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_013_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_023_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_13_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_123_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_23_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_4_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_45_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_045_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0145_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01245_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 1;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 2;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0245_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_145_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1245_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_245_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_04_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_014_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_024_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_14_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_124_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_24_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_5_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_05_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_015_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0125_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 2;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 4;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_025_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_15_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_125_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_25_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_0_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_01_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_012_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 3;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 8;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_02_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_1_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_12_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 4;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 16;
}

GKYL_CU_DH void gkyl_dg_eval_at_coord_proj_6x_ser_p1_eval_dirs_2_target_basis(int *cdim, int *ndim, enum gkyl_basis_type *btype, int *poly_order, int *num_basis) 
{

  *cdim = 1;
  *ndim = 5;
  *btype = GKYL_BASIS_MODAL_SERENDIPITY;
  *poly_order = 1;
  *num_basis = 32;
}

