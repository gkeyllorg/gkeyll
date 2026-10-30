#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -(0.0625*(25.980762113533157*coeff[0]*qr[1]-25.980762113533157*coeff[0]*ql[1]-15.0*coeff[0]*qr[0]-15.0*coeff[0]*ql[0]+30.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0625*(33.0*coeff[0]*qr[1]+33.0*coeff[0]*ql[1]+114.0*coeff[0]*qc[1]-25.980762113533157*coeff[0]*qr[0]+25.980762113533157*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0625*(25.980762113533157*coeff[0]*qr[5]-25.980762113533157*coeff[0]*ql[5]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0625*(25.980762113533157*coeff[0]*qr[6]-25.980762113533157*coeff[0]*ql[6]-15.0*coeff[0]*qr[3]-15.0*coeff[0]*ql[3]+30.0*coeff[0]*qc[3])*Jfac); 
  out[4] += -(0.0625*(25.980762113533157*coeff[0]*qr[8]-25.980762113533157*coeff[0]*ql[8]-15.0*coeff[0]*qr[4]-15.0*coeff[0]*ql[4]+30.0*coeff[0]*qc[4])*Jfac); 
  out[5] += -(0.0625*(33.0*coeff[0]*qr[5]+33.0*coeff[0]*ql[5]+114.0*coeff[0]*qc[5]-25.980762113533157*coeff[0]*qr[2]+25.980762113533157*coeff[0]*ql[2])*Jfac); 
  out[6] += -(0.0625*(33.0*coeff[0]*qr[6]+33.0*coeff[0]*ql[6]+114.0*coeff[0]*qc[6]-25.980762113533157*coeff[0]*qr[3]+25.980762113533157*coeff[0]*ql[3])*Jfac); 
  out[7] += -(0.0625*(25.980762113533157*coeff[0]*qr[13]-25.980762113533157*coeff[0]*ql[13]-15.0*coeff[0]*qr[7]-15.0*coeff[0]*ql[7]+30.0*coeff[0]*qc[7])*Jfac); 
  out[8] += -(0.0625*(33.0*coeff[0]*qr[8]+33.0*coeff[0]*ql[8]+114.0*coeff[0]*qc[8]-25.980762113533157*coeff[0]*qr[4]+25.980762113533157*coeff[0]*ql[4])*Jfac); 
  out[9] += -(0.0625*(25.980762113533157*coeff[0]*qr[14]-25.980762113533157*coeff[0]*ql[14]-15.0*coeff[0]*qr[9]-15.0*coeff[0]*ql[9]+30.0*coeff[0]*qc[9])*Jfac); 
  out[10] += -(0.0625*(25.980762113533157*coeff[0]*qr[15]-25.980762113533157*coeff[0]*ql[15]-15.0*coeff[0]*qr[10]-15.0*coeff[0]*ql[10]+30.0*coeff[0]*qc[10])*Jfac); 
  out[11] += -(0.0625*(25.98076211353316*coeff[0]*qr[17]-25.98076211353316*coeff[0]*ql[17]-15.0*coeff[0]*qr[11]-15.0*coeff[0]*ql[11]+30.0*coeff[0]*qc[11])*Jfac); 
  out[12] += -(0.0625*(25.98076211353316*coeff[0]*qr[20]-25.98076211353316*coeff[0]*ql[20]-15.0*coeff[0]*qr[12]-15.0*coeff[0]*ql[12]+30.0*coeff[0]*qc[12])*Jfac); 
  out[13] += -(0.0625*(33.0*coeff[0]*qr[13]+33.0*coeff[0]*ql[13]+114.0*coeff[0]*qc[13]-25.980762113533157*coeff[0]*qr[7]+25.980762113533157*coeff[0]*ql[7])*Jfac); 
  out[14] += -(0.0625*(33.0*coeff[0]*qr[14]+33.0*coeff[0]*ql[14]+114.0*coeff[0]*qc[14]-25.980762113533157*coeff[0]*qr[9]+25.980762113533157*coeff[0]*ql[9])*Jfac); 
  out[15] += -(0.0625*(33.0*coeff[0]*qr[15]+33.0*coeff[0]*ql[15]+114.0*coeff[0]*qc[15]-25.980762113533157*coeff[0]*qr[10]+25.980762113533157*coeff[0]*ql[10])*Jfac); 
  out[16] += -(0.0625*(25.980762113533157*coeff[0]*qr[23]-25.980762113533157*coeff[0]*ql[23]-15.0*coeff[0]*qr[16]-15.0*coeff[0]*ql[16]+30.0*coeff[0]*qc[16])*Jfac); 
  out[17] += -(0.0625*(33.0*coeff[0]*qr[17]+33.0*coeff[0]*ql[17]+114.0*coeff[0]*qc[17]-25.98076211353316*coeff[0]*qr[11]+25.98076211353316*coeff[0]*ql[11])*Jfac); 
  out[18] += -(0.0625*(25.98076211353316*coeff[0]*qr[24]-25.98076211353316*coeff[0]*ql[24]-15.0*coeff[0]*qr[18]-15.0*coeff[0]*ql[18]+30.0*coeff[0]*qc[18])*Jfac); 
  out[19] += -(0.0625*(25.98076211353316*coeff[0]*qr[25]-25.98076211353316*coeff[0]*ql[25]-15.0*coeff[0]*qr[19]-15.0*coeff[0]*ql[19]+30.0*coeff[0]*qc[19])*Jfac); 
  out[20] += -(0.0625*(33.0*coeff[0]*qr[20]+33.0*coeff[0]*ql[20]+114.0*coeff[0]*qc[20]-25.98076211353316*coeff[0]*qr[12]+25.98076211353316*coeff[0]*ql[12])*Jfac); 
  out[21] += -(0.0625*(25.98076211353316*coeff[0]*qr[27]-25.98076211353316*coeff[0]*ql[27]-15.0*coeff[0]*qr[21]-15.0*coeff[0]*ql[21]+30.0*coeff[0]*qc[21])*Jfac); 
  out[22] += -(0.0625*(25.98076211353316*coeff[0]*qr[28]-25.98076211353316*coeff[0]*ql[28]-15.0*coeff[0]*qr[22]-15.0*coeff[0]*ql[22]+30.0*coeff[0]*qc[22])*Jfac); 
  out[23] += -(0.0625*(33.0*coeff[0]*qr[23]+33.0*coeff[0]*ql[23]+114.0*coeff[0]*qc[23]-25.980762113533157*coeff[0]*qr[16]+25.980762113533157*coeff[0]*ql[16])*Jfac); 
  out[24] += -(0.0625*(33.0*coeff[0]*qr[24]+33.0*coeff[0]*ql[24]+114.0*coeff[0]*qc[24]-25.98076211353316*coeff[0]*qr[18]+25.98076211353316*coeff[0]*ql[18])*Jfac); 
  out[25] += -(0.0625*(33.0*coeff[0]*qr[25]+33.0*coeff[0]*ql[25]+114.0*coeff[0]*qc[25]-25.98076211353316*coeff[0]*qr[19]+25.98076211353316*coeff[0]*ql[19])*Jfac); 
  out[26] += -(0.0625*(25.98076211353316*coeff[0]*qr[31]-25.98076211353316*coeff[0]*ql[31]-15.0*coeff[0]*qr[26]-15.0*coeff[0]*ql[26]+30.0*coeff[0]*qc[26])*Jfac); 
  out[27] += -(0.0625*(33.0*coeff[0]*qr[27]+33.0*coeff[0]*ql[27]+114.0*coeff[0]*qc[27]-25.98076211353316*coeff[0]*qr[21]+25.98076211353316*coeff[0]*ql[21])*Jfac); 
  out[28] += -(0.0625*(33.0*coeff[0]*qr[28]+33.0*coeff[0]*ql[28]+114.0*coeff[0]*qc[28]-25.98076211353316*coeff[0]*qr[22]+25.98076211353316*coeff[0]*ql[22])*Jfac); 
  out[29] += -(0.0625*(25.98076211353316*coeff[0]*qr[32]-25.98076211353316*coeff[0]*ql[32]-15.0*coeff[0]*qr[29]-15.0*coeff[0]*ql[29]+30.0*coeff[0]*qc[29])*Jfac); 
  out[30] += -(0.0625*(25.980762113533157*coeff[0]*qr[33]-25.980762113533157*coeff[0]*ql[33]-15.0*coeff[0]*qr[30]-15.0*coeff[0]*ql[30]+30.0*coeff[0]*qc[30])*Jfac); 
  out[31] += -(0.0625*(33.0*coeff[0]*qr[31]+33.0*coeff[0]*ql[31]+114.0*coeff[0]*qc[31]-25.98076211353316*coeff[0]*qr[26]+25.98076211353316*coeff[0]*ql[26])*Jfac); 
  out[32] += -(0.0625*(33.0*coeff[0]*qr[32]+33.0*coeff[0]*ql[32]+114.0*coeff[0]*qc[32]-25.98076211353316*coeff[0]*qr[29]+25.98076211353316*coeff[0]*ql[29])*Jfac); 
  out[33] += -(0.0625*(33.0*coeff[0]*qr[33]+33.0*coeff[0]*ql[33]+114.0*coeff[0]*qc[33]-25.980762113533157*coeff[0]*qr[30]+25.980762113533157*coeff[0]*ql[30])*Jfac); 
  out[34] += -(0.0625*(25.980762113533157*coeff[0]*qr[35]-25.980762113533157*coeff[0]*ql[35]-15.0*coeff[0]*qr[34]-15.0*coeff[0]*ql[34]+30.0*coeff[0]*qc[34])*Jfac); 
  out[35] += -(0.0625*(33.0*coeff[0]*qr[35]+33.0*coeff[0]*ql[35]+114.0*coeff[0]*qc[35]-25.980762113533157*coeff[0]*qr[34]+25.980762113533157*coeff[0]*ql[34])*Jfac); 

  return 0.;

}

