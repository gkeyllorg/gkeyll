#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.0625*(8.660254037844386*coeff[0]*qr[1]-8.660254037844386*coeff[0]*ql[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0625*(7.0*coeff[0]*qr[1]+7.0*coeff[0]*ql[1]+46.0*coeff[0]*qc[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0625*(8.660254037844386*coeff[0]*qr[4]-8.660254037844386*coeff[0]*ql[4]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0625*(8.660254037844386*coeff[0]*qr[5]-8.660254037844386*coeff[0]*ql[5]-9.0*coeff[0]*qr[3]-9.0*coeff[0]*ql[3]+18.0*coeff[0]*qc[3])*Jfac); 
  out[4] += -(0.0625*(7.0*coeff[0]*qr[4]+7.0*coeff[0]*ql[4]+46.0*coeff[0]*qc[4]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*Jfac); 
  out[5] += -(0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844386*coeff[0]*qr[3]+8.660254037844386*coeff[0]*ql[3])*Jfac); 
  out[6] += -(0.0625*(8.660254037844386*coeff[0]*qr[9]-8.660254037844386*coeff[0]*ql[9]-9.0*coeff[0]*qr[6]-9.0*coeff[0]*ql[6]+18.0*coeff[0]*qc[6])*Jfac); 
  out[7] += -(0.0625*(8.660254037844387*coeff[0]*qr[10]-8.660254037844387*coeff[0]*ql[10]-9.0*coeff[0]*qr[7]-9.0*coeff[0]*ql[7]+18.0*coeff[0]*qc[7])*Jfac); 
  out[8] += -(0.0625*(8.660254037844387*coeff[0]*qr[12]-8.660254037844387*coeff[0]*ql[12]-9.0*coeff[0]*qr[8]-9.0*coeff[0]*ql[8]+18.0*coeff[0]*qc[8])*Jfac); 
  out[9] += -(0.0625*(7.0*coeff[0]*qr[9]+7.0*coeff[0]*ql[9]+46.0*coeff[0]*qc[9]-8.660254037844386*coeff[0]*qr[6]+8.660254037844386*coeff[0]*ql[6])*Jfac); 
  out[10] += -(0.0625*(7.0*coeff[0]*qr[10]+7.0*coeff[0]*ql[10]+46.0*coeff[0]*qc[10]-8.660254037844387*coeff[0]*qr[7]+8.660254037844387*coeff[0]*ql[7])*Jfac); 
  out[11] += -(0.0625*(8.660254037844387*coeff[0]*qr[14]-8.660254037844387*coeff[0]*ql[14]-9.0*coeff[0]*qr[11]-9.0*coeff[0]*ql[11]+18.0*coeff[0]*qc[11])*Jfac); 
  out[12] += -(0.0625*(7.0*coeff[0]*qr[12]+7.0*coeff[0]*ql[12]+46.0*coeff[0]*qc[12]-8.660254037844387*coeff[0]*qr[8]+8.660254037844387*coeff[0]*ql[8])*Jfac); 
  out[13] += -(0.0625*(8.660254037844387*coeff[0]*qr[15]-8.660254037844387*coeff[0]*ql[15]-9.0*coeff[0]*qr[13]-9.0*coeff[0]*ql[13]+18.0*coeff[0]*qc[13])*Jfac); 
  out[14] += -(0.0625*(7.0*coeff[0]*qr[14]+7.0*coeff[0]*ql[14]+46.0*coeff[0]*qc[14]-8.660254037844387*coeff[0]*qr[11]+8.660254037844387*coeff[0]*ql[11])*Jfac); 
  out[15] += -(0.0625*(7.0*coeff[0]*qr[15]+7.0*coeff[0]*ql[15]+46.0*coeff[0]*qc[15]-8.660254037844387*coeff[0]*qr[13]+8.660254037844387*coeff[0]*ql[13])*Jfac); 
  out[16] += -(0.0625*(8.660254037844386*coeff[0]*qr[17]-8.660254037844386*coeff[0]*ql[17]-9.0*coeff[0]*qr[16]-9.0*coeff[0]*ql[16]+18.0*coeff[0]*qc[16])*Jfac); 
  out[17] += -(0.0625*(7.0*coeff[0]*qr[17]+7.0*coeff[0]*ql[17]+46.0*coeff[0]*qc[17]-8.660254037844386*coeff[0]*qr[16]+8.660254037844386*coeff[0]*ql[16])*Jfac); 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_vlasov_order2_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.044194173824159196*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[1]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[1]+30.0*coeff[1]*qc[1]+(15.58845726811989*ql[0]-15.58845726811989*qr[0])*coeff[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[1]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[1]+46.0*coeff[0]*qc[1]+(-(15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.044194173824159196*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[4]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[4]+30.0*coeff[1]*qc[4]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[2]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[2]+18.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.044194173824159196*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[5]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[5]+30.0*coeff[1]*qc[5]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[3]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[3]+18.0*coeff[0]*qc[3])*Jfac); 
  out[4] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[4]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[4]+46.0*coeff[0]*qc[4]+(-(15.0*coeff[1])-8.660254037844386*coeff[0])*qr[2]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[2]+30.0*coeff[1]*qc[2])*Jfac); 
  out[5] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[5]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[5]+46.0*coeff[0]*qc[5]+(-(15.0*coeff[1])-8.660254037844386*coeff[0])*qr[3]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[3]+30.0*coeff[1]*qc[3])*Jfac); 
  out[6] += -(0.044194173824159196*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[9]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[9]+30.0*coeff[1]*qc[9]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[6]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[6]+18.0*coeff[0]*qc[6])*Jfac); 
  out[7] += -(0.044194173824159196*((15.000000000000002*coeff[1]+8.660254037844387*coeff[0])*qr[10]+(15.000000000000002*coeff[1]-8.660254037844387*coeff[0])*ql[10]+30.000000000000004*coeff[1]*qc[10]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[7]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[7]+18.0*coeff[0]*qc[7])*Jfac); 
  out[8] += -(0.044194173824159196*((15.000000000000002*coeff[1]+8.660254037844387*coeff[0])*qr[12]+(15.000000000000002*coeff[1]-8.660254037844387*coeff[0])*ql[12]+30.000000000000004*coeff[1]*qc[12]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[8]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[8]+18.0*coeff[0]*qc[8])*Jfac); 
  out[9] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[9]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[9]+46.0*coeff[0]*qc[9]+(-(15.0*coeff[1])-8.660254037844386*coeff[0])*qr[6]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[6]+30.0*coeff[1]*qc[6])*Jfac); 
  out[10] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[10]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[10]+46.0*coeff[0]*qc[10]+(-(15.000000000000002*coeff[1])-8.660254037844387*coeff[0])*qr[7]+(8.660254037844387*coeff[0]-15.000000000000002*coeff[1])*ql[7]+30.000000000000004*coeff[1]*qc[7])*Jfac); 
  out[11] += -(0.044194173824159196*((15.000000000000002*coeff[1]+8.660254037844387*coeff[0])*qr[14]+(15.000000000000002*coeff[1]-8.660254037844387*coeff[0])*ql[14]+30.000000000000004*coeff[1]*qc[14]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[11]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[11]+18.0*coeff[0]*qc[11])*Jfac); 
  out[12] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[12]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[12]+46.0*coeff[0]*qc[12]+(-(15.000000000000002*coeff[1])-8.660254037844387*coeff[0])*qr[8]+(8.660254037844387*coeff[0]-15.000000000000002*coeff[1])*ql[8]+30.000000000000004*coeff[1]*qc[8])*Jfac); 
  out[13] += -(0.044194173824159196*((15.000000000000002*coeff[1]+8.660254037844387*coeff[0])*qr[15]+(15.000000000000002*coeff[1]-8.660254037844387*coeff[0])*ql[15]+30.000000000000004*coeff[1]*qc[15]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[13]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[13]+18.0*coeff[0]*qc[13])*Jfac); 
  out[14] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[14]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[14]+46.0*coeff[0]*qc[14]+(-(15.000000000000002*coeff[1])-8.660254037844387*coeff[0])*qr[11]+(8.660254037844387*coeff[0]-15.000000000000002*coeff[1])*ql[11]+30.000000000000004*coeff[1]*qc[11])*Jfac); 
  out[15] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[15]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[15]+46.0*coeff[0]*qc[15]+(-(15.000000000000002*coeff[1])-8.660254037844387*coeff[0])*qr[13]+(8.660254037844387*coeff[0]-15.000000000000002*coeff[1])*ql[13]+30.000000000000004*coeff[1]*qc[13])*Jfac); 
  out[16] += -(0.044194173824159196*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[17]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[17]+30.0*coeff[1]*qc[17]+(-(15.58845726811989*coeff[1])-9.0*coeff[0])*qr[16]+(15.58845726811989*coeff[1]-9.0*coeff[0])*ql[16]+18.0*coeff[0]*qc[16])*Jfac); 
  out[17] += -(0.044194173824159196*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[17]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[17]+46.0*coeff[0]*qc[17]+(-(15.0*coeff[1])-8.660254037844386*coeff[0])*qr[16]+(8.660254037844386*coeff[0]-15.0*coeff[1])*ql[16]+30.0*coeff[1]*qc[16])*Jfac); 

  return 0.;

}

