#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[0]*qr[7]+53.66563145999496*coeff[0]*ql[7]-107.3312629199899*coeff[0]*qc[7]-95.26279441628824*coeff[0]*qr[1]+95.26279441628824*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.003125*(236.2519841186524*coeff[0]*qr[7]-236.2519841186524*coeff[0]*ql[7]-465.0*coeff[0]*qr[1]-465.0*coeff[0]*ql[1]-1710.0*coeff[0]*qc[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0125*(53.66563145999495*coeff[0]*qr[11]+53.66563145999495*coeff[0]*ql[11]-107.3312629199899*coeff[0]*qc[11]-95.26279441628824*coeff[0]*qr[4]+95.26279441628824*coeff[0]*ql[4]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.0125*(53.66563145999495*coeff[0]*qr[13]+53.66563145999495*coeff[0]*ql[13]-107.3312629199899*coeff[0]*qc[13]-95.26279441628824*coeff[0]*qr[5]+95.26279441628824*coeff[0]*ql[5]+75.0*coeff[0]*qr[3]+75.0*coeff[0]*ql[3]-150.0*coeff[0]*qc[3])*Jfac; 
  out[4] += 0.003125*(236.2519841186524*coeff[0]*qr[11]-236.2519841186524*coeff[0]*ql[11]-465.0*coeff[0]*qr[4]-465.0*coeff[0]*ql[4]-1710.0*coeff[0]*qc[4]+381.051177665153*coeff[0]*qr[2]-381.051177665153*coeff[0]*ql[2])*Jfac; 
  out[5] += 0.003125*(236.2519841186524*coeff[0]*qr[13]-236.2519841186524*coeff[0]*ql[13]-465.0*coeff[0]*qr[5]-465.0*coeff[0]*ql[5]-1710.0*coeff[0]*qc[5]+381.051177665153*coeff[0]*qr[3]-381.051177665153*coeff[0]*ql[3])*Jfac; 
  out[6] += 0.0125*(53.66563145999496*coeff[0]*qr[17]+53.66563145999496*coeff[0]*ql[17]-107.3312629199899*coeff[0]*qc[17]-95.26279441628824*coeff[0]*qr[10]+95.26279441628824*coeff[0]*ql[10]+75.0*coeff[0]*qr[6]+75.0*coeff[0]*ql[6]-150.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.015625*(9.0*coeff[0]*qr[7]+9.0*coeff[0]*ql[7]+402.0*coeff[0]*qc[7]+19.36491673103709*coeff[0]*qr[1]-19.36491673103709*coeff[0]*ql[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 
  out[8] += -0.0125*(95.26279441628826*coeff[0]*qr[12]-95.26279441628826*coeff[0]*ql[12]-75.0*coeff[0]*qr[8]-75.0*coeff[0]*ql[8]+150.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0125*(95.26279441628826*coeff[0]*qr[15]-95.26279441628826*coeff[0]*ql[15]-75.0*coeff[0]*qr[9]-75.0*coeff[0]*ql[9]+150.0*coeff[0]*qc[9])*Jfac; 
  out[10] += 0.003125*(236.2519841186524*coeff[0]*qr[17]-236.2519841186524*coeff[0]*ql[17]-465.0*coeff[0]*qr[10]-465.0*coeff[0]*ql[10]-1710.0*coeff[0]*qc[10]+381.051177665153*coeff[0]*qr[6]-381.051177665153*coeff[0]*ql[6])*Jfac; 
  out[11] += -0.015625*(9.0*coeff[0]*qr[11]+9.0*coeff[0]*ql[11]+402.0*coeff[0]*qc[11]+19.36491673103708*coeff[0]*qr[4]-19.36491673103708*coeff[0]*ql[4]-26.83281572999747*coeff[0]*qr[2]-26.83281572999747*coeff[0]*ql[2]+53.66563145999495*coeff[0]*qc[2])*Jfac; 
  out[12] += -0.003125*(465.0*coeff[0]*qr[12]+465.0*coeff[0]*ql[12]+1710.0*coeff[0]*qc[12]-381.051177665153*coeff[0]*qr[8]+381.051177665153*coeff[0]*ql[8])*Jfac; 
  out[13] += -0.015625*(9.0*coeff[0]*qr[13]+9.0*coeff[0]*ql[13]+402.0*coeff[0]*qc[13]+19.36491673103708*coeff[0]*qr[5]-19.36491673103708*coeff[0]*ql[5]-26.83281572999747*coeff[0]*qr[3]-26.83281572999747*coeff[0]*ql[3]+53.66563145999495*coeff[0]*qc[3])*Jfac; 
  out[14] += -0.0125*(95.26279441628826*coeff[0]*qr[18]-95.26279441628826*coeff[0]*ql[18]-75.0*coeff[0]*qr[14]-75.0*coeff[0]*ql[14]+150.0*coeff[0]*qc[14])*Jfac; 
  out[15] += -0.003125*(465.0*coeff[0]*qr[15]+465.0*coeff[0]*ql[15]+1710.0*coeff[0]*qc[15]-381.051177665153*coeff[0]*qr[9]+381.051177665153*coeff[0]*ql[9])*Jfac; 
  out[16] += -0.0125*(95.26279441628826*coeff[0]*qr[19]-95.26279441628826*coeff[0]*ql[19]-75.0*coeff[0]*qr[16]-75.0*coeff[0]*ql[16]+150.0*coeff[0]*qc[16])*Jfac; 
  out[17] += -0.015625*(9.0*coeff[0]*qr[17]+9.0*coeff[0]*ql[17]+402.0*coeff[0]*qc[17]+19.36491673103709*coeff[0]*qr[10]-19.36491673103709*coeff[0]*ql[10]-26.83281572999748*coeff[0]*qr[6]-26.83281572999748*coeff[0]*ql[6]+53.66563145999496*coeff[0]*qc[6])*Jfac; 
  out[18] += -0.003125*(465.0*coeff[0]*qr[18]+465.0*coeff[0]*ql[18]+1710.0*coeff[0]*qc[18]-381.051177665153*coeff[0]*qr[14]+381.051177665153*coeff[0]*ql[14])*Jfac; 
  out[19] += -0.003125*(465.0*coeff[0]*qr[19]+465.0*coeff[0]*ql[19]+1710.0*coeff[0]*qc[19]-381.051177665153*coeff[0]*qr[16]+381.051177665153*coeff[0]*ql[16])*Jfac; 

  return 0.;

}

