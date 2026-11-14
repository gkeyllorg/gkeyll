#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[1]*qr[8]+53.66563145999496*coeff[1]*ql[8]-107.3312629199899*coeff[1]*qc[8]-95.26279441628824*coeff[1]*qr[2]+95.26279441628824*coeff[1]*ql[2]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0125*(53.66563145999495*coeff[1]*qr[12]+53.66563145999495*coeff[1]*ql[12]-107.3312629199899*coeff[1]*qc[12]-95.26279441628824*coeff[1]*qr[4]+95.26279441628824*coeff[1]*ql[4]+75.0*coeff[1]*qr[1]+75.0*coeff[1]*ql[1]-150.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.003125*(236.2519841186524*coeff[1]*qr[8]-236.2519841186524*coeff[1]*ql[8]-465.0*coeff[1]*qr[2]-465.0*coeff[1]*ql[2]-1710.0*coeff[1]*qc[2]+(381.051177665153*qr[0]-381.051177665153*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.0125*(53.66563145999495*coeff[1]*qr[14]+53.66563145999495*coeff[1]*ql[14]-107.3312629199899*coeff[1]*qc[14]-95.26279441628824*coeff[1]*qr[6]+95.26279441628824*coeff[1]*ql[6]+75.0*coeff[1]*qr[3]+75.0*coeff[1]*ql[3]-150.0*coeff[1]*qc[3])*Jfac; 
  out[4] += 0.003125*(236.2519841186524*coeff[1]*qr[12]-236.2519841186524*coeff[1]*ql[12]-465.0*coeff[1]*qr[4]-465.0*coeff[1]*ql[4]-1710.0*coeff[1]*qc[4]+381.051177665153*coeff[1]*qr[1]-381.051177665153*coeff[1]*ql[1])*Jfac; 
  out[5] += 0.0125*(53.66563145999496*coeff[1]*qr[18]+53.66563145999496*coeff[1]*ql[18]-107.3312629199899*coeff[1]*qc[18]-95.26279441628824*coeff[1]*qr[10]+95.26279441628824*coeff[1]*ql[10]+75.0*coeff[1]*qr[5]+75.0*coeff[1]*ql[5]-150.0*coeff[1]*qc[5])*Jfac; 
  out[6] += 0.003125*(236.2519841186524*coeff[1]*qr[14]-236.2519841186524*coeff[1]*ql[14]-465.0*coeff[1]*qr[6]-465.0*coeff[1]*ql[6]-1710.0*coeff[1]*qc[6]+381.051177665153*coeff[1]*qr[3]-381.051177665153*coeff[1]*ql[3])*Jfac; 
  out[7] += -0.0125*(95.26279441628826*coeff[1]*qr[11]-95.26279441628826*coeff[1]*ql[11]-75.0*coeff[1]*qr[7]-75.0*coeff[1]*ql[7]+150.0*coeff[1]*qc[7])*Jfac; 
  out[8] += -0.015625*(9.0*coeff[1]*qr[8]+9.0*coeff[1]*ql[8]+402.0*coeff[1]*qc[8]+19.36491673103709*coeff[1]*qr[2]-19.36491673103709*coeff[1]*ql[2]+((-26.83281572999748*qr[0])-26.83281572999748*ql[0]+53.66563145999496*qc[0])*coeff[1])*Jfac; 
  out[9] += -0.0125*(95.26279441628826*coeff[1]*qr[16]-95.26279441628826*coeff[1]*ql[16]-75.0*coeff[1]*qr[9]-75.0*coeff[1]*ql[9]+150.0*coeff[1]*qc[9])*Jfac; 
  out[10] += 0.003125*(236.2519841186524*coeff[1]*qr[18]-236.2519841186524*coeff[1]*ql[18]-465.0*coeff[1]*qr[10]-465.0*coeff[1]*ql[10]-1710.0*coeff[1]*qc[10]+381.051177665153*coeff[1]*qr[5]-381.051177665153*coeff[1]*ql[5])*Jfac; 
  out[11] += -0.003125*(465.0*coeff[1]*qr[11]+465.0*coeff[1]*ql[11]+1710.0*coeff[1]*qc[11]-381.051177665153*coeff[1]*qr[7]+381.051177665153*coeff[1]*ql[7])*Jfac; 
  out[12] += -0.015625*(9.0*coeff[1]*qr[12]+9.0*coeff[1]*ql[12]+402.0*coeff[1]*qc[12]+19.36491673103708*coeff[1]*qr[4]-19.36491673103708*coeff[1]*ql[4]-26.83281572999747*coeff[1]*qr[1]-26.83281572999747*coeff[1]*ql[1]+53.66563145999495*coeff[1]*qc[1])*Jfac; 
  out[13] += -0.0125*(95.26279441628826*coeff[1]*qr[17]-95.26279441628826*coeff[1]*ql[17]-75.0*coeff[1]*qr[13]-75.0*coeff[1]*ql[13]+150.0*coeff[1]*qc[13])*Jfac; 
  out[14] += -0.015625*(9.0*coeff[1]*qr[14]+9.0*coeff[1]*ql[14]+402.0*coeff[1]*qc[14]+19.36491673103708*coeff[1]*qr[6]-19.36491673103708*coeff[1]*ql[6]-26.83281572999747*coeff[1]*qr[3]-26.83281572999747*coeff[1]*ql[3]+53.66563145999495*coeff[1]*qc[3])*Jfac; 
  out[15] += -0.0125*(95.26279441628826*coeff[1]*qr[19]-95.26279441628826*coeff[1]*ql[19]-75.0*coeff[1]*qr[15]-75.0*coeff[1]*ql[15]+150.0*coeff[1]*qc[15])*Jfac; 
  out[16] += -0.003125*(465.0*coeff[1]*qr[16]+465.0*coeff[1]*ql[16]+1710.0*coeff[1]*qc[16]-381.051177665153*coeff[1]*qr[9]+381.051177665153*coeff[1]*ql[9])*Jfac; 
  out[17] += -0.003125*(465.0*coeff[1]*qr[17]+465.0*coeff[1]*ql[17]+1710.0*coeff[1]*qc[17]-381.051177665153*coeff[1]*qr[13]+381.051177665153*coeff[1]*ql[13])*Jfac; 
  out[18] += -0.015625*(9.0*coeff[1]*qr[18]+9.0*coeff[1]*ql[18]+402.0*coeff[1]*qc[18]+19.36491673103709*coeff[1]*qr[10]-19.36491673103709*coeff[1]*ql[10]-26.83281572999748*coeff[1]*qr[5]-26.83281572999748*coeff[1]*ql[5]+53.66563145999496*coeff[1]*qc[5])*Jfac; 
  out[19] += -0.003125*(465.0*coeff[1]*qr[19]+465.0*coeff[1]*ql[19]+1710.0*coeff[1]*qc[19]-381.051177665153*coeff[1]*qr[15]+381.051177665153*coeff[1]*ql[15])*Jfac; 

  return 0.;

}

