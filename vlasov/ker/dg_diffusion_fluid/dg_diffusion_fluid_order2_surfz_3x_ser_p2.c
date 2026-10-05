#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[2]*qr[9]+53.66563145999496*coeff[2]*ql[9]-107.3312629199899*coeff[2]*qc[9]-95.26279441628824*coeff[2]*qr[3]+95.26279441628824*coeff[2]*ql[3]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[2])*Jfac; 
  out[1] += 0.0125*(53.66563145999495*coeff[2]*qr[15]+53.66563145999495*coeff[2]*ql[15]-107.3312629199899*coeff[2]*qc[15]-95.26279441628824*coeff[2]*qr[5]+95.26279441628824*coeff[2]*ql[5]+(75.0*qr[1]+75.0*ql[1]-150.0*qc[1])*coeff[2])*Jfac; 
  out[2] += 0.0125*(53.66563145999495*coeff[2]*qr[16]+53.66563145999495*coeff[2]*ql[16]-107.3312629199899*coeff[2]*qc[16]-95.26279441628824*coeff[2]*qr[6]+95.26279441628824*coeff[2]*ql[6]+75.0*coeff[2]*qr[2]+75.0*coeff[2]*ql[2]-150.0*coeff[2]*qc[2])*Jfac; 
  out[3] += 0.003125*(236.2519841186524*coeff[2]*qr[9]-236.2519841186524*coeff[2]*ql[9]-465.0*coeff[2]*qr[3]-465.0*coeff[2]*ql[3]-1710.0*coeff[2]*qc[3]+(381.051177665153*qr[0]-381.051177665153*ql[0])*coeff[2])*Jfac; 
  out[4] += 0.0125*(53.66563145999496*coeff[2]*qr[19]+53.66563145999496*coeff[2]*ql[19]-107.3312629199899*coeff[2]*qc[19]-95.26279441628824*coeff[2]*qr[10]+95.26279441628824*coeff[2]*ql[10]+75.0*coeff[2]*qr[4]+75.0*coeff[2]*ql[4]-150.0*coeff[2]*qc[4])*Jfac; 
  out[5] += 0.003125*(236.2519841186524*coeff[2]*qr[15]-236.2519841186524*coeff[2]*ql[15]-465.0*coeff[2]*qr[5]-465.0*coeff[2]*ql[5]-1710.0*coeff[2]*qc[5]+(381.051177665153*qr[1]-381.051177665153*ql[1])*coeff[2])*Jfac; 
  out[6] += 0.003125*(236.2519841186524*coeff[2]*qr[16]-236.2519841186524*coeff[2]*ql[16]-465.0*coeff[2]*qr[6]-465.0*coeff[2]*ql[6]-1710.0*coeff[2]*qc[6]+381.051177665153*coeff[2]*qr[2]-381.051177665153*coeff[2]*ql[2])*Jfac; 
  out[7] += -0.0125*(95.26279441628826*coeff[2]*qr[13]-95.26279441628826*coeff[2]*ql[13]-75.0*coeff[2]*qr[7]-75.0*coeff[2]*ql[7]+150.0*coeff[2]*qc[7])*Jfac; 
  out[8] += -0.0125*(95.26279441628826*coeff[2]*qr[14]-95.26279441628826*coeff[2]*ql[14]-75.0*coeff[2]*qr[8]-75.0*coeff[2]*ql[8]+150.0*coeff[2]*qc[8])*Jfac; 
  out[9] += -0.015625*(9.0*coeff[2]*qr[9]+9.0*coeff[2]*ql[9]+402.0*coeff[2]*qc[9]+19.36491673103709*coeff[2]*qr[3]-19.36491673103709*coeff[2]*ql[3]+((-26.83281572999748*qr[0])-26.83281572999748*ql[0]+53.66563145999496*qc[0])*coeff[2])*Jfac; 
  out[10] += 0.003125*(236.2519841186524*coeff[2]*qr[19]-236.2519841186524*coeff[2]*ql[19]-465.0*coeff[2]*qr[10]-465.0*coeff[2]*ql[10]-1710.0*coeff[2]*qc[10]+381.051177665153*coeff[2]*qr[4]-381.051177665153*coeff[2]*ql[4])*Jfac; 
  out[11] += -0.0125*(95.26279441628826*coeff[2]*qr[17]-95.26279441628826*coeff[2]*ql[17]-75.0*coeff[2]*qr[11]-75.0*coeff[2]*ql[11]+150.0*coeff[2]*qc[11])*Jfac; 
  out[12] += -0.0125*(95.26279441628826*coeff[2]*qr[18]-95.26279441628826*coeff[2]*ql[18]-75.0*coeff[2]*qr[12]-75.0*coeff[2]*ql[12]+150.0*coeff[2]*qc[12])*Jfac; 
  out[13] += -0.003125*(465.0*coeff[2]*qr[13]+465.0*coeff[2]*ql[13]+1710.0*coeff[2]*qc[13]-381.051177665153*coeff[2]*qr[7]+381.051177665153*coeff[2]*ql[7])*Jfac; 
  out[14] += -0.003125*(465.0*coeff[2]*qr[14]+465.0*coeff[2]*ql[14]+1710.0*coeff[2]*qc[14]-381.051177665153*coeff[2]*qr[8]+381.051177665153*coeff[2]*ql[8])*Jfac; 
  out[15] += -0.015625*(9.0*coeff[2]*qr[15]+9.0*coeff[2]*ql[15]+402.0*coeff[2]*qc[15]+19.36491673103708*coeff[2]*qr[5]-19.36491673103708*coeff[2]*ql[5]+((-26.83281572999747*qr[1])-26.83281572999747*ql[1]+53.66563145999495*qc[1])*coeff[2])*Jfac; 
  out[16] += -0.015625*(9.0*coeff[2]*qr[16]+9.0*coeff[2]*ql[16]+402.0*coeff[2]*qc[16]+19.36491673103708*coeff[2]*qr[6]-19.36491673103708*coeff[2]*ql[6]-26.83281572999747*coeff[2]*qr[2]-26.83281572999747*coeff[2]*ql[2]+53.66563145999495*coeff[2]*qc[2])*Jfac; 
  out[17] += -0.003125*(465.0*coeff[2]*qr[17]+465.0*coeff[2]*ql[17]+1710.0*coeff[2]*qc[17]-381.051177665153*coeff[2]*qr[11]+381.051177665153*coeff[2]*ql[11])*Jfac; 
  out[18] += -0.003125*(465.0*coeff[2]*qr[18]+465.0*coeff[2]*ql[18]+1710.0*coeff[2]*qc[18]-381.051177665153*coeff[2]*qr[12]+381.051177665153*coeff[2]*ql[12])*Jfac; 
  out[19] += -0.015625*(9.0*coeff[2]*qr[19]+9.0*coeff[2]*ql[19]+402.0*coeff[2]*qc[19]+19.36491673103709*coeff[2]*qr[10]-19.36491673103709*coeff[2]*ql[10]-26.83281572999748*coeff[2]*qr[4]-26.83281572999748*coeff[2]*ql[4]+53.66563145999496*coeff[2]*qc[4])*Jfac; 

  return 0.;

}

