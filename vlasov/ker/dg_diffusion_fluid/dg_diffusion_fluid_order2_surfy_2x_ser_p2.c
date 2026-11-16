#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[1]*qr[5]+53.66563145999496*coeff[1]*ql[5]-107.3312629199899*coeff[1]*qc[5]-95.26279441628824*coeff[1]*qr[2]+95.26279441628824*coeff[1]*ql[2]+(75.0*qr[0]+75.0*ql[0]-150.0*qc[0])*coeff[1])*Jfac; 
  out[1] += 0.0125*(53.66563145999495*coeff[1]*qr[7]+53.66563145999495*coeff[1]*ql[7]-107.3312629199899*coeff[1]*qc[7]-95.26279441628824*coeff[1]*qr[3]+95.26279441628824*coeff[1]*ql[3]+75.0*coeff[1]*qr[1]+75.0*coeff[1]*ql[1]-150.0*coeff[1]*qc[1])*Jfac; 
  out[2] += 0.003125*(236.2519841186524*coeff[1]*qr[5]-236.2519841186524*coeff[1]*ql[5]-465.0*coeff[1]*qr[2]-465.0*coeff[1]*ql[2]-1710.0*coeff[1]*qc[2]+(381.051177665153*qr[0]-381.051177665153*ql[0])*coeff[1])*Jfac; 
  out[3] += 0.003125*(236.2519841186524*coeff[1]*qr[7]-236.2519841186524*coeff[1]*ql[7]-465.0*coeff[1]*qr[3]-465.0*coeff[1]*ql[3]-1710.0*coeff[1]*qc[3]+381.051177665153*coeff[1]*qr[1]-381.051177665153*coeff[1]*ql[1])*Jfac; 
  out[4] += -0.0125*(95.26279441628826*coeff[1]*qr[6]-95.26279441628826*coeff[1]*ql[6]-75.0*coeff[1]*qr[4]-75.0*coeff[1]*ql[4]+150.0*coeff[1]*qc[4])*Jfac; 
  out[5] += -0.015625*(9.0*coeff[1]*qr[5]+9.0*coeff[1]*ql[5]+402.0*coeff[1]*qc[5]+19.36491673103709*coeff[1]*qr[2]-19.36491673103709*coeff[1]*ql[2]+((-26.83281572999748*qr[0])-26.83281572999748*ql[0]+53.66563145999496*qc[0])*coeff[1])*Jfac; 
  out[6] += -0.003125*(465.0*coeff[1]*qr[6]+465.0*coeff[1]*ql[6]+1710.0*coeff[1]*qc[6]-381.051177665153*coeff[1]*qr[4]+381.051177665153*coeff[1]*ql[4])*Jfac; 
  out[7] += -0.015625*(9.0*coeff[1]*qr[7]+9.0*coeff[1]*ql[7]+402.0*coeff[1]*qc[7]+19.36491673103708*coeff[1]*qr[3]-19.36491673103708*coeff[1]*ql[3]-26.83281572999747*coeff[1]*qr[1]-26.83281572999747*coeff[1]*ql[1]+53.66563145999495*coeff[1]*qc[1])*Jfac; 

  return 0.;

}

