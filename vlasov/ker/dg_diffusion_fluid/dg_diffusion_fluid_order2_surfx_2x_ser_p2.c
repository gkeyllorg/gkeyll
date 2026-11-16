#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[0]*qr[4]+53.66563145999496*coeff[0]*ql[4]-107.3312629199899*coeff[0]*qc[4]-95.26279441628824*coeff[0]*qr[1]+95.26279441628824*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.003125*(236.2519841186524*coeff[0]*qr[4]-236.2519841186524*coeff[0]*ql[4]-465.0*coeff[0]*qr[1]-465.0*coeff[0]*ql[1]-1710.0*coeff[0]*qc[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += 0.0125*(53.66563145999495*coeff[0]*qr[6]+53.66563145999495*coeff[0]*ql[6]-107.3312629199899*coeff[0]*qc[6]-95.26279441628824*coeff[0]*qr[3]+95.26279441628824*coeff[0]*ql[3]+75.0*coeff[0]*qr[2]+75.0*coeff[0]*ql[2]-150.0*coeff[0]*qc[2])*Jfac; 
  out[3] += 0.003125*(236.2519841186524*coeff[0]*qr[6]-236.2519841186524*coeff[0]*ql[6]-465.0*coeff[0]*qr[3]-465.0*coeff[0]*ql[3]-1710.0*coeff[0]*qc[3]+381.051177665153*coeff[0]*qr[2]-381.051177665153*coeff[0]*ql[2])*Jfac; 
  out[4] += -0.015625*(9.0*coeff[0]*qr[4]+9.0*coeff[0]*ql[4]+402.0*coeff[0]*qc[4]+19.36491673103709*coeff[0]*qr[1]-19.36491673103709*coeff[0]*ql[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 
  out[5] += -0.0125*(95.26279441628826*coeff[0]*qr[7]-95.26279441628826*coeff[0]*ql[7]-75.0*coeff[0]*qr[5]-75.0*coeff[0]*ql[5]+150.0*coeff[0]*qc[5])*Jfac; 
  out[6] += -0.015625*(9.0*coeff[0]*qr[6]+9.0*coeff[0]*ql[6]+402.0*coeff[0]*qc[6]+19.36491673103708*coeff[0]*qr[3]-19.36491673103708*coeff[0]*ql[3]-26.83281572999747*coeff[0]*qr[2]-26.83281572999747*coeff[0]*ql[2]+53.66563145999495*coeff[0]*qc[2])*Jfac; 
  out[7] += -0.003125*(465.0*coeff[0]*qr[7]+465.0*coeff[0]*ql[7]+1710.0*coeff[0]*qc[7]-381.051177665153*coeff[0]*qr[5]+381.051177665153*coeff[0]*ql[5])*Jfac; 

  return 0.;

}

