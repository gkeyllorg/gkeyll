#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += 0.0125*(53.66563145999496*coeff[0]*qr[2]+53.66563145999496*coeff[0]*ql[2]-107.3312629199899*coeff[0]*qc[2]-95.26279441628824*coeff[0]*qr[1]+95.26279441628824*coeff[0]*ql[1]+75.0*coeff[0]*qr[0]+75.0*coeff[0]*ql[0]-150.0*coeff[0]*qc[0])*Jfac; 
  out[1] += 0.003125*(236.2519841186524*coeff[0]*qr[2]-236.2519841186524*coeff[0]*ql[2]-465.0*coeff[0]*qr[1]-465.0*coeff[0]*ql[1]-1710.0*coeff[0]*qc[1]+381.051177665153*coeff[0]*qr[0]-381.051177665153*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.015625*(9.0*coeff[0]*qr[2]+9.0*coeff[0]*ql[2]+402.0*coeff[0]*qc[2]+19.36491673103709*coeff[0]*qr[1]-19.36491673103709*coeff[0]*ql[1]-26.83281572999748*coeff[0]*qr[0]-26.83281572999748*coeff[0]*ql[0]+53.66563145999496*coeff[0]*qc[0])*Jfac; 

  return 0.;

}

