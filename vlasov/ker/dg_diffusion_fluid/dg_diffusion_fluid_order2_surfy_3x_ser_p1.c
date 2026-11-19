#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += -0.0625*(8.660254037844386*coeff[1]*qr[2]-8.660254037844386*coeff[1]*ql[2]+((-9.0*qr[0])-9.0*ql[0]+18.0*qc[0])*coeff[1])*Jfac; 
  out[1] += -0.0625*(8.660254037844386*coeff[1]*qr[4]-8.660254037844386*coeff[1]*ql[4]-9.0*coeff[1]*qr[1]-9.0*coeff[1]*ql[1]+18.0*coeff[1]*qc[1])*Jfac; 
  out[2] += -0.0625*(7.0*coeff[1]*qr[2]+7.0*coeff[1]*ql[2]+46.0*coeff[1]*qc[2]+(8.660254037844386*ql[0]-8.660254037844386*qr[0])*coeff[1])*Jfac; 
  out[3] += -0.0625*(8.660254037844386*coeff[1]*qr[6]-8.660254037844386*coeff[1]*ql[6]-9.0*coeff[1]*qr[3]-9.0*coeff[1]*ql[3]+18.0*coeff[1]*qc[3])*Jfac; 
  out[4] += -0.0625*(7.0*coeff[1]*qr[4]+7.0*coeff[1]*ql[4]+46.0*coeff[1]*qc[4]-8.660254037844386*coeff[1]*qr[1]+8.660254037844386*coeff[1]*ql[1])*Jfac; 
  out[5] += -0.0625*(8.660254037844386*coeff[1]*qr[7]-8.660254037844386*coeff[1]*ql[7]-9.0*coeff[1]*qr[5]-9.0*coeff[1]*ql[5]+18.0*coeff[1]*qc[5])*Jfac; 
  out[6] += -0.0625*(7.0*coeff[1]*qr[6]+7.0*coeff[1]*ql[6]+46.0*coeff[1]*qc[6]-8.660254037844386*coeff[1]*qr[3]+8.660254037844386*coeff[1]*ql[3])*Jfac; 
  out[7] += -0.0625*(7.0*coeff[1]*qr[7]+7.0*coeff[1]*ql[7]+46.0*coeff[1]*qc[7]-8.660254037844386*coeff[1]*qr[5]+8.660254037844386*coeff[1]*ql[5])*Jfac; 

  return 0.;

}

