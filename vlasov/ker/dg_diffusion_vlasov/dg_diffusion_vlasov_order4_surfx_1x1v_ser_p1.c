#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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
  out[2] += -(0.0625*(25.980762113533157*coeff[0]*qr[3]-25.980762113533157*coeff[0]*ql[3]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0625*(33.0*coeff[0]*qr[3]+33.0*coeff[0]*ql[3]+114.0*coeff[0]*qc[3]-25.980762113533157*coeff[0]*qr[2]+25.980762113533157*coeff[0]*ql[2])*Jfac); 
  out[4] += -(0.0625*(25.98076211353316*coeff[0]*qr[5]-25.98076211353316*coeff[0]*ql[5]-15.0*coeff[0]*qr[4]-15.0*coeff[0]*ql[4]+30.0*coeff[0]*qc[4])*Jfac); 
  out[5] += -(0.0625*(33.0*coeff[0]*qr[5]+33.0*coeff[0]*ql[5]+114.0*coeff[0]*qc[5]-25.98076211353316*coeff[0]*qr[4]+25.98076211353316*coeff[0]*ql[4])*Jfac); 

  return 0.;

}

