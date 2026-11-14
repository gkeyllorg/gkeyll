#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += -(0.0078125*(27502.584878516423*coeff[0]*qr[3]-27502.584878516423*coeff[0]*ql[3]-40148.60053600873*coeff[0]*qr[2]-40148.60053600873*coeff[0]*ql[2]+80297.20107201746*coeff[0]*qc[2]+34372.54827620435*coeff[0]*qr[1]-34372.54827620435*coeff[0]*ql[1]-19845.0*coeff[0]*qr[0]-19845.0*coeff[0]*ql[0]+39690.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0078125*(41036.96534832954*coeff[0]*qr[3]+41036.96534832954*coeff[0]*ql[3]+108469.5666996047*coeff[0]*qc[3]-64659.45696493282*coeff[0]*qr[2]+64659.45696493282*coeff[0]*ql[2]+58275.0*coeff[0]*qr[1]+58275.0*coeff[0]*ql[1]+121590.0*coeff[0]*qc[1]-34372.54827620435*coeff[0]*qr[0]+34372.54827620435*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.00390625*(45524.23393095156*coeff[0]*qr[3]-45524.23393095156*coeff[0]*ql[3]-96075.0*coeff[0]*qr[2]-96075.0*coeff[0]*ql[2]+343350.0*coeff[0]*qc[2]+101259.1495865929*coeff[0]*qr[1]-101259.1495865929*coeff[0]*ql[1]-64096.88857503149*coeff[0]*qr[0]-64096.88857503149*coeff[0]*ql[0]+128193.77715006298*coeff[0]*qc[0])*Jfac); 
  out[3] += 0.00390625*(102555.0*coeff[0]*qr[3]+102555.0*coeff[0]*ql[3]-163530.0*coeff[0]*qc[3]-116339.70893465397*coeff[0]*qr[2]+116339.70893465397*coeff[0]*ql[2]+81455.28297784003*coeff[0]*qr[1]+81455.28297784003*coeff[0]*ql[1]+120017.65745089341*coeff[0]*qc[1]-40837.17148628197*coeff[0]*qr[0]+40837.17148628197*coeff[0]*ql[0])*Jfac; 

  return 0.;

}

