#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.0078125*(87.30979326513149*coeff[0]*qr[3]-87.30979326513149*coeff[0]*ql[3]-199.0100499974813*coeff[0]*qr[2]-199.0100499974813*coeff[0]*ql[2]+398.0200999949626*coeff[0]*qc[2]+247.68326548234944*coeff[0]*qr[1]-247.68326548234944*coeff[0]*ql[1]-175.0*coeff[0]*qr[0]-175.0*coeff[0]*ql[0]+350.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0011160714285714285*(765.2901410576253*coeff[0]*qr[3]+765.2901410576253*coeff[0]*ql[3]+2703.7196600239454*coeff[0]*qc[3]-1870.6509562181823*coeff[0]*qr[2]+1870.6509562181823*coeff[0]*ql[2]+2415.0*coeff[0]*qr[1]+2415.0*coeff[0]*ql[1]+7182.0*coeff[0]*qc[1]-1733.782858376446*coeff[0]*qr[0]+1733.782858376446*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0011160714285714285*(230.72711154088506*coeff[0]*qr[3]-230.72711154088506*coeff[0]*ql[3]-1015.0*coeff[0]*qr[2]-1015.0*coeff[0]*ql[2]+10430.0*coeff[0]*qc[2]+1599.542121983663*coeff[0]*qr[1]-1599.542121983663*coeff[0]*ql[1]-1236.5455915573837*coeff[0]*qr[0]-1236.5455915573837*coeff[0]*ql[0]+2473.0911831147673*coeff[0]*qc[0])*Jfac); 
  out[3] += 0.0078125*(153.0*coeff[0]*qr[3]+153.0*coeff[0]*ql[3]-1230.0*coeff[0]*qc[3]-183.3984732760881*coeff[0]*qr[2]+183.3984732760881*coeff[0]*ql[2]+114.564392373896*coeff[0]*qr[1]+114.564392373896*coeff[0]*ql[1]+82.4863625092051*coeff[0]*qc[1]-44.977772288098045*coeff[0]*qr[0]+44.977772288098045*coeff[0]*ql[0])*Jfac; 

  return 0.;

}

