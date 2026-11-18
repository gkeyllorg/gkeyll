#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.0078125*(87.30979326513149*coeff[0]*qr[8]-87.30979326513149*coeff[0]*ql[8]-199.0100499974813*coeff[0]*qr[4]-199.0100499974813*coeff[0]*ql[4]+398.0200999949626*coeff[0]*qc[4]+247.68326548234944*coeff[0]*qr[1]-247.68326548234944*coeff[0]*ql[1]-175.0*coeff[0]*qr[0]-175.0*coeff[0]*ql[0]+350.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0011160714285714285*(765.2901410576253*coeff[0]*qr[8]+765.2901410576253*coeff[0]*ql[8]+2703.7196600239454*coeff[0]*qc[8]-1870.6509562181823*coeff[0]*qr[4]+1870.6509562181823*coeff[0]*ql[4]+2415.0*coeff[0]*qr[1]+2415.0*coeff[0]*ql[1]+7182.0*coeff[0]*qc[1]-1733.782858376446*coeff[0]*qr[0]+1733.782858376446*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0026041666666666665*(261.9293797953944*coeff[0]*qr[10]-261.9293797953944*coeff[0]*ql[10]-597.0301499924439*coeff[0]*qr[6]-597.0301499924439*coeff[0]*ql[6]+1194.0602999848877*coeff[0]*qc[6]+743.0497964470483*coeff[0]*qr[3]-743.0497964470483*coeff[0]*ql[3]-525.0*coeff[0]*qr[2]-525.0*coeff[0]*ql[2]+1050.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0011160714285714285*(765.2901410576253*coeff[0]*qr[10]+765.2901410576253*coeff[0]*ql[10]+2703.7196600239454*coeff[0]*qc[10]-1870.6509562181825*coeff[0]*qr[6]+1870.6509562181825*coeff[0]*ql[6]+2415.0*coeff[0]*qr[3]+2415.0*coeff[0]*ql[3]+7182.0*coeff[0]*qc[3]-1733.782858376446*coeff[0]*qr[2]+1733.782858376446*coeff[0]*ql[2])*Jfac); 
  out[4] += -(0.0011160714285714285*(230.72711154088506*coeff[0]*qr[8]-230.72711154088506*coeff[0]*ql[8]-1015.0*coeff[0]*qr[4]-1015.0*coeff[0]*ql[4]+10430.0*coeff[0]*qc[4]+1599.542121983663*coeff[0]*qr[1]-1599.542121983663*coeff[0]*ql[1]-1236.5455915573837*coeff[0]*qr[0]-1236.5455915573837*coeff[0]*ql[0]+2473.0911831147673*coeff[0]*qc[0])*Jfac); 
  out[5] += -(0.0015625*(1238.4163274117475*coeff[0]*qr[7]-1238.4163274117475*coeff[0]*ql[7]-875.0*coeff[0]*qr[5]-875.0*coeff[0]*ql[5]+1750.0*coeff[0]*qc[5])*Jfac); 
  out[6] += -(3.720238095238095e-4*(692.1813346226551*coeff[0]*qr[10]-692.1813346226551*coeff[0]*ql[10]-3045.0*coeff[0]*qr[6]-3045.0*coeff[0]*ql[6]+31290.0*coeff[0]*qc[6]+4798.62636595099*coeff[0]*qr[3]-4798.62636595099*coeff[0]*ql[3]-3709.636774672151*coeff[0]*qr[2]-3709.636774672151*coeff[0]*ql[2]+7419.273549344302*coeff[0]*qc[2])*Jfac); 
  out[7] += -(0.0015625*(1725.0*coeff[0]*qr[7]+1725.0*coeff[0]*ql[7]+5130.0*coeff[0]*qc[7]-1238.4163274117475*coeff[0]*qr[5]+1238.4163274117475*coeff[0]*ql[5])*Jfac); 
  out[8] += 0.0078125*(153.0*coeff[0]*qr[8]+153.0*coeff[0]*ql[8]-1230.0*coeff[0]*qc[8]-183.3984732760881*coeff[0]*qr[4]+183.3984732760881*coeff[0]*ql[4]+114.564392373896*coeff[0]*qr[1]+114.564392373896*coeff[0]*ql[1]+82.4863625092051*coeff[0]*qc[1]-44.977772288098045*coeff[0]*qr[0]+44.977772288098045*coeff[0]*ql[0])*Jfac; 
  out[9] += -(0.0011160714285714285*(1733.782858376446*coeff[0]*qr[11]-1733.782858376446*coeff[0]*ql[11]-1225.0*coeff[0]*qr[9]-1225.0*coeff[0]*ql[9]+2450.0*coeff[0]*qc[9])*Jfac); 
  out[10] += 0.0026041666666666665*(459.0*coeff[0]*qr[10]+459.0*coeff[0]*ql[10]-3690.0*coeff[0]*qc[10]-550.1954198282642*coeff[0]*qr[6]+550.1954198282642*coeff[0]*ql[6]+343.693177121688*coeff[0]*qr[3]+343.693177121688*coeff[0]*ql[3]+247.45908752761534*coeff[0]*qc[3]-134.9333168642941*coeff[0]*qr[2]+134.9333168642941*coeff[0]*ql[2])*Jfac; 
  out[11] += -(0.0011160714285714285*(2415.0*coeff[0]*qr[11]+2415.0*coeff[0]*ql[11]+7182.0*coeff[0]*qc[11]-1733.782858376446*coeff[0]*qr[9]+1733.782858376446*coeff[0]*ql[9])*Jfac); 

  return 0.;

}

