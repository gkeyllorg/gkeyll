#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  out[0] += -(0.0078125*(87.30979326513149*coeff[1]*qr[9]-87.30979326513149*coeff[1]*ql[9]-199.0100499974813*coeff[1]*qr[5]-199.0100499974813*coeff[1]*ql[5]+398.0200999949626*coeff[1]*qc[5]+247.68326548234944*coeff[1]*qr[2]-247.68326548234944*coeff[1]*ql[2]+(-(175.0*qr[0])-175.0*ql[0]+350.0*qc[0])*coeff[1])*Jfac); 
  out[1] += -(0.0026041666666666665*(261.9293797953944*coeff[1]*qr[11]-261.9293797953944*coeff[1]*ql[11]-597.0301499924439*coeff[1]*qr[7]-597.0301499924439*coeff[1]*ql[7]+1194.0602999848877*coeff[1]*qc[7]+743.0497964470483*coeff[1]*qr[3]-743.0497964470483*coeff[1]*ql[3]-525.0*coeff[1]*qr[1]-525.0*coeff[1]*ql[1]+1050.0*coeff[1]*qc[1])*Jfac); 
  out[2] += -(0.0011160714285714285*(765.2901410576253*coeff[1]*qr[9]+765.2901410576253*coeff[1]*ql[9]+2703.7196600239454*coeff[1]*qc[9]-1870.6509562181823*coeff[1]*qr[5]+1870.6509562181823*coeff[1]*ql[5]+2415.0*coeff[1]*qr[2]+2415.0*coeff[1]*ql[2]+7182.0*coeff[1]*qc[2]+(1733.782858376446*ql[0]-1733.782858376446*qr[0])*coeff[1])*Jfac); 
  out[3] += -(0.0011160714285714285*(765.2901410576253*coeff[1]*qr[11]+765.2901410576253*coeff[1]*ql[11]+2703.7196600239454*coeff[1]*qc[11]-1870.6509562181825*coeff[1]*qr[7]+1870.6509562181825*coeff[1]*ql[7]+2415.0*coeff[1]*qr[3]+2415.0*coeff[1]*ql[3]+7182.0*coeff[1]*qc[3]-1733.782858376446*coeff[1]*qr[1]+1733.782858376446*coeff[1]*ql[1])*Jfac); 
  out[4] += -(0.0015625*(1238.4163274117475*coeff[1]*qr[6]-1238.4163274117475*coeff[1]*ql[6]-875.0*coeff[1]*qr[4]-875.0*coeff[1]*ql[4]+1750.0*coeff[1]*qc[4])*Jfac); 
  out[5] += -(0.0011160714285714285*(230.72711154088506*coeff[1]*qr[9]-230.72711154088506*coeff[1]*ql[9]-1015.0*coeff[1]*qr[5]-1015.0*coeff[1]*ql[5]+10430.0*coeff[1]*qc[5]+1599.542121983663*coeff[1]*qr[2]-1599.542121983663*coeff[1]*ql[2]+(-(1236.5455915573837*qr[0])-1236.5455915573837*ql[0]+2473.0911831147673*qc[0])*coeff[1])*Jfac); 
  out[6] += -(0.0015625*(1725.0*coeff[1]*qr[6]+1725.0*coeff[1]*ql[6]+5130.0*coeff[1]*qc[6]-1238.4163274117475*coeff[1]*qr[4]+1238.4163274117475*coeff[1]*ql[4])*Jfac); 
  out[7] += -(3.720238095238095e-4*(692.1813346226551*coeff[1]*qr[11]-692.1813346226551*coeff[1]*ql[11]-3045.0*coeff[1]*qr[7]-3045.0*coeff[1]*ql[7]+31290.0*coeff[1]*qc[7]+4798.62636595099*coeff[1]*qr[3]-4798.62636595099*coeff[1]*ql[3]-3709.636774672151*coeff[1]*qr[1]-3709.636774672151*coeff[1]*ql[1]+7419.273549344302*coeff[1]*qc[1])*Jfac); 
  out[8] += -(0.0011160714285714285*(1733.782858376446*coeff[1]*qr[10]-1733.782858376446*coeff[1]*ql[10]-1225.0*coeff[1]*qr[8]-1225.0*coeff[1]*ql[8]+2450.0*coeff[1]*qc[8])*Jfac); 
  out[9] += 0.0078125*(153.0*coeff[1]*qr[9]+153.0*coeff[1]*ql[9]-1230.0*coeff[1]*qc[9]-183.3984732760881*coeff[1]*qr[5]+183.3984732760881*coeff[1]*ql[5]+114.564392373896*coeff[1]*qr[2]+114.564392373896*coeff[1]*ql[2]+82.4863625092051*coeff[1]*qc[2]+(44.977772288098045*ql[0]-44.977772288098045*qr[0])*coeff[1])*Jfac; 
  out[10] += -(0.0011160714285714285*(2415.0*coeff[1]*qr[10]+2415.0*coeff[1]*ql[10]+7182.0*coeff[1]*qc[10]-1733.782858376446*coeff[1]*qr[8]+1733.782858376446*coeff[1]*ql[8])*Jfac); 
  out[11] += 0.0026041666666666665*(459.0*coeff[1]*qr[11]+459.0*coeff[1]*ql[11]-3690.0*coeff[1]*qc[11]-550.1954198282642*coeff[1]*qr[7]+550.1954198282642*coeff[1]*ql[7]+343.693177121688*coeff[1]*qr[3]+343.693177121688*coeff[1]*ql[3]+247.45908752761534*coeff[1]*qc[3]-134.9333168642941*coeff[1]*qr[1]+134.9333168642941*coeff[1]*ql[1])*Jfac; 

  return 0.;

}

