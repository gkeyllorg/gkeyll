#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  out[0] += -(0.0078125*(27502.584878516423*coeff[0]*qr[8]-27502.584878516423*coeff[0]*ql[8]-40148.60053600873*coeff[0]*qr[4]-40148.60053600873*coeff[0]*ql[4]+80297.20107201746*coeff[0]*qc[4]+34372.54827620435*coeff[0]*qr[1]-34372.54827620435*coeff[0]*ql[1]-19845.0*coeff[0]*qr[0]-19845.0*coeff[0]*ql[0]+39690.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0078125*(41036.96534832954*coeff[0]*qr[8]+41036.96534832954*coeff[0]*ql[8]+108469.5666996047*coeff[0]*qc[8]-64659.45696493282*coeff[0]*qr[4]+64659.45696493282*coeff[0]*ql[4]+58275.0*coeff[0]*qr[1]+58275.0*coeff[0]*ql[1]+121590.0*coeff[0]*qc[1]-34372.54827620435*coeff[0]*qr[0]+34372.54827620435*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0078125*(27502.58487851641*coeff[0]*qr[10]-27502.58487851641*coeff[0]*ql[10]-40148.600536008715*coeff[0]*qr[6]-40148.600536008715*coeff[0]*ql[6]+80297.20107201743*coeff[0]*qc[6]+34372.54827620435*coeff[0]*qr[3]-34372.54827620435*coeff[0]*ql[3]-19845.0*coeff[0]*qr[2]-19845.0*coeff[0]*ql[2]+39690.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0078125*(41036.96534832955*coeff[0]*qr[10]+41036.96534832955*coeff[0]*ql[10]+108469.56669960472*coeff[0]*qc[10]-64659.45696493283*coeff[0]*qr[6]+64659.45696493283*coeff[0]*ql[6]+58275.0*coeff[0]*qr[3]+58275.0*coeff[0]*ql[3]+121590.0*coeff[0]*qc[3]-34372.54827620435*coeff[0]*qr[2]+34372.54827620435*coeff[0]*ql[2])*Jfac); 
  out[4] += -(0.00390625*(45524.23393095156*coeff[0]*qr[8]-45524.23393095156*coeff[0]*ql[8]-96075.0*coeff[0]*qr[4]-96075.0*coeff[0]*ql[4]+343350.0*coeff[0]*qc[4]+101259.1495865929*coeff[0]*qr[1]-101259.1495865929*coeff[0]*ql[1]-64096.88857503149*coeff[0]*qr[0]-64096.88857503149*coeff[0]*ql[0]+128193.77715006298*coeff[0]*qc[0])*Jfac); 
  out[5] += -(0.0078125*(34372.548276204376*coeff[0]*qr[7]-34372.548276204376*coeff[0]*ql[7]-19845.0*coeff[0]*qr[5]-19845.0*coeff[0]*ql[5]+39690.0*coeff[0]*qc[5])*Jfac); 
  out[6] += -(0.00390625*(45524.233930951545*coeff[0]*qr[10]-45524.233930951545*coeff[0]*ql[10]-96075.0*coeff[0]*qr[6]-96075.0*coeff[0]*ql[6]+343350.0*coeff[0]*qc[6]+101259.14958659293*coeff[0]*qr[3]-101259.14958659293*coeff[0]*ql[3]-64096.88857503146*coeff[0]*qr[2]-64096.88857503146*coeff[0]*ql[2]+128193.77715006292*coeff[0]*qc[2])*Jfac); 
  out[7] += -(0.0078125*(58275.0*coeff[0]*qr[7]+58275.0*coeff[0]*ql[7]+121590.0*coeff[0]*qc[7]-34372.548276204376*coeff[0]*qr[5]+34372.548276204376*coeff[0]*ql[5])*Jfac); 
  out[8] += 0.00390625*(102555.0*coeff[0]*qr[8]+102555.0*coeff[0]*ql[8]-163530.0*coeff[0]*qc[8]-116339.70893465397*coeff[0]*qr[4]+116339.70893465397*coeff[0]*ql[4]+81455.28297784003*coeff[0]*qr[1]+81455.28297784003*coeff[0]*ql[1]+120017.65745089341*coeff[0]*qc[1]-40837.17148628197*coeff[0]*qr[0]+40837.17148628197*coeff[0]*ql[0])*Jfac; 
  out[9] += -(0.0078125*(34372.54827620437*coeff[0]*qr[11]-34372.54827620437*coeff[0]*ql[11]-19845.0*coeff[0]*qr[9]-19845.0*coeff[0]*ql[9]+39690.0*coeff[0]*qc[9])*Jfac); 
  out[10] += 0.00390625*(102555.0*coeff[0]*qr[10]+102555.0*coeff[0]*ql[10]-163530.0*coeff[0]*qc[10]-116339.70893465396*coeff[0]*qr[6]+116339.70893465396*coeff[0]*ql[6]+81455.28297784005*coeff[0]*qr[3]+81455.28297784005*coeff[0]*ql[3]+120017.65745089344*coeff[0]*qc[3]-40837.17148628195*coeff[0]*qr[2]+40837.17148628195*coeff[0]*ql[2])*Jfac; 
  out[11] += -(0.0078125*(58275.0*coeff[0]*qr[11]+58275.0*coeff[0]*ql[11]+121590.0*coeff[0]*qc[11]-34372.54827620437*coeff[0]*qr[9]+34372.54827620437*coeff[0]*ql[9])*Jfac); 

  return 0.;

}

