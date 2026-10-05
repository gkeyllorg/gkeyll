#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  out[0] += -(0.0078125*(27502.584878516423*coeff[1]*qr[9]-27502.584878516423*coeff[1]*ql[9]-40148.60053600873*coeff[1]*qr[5]-40148.60053600873*coeff[1]*ql[5]+80297.20107201746*coeff[1]*qc[5]+34372.54827620435*coeff[1]*qr[2]-34372.54827620435*coeff[1]*ql[2]+(-(19845.0*qr[0])-19845.0*ql[0]+39690.0*qc[0])*coeff[1])*Jfac); 
  out[1] += -(0.0078125*(27502.58487851641*coeff[1]*qr[11]-27502.58487851641*coeff[1]*ql[11]-40148.600536008715*coeff[1]*qr[7]-40148.600536008715*coeff[1]*ql[7]+80297.20107201743*coeff[1]*qc[7]+34372.54827620435*coeff[1]*qr[3]-34372.54827620435*coeff[1]*ql[3]-19845.0*coeff[1]*qr[1]-19845.0*coeff[1]*ql[1]+39690.0*coeff[1]*qc[1])*Jfac); 
  out[2] += -(0.0078125*(41036.96534832954*coeff[1]*qr[9]+41036.96534832954*coeff[1]*ql[9]+108469.5666996047*coeff[1]*qc[9]-64659.45696493282*coeff[1]*qr[5]+64659.45696493282*coeff[1]*ql[5]+58275.0*coeff[1]*qr[2]+58275.0*coeff[1]*ql[2]+121590.0*coeff[1]*qc[2]+(34372.54827620435*ql[0]-34372.54827620435*qr[0])*coeff[1])*Jfac); 
  out[3] += -(0.0078125*(41036.96534832955*coeff[1]*qr[11]+41036.96534832955*coeff[1]*ql[11]+108469.56669960472*coeff[1]*qc[11]-64659.45696493283*coeff[1]*qr[7]+64659.45696493283*coeff[1]*ql[7]+58275.0*coeff[1]*qr[3]+58275.0*coeff[1]*ql[3]+121590.0*coeff[1]*qc[3]-34372.54827620435*coeff[1]*qr[1]+34372.54827620435*coeff[1]*ql[1])*Jfac); 
  out[4] += -(0.0078125*(34372.548276204376*coeff[1]*qr[6]-34372.548276204376*coeff[1]*ql[6]-19845.0*coeff[1]*qr[4]-19845.0*coeff[1]*ql[4]+39690.0*coeff[1]*qc[4])*Jfac); 
  out[5] += -(0.00390625*(45524.23393095156*coeff[1]*qr[9]-45524.23393095156*coeff[1]*ql[9]-96075.0*coeff[1]*qr[5]-96075.0*coeff[1]*ql[5]+343350.0*coeff[1]*qc[5]+101259.1495865929*coeff[1]*qr[2]-101259.1495865929*coeff[1]*ql[2]+(-(64096.88857503149*qr[0])-64096.88857503149*ql[0]+128193.77715006298*qc[0])*coeff[1])*Jfac); 
  out[6] += -(0.0078125*(58275.0*coeff[1]*qr[6]+58275.0*coeff[1]*ql[6]+121590.0*coeff[1]*qc[6]-34372.548276204376*coeff[1]*qr[4]+34372.548276204376*coeff[1]*ql[4])*Jfac); 
  out[7] += -(0.00390625*(45524.233930951545*coeff[1]*qr[11]-45524.233930951545*coeff[1]*ql[11]-96075.0*coeff[1]*qr[7]-96075.0*coeff[1]*ql[7]+343350.0*coeff[1]*qc[7]+101259.14958659293*coeff[1]*qr[3]-101259.14958659293*coeff[1]*ql[3]-64096.88857503146*coeff[1]*qr[1]-64096.88857503146*coeff[1]*ql[1]+128193.77715006292*coeff[1]*qc[1])*Jfac); 
  out[8] += -(0.0078125*(34372.54827620437*coeff[1]*qr[10]-34372.54827620437*coeff[1]*ql[10]-19845.0*coeff[1]*qr[8]-19845.0*coeff[1]*ql[8]+39690.0*coeff[1]*qc[8])*Jfac); 
  out[9] += 0.00390625*(102555.0*coeff[1]*qr[9]+102555.0*coeff[1]*ql[9]-163530.0*coeff[1]*qc[9]-116339.70893465397*coeff[1]*qr[5]+116339.70893465397*coeff[1]*ql[5]+81455.28297784003*coeff[1]*qr[2]+81455.28297784003*coeff[1]*ql[2]+120017.65745089341*coeff[1]*qc[2]+(40837.17148628197*ql[0]-40837.17148628197*qr[0])*coeff[1])*Jfac; 
  out[10] += -(0.0078125*(58275.0*coeff[1]*qr[10]+58275.0*coeff[1]*ql[10]+121590.0*coeff[1]*qc[10]-34372.54827620437*coeff[1]*qr[8]+34372.54827620437*coeff[1]*ql[8])*Jfac); 
  out[11] += 0.00390625*(102555.0*coeff[1]*qr[11]+102555.0*coeff[1]*ql[11]-163530.0*coeff[1]*qc[11]-116339.70893465396*coeff[1]*qr[7]+116339.70893465396*coeff[1]*ql[7]+81455.28297784005*coeff[1]*qr[3]+81455.28297784005*coeff[1]*ql[3]+120017.65745089344*coeff[1]*qc[3]-40837.17148628195*coeff[1]*qr[1]+40837.17148628195*coeff[1]*ql[1])*Jfac; 

  return 0.;

}

