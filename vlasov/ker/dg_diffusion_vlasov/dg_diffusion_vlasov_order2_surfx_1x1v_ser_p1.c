#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.0625*(8.660254037844386*coeff[0]*qr[1]-8.660254037844386*coeff[0]*ql[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.0625*(7.0*coeff[0]*qr[1]+7.0*coeff[0]*ql[1]+46.0*coeff[0]*qc[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.0625*(8.660254037844386*coeff[0]*qr[3]-8.660254037844386*coeff[0]*ql[3]-9.0*coeff[0]*qr[2]-9.0*coeff[0]*ql[2]+18.0*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.0625*(7.0*coeff[0]*qr[3]+7.0*coeff[0]*ql[3]+46.0*coeff[0]*qc[3]-8.660254037844386*coeff[0]*qr[2]+8.660254037844386*coeff[0]*ql[2])*Jfac); 
  out[4] += -(0.0625*(8.660254037844387*coeff[0]*qr[5]-8.660254037844387*coeff[0]*ql[5]-9.0*coeff[0]*qr[4]-9.0*coeff[0]*ql[4]+18.0*coeff[0]*qc[4])*Jfac); 
  out[5] += -(0.0625*(7.0*coeff[0]*qr[5]+7.0*coeff[0]*ql[5]+46.0*coeff[0]*qc[5]-8.660254037844387*coeff[0]*qr[4]+8.660254037844387*coeff[0]*ql[4])*Jfac); 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_vlasov_order2_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -(0.03125*((21.213203435596427*coeff[1]+12.24744871391589*coeff[0])*qr[1]+(21.213203435596427*coeff[1]-12.24744871391589*coeff[0])*ql[1]+42.42640687119286*coeff[1]*qc[1]+(22.0454076850486*ql[0]-22.0454076850486*qr[0])*coeff[1]-12.727922061357857*coeff[0]*qr[0]-12.727922061357857*coeff[0]*ql[0]+25.455844122715718*coeff[0]*qc[0])*Jfac); 
  out[1] += -(0.03125*((17.146428199482244*coeff[1]+9.899494936611665*coeff[0])*qr[1]+(9.899494936611665*coeff[0]-17.146428199482244*coeff[1])*ql[1]+65.05382386916239*coeff[0]*qc[1]+(-(21.213203435596427*qr[0])-21.213203435596427*ql[0]+42.42640687119286*qc[0])*coeff[1]-12.24744871391589*coeff[0]*qr[0]+12.24744871391589*coeff[0]*ql[0])*Jfac); 
  out[2] += -(0.03125*((21.213203435596427*coeff[1]+12.24744871391589*coeff[0])*qr[3]+(21.213203435596427*coeff[1]-12.24744871391589*coeff[0])*ql[3]+42.42640687119286*coeff[1]*qc[3]+(-(22.0454076850486*coeff[1])-12.727922061357857*coeff[0])*qr[2]+(22.0454076850486*coeff[1]-12.727922061357857*coeff[0])*ql[2]+25.455844122715718*coeff[0]*qc[2])*Jfac); 
  out[3] += -(0.03125*((17.146428199482244*coeff[1]+9.899494936611665*coeff[0])*qr[3]+(9.899494936611665*coeff[0]-17.146428199482244*coeff[1])*ql[3]+65.05382386916239*coeff[0]*qc[3]+(-(21.213203435596427*coeff[1])-12.24744871391589*coeff[0])*qr[2]+(12.24744871391589*coeff[0]-21.213203435596427*coeff[1])*ql[2]+42.42640687119286*coeff[1]*qc[2])*Jfac); 
  out[4] += -(0.03125*((21.21320343559643*coeff[1]+12.247448713915892*coeff[0])*qr[5]+(21.21320343559643*coeff[1]-12.247448713915892*coeff[0])*ql[5]+42.42640687119286*coeff[1]*qc[5]+(-(22.0454076850486*coeff[1])-12.727922061357857*coeff[0])*qr[4]+(22.0454076850486*coeff[1]-12.727922061357857*coeff[0])*ql[4]+25.455844122715718*coeff[0]*qc[4])*Jfac); 
  out[5] += -(0.03125*((17.146428199482244*coeff[1]+9.899494936611665*coeff[0])*qr[5]+(9.899494936611665*coeff[0]-17.146428199482244*coeff[1])*ql[5]+65.05382386916239*coeff[0]*qc[5]+(-(21.21320343559643*coeff[1])-12.247448713915892*coeff[0])*qr[4]+(12.247448713915892*coeff[0]-21.21320343559643*coeff[1])*ql[4]+42.42640687119286*coeff[1]*qc[4])*Jfac); 

  return 0.;

}

