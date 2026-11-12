#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfy_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],6.);

  double vol_incr[12] = {0.0}; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(214.86394436340956*coeff[1]*fSkin[9])-214.86394436340956*coeff[1]*fEdge[9]-313.6609416875682*coeff[1]*fSkin[5]+313.6609416875682*coeff[1]*fEdge[5]-268.53553340784646*coeff[1]*fSkin[2]-268.53553340784646*coeff[1]*fEdge[2]-155.0390625*fSkin[0]*coeff[1]+155.0390625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(214.86394436340944*coeff[1]*fSkin[11])-214.86394436340944*coeff[1]*fEdge[11]-313.6609416875681*coeff[1]*fSkin[7]+313.6609416875681*coeff[1]*fEdge[7]-268.53553340784646*coeff[1]*fSkin[3]-268.53553340784646*coeff[1]*fEdge[3]-155.0390625*coeff[1]*fSkin[1]+155.0390625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(423.7092449203309*coeff[1]*fSkin[9])-320.60129178382454*coeff[1]*fEdge[9]-581.4013671669961*coeff[1]*fSkin[5]+505.15200753853765*coeff[1]*fEdge[5]-474.9609375*coeff[1]*fSkin[2]-455.2734375*coeff[1]*fEdge[2]-268.53553340784646*fSkin[0]*coeff[1]+268.53553340784646*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(423.70924492033095*coeff[1]*fSkin[11])-320.6012917838246*coeff[1]*fEdge[11]-581.4013671669962*coeff[1]*fSkin[7]+505.15200753853776*coeff[1]*fEdge[7]-474.9609375*coeff[1]*fSkin[3]-455.2734375*coeff[1]*fEdge[3]-268.53553340784646*coeff[1]*fSkin[1]+268.53553340784646*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = -(268.5355334078467*coeff[1]*fSkin[6])-268.5355334078467*coeff[1]*fEdge[6]-155.0390625*coeff[1]*fSkin[4]+155.0390625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(577.1644241520037*coeff[1]*fSkin[9])-177.82903879277956*coeff[1]*fEdge[9]-670.60546875*coeff[1]*fSkin[5]+375.29296875*coeff[1]*fEdge[5]-471.79291270108683*coeff[1]*fSkin[2]-395.54355307262836*coeff[1]*fEdge[2]-250.37847099621678*fSkin[0]*coeff[1]+250.37847099621678*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(474.9609375*coeff[1]*fSkin[6])-455.2734375*coeff[1]*fEdge[6]-268.5355334078467*coeff[1]*fSkin[4]+268.5355334078467*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = -(577.1644241520036*coeff[1]*fSkin[11])-177.8290387927795*coeff[1]*fEdge[11]-670.60546875*coeff[1]*fSkin[7]+375.29296875*coeff[1]*fEdge[7]-471.7929127010872*coeff[1]*fSkin[3]-395.5435530726286*coeff[1]*fEdge[3]-250.37847099621658*coeff[1]*fSkin[1]+250.37847099621658*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = -(268.53553340784663*coeff[1]*fSkin[10])-268.53553340784663*coeff[1]*fEdge[10]-155.0390625*coeff[1]*fSkin[8]+155.0390625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(319.39453125*coeff[1]*fSkin[9])+400.60546875*coeff[1]*fEdge[9]-3.119807698118848*coeff[1]*fSkin[5]-454.45198802599225*coeff[1]*fEdge[5]+234.40948720877657*coeff[1]*fSkin[2]+318.18469913218803*coeff[1]*fEdge[2]+159.52020111828895*fSkin[0]*coeff[1]-159.52020111828895*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = -(474.9609375*coeff[1]*fSkin[10])-455.2734375*coeff[1]*fEdge[10]-268.53553340784663*coeff[1]*fSkin[8]+268.53553340784663*coeff[1]*fEdge[8]; 
  edgeSurf_incr[11] = -(319.39453125*coeff[1]*fSkin[11])+400.60546875*coeff[1]*fEdge[11]-3.119807698118848*coeff[1]*fSkin[7]-454.4519880259921*coeff[1]*fEdge[7]+234.40948720877623*coeff[1]*fSkin[3]+318.1846991321877*coeff[1]*fEdge[3]+159.52020111828892*coeff[1]*fSkin[1]-159.52020111828892*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = -(103.10795313650637*coeff[1]*fSkin[9])+76.2493596284585*coeff[1]*fSkin[5]-19.6875*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = -(103.1079531365064*coeff[1]*fSkin[11])+76.24935962845854*coeff[1]*fSkin[7]-19.6875*coeff[1]*fSkin[3]; 
  boundSurf_incr[5] = 399.3353853592242*coeff[1]*fSkin[9]-295.3125*coeff[1]*fSkin[5]+76.2493596284585*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = -(19.6875*coeff[1]*fSkin[6]); 
  boundSurf_incr[7] = 399.3353853592241*coeff[1]*fSkin[11]-295.3125*coeff[1]*fSkin[7]+76.24935962845854*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = -(720.0*coeff[1]*fSkin[9])+457.5717957241111*coeff[1]*fSkin[5]-83.77521192341143*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = -(19.6875*coeff[1]*fSkin[10]); 
  boundSurf_incr[11] = -(720.0*coeff[1]*fSkin[11])+457.57179572411087*coeff[1]*fSkin[7]-83.77521192341143*coeff[1]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = 214.86394436340956*coeff[1]*fSkin[9]+214.86394436340956*coeff[1]*fEdge[9]-313.6609416875682*coeff[1]*fSkin[5]+313.6609416875682*coeff[1]*fEdge[5]+268.53553340784646*coeff[1]*fSkin[2]+268.53553340784646*coeff[1]*fEdge[2]-155.0390625*fSkin[0]*coeff[1]+155.0390625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 214.86394436340944*coeff[1]*fSkin[11]+214.86394436340944*coeff[1]*fEdge[11]-313.6609416875681*coeff[1]*fSkin[7]+313.6609416875681*coeff[1]*fEdge[7]+268.53553340784646*coeff[1]*fSkin[3]+268.53553340784646*coeff[1]*fEdge[3]-155.0390625*coeff[1]*fSkin[1]+155.0390625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(423.7092449203309*coeff[1]*fSkin[9])-320.60129178382454*coeff[1]*fEdge[9]+581.4013671669961*coeff[1]*fSkin[5]-505.15200753853765*coeff[1]*fEdge[5]-474.9609375*coeff[1]*fSkin[2]-455.2734375*coeff[1]*fEdge[2]+268.53553340784646*fSkin[0]*coeff[1]-268.53553340784646*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(423.70924492033095*coeff[1]*fSkin[11])-320.6012917838246*coeff[1]*fEdge[11]+581.4013671669962*coeff[1]*fSkin[7]-505.15200753853776*coeff[1]*fEdge[7]-474.9609375*coeff[1]*fSkin[3]-455.2734375*coeff[1]*fEdge[3]+268.53553340784646*coeff[1]*fSkin[1]-268.53553340784646*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = 268.5355334078467*coeff[1]*fSkin[6]+268.5355334078467*coeff[1]*fEdge[6]-155.0390625*coeff[1]*fSkin[4]+155.0390625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 577.1644241520037*coeff[1]*fSkin[9]+177.82903879277956*coeff[1]*fEdge[9]-670.60546875*coeff[1]*fSkin[5]+375.29296875*coeff[1]*fEdge[5]+471.79291270108683*coeff[1]*fSkin[2]+395.54355307262836*coeff[1]*fEdge[2]-250.37847099621678*fSkin[0]*coeff[1]+250.37847099621678*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(474.9609375*coeff[1]*fSkin[6])-455.2734375*coeff[1]*fEdge[6]+268.5355334078467*coeff[1]*fSkin[4]-268.5355334078467*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = 577.1644241520036*coeff[1]*fSkin[11]+177.8290387927795*coeff[1]*fEdge[11]-670.60546875*coeff[1]*fSkin[7]+375.29296875*coeff[1]*fEdge[7]+471.7929127010872*coeff[1]*fSkin[3]+395.5435530726286*coeff[1]*fEdge[3]-250.37847099621658*coeff[1]*fSkin[1]+250.37847099621658*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = 268.53553340784663*coeff[1]*fSkin[10]+268.53553340784663*coeff[1]*fEdge[10]-155.0390625*coeff[1]*fSkin[8]+155.0390625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(319.39453125*coeff[1]*fSkin[9])+400.60546875*coeff[1]*fEdge[9]+3.119807698118848*coeff[1]*fSkin[5]+454.45198802599225*coeff[1]*fEdge[5]+234.40948720877657*coeff[1]*fSkin[2]+318.18469913218803*coeff[1]*fEdge[2]-159.52020111828895*fSkin[0]*coeff[1]+159.52020111828895*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = -(474.9609375*coeff[1]*fSkin[10])-455.2734375*coeff[1]*fEdge[10]+268.53553340784663*coeff[1]*fSkin[8]-268.53553340784663*coeff[1]*fEdge[8]; 
  edgeSurf_incr[11] = -(319.39453125*coeff[1]*fSkin[11])+400.60546875*coeff[1]*fEdge[11]+3.119807698118848*coeff[1]*fSkin[7]+454.4519880259921*coeff[1]*fEdge[7]+234.40948720877623*coeff[1]*fSkin[3]+318.1846991321877*coeff[1]*fEdge[3]-159.52020111828892*coeff[1]*fSkin[1]+159.52020111828892*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = -(103.10795313650637*coeff[1]*fSkin[9])-76.2493596284585*coeff[1]*fSkin[5]-19.6875*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = -(103.1079531365064*coeff[1]*fSkin[11])-76.24935962845854*coeff[1]*fSkin[7]-19.6875*coeff[1]*fSkin[3]; 
  boundSurf_incr[5] = -(399.3353853592242*coeff[1]*fSkin[9])-295.3125*coeff[1]*fSkin[5]-76.2493596284585*coeff[1]*fSkin[2]; 
  boundSurf_incr[6] = -(19.6875*coeff[1]*fSkin[6]); 
  boundSurf_incr[7] = -(399.3353853592241*coeff[1]*fSkin[11])-295.3125*coeff[1]*fSkin[7]-76.24935962845854*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = -(720.0*coeff[1]*fSkin[9])-457.5717957241111*coeff[1]*fSkin[5]-83.77521192341143*coeff[1]*fSkin[2]; 
  boundSurf_incr[10] = -(19.6875*coeff[1]*fSkin[10]); 
  boundSurf_incr[11] = -(720.0*coeff[1]*fSkin[11])-457.57179572411087*coeff[1]*fSkin[7]-83.77521192341143*coeff[1]*fSkin[3]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 

  return 0.;
}

