#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[16] = {0.0}; 

  double edgeSurf_incr[16] = {0.0}; 
  double boundSurf_incr[16] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(214.86394436340956*coeff[0]*fSkin[8])-214.86394436340956*coeff[0]*fEdge[8]-313.6609416875682*coeff[0]*fSkin[4]+313.6609416875682*coeff[0]*fEdge[4]-268.53553340784646*coeff[0]*fSkin[1]-268.53553340784646*coeff[0]*fEdge[1]-155.0390625*coeff[0]*fSkin[0]+155.0390625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(423.7092449203309*coeff[0]*fSkin[8])-320.60129178382454*coeff[0]*fEdge[8]-581.4013671669961*coeff[0]*fSkin[4]+505.15200753853765*coeff[0]*fEdge[4]-474.9609375*coeff[0]*fSkin[1]-455.2734375*coeff[0]*fEdge[1]-268.53553340784646*coeff[0]*fSkin[0]+268.53553340784646*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(214.86394436340944*coeff[0]*fSkin[11])-214.86394436340944*coeff[0]*fEdge[11]-313.6609416875681*coeff[0]*fSkin[6]+313.6609416875681*coeff[0]*fEdge[6]-268.53553340784646*coeff[0]*fSkin[3]-268.53553340784646*coeff[0]*fEdge[3]-155.0390625*coeff[0]*fSkin[2]+155.0390625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(423.70924492033095*coeff[0]*fSkin[11])-320.6012917838246*coeff[0]*fEdge[11]-581.4013671669962*coeff[0]*fSkin[6]+505.15200753853776*coeff[0]*fEdge[6]-474.9609375*coeff[0]*fSkin[3]-455.2734375*coeff[0]*fEdge[3]-268.53553340784646*coeff[0]*fSkin[2]+268.53553340784646*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(577.1644241520037*coeff[0]*fSkin[8])-177.82903879277956*coeff[0]*fEdge[8]-670.60546875*coeff[0]*fSkin[4]+375.29296875*coeff[0]*fEdge[4]-471.79291270108683*coeff[0]*fSkin[1]-395.54355307262836*coeff[0]*fEdge[1]-250.37847099621678*coeff[0]*fSkin[0]+250.37847099621678*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(214.86394436340953*coeff[0]*fSkin[13])-214.86394436340953*coeff[0]*fEdge[13]-313.6609416875682*coeff[0]*fSkin[10]+313.6609416875682*coeff[0]*fEdge[10]-268.5355334078467*coeff[0]*fSkin[7]-268.5355334078467*coeff[0]*fEdge[7]-155.0390625*coeff[0]*fSkin[5]+155.0390625*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(577.1644241520036*coeff[0]*fSkin[11])-177.8290387927795*coeff[0]*fEdge[11]-670.60546875*coeff[0]*fSkin[6]+375.29296875*coeff[0]*fEdge[6]-471.7929127010872*coeff[0]*fSkin[3]-395.5435530726286*coeff[0]*fEdge[3]-250.37847099621658*coeff[0]*fSkin[2]+250.37847099621658*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(423.709244920331*coeff[0]*fSkin[13])-320.6012917838246*coeff[0]*fEdge[13]-581.4013671669962*coeff[0]*fSkin[10]+505.15200753853776*coeff[0]*fEdge[10]-474.9609375*coeff[0]*fSkin[7]-455.2734375*coeff[0]*fEdge[7]-268.5355334078467*coeff[0]*fSkin[5]+268.5355334078467*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(319.39453125*coeff[0]*fSkin[8])+400.60546875*coeff[0]*fEdge[8]-3.119807698118848*coeff[0]*fSkin[4]-454.45198802599225*coeff[0]*fEdge[4]+234.40948720877657*coeff[0]*fSkin[1]+318.18469913218803*coeff[0]*fEdge[1]+159.52020111828895*coeff[0]*fSkin[0]-159.52020111828895*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = -(214.86394436340956*coeff[0]*fSkin[15])-214.86394436340956*coeff[0]*fEdge[15]-313.6609416875682*coeff[0]*fSkin[14]+313.6609416875682*coeff[0]*fEdge[14]-268.53553340784663*coeff[0]*fSkin[12]-268.53553340784663*coeff[0]*fEdge[12]-155.0390625*coeff[0]*fSkin[9]+155.0390625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(577.1644241520036*coeff[0]*fSkin[13])-177.8290387927795*coeff[0]*fEdge[13]-670.60546875*coeff[0]*fSkin[10]+375.29296875*coeff[0]*fEdge[10]-471.7929127010872*coeff[0]*fSkin[7]-395.5435530726286*coeff[0]*fEdge[7]-250.37847099621678*coeff[0]*fSkin[5]+250.37847099621678*coeff[0]*fEdge[5]; 
  edgeSurf_incr[11] = -(319.39453125*coeff[0]*fSkin[11])+400.60546875*coeff[0]*fEdge[11]-3.119807698118848*coeff[0]*fSkin[6]-454.4519880259921*coeff[0]*fEdge[6]+234.40948720877623*coeff[0]*fSkin[3]+318.1846991321877*coeff[0]*fEdge[3]+159.52020111828892*coeff[0]*fSkin[2]-159.52020111828892*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = -(423.70924492033095*coeff[0]*fSkin[15])-320.6012917838246*coeff[0]*fEdge[15]-581.4013671669962*coeff[0]*fSkin[14]+505.1520075385377*coeff[0]*fEdge[14]-474.9609375*coeff[0]*fSkin[12]-455.2734375*coeff[0]*fEdge[12]-268.53553340784663*coeff[0]*fSkin[9]+268.53553340784663*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = -(319.39453125*coeff[0]*fSkin[13])+400.60546875*coeff[0]*fEdge[13]-3.119807698118933*coeff[0]*fSkin[10]-454.4519880259919*coeff[0]*fEdge[10]+234.40948720877634*coeff[0]*fSkin[7]+318.1846991321877*coeff[0]*fEdge[7]+159.52020111828898*coeff[0]*fSkin[5]-159.52020111828898*coeff[0]*fEdge[5]; 
  edgeSurf_incr[14] = -(577.1644241520036*coeff[0]*fSkin[15])-177.8290387927795*coeff[0]*fEdge[15]-670.60546875*coeff[0]*fSkin[14]+375.29296875*coeff[0]*fEdge[14]-471.79291270108706*coeff[0]*fSkin[12]-395.5435530726286*coeff[0]*fEdge[12]-250.37847099621672*coeff[0]*fSkin[9]+250.37847099621672*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = -(319.39453125*coeff[0]*fSkin[15])+400.60546875*coeff[0]*fEdge[15]-3.119807698118933*coeff[0]*fSkin[14]-454.4519880259919*coeff[0]*fEdge[14]+234.40948720877623*coeff[0]*fSkin[12]+318.1846991321877*coeff[0]*fEdge[12]+159.52020111828895*coeff[0]*fSkin[9]-159.52020111828895*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = -(103.10795313650637*coeff[0]*fSkin[8])+76.2493596284585*coeff[0]*fSkin[4]-19.6875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = -(103.1079531365064*coeff[0]*fSkin[11])+76.24935962845854*coeff[0]*fSkin[6]-19.6875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = 399.3353853592242*coeff[0]*fSkin[8]-295.3125*coeff[0]*fSkin[4]+76.2493596284585*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = 399.3353853592241*coeff[0]*fSkin[11]-295.3125*coeff[0]*fSkin[6]+76.24935962845854*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(103.10795313650641*coeff[0]*fSkin[13])+76.24935962845854*coeff[0]*fSkin[10]-19.6875*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = -(720.0*coeff[0]*fSkin[8])+457.5717957241111*coeff[0]*fSkin[4]-83.77521192341143*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 399.33538535922406*coeff[0]*fSkin[13]-295.3125*coeff[0]*fSkin[10]+76.24935962845854*coeff[0]*fSkin[7]; 
  boundSurf_incr[11] = -(720.0*coeff[0]*fSkin[11])+457.57179572411087*coeff[0]*fSkin[6]-83.77521192341143*coeff[0]*fSkin[3]; 
  boundSurf_incr[12] = -(103.1079531365064*coeff[0]*fSkin[15])+76.24935962845852*coeff[0]*fSkin[14]-19.6875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = -(720.0*coeff[0]*fSkin[13])+457.5717957241109*coeff[0]*fSkin[10]-83.77521192341148*coeff[0]*fSkin[7]; 
  boundSurf_incr[14] = 399.33538535922406*coeff[0]*fSkin[15]-295.3125*coeff[0]*fSkin[14]+76.24935962845852*coeff[0]*fSkin[12]; 
  boundSurf_incr[15] = -(720.0*coeff[0]*fSkin[15])+457.5717957241109*coeff[0]*fSkin[14]-83.77521192341143*coeff[0]*fSkin[12]; 

  } else { 

  edgeSurf_incr[0] = 214.86394436340956*coeff[0]*fSkin[8]+214.86394436340956*coeff[0]*fEdge[8]-313.6609416875682*coeff[0]*fSkin[4]+313.6609416875682*coeff[0]*fEdge[4]+268.53553340784646*coeff[0]*fSkin[1]+268.53553340784646*coeff[0]*fEdge[1]-155.0390625*coeff[0]*fSkin[0]+155.0390625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(423.7092449203309*coeff[0]*fSkin[8])-320.60129178382454*coeff[0]*fEdge[8]+581.4013671669961*coeff[0]*fSkin[4]-505.15200753853765*coeff[0]*fEdge[4]-474.9609375*coeff[0]*fSkin[1]-455.2734375*coeff[0]*fEdge[1]+268.53553340784646*coeff[0]*fSkin[0]-268.53553340784646*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 214.86394436340944*coeff[0]*fSkin[11]+214.86394436340944*coeff[0]*fEdge[11]-313.6609416875681*coeff[0]*fSkin[6]+313.6609416875681*coeff[0]*fEdge[6]+268.53553340784646*coeff[0]*fSkin[3]+268.53553340784646*coeff[0]*fEdge[3]-155.0390625*coeff[0]*fSkin[2]+155.0390625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(423.70924492033095*coeff[0]*fSkin[11])-320.6012917838246*coeff[0]*fEdge[11]+581.4013671669962*coeff[0]*fSkin[6]-505.15200753853776*coeff[0]*fEdge[6]-474.9609375*coeff[0]*fSkin[3]-455.2734375*coeff[0]*fEdge[3]+268.53553340784646*coeff[0]*fSkin[2]-268.53553340784646*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 577.1644241520037*coeff[0]*fSkin[8]+177.82903879277956*coeff[0]*fEdge[8]-670.60546875*coeff[0]*fSkin[4]+375.29296875*coeff[0]*fEdge[4]+471.79291270108683*coeff[0]*fSkin[1]+395.54355307262836*coeff[0]*fEdge[1]-250.37847099621678*coeff[0]*fSkin[0]+250.37847099621678*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 214.86394436340953*coeff[0]*fSkin[13]+214.86394436340953*coeff[0]*fEdge[13]-313.6609416875682*coeff[0]*fSkin[10]+313.6609416875682*coeff[0]*fEdge[10]+268.5355334078467*coeff[0]*fSkin[7]+268.5355334078467*coeff[0]*fEdge[7]-155.0390625*coeff[0]*fSkin[5]+155.0390625*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 577.1644241520036*coeff[0]*fSkin[11]+177.8290387927795*coeff[0]*fEdge[11]-670.60546875*coeff[0]*fSkin[6]+375.29296875*coeff[0]*fEdge[6]+471.7929127010872*coeff[0]*fSkin[3]+395.5435530726286*coeff[0]*fEdge[3]-250.37847099621658*coeff[0]*fSkin[2]+250.37847099621658*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(423.709244920331*coeff[0]*fSkin[13])-320.6012917838246*coeff[0]*fEdge[13]+581.4013671669962*coeff[0]*fSkin[10]-505.15200753853776*coeff[0]*fEdge[10]-474.9609375*coeff[0]*fSkin[7]-455.2734375*coeff[0]*fEdge[7]+268.5355334078467*coeff[0]*fSkin[5]-268.5355334078467*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(319.39453125*coeff[0]*fSkin[8])+400.60546875*coeff[0]*fEdge[8]+3.119807698118848*coeff[0]*fSkin[4]+454.45198802599225*coeff[0]*fEdge[4]+234.40948720877657*coeff[0]*fSkin[1]+318.18469913218803*coeff[0]*fEdge[1]-159.52020111828895*coeff[0]*fSkin[0]+159.52020111828895*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = 214.86394436340956*coeff[0]*fSkin[15]+214.86394436340956*coeff[0]*fEdge[15]-313.6609416875682*coeff[0]*fSkin[14]+313.6609416875682*coeff[0]*fEdge[14]+268.53553340784663*coeff[0]*fSkin[12]+268.53553340784663*coeff[0]*fEdge[12]-155.0390625*coeff[0]*fSkin[9]+155.0390625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 577.1644241520036*coeff[0]*fSkin[13]+177.8290387927795*coeff[0]*fEdge[13]-670.60546875*coeff[0]*fSkin[10]+375.29296875*coeff[0]*fEdge[10]+471.7929127010872*coeff[0]*fSkin[7]+395.5435530726286*coeff[0]*fEdge[7]-250.37847099621678*coeff[0]*fSkin[5]+250.37847099621678*coeff[0]*fEdge[5]; 
  edgeSurf_incr[11] = -(319.39453125*coeff[0]*fSkin[11])+400.60546875*coeff[0]*fEdge[11]+3.119807698118848*coeff[0]*fSkin[6]+454.4519880259921*coeff[0]*fEdge[6]+234.40948720877623*coeff[0]*fSkin[3]+318.1846991321877*coeff[0]*fEdge[3]-159.52020111828892*coeff[0]*fSkin[2]+159.52020111828892*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = -(423.70924492033095*coeff[0]*fSkin[15])-320.6012917838246*coeff[0]*fEdge[15]+581.4013671669962*coeff[0]*fSkin[14]-505.1520075385377*coeff[0]*fEdge[14]-474.9609375*coeff[0]*fSkin[12]-455.2734375*coeff[0]*fEdge[12]+268.53553340784663*coeff[0]*fSkin[9]-268.53553340784663*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = -(319.39453125*coeff[0]*fSkin[13])+400.60546875*coeff[0]*fEdge[13]+3.119807698118933*coeff[0]*fSkin[10]+454.4519880259919*coeff[0]*fEdge[10]+234.40948720877634*coeff[0]*fSkin[7]+318.1846991321877*coeff[0]*fEdge[7]-159.52020111828898*coeff[0]*fSkin[5]+159.52020111828898*coeff[0]*fEdge[5]; 
  edgeSurf_incr[14] = 577.1644241520036*coeff[0]*fSkin[15]+177.8290387927795*coeff[0]*fEdge[15]-670.60546875*coeff[0]*fSkin[14]+375.29296875*coeff[0]*fEdge[14]+471.79291270108706*coeff[0]*fSkin[12]+395.5435530726286*coeff[0]*fEdge[12]-250.37847099621672*coeff[0]*fSkin[9]+250.37847099621672*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = -(319.39453125*coeff[0]*fSkin[15])+400.60546875*coeff[0]*fEdge[15]+3.119807698118933*coeff[0]*fSkin[14]+454.4519880259919*coeff[0]*fEdge[14]+234.40948720877623*coeff[0]*fSkin[12]+318.1846991321877*coeff[0]*fEdge[12]-159.52020111828895*coeff[0]*fSkin[9]+159.52020111828895*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = -(103.10795313650637*coeff[0]*fSkin[8])-76.2493596284585*coeff[0]*fSkin[4]-19.6875*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = -(103.1079531365064*coeff[0]*fSkin[11])-76.24935962845854*coeff[0]*fSkin[6]-19.6875*coeff[0]*fSkin[3]; 
  boundSurf_incr[4] = -(399.3353853592242*coeff[0]*fSkin[8])-295.3125*coeff[0]*fSkin[4]-76.2493596284585*coeff[0]*fSkin[1]; 
  boundSurf_incr[6] = -(399.3353853592241*coeff[0]*fSkin[11])-295.3125*coeff[0]*fSkin[6]-76.24935962845854*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(103.10795313650641*coeff[0]*fSkin[13])-76.24935962845854*coeff[0]*fSkin[10]-19.6875*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = -(720.0*coeff[0]*fSkin[8])-457.5717957241111*coeff[0]*fSkin[4]-83.77521192341143*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = -(399.33538535922406*coeff[0]*fSkin[13])-295.3125*coeff[0]*fSkin[10]-76.24935962845854*coeff[0]*fSkin[7]; 
  boundSurf_incr[11] = -(720.0*coeff[0]*fSkin[11])-457.57179572411087*coeff[0]*fSkin[6]-83.77521192341143*coeff[0]*fSkin[3]; 
  boundSurf_incr[12] = -(103.1079531365064*coeff[0]*fSkin[15])-76.24935962845852*coeff[0]*fSkin[14]-19.6875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = -(720.0*coeff[0]*fSkin[13])-457.5717957241109*coeff[0]*fSkin[10]-83.77521192341148*coeff[0]*fSkin[7]; 
  boundSurf_incr[14] = -(399.33538535922406*coeff[0]*fSkin[15])-295.3125*coeff[0]*fSkin[14]-76.24935962845852*coeff[0]*fSkin[12]; 
  boundSurf_incr[15] = -(720.0*coeff[0]*fSkin[15])-457.5717957241109*coeff[0]*fSkin[14]-83.77521192341143*coeff[0]*fSkin[12]; 

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
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 

  return 0.;
}

