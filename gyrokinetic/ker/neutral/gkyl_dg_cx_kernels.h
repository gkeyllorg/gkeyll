#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
 
GKYL_CU_DH static inline 
double sigma_cx_1x_ser_p1(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a: constant in fitting function.
  // b: constant in fitting function.
  // maxwellian_moms_ion[6]: ion prim vars.
  // maxwellian_moms_neut[10]: neut prim vars.
  // u_ion[6]: ion drift velocity vector (upar_i b_1, upar_i b_2, upar_i b_3).
  // v_sigma_cx: cell ave cross section fitting eqn.

  double m0_neut_av = 0.7071067811865476*maxwellian_moms_neut[0]; 
 
  const double *vt_sq_ion = &maxwellian_moms_ion[4]; 
  const double *u_neut = &maxwellian_moms_neut[2]; 
  const double *vt_sq_neut = &maxwellian_moms_neut[8]; 
 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.5*pow(u_neut[4],2)-1.0*u_ion[4]*u_neut[4]+0.5*pow(u_ion[4],2)+0.5*pow(u_neut[2],2)-1.0*u_ion[2]*u_neut[2]+0.5*pow(u_ion[2],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
  double v_cx = sqrt(2.5464790894703255*vt_sq_neut_av+2.5464790894703255*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 1.4142135623730951*a*v_cx-1.4142135623730951*b*v_cx*log(v_cx); 
 
  return 0.2357022603955158*v_sigma_cx[0]*m0_neut_av; 
  }
} 
 
GKYL_CU_DH static inline 
double sigma_cx_1x_ser_p2(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a: constant in fitting function.
  // b: constant in fitting function.
  // maxwellian_moms_ion[9]: ion prim vars.
  // maxwellian_moms_neut[15]: neut prim vars.
  // u_ion[9]: ion drift velocity vector (upar_i b_1, upar_i b_2, upar_i b_3).
  // v_sigma_cx: cell ave cross section fitting eqn.

  double m0_neut_av = 0.7071067811865476*maxwellian_moms_neut[0]; 
 
  const double *vt_sq_ion = &maxwellian_moms_ion[6]; 
  const double *u_neut = &maxwellian_moms_neut[3]; 
  const double *vt_sq_neut = &maxwellian_moms_neut[12]; 
 
  double vt_sq_ion_av = 0.7071067811865476*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.7071067811865476*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.5*pow(u_neut[6],2)-1.0*u_ion[6]*u_neut[6]+0.5*pow(u_ion[6],2)+0.5*pow(u_neut[3],2)-1.0*u_ion[3]*u_neut[3]+0.5*pow(u_ion[3],2)+0.5*pow(u_neut[0],2)-1.0*u_ion[0]*u_neut[0]+0.5*pow(u_ion[0],2); 
 
  double v_cx = sqrt(2.5464790894703255*vt_sq_neut_av+2.5464790894703255*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 1.4142135623730951*a*v_cx-1.4142135623730951*b*v_cx*log(v_cx); 
 
  return 0.1414213562373095*v_sigma_cx[0]*m0_neut_av; 
  }
} 
 
GKYL_CU_DH static inline 
double sigma_cx_2x_ser_p1(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a: constant in fitting function.
  // b: constant in fitting function.
  // maxwellian_moms_ion[12]: ion prim vars.
  // maxwellian_moms_neut[20]: neut prim vars.
  // u_ion[12]: ion drift velocity vector (upar_i b_1, upar_i b_2, upar_i b_3).
  // v_sigma_cx: cell ave cross section fitting eqn.

  double m0_neut_av = 0.5*maxwellian_moms_neut[0]; 
 
  const double *vt_sq_ion = &maxwellian_moms_ion[8]; 
  const double *u_neut = &maxwellian_moms_neut[4]; 
  const double *vt_sq_neut = &maxwellian_moms_neut[16]; 
 
  double vt_sq_ion_av = 0.5*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.5*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.25*pow(u_neut[8],2)-0.5*u_ion[8]*u_neut[8]+0.25*pow(u_ion[8],2)+0.25*pow(u_neut[4],2)-0.5*u_ion[4]*u_neut[4]+0.25*pow(u_ion[4],2)+0.25*pow(u_neut[0],2)-0.5*u_ion[0]*u_neut[0]+0.25*pow(u_ion[0],2); 
 
  double v_cx = sqrt(2.5464790894703255*vt_sq_neut_av+2.5464790894703255*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 2.0*a*v_cx-2.0*b*v_cx*log(v_cx); 
 
  return 0.16666666666666666*v_sigma_cx[0]*m0_neut_av; 
  }
} 
 
GKYL_CU_DH static inline 
double sigma_cx_2x_ser_p2(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a: constant in fitting function.
  // b: constant in fitting function.
  // maxwellian_moms_ion[24]: ion prim vars.
  // maxwellian_moms_neut[40]: neut prim vars.
  // u_ion[24]: ion drift velocity vector (upar_i b_1, upar_i b_2, upar_i b_3).
  // v_sigma_cx: cell ave cross section fitting eqn.

  double m0_neut_av = 0.5*maxwellian_moms_neut[0]; 
 
  const double *vt_sq_ion = &maxwellian_moms_ion[16]; 
  const double *u_neut = &maxwellian_moms_neut[8]; 
  const double *vt_sq_neut = &maxwellian_moms_neut[32]; 
 
  double vt_sq_ion_av = 0.5*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.5*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.25*pow(u_neut[16],2)-0.5*u_ion[16]*u_neut[16]+0.25*pow(u_ion[16],2)+0.25*pow(u_neut[8],2)-0.5*u_ion[8]*u_neut[8]+0.25*pow(u_ion[8],2)+0.25*pow(u_neut[0],2)-0.5*u_ion[0]*u_neut[0]+0.25*pow(u_ion[0],2); 
 
  double v_cx = sqrt(2.5464790894703255*vt_sq_neut_av+2.5464790894703255*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 2.0*a*v_cx-2.0*b*v_cx*log(v_cx); 
 
  return 0.1*v_sigma_cx[0]*m0_neut_av; 
  }
} 
 
GKYL_CU_DH static inline 
double sigma_cx_3x_ser_p1(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *maxwellian_moms_ion, const double *maxwellian_moms_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) 
{ 
  // a: constant in fitting function.
  // b: constant in fitting function.
  // maxwellian_moms_ion[24]: ion prim vars.
  // maxwellian_moms_neut[40]: neut prim vars.
  // u_ion[24]: ion drift velocity vector (upar_i b_1, upar_i b_2, upar_i b_3).
  // v_sigma_cx: cell ave cross section fitting eqn.

  double m0_neut_av = 0.35355339059327384*maxwellian_moms_neut[0]; 
 
  const double *vt_sq_ion = &maxwellian_moms_ion[16]; 
  const double *u_neut = &maxwellian_moms_neut[8]; 
  const double *vt_sq_neut = &maxwellian_moms_neut[32]; 
 
  double vt_sq_ion_av = 0.35355339059327384*vt_sq_ion[0]; 
  double vt_sq_neut_av = 0.35355339059327384*vt_sq_neut[0]; 
  if ((vt_sq_ion_av > 0.) && (vt_sq_ion_av < vt_sq_ion_min)) vt_sq_ion_av = vt_sq_ion_min;
  if ((vt_sq_neut_av > 0.) && (vt_sq_neut_av < vt_sq_neut_min)) vt_sq_neut_av = vt_sq_neut_min;
  
  if (m0_neut_av <= 0 || vt_sq_neut_av <= 0 || vt_sq_ion_av <= 0) { 
    v_sigma_cx[0] = 0.0;
    return 0.0; 
  } else {
  double v_in_sq_av = 0.125*pow(u_neut[16],2)-0.25*u_ion[16]*u_neut[16]+0.125*pow(u_ion[16],2)+0.125*pow(u_neut[8],2)-0.25*u_ion[8]*u_neut[8]+0.125*pow(u_ion[8],2)+0.125*pow(u_neut[0],2)-0.25*u_ion[0]*u_neut[0]+0.125*pow(u_ion[0],2); 
 
  double v_cx = sqrt(2.5464790894703255*vt_sq_neut_av+2.5464790894703255*vt_sq_ion_av+v_in_sq_av);
  v_sigma_cx[0] = 2.8284271247461907*a*v_cx-2.8284271247461907*b*v_cx*log(v_cx); 
 
  return 0.11785113019775789*v_sigma_cx[0]*m0_neut_av; 
  }
} 
