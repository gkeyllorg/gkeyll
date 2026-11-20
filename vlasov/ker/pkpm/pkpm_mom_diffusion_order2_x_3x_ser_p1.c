#include <gkyl_pkpm_vars_kernels.h> 

// The derivatives of orders 0 and 1 (in units of dxi/2) at a face of the cubic whose p=1 projections 
// on the cells of half-widths hl (left) and hr (right) are (al0, al1) and (ar0, ar1). 
GKYL_CU_DH static inline void pkpm_mom_diffusion_order2_x_3x_ser_p1_face(double hl, double hr, double al0, double al1, double ar0, double ar1,
  double* GKYL_RESTRICT dq) 
{ 
  const double s = hl + hr, rd = 1.0/(hl*hr*s*s*s); 
  dq[0] = (1.224744871391589*al1*hl*pow(hr,4)+0.7071067811865475*al0*hl*pow(hr,4)+2.041241452319315*al1*pow(hl,2)*pow(hr,3)+2.1213203435596424*al0*pow(hl,2)*pow(hr,3)-2.041241452319315*ar1*pow(hl,3)*pow(hr,2)+2.1213203435596424*ar0*pow(hl,3)*pow(hr,2)-1.224744871391589*ar1*pow(hl,4)*hr+0.7071067811865475*ar0*pow(hl,4)*hr)*rd; 
  dq[1] = (1.224744871391589*al1*pow(hr,4)-1.224744871391589*al1*hl*pow(hr,3)-6.123724356957944*ar1*pow(hl,2)*pow(hr,2)+6.363961030678928*ar0*pow(hl,2)*pow(hr,2)-6.123724356957944*al1*pow(hl,2)*pow(hr,2)-6.363961030678928*al0*pow(hl,2)*pow(hr,2)-1.224744871391589*ar1*pow(hl,3)*hr+1.224744871391589*ar1*pow(hl,4))*rd; 
} 

// Recovery-DG diffusion of order 2 of one p=1 mode pair along the direction: increments of the 
// current cell from the data m = (J M)/J_d of the left, current, right cells of half-widths jl, jc, jr 
// (in units of dxi/2), without the factor (-1)^(n/2+1) D (2/dxi)^n. 
GKYL_CU_DH static inline void pkpm_mom_diffusion_order2_x_3x_ser_p1_1d(double jl, double jc, double jr, double ml0, double ml1, double mc0,
  double mc1, double mr0, double mr1, double* GKYL_RESTRICT o) 
{ 
  double dR[2], dL[2]; 
  pkpm_mom_diffusion_order2_x_3x_ser_p1_face(jc, jr, mc0, mc1, mr0, mr1, dR); 
  pkpm_mom_diffusion_order2_x_3x_ser_p1_face(jl, jc, ml0, ml1, mc0, mc1, dL); 
  o[0] = 0.7071067811865475*(dR[1] - dL[1]); 
  o[1] = 1.224744871391589*(dR[1] + dL[1]) - (1.224744871391589/jc)*(dR[0] - dL[0]); 
} 

GKYL_CU_DH double pkpm_mom_diffusion_order2_x_3x_ser_p1(const double *dxv, double D, const double *jac_l, const double *jac_c,
  const double *jac_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]:  Cell spacing (computational). 
  // D:          Diffusion coefficient of direction x. 
  // jac_l/c/r:  Position-map Jacobian per direction of the left/current/right cells (cdim blocks of 2 equal entries). 
  // fl/fc/fr:   J-weighted fluid state [rho, M_x, M_y, M_z, E] of the left/current/right cells. 
  // out:        Incremented right-hand side (momentum blocks). 
  // Returns the time-step rate. 
  const double jl = jac_l[0], jc = jac_c[0], jr = jac_r[0]; 
  const double rjl = 1.0/jl, rjc = 1.0/jc, rjr = 1.0/jr; 
  const double fac = 1.0*D*pow(2.0/dxv[0], 2); 
  double o[2]; 
  for (int c = 1; c < 4; ++c) { 
    const double *ql = &fl[c*8], *qc = &fc[c*8], *qr = &fr[c*8]; 
    double *oc = &out[c*8]; 
    pkpm_mom_diffusion_order2_x_3x_ser_p1_1d(jl, jc, jr, ql[0]*rjl, ql[1]*rjl, qc[0]*rjc, qc[1]*rjc, qr[0]*rjr, qr[1]*rjr, o); 
    oc[0] += fac*o[0]; oc[1] += fac*o[1]; 
    pkpm_mom_diffusion_order2_x_3x_ser_p1_1d(jl, jc, jr, ql[2]*rjl, ql[4]*rjl, qc[2]*rjc, qc[4]*rjc, qr[2]*rjr, qr[4]*rjr, o); 
    oc[2] += fac*o[0]; oc[4] += fac*o[1]; 
    pkpm_mom_diffusion_order2_x_3x_ser_p1_1d(jl, jc, jr, ql[3]*rjl, ql[5]*rjl, qc[3]*rjc, qc[5]*rjc, qr[3]*rjr, qr[5]*rjr, o); 
    oc[3] += fac*o[0]; oc[5] += fac*o[1]; 
    pkpm_mom_diffusion_order2_x_3x_ser_p1_1d(jl, jc, jr, ql[6]*rjl, ql[7]*rjl, qc[6]*rjc, qc[7]*rjc, qr[6]*rjr, qr[7]*rjr, o); 
    oc[6] += fac*o[0]; oc[7] += fac*o[1]; 
  } 
  return 4*D*pow(2.0/(dxv[0]*jc), 2); 
} 
