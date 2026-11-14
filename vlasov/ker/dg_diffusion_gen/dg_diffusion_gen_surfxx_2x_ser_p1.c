#include <gkyl_dg_diffusion_gen_kernels.h>

GKYL_CU_DH double
dg_diffusion_gen_surfxx_2x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* q[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jxx = 4/dx[0]/dx[0];

  const double* Dxx = &Dij[0];
  const double* qlc = q[0];
  const double* qcc = q[1];
  const double* quc = q[2];

  out[0] += Jxx*((-0.46875*Dxx[3]*quc[3])-0.270632938682637*Dxx[2]*quc[3]-0.46875*Dxx[3]*qlc[3]+0.270632938682637*Dxx[2]*qlc[3]-0.9375*Dxx[3]*qcc[3]+0.4871392896287466*quc[2]*Dxx[3]-0.4871392896287466*qlc[2]*Dxx[3]+0.28125*Dxx[2]*quc[2]+0.28125*Dxx[2]*qlc[2]-0.5625*Dxx[2]*qcc[2]-0.46875*Dxx[1]*quc[1]-0.270632938682637*Dxx[0]*quc[1]-0.46875*Dxx[1]*qlc[1]+0.270632938682637*Dxx[0]*qlc[1]-0.9375*Dxx[1]*qcc[1]+0.4871392896287466*quc[0]*Dxx[1]-0.4871392896287466*qlc[0]*Dxx[1]+0.28125*Dxx[0]*quc[0]+0.28125*Dxx[0]*qlc[0]-0.5625*Dxx[0]*qcc[0]);
  out[1] += Jxx*((-0.3788861141556919*Dxx[3]*quc[3])-0.21875*Dxx[2]*quc[3]+0.3788861141556919*Dxx[3]*qlc[3]-0.21875*Dxx[2]*qlc[3]-1.4375*Dxx[2]*qcc[3]+0.46875*quc[2]*Dxx[3]+0.46875*qlc[2]*Dxx[3]-0.9375*qcc[2]*Dxx[3]+0.270632938682637*Dxx[2]*quc[2]-0.270632938682637*Dxx[2]*qlc[2]-0.3788861141556919*Dxx[1]*quc[1]-0.21875*Dxx[0]*quc[1]+0.3788861141556919*Dxx[1]*qlc[1]-0.21875*Dxx[0]*qlc[1]-1.4375*Dxx[0]*qcc[1]+0.46875*quc[0]*Dxx[1]+0.46875*qlc[0]*Dxx[1]-0.9375*qcc[0]*Dxx[1]+0.270632938682637*Dxx[0]*quc[0]-0.270632938682637*Dxx[0]*qlc[0]);
  out[2] += Jxx*((-0.46875*Dxx[1]*quc[3])-0.270632938682637*Dxx[0]*quc[3]-0.46875*Dxx[1]*qlc[3]+0.270632938682637*Dxx[0]*qlc[3]-0.9375*Dxx[1]*qcc[3]-0.46875*quc[1]*Dxx[3]-0.46875*qlc[1]*Dxx[3]-0.9375*qcc[1]*Dxx[3]+0.4871392896287466*quc[0]*Dxx[3]-0.4871392896287466*qlc[0]*Dxx[3]+0.4871392896287466*Dxx[1]*quc[2]+0.28125*Dxx[0]*quc[2]-0.4871392896287466*Dxx[1]*qlc[2]+0.28125*Dxx[0]*qlc[2]-0.5625*Dxx[0]*qcc[2]-0.270632938682637*quc[1]*Dxx[2]+0.270632938682637*qlc[1]*Dxx[2]+0.28125*quc[0]*Dxx[2]+0.28125*qlc[0]*Dxx[2]-0.5625*qcc[0]*Dxx[2]);
  out[3] += Jxx*((-0.3788861141556919*Dxx[1]*quc[3])-0.21875*Dxx[0]*quc[3]+0.3788861141556919*Dxx[1]*qlc[3]-0.21875*Dxx[0]*qlc[3]-1.4375*Dxx[0]*qcc[3]-0.3788861141556919*quc[1]*Dxx[3]+0.3788861141556919*qlc[1]*Dxx[3]+0.46875*quc[0]*Dxx[3]+0.46875*qlc[0]*Dxx[3]-0.9375*qcc[0]*Dxx[3]+0.46875*Dxx[1]*quc[2]+0.270632938682637*Dxx[0]*quc[2]+0.46875*Dxx[1]*qlc[2]-0.270632938682637*Dxx[0]*qlc[2]-0.9375*Dxx[1]*qcc[2]-0.21875*quc[1]*Dxx[2]-0.21875*qlc[1]*Dxx[2]-1.4375*qcc[1]*Dxx[2]+0.270632938682637*quc[0]*Dxx[2]-0.270632938682637*qlc[0]*Dxx[2]);
  return 0.;
}
