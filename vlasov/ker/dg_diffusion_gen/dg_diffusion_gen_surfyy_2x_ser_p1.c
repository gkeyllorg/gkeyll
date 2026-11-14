#include <gkyl_dg_diffusion_gen_kernels.h>

GKYL_CU_DH double
dg_diffusion_gen_surfyy_2x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* q[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jyy = 4/dx[1]/dx[1];

  const double* Dyy = &Dij[8];

  const double* qcl = q[0];
  const double* qcc = q[1];
  const double* qcu = q[2];

  out[0] += Jyy*((-0.46875*Dyy[3]*qcu[3])-0.270632938682637*Dyy[1]*qcu[3]-0.46875*Dyy[3]*qcl[3]+0.270632938682637*Dyy[1]*qcl[3]-0.9375*Dyy[3]*qcc[3]+0.4871392896287466*qcu[1]*Dyy[3]-0.4871392896287466*qcl[1]*Dyy[3]-0.46875*Dyy[2]*qcu[2]-0.270632938682637*Dyy[0]*qcu[2]-0.46875*Dyy[2]*qcl[2]+0.270632938682637*Dyy[0]*qcl[2]-0.9375*Dyy[2]*qcc[2]+0.4871392896287466*qcu[0]*Dyy[2]-0.4871392896287466*qcl[0]*Dyy[2]+0.28125*Dyy[1]*qcu[1]+0.28125*Dyy[1]*qcl[1]-0.5625*Dyy[1]*qcc[1]+0.28125*Dyy[0]*qcu[0]+0.28125*Dyy[0]*qcl[0]-0.5625*Dyy[0]*qcc[0]);
  out[1] += Jyy*((-0.46875*Dyy[2]*qcu[3])-0.270632938682637*Dyy[0]*qcu[3]-0.46875*Dyy[2]*qcl[3]+0.270632938682637*Dyy[0]*qcl[3]-0.9375*Dyy[2]*qcc[3]-0.46875*qcu[2]*Dyy[3]-0.46875*qcl[2]*Dyy[3]-0.9375*qcc[2]*Dyy[3]+0.4871392896287466*qcu[0]*Dyy[3]-0.4871392896287466*qcl[0]*Dyy[3]-0.270632938682637*Dyy[1]*qcu[2]+0.270632938682637*Dyy[1]*qcl[2]+0.4871392896287466*qcu[1]*Dyy[2]-0.4871392896287466*qcl[1]*Dyy[2]+0.28125*Dyy[0]*qcu[1]+0.28125*Dyy[0]*qcl[1]-0.5625*Dyy[0]*qcc[1]+0.28125*qcu[0]*Dyy[1]+0.28125*qcl[0]*Dyy[1]-0.5625*qcc[0]*Dyy[1]);
  out[2] += Jyy*((-0.3788861141556919*Dyy[3]*qcu[3])-0.21875*Dyy[1]*qcu[3]+0.3788861141556919*Dyy[3]*qcl[3]-0.21875*Dyy[1]*qcl[3]-1.4375*Dyy[1]*qcc[3]+0.46875*qcu[1]*Dyy[3]+0.46875*qcl[1]*Dyy[3]-0.9375*qcc[1]*Dyy[3]-0.3788861141556919*Dyy[2]*qcu[2]-0.21875*Dyy[0]*qcu[2]+0.3788861141556919*Dyy[2]*qcl[2]-0.21875*Dyy[0]*qcl[2]-1.4375*Dyy[0]*qcc[2]+0.46875*qcu[0]*Dyy[2]+0.46875*qcl[0]*Dyy[2]-0.9375*qcc[0]*Dyy[2]+0.270632938682637*Dyy[1]*qcu[1]-0.270632938682637*Dyy[1]*qcl[1]+0.270632938682637*Dyy[0]*qcu[0]-0.270632938682637*Dyy[0]*qcl[0]);
  out[3] += Jyy*((-0.3788861141556919*Dyy[2]*qcu[3])-0.21875*Dyy[0]*qcu[3]+0.3788861141556919*Dyy[2]*qcl[3]-0.21875*Dyy[0]*qcl[3]-1.4375*Dyy[0]*qcc[3]-0.3788861141556919*qcu[2]*Dyy[3]+0.3788861141556919*qcl[2]*Dyy[3]+0.46875*qcu[0]*Dyy[3]+0.46875*qcl[0]*Dyy[3]-0.9375*qcc[0]*Dyy[3]-0.21875*Dyy[1]*qcu[2]-0.21875*Dyy[1]*qcl[2]-1.4375*Dyy[1]*qcc[2]+0.46875*qcu[1]*Dyy[2]+0.46875*qcl[1]*Dyy[2]-0.9375*qcc[1]*Dyy[2]+0.270632938682637*Dyy[0]*qcu[1]-0.270632938682637*Dyy[0]*qcl[1]+0.270632938682637*qcu[0]*Dyy[1]-0.270632938682637*qcl[0]*Dyy[1]);
  return 0.;
}
