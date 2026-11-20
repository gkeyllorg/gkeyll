#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_bvar_grad_x_2x_ser_p1(const double *dxv, const double *jacob_pos, const double *b_l, const double *b_c, const double *b_r, double* GKYL_RESTRICT grad_b) 
{ 
  // dxv[NDIM]: Cell spacing (computational). 
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 2 equal entries). 
  // b_l/b_c/b_r: Magnetic field unit vector [bx, by, bz] in the left/current/right cell in direction x. 
  // grad_b:    Output weak gradient of b, block of direction x: [d_x bx, d_x by, d_x bz]. 
  const double dx1 = 2.0/(dxv[0]*jacob_pos[0]); 
  const double *bx_l = &b_l[0]; 
  const double *bx_c = &b_c[0]; 
  const double *bx_r = &b_r[0]; 
  double *grad_bx = &grad_b[0]; 
  const double *by_l = &b_l[4]; 
  const double *by_c = &b_c[4]; 
  const double *by_r = &b_r[4]; 
  double *grad_by = &grad_b[4]; 
  const double *bz_l = &b_l[8]; 
  const double *bz_c = &b_c[8]; 
  const double *bz_r = &b_r[8]; 
  double *grad_bz = &grad_b[8]; 
 
  grad_bx[0] = -(0.4330127018922193*bx_r[1]*dx1)-0.4330127018922193*bx_l[1]*dx1+0.8660254037844386*bx_c[1]*dx1+0.25*bx_r[0]*dx1-0.25*bx_l[0]*dx1; 
  grad_bx[1] = -(0.75*bx_r[1]*dx1)+0.75*bx_l[1]*dx1+0.4330127018922193*bx_r[0]*dx1+0.4330127018922193*bx_l[0]*dx1-0.8660254037844385*bx_c[0]*dx1; 
  grad_bx[2] = -(0.4330127018922193*bx_r[3]*dx1)-0.4330127018922193*bx_l[3]*dx1+0.8660254037844386*bx_c[3]*dx1+0.25*bx_r[2]*dx1-0.25*bx_l[2]*dx1; 
  grad_bx[3] = -(0.75*bx_r[3]*dx1)+0.75*bx_l[3]*dx1+0.4330127018922193*bx_r[2]*dx1+0.4330127018922193*bx_l[2]*dx1-0.8660254037844385*bx_c[2]*dx1; 
  grad_by[0] = -(0.4330127018922193*by_r[1]*dx1)-0.4330127018922193*by_l[1]*dx1+0.8660254037844386*by_c[1]*dx1+0.25*by_r[0]*dx1-0.25*by_l[0]*dx1; 
  grad_by[1] = -(0.75*by_r[1]*dx1)+0.75*by_l[1]*dx1+0.4330127018922193*by_r[0]*dx1+0.4330127018922193*by_l[0]*dx1-0.8660254037844385*by_c[0]*dx1; 
  grad_by[2] = -(0.4330127018922193*by_r[3]*dx1)-0.4330127018922193*by_l[3]*dx1+0.8660254037844386*by_c[3]*dx1+0.25*by_r[2]*dx1-0.25*by_l[2]*dx1; 
  grad_by[3] = -(0.75*by_r[3]*dx1)+0.75*by_l[3]*dx1+0.4330127018922193*by_r[2]*dx1+0.4330127018922193*by_l[2]*dx1-0.8660254037844385*by_c[2]*dx1; 
  grad_bz[0] = -(0.4330127018922193*bz_r[1]*dx1)-0.4330127018922193*bz_l[1]*dx1+0.8660254037844386*bz_c[1]*dx1+0.25*bz_r[0]*dx1-0.25*bz_l[0]*dx1; 
  grad_bz[1] = -(0.75*bz_r[1]*dx1)+0.75*bz_l[1]*dx1+0.4330127018922193*bz_r[0]*dx1+0.4330127018922193*bz_l[0]*dx1-0.8660254037844385*bz_c[0]*dx1; 
  grad_bz[2] = -(0.4330127018922193*bz_r[3]*dx1)-0.4330127018922193*bz_l[3]*dx1+0.8660254037844386*bz_c[3]*dx1+0.25*bz_r[2]*dx1-0.25*bz_l[2]*dx1; 
  grad_bz[3] = -(0.75*bz_r[3]*dx1)+0.75*bz_l[3]*dx1+0.4330127018922193*bz_r[2]*dx1+0.4330127018922193*bz_l[2]*dx1-0.8660254037844385*bz_c[2]*dx1; 
} 
