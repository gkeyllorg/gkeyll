#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_in_cell_grad_x_1x_ser_p1(const double *dxv, const double *jacob_pos, double blend_lo, double blend_up,
  const double *q_l, const double *q_c, const double *q_r, double* GKYL_RESTRICT grad) 
{ 
  // dxv[NDIM]: Cell spacing (computational). 
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 2 equal entries). 
  // blend_lo/up: Weight of the label-paired flux in the lower/upper face of direction x. 
  // q_l/q_c/q_r: Three-component field in the left/current/right cell in direction x. 
  // grad:      Output gradient, block of direction x: [d_x qx, d_x qy, d_x qz]. 
  const double dx1 = 2.0/(dxv[0]*jacob_pos[0]); 
  const double *qx_l = &q_l[0]; 
  const double *qx_c = &q_c[0]; 
  const double *qx_r = &q_r[0]; 
  double *grad_qx = &grad[0]; 
  const double *qy_l = &q_l[2]; 
  const double *qy_c = &q_c[2]; 
  const double *qy_r = &q_r[2]; 
  double *grad_qy = &grad[2]; 
  const double *qz_l = &q_l[4]; 
  const double *qz_c = &q_c[4]; 
  const double *qz_r = &q_r[4]; 
  double *grad_qz = &grad[4]; 
 
  grad_qx[0] = -(0.4330127018922193*qx_r[1]*blend_up*dx1)-0.4330127018922193*qx_c[1]*blend_up*dx1+0.25*qx_r[0]*blend_up*dx1-0.25*qx_c[0]*blend_up*dx1-0.4330127018922193*qx_l[1]*blend_lo*dx1-0.4330127018922193*qx_c[1]*blend_lo*dx1-0.25*qx_l[0]*blend_lo*dx1+0.25*qx_c[0]*blend_lo*dx1+1.7320508075688772*qx_c[1]*dx1; 
  grad_qx[1] = -(0.75*qx_r[1]*blend_up*dx1)-0.75*qx_c[1]*blend_up*dx1+0.4330127018922193*qx_r[0]*blend_up*dx1-0.4330127018922193*qx_c[0]*blend_up*dx1+0.75*qx_l[1]*blend_lo*dx1+0.75*qx_c[1]*blend_lo*dx1+0.4330127018922193*qx_l[0]*blend_lo*dx1-0.4330127018922193*qx_c[0]*blend_lo*dx1; 
  grad_qy[0] = -(0.4330127018922193*qy_r[1]*blend_up*dx1)-0.4330127018922193*qy_c[1]*blend_up*dx1+0.25*qy_r[0]*blend_up*dx1-0.25*qy_c[0]*blend_up*dx1-0.4330127018922193*qy_l[1]*blend_lo*dx1-0.4330127018922193*qy_c[1]*blend_lo*dx1-0.25*qy_l[0]*blend_lo*dx1+0.25*qy_c[0]*blend_lo*dx1+1.7320508075688772*qy_c[1]*dx1; 
  grad_qy[1] = -(0.75*qy_r[1]*blend_up*dx1)-0.75*qy_c[1]*blend_up*dx1+0.4330127018922193*qy_r[0]*blend_up*dx1-0.4330127018922193*qy_c[0]*blend_up*dx1+0.75*qy_l[1]*blend_lo*dx1+0.75*qy_c[1]*blend_lo*dx1+0.4330127018922193*qy_l[0]*blend_lo*dx1-0.4330127018922193*qy_c[0]*blend_lo*dx1; 
  grad_qz[0] = -(0.4330127018922193*qz_r[1]*blend_up*dx1)-0.4330127018922193*qz_c[1]*blend_up*dx1+0.25*qz_r[0]*blend_up*dx1-0.25*qz_c[0]*blend_up*dx1-0.4330127018922193*qz_l[1]*blend_lo*dx1-0.4330127018922193*qz_c[1]*blend_lo*dx1-0.25*qz_l[0]*blend_lo*dx1+0.25*qz_c[0]*blend_lo*dx1+1.7320508075688772*qz_c[1]*dx1; 
  grad_qz[1] = -(0.75*qz_r[1]*blend_up*dx1)-0.75*qz_c[1]*blend_up*dx1+0.4330127018922193*qz_r[0]*blend_up*dx1-0.4330127018922193*qz_c[0]*blend_up*dx1+0.75*qz_l[1]*blend_lo*dx1+0.75*qz_c[1]*blend_lo*dx1+0.4330127018922193*qz_l[0]*blend_lo*dx1-0.4330127018922193*qz_c[0]*blend_lo*dx1; 
} 
