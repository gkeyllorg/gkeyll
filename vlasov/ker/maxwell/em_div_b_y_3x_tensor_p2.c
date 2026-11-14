#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void em_div_b_y_3x_tensor_p2(const double *dxv, 
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  const double *bvar_c, double* GKYL_RESTRICT max_b, double* GKYL_RESTRICT div_b) 
{ 
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector in left/center/right cells in each direction. 
  // bvar_c:          Input volume expansion of magnetic field unit vector in center cell. 
  // max_b:           Output surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|. 
  // div_b:           Output volume expansion of div(b).

  const double dx1 = 2.0/dxv[1]; 
  const double *b_c = &bvar_c[27]; 
  const double *b_surf_lr = &bvar_surf_l[81]; 
  const double *b_surf_cl = &bvar_surf_c[72]; 
  const double *b_surf_cr = &bvar_surf_c[81]; 
  const double *b_surf_rl = &bvar_surf_r[72]; 

  double *max_b_l = &max_b[18]; 
  double *max_b_r = &max_b[27]; 

  double bl_r = 0.0; 
  double bc_l = 0.0; 
  double bc_r = 0.0; 
  double br_l = 0.0; 
  double max_b_quad_l[9] = {0.0}; 
  double max_b_quad_r[9] = {0.0}; 

  bl_r = 0.4*b_surf_lr[8]-0.5999999999999995*b_surf_lr[7]-0.5999999999999999*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]+0.9*b_surf_lr[3]-0.6708203932499369*b_surf_lr[2]-0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = 0.4*b_surf_cl[8]-0.5999999999999995*b_surf_cl[7]-0.5999999999999999*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]+0.9*b_surf_cl[3]-0.6708203932499369*b_surf_cl[2]-0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = 0.4*b_surf_cr[8]-0.5999999999999995*b_surf_cr[7]-0.5999999999999999*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]+0.9*b_surf_cr[3]-0.6708203932499369*b_surf_cr[2]-0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = 0.4*b_surf_rl[8]-0.5999999999999995*b_surf_rl[7]-0.5999999999999999*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]+0.9*b_surf_rl[3]-0.6708203932499369*b_surf_rl[2]-0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[0] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[0] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = (-0.5*b_surf_lr[8])+0.75*b_surf_lr[7]-0.5590169943749475*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]-0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = (-0.5*b_surf_cl[8])+0.75*b_surf_cl[7]-0.5590169943749475*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]-0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = (-0.5*b_surf_cr[8])+0.75*b_surf_cr[7]-0.5590169943749475*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]-0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = (-0.5*b_surf_rl[8])+0.75*b_surf_rl[7]-0.5590169943749475*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]-0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[1] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[1] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.4*b_surf_lr[8]-0.5999999999999995*b_surf_lr[7]+0.5999999999999999*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]-0.9*b_surf_lr[3]+0.6708203932499369*b_surf_lr[2]-0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = 0.4*b_surf_cl[8]-0.5999999999999995*b_surf_cl[7]+0.5999999999999999*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]-0.9*b_surf_cl[3]+0.6708203932499369*b_surf_cl[2]-0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = 0.4*b_surf_cr[8]-0.5999999999999995*b_surf_cr[7]+0.5999999999999999*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]-0.9*b_surf_cr[3]+0.6708203932499369*b_surf_cr[2]-0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = 0.4*b_surf_rl[8]-0.5999999999999995*b_surf_rl[7]+0.5999999999999999*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]-0.9*b_surf_rl[3]+0.6708203932499369*b_surf_rl[2]-0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[2] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[2] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = (-0.5*b_surf_lr[8])+0.75*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]-0.5590169943749475*b_surf_lr[4]-0.6708203932499369*b_surf_lr[2]+0.5*b_surf_lr[0]; 
  bc_l = (-0.5*b_surf_cl[8])+0.75*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]-0.5590169943749475*b_surf_cl[4]-0.6708203932499369*b_surf_cl[2]+0.5*b_surf_cl[0]; 
  bc_r = (-0.5*b_surf_cr[8])+0.75*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]-0.5590169943749475*b_surf_cr[4]-0.6708203932499369*b_surf_cr[2]+0.5*b_surf_cr[0]; 
  br_l = (-0.5*b_surf_rl[8])+0.75*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]-0.5590169943749475*b_surf_rl[4]-0.6708203932499369*b_surf_rl[2]+0.5*b_surf_rl[0]; 
  max_b_quad_l[3] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[3] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.625*b_surf_lr[8]-0.5590169943749475*b_surf_lr[5]-0.5590169943749475*b_surf_lr[4]+0.5*b_surf_lr[0]; 
  bc_l = 0.625*b_surf_cl[8]-0.5590169943749475*b_surf_cl[5]-0.5590169943749475*b_surf_cl[4]+0.5*b_surf_cl[0]; 
  bc_r = 0.625*b_surf_cr[8]-0.5590169943749475*b_surf_cr[5]-0.5590169943749475*b_surf_cr[4]+0.5*b_surf_cr[0]; 
  br_l = 0.625*b_surf_rl[8]-0.5590169943749475*b_surf_rl[5]-0.5590169943749475*b_surf_rl[4]+0.5*b_surf_rl[0]; 
  max_b_quad_l[4] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[4] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = (-0.5*b_surf_lr[8])-0.75*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]-0.5590169943749475*b_surf_lr[4]+0.6708203932499369*b_surf_lr[2]+0.5*b_surf_lr[0]; 
  bc_l = (-0.5*b_surf_cl[8])-0.75*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]-0.5590169943749475*b_surf_cl[4]+0.6708203932499369*b_surf_cl[2]+0.5*b_surf_cl[0]; 
  bc_r = (-0.5*b_surf_cr[8])-0.75*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]-0.5590169943749475*b_surf_cr[4]+0.6708203932499369*b_surf_cr[2]+0.5*b_surf_cr[0]; 
  br_l = (-0.5*b_surf_rl[8])-0.75*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]-0.5590169943749475*b_surf_rl[4]+0.6708203932499369*b_surf_rl[2]+0.5*b_surf_rl[0]; 
  max_b_quad_l[5] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[5] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.4*b_surf_lr[8]+0.5999999999999995*b_surf_lr[7]-0.5999999999999999*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]-0.9*b_surf_lr[3]-0.6708203932499369*b_surf_lr[2]+0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = 0.4*b_surf_cl[8]+0.5999999999999995*b_surf_cl[7]-0.5999999999999999*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]-0.9*b_surf_cl[3]-0.6708203932499369*b_surf_cl[2]+0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = 0.4*b_surf_cr[8]+0.5999999999999995*b_surf_cr[7]-0.5999999999999999*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]-0.9*b_surf_cr[3]-0.6708203932499369*b_surf_cr[2]+0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = 0.4*b_surf_rl[8]+0.5999999999999995*b_surf_rl[7]-0.5999999999999999*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]-0.9*b_surf_rl[3]-0.6708203932499369*b_surf_rl[2]+0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[6] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[6] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = (-0.5*b_surf_lr[8])-0.75*b_surf_lr[7]-0.5590169943749475*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]+0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = (-0.5*b_surf_cl[8])-0.75*b_surf_cl[7]-0.5590169943749475*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]+0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = (-0.5*b_surf_cr[8])-0.75*b_surf_cr[7]-0.5590169943749475*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]+0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = (-0.5*b_surf_rl[8])-0.75*b_surf_rl[7]-0.5590169943749475*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]+0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[7] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[7] = fmax(fabs(bc_r), fabs(br_l)); 

  bl_r = 0.4*b_surf_lr[8]+0.5999999999999995*b_surf_lr[7]+0.5999999999999999*b_surf_lr[6]+0.4472135954999579*b_surf_lr[5]+0.4472135954999579*b_surf_lr[4]+0.9*b_surf_lr[3]+0.6708203932499369*b_surf_lr[2]+0.6708203932499369*b_surf_lr[1]+0.5*b_surf_lr[0]; 
  bc_l = 0.4*b_surf_cl[8]+0.5999999999999995*b_surf_cl[7]+0.5999999999999999*b_surf_cl[6]+0.4472135954999579*b_surf_cl[5]+0.4472135954999579*b_surf_cl[4]+0.9*b_surf_cl[3]+0.6708203932499369*b_surf_cl[2]+0.6708203932499369*b_surf_cl[1]+0.5*b_surf_cl[0]; 
  bc_r = 0.4*b_surf_cr[8]+0.5999999999999995*b_surf_cr[7]+0.5999999999999999*b_surf_cr[6]+0.4472135954999579*b_surf_cr[5]+0.4472135954999579*b_surf_cr[4]+0.9*b_surf_cr[3]+0.6708203932499369*b_surf_cr[2]+0.6708203932499369*b_surf_cr[1]+0.5*b_surf_cr[0]; 
  br_l = 0.4*b_surf_rl[8]+0.5999999999999995*b_surf_rl[7]+0.5999999999999999*b_surf_rl[6]+0.4472135954999579*b_surf_rl[5]+0.4472135954999579*b_surf_rl[4]+0.9*b_surf_rl[3]+0.6708203932499369*b_surf_rl[2]+0.6708203932499369*b_surf_rl[1]+0.5*b_surf_rl[0]; 
  max_b_quad_l[8] = fmax(fabs(bl_r), fabs(bc_l)); 
  max_b_quad_r[8] = fmax(fabs(bc_r), fabs(br_l)); 

  tensor_3x_p2_upwind_quad_to_modal(max_b_quad_l, max_b_l); 
  tensor_3x_p2_upwind_quad_to_modal(max_b_quad_r, max_b_r); 

  div_b[0] += (0.3535533905932737*b_surf_rl[0]-0.3535533905932737*b_surf_lr[0]+0.3535533905932737*b_surf_cr[0]-0.3535533905932737*b_surf_cl[0])*dx1; 
  div_b[1] += (0.3535533905932737*b_surf_rl[1]-0.3535533905932737*b_surf_lr[1]+0.3535533905932737*b_surf_cr[1]-0.3535533905932737*b_surf_cl[1])*dx1; 
  div_b[2] += (0.6123724356957944*(b_surf_rl[0]+b_surf_lr[0]+b_surf_cr[0]+b_surf_cl[0])-1.732050807568877*b_c[0])*dx1; 
  div_b[3] += (0.3535533905932737*b_surf_rl[2]-0.3535533905932737*b_surf_lr[2]+0.3535533905932737*b_surf_cr[2]-0.3535533905932737*b_surf_cl[2])*dx1; 
  div_b[4] += (0.6123724356957944*(b_surf_rl[1]+b_surf_lr[1]+b_surf_cr[1]+b_surf_cl[1])-1.732050807568877*b_c[1])*dx1; 
  div_b[5] += (0.3535533905932737*b_surf_rl[3]-0.3535533905932737*b_surf_lr[3]+0.3535533905932737*b_surf_cr[3]-0.3535533905932737*b_surf_cl[3])*dx1; 
  div_b[6] += (0.6123724356957944*(b_surf_rl[2]+b_surf_lr[2]+b_surf_cr[2]+b_surf_cl[2])-1.732050807568877*b_c[3])*dx1; 
  div_b[7] += (0.3535533905932737*b_surf_rl[4]-0.3535533905932737*b_surf_lr[4]+0.3535533905932737*b_surf_cr[4]-0.3535533905932737*b_surf_cl[4])*dx1; 
  div_b[8] += ((-3.872983346207417*b_c[2])+0.7905694150420947*b_surf_rl[0]-0.7905694150420947*b_surf_lr[0]+0.7905694150420947*b_surf_cr[0]-0.7905694150420947*b_surf_cl[0])*dx1; 
  div_b[9] += (0.3535533905932737*b_surf_rl[5]-0.3535533905932737*b_surf_lr[5]+0.3535533905932737*b_surf_cr[5]-0.3535533905932737*b_surf_cl[5])*dx1; 
  div_b[10] += (0.6123724356957944*(b_surf_rl[3]+b_surf_lr[3]+b_surf_cr[3]+b_surf_cl[3])-1.732050807568877*b_c[5])*dx1; 
  div_b[11] += (0.6123724356957944*(b_surf_rl[4]+b_surf_lr[4]+b_surf_cr[4]+b_surf_cl[4])-1.732050807568877*b_c[7])*dx1; 
  div_b[12] += ((-3.872983346207417*b_c[4])+0.7905694150420948*b_surf_rl[1]-0.7905694150420948*b_surf_lr[1]+0.7905694150420948*b_surf_cr[1]-0.7905694150420948*b_surf_cl[1])*dx1; 
  div_b[13] += (0.3535533905932737*b_surf_rl[6]-0.3535533905932737*b_surf_lr[6]+0.3535533905932737*b_surf_cr[6]-0.3535533905932737*b_surf_cl[6])*dx1; 
  div_b[14] += ((-3.872983346207417*b_c[6])+0.7905694150420948*b_surf_rl[2]-0.7905694150420948*b_surf_lr[2]+0.7905694150420948*b_surf_cr[2]-0.7905694150420948*b_surf_cl[2])*dx1; 
  div_b[15] += (0.3535533905932737*b_surf_rl[7]-0.3535533905932737*b_surf_lr[7]+0.3535533905932737*b_surf_cr[7]-0.3535533905932737*b_surf_cl[7])*dx1; 
  div_b[16] += (0.6123724356957944*(b_surf_rl[5]+b_surf_lr[5]+b_surf_cr[5]+b_surf_cl[5])-1.732050807568877*b_c[9])*dx1; 
  div_b[17] += (0.6123724356957944*(b_surf_rl[6]+b_surf_lr[6]+b_surf_cr[6]+b_surf_cl[6])-1.732050807568877*b_c[13])*dx1; 
  div_b[18] += ((-3.872983346207417*b_c[10])+0.7905694150420947*b_surf_rl[3]-0.7905694150420947*b_surf_lr[3]+0.7905694150420947*b_surf_cr[3]-0.7905694150420947*b_surf_cl[3])*dx1; 
  div_b[19] += (0.6123724356957944*(b_surf_rl[7]+b_surf_lr[7]+b_surf_cr[7]+b_surf_cl[7])-1.732050807568877*b_c[15])*dx1; 
  div_b[20] += ((-3.872983346207417*b_c[11])+0.7905694150420947*b_surf_rl[4]-0.7905694150420947*b_surf_lr[4]+0.7905694150420947*b_surf_cr[4]-0.7905694150420947*b_surf_cl[4])*dx1; 
  div_b[21] += (0.3535533905932737*b_surf_rl[8]-0.3535533905932737*b_surf_lr[8]+0.3535533905932737*b_surf_cr[8]-0.3535533905932737*b_surf_cl[8])*dx1; 
  div_b[22] += ((-3.872983346207417*b_c[16])+0.7905694150420947*b_surf_rl[5]-0.7905694150420947*b_surf_lr[5]+0.7905694150420947*b_surf_cr[5]-0.7905694150420947*b_surf_cl[5])*dx1; 
  div_b[23] += ((-3.872983346207417*b_c[17])+0.7905694150420948*b_surf_rl[6]-0.7905694150420948*b_surf_lr[6]+0.7905694150420948*b_surf_cr[6]-0.7905694150420948*b_surf_cl[6])*dx1; 
  div_b[24] += (0.6123724356957944*(b_surf_rl[8]+b_surf_lr[8]+b_surf_cr[8]+b_surf_cl[8])-1.732050807568877*b_c[21])*dx1; 
  div_b[25] += ((-3.872983346207417*b_c[19])+0.7905694150420948*b_surf_rl[7]-0.7905694150420948*b_surf_lr[7]+0.7905694150420948*b_surf_cr[7]-0.7905694150420948*b_surf_cl[7])*dx1; 
  div_b[26] += ((-3.872983346207417*b_c[24])+0.7905694150420947*b_surf_rl[8]-0.7905694150420947*b_surf_lr[8]+0.7905694150420947*b_surf_cr[8]-0.7905694150420947*b_surf_cl[8])*dx1; 

} 
