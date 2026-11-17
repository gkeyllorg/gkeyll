#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
GKYL_CU_DH double iso_euler_surfx_1x_ser_p2(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn,
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r,
    const double *fluid_l, const double *fluid_c, const double *fluid_r,
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:      Cell-center coordinates.
  // dxv[NDIM]:    Cell spacing.
  // wv_eqn:       Wave equation: the Riemann solver for the flux at the interfaces.
  // geom_l:       Geometry for the left surface update.
  // geom_r:       Geometry for the right surface update.
  // fluid_l/c/r:  [rho, rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // out:          Incremented output.

  const double dx1 = 2.0/dxv[0]; 
  const double vt = container_of(wv_eqn, struct wv_iso_euler, eqn)->vt; 

  const double *rho_l = &fluid_l[0]; 
  const double *rhoux_l = &fluid_l[3]; 
  const double *rhouy_l = &fluid_l[6]; 
  const double *rhouz_l = &fluid_l[9]; 

  const double *rho_c = &fluid_c[0]; 
  const double *rhoux_c = &fluid_c[3]; 
  const double *rhouy_c = &fluid_c[6]; 
  const double *rhouz_c = &fluid_c[9]; 

  const double *rho_r = &fluid_r[0]; 
  const double *rhoux_r = &fluid_r[3]; 
  const double *rhouy_r = &fluid_r[6]; 
  const double *rhouz_r = &fluid_r[9]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[3]; 
  double *outrhouy = &out[6]; 
  double *outrhouz = &out[9]; 

  double q_lr[4] = {0.0}; 
  double q_lr_local[4] = {0.0}; 
  double q_cl[4] = {0.0}; 
  double q_cl_local[4] = {0.0}; 
  double q_cr[4] = {0.0}; 
  double q_cr_local[4] = {0.0}; 
  double q_rl[4] = {0.0}; 
  double q_rl_local[4] = {0.0}; 
  double delta_l[4] = {0.0}; 
  double waves_l[12] = {0.0}; 
  double speeds_l[3] = {0.0}; 
  double amdq_l_local[4] = {0.0}; 
  double apdq_l_local[4] = {0.0}; 
  double fhat_l_local[4] = {0.0}; 
  double fhat_l[4] = {0.0}; 
  double lenr_l = geom_l->lenr[0]; 
  double delta_r[4] = {0.0}; 
  double waves_r[12] = {0.0}; 
  double speeds_r[3] = {0.0}; 
  double amdq_r_local[4] = {0.0}; 
  double apdq_r_local[4] = {0.0}; 
  double fhat_r_local[4] = {0.0}; 
  double fhat_r[4] = {0.0}; 
  double lenr_r = geom_r->lenr[0]; 
  double flux_rho_l[1] = {0.0}; 
  double flux_rho_r[1] = {0.0}; 
  double flux_rhoux_l[1] = {0.0}; 
  double flux_rhoux_r[1] = {0.0}; 
  double flux_rhouy_l[1] = {0.0}; 
  double flux_rhouy_r[1] = {0.0}; 
  double flux_rhouz_l[1] = {0.0}; 
  double flux_rhouz_r[1] = {0.0}; 

  q_lr[0] = 1.5811388300841895*rho_l[2]+1.224744871391589*rho_l[1]+0.7071067811865475*rho_l[0]; 
  q_lr[1] = 1.5811388300841895*rhoux_l[2]+1.224744871391589*rhoux_l[1]+0.7071067811865475*rhoux_l[0]; 
  q_lr[2] = 1.5811388300841895*rhouy_l[2]+1.224744871391589*rhouy_l[1]+0.7071067811865475*rhouy_l[0]; 
  q_lr[3] = 1.5811388300841895*rhouz_l[2]+1.224744871391589*rhouz_l[1]+0.7071067811865475*rhouz_l[0]; 
  q_cl[0] = 1.5811388300841895*rho_c[2]-1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0]; 
  q_cl[1] = 1.5811388300841895*rhoux_c[2]-1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0]; 
  q_cl[2] = 1.5811388300841895*rhouy_c[2]-1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0]; 
  q_cl[3] = 1.5811388300841895*rhouz_c[2]-1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0]; 
  q_cr[0] = 1.5811388300841895*rho_c[2]+1.224744871391589*rho_c[1]+0.7071067811865475*rho_c[0]; 
  q_cr[1] = 1.5811388300841895*rhoux_c[2]+1.224744871391589*rhoux_c[1]+0.7071067811865475*rhoux_c[0]; 
  q_cr[2] = 1.5811388300841895*rhouy_c[2]+1.224744871391589*rhouy_c[1]+0.7071067811865475*rhouy_c[0]; 
  q_cr[3] = 1.5811388300841895*rhouz_c[2]+1.224744871391589*rhouz_c[1]+0.7071067811865475*rhouz_c[0]; 
  q_rl[0] = 1.5811388300841895*rho_r[2]-1.224744871391589*rho_r[1]+0.7071067811865475*rho_r[0]; 
  q_rl[1] = 1.5811388300841895*rhoux_r[2]-1.224744871391589*rhoux_r[1]+0.7071067811865475*rhoux_r[0]; 
  q_rl[2] = 1.5811388300841895*rhouy_r[2]-1.224744871391589*rhouy_r[1]+0.7071067811865475*rhouy_r[0]; 
  q_rl[3] = 1.5811388300841895*rhouz_r[2]-1.224744871391589*rhouz_r[1]+0.7071067811865475*rhouz_r[0]; 

  iso_euler_rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  iso_euler_rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  iso_euler_rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  iso_euler_rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  for (int i = 0; i < 4; ++i) { 
    delta_l[i] = q_cl_local[i] - q_lr_local[i]; 
    delta_r[i] = q_rl_local[i] - q_cr_local[i]; 
  } 
  gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, 1.0, 1.0, waves_l, speeds_l); 
  gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, 1.0, 1.0, waves_r, speeds_r); 
  for (int i = 0; i < 3; ++i) { 
    speeds_l[i] *= lenr_l; 
    speeds_r[i] *= lenr_r; 
  } 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, 1.0, 1.0, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, 1.0, 1.0, waves_r, speeds_r, amdq_r_local, apdq_r_local); 
  iso_euler_flux(vt, q_lr_local, fhat_l_local); 
  iso_euler_flux(vt, q_cr_local, fhat_r_local); 
  for (int i = 0; i < 4; ++i) { 
    fhat_l_local[i] = lenr_l*fhat_l_local[i] + amdq_l_local[i]; 
    fhat_r_local[i] = lenr_r*fhat_r_local[i] + amdq_r_local[i]; 
  } 
  iso_euler_rot_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], fhat_l_local, fhat_l); 
  iso_euler_rot_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], fhat_r_local, fhat_r); 
  flux_rho_l[0] = fhat_l[0]; 
  flux_rho_r[0] = fhat_r[0]; 
  flux_rhoux_l[0] = fhat_l[1]; 
  flux_rhoux_r[0] = fhat_r[1]; 
  flux_rhouy_l[0] = fhat_l[2]; 
  flux_rhouy_r[0] = fhat_r[2]; 
  flux_rhouz_l[0] = fhat_l[3]; 
  flux_rhouz_r[0] = fhat_r[3]; 

  outrho[0] += (0.7071067811865475*flux_rho_l[0]-0.7071067811865475*flux_rho_r[0])*dx1; 
  outrho[1] += -(1.224744871391589*(flux_rho_r[0]+flux_rho_l[0])*dx1); 
  outrho[2] += (1.5811388300841895*flux_rho_l[0]-1.5811388300841895*flux_rho_r[0])*dx1; 

  outrhoux[0] += (0.7071067811865475*flux_rhoux_l[0]-0.7071067811865475*flux_rhoux_r[0])*dx1; 
  outrhoux[1] += -(1.224744871391589*(flux_rhoux_r[0]+flux_rhoux_l[0])*dx1); 
  outrhoux[2] += (1.5811388300841895*flux_rhoux_l[0]-1.5811388300841895*flux_rhoux_r[0])*dx1; 

  outrhouy[0] += (0.7071067811865475*flux_rhouy_l[0]-0.7071067811865475*flux_rhouy_r[0])*dx1; 
  outrhouy[1] += -(1.224744871391589*(flux_rhouy_r[0]+flux_rhouy_l[0])*dx1); 
  outrhouy[2] += (1.5811388300841895*flux_rhouy_l[0]-1.5811388300841895*flux_rhouy_r[0])*dx1; 

  outrhouz[0] += (0.7071067811865475*flux_rhouz_l[0]-0.7071067811865475*flux_rhouz_r[0])*dx1; 
  outrhouz[1] += -(1.224744871391589*(flux_rhouz_r[0]+flux_rhouz_l[0])*dx1); 
  outrhouz[2] += (1.5811388300841895*flux_rhouz_l[0]-1.5811388300841895*flux_rhouz_r[0])*dx1; 

  return 0.;

} 
