#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_surfy_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
    const struct gkyl_wv_eqn *wv_eqn,
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r,
    const double *fluid_l, const double *fluid_c, const double *fluid_r,
    double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:      Cell-center coordinates.
  // dxv[NDIM]:    Cell spacing (computational).
  // jacob_pos:    Position-map Jacobian per configuration direction of the cell (cdim blocks of
  //               2 equal entries): its physical size in direction d is J_d dxv[d].
  // wv_eqn:       Wave equation: the Riemann solver for the flux at the interfaces.
  // geom_l:       Geometry for the left surface update.
  // geom_r:       Geometry for the right surface update.
  // fluid_l/c/r:  [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.
  // out:          Incremented output.

  const double dx1 = 2.0/(dxv[1]*jacob_pos[2]); 
  const double gas_gamma = container_of(wv_eqn, struct wv_euler, eqn)->gas_gamma; 

  const double *rho_l = &fluid_l[0]; 
  const double *rhoux_l = &fluid_l[4]; 
  const double *rhouy_l = &fluid_l[8]; 
  const double *rhouz_l = &fluid_l[12]; 
  const double *energy_l = &fluid_l[16]; 

  const double *rho_c = &fluid_c[0]; 
  const double *rhoux_c = &fluid_c[4]; 
  const double *rhouy_c = &fluid_c[8]; 
  const double *rhouz_c = &fluid_c[12]; 
  const double *energy_c = &fluid_c[16]; 

  const double *rho_r = &fluid_r[0]; 
  const double *rhoux_r = &fluid_r[4]; 
  const double *rhouy_r = &fluid_r[8]; 
  const double *rhouz_r = &fluid_r[12]; 
  const double *energy_r = &fluid_r[16]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[4]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[12]; 
  double *outenergy = &out[16]; 

  double q_lr[5] = {0.0}; 
  double q_lr_local[5] = {0.0}; 
  double q_cl[5] = {0.0}; 
  double q_cl_local[5] = {0.0}; 
  double q_cr[5] = {0.0}; 
  double q_cr_local[5] = {0.0}; 
  double q_rl[5] = {0.0}; 
  double q_rl_local[5] = {0.0}; 
  double delta_l[5] = {0.0}; 
  double waves_l[15] = {0.0}; 
  double speeds_l[3] = {0.0}; 
  double amdq_l_local[5] = {0.0}; 
  double apdq_l_local[5] = {0.0}; 
  double fhat_l_local[5] = {0.0}; 
  double fhat_l[5] = {0.0}; 
  double lenr_l = geom_l->lenr[1]; 
  double delta_r[5] = {0.0}; 
  double waves_r[15] = {0.0}; 
  double speeds_r[3] = {0.0}; 
  double amdq_r_local[5] = {0.0}; 
  double apdq_r_local[5] = {0.0}; 
  double fhat_r_local[5] = {0.0}; 
  double fhat_r[5] = {0.0}; 
  double lenr_r = geom_r->lenr[1]; 
  double flux_rho_quad_l[2] = {0.0}; 
  double flux_rho_l[2] = {0.0}; 
  double flux_rho_quad_r[2] = {0.0}; 
  double flux_rho_r[2] = {0.0}; 
  double flux_rhoux_quad_l[2] = {0.0}; 
  double flux_rhoux_l[2] = {0.0}; 
  double flux_rhoux_quad_r[2] = {0.0}; 
  double flux_rhoux_r[2] = {0.0}; 
  double flux_rhouy_quad_l[2] = {0.0}; 
  double flux_rhouy_l[2] = {0.0}; 
  double flux_rhouy_quad_r[2] = {0.0}; 
  double flux_rhouy_r[2] = {0.0}; 
  double flux_rhouz_quad_l[2] = {0.0}; 
  double flux_rhouz_l[2] = {0.0}; 
  double flux_rhouz_quad_r[2] = {0.0}; 
  double flux_rhouz_r[2] = {0.0}; 
  double flux_energy_quad_l[2] = {0.0}; 
  double flux_energy_l[2] = {0.0}; 
  double flux_energy_quad_r[2] = {0.0}; 
  double flux_energy_r[2] = {0.0}; 

  q_lr[0] = ser_2x_p1_surfx2_eval_quad_node_0_r(rho_l); 
  q_lr[1] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_l); 
  q_lr[2] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_l); 
  q_lr[3] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_l); 
  q_lr[4] = ser_2x_p1_surfx2_eval_quad_node_0_r(energy_l); 
  q_cl[0] = ser_2x_p1_surfx2_eval_quad_node_0_l(rho_c); 
  q_cl[1] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_c); 
  q_cl[2] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_c); 
  q_cl[3] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_c); 
  q_cl[4] = ser_2x_p1_surfx2_eval_quad_node_0_l(energy_c); 
  q_cr[0] = ser_2x_p1_surfx2_eval_quad_node_0_r(rho_c); 
  q_cr[1] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_c); 
  q_cr[2] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_c); 
  q_cr[3] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_c); 
  q_cr[4] = ser_2x_p1_surfx2_eval_quad_node_0_r(energy_c); 
  q_rl[0] = ser_2x_p1_surfx2_eval_quad_node_0_l(rho_r); 
  q_rl[1] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_r); 
  q_rl[2] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_r); 
  q_rl[3] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_r); 
  q_rl[4] = ser_2x_p1_surfx2_eval_quad_node_0_l(energy_r); 

  rot_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  rot_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  rot_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  rot_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  for (int i = 0; i < 5; ++i) { 
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
  gkyl_euler_flux(gas_gamma, q_lr_local, fhat_l_local); 
  gkyl_euler_flux(gas_gamma, q_cr_local, fhat_r_local); 
  for (int i = 0; i < 5; ++i) { 
    fhat_l_local[i] = lenr_l*fhat_l_local[i] + amdq_l_local[i]; 
    fhat_r_local[i] = lenr_r*fhat_r_local[i] + amdq_r_local[i]; 
  } 
  rot_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], fhat_l_local, fhat_l); 
  rot_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], fhat_r_local, fhat_r); 
  flux_rho_quad_l[0] = fhat_l[0]; 
  flux_rho_quad_r[0] = fhat_r[0]; 
  flux_rhoux_quad_l[0] = fhat_l[1]; 
  flux_rhoux_quad_r[0] = fhat_r[1]; 
  flux_rhouy_quad_l[0] = fhat_l[2]; 
  flux_rhouy_quad_r[0] = fhat_r[2]; 
  flux_rhouz_quad_l[0] = fhat_l[3]; 
  flux_rhouz_quad_r[0] = fhat_r[3]; 
  flux_energy_quad_l[0] = fhat_l[4]; 
  flux_energy_quad_r[0] = fhat_r[4]; 

  q_lr[0] = ser_2x_p1_surfx2_eval_quad_node_1_r(rho_l); 
  q_lr[1] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_l); 
  q_lr[2] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_l); 
  q_lr[3] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_l); 
  q_lr[4] = ser_2x_p1_surfx2_eval_quad_node_1_r(energy_l); 
  q_cl[0] = ser_2x_p1_surfx2_eval_quad_node_1_l(rho_c); 
  q_cl[1] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_c); 
  q_cl[2] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_c); 
  q_cl[3] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_c); 
  q_cl[4] = ser_2x_p1_surfx2_eval_quad_node_1_l(energy_c); 
  q_cr[0] = ser_2x_p1_surfx2_eval_quad_node_1_r(rho_c); 
  q_cr[1] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_c); 
  q_cr[2] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_c); 
  q_cr[3] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_c); 
  q_cr[4] = ser_2x_p1_surfx2_eval_quad_node_1_r(energy_c); 
  q_rl[0] = ser_2x_p1_surfx2_eval_quad_node_1_l(rho_r); 
  q_rl[1] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_r); 
  q_rl[2] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_r); 
  q_rl[3] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_r); 
  q_rl[4] = ser_2x_p1_surfx2_eval_quad_node_1_l(energy_r); 

  rot_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_lr, q_lr_local); 
  rot_to_local(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], q_cl, q_cl_local); 
  rot_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_cr, q_cr_local); 
  rot_to_local(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], q_rl, q_rl_local); 

  for (int i = 0; i < 5; ++i) { 
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
  gkyl_euler_flux(gas_gamma, q_lr_local, fhat_l_local); 
  gkyl_euler_flux(gas_gamma, q_cr_local, fhat_r_local); 
  for (int i = 0; i < 5; ++i) { 
    fhat_l_local[i] = lenr_l*fhat_l_local[i] + amdq_l_local[i]; 
    fhat_r_local[i] = lenr_r*fhat_r_local[i] + amdq_r_local[i]; 
  } 
  rot_to_global(wv_eqn, geom_l->tau1[1], geom_l->tau2[1], geom_l->norm[1], fhat_l_local, fhat_l); 
  rot_to_global(wv_eqn, geom_r->tau1[1], geom_r->tau2[1], geom_r->norm[1], fhat_r_local, fhat_r); 
  flux_rho_quad_l[1] = fhat_l[0]; 
  flux_rho_quad_r[1] = fhat_r[0]; 
  flux_rhoux_quad_l[1] = fhat_l[1]; 
  flux_rhoux_quad_r[1] = fhat_r[1]; 
  flux_rhouy_quad_l[1] = fhat_l[2]; 
  flux_rhouy_quad_r[1] = fhat_r[2]; 
  flux_rhouz_quad_l[1] = fhat_l[3]; 
  flux_rhouz_quad_r[1] = fhat_r[3]; 
  flux_energy_quad_l[1] = fhat_l[4]; 
  flux_energy_quad_r[1] = fhat_r[4]; 

  ser_2x_p1_upwind_quad_to_modal(flux_rho_quad_l, flux_rho_l); 
  ser_2x_p1_upwind_quad_to_modal(flux_rho_quad_r, flux_rho_r); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhoux_quad_l, flux_rhoux_l); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhoux_quad_r, flux_rhoux_r); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouy_quad_l, flux_rhouy_l); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouy_quad_r, flux_rhouy_r); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouz_quad_l, flux_rhouz_l); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouz_quad_r, flux_rhouz_r); 
  ser_2x_p1_upwind_quad_to_modal(flux_energy_quad_l, flux_energy_l); 
  ser_2x_p1_upwind_quad_to_modal(flux_energy_quad_r, flux_energy_r); 

  outrho[0] += (0.7071067811865475*flux_rho_l[0]-0.7071067811865475*flux_rho_r[0])*dx1; 
  outrho[1] += (0.7071067811865475*flux_rho_l[1]-0.7071067811865475*flux_rho_r[1])*dx1; 
  outrho[2] += -(1.224744871391589*(flux_rho_r[0]+flux_rho_l[0])*dx1); 
  outrho[3] += -(1.224744871391589*(flux_rho_r[1]+flux_rho_l[1])*dx1); 

  outrhoux[0] += (0.7071067811865475*flux_rhoux_l[0]-0.7071067811865475*flux_rhoux_r[0])*dx1; 
  outrhoux[1] += (0.7071067811865475*flux_rhoux_l[1]-0.7071067811865475*flux_rhoux_r[1])*dx1; 
  outrhoux[2] += -(1.224744871391589*(flux_rhoux_r[0]+flux_rhoux_l[0])*dx1); 
  outrhoux[3] += -(1.224744871391589*(flux_rhoux_r[1]+flux_rhoux_l[1])*dx1); 

  outrhouy[0] += (0.7071067811865475*flux_rhouy_l[0]-0.7071067811865475*flux_rhouy_r[0])*dx1; 
  outrhouy[1] += (0.7071067811865475*flux_rhouy_l[1]-0.7071067811865475*flux_rhouy_r[1])*dx1; 
  outrhouy[2] += -(1.224744871391589*(flux_rhouy_r[0]+flux_rhouy_l[0])*dx1); 
  outrhouy[3] += -(1.224744871391589*(flux_rhouy_r[1]+flux_rhouy_l[1])*dx1); 

  outrhouz[0] += (0.7071067811865475*flux_rhouz_l[0]-0.7071067811865475*flux_rhouz_r[0])*dx1; 
  outrhouz[1] += (0.7071067811865475*flux_rhouz_l[1]-0.7071067811865475*flux_rhouz_r[1])*dx1; 
  outrhouz[2] += -(1.224744871391589*(flux_rhouz_r[0]+flux_rhouz_l[0])*dx1); 
  outrhouz[3] += -(1.224744871391589*(flux_rhouz_r[1]+flux_rhouz_l[1])*dx1); 

  outenergy[0] += (0.7071067811865475*flux_energy_l[0]-0.7071067811865475*flux_energy_r[0])*dx1; 
  outenergy[1] += (0.7071067811865475*flux_energy_l[1]-0.7071067811865475*flux_energy_r[1])*dx1; 
  outenergy[2] += -(1.224744871391589*(flux_energy_r[0]+flux_energy_l[0])*dx1); 
  outenergy[3] += -(1.224744871391589*(flux_energy_r[1]+flux_energy_l[1])*dx1); 

  return 0.;

} 
