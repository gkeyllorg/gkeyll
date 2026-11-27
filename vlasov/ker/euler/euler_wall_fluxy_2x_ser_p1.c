#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void euler_wall_fluxy_2x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall) 
{ 
  // wv_eqn:     Wave equation (its adiabatic index or thermal speed).
  // geom:       Geometry of the wall face.
  // edge:       -1 if the wall is the lower face of the boundary cell, +1 if its upper face.
  // fluid_skin: [rho, rho ux, rho uy, rho uz, energy], Fluid state of the boundary cell.
  // flux_wall:  Output surface expansion of the flux through the wall, along the coordinate
  //             direction ([F_rho, F_rhoux, F_rhouy, F_rhouz, F_energy], 2 modes each).

  const double gas_gamma = container_of(wv_eqn, struct wv_euler, eqn)->gas_gamma; 
  const double sgn = edge; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[4]; 
  const double *rhouy_skin = &fluid_skin[8]; 
  const double *rhouz_skin = &fluid_skin[12]; 
  const double *energy_skin = &fluid_skin[16]; 

  double q[5] = {0.0}; 
  double q_local[5] = {0.0}; 
  double fhat_local[5] = {0.0}; 
  double fhat[5] = {0.0}; 
  double flux_rho_quad[2] = {0.0}; 
  double flux_rhoux_quad[2] = {0.0}; 
  double flux_rhouy_quad[2] = {0.0}; 
  double flux_rhouz_quad[2] = {0.0}; 
  double flux_energy_quad[2] = {0.0}; 

  if (edge == -1) { 
    q[0] = ser_2x_p1_surfx2_eval_quad_node_0_l(rho_skin); 
    q[1] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhoux_skin); 
    q[2] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouy_skin); 
    q[3] = ser_2x_p1_surfx2_eval_quad_node_0_l(rhouz_skin); 
    q[4] = ser_2x_p1_surfx2_eval_quad_node_0_l(energy_skin); 
  } else { 
    q[0] = ser_2x_p1_surfx2_eval_quad_node_0_r(rho_skin); 
    q[1] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhoux_skin); 
    q[2] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouy_skin); 
    q[3] = ser_2x_p1_surfx2_eval_quad_node_0_r(rhouz_skin); 
    q[4] = ser_2x_p1_surfx2_eval_quad_node_0_r(energy_skin); 
  } 
  rot_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q, q_local); 
  euler_wall_flux(gas_gamma, sgn, q_local, fhat_local); 
  rot_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], fhat_local, fhat); 
  flux_rho_quad[0] = fhat[0]; 
  flux_rhoux_quad[0] = fhat[1]; 
  flux_rhouy_quad[0] = fhat[2]; 
  flux_rhouz_quad[0] = fhat[3]; 
  flux_energy_quad[0] = fhat[4]; 

  if (edge == -1) { 
    q[0] = ser_2x_p1_surfx2_eval_quad_node_1_l(rho_skin); 
    q[1] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhoux_skin); 
    q[2] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouy_skin); 
    q[3] = ser_2x_p1_surfx2_eval_quad_node_1_l(rhouz_skin); 
    q[4] = ser_2x_p1_surfx2_eval_quad_node_1_l(energy_skin); 
  } else { 
    q[0] = ser_2x_p1_surfx2_eval_quad_node_1_r(rho_skin); 
    q[1] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhoux_skin); 
    q[2] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouy_skin); 
    q[3] = ser_2x_p1_surfx2_eval_quad_node_1_r(rhouz_skin); 
    q[4] = ser_2x_p1_surfx2_eval_quad_node_1_r(energy_skin); 
  } 
  rot_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q, q_local); 
  euler_wall_flux(gas_gamma, sgn, q_local, fhat_local); 
  rot_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], fhat_local, fhat); 
  flux_rho_quad[1] = fhat[0]; 
  flux_rhoux_quad[1] = fhat[1]; 
  flux_rhouy_quad[1] = fhat[2]; 
  flux_rhouz_quad[1] = fhat[3]; 
  flux_energy_quad[1] = fhat[4]; 

  ser_2x_p1_upwind_quad_to_modal(flux_rho_quad, &flux_wall[0]); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhoux_quad, &flux_wall[2]); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouy_quad, &flux_wall[4]); 
  ser_2x_p1_upwind_quad_to_modal(flux_rhouz_quad, &flux_wall[6]); 
  ser_2x_p1_upwind_quad_to_modal(flux_energy_quad, &flux_wall[8]); 

} 
