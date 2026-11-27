#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
#include <gkyl_basis_tensor_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void iso_euler_wall_fluxy_2x_tensor_p2(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall) 
{ 
  // wv_eqn:     Wave equation (its adiabatic index or thermal speed).
  // geom:       Geometry of the wall face.
  // edge:       -1 if the wall is the lower face of the boundary cell, +1 if its upper face.
  // fluid_skin: [rho, rho ux, rho uy, rho uz], Fluid state of the boundary cell.
  // flux_wall:  Output surface expansion of the flux through the wall, along the coordinate
  //             direction ([F_rho, F_rhoux, F_rhouy, F_rhouz], 3 modes each).

  const double vt = container_of(wv_eqn, struct wv_iso_euler, eqn)->vt; 
  const double sgn = edge; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[9]; 
  const double *rhouy_skin = &fluid_skin[18]; 
  const double *rhouz_skin = &fluid_skin[27]; 

  double q[4] = {0.0}; 
  double q_local[4] = {0.0}; 
  double fhat_local[4] = {0.0}; 
  double fhat[4] = {0.0}; 
  double flux_rho_quad[3] = {0.0}; 
  double flux_rhoux_quad[3] = {0.0}; 
  double flux_rhouy_quad[3] = {0.0}; 
  double flux_rhouz_quad[3] = {0.0}; 

  if (edge == -1) { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_0_l(rhouz_skin); 
  } else { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_0_r(rhouz_skin); 
  } 
  iso_euler_rot_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q, q_local); 
  iso_euler_wall_flux(vt, sgn, q_local, fhat_local); 
  iso_euler_rot_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], fhat_local, fhat); 
  flux_rho_quad[0] = fhat[0]; 
  flux_rhoux_quad[0] = fhat[1]; 
  flux_rhouy_quad[0] = fhat[2]; 
  flux_rhouz_quad[0] = fhat[3]; 

  if (edge == -1) { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_1_l(rhouz_skin); 
  } else { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_1_r(rhouz_skin); 
  } 
  iso_euler_rot_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q, q_local); 
  iso_euler_wall_flux(vt, sgn, q_local, fhat_local); 
  iso_euler_rot_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], fhat_local, fhat); 
  flux_rho_quad[1] = fhat[0]; 
  flux_rhoux_quad[1] = fhat[1]; 
  flux_rhouy_quad[1] = fhat[2]; 
  flux_rhouz_quad[1] = fhat[3]; 

  if (edge == -1) { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_2_l(rhouz_skin); 
  } else { 
    q[0] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rho_skin); 
    q[1] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhoux_skin); 
    q[2] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouy_skin); 
    q[3] = tensor_2x_p2_surfx2_eval_quad_node_2_r(rhouz_skin); 
  } 
  iso_euler_rot_to_local(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], q, q_local); 
  iso_euler_wall_flux(vt, sgn, q_local, fhat_local); 
  iso_euler_rot_to_global(wv_eqn, geom->tau1[1], geom->tau2[1], geom->norm[1], fhat_local, fhat); 
  flux_rho_quad[2] = fhat[0]; 
  flux_rhoux_quad[2] = fhat[1]; 
  flux_rhouy_quad[2] = fhat[2]; 
  flux_rhouz_quad[2] = fhat[3]; 

  tensor_2x_p2_upwind_quad_to_modal(flux_rho_quad, &flux_wall[0]); 
  tensor_2x_p2_upwind_quad_to_modal(flux_rhoux_quad, &flux_wall[3]); 
  tensor_2x_p2_upwind_quad_to_modal(flux_rhouy_quad, &flux_wall[6]); 
  tensor_2x_p2_upwind_quad_to_modal(flux_rhouz_quad, &flux_wall[9]); 

} 
