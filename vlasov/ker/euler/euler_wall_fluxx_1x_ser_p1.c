#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_wall_fluxx_1x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall) 
{ 
  // wv_eqn:     Wave equation (its adiabatic index or thermal speed).
  // geom:       Geometry of the wall face.
  // edge:       -1 if the wall is the lower face of the boundary cell, +1 if its upper face.
  // fluid_skin: [rho, rho ux, rho uy, rho uz, energy], Fluid state of the boundary cell.
  // flux_wall:  Output surface expansion of the flux through the wall, along the coordinate
  //             direction ([F_rho, F_rhoux, F_rhouy, F_rhouz, F_energy], 1 modes each).

  const double gas_gamma = container_of(wv_eqn, struct wv_euler, eqn)->gas_gamma; 
  const double sgn = edge; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[2]; 
  const double *rhouy_skin = &fluid_skin[4]; 
  const double *rhouz_skin = &fluid_skin[6]; 
  const double *energy_skin = &fluid_skin[8]; 

  double q[5] = {0.0}; 
  double q_local[5] = {0.0}; 
  double fhat_local[5] = {0.0}; 
  double fhat[5] = {0.0}; 

  if (edge == -1) { 
    q[0] = 0.7071067811865475*rho_skin[0]-1.224744871391589*rho_skin[1]; 
    q[1] = 0.7071067811865475*rhoux_skin[0]-1.224744871391589*rhoux_skin[1]; 
    q[2] = 0.7071067811865475*rhouy_skin[0]-1.224744871391589*rhouy_skin[1]; 
    q[3] = 0.7071067811865475*rhouz_skin[0]-1.224744871391589*rhouz_skin[1]; 
    q[4] = 0.7071067811865475*energy_skin[0]-1.224744871391589*energy_skin[1]; 
  } else { 
    q[0] = 1.224744871391589*rho_skin[1]+0.7071067811865475*rho_skin[0]; 
    q[1] = 1.224744871391589*rhoux_skin[1]+0.7071067811865475*rhoux_skin[0]; 
    q[2] = 1.224744871391589*rhouy_skin[1]+0.7071067811865475*rhouy_skin[0]; 
    q[3] = 1.224744871391589*rhouz_skin[1]+0.7071067811865475*rhouz_skin[0]; 
    q[4] = 1.224744871391589*energy_skin[1]+0.7071067811865475*energy_skin[0]; 
  } 
  rot_to_local(wv_eqn, geom->tau1[0], geom->tau2[0], geom->norm[0], q, q_local); 
  euler_wall_flux(gas_gamma, sgn, q_local, fhat_local); 
  rot_to_global(wv_eqn, geom->tau1[0], geom->tau2[0], geom->norm[0], fhat_local, fhat); 
  flux_wall[0] = fhat[0]; 
  flux_wall[1] = fhat[1]; 
  flux_wall[2] = fhat[2]; 
  flux_wall[3] = fhat[3]; 
  flux_wall[4] = fhat[4]; 

} 
