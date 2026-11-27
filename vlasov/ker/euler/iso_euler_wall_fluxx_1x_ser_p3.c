#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
GKYL_CU_DH void iso_euler_wall_fluxx_1x_ser_p3(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom,
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall) 
{ 
  // wv_eqn:     Wave equation (its adiabatic index or thermal speed).
  // geom:       Geometry of the wall face.
  // edge:       -1 if the wall is the lower face of the boundary cell, +1 if its upper face.
  // fluid_skin: [rho, rho ux, rho uy, rho uz], Fluid state of the boundary cell.
  // flux_wall:  Output surface expansion of the flux through the wall, along the coordinate
  //             direction ([F_rho, F_rhoux, F_rhouy, F_rhouz], 1 modes each).

  const double vt = container_of(wv_eqn, struct wv_iso_euler, eqn)->vt; 
  const double sgn = edge; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[4]; 
  const double *rhouy_skin = &fluid_skin[8]; 
  const double *rhouz_skin = &fluid_skin[12]; 

  double q[4] = {0.0}; 
  double q_local[4] = {0.0}; 
  double fhat_local[4] = {0.0}; 
  double fhat[4] = {0.0}; 

  if (edge == -1) { 
    q[0] = -(1.8708286933869707*rho_skin[3])+1.5811388300841895*rho_skin[2]-1.224744871391589*rho_skin[1]+0.7071067811865475*rho_skin[0]; 
    q[1] = -(1.8708286933869707*rhoux_skin[3])+1.5811388300841895*rhoux_skin[2]-1.224744871391589*rhoux_skin[1]+0.7071067811865475*rhoux_skin[0]; 
    q[2] = -(1.8708286933869707*rhouy_skin[3])+1.5811388300841895*rhouy_skin[2]-1.224744871391589*rhouy_skin[1]+0.7071067811865475*rhouy_skin[0]; 
    q[3] = -(1.8708286933869707*rhouz_skin[3])+1.5811388300841895*rhouz_skin[2]-1.224744871391589*rhouz_skin[1]+0.7071067811865475*rhouz_skin[0]; 
  } else { 
    q[0] = 1.8708286933869707*rho_skin[3]+1.5811388300841895*rho_skin[2]+1.224744871391589*rho_skin[1]+0.7071067811865475*rho_skin[0]; 
    q[1] = 1.8708286933869707*rhoux_skin[3]+1.5811388300841895*rhoux_skin[2]+1.224744871391589*rhoux_skin[1]+0.7071067811865475*rhoux_skin[0]; 
    q[2] = 1.8708286933869707*rhouy_skin[3]+1.5811388300841895*rhouy_skin[2]+1.224744871391589*rhouy_skin[1]+0.7071067811865475*rhouy_skin[0]; 
    q[3] = 1.8708286933869707*rhouz_skin[3]+1.5811388300841895*rhouz_skin[2]+1.224744871391589*rhouz_skin[1]+0.7071067811865475*rhouz_skin[0]; 
  } 
  iso_euler_rot_to_local(wv_eqn, geom->tau1[0], geom->tau2[0], geom->norm[0], q, q_local); 
  iso_euler_wall_flux(vt, sgn, q_local, fhat_local); 
  iso_euler_rot_to_global(wv_eqn, geom->tau1[0], geom->tau2[0], geom->norm[0], fhat_local, fhat); 
  flux_wall[0] = fhat[0]; 
  flux_wall[1] = fhat[1]; 
  flux_wall[2] = fhat[2]; 
  flux_wall[3] = fhat[3]; 

} 
