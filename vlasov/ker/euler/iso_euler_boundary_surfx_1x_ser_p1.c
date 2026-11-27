#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
GKYL_CU_DH double iso_euler_boundary_surfx_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
    const struct gkyl_wv_eqn *wv_eqn,
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge,
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing (computational).
  // jacob_pos:  Position-map Jacobian per configuration direction of the boundary cell (cdim
  //             blocks of 2 equal entries): its physical size in direction d is J_d dxv[d].
  // wv_eqn:     Wave equation: the Riemann solver for the flux at the interior face.
  // geom_l:     Geometry of the lower face of the boundary cell.
  // geom_r:     Geometry of the upper face of the boundary cell.
  // edge:       -1 if the wall is the lower face of the boundary cell, +1 if its upper face.
  // flux_wall:  Surface expansion of the flux through the wall (iso_euler_wall_flux).
  // fluid_edge/skin: [rho, rho ux, rho uy, rho uz], Fluid state of the interior neighbor/the boundary cell.
  // out:        Incremented output of the boundary cell.

  const double dx1 = 2.0/(dxv[0]*jacob_pos[0]); 
  const double vt = container_of(wv_eqn, struct wv_iso_euler, eqn)->vt; 

  const double *rho_edge = &fluid_edge[0]; 
  const double *rhoux_edge = &fluid_edge[2]; 
  const double *rhouy_edge = &fluid_edge[4]; 
  const double *rhouz_edge = &fluid_edge[6]; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[2]; 
  const double *rhouy_skin = &fluid_skin[4]; 
  const double *rhouz_skin = &fluid_skin[6]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[2]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[6]; 

  double q_lft[4] = {0.0}; 
  double q_lft_local[4] = {0.0}; 
  double q_rgt[4] = {0.0}; 
  double q_rgt_local[4] = {0.0}; 
  double delta[4] = {0.0}; 
  double waves[12] = {0.0}; 
  double speeds[3] = {0.0}; 
  double amdq_local[4] = {0.0}; 
  double apdq_local[4] = {0.0}; 
  double fhat_local[4] = {0.0}; 
  double fhat[4] = {0.0}; 
  double flux_rho_l[1] = {0.0}; 
  double flux_rho_r[1] = {0.0}; 
  double flux_rhoux_l[1] = {0.0}; 
  double flux_rhoux_r[1] = {0.0}; 
  double flux_rhouy_l[1] = {0.0}; 
  double flux_rhouy_r[1] = {0.0}; 
  double flux_rhouz_l[1] = {0.0}; 
  double flux_rhouz_r[1] = {0.0}; 

  if (edge == -1) { 
    q_lft[0] = 1.224744871391589*rho_skin[1]+0.7071067811865475*rho_skin[0]; 
    q_rgt[0] = 0.7071067811865475*rho_edge[0]-1.224744871391589*rho_edge[1]; 
    q_lft[1] = 1.224744871391589*rhoux_skin[1]+0.7071067811865475*rhoux_skin[0]; 
    q_rgt[1] = 0.7071067811865475*rhoux_edge[0]-1.224744871391589*rhoux_edge[1]; 
    q_lft[2] = 1.224744871391589*rhouy_skin[1]+0.7071067811865475*rhouy_skin[0]; 
    q_rgt[2] = 0.7071067811865475*rhouy_edge[0]-1.224744871391589*rhouy_edge[1]; 
    q_lft[3] = 1.224744871391589*rhouz_skin[1]+0.7071067811865475*rhouz_skin[0]; 
    q_rgt[3] = 0.7071067811865475*rhouz_edge[0]-1.224744871391589*rhouz_edge[1]; 
    iso_euler_rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_lft, q_lft_local); 
    iso_euler_rot_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rgt, q_rgt_local); 
    for (int i = 0; i < 4; ++i) { 
      delta[i] = q_rgt_local[i] - q_lft_local[i]; 
    } 
    gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q_lft_local, q_rgt_local, 1.0, 1.0, waves, speeds); 
    for (int i = 0; i < 3; ++i) { 
      speeds[i] *= geom_r->lenr[0]; 
    } 
    gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lft_local, q_rgt_local, 1.0, 1.0, waves, speeds, amdq_local, apdq_local); 
    iso_euler_flux(vt, q_lft_local, fhat_local); 
    for (int i = 0; i < 4; ++i) { 
      fhat_local[i] = geom_r->lenr[0]*fhat_local[i] + amdq_local[i]; 
    } 
    iso_euler_rot_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], fhat_local, fhat); 
    flux_rho_r[0] = fhat[0]; 
    flux_rhoux_r[0] = fhat[1]; 
    flux_rhouy_r[0] = fhat[2]; 
    flux_rhouz_r[0] = fhat[3]; 

    flux_rho_l[0] = flux_wall[0]; 
    flux_rhoux_l[0] = flux_wall[1]; 
    flux_rhouy_l[0] = flux_wall[2]; 
    flux_rhouz_l[0] = flux_wall[3]; 
  } else { 
    q_lft[0] = 1.224744871391589*rho_edge[1]+0.7071067811865475*rho_edge[0]; 
    q_rgt[0] = 0.7071067811865475*rho_skin[0]-1.224744871391589*rho_skin[1]; 
    q_lft[1] = 1.224744871391589*rhoux_edge[1]+0.7071067811865475*rhoux_edge[0]; 
    q_rgt[1] = 0.7071067811865475*rhoux_skin[0]-1.224744871391589*rhoux_skin[1]; 
    q_lft[2] = 1.224744871391589*rhouy_edge[1]+0.7071067811865475*rhouy_edge[0]; 
    q_rgt[2] = 0.7071067811865475*rhouy_skin[0]-1.224744871391589*rhouy_skin[1]; 
    q_lft[3] = 1.224744871391589*rhouz_edge[1]+0.7071067811865475*rhouz_edge[0]; 
    q_rgt[3] = 0.7071067811865475*rhouz_skin[0]-1.224744871391589*rhouz_skin[1]; 
    iso_euler_rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lft, q_lft_local); 
    iso_euler_rot_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_rgt, q_rgt_local); 
    for (int i = 0; i < 4; ++i) { 
      delta[i] = q_rgt_local[i] - q_lft_local[i]; 
    } 
    gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta, q_lft_local, q_rgt_local, 1.0, 1.0, waves, speeds); 
    for (int i = 0; i < 3; ++i) { 
      speeds[i] *= geom_l->lenr[0]; 
    } 
    gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lft_local, q_rgt_local, 1.0, 1.0, waves, speeds, amdq_local, apdq_local); 
    iso_euler_flux(vt, q_lft_local, fhat_local); 
    for (int i = 0; i < 4; ++i) { 
      fhat_local[i] = geom_l->lenr[0]*fhat_local[i] + amdq_local[i]; 
    } 
    iso_euler_rot_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], fhat_local, fhat); 
    flux_rho_l[0] = fhat[0]; 
    flux_rhoux_l[0] = fhat[1]; 
    flux_rhouy_l[0] = fhat[2]; 
    flux_rhouz_l[0] = fhat[3]; 

    flux_rho_r[0] = flux_wall[0]; 
    flux_rhoux_r[0] = flux_wall[1]; 
    flux_rhouy_r[0] = flux_wall[2]; 
    flux_rhouz_r[0] = flux_wall[3]; 
  } 

  outrho[0] += (0.7071067811865475*flux_rho_l[0]-0.7071067811865475*flux_rho_r[0])*dx1; 
  outrho[1] += -(1.224744871391589*(flux_rho_r[0]+flux_rho_l[0])*dx1); 

  outrhoux[0] += (0.7071067811865475*flux_rhoux_l[0]-0.7071067811865475*flux_rhoux_r[0])*dx1; 
  outrhoux[1] += -(1.224744871391589*(flux_rhoux_r[0]+flux_rhoux_l[0])*dx1); 

  outrhouy[0] += (0.7071067811865475*flux_rhouy_l[0]-0.7071067811865475*flux_rhouy_r[0])*dx1; 
  outrhouy[1] += -(1.224744871391589*(flux_rhouy_r[0]+flux_rhouy_l[0])*dx1); 

  outrhouz[0] += (0.7071067811865475*flux_rhouz_l[0]-0.7071067811865475*flux_rhouz_r[0])*dx1; 
  outrhouz[1] += -(1.224744871391589*(flux_rhouz_r[0]+flux_rhouz_l[0])*dx1); 

  return 0.;

} 
