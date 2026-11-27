#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
#include <gkyl_basis_tensor_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double iso_euler_boundary_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos,
    const struct gkyl_wv_eqn *wv_eqn,
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge,
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing (computational).
  // jacob_pos:  Position-map Jacobian per configuration direction of the boundary cell (cdim
  //             blocks of 3 equal entries): its physical size in direction d is J_d dxv[d].
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
  const double *rhoux_edge = &fluid_edge[9]; 
  const double *rhouy_edge = &fluid_edge[18]; 
  const double *rhouz_edge = &fluid_edge[27]; 

  const double *rho_skin = &fluid_skin[0]; 
  const double *rhoux_skin = &fluid_skin[9]; 
  const double *rhouy_skin = &fluid_skin[18]; 
  const double *rhouz_skin = &fluid_skin[27]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[9]; 
  double *outrhouy = &out[18]; 
  double *outrhouz = &out[27]; 

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
  double flux_rho_quad_l[3] = {0.0}; 
  double flux_rho_l[3] = {0.0}; 
  double flux_rho_quad_r[3] = {0.0}; 
  double flux_rho_r[3] = {0.0}; 
  double flux_rhoux_quad_l[3] = {0.0}; 
  double flux_rhoux_l[3] = {0.0}; 
  double flux_rhoux_quad_r[3] = {0.0}; 
  double flux_rhoux_r[3] = {0.0}; 
  double flux_rhouy_quad_l[3] = {0.0}; 
  double flux_rhouy_l[3] = {0.0}; 
  double flux_rhouy_quad_r[3] = {0.0}; 
  double flux_rhouy_r[3] = {0.0}; 
  double flux_rhouz_quad_l[3] = {0.0}; 
  double flux_rhouz_l[3] = {0.0}; 
  double flux_rhouz_quad_r[3] = {0.0}; 
  double flux_rhouz_r[3] = {0.0}; 

  if (edge == -1) { 
    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rho_skin); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rho_edge); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhoux_skin); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhoux_edge); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhouy_skin); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhouy_edge); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhouz_skin); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhouz_edge); 
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
    flux_rho_quad_r[0] = fhat[0]; 
    flux_rhoux_quad_r[0] = fhat[1]; 
    flux_rhouy_quad_r[0] = fhat[2]; 
    flux_rhouz_quad_r[0] = fhat[3]; 

    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rho_skin); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rho_edge); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhoux_skin); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhoux_edge); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhouy_skin); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhouy_edge); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhouz_skin); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhouz_edge); 
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
    flux_rho_quad_r[1] = fhat[0]; 
    flux_rhoux_quad_r[1] = fhat[1]; 
    flux_rhouy_quad_r[1] = fhat[2]; 
    flux_rhouz_quad_r[1] = fhat[3]; 

    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rho_skin); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rho_edge); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhoux_skin); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhoux_edge); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhouy_skin); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhouy_edge); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhouz_skin); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhouz_edge); 
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
    flux_rho_quad_r[2] = fhat[0]; 
    flux_rhoux_quad_r[2] = fhat[1]; 
    flux_rhouy_quad_r[2] = fhat[2]; 
    flux_rhouz_quad_r[2] = fhat[3]; 

    tensor_2x_p2_upwind_quad_to_modal(flux_rho_quad_r, flux_rho_r); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhoux_quad_r, flux_rhoux_r); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhouy_quad_r, flux_rhouy_r); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhouz_quad_r, flux_rhouz_r); 
    flux_rho_l[0] = flux_wall[0]; 
    flux_rho_l[1] = flux_wall[1]; 
    flux_rho_l[2] = flux_wall[2]; 
    flux_rhoux_l[0] = flux_wall[3]; 
    flux_rhoux_l[1] = flux_wall[4]; 
    flux_rhoux_l[2] = flux_wall[5]; 
    flux_rhouy_l[0] = flux_wall[6]; 
    flux_rhouy_l[1] = flux_wall[7]; 
    flux_rhouy_l[2] = flux_wall[8]; 
    flux_rhouz_l[0] = flux_wall[9]; 
    flux_rhouz_l[1] = flux_wall[10]; 
    flux_rhouz_l[2] = flux_wall[11]; 
  } else { 
    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rho_edge); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rho_skin); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhoux_edge); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhoux_skin); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhouy_edge); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhouy_skin); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rhouz_edge); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rhouz_skin); 
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
    flux_rho_quad_l[0] = fhat[0]; 
    flux_rhoux_quad_l[0] = fhat[1]; 
    flux_rhouy_quad_l[0] = fhat[2]; 
    flux_rhouz_quad_l[0] = fhat[3]; 

    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rho_edge); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rho_skin); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhoux_edge); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhoux_skin); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhouy_edge); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhouy_skin); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rhouz_edge); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rhouz_skin); 
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
    flux_rho_quad_l[1] = fhat[0]; 
    flux_rhoux_quad_l[1] = fhat[1]; 
    flux_rhouy_quad_l[1] = fhat[2]; 
    flux_rhouz_quad_l[1] = fhat[3]; 

    q_lft[0] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rho_edge); 
    q_rgt[0] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rho_skin); 
    q_lft[1] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhoux_edge); 
    q_rgt[1] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhoux_skin); 
    q_lft[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhouy_edge); 
    q_rgt[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhouy_skin); 
    q_lft[3] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rhouz_edge); 
    q_rgt[3] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rhouz_skin); 
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
    flux_rho_quad_l[2] = fhat[0]; 
    flux_rhoux_quad_l[2] = fhat[1]; 
    flux_rhouy_quad_l[2] = fhat[2]; 
    flux_rhouz_quad_l[2] = fhat[3]; 

    tensor_2x_p2_upwind_quad_to_modal(flux_rho_quad_l, flux_rho_l); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhoux_quad_l, flux_rhoux_l); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhouy_quad_l, flux_rhouy_l); 
    tensor_2x_p2_upwind_quad_to_modal(flux_rhouz_quad_l, flux_rhouz_l); 
    flux_rho_r[0] = flux_wall[0]; 
    flux_rho_r[1] = flux_wall[1]; 
    flux_rho_r[2] = flux_wall[2]; 
    flux_rhoux_r[0] = flux_wall[3]; 
    flux_rhoux_r[1] = flux_wall[4]; 
    flux_rhoux_r[2] = flux_wall[5]; 
    flux_rhouy_r[0] = flux_wall[6]; 
    flux_rhouy_r[1] = flux_wall[7]; 
    flux_rhouy_r[2] = flux_wall[8]; 
    flux_rhouz_r[0] = flux_wall[9]; 
    flux_rhouz_r[1] = flux_wall[10]; 
    flux_rhouz_r[2] = flux_wall[11]; 
  } 

  outrho[0] += (0.7071067811865475*flux_rho_l[0]-0.7071067811865475*flux_rho_r[0])*dx1; 
  outrho[1] += -(1.224744871391589*(flux_rho_r[0]+flux_rho_l[0])*dx1); 
  outrho[2] += (0.7071067811865475*flux_rho_l[1]-0.7071067811865475*flux_rho_r[1])*dx1; 
  outrho[3] += -(1.224744871391589*(flux_rho_r[1]+flux_rho_l[1])*dx1); 
  outrho[4] += (1.5811388300841895*flux_rho_l[0]-1.5811388300841895*flux_rho_r[0])*dx1; 
  outrho[5] += (0.7071067811865475*flux_rho_l[2]-0.7071067811865475*flux_rho_r[2])*dx1; 
  outrho[6] += (1.5811388300841898*flux_rho_l[1]-1.5811388300841898*flux_rho_r[1])*dx1; 
  outrho[7] += -(1.224744871391589*(flux_rho_r[2]+flux_rho_l[2])*dx1); 
  outrho[8] += (1.5811388300841895*flux_rho_l[2]-1.5811388300841895*flux_rho_r[2])*dx1; 

  outrhoux[0] += (0.7071067811865475*flux_rhoux_l[0]-0.7071067811865475*flux_rhoux_r[0])*dx1; 
  outrhoux[1] += -(1.224744871391589*(flux_rhoux_r[0]+flux_rhoux_l[0])*dx1); 
  outrhoux[2] += (0.7071067811865475*flux_rhoux_l[1]-0.7071067811865475*flux_rhoux_r[1])*dx1; 
  outrhoux[3] += -(1.224744871391589*(flux_rhoux_r[1]+flux_rhoux_l[1])*dx1); 
  outrhoux[4] += (1.5811388300841895*flux_rhoux_l[0]-1.5811388300841895*flux_rhoux_r[0])*dx1; 
  outrhoux[5] += (0.7071067811865475*flux_rhoux_l[2]-0.7071067811865475*flux_rhoux_r[2])*dx1; 
  outrhoux[6] += (1.5811388300841898*flux_rhoux_l[1]-1.5811388300841898*flux_rhoux_r[1])*dx1; 
  outrhoux[7] += -(1.224744871391589*(flux_rhoux_r[2]+flux_rhoux_l[2])*dx1); 
  outrhoux[8] += (1.5811388300841895*flux_rhoux_l[2]-1.5811388300841895*flux_rhoux_r[2])*dx1; 

  outrhouy[0] += (0.7071067811865475*flux_rhouy_l[0]-0.7071067811865475*flux_rhouy_r[0])*dx1; 
  outrhouy[1] += -(1.224744871391589*(flux_rhouy_r[0]+flux_rhouy_l[0])*dx1); 
  outrhouy[2] += (0.7071067811865475*flux_rhouy_l[1]-0.7071067811865475*flux_rhouy_r[1])*dx1; 
  outrhouy[3] += -(1.224744871391589*(flux_rhouy_r[1]+flux_rhouy_l[1])*dx1); 
  outrhouy[4] += (1.5811388300841895*flux_rhouy_l[0]-1.5811388300841895*flux_rhouy_r[0])*dx1; 
  outrhouy[5] += (0.7071067811865475*flux_rhouy_l[2]-0.7071067811865475*flux_rhouy_r[2])*dx1; 
  outrhouy[6] += (1.5811388300841898*flux_rhouy_l[1]-1.5811388300841898*flux_rhouy_r[1])*dx1; 
  outrhouy[7] += -(1.224744871391589*(flux_rhouy_r[2]+flux_rhouy_l[2])*dx1); 
  outrhouy[8] += (1.5811388300841895*flux_rhouy_l[2]-1.5811388300841895*flux_rhouy_r[2])*dx1; 

  outrhouz[0] += (0.7071067811865475*flux_rhouz_l[0]-0.7071067811865475*flux_rhouz_r[0])*dx1; 
  outrhouz[1] += -(1.224744871391589*(flux_rhouz_r[0]+flux_rhouz_l[0])*dx1); 
  outrhouz[2] += (0.7071067811865475*flux_rhouz_l[1]-0.7071067811865475*flux_rhouz_r[1])*dx1; 
  outrhouz[3] += -(1.224744871391589*(flux_rhouz_r[1]+flux_rhouz_l[1])*dx1); 
  outrhouz[4] += (1.5811388300841895*flux_rhouz_l[0]-1.5811388300841895*flux_rhouz_r[0])*dx1; 
  outrhouz[5] += (0.7071067811865475*flux_rhouz_l[2]-0.7071067811865475*flux_rhouz_r[2])*dx1; 
  outrhouz[6] += (1.5811388300841898*flux_rhouz_l[1]-1.5811388300841898*flux_rhouz_r[1])*dx1; 
  outrhouz[7] += -(1.224744871391589*(flux_rhouz_r[2]+flux_rhouz_l[2])*dx1); 
  outrhouz[8] += (1.5811388300841895*flux_rhouz_l[2]-1.5811388300841895*flux_rhouz_r[2])*dx1; 

  return 0.;

} 
