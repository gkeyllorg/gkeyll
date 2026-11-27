#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_wv_euler_priv.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

// Flux through an absorbing wall of the state q = [rho, rho un, rho ut1, rho ut2, E] in the frame 
// of the face: that of the particles of the Maxwellian with the density, flow and temperature 
// T/m = p/rho of q that move toward the wall (none come back). The normal translational degree 
// of freedom is kinetic (half-range moments of the Maxwellian); the rest of the internal energy, 
// T (1/(gas_gamma - 1) - 1/2) per particle, is carried at its mean. sgn = +1 if the wall is on 
// the side the normal points to, -1 if on the other side; the flux is along the normal. 
GKYL_CU_DH static inline void 
euler_wall_flux(double gas_gamma, double sgn, const double *q, double *flux) 
{ 
  const double rho = q[0]; 
  const double un = sgn*q[1]/rho, ut1 = q[2]/rho, ut2 = q[3]/rho; // un: toward the wall 
  const double vtsq = (gas_gamma - 1.0)*(q[4] - 0.5*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])/rho)/rho; 
  const double vt = sqrt(vtsq); 
  const double s = 0.7071067811865476*un/vt; 
  const double hf = 0.5*erfc(-s); // fraction of the particles moving toward the wall 
  const double ex = 0.3989422804014327*vt*exp(-s*s); 
  // half-range moments of the normal velocity per unit mass density 
  const double m1 = un*hf + ex; 
  const double m2 = (un*un + vtsq)*hf + un*ex; 
  const double m3 = un*(un*un + 3.0*vtsq)*hf + (un*un + 2.0*vtsq)*ex; 
  flux[0] = sgn*rho*m1; 
  flux[1] = rho*m2; 
  flux[2] = ut1*flux[0]; 
  flux[3] = ut2*flux[0]; 
  flux[4] = sgn*0.5*rho*m3 + flux[0]*(0.5*(ut1*ut1 + ut2*ut2) + vtsq*(1.0/(gas_gamma - 1.0) - 0.5)); 
} 

GKYL_CU_DH void fluid_vars_pressure_1x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_1x_ser_p1(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_1x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_1x_ser_p1(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_1x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_1x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_positivity_1x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 

GKYL_CU_DH void fluid_vars_pressure_1x_ser_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_1x_ser_p2(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_1x_ser_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_1x_ser_p2(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_1x_ser_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_1x_ser_p2(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_1x_ser_p2(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_1x_ser_p2(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_1x_ser_p3(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_1x_ser_p3(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_1x_ser_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_1x_ser_p3(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_1x_ser_p3(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_1x_ser_p3(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_1x_ser_p3(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_1x_ser_p3(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_1x_ser_p3(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_1x_ser_p3(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_1x_ser_p3(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 

GKYL_CU_DH void fluid_vars_pressure_2x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_2x_ser_p1(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_2x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limitery_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfy_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxy_2x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfy_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_2x_ser_p1(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_2x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_2x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_positivity_2x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 

GKYL_CU_DH void fluid_vars_pressure_3x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_3x_ser_p1(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_3x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limitery_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfy_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxy_3x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfy_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterz_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfz_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxz_3x_ser_p1(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfz_3x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_3x_ser_p1(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_3x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_3x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_positivity_3x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 

GKYL_CU_DH void fluid_vars_pressure_2x_tensor_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_ke_2x_tensor_p2(const double *fluid, const double *u, double* GKYL_RESTRICT ke); 
GKYL_CU_DH int fluid_vars_u_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *fluid); 
GKYL_CU_DH void fluid_vars_u_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf); 
GKYL_CU_DH double euler_vol_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limiterx_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxx_2x_tensor_p2(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfx_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_limitery_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double euler_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_wall_fluxy_2x_tensor_p2(const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, 
    const int edge, const double *fluid_skin, double* GKYL_RESTRICT flux_wall); 
GKYL_CU_DH double euler_boundary_surfy_2x_tensor_p2(const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, const int edge, 
    const double *flux_wall, const double *fluid_edge, const double *fluid_skin, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_integrated_2x_tensor_p2(const double *fluid, const double* u_i, const double* p_ij, double* GKYL_RESTRICT int_fluid_vars); 
GKYL_CU_DH void fluid_vars_source_2x_tensor_p2(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_em_source_2x_tensor_p2(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 

EXTERN_C_END 
