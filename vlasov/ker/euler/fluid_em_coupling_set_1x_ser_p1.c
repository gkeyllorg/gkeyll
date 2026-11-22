#include <gkyl_mat.h> 
#include <gkyl_fluid_em_coupling_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_set_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  double epsilon0, bool is_field_static, double dt, struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n,
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current,
  const double *sigma, const double *fluid[GKYL_MAX_SPECIES], const double *em) 
{ 
  // count:             index of this cell; node n has the system 2*count + n of the batch. 
  // num_species:       number of species (momentum equations). 
  // qbym:              charge/mass ratio of each species. 
  // epsilon0:          permittivity of free space. 
  // is_field_static:   true if the self-consistent electric field is not advanced. 
  // dt:                time step. 
  // A_n, rhs_n:        batched left-hand sides and right-hand sides, one system per node in the 
  //                    unknowns [Mbar_1, ..., Mbar_S, epsilon0 Ebar] (3 S + 3). 
  // app_accel:         applied accelerations [ax, ay, az] per species. 
  // ext_em:            external fields [Ex, Ey, Ez, Bx, By, Bz]. 
  // app_current:       applied current [Jx, Jy, Jz]. 
  // sigma:             resistive layer, the damping rate of E at the 2 nodes (NULL: none). 
  // fluid:             fluid state [rho, rho ux, rho uy, rho uz, E] per species at the old time. 
  // em:                [Ex, Ey, Ez, Bx, By, Bz] at the old time. 
  const int rE = 3*num_species; 
  double Eq[3][2], Bq[3][2], extEq[3][2], Jq[3][2], sigq[2], phi2q[2], tmp[2]; 
  for (int c = 0; c < 3; ++c) { 
    fluid_em_coupling_1x_ser_p1_nodal(&em[c*2], Eq[c]); 
    fluid_em_coupling_1x_ser_p1_nodal(&em[(3 + c)*2], Bq[c]); 
    fluid_em_coupling_1x_ser_p1_nodal(&ext_em[(3 + c)*2], tmp); 
    for (int n = 0; n < 2; ++n) Bq[c][n] += tmp[n]; 
    fluid_em_coupling_1x_ser_p1_nodal(&ext_em[c*2], extEq[c]); 
    fluid_em_coupling_1x_ser_p1_nodal(&app_current[c*2], Jq[c]); 
  } 
  for (int n = 0; n < 2; ++n) sigq[n] = sigma ? sigma[n] : 0.0; 
  // Ampere's law: epsilon0 Ebar + dt phi2 sum_s qbym_s Mbar_s = epsilon0 phi E^n - dt phi2 J_app. 
  for (int n = 0; n < 2; ++n) { 
    struct gkyl_mat lhs = gkyl_nmat_get(A_n, 2*count + n); 
    struct gkyl_mat rhs = gkyl_nmat_get(rhs_n, 2*count + n); 
    gkyl_mat_clear(&lhs, 0.0); 
    gkyl_mat_clear(&rhs, 0.0); 
    double g2, phi; 
    fluid_em_coupling_layer_coeffs(sigq[n]*dt, &g2, &phi, &phi2q[n]); 
    for (int c = 0; c < 3; ++c) { 
      gkyl_mat_set(&lhs, rE + c, rE + c, 1.0); 
      gkyl_mat_set(&rhs, rE + c, 0, epsilon0*phi*Eq[c][n] - dt*phi2q[n]*Jq[c][n]); 
    } 
  } 
  // Momentum: Mbar_s - dt/2 qbym_s (rho_s Ebar + Mbar_s x B) = M_s^n + dt/2 (qbym_s rho_s E_ext + rho_s a_s). 
  for (int s = 0; s < num_species; ++s) { 
    double rhoq[2], Mq[3][2], aq[3][2]; 
    fluid_em_coupling_1x_ser_p1_nodal(fluid[s], rhoq); 
    for (int c = 0; c < 3; ++c) { 
      fluid_em_coupling_1x_ser_p1_nodal(&fluid[s][(1 + c)*2], Mq[c]); 
      fluid_em_coupling_1x_ser_p1_nodal(&app_accel[s][c*2], aq[c]); 
    } 
    const int r = 3*s; 
    const double f = -0.5*dt*qbym[s]; 
    for (int n = 0; n < 2; ++n) { 
      struct gkyl_mat lhs = gkyl_nmat_get(A_n, 2*count + n); 
      struct gkyl_mat rhs = gkyl_nmat_get(rhs_n, 2*count + n); 
      for (int c = 0; c < 3; ++c) { 
        gkyl_mat_set(&lhs, r + c, r + c, 1.0); 
        gkyl_mat_set(&lhs, r + c, rE + c, f*rhoq[n]/epsilon0); 
        gkyl_mat_set(&rhs, r + c, 0, Mq[c][n] + 0.5*dt*rhoq[n]*(qbym[s]*extEq[c][n] + aq[c][n])); 
        if (!is_field_static) gkyl_mat_set(&lhs, rE + c, r + c, dt*phi2q[n]*qbym[s]); 
      } 
      // -dt/2 qbym (Mbar x B) 
      gkyl_mat_set(&lhs, r + 0, r + 1, f*Bq[2][n]); 
      gkyl_mat_set(&lhs, r + 0, r + 2, -f*Bq[1][n]); 
      gkyl_mat_set(&lhs, r + 1, r + 2, f*Bq[0][n]); 
      gkyl_mat_set(&lhs, r + 1, r + 0, -f*Bq[2][n]); 
      gkyl_mat_set(&lhs, r + 2, r + 0, f*Bq[1][n]); 
      gkyl_mat_set(&lhs, r + 2, r + 1, -f*Bq[0][n]); 
    } 
  } 
} 
