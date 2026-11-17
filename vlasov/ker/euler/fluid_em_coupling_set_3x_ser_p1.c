#include <gkyl_mat.h> 
#include <gkyl_fluid_em_coupling_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_set_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  double epsilon0, bool is_field_static, double dt, struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n,
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current,
  const double *fluid[GKYL_MAX_SPECIES], const double *em) 
{ 
  // count:             index of the matrix of this cell in the batch. 
  // num_species:       number of species (momentum equations). 
  // qbym:              charge/mass ratio of each species. 
  // epsilon0:          permittivity of free space. 
  // is_field_static:   true if the self-consistent electric field is not advanced. 
  // dt:                time step. 
  // A_n, rhs_n:        batched left-hand side and right-hand side; block diagonal over the 8 nodes, 
  //                    unknowns per node [Mbar_1, ..., Mbar_S, epsilon0 Ebar] (3 S + 3). 
  // app_accel:         applied accelerations [ax, ay, az] per species. 
  // ext_em:            external fields [Ex, Ey, Ez, Bx, By, Bz]. 
  // app_current:       applied current [Jx, Jy, Jz]. 
  // fluid:             fluid state [rho, rho ux, rho uy, rho uz, E] per species at the old time. 
  // em:                [Ex, Ey, Ez, Bx, By, Bz] at the old time. 
  struct gkyl_mat lhs = gkyl_nmat_get(A_n, count); 
  struct gkyl_mat rhs = gkyl_nmat_get(rhs_n, count); 
  gkyl_mat_clear(&lhs, 0.0); 
  gkyl_mat_clear(&rhs, 0.0); 
  const int NU = 3*num_species + 3; 
  double Eq[3][8], Bq[3][8], extEq[3][8], Jq[3][8], tmp[8]; 
  for (int c = 0; c < 3; ++c) { 
    fluid_em_coupling_3x_ser_p1_nodal(&em[c*8], Eq[c]); 
    fluid_em_coupling_3x_ser_p1_nodal(&em[(3 + c)*8], Bq[c]); 
    fluid_em_coupling_3x_ser_p1_nodal(&ext_em[(3 + c)*8], tmp); 
    for (int n = 0; n < 8; ++n) Bq[c][n] += tmp[n]; 
    fluid_em_coupling_3x_ser_p1_nodal(&ext_em[c*8], extEq[c]); 
    fluid_em_coupling_3x_ser_p1_nodal(&app_current[c*8], Jq[c]); 
  } 
  // Ampere's law: epsilon0 Ebar + dt/2 sum_s qbym_s Mbar_s = epsilon0 E^n - dt/2 J_app. 
  for (int n = 0; n < 8; ++n) { 
    for (int c = 0; c < 3; ++c) { 
      const int rE = n*NU + 3*num_species + c; 
      gkyl_mat_set(&lhs, rE, rE, 1.0); 
      gkyl_mat_set(&rhs, rE, 0, epsilon0*Eq[c][n] - 0.5*dt*Jq[c][n]); 
    } 
  } 
  // Momentum: Mbar_s - dt/2 qbym_s (rho_s Ebar + Mbar_s x B) = M_s^n + dt/2 (qbym_s rho_s E_ext + rho_s a_s). 
  for (int s = 0; s < num_species; ++s) { 
    double rhoq[8], Mq[3][8], aq[3][8]; 
    fluid_em_coupling_3x_ser_p1_nodal(fluid[s], rhoq); 
    for (int c = 0; c < 3; ++c) { 
      fluid_em_coupling_3x_ser_p1_nodal(&fluid[s][(1 + c)*8], Mq[c]); 
      fluid_em_coupling_3x_ser_p1_nodal(&app_accel[s][c*8], aq[c]); 
    } 
    const double f = -0.5*dt*qbym[s]; 
    for (int n = 0; n < 8; ++n) { 
      const int r = n*NU + 3*s, rE = n*NU + 3*num_species; 
      for (int c = 0; c < 3; ++c) { 
        gkyl_mat_set(&lhs, r + c, r + c, 1.0); 
        gkyl_mat_set(&lhs, r + c, rE + c, f*rhoq[n]/epsilon0); 
        gkyl_mat_set(&rhs, r + c, 0, Mq[c][n] + 0.5*dt*rhoq[n]*(qbym[s]*extEq[c][n] + aq[c][n])); 
        if (!is_field_static) gkyl_mat_set(&lhs, rE + c, r + c, 0.5*dt*qbym[s]); 
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
