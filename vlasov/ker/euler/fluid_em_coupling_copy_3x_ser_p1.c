#include <gkyl_mat.h> 
#include <gkyl_fluid_em_coupling_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_copy_3x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, struct gkyl_nmat *x,
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:       index of the solution of this cell in the batch. 
  // has_energy:  whether the state of each species carries an energy (Euler; not isothermal Euler). 
  // x:           solution [Mbar_1, ..., Mbar_S, epsilon0 Ebar] per node. 
  // fluid:       fluid state [rho, rho ux, rho uy, rho uz, E] per species: old time in, new time out 
  //              (M^{n+1} = 2 Mbar - M^n; E, where the state has one, gains the change of the bulk 
  //              kinetic energy |M|^2/(2 rho) node by node, the work of the field on the fluid; rho 
  //              is unchanged). 
  // em:          [Ex, Ey, Ez, Bx, By, Bz]: E at the old time in, new time out (E^{n+1} = 2 Ebar - E^n). 
  struct gkyl_mat sol = gkyl_nmat_get(x, count); 
  const int NU = 3*num_species + 3; 
  double q[8], qnew[8]; 
  for (int c = 0; c < 3; ++c) { 
    fluid_em_coupling_3x_ser_p1_nodal(&em[c*8], q); 
    for (int n = 0; n < 8; ++n) qnew[n] = 2.0*gkyl_mat_get(&sol, n*NU + 3*num_species + c, 0)/epsilon0 - q[n]; 
    fluid_em_coupling_3x_ser_p1_modal(qnew, &em[c*8]); 
  } 
  double rhoq[8], dke[8]; 
  for (int s = 0; s < num_species; ++s) { 
    fluid_em_coupling_3x_ser_p1_nodal(fluid[s], rhoq); 
    for (int n = 0; n < 8; ++n) dke[n] = 0.0; 
    for (int c = 0; c < 3; ++c) { 
      fluid_em_coupling_3x_ser_p1_nodal(&fluid[s][(1 + c)*8], q); 
      for (int n = 0; n < 8; ++n) { 
        qnew[n] = 2.0*gkyl_mat_get(&sol, n*NU + 3*s + c, 0) - q[n]; 
        dke[n] += qnew[n]*qnew[n] - q[n]*q[n]; 
      } 
      fluid_em_coupling_3x_ser_p1_modal(qnew, &fluid[s][(1 + c)*8]); 
    } 
    // the work of the field on the fluid is the change of its bulk kinetic energy 
    if (has_energy[s]) { 
      fluid_em_coupling_3x_ser_p1_nodal(&fluid[s][4*8], q); 
      for (int n = 0; n < 8; ++n) qnew[n] = q[n] + (rhoq[n] > 0.0 ? 0.5*dke[n]/rhoq[n] : 0.0); 
      fluid_em_coupling_3x_ser_p1_modal(qnew, &fluid[s][4*8]); 
    } 
  } 
} 
