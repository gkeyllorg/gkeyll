#include <gkyl_mat.h> 
#include <gkyl_fluid_em_coupling_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_copy_1x_ser_p1(int count, int num_species, double qbym[GKYL_MAX_SPECIES],
  bool has_energy[GKYL_MAX_SPECIES], double epsilon0, double dt, const double *sigma,
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:       index of this cell; node n has the solution 2*count + n of the batch. 
  // has_energy:  whether the state of each species carries an energy (Euler; not isothermal Euler). 
  // dt, sigma:   time step and resistive layer (NULL: none), as in the set kernel. 
  // x:           solution [Mbar_1, ..., Mbar_S, epsilon0 Ebar] per node. 
  // fluid:       fluid state [rho, rho ux, rho uy, rho uz, E] per species: old time in, new time out 
  //              (M^{n+1} = 2 Mbar - M^n; E, where the state has one, gains the change of the bulk 
  //              kinetic energy |M|^2/(2 rho) node by node, the work of the field on the fluid; rho 
  //              is unchanged). 
  // em:          [Ex, Ey, Ez, Bx, By, Bz]: E at the old time in, new time out 
  //              (E^{n+1} = (phi/phi2) Ebar - (phi^2/phi2 - g2) E^n, 2 Ebar - E^n without a layer). 
  const int rE = 3*num_species; 
  double sigq[2], cEbar[2], cE[2], q[2], qnew[2]; 
  for (int n = 0; n < 2; ++n) sigq[n] = sigma ? sigma[n] : 0.0; 
  for (int n = 0; n < 2; ++n) { 
    double g2, phi, phi2; 
    fluid_em_coupling_layer_coeffs(sigq[n]*dt, &g2, &phi, &phi2); 
    cEbar[n] = phi/phi2; 
    cE[n] = phi*phi/phi2 - g2; 
  } 
  for (int c = 0; c < 3; ++c) { 
    fluid_em_coupling_1x_ser_p1_nodal(&em[c*2], q); 
    for (int n = 0; n < 2; ++n) { 
      const struct gkyl_mat sol = gkyl_nmat_get(x, 2*count + n); 
      qnew[n] = cEbar[n]*gkyl_mat_get(&sol, rE + c, 0)/epsilon0 - cE[n]*q[n]; 
    } 
    fluid_em_coupling_1x_ser_p1_modal(qnew, &em[c*2]); 
  } 
  double rhoq[2], dke[2]; 
  for (int s = 0; s < num_species; ++s) { 
    fluid_em_coupling_1x_ser_p1_nodal(fluid[s], rhoq); 
    for (int n = 0; n < 2; ++n) dke[n] = 0.0; 
    for (int c = 0; c < 3; ++c) { 
      fluid_em_coupling_1x_ser_p1_nodal(&fluid[s][(1 + c)*2], q); 
      for (int n = 0; n < 2; ++n) { 
        const struct gkyl_mat sol = gkyl_nmat_get(x, 2*count + n); 
        qnew[n] = 2.0*gkyl_mat_get(&sol, 3*s + c, 0) - q[n]; 
        dke[n] += qnew[n]*qnew[n] - q[n]*q[n]; 
      } 
      fluid_em_coupling_1x_ser_p1_modal(qnew, &fluid[s][(1 + c)*2]); 
    } 
    // the work of the field on the fluid is the change of its bulk kinetic energy 
    if (has_energy[s]) { 
      fluid_em_coupling_1x_ser_p1_nodal(&fluid[s][4*2], q); 
      for (int n = 0; n < 2; ++n) qnew[n] = q[n] + (rhoq[n] > 0.0 ? 0.5*dke[n]/rhoq[n] : 0.0); 
      fluid_em_coupling_1x_ser_p1_modal(qnew, &fluid[s][4*2]); 
    } 
  } 
} 
