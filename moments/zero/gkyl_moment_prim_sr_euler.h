#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Compute primitive variables given conserved variables.
 * 
 * See "Grid-based Methods in Relativistic Hydrodynamics and
 * Magnetohydrodynamics", Marti and Muller, Living Reviews in
 * Comp. Astro. vol 1 (3), 2015. URL:
 * https://link.springer.com/article/10.1007/lrca-2015-3
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param v Primitive variables (output)
 */
static inline void
gkyl_sr_euler_prim_vars(double gas_gamma, const double q[5], double v[5])
{
  double us=0., vs=0., ws=0., q2s = 0., cs2 = 0.;
  double gammas=0., rhos=0., rhoEpss = 0., fs = 0., dfs = 0., fac0 = 1.;
  double g1 = gas_gamma - 1;
  double ps = 0., ps2 = 1.; 
  double tol = 1.e-6;
  size_t iter = 0;
  
  while (fabs(ps2 - ps) > tol) {
    iter += 1;
    ps = ps2;
    fac0 = q[1] + ps;
    us = q[2] / fac0;
    vs = q[3] / fac0;
    ws = q[4] / fac0;
    q2s = us*us + vs*vs + ws*ws;
    //gammas = 1. / sqrt(1. - q2s);
    gammas = pow(fabs(1. - q2s), -0.5);
    
    rhos = q[0] / gammas;
    rhoEpss = (q[1]  - gammas*q[0] + ps*(1 - gammas*gammas)) / (gammas*gammas);
    fs = g1*rhoEpss - ps; // (gas_gamma - 1)*rhos*eps - ps = p(rhos,epss) - ps, eqn 55
    cs2 = gas_gamma*gammas*gammas*ps / fac0;
    dfs = q2s*cs2 - 1; // eqn 60 for df / dp
    ps2 = ps - fs / dfs;

    //printf("---> Iteration %ld (%g) \n", iter, ps2);
    //printf(" %lg %lg %lg %lg %lg\n", q[0], q[1], q[2], q[3], q[4]);
  }
  //printf("Iterations %ld. Error %lg \n", iter, fabs(ps2-ps));
  
  fac0 = q[1] + ps2;
  us = q[2] / fac0;  
  vs = q[3] / fac0;
  ws = q[4] / fac0;
  q2s = us*us + vs*vs + ws*ws;
  //gammas = 1 / sqrt(1 - q2s);
  gammas = pow(fabs(1. - q2s), -0.5);
  
  v[0] = q[0] / gammas; // rho
  v[1] = ps2; // p
  v[2] = us; 
  v[3] = vs;
  v[4] = ws;
}

/**
 * Compute maximum absolute speed.
 *  
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_sr_euler_max_abs_speed(double gas_gamma, const double q[5]);

/**
 * Compute flux. Assumes rotation to local coordinate system.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_sr_euler_flux(double gas_gamma, const double q[5], double flux[5]);
