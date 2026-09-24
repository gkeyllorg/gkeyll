#pragma once

/**
 * Compute the fluid pressure given the conserved variables.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
static inline double
gkyl_mhd_pressure(double gas_gamma, const double q[8])
{
  return (gas_gamma-1) *
    (q[4] - 0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]
      - 0.5*(q[5]*q[5]+q[6]*q[6]+q[7]*q[7]));
}

/**
 * Compute MHD fast Alfven speed.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Fast Alfven speed for given q
 */
double gkyl_mhd_fast_speed(double gas_gamma, const double q[8]);

/**
 * Compute maximum absolute speed.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_mhd_max_abs_speed(double gas_gamma, const double q[8]);

/**
 * Compute eigen speeds due to roe.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @param buf Array with four elements to hold computed u, v, w, cfast.
 * @return buf
 */
void gkyl_mhd_eigen_speeds_roe(const double gamma, const double *ql,
    const double *qr, double buf[]);

/**
 * Compute maximum absolute speed using MHD flow and fast Alfven speed due to a
 * Roe scheme.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param ql Conserved variables of the left stage
 * @param qr Conserved variables of the right stage
 * @return Maximum absolute speed for given ql and qr
 */
double gkyl_mhd_max_abs_speed_roe(const double gamma, const double *ql, const double *qr);

/*
 * Compute conserved variables from primitive variables. Prim vars in
 * the order [rho, vx, vy, vz, pr, Bx, By, Bz]
 *
 * @param gas_gamma Gas adiabatic constant
 * @param pv Primitive variables
 * @param q Conserved variables
 */
void gkyl_mhd_cons_vars(double gas_gamma, const double pv[8], double q[8]);

/**
 * Compute flux. Assumes rotation to local coordinate system.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_mhd_flux(double gas_gamma, const double q[8], double flux[8]);

/**
 * Compute flux for the GLM-MHD equations. Assumes rotation to local coordinate system.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param ch The propagation speed of div(B) errors in the GLM method.
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_glm_mhd_flux(double gas_gamma, double ch, const double q[9], double flux[9]);
