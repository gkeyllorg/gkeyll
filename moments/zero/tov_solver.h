#pragma once

#include <stdbool.h>

#define TOV_MAX_POINTS 500000

// Output struct for BL evaluation
struct tov_eval_bl {
    double rho;
    double P;
    double eps;
    double e;
    double m;
    double Phi;
};

// Output struct for fluid evaluation (CKS/SKS)
struct tov_eval_fluid {
    double rho;
    double P;
    double eps;
    double e;
};

// TOV solution object
struct gkyl_tov;

/**
 * Create a new TOV solution object.
 * Integrates the TOV equations for a polytropic EOS P = K * rho^Gamma.
 *
 * @param K Polytropic constant
 * @param Gamma Adiabatic index
 * @param rho_c Central density
 * @param dr Radial step size
 * @return TOV solution object
 */
struct gkyl_tov* gkyl_tov_new(double K, double Gamma, double rho_c, double dr);


/**
 * Evaluate TOV solution in Boyer-Lindquist coordinates at areal radius r.
 *
 * @param tov TOV solution object
 * @param r Areal radius
 * @param out Output struct
 * @return true on success
 */
bool gkyl_tov_eval_bl(const struct gkyl_tov *tov, double r, struct tov_eval_bl *out);

/** Evaluate m' and m'' from the TOV equations and polytropic EOS, not by
 * differentiating the piecewise-linear radial table. */
bool gkyl_tov_eval_mass_derivs(const struct gkyl_tov *tov, double r,
    double *dm, double *d2m);

/**
 * Evaluate fluid variables at position (x,y,z) relative to star center.
 * r = sqrt(x^2+y^2+z^2) is the areal radius shared between Schwarzschild and CKS.
 *
 * @param tov TOV solution object
 * @param x, y, z Position relative to star center
 * @param out Output struct
 * @return true on success
 */
bool gkyl_tov_eval(const struct gkyl_tov *tov, double x, double y, double z,
    struct tov_eval_fluid *out);

/**
 * Get gravitational mass of the star.
 */
double gkyl_tov_star_mass(const struct gkyl_tov *tov);

/**
 * Get areal radius of the stellar surface.
 */
double gkyl_tov_star_radius(const struct gkyl_tov *tov);

/**
 * Release TOV solution object and free memory.
 */
void gkyl_tov_solution_release(struct gkyl_tov *tov);

/**
 * Refresh lapse potential Phi and enclosed mass m from the current
 * one-dimensional GR TOV state.
 *
 * The state q is assumed to be contiguous with 8 components per cell.
 */
void gkyl_gr_tov_refresh_geometry_from_state(double gas_gamma, double p_atm, int ncells, double *q);
