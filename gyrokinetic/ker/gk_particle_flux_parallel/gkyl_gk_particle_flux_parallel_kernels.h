#pragma once
#include <gkyl_util.h>

EXTERN_C_BEG

/**
 * Integrate the surface nodal flux over reference velocity coordinates.
 *
 * @param factor Signed velocity-cell measure and gyroangle prefactor.
 * @param flux Surface nodal flux.
 * @param out Integrated flux at configuration-space surface nodes.
 */
GKYL_CU_DH void gk_particle_flux_parallel_integral_1x1v_ser_p1(double factor, const double *flux, double *out);

/**
 * Integrate the surface nodal flux over reference velocity coordinates.
 *
 * @param factor Signed velocity-cell measure and gyroangle prefactor.
 * @param flux Surface nodal flux.
 * @param out Integrated flux at configuration-space surface nodes.
 */
GKYL_CU_DH void gk_particle_flux_parallel_integral_1x2v_ser_p1(double factor, const double *flux, double *out);

/**
 * Integrate the surface nodal flux over reference velocity coordinates.
 *
 * @param factor Signed velocity-cell measure and gyroangle prefactor.
 * @param flux Surface nodal flux.
 * @param out Integrated flux at configuration-space surface nodes.
 */
GKYL_CU_DH void gk_particle_flux_parallel_integral_2x2v_ser_p1(double factor, const double *flux, double *out);

/**
 * Integrate the surface nodal flux over reference velocity coordinates.
 *
 * @param factor Signed velocity-cell measure and gyroangle prefactor.
 * @param flux Surface nodal flux.
 * @param out Integrated flux at configuration-space surface nodes.
 */
GKYL_CU_DH void gk_particle_flux_parallel_integral_3x2v_ser_p1(double factor, const double *flux, double *out);

EXTERN_C_END
