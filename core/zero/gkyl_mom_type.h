#pragma once

#include <gkyl_ref_count.h>
#include <stdint.h>

// Define options for moments of a distribution function.
enum gkyl_distribution_moments {
  GKYL_F_MOMENT_M0 = 0, // Number density.
  GKYL_F_MOMENT_M1, // Momentum density.
  GKYL_F_MOMENT_M2, // Kinetic energy density.
  GKYL_F_MOMENT_M2PAR, // Parallel kinetic energy density.
  GKYL_F_MOMENT_M2PERP, // Perpendicular kinetic energy density.
  GKYL_F_MOMENT_M2IJ, // Kinetic energy tensor..
  GKYL_F_MOMENT_M3, // Heat flux.
  GKYL_F_MOMENT_M3PAR, // Parallel energy flux.
  GKYL_F_MOMENT_M3PERP, // Perpendicular energy flux.
  GKYL_F_MOMENT_M3IJK, // Heat flux in lab frame.
  GKYL_F_MOMENT_MAXWELLIAN, // M0, drift speed, T/m.
  GKYL_F_MOMENT_BIMAXWELLIAN, // M0, drift speed, Tpar/m, Tperp/m.
  GKYL_F_MOMENT_LTE, // Maxwellian or Maxwell-Juttner moments.
  GKYL_F_MOMENT_M0M1M2,  // M0, M1, M2.
  GKYL_F_MOMENT_M0M1M2PARM2PERP,  // M0, M1, M2par, M2perp.
  GKYL_F_MOMENT_HAMILTONIAN,  // M0, mass*M1, H moments.
  GKYL_F_MOMENT_M1_FROM_H, // dH/dv / m moment.
  GKYL_F_MOMENT_ENERGY, // H moment.
  GKYL_F_MOMENT_M0ENERGYM3, // M0, Energy (H) and M3 moments.
  GKYL_F_MOMENT_NI, // M0, M1i for-vector.
  GKYL_F_MOMENT_TIJ, // Stress-energy tensor.
};

// String names corresponding to the enum options above.
static const char *gkyl_distribution_moments_strs[] = {
  "M0",
  "M1",
  "M2",
  "M2par",
  "M2perp",
  "M2ij",
  "M3",
  "M3par",
  "M3perp",
  "M3ijk",
  "MaxwellianMoments",
  "BiMaxwellianMoments",
  "LTEMoments",
  "M0M1M2",
  "M0M1M2parM2perp",
  "HamiltonianMoments",
  "M1_from_H",
  "EnergyMoment",
  "M0EnergyM3",
  "Ni",
  "Tij",
};

// Descriptions for each distribution moment.
static const char *gkyl_distribution_moments_descriptions[] = {
  "Number density.",
  "Momentum density divided by the mass.",
  "Kinetic energy density times 2/mass.",
  "Parallel kinetic energy density times 2/mass.",
  "Perpendicular kinetic energy density times 2/mass.",
  "Second velocity moment (tensor) of the distribution.",
  "Heat flux density times 2/mass.",
  "Parallel heat flux density times 2/mass.",
  "Perpendicular heat flux density times 2/mass.",
  "Third velocity moment (tensor) times 2/mass.",
  "Number density, drift velocity, and thermal speed squared.",
  "Number density, drift velocity, parallel and perpendicular thermal speed squared.",
  "Number density, drift velocity, and thermal speed squared.",
  "Number density, momentum density divided by mass, and kinetic energy density times 2/mass.",
  "Number density, momentum density divided by mass, and parallel and perpendicular kinetic energy density times 2/mass.",
  "Number density, momentum density, and Hamiltonian energy density.",
  "Momentum density (divided by mass) from the velocity gradient of the Hamiltonian.",
  "Hamiltonian energy density (for neutrals).",
  "Number density, Hamiltonian energy density, and heat flux density times 2/mass.",
  "4-momentum (M0, M1).",
  "Stress-energy tensor.",
};

// Forward declare for use in function pointers
struct gkyl_mom_type;

/**
 * Function pointer type to compute the needed moment.
 */
typedef void (*momf_t)(const struct gkyl_mom_type *momt,
  const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param);

struct gkyl_mom_type {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of components in moment
  momf_t kernel; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flags;
  struct gkyl_mom_type *on_dev; // pointer to itself or device data
};

/**
 * Check if moment type is on device.
 *
 * @param momt Moment type to check
 * @return true if momt on device, false otherwise
 */
bool gkyl_mom_type_is_cu_dev(const struct gkyl_mom_type *momt);

/**
 * Acquire pointer to moment object. Delete using the release() method
 *
 * @param momt Moment object to get pointer from.
 * @return acquired object
 */
struct gkyl_mom_type* gkyl_mom_type_acquire(const struct gkyl_mom_type* momt);

/**
 * Delete moment object
 *
 * @param momt Moment object to delete.
 */
void gkyl_mom_type_release(const struct gkyl_mom_type* momt);

/**
 * Calculate moment specified by mom_type object.
 *
 * @param momt Moment type object
 * @param xc Cell center coordinates
 * @param dx Cell size in each direction
 * @param idx Index into phase-space cell
 * @param f Input pointer to distribution function in cell
 * @param out On output, contribution to moment from phase-space cell
 */
void gkyl_mom_type_calc(const struct gkyl_mom_type* momt,
  const double *xc, const double *dx, const int *idx,
  const double *f, double* GKYL_RESTRICT out, void *param);

/**
 * Get number of moments specified by mom_type object
 *
 * @param momt Moment type object
 * returns int Number of moments
 */
int gkyl_mom_type_num_mom(const struct gkyl_mom_type* momt);
