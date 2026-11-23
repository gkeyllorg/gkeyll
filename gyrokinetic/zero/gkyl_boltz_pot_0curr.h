#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

typedef struct gkyl_boltz_pot_0curr gkyl_boltz_pot_0curr;

/**
 * Create an isothermal Boltzmann electron potential updater imposing
 * zero net outward current summed over both ends of each field line.
 * The parallel direction is the last configuration coordinate.
 *
 * @param basis Configuration-space serendipity basis (1-3 dimensions, p=1 or p=2).
 * @param mass_e Electron mass.
 * @param charge_e Signed electron charge (negative).
 * @param temp_e Electron temperature in energy units.
 * @param use_gpu Whether to operate on GPU arrays.
 * @return Newly allocated updater; release with gkyl_boltz_pot_0curr_release.
 */
gkyl_boltz_pot_0curr *gkyl_boltz_pot_0curr_new(
  const struct gkyl_basis *basis, double mass_e, double charge_e, double temp_e, bool use_gpu
);

/**
 * Assemble values at one physical end of a field line. Evaluate density and
 * prescribed wall potential on the skin face; copy the supplied nodal flux and
 * geometry factor. This does not compute the sheath plasma potential.
 * Output entries contain four blocks of nq surface Gauss-Legendre nodal values:
 * (F, K, m0_elc, phi_wall), with nq=(p+1)^(cdim-1).
 * F=sum_i (q_i/|q_e|)*Gamma_i is outward ion flux per transverse coordinate area.
 * K is the physical face-area Jacobian times |b dot n|.
 * All arrays must use the updater's CPU/GPU backend.
 *
 * @param up Updater.
 * @param edge GKYL_LOWER_EDGE or GKYL_UPPER_EDGE in the parallel direction.
 * @param skin Configuration-space skin range retaining density/wall array strides.
 * @param flux F at surface nodes; skin.volume entries of nq components.
 * @param area K at the same nodes, with the same layout as flux.
 * @param m0_elc Physical quasineutral electron density, without J, in the volume modal basis.
 * @param wall Prescribed wall potential in the volume modal basis; NULL means zero.
 * @param end_vals Output boundary values; skin.volume entries of 4*nq components.
 */
void gkyl_boltz_pot_0curr_boundary(
  const gkyl_boltz_pot_0curr *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
);

/**
 * Calculate the reference density and potential for zero net outward current summed over both ends of each field line.
 * Boundary values must already be available on the calling rank. The output
 * has two blocks of basis.num_basis modal coefficients, (n_ref, phi_ref),
 * constant along the parallel coordinate. The density clamp is applied later
 * by _phi_calc. All arrays must use the updater's CPU/GPU backend.
 *
 * @param up Updater.
 * @param ghost Output ghost range retaining nphi_ref strides; _phi_calc reads
 *   the local lower parallel ghost slab.
 * @param lower Lower-end values from _boundary, ordered by transverse cell.
 * @param upper Upper-end values, with the same transverse order as lower.
 * @param nphi_ref Output reference density and potential on the ghost range.
 */
void gkyl_boltz_pot_0curr_reference(
  const gkyl_boltz_pot_0curr *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
);

/**
 * Compute the DG potential at volume quadrature nodes and project onto the basis:
 * phi = phi_ref + (T_e/|q_e|)*log(max(m0_elc,n_ref)/n_ref).
 * A nonfinite or nonpositive reference density produces zero nodal potential.
 * Reference data are read from the local lower parallel ghost slab. This does
 * not perform the app's subsequent FEM projection. All arrays must use the
 * updater's CPU/GPU backend.
 *
 * @param up Updater.
 * @param local Interior configuration range retaining density/output strides.
 * @param local_ext Extended configuration range including the lower ghost slab.
 * @param m0_elc Physical electron density without J, in the volume modal basis.
 * @param nphi_ref Two modal blocks (n_ref, phi_ref) in the lower ghost slab.
 * @param phi Output volume-modal potential on local.
 */
void gkyl_boltz_pot_0curr_phi_calc(
  const gkyl_boltz_pot_0curr *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
);

/**
 * Release the updater and its device resources.
 *
 * @param up Updater to release.
 */
void gkyl_boltz_pot_0curr_release(gkyl_boltz_pot_0curr *up);
