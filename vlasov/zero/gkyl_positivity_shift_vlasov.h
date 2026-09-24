#pragma once

#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_positivity_shift_vlasov gkyl_positivity_shift_vlasov;

/**
 * Create a new upater which shifts f up to a floor wherever f<0, and computes
 * contribution to the integrated moments of the shift in f in each cell (the
 * final reduction needs to be performed outside of this updater, as it may
 * involve a communicator).
 *
 * @param cbasis Conf-space basis.
 * @param pbasis Phase-space basis.
 * @param pgrid Phase-space grid.
 * @param conf_rng_ext Extended configuration space range.
 * @param use_gpu bool to determine if on GPU.
 * @return New positivity shift updater pointer.
 */
struct gkyl_positivity_shift_vlasov*
gkyl_positivity_shift_vlasov_new(struct gkyl_basis cbasis, struct gkyl_basis pbasis,
  struct gkyl_rect_grid pgrid, const struct gkyl_range *conf_rng_ext, bool use_gpu);

/**
 * Run the positivity shift updater in the indicated range.
 *
 * @param up Positivity shift updater.
 * @param conf_rng Config-space range.
 * @param phase_rng Phase-space range.
 * @param distf Distribution function array.
 * @param m0 Output M0 moment array.
 * @param delta_m0 M0 moment of the shift in f.
 */
void
gkyl_positivity_shift_vlasov_advance(gkyl_positivity_shift_vlasov* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT m0,
  struct gkyl_array *GKYL_RESTRICT delta_m0);

/**
 * Release the memory associated with this positivity shift updater.
 *
 * @param up Positivity shift updater.
 */
void
gkyl_positivity_shift_vlasov_release(gkyl_positivity_shift_vlasov* up);
