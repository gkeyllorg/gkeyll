#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_coll_nu gkyl_pkpm_coll_nu;

/**
 * Create an updater for the collision frequency at the 2^cdim interior Gauss-Legendre nodes of
 * configuration space. Nodal arrays have the same number of components as the p=1 modal ones.
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param spitzer True to scale the reference frequency by n / (vtsq_self + vtsq_other)^(3/2),
 *   false to use the reference frequency as it is.
 * @param vtsq_min Floor on each thermal speed squared (used only if spitzer is true).
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_pkpm_coll_nu *gkyl_pkpm_coll_nu_new(
  const struct gkyl_basis *conf_basis, bool spitzer, double vtsq_min, bool use_gpu
);

/**
 * Compute the collision frequency at the nodes. All arrays hold node values. The result is zero
 * at nodes where n_other is not positive. n_other, vtsq_self and vtsq_other are not read unless
 * the updater was created with spitzer = true.
 *
 * @param up Updater.
 * @param range Configuration-space range.
 * @param nu_ref Reference collision frequency.
 * @param n_other Number density of the colliding partner.
 * @param vtsq_self Thermal speed squared of this species.
 * @param vtsq_other Thermal speed squared of the colliding partner.
 * @param nu Output collision frequency.
 */
void gkyl_pkpm_coll_nu_advance(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *n_other,
  const struct gkyl_array *vtsq_self, const struct gkyl_array *vtsq_other, struct gkyl_array *nu
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_coll_nu_release(struct gkyl_pkpm_coll_nu *up);
