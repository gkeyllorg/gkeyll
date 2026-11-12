#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_pkpm_coll_nu gkyl_pkpm_coll_nu;

/**
 * Create an updater for the collision frequency of a species s colliding with a species r
 * (r = s for self collisions) at the 2^cdim interior Gauss-Legendre nodes of configuration
 * space. Nodal arrays have the same number of components as the p=1 modal ones.
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param spitzer True to scale the reference frequency by n_r / (vtsq_s + vtsq_r)^(3/2),
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
 * at nodes where the density of the partner is not positive. The node fields are not read unless
 * the updater was created with spitzer = true.
 *
 * @param up Updater.
 * @param range Configuration-space range.
 * @param nu_ref Reference collision frequency.
 * @param prim_self Node fields [n, u (3), vtsq] of this species.
 * @param prim_other Node fields [n, u (3), vtsq] of the colliding partner.
 * @param nu Output collision frequency.
 */
void gkyl_pkpm_coll_nu_advance(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *prim_self,
  const struct gkyl_array *prim_other, struct gkyl_array *nu
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_pkpm_coll_nu_release(struct gkyl_pkpm_coll_nu *up);
