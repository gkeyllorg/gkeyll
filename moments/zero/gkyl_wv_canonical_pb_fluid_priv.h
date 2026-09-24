#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_can_pb_incompress_euler {
  struct gkyl_wv_eqn eqn; // base object
};

/**
 * Free incompressible Euler eqn object.
 *
 * @param ref Reference counter for incompressible Euler eqn
 */
void gkyl_wv_can_pb_incompress_euler_free(const struct gkyl_ref_count *ref);

struct wv_can_pb_hasegawa_mima {
  struct gkyl_wv_eqn eqn; // base object
};

/**
 * Free Hasegawa-Mima eqn object.
 *
 * @param ref Reference counter for Hasegawa-Mima eqn
 */
void gkyl_wv_can_pb_hasegawa_mima_free(const struct gkyl_ref_count *ref);

struct wv_can_pb_hasegawa_wakatani {
  struct gkyl_wv_eqn eqn; // base object
  double alpha; // Adiabaticity parameter for adiabatic coupling of vorticity and density.
  bool is_modified; // is_modified Boolean parameter for if we are doing the modified Hasegawa-Wakatani
};

/**
 * Free Hasegawa-Wakatani eqn object.
 *
 * @param ref Reference counter for Hasegawa-Wakatani eqn
 */
void gkyl_wv_can_pb_hasegawa_wakatani_free(const struct gkyl_ref_count *ref);
