// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_pkpm_layout.h>

typedef void (*pkpm_coll_nu_advance_t)(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *prim_self,
  const struct gkyl_array *prim_other, struct gkyl_array *nu
);

struct gkyl_pkpm_coll_nu {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double vtsq_min; // floor on each thermal speed squared
  pkpm_coll_nu_advance_t advance; // the reference frequency or its Spitzer scaling
};
