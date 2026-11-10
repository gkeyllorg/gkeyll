// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_pkpm_coll_nodal.h>
#include <gkyl_range.h>

typedef void (*pkpm_coll_nu_conf_to_nodal_t)(const double *c, double *GKYL_RESTRICT n);

GKYL_CU_D static const pkpm_coll_nu_conf_to_nodal_t pkpm_coll_nu_conf_to_nodal_list[] = {
  pkpm_coll_conf_to_nodal_1x, pkpm_coll_conf_to_nodal_2x, pkpm_coll_conf_to_nodal_3x
};

typedef void (*pkpm_coll_nu_advance_t)(
  const struct gkyl_pkpm_coll_nu *up, const struct gkyl_range *range,
  const struct gkyl_array *nu_ref, const struct gkyl_array *n_other,
  const struct gkyl_array *vtsq_self, const struct gkyl_array *vtsq_other, struct gkyl_array *nu
);

struct gkyl_pkpm_coll_nu {
  int num_nodes; // number of nodes (and of p=1 basis functions) per cell
  double vtsq_min; // floor on each thermal speed squared
  pkpm_coll_nu_conf_to_nodal_t conf_to_nodal; // modal-to-nodal conversion of this cdim
  pkpm_coll_nu_advance_t advance; // the reference frequency or its Spitzer scaling
};
