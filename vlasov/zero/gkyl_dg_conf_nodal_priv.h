// Private header: not for direct use
#pragma once

#include <gkyl_dg_conf_nodal_kernels.h>
#include <gkyl_util.h>

typedef void (*dg_conf_nodal_kern_t)(const double *in, double *GKYL_RESTRICT out);

// for use in kernel tables
typedef struct {
  dg_conf_nodal_kern_t kernels[3];
} gkyl_dg_conf_nodal_kern_list;

// Configuration-space kernels (p=1 serendipity), by cdim.
GKYL_CU_D static const gkyl_dg_conf_nodal_kern_list dg_conf_nodal_conf_to_nodal_kernels = {
  dg_conf_nodal_conf_to_nodal_1x_ser_p1, dg_conf_nodal_conf_to_nodal_2x_ser_p1,
  dg_conf_nodal_conf_to_nodal_3x_ser_p1
};

GKYL_CU_D static const gkyl_dg_conf_nodal_kern_list dg_conf_nodal_conf_to_modal_kernels = {
  dg_conf_nodal_conf_to_modal_1x_ser_p1, dg_conf_nodal_conf_to_modal_2x_ser_p1,
  dg_conf_nodal_conf_to_modal_3x_ser_p1
};

// Phase-space kernels (hybrid basis), by cdim and then by vdim.
GKYL_CU_D static const gkyl_dg_conf_nodal_kern_list dg_conf_nodal_phase_to_nodal_kernels[] = {
  // 1x kernels
  {dg_conf_nodal_phase_to_nodal_1x1v_hyb_p1, NULL, NULL}, // 0
  // 2x kernels
  {dg_conf_nodal_phase_to_nodal_2x1v_hyb_p1, NULL, NULL}, // 1
  // 3x kernels
  {dg_conf_nodal_phase_to_nodal_3x1v_hyb_p1, NULL, NULL} // 2
};

GKYL_CU_D static const gkyl_dg_conf_nodal_kern_list dg_conf_nodal_phase_to_modal_kernels[] = {
  // 1x kernels
  {dg_conf_nodal_phase_to_modal_accumulate_1x1v_hyb_p1, NULL, NULL}, // 0
  // 2x kernels
  {dg_conf_nodal_phase_to_modal_accumulate_2x1v_hyb_p1, NULL, NULL}, // 1
  // 3x kernels
  {dg_conf_nodal_phase_to_modal_accumulate_3x1v_hyb_p1, NULL, NULL} // 2
};

struct gkyl_dg_conf_nodal {
  int num_conf_basis; // number of configuration-space basis functions (and of nodes)
  int num_phase_basis; // number of phase-space basis functions (0 without a phase basis)
  dg_conf_nodal_kern_t conf_to_nodal; // configuration-space expansion to node values
  dg_conf_nodal_kern_t conf_to_modal; // configuration-space node values to expansion
  dg_conf_nodal_kern_t phase_to_nodal; // phase-space expansion to the nodal form
  dg_conf_nodal_kern_t phase_to_modal; // accumulate the nodal form onto a phase-space expansion
};
