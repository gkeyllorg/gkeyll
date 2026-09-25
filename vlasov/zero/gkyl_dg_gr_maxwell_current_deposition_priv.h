// Private header: not for direct use.
#pragma once

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_util.h>

typedef void (*gr_maxwell_current_deposition_t)(
  const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
  const double *vierb_con_nodal, const double *m0, const double *m1i, double *GKYL_RESTRICT out
);

// The cv_index[cd].vdim[vd] is used to index the kernel tables below.
GKYL_CU_TABLE
static const struct {
  int vdim[4];
} current_dep_cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense.
  {-1, 0, 1, 2}, // 1x kernel indices.
  {-1, -1, 3, 4}, // 2x kernel indices.
  {-1, -1, -1, 5} // 3x kernel indices.
};

typedef struct {
  gr_maxwell_current_deposition_t kernels[4];
} gkyl_gr_maxwell_current_deposition_kern_list;

struct gkyl_dg_gr_maxwell_current_deposition {
  int cdim; // Configuration-space dimensions.
  int vdim; // Velocity-space dimensions.
  int num_basis; // Number of configuration-space basis functions.
  gr_maxwell_current_deposition_t current_deposition; // Current deposition kernel.

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_current_deposition *on_dev; // Pointer to itself or device data.
};

GKYL_CU_D static const gkyl_gr_maxwell_current_deposition_kern_list
  ser_current_deposition_kernels[] = {
    {NULL, gr_maxwell_current_dep_1x_1v_ser_p1, gr_maxwell_current_dep_1x_1v_ser_p2, NULL}, // 0
    {NULL, gr_maxwell_current_dep_1x_2v_ser_p1, gr_maxwell_current_dep_1x_2v_ser_p2, NULL}, // 1
    {NULL, gr_maxwell_current_dep_1x_3v_ser_p1, gr_maxwell_current_dep_1x_3v_ser_p2, NULL}, // 2
    {NULL, gr_maxwell_current_dep_2x_2v_ser_p1, gr_maxwell_current_dep_2x_2v_ser_p2, NULL}, // 3
    {NULL, gr_maxwell_current_dep_2x_3v_ser_p1, gr_maxwell_current_dep_2x_3v_ser_p2, NULL}, // 4
    {NULL, gr_maxwell_current_dep_3x_3v_ser_p1, gr_maxwell_current_dep_3x_3v_ser_p2, NULL} // 5
};

GKYL_CU_D static const gkyl_gr_maxwell_current_deposition_kern_list
  tensor_current_deposition_kernels[] = {
    {NULL, gr_maxwell_current_dep_1x_1v_tensor_p1, gr_maxwell_current_dep_1x_1v_tensor_p2, NULL
    }, // 0
    {NULL, gr_maxwell_current_dep_1x_2v_tensor_p1, gr_maxwell_current_dep_1x_2v_tensor_p2, NULL
    }, // 1
    {NULL, gr_maxwell_current_dep_1x_3v_tensor_p1, gr_maxwell_current_dep_1x_3v_tensor_p2, NULL
    }, // 2
    {NULL, gr_maxwell_current_dep_2x_2v_tensor_p1, gr_maxwell_current_dep_2x_2v_tensor_p2, NULL
    }, // 3
    {NULL, gr_maxwell_current_dep_2x_3v_tensor_p1, gr_maxwell_current_dep_2x_3v_tensor_p2, NULL
    }, // 4
    {NULL, gr_maxwell_current_dep_3x_3v_tensor_p1, gr_maxwell_current_dep_3x_3v_tensor_p2, NULL
    } // 5
};
