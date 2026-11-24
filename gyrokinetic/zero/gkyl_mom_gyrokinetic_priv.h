#pragma once

// Private header for mom_gyrokinetic updater, not for direct use in user code.

#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D static struct {
  int vdim[3];
} cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1, 0, 1}, // 1x kernel indices
  {-1, -1, 2}, // 2x kernel indices
  {-1, -1, 3} // 3x kernel indices
};

// Kernel signature.
typedef void (*wmomf_t)(
  const double *dxv, const double *vmap, double m_, double q_, const double *bmag,
  const double *phi, const double *f, double *GKYL_RESTRICT out
);

// for use in kernel tables
typedef struct {
  wmomf_t kernels[3]; // One for each poly_order.
} gkyl_gyrokinetic_mom_kern_list;

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M0_1x1v_ser_p1, mom_gyrokinetic_M0_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M0_1x2v_ser_p1, mom_gyrokinetic_M0_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M0_2x2v_ser_p1, mom_gyrokinetic_M0_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M0_3x2v_ser_p1, NULL} // 3
};

// M1 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m1_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M1_1x1v_ser_p1, mom_gyrokinetic_M1_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M1_1x2v_ser_p1, mom_gyrokinetic_M1_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M1_2x2v_ser_p1, mom_gyrokinetic_M1_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M1_3x2v_ser_p1, NULL} // 3
};

// M2 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M2_1x1v_ser_p1, mom_gyrokinetic_M2_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M2_1x2v_ser_p1, mom_gyrokinetic_M2_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M2_2x2v_ser_p1, mom_gyrokinetic_M2_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M2_3x2v_ser_p1, NULL} // 3
};

// M3 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m3_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M3_1x1v_ser_p1, mom_gyrokinetic_M3_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M3_1x2v_ser_p1, mom_gyrokinetic_M3_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M3_2x2v_ser_p1, mom_gyrokinetic_M3_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M3_3x2v_ser_p1, NULL} // 3
};

// M2 parallel kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m2par_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M2PAR_1x1v_ser_p1, mom_gyrokinetic_M2PAR_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M2PAR_1x2v_ser_p1, mom_gyrokinetic_M2PAR_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M2PAR_2x2v_ser_p1, mom_gyrokinetic_M2PAR_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M2PAR_3x2v_ser_p1, NULL} // 3
};

// M2 perpendicular kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m2perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, mom_gyrokinetic_M2PERP_1x2v_ser_p1, mom_gyrokinetic_M2PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M2PERP_2x2v_ser_p1, mom_gyrokinetic_M2PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M2PERP_3x2v_ser_p1, NULL} // 3
};

// M3 parallel kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m3par_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M3PAR_1x1v_ser_p1, mom_gyrokinetic_M3PAR_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M3PAR_1x2v_ser_p1, mom_gyrokinetic_M3PAR_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M3PAR_2x2v_ser_p1, mom_gyrokinetic_M3PAR_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M3PAR_3x2v_ser_p1, NULL} // 3
};

// M3 perpendicular kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m3perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, mom_gyrokinetic_M3PERP_1x2v_ser_p1, mom_gyrokinetic_M3PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M3PERP_2x2v_ser_p1, mom_gyrokinetic_M3PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M3PERP_3x2v_ser_p1, NULL} // 3
};

// Density, parallel momentum and kinetic energy computed together.
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M0M1M2_1x1v_ser_p1, mom_gyrokinetic_M0M1M2_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M0M1M2_1x2v_ser_p1, mom_gyrokinetic_M0M1M2_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M0M1M2_2x2v_ser_p1, mom_gyrokinetic_M0M1M2_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M0M1M2_3x2v_ser_p1, NULL} // 3
};

// Density, parallel momentum, kinetic energy and energy flux computed together.
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2m3_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1, mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1, mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1, mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1, NULL} // 3
};

// Density, parallel momentum, and parallel and perpendicular kinetic energy computed together.
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2parm2perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1, mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2
  }, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1, mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2
  }, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1, NULL} // 3
};

// Density, parallel momentum, and parallel and perpendicular kinetic energy,
// and parallel and perpendicular energy flux computed together.
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2parm2perpm3parm3perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1,
   mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1,
   mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1, NULL} // 3
};

// Hamiltonian moment kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_hamiltonian_kernels[] = {
  // 1x kernels
  {NULL, mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1, mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2}, // 0
  {NULL, mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1, mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1, mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1, NULL} // 3
};

// Integrated M0 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m0_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M0_1x1v_ser_p1, int_mom_gyrokinetic_M0_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M0_1x2v_ser_p1, int_mom_gyrokinetic_M0_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M0_2x2v_ser_p1, int_mom_gyrokinetic_M0_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M0_3x2v_ser_p1, NULL} // 3
};

// Integrated M1 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m1_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M1_1x1v_ser_p1, int_mom_gyrokinetic_M1_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M1_1x2v_ser_p1, int_mom_gyrokinetic_M1_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M1_2x2v_ser_p1, int_mom_gyrokinetic_M1_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M1_3x2v_ser_p1, NULL} // 3
};

// Integrated M2 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m2_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M2_1x1v_ser_p1, int_mom_gyrokinetic_M2_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M2_1x2v_ser_p1, int_mom_gyrokinetic_M2_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M2_2x2v_ser_p1, int_mom_gyrokinetic_M2_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M2_3x2v_ser_p1, NULL} // 3
};

// Integrated M3 kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m3_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M3_1x1v_ser_p1, int_mom_gyrokinetic_M3_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M3_1x2v_ser_p1, int_mom_gyrokinetic_M3_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M3_2x2v_ser_p1, int_mom_gyrokinetic_M3_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M3_3x2v_ser_p1, NULL} // 3
};

// Integrated M2 parallel kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m2par_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M2PAR_1x1v_ser_p1, int_mom_gyrokinetic_M2PAR_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M2PAR_1x2v_ser_p1, int_mom_gyrokinetic_M2PAR_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M2PAR_2x2v_ser_p1, int_mom_gyrokinetic_M2PAR_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M2PAR_3x2v_ser_p1, NULL} // 3
};

// Integrated M2 perpendicular kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m2perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, int_mom_gyrokinetic_M2PERP_1x2v_ser_p1, int_mom_gyrokinetic_M2PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M2PERP_2x2v_ser_p1, int_mom_gyrokinetic_M2PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M2PERP_3x2v_ser_p1, NULL} // 3
};

// Integrated M3 parallel kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m3par_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M3PAR_1x1v_ser_p1, int_mom_gyrokinetic_M3PAR_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M3PAR_1x2v_ser_p1, int_mom_gyrokinetic_M3PAR_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M3PAR_2x2v_ser_p1, int_mom_gyrokinetic_M3PAR_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M3PAR_3x2v_ser_p1, NULL} // 3
};

// Integrated M3 perpendicular kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m3perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, int_mom_gyrokinetic_M3PERP_1x2v_ser_p1, int_mom_gyrokinetic_M3PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M3PERP_2x2v_ser_p1, int_mom_gyrokinetic_M3PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M3PERP_3x2v_ser_p1, NULL} // 3
};

// Integrated M0M1M2_moments kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2_1x1v_ser_p1, int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M0M1M2_1x2v_ser_p1, int_mom_gyrokinetic_M0M1M2_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2_2x2v_ser_p1, int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1, NULL} // 3
};

// Integrated M0M1M2M3_moments kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2m3_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1, int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2}, // 0
  {NULL, int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1, int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1, int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1, NULL} // 3
};

// Integrated M0M1M2PARM2PERP_moments kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2parm2perp_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL}, // 0
  {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1,
   int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2}, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1,
   int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2}, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1, NULL} // 3
};

// Integrated M0M1M2PARM2PERPM3PARM3PERP_moments kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2parm2perpm3parm3perp_kernels[] =
  {
    // 1x kernels
    {NULL, NULL, NULL}, // 0
    {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1,
     int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2}, // 1
    // 2x kernels
    {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1,
     int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2}, // 2
    // 3x kernels
    {NULL, int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1, NULL} // 3
};

// Integrated Hamiltonian kernel list
GKYL_CU_D static const gkyl_gyrokinetic_mom_kern_list ser_int_hamiltonian_kernels[] = {
  // 1x kernels
  {NULL, int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1, int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2
  }, // 0
  {NULL, int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1, int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2
  }, // 1
  // 2x kernels
  {NULL, int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1, int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2
  }, // 2
  // 3x kernels
  {NULL, int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1, NULL} // 3
};

struct gkyl_mom_gyrokinetic_kernels {
  wmomf_t wmom; // Kernels that computes moment.
};

// Object definition.
struct gkyl_mom_gyrokinetic {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  double mass; // Species mass.
  double charge; // Species charge.
  const struct gk_geometry *gk_geom; // Pointer to geometry object.
  const struct gkyl_velocity_map *vel_map; // Pointer to velocity mapping object.
  bool use_gpu; // Whether to run on the GPU.
  struct gkyl_mom_gyrokinetic_kernels *kernels;
  int num_mom; // Number of moments to compute.
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.

void mom_gk_choose_kernel_cu(
  struct gkyl_mom_gyrokinetic_kernels *kernels, int *num_mom, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type, bool is_integrated
);

void gkyl_mom_gyrokinetic_advance_cu(
  struct gkyl_mom_gyrokinetic *up, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT phi,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout
);
#endif

GKYL_CU_D static void
mom_gk_choose_kernel(
  struct gkyl_mom_gyrokinetic_kernels *kernels, int *num_mom, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, enum gkyl_distribution_moments mom_type, bool is_integrated
)
{
  enum gkyl_basis_type cbasis_type = cbasis->b_type;
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = pbasis->poly_order;

  // Choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, *m3_kernels,
    *m2par_kernels, *m2perp_kernels, *m3par_kernels, *m3perp_kernels, *m0m1m2_kernels,
    *m0m1m2m3_kernels, *m0m1m2parm2perp_kernels, *m0m1m2parm2perpm3parm3perp_kernels,
    *hamiltonian_kernels;

  switch (cbasis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (is_integrated) {
        m0_kernels = ser_int_m0_kernels;
        m1_kernels = ser_int_m1_kernels;
        m2_kernels = ser_int_m2_kernels;
        m3_kernels = ser_int_m3_kernels;
        m2par_kernels = ser_int_m2par_kernels;
        m2perp_kernels = ser_int_m2perp_kernels;
        m3par_kernels = ser_int_m3par_kernels;
        m3perp_kernels = ser_int_m3perp_kernels;
        m0m1m2_kernels = ser_int_m0m1m2_kernels;
        m0m1m2m3_kernels = ser_int_m0m1m2m3_kernels;
        m0m1m2parm2perp_kernels = ser_int_m0m1m2parm2perp_kernels;
        m0m1m2parm2perpm3parm3perp_kernels = ser_int_m0m1m2parm2perpm3parm3perp_kernels;
        hamiltonian_kernels = ser_int_hamiltonian_kernels;
      } else {
        m0_kernels = ser_m0_kernels;
        m1_kernels = ser_m1_kernels;
        m2_kernels = ser_m2_kernels;
        m3_kernels = ser_m3_kernels;
        m2par_kernels = ser_m2par_kernels;
        m2perp_kernels = ser_m2perp_kernels;
        m3par_kernels = ser_m3par_kernels;
        m3perp_kernels = ser_m3perp_kernels;
        m0m1m2_kernels = ser_m0m1m2_kernels;
        m0m1m2m3_kernels = ser_m0m1m2m3_kernels;
        m0m1m2parm2perp_kernels = ser_m0m1m2parm2perp_kernels;
        m0m1m2parm2perpm3parm3perp_kernels = ser_m0m1m2parm2perpm3parm3perp_kernels;
        hamiltonian_kernels = ser_hamiltonian_kernels;
      }
      break;

    default:
      assert(false);
      break;
  }

  if (mom_type == GKYL_F_MOMENT_M0) { // Density.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M1) { // Parallel momentum.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M2) { // Total kinetic energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M3) { // Total energy flux.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M2PAR) { // Parallel energy.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m2par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M2PERP) { // Perpendicular energy.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M3PAR) { // Parallel energy flux.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m3par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M3PERP) { // Perpendicular energy flux.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 1;
  } else if (mom_type == GKYL_F_MOMENT_M0M1M2) {
    // Density, parallel momentum, and total energy computed together.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m0m1m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 3;
  } else if (mom_type == GKYL_F_MOMENT_M0M1M2M3) {
    // Density, parallel momentum, total energy and total energy flux computed together.
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m0m1m2m3_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 4;
  } else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP) {
    // Density, parallel momentum, parallel and perpendicular
    // kinetic energy computed together.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = m0m1m2parm2perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = vdim + 2;
  } else if (mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERPM3PARM3PERP) {
    // Density, parallel momentum, parallel and perpendicular
    // kinetic energy, and parallel and perpendicular energy flux computed together.
    assert(vdim == 2);
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(
      NULL != m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]
    );

    kernels->wmom =
      m0m1m2parm2perpm3parm3perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 2 * vdim + 2;
  } else if (mom_type == GKYL_F_MOMENT_HAMILTONIAN) {
    // M0, mass*M0 and total particle energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

    kernels->wmom = hamiltonian_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    *num_mom = 3;
  } else {
    // string not recognized
    printf("Error: requested moment %d.\n", mom_type);
    assert(false);
  }
}
