#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_kernels.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*sheath_reflectedf_t)(const double *vmap, const double q2Dm,
  const double *phi, const double *phiWall, const double *f, double *fRefl);

typedef struct { sheath_reflectedf_t kernels[3]; } sheath_reflectedf_kern_list;  // For use in kernel tables.
typedef struct { sheath_reflectedf_kern_list list[4]; } edged_sheath_reflectedf_kern_list;

// Serendipity  kernels.
GKYL_CU_D
static const edged_sheath_reflectedf_kern_list ser_sheath_reflect_list[] = {
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1, NULL },
          },
  },
  { .list={
           { bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1, NULL },
           { bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1, NULL },
          },
  },
};

struct gkyl_bc_sheath_gyrokinetic_kernels {
  sheath_reflectedf_t reflectedf;  // reflectedf kernel.
};

// Primary struct in this updater.
struct gkyl_bc_sheath_gyrokinetic {
  int dir; // Direction perpendicular to the sheath boundary.
  int cdim; // Conf-space dimensionality.
  enum gkyl_edge_loc edge; // Lower or upper boundary.
  const struct gkyl_basis *basis; // Phase-space basis.
  bool use_gpu; // Whether to run on GPU.
  double q2Dm; // charge-to-mass ratio times 2.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels;  // reflectedf kernel.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels_cu;  // device copy.
  const struct gkyl_range *skin_r, *ghost_r; // Skin and ghost ranges.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping.
};

void
gkyl_bc_gksheath_choose_reflectedf_kernel_cu(const struct gkyl_basis *basis, enum gkyl_edge_loc edge, struct gkyl_bc_sheath_gyrokinetic_kernels *kers);

GKYL_CU_D
static sheath_reflectedf_t
bc_gksheath_choose_reflectedf_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_sheath_reflect_list[edge].list[dim-2].kernels[poly_order-1];
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static void
bc_gksheath_reflect(int dir, const struct gkyl_basis *basis, int cdim, double *out, const double *inp)
{
  basis->flip_odd_sign(dir, inp, out);
  basis->flip_odd_sign(cdim, out, out); // cdim is the vpar direction.
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to apply the sheath BC.

 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param distf Distribution function array to apply BC to.
 * @param conf_r Configuration space range (to index phi).
 */
void gkyl_bc_sheath_gyrokinetic_advance_cu(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r);

#endif
