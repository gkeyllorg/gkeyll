#pragma once

// Private header for bc_sheath_gyrokinetic updater, not for direct use in user code.

#include "gkyl_bc_sheath_gyrokinetic.h"
#include "gkyl_bc_sheath_gyrokinetic_kernels.h"
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*sheath_reflectedf_t)(
  const double *vmap, const double *vcutsq, const double *f, double *fRefl
);

typedef struct {
  sheath_reflectedf_t kernels[3];
} sheath_reflectedf_kern_list; // For use in kernel tables.
typedef struct {
  sheath_reflectedf_kern_list list[4];
} edged_sheath_reflectedf_kern_list;

// Serendipity  kernels.
GKYL_CU_D static const edged_sheath_reflectedf_kern_list ser_sheath_reflect_list[] = {
  {
    .list =
      {{bc_sheath_gyrokinetic_reflectedf_lower_1x1v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_lower_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_lower_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_lower_3x2v_ser_p1, NULL}},
  },
  {
    .list =
      {{bc_sheath_gyrokinetic_reflectedf_upper_1x1v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_upper_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_upper_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_reflectedf_upper_3x2v_ser_p1, NULL}},
  }
};

typedef void (*sheath_vcutsq_const_t)(
  const double *phi, const double *phi_wall, double q2Dm, double *vcutsq_out
);
typedef struct {
  sheath_vcutsq_const_t kernels[3];
} sheath_vcutsq_const_kern_list; // For use in kernel tables.
typedef struct {
  sheath_vcutsq_const_kern_list dim_list[4];
} edged_sheath_vcutsq_const_kern_list;

// Serendipity kernels to calculate vcut using conducting sheath model.
GKYL_CU_D static const edged_sheath_vcutsq_const_kern_list ser_sheath_vcutsq_const_list[] = {
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_lower_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_lower_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_lower_3x2v_ser_p1, NULL}},
  },
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_upper_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_upper_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_const_upper_3x2v_ser_p1, NULL}},
  }
};

// Function pointer type for surrogate sheath BC to determine vcut factor from pre-computed NN output.
typedef void (*sheath_vcutsq_surr_t)(
  const double *vmap, const float *nn_out, int n_out, const double *temperature, const double *bmag,
  double *vcutsq_out
);

typedef struct {
  sheath_vcutsq_surr_t kernels[3];
} sheath_vcutsq_surr_kern_list; // For use in kernel tables.
typedef struct {
  sheath_vcutsq_surr_kern_list dim_list[4];
} edged_sheath_vcutsq_surr_kern_list;

// Serendipity to evaluate vcut from the KANN surrogate model output.
GKYL_CU_D static const edged_sheath_vcutsq_surr_kern_list ser_sheath_vcutsq_surr_list[] = {
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_lower_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_lower_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_lower_3x2v_ser_p1, NULL}},
  },
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_upper_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_upper_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_vcutsq_surr_upper_3x2v_ser_p1, NULL}},
  }
};

// Function pointer type for the build the KANN input.
typedef void (*sheath_surrogate_inp_t)(
  const double *phi, const double *phi_wall, const double *density, const double *temperature,
  const double *bmag, const double *bimpact_angle, int n_inp, float *inp_out
);
typedef struct {
  sheath_surrogate_inp_t kernels[3];
} sheath_inp_kern_list;
typedef struct {
  sheath_inp_kern_list dim_list[4];
} edged_sheath_inp_kern_list;

// Serendipity kernels to build the KANN input.
GKYL_CU_D static const edged_sheath_inp_kern_list ser_sheath_vcutsq_input_list[] = {
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_build_input_lower_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_build_input_lower_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_build_input_lower_3x2v_ser_p1, NULL}},
  },
  {
    .dim_list =
      {{NULL, NULL},
       {bc_sheath_gyrokinetic_build_input_upper_1x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_build_input_upper_2x2v_ser_p1, NULL},
       {bc_sheath_gyrokinetic_build_input_upper_3x2v_ser_p1, NULL}},
  }
};

struct gkyl_bc_sheath_gyrokinetic_kernels {
  sheath_reflectedf_t reflectedf; // reflectedf kernel.
  sheath_surrogate_inp_t vcutsq_input; // NN input construction kernel.
  sheath_vcutsq_surr_t vcutsq_surr; // vcut_fact DG projection kernel (uses pre-computed NN output).
  sheath_vcutsq_const_t vcutsq_const; // vcut calculation using conducting sheath model.
};

// Primary struct in this updater.
struct gkyl_bc_sheath_gyrokinetic {
  int dir; // Direction perpendicular to the sheath boundary.
  int cdim; // Conf-space dimensionality.
  enum gkyl_gyrokinetic_bc_type type; // Type of gyrokinetic BC.
  enum gkyl_edge_loc edge; // Lower or upper boundary.
  const struct gkyl_basis *basis; // Phase-space basis.
  bool use_gpu; // Whether to run on GPU.
  double q2Dm; // charge-to-mass ratio times 2.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels; // reflectedf kernel.
  struct gkyl_bc_sheath_gyrokinetic_kernels *kernels_cu; // device copy.
  const struct gkyl_range *skin_r, *ghost_r; // Skin and ghost ranges.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping.
  int vcutsq_dim; // Dimensionality of vcutsq array.
  struct gkyl_array *vcutsq; // factor for mu dependent vcut in sheath BC.
  struct gkyl_basis
    vcutsq_basis; // Basis for vcutsq array (expansion in perpendicular config space and mu).
  struct gkyl_range vcutsq_local; // Range for vcutsq array.
  struct gkyl_range perp_local; // Peprendicular conf space range.
  int perp_node_per_cell; // Number of nodes per cell in the perpendicular conf. space.
  void (*update_vcutsq)(
    const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
    const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
    const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
    const struct gkyl_range *conf_r
  ); // Function pointer to update vcutsq array.
  struct gkyl_kann_net *kann_net; // KANN sheath surrogate.
  struct gkyl_kn_vec *kann_inp; // Input vector for KANN surrogate (3 x number of perp nodes).
  struct gkyl_kn_vec *
    kann_out; // Output vector for KANN surrogate: dim_out x number of perp nodes, where each node block is dim_out/2 vcut values followed by dim_out/2 mu-grid values.
  struct gkyl_array *kann_infer_xy_out; // Output array for KANN on a perpendicular plane.
  double *nn_out; // Buffer to hold raw NN outputs for all perp nodes (device memory when using GPU).

  // Diagnostic output of the vcutsq array.
  bool vcutsq_diag_on; // Whether the vcutsq diagnostic output is set up.
  bool vcutsq_write; // Whether this rank abuts the boundary and writes vcutsq.
  struct gkyl_rect_grid vcutsq_grid; // Grid for vcutsq output (perp conf-space + mu).
  struct gkyl_range vcutsq_diag_range,
    vcutsq_diag_range_ext; // Output ranges (matching vcutsq_grid).
  struct gkyl_array *vcutsq_ho; // Host buffer for vcutsq output.
};

void gkyl_bc_gksheath_choose_reflectedf_kernel_cu(
  const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers
);

GKYL_CU_D static sheath_reflectedf_t
bc_gksheath_choose_reflectedf_kernel(const struct gkyl_basis *basis, enum gkyl_edge_loc edge)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  sheath_reflectedf_t kern = NULL;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
      kern = ser_sheath_reflect_list[edge].list[dim - 2].kernels[poly_order - 1];
      break;
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kern = ser_sheath_reflect_list[edge].list[dim - 2].kernels[poly_order - 1];
      break;
    default:
      kern = NULL;
      break;
  }
  assert(kern
  ); // Reflection kernel must be non-null since we always use reflection BCs in the sheath updater.
  return kern;
}

GKYL_CU_D static void
bc_gksheath_reflect(
  int dir, const struct gkyl_basis *basis, int cdim, double *out, const double *inp
)
{
  basis->flip_odd_sign(dir, inp, out);
  basis->flip_odd_sign(cdim, out, out); // cdim is the vpar direction.
}

void gkyl_bc_gksheath_choose_surrogate_kernels_cu(
  const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers
);

GKYL_CU_D static void
bc_gksheath_choose_vcutsq_surr_kernels(
  const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers
)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->vcutsq_surr =
        ser_sheath_vcutsq_surr_list[edge].dim_list[dim - 2].kernels[poly_order - 1];
      kers->vcutsq_input =
        ser_sheath_vcutsq_input_list[edge].dim_list[dim - 2].kernels[poly_order - 1];
      break;
    default:
      assert(false);
  }
}

void gkyl_bc_gksheath_choose_const_kernels_cu(
  const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers
);

GKYL_CU_D static void
bc_gksheath_choose_vcutsq_const_kernels(
  const struct gkyl_basis *basis, enum gkyl_edge_loc edge,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers
)
{
  int dim = basis->ndim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->vcutsq_const =
        ser_sheath_vcutsq_const_list[edge].dim_list[dim - 2].kernels[poly_order - 1];
      break;
    default:
      assert(false);
  }
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
void gkyl_bc_sheath_gyrokinetic_advance_cu(
  const struct gkyl_bc_sheath_gyrokinetic *up, struct gkyl_array *distf,
  const struct gkyl_range *conf_r
);

/**
 * CUDA device function to update the vcutsq array with surrogate.
 * 
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param dens Density array.
 * @param temp Temperature array.
 * @param bmag Magnetic field magnitude array.
 * @param bimpact_angle Magnetic field impact angle array.
 * @param conf_r Configuration space range (to index phi, dens, temp, bmag, and bimpact_angle).
 * 
 */
void bc_gksheath_update_vcutsq_surrogate_cu(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
  const struct gkyl_range *conf_r
);

/**
 * CUDA device function to update the vcutsq array with conducting sheath model.
 * 
 * @param up BC updater.
 * @param phi Electrostatic potential.
 * @param phi_wall Wall potential.
 * @param dens Density array.
 * @param temp Temperature array.
 * @param bmag Magnetic field magnitude array.
 * @param bimpact_angle Magnetic field impact angle array.
 * @param conf_r Configuration space range (to index phi, dens, temp, bmag, and bimpact_angle).
 * 
 */
void bc_gksheath_update_vcutsq_const_cu(
  const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, const struct gkyl_array *dens, const struct gkyl_array *temp,
  const struct gkyl_array *bmag, const struct gkyl_array *bimpact_angle,
  const struct gkyl_range *conf_r
);
#endif
