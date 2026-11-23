#pragma once
#include <assert.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_boltz_potential_kernels.h>

// Both reference kernels accept two ends; the ambipolar kernel ignores upper.
typedef void (*boltz_reference_t)(double, double, double, const double *, const double *, double *);
typedef void (*boltz_boundary_t)(const double *, const double *, const double *, const double *, double *);
typedef void (*boltz_phi_t)(double, double, const double *, const double *, double *);

struct boltz_potential {
  int cdim; // Configuration dimension; last direction is parallel.
  double mass_e, charge_e, temp_e; // Electron mass, signed charge, and temperature.
  bool use_gpu; // Whether arrays reside on the GPU.
  boltz_boundary_t boundary[2]; // Lower/upper face packing kernels.
  boltz_reference_t reference; // Model-specific reference kernel.
  boltz_phi_t phi; // Shared volume response kernel.
  struct boltz_potential *on_dev; // Device copy, or self on CPU.
};

GKYL_CU_TABLE static const boltz_boundary_t boltz_boundary_lower_kernels[3][2] = {
  {boltz_potential_boundary_lower_1x_ser_p1, boltz_potential_boundary_lower_1x_ser_p2},
  {boltz_potential_boundary_lower_2x_ser_p1, boltz_potential_boundary_lower_2x_ser_p2},
  {boltz_potential_boundary_lower_3x_ser_p1, boltz_potential_boundary_lower_3x_ser_p2}
};

GKYL_CU_TABLE static const boltz_boundary_t boltz_boundary_upper_kernels[3][2] = {
  {boltz_potential_boundary_upper_1x_ser_p1, boltz_potential_boundary_upper_1x_ser_p2},
  {boltz_potential_boundary_upper_2x_ser_p1, boltz_potential_boundary_upper_2x_ser_p2},
  {boltz_potential_boundary_upper_3x_ser_p1, boltz_potential_boundary_upper_3x_ser_p2}
};

GKYL_CU_TABLE static const boltz_phi_t boltz_phi_kernels[3][2] = {
  {boltz_potential_phi_1x_ser_p1, boltz_potential_phi_1x_ser_p2},
  {boltz_potential_phi_2x_ser_p1, boltz_potential_phi_2x_ser_p2},
  {boltz_potential_phi_3x_ser_p1, boltz_potential_phi_3x_ser_p2}
};

// Select common and model-specific kernels once during construction.
GKYL_CU_DH static inline void
boltz_potential_set_kernels(struct boltz_potential *up, int poly_order, boltz_reference_t reference)
{
  up->boundary[0] = boltz_boundary_lower_kernels[up->cdim - 1][poly_order - 1];
  up->boundary[1] = boltz_boundary_upper_kernels[up->cdim - 1][poly_order - 1];
  up->phi = boltz_phi_kernels[up->cdim - 1][poly_order - 1];
  up->reference = reference;
}

/**
 * Initialize shared host state; the model constructor installs kernel pointers.
 *
 * @param up Shared state to initialize.
 * @param basis Configuration basis.
 * @param mass_e Electron mass.
 * @param charge_e Signed electron charge.
 * @param temp_e Electron temperature in energy units.
 * @param use_gpu Whether to operate on device arrays.
 */
void boltz_potential_init(
  struct boltz_potential *up, const struct gkyl_basis *basis, double mass_e, double charge_e,
  double temp_e, bool use_gpu
);

/**
 * Release the device copy, if present; the enclosing updater owns the host state.
 *
 * @param up Shared state.
 */
void boltz_potential_release(struct boltz_potential *up);

/**
 * Apply the shared boundary traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param edge Parallel boundary edge.
 * @param skin Skin range retaining volume-array strides.
 * @param flux Nodal charge-number flux.
 * @param area Nodal projected face-area Jacobian.
 * @param m0_elc Physical electron density in the volume basis.
 * @param wall Volume-modal wall potential, or NULL for zero.
 * @param end_vals Packed nodal boundary output.
 */
void boltz_potential_boundary(
  const struct boltz_potential *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
);

#ifdef GKYL_HAVE_CUDA
/**
 * Apply the shared boundary traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param edge Parallel boundary edge.
 * @param skin Skin range retaining volume-array strides.
 * @param flux Nodal charge-number flux.
 * @param area Nodal projected face-area Jacobian.
 * @param m0_elc Physical electron density in the volume basis.
 * @param wall Volume-modal wall potential, or NULL for zero.
 * @param end_vals Packed nodal boundary output.
 */
void boltz_potential_boundary_cu(
  const struct boltz_potential *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
);
#endif

/**
 * Apply the shared reference traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param ghost Output lower ghost slab retaining array strides.
 * @param lower Packed values at the lower or selected end.
 * @param upper Packed upper values; pass lower again for the single-end model.
 * @param nphi_ref Output modal reference density and potential.
 */
void boltz_potential_reference(
  const struct boltz_potential *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
);

#ifdef GKYL_HAVE_CUDA
/**
 * Apply the shared reference traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param ghost Output lower ghost slab retaining array strides.
 * @param lower Packed values at the lower or selected end.
 * @param upper Packed upper values; pass lower again for the single-end model.
 * @param nphi_ref Output modal reference density and potential.
 */
void boltz_potential_reference_cu(
  const struct boltz_potential *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
);
#endif

/**
 * Apply the shared phi_calc traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param local Interior range retaining array strides.
 * @param local_ext Extended range including the reference ghost slab.
 * @param m0_elc Physical electron density in the volume basis.
 * @param nphi_ref Reference data in the lower parallel ghost slab.
 * @param phi Output DG potential.
 */
void boltz_potential_phi_calc(
  const struct boltz_potential *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
);

#ifdef GKYL_HAVE_CUDA
/**
 * Apply the shared phi_calc traversal using the selected kernels.
 *
 * @param up Shared state.
 * @param local Interior range retaining array strides.
 * @param local_ext Extended range including the reference ghost slab.
 * @param m0_elc Physical electron density in the volume basis.
 * @param nphi_ref Reference data in the lower parallel ghost slab.
 * @param phi Output DG potential.
 */
void boltz_potential_phi_calc_cu(
  const struct boltz_potential *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
);
#endif
