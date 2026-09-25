// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH static inline void
comp_to_phys(
  int ndim, const double *eta, const double *GKYL_RESTRICT dx, const double *GKYL_RESTRICT xc,
  double *GKYL_RESTRICT xout
)
{
  for (int d = 0; d < ndim; ++d) {
    xout[d] = 0.5 * dx[d] * eta[d] + xc[d];
  }
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i = 0; i < cdim; ++i) {
    out[i] = cidx[i];
  }
  for (int i = cdim; i < pdim; ++i) {
    out[i] = vidx[i - cdim];
  }
}

struct gkyl_vlasov_lte_proj_on_basis {
  struct gkyl_rect_grid phase_grid;
  int cdim; // Configuration-space dimension
  int pdim; // Phase-space dimension
  int vdim; // Velocity-space dimension

  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis vmap_basis; // Velocity-space basis for vmap.
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // number of Configuration-space basis functions
  int num_phase_basis; // number of Phase-space basis functions

  struct gkyl_basis *vmap_basis_on_dev; // Velocity-space basis for vmap on GPUs.

  int num_comp_V_drift; // Number of components in V_drift (in relativity this is the four-velocity).
  bool is_relativistic; // Boolean for if we are projecting the relativistic LTE
  bool is_canonical_pb; // Boolean for if we are projecting the canonical-pb LTE
  bool is_bimaxwellian; // Boolean for if we are projecting a bi-Maxwellian
  bool use_gpu; // Boolean if we are performing projection on device.

  struct gkyl_range conf_qrange; // Range of Configuration-space ordinates.
  struct gkyl_range vel_qrange; // Range of Velocity-space ordinates.
  struct gkyl_range phase_qrange; // Range of Phase-space ordinates.

  struct gkyl_range vel_range;

  // for quadrature in Phase-space
  int tot_quad; // total number of Phase-space quadrature points
  struct gkyl_array *ordinates; // Phase-space ordinates for quadrature
  struct gkyl_array *weights; // Phase-space weights for quadrature
  struct gkyl_array *basis_at_ords; // Phase-space basis functions at ordinates

  // for quadrature in Configuration-space
  int tot_conf_quad; // total number of Configuration-space quadrature points
  struct gkyl_array *conf_ordinates; // Configuration-space ordinates for quadrature
  struct gkyl_array *conf_weights; // weights for Configuration-space quadrature
  struct gkyl_array *conf_basis_at_ords; // Configuration-space basis functions at ordinates

  struct gkyl_array *fun_at_ords; // function LTE distribution evaluated at
    // ordinates in a cell.

  int *p2c_qidx; // Mapping between Configuration-space and Phase-space ordinates.
  int *p2v_qidx; // Mapping between Velocity-space and Phase-space ordinates.
  struct gkyl_array *f_lte_quad; // Array keeping f_lte at phase-space quadrature nodes
  struct gkyl_array *moms_lte_quad; // Array keeping moms_lte (n, V_drift, T/m)
    // at configuration-space quadrature nodes
  struct gkyl_array *expamp_quad; // Array keeping exponential pre-factor
    // at configuration-space quadrature nodes
  struct gkyl_array *h_ij_quad; // metric configuration-space quadrature nodes, covariant components
  struct gkyl_array *
    h_ij_inv_quad; // metric inverse at configuration-space quadrature nodes, contravariant components
  struct gkyl_array *det_h_quad; // metric determinant at configuration-space quadrature nodes

  bool
    use_extended_hamil_def; // (Can-pb) quantity, bool to determine if we wish to use an extended can-pb Hamil
  struct gkyl_array *
    background_flows_quad; // (Can-bp quantity) background flows from the geometry (i.e. constant rotation)
  struct gkyl_array *
    effective_potential_quad; // (Can-bp quantity) background effective potential from the geometry (i.e. constant rotation effective potential)

  bool use_vmap;
  const struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object (acquired; 0 if not given).
  struct gkyl_array *vmap; // Borrowed from vel_map; kept alive by the acquired struct.
  struct gkyl_array *jacob_vel_gauss; // Borrowed from vel_map; kept alive by the acquired struct.

  struct gkyl_vlasov_lte_moments *moments_up; // LTE moment calculation routine for computing density
  struct gkyl_array *num_ratio; // Number density ratio: num_ratio = n_target/n0
  struct gkyl_dg_bin_op_mem *mem; // bin_op memory to compute ratio and rescale distribution function

  struct gkyl_mat_mm_array_mem *phase_nodal_to_modal_mem; // structure of data which converts
    // stores the info to convert phase
    // space nodal to modal gkyl arrays
};
