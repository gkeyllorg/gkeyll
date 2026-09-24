// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mat.h>
#include <gkyl_mat_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 
#include <gkyl_util.h>
#include <assert.h>

GKYL_CU_DH
static inline void
log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}

struct gkyl_loss_cone_mask_gyrokinetic {
  int cdim; // Configuration-space dimension.
  int pdim; // Phase-space dimension.
  int vdim; // Velocity-space dimension.

  const struct gkyl_rect_grid *grid_phase;
  int num_basis_conf; // Number of configuration-space basis functions.
  int num_basis_phase; // Number of phase-space basis functions.
  double norm_fac; // Normalization factor.

  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  double mass; // Species mass.
  double charge; // Species charge.
  double *bmag_max; // Maximum magnetic field amplitude.
  double *bmag_max_loc; // Location of bmag_max.
  bool use_gpu; // Boolean if we are performing projection on device.

  loss_cone_mask_gyrokinetic_c2p_t c2p_pos; // Function transforming position comp to phys coords.
  void *c2p_pos_ctx; // Context for the c2p_pos mapping.

  bool cellwise_trap_loss; // Whether a whole cell is trapped/lost, or whether
                           // high-order distinction within a cell is allowed.
  struct gkyl_range conf_qrange; // Range of Configuration-space ordinates.
  struct gkyl_range phase_qrange; // Range of Phase-space ordinates.

  // For quadrature in phase-space.
  int tot_quad_phase; // Total number of quadrature points.
  struct gkyl_array *ordinates_phase; // Ordinates.
  struct gkyl_array *weights_phase; // Weights.
  struct gkyl_array *basis_at_ords_phase; // Basis functions at ordinates.

  // For quadrature in configuration-space.
  int tot_quad_conf; // Total number of quadrature points.
  struct gkyl_array *ordinates_conf; // Ordinates.
  struct gkyl_array *weights_conf; // Weights.
  struct gkyl_array *basis_at_ords_conf; // Basis functions at ordinates.

  struct gkyl_array *fun_at_ords; // Mask we are projecting at ordinates in a cell.

  int *p2c_qidx;  // Mapping between configuration-space and phase-space ordinates.
  struct gkyl_array *mask_out_quad; // Array keeping f_lte at phase-space quadrature nodes.
  struct gkyl_array *qDphiDbmag_quad; // Array keeping q*(phi-phi_m)/(B_max-B)
                                      // at configuration-space quadrature nodes.
  struct gkyl_array *Dbmag_quad; // B_max-B at configuration-space quadrature nodes.

  struct gkyl_mat_mm_array_mem *phase_nodal_to_modal_mem; // Structure of data which converts  
                                                          // stores the info to convert phase
                                                          // space nodal to modal gkyl arrays.
};

#ifdef GKYL_HAVE_CUDA
/**
 * Obtain bmag_max-bmag at conf-space quadrature nodes and store it in up->Dbmag_quad.
 *
 * @param up Project on basis updater to run.
 * @param conf_rng Configuration-space range.
 * @param bmag Magnetic field magnitude.
 * @param bmag_max Maximum bmag.
 */
void 
gkyl_loss_cone_mask_gyrokinetic_Dbmag_quad_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *conf_range, const struct gkyl_array *bmag, const double *bmag_max);

/**
 * Compute projection of the loss cone masking function on the phase-space basis
 * on the GPU.
 *
 * @param up Project on basis updater to run.
 * @param phase_rng Phase-space range.
 * @param conf_rng Configuration-space range.
 * @param phi Electrostatic potential.
 * @param phi_m Electrostatic potential at the mirror throat (on GPU).
 * @param mask_out Output masking function.
 */
void
gkyl_loss_cone_mask_gyrokinetic_advance_cu(gkyl_loss_cone_mask_gyrokinetic *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *phi, const double *phi_m, struct gkyl_array *mask_out);
#endif
