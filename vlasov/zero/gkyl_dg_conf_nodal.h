#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_dg_conf_nodal gkyl_dg_conf_nodal;

/**
 * Create an updater that converts DG expansions to and from a form that is nodal in
 * configuration space and modal in velocity space, for operators that act in velocity space
 * alone at each configuration-space node (collisions).
 *
 * The configuration basis is p=1 and the nodes are its 2^cdim interior Gauss-Legendre points
 * x_d = -+1/sqrt(3), in the row-major order of gkyl_range (the last direction fastest): node
 * q = sum_d i_d 2^(cdim-1-d), with i_d = 0 (1) at the lower (upper) ordinate of direction d.
 * This is the order of the quadrature nodes of the projection updaters and of the nodes of every
 * kernel and updater that stores values per configuration node.
 *
 * A configuration-space expansion becomes its 2^cdim node values. A phase-space expansion, in a
 * basis that is the tensor product of the configuration basis and a velocity basis of NV
 * functions (the hybrid basis), becomes fnodal[q*NV + j]: the velocity mode j at node q, so each
 * node holds a contiguous velocity-space expansion. Nodal arrays have the same number of
 * components as the modal ones, and arrays with several expansions per cell (blocks) are
 * converted block by block.
 *
 * At p=1 the weak product (or division) of configuration-space expansions has the pointwise
 * product (or quotient) of their node values, so operators with configuration-dependent
 * coefficients can be applied node by node.
 *
 * @param conf_basis Configuration-space basis (p=1 serendipity).
 * @param phase_basis Phase-space basis (hybrid), or 0 if only configuration-space arrays are
 *   converted.
 * @param use_gpu Whether to run on the GPU (not available yet).
 * @return New updater pointer.
 */
struct gkyl_dg_conf_nodal *gkyl_dg_conf_nodal_new(
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu
);

/**
 * Convert configuration-space expansions to their node values.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param fmodal Input expansions (any number of blocks).
 * @param fnodal Output node values.
 */
void gkyl_dg_conf_nodal_conf_to_nodal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fmodal, struct gkyl_array *fnodal
);

/**
 * Convert configuration-space node values to expansions.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range.
 * @param fnodal Input node values (any number of blocks).
 * @param fmodal Output expansions.
 */
void gkyl_dg_conf_nodal_conf_to_modal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fnodal, struct gkyl_array *fmodal
);

/**
 * Convert phase-space expansions to the configuration-nodal, velocity-modal form.
 *
 * @param up Updater.
 * @param phase_range Phase-space range.
 * @param fmodal Input expansions (any number of blocks).
 * @param fnodal Output nodal form.
 */
void gkyl_dg_conf_nodal_phase_to_nodal(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *phase_range,
  const struct gkyl_array *fmodal, struct gkyl_array *fnodal
);

/**
 * Accumulate the expansions of a configuration-nodal, velocity-modal array onto phase-space
 * expansions.
 *
 * @param up Updater.
 * @param phase_range Phase-space range.
 * @param fnodal Input nodal form (any number of blocks).
 * @param fmodal Expansions to increment.
 */
void gkyl_dg_conf_nodal_phase_to_modal_accumulate(
  const struct gkyl_dg_conf_nodal *up, const struct gkyl_range *phase_range,
  const struct gkyl_array *fnodal, struct gkyl_array *fmodal
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_conf_nodal_release(struct gkyl_dg_conf_nodal *up);
