#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

// Forward declaration to avoid pulling in comm/IO headers.
struct gkyl_comm;

struct gkyl_vlasov_position_map_inp {
  evalf_t eval_pmap; // Position mapping in this configuration-space dimension. NULL => identity.
  void *ctx; // Context for function evaluation. Can be NULL.
};

// Configuration-space mapping object for the Vlasov app. Packages a
// (possibly non-uniform) computational-to-physical position mapping in each
// configuration dimension, along with the configuration-space Jacobian
// evaluated at the Gauss-Legendre quadrature points needed by the streaming
// (spatial advection) volume and surface operations. This object is the
// configuration-space analogue of gkyl_vlasov_velocity_map and is *always*
// created by the app: when no user mapping is supplied the identity map is
// used, and there is a single code path — the map is always populated, always
// written at frame 0, and the Jacobian is always divided/rescaled — so
// consumers handle uniform and non-uniform configuration grids transparently.
//
// Unlike the velocity map, the position map supports only a C^0 piecewise
// LINEAR representation (for both Serendipity and tensor bases): only the
// tensor basis can divide out a continuous quadratic Jacobian without aliasing
// errors, so a higher-order conf map is deferred. The linear map is stored as a
// degenerate cubic in the same cdim*4 layout the velocity map uses (quadratic
// and cubic coefficients identically zero), so the post-processing/I/O tooling
// is shared. The conf Jacobian J_x is therefore piecewise constant and
// discontinuous at cell interfaces: the Jacobian arrays hold that constant
// replicated at every quadrature point, and consumers needing the Jacobian at
// an interface must take the minimum of the two adjacent cells when estimating
// maximum frequencies. Because J_x is a per-cell constant, dividing/rescaling
// the (J_x J_v f)-weighted distribution by J_x is an exact scalar operation,
// performed directly in this file with no generated nodal kernels.
//
// This is a host-side, reference-counted container; like the velocity map there
// is no device-side mirror of the struct. When created with use_gpu=true the
// solver arrays (pmap, jacob_pos, jacob_pos_surf, jacob_pos_gauss) are
// device-resident gkyl_arrays, and consumers acquire the struct (for lifetime
// safety) and unpack raw device pointers inside their _cu.cu initialization.
struct gkyl_vlasov_position_map {
  struct gkyl_rect_grid grid_pos; // Configuration-space grid.
  struct gkyl_range
    local_pos; // Interior configuration-space range (I/O sub-range; conf indexing of the local phase range).
  struct gkyl_range
    local_ext_pos; // Extended configuration-space range the arrays are actually defined on (includes ghost cells).
  struct gkyl_basis basis_pos; // Configuration-space basis (b_type and poly_order drive consumers).
  struct gkyl_basis
    basis_pgkyl; // Basis of the I/O representation pmap_pgkyl_host (p=3 serendipity); used to build the write metadata.
  bool is_identity; // True if no user mapping was given in any direction.

  // Solver arrays; device-resident when created with use_gpu=true.
  struct gkyl_array *
    pmap; // Linear map (degenerate-cubic layout) in each config dimension (cdim*4 components: x, then y, then z).
  struct gkyl_array *
    jacob_pos; // Jacobian (derivative of pmap) in each config dimension at 1D Gauss-Legendre quadrature points (cdim*(p+1) components).
  struct gkyl_array *
    jacob_pos_surf; // Jacobian in each config dimension at the (higher order) 1D Gauss-Legendre quadrature points used by surface updates (cdim*(p+2) components).
  struct gkyl_array *
    jacob_pos_gauss; // Total configuration-space Jacobian at the full Gauss-Legendre quadrature points (tensor(cdim,p).num_basis components). Constant per cell.

  // Host mirrors of the solver arrays (acquired aliases when !use_gpu).
  struct gkyl_array *pmap_host;
  struct gkyl_array *jacob_pos_host;
  struct gkyl_array *jacob_pos_surf_host;
  struct gkyl_array *jacob_pos_gauss_host;

  // I/O-only arrays; always host-resident.
  struct gkyl_array
    *pmap_pgkyl_host; // Map for I/O (defined in the full 1X, 2X, or 3X p=3 serendipity layout).
  struct gkyl_array *pmap_avg_pgkyl_host; // Cell average of the mapping for I/O.

  uint32_t flags;
  struct gkyl_ref_count ref_count;
};

/**
 * Create a new configuration-space position map based on the input mapping
 * functions, using a C^0 piecewise linear representation in each configuration
 * dimension (for both Serendipity and tensor bases). Directions with
 * eval_pmap==NULL get the identity map, so this constructor is valid for every
 * app. All arrays (including host mirrors and I/O arrays) are allocated and
 * populated here; when use_gpu=true the solver arrays are copied onto device.
 *
 * @param cgrid Configuration-space grid object.
 * @param crange Configuration-space interior range (used for I/O and conf indexing).
 * @param crange_ext Extended configuration-space range the arrays are defined on (ghost cells copy the adjacent skin cell).
 * @param conf_basis Configuration-space basis.
 * @param inp_pmap Position mapping input (function and context) in each configuration-space dimension.
 * @param use_gpu Whether to place the solver arrays on device.
 * @return New position map object.
 */
struct gkyl_vlasov_position_map *gkyl_vlasov_position_map_new(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange,
  const struct gkyl_range *crange_ext, const struct gkyl_basis *conf_basis,
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM], bool use_gpu
);

/**
 * Acquire a pointer to the position map object, incrementing its reference
 * count. Release with gkyl_vlasov_position_map_release.
 *
 * @param vpm Position map object.
 * @return Acquired position map object.
 */
struct gkyl_vlasov_position_map *gkyl_vlasov_position_map_acquire(
  const struct gkyl_vlasov_position_map *vpm
);

/**
 * Release the position map object, decrementing its reference count and
 * freeing all contained arrays when the count reaches zero.
 *
 * @param vpm Position map object.
 */
void gkyl_vlasov_position_map_release(const struct gkyl_vlasov_position_map *vpm);

/**
 * Check whether the solver arrays in this position map are device-resident.
 *
 * @param vpm Position map object.
 * @return True if the solver arrays live on device.
 */
bool gkyl_vlasov_position_map_is_cu_dev(const struct gkyl_vlasov_position_map *vpm);

/**
 * Evaluate the (DG) computational-to-physical configuration-space map at a
 * computational coordinate. Used as the c2p coordinate transform when
 * projecting initial conditions / sources on a non-uniform conf mesh, so the
 * projection samples the user function at the same physical coordinates the
 * solver kernels assume. Host-side only (uses the host map array).
 *
 * @param vpm Position map object.
 * @param xc Computational configuration-space coordinate (cdim components).
 * @param xp On output, the physical configuration-space coordinate (cdim components).
 */
void gkyl_vlasov_position_map_eval_mc2p(
  const struct gkyl_vlasov_position_map *vpm, const double *xc, double *xp
);

/**
 * Write the position map to %s-%s_pmap.gkyl and its cell average to
 * %s-%s_pmap_avg.gkyl, on rank 0 of comm. The map is static in time so this is
 * intended to be called once (e.g. at frame 0); it is always written, uniform
 * (identity) grids included, and both files are always written so cell-average
 * post-processing has a matching map. The metadata is constructed internally
 * from the stored I/O basis — the files hold the p=3 serendipity pgkyl
 * representation (p=0 for the cell average), not the app's conf basis, so
 * caller-supplied app metadata would mislabel the data.
 *
 * @param vpm Position map object.
 * @param comm Communicator the app belongs to.
 * @param app_name Name of the app.
 * @param name Name component used in the file name (e.g. "position-map").
 */
void gkyl_vlasov_position_map_write(
  const struct gkyl_vlasov_position_map *vpm, struct gkyl_comm *comm, const char *app_name,
  const char *name
);

/**
 * Divide out the configuration-space Jacobian from Jf to obtain f. The conf
 * Jacobian is a per-cell constant (C^0 linear map), so the division is exact
 * both modally and nodally: each phase-space coefficient is divided by the
 * total conf Jacobian of its configuration cell. Runs on device when the map
 * was created with use_gpu=true.
 *
 * @param vpm Position map object.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param Jf Input array carrying the conf Jacobian weight.
 * @param f_no_J Output array with configuration-space Jacobian divided out.
 */
void gkyl_vlasov_position_map_divide_jacobpos(
  const struct gkyl_vlasov_position_map *vpm, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range, const struct gkyl_array *Jf, struct gkyl_array *f_no_J
);

/**
 * Multiply f by the configuration-space Jacobian to obtain Jf; the inverse of
 * gkyl_vlasov_position_map_divide_jacobpos.
 *
 * @param vpm Position map object.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param f_no_J Input array with configuration-space Jacobian divided out.
 * @param Jf Output array carrying the conf Jacobian weight.
 */
void gkyl_vlasov_position_map_rescale_jacobpos(
  const struct gkyl_vlasov_position_map *vpm, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range, const struct gkyl_array *f_no_J, struct gkyl_array *Jf
);

/**
 * Divide the leading components of a configuration-space field (e.g. a velocity
 * moment) by the per-cell constant total conf Jacobian J. Moments of the stored
 * J_x J_v f carry J (the J_x factors out of the velocity integral), so this
 * converts them to physical moments for I/O. Operates on HOST arrays over the
 * configuration range (moment I/O is host-side).
 *
 * Only the first num_coeff_divide coefficients of each cell are scaled by 1/J;
 * the remaining coefficients are copied through unchanged. This supports
 * derived moment diagnostics (e.g. the LTE moments n, V_drift, T/m) where only
 * the extensive density n (the first num_basis coefficients) carries J while
 * the intensive V_drift and T/m are velocity ratios in which J cancels. For
 * pure moments (M0, M1i, M2, ...) every component carries J, so pass the full
 * component count (mom->ncomp).
 *
 * @param vpm Position map object.
 * @param conf_range Configuration-space range to update.
 * @param num_coeff_divide Number of leading coefficients per cell to divide by J.
 * @param Jmom Input configuration-space field carrying the conf Jacobian weight.
 * @param mom_no_J Output configuration-space field with the conf Jacobian divided out.
 */
void gkyl_vlasov_position_map_divide_jacobpos_conf(
  const struct gkyl_vlasov_position_map *vpm, const struct gkyl_range *conf_range,
  int num_coeff_divide, const struct gkyl_array *Jmom, struct gkyl_array *mom_no_J
);

/**
 * Multiply every component of a configuration-space field by the per-cell
 * constant total conf Jacobian J (host-side). Inverse of the full-component
 * divide; used to convert a physical conf field to the J-weighted form used for
 * evolution, e.g. E -> J*E for the Vlasov-Maxwell field at IC/restart.
 *
 * @param vpm Position map object.
 * @param conf_range Configuration-space range to update.
 * @param a_no_J Input physical configuration-space field.
 * @param Ja Output J-weighted configuration-space field.
 */
void gkyl_vlasov_position_map_rescale_jacobpos_conf(
  const struct gkyl_vlasov_position_map *vpm, const struct gkyl_range *conf_range,
  const struct gkyl_array *a_no_J, struct gkyl_array *Ja
);

/**
 * Device versions of the configuration-space divide/rescale by the per-cell
 * position-map Jacobian. jac_range is the range the Jacobian array is indexed
 * with (the map's local_pos); only cells in conf_range are written.
 */
void gkyl_vlasov_position_map_divide_jacobpos_conf_cu(
  const struct gkyl_range *jac_range, const struct gkyl_range *conf_range, int num_coeff_divide,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *Jmom,
  struct gkyl_array *mom_no_J
);

void gkyl_vlasov_position_map_rescale_jacobpos_conf_cu(
  const struct gkyl_range *jac_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *jacob_pos_gauss, const struct gkyl_array *a_no_J, struct gkyl_array *Ja
);
