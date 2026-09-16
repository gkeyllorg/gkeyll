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

struct gkyl_vlasov_velocity_map_inp {
  evalf_t eval_vmap; // Velocity mapping in this velocity-space dimension used by dg_basis_ops to construct C^1 mapping. NULL => identity.
  void *ctx; // Context for function evaluation. Can be NULL.
};

// Representation of the stored velocity map, chosen by the velocity basis:
// tensor bases use the C^1 cubic (the tensor p=1 hybrid phase basis hands
// the map its p=2 tensor velocity basis, so it lands here too); Serendipity
// bases use the C^0 piecewise linear, stored as a degenerate cubic in the
// same vdim*4 layout (quadratic and cubic coefficients zero). The C^0 map's
// Jacobian is piecewise constant and
// discontinuous at cell interfaces: the Jacobian arrays hold that constant
// replicated at every quadrature point, and consumers needing the Jacobian
// at an interface must take the minimum of the two adjacent cells when
// estimating maximum frequencies.
enum gkyl_vlasov_vmap_rep {
  GKYL_VLASOV_VMAP_C1_CUBIC = 0,
  GKYL_VLASOV_VMAP_C0_LINEAR,
};

// Velocity-space mapping object for the Vlasov app. Packages the
// representation (C^1 cubic or C^0 linear, by velocity basis) of the
// (possibly non-uniform) computational-to-physical velocity mapping in each
// velocity dimension, along with the velocity-space Jacobian evaluated at the
// Gauss-Legendre quadrature points needed by volume, surface, and projection
// operations. This object is *always* created by the app: when no user
// mapping is supplied the identity map is used, and there is a single code
// path — the map is always populated, always written at frame 0, and the
// Jacobian is always divided/rescaled — so consumers handle uniform and
// non-uniform velocity grids transparently.
//
// This is a host-side, reference-counted container. There is deliberately no
// device-side mirror of the struct itself: when created with use_gpu=true the
// solver arrays (vmap, jacob_vel, jacob_vel_surf, jacob_vel_gauss) are
// device-resident gkyl_arrays, and consumers acquire the struct (for lifetime
// safety) and unpack raw device pointers (e.g. vvm->jacob_vel->on_dev) inside
// their _cu.cu initialization.
struct gkyl_vlasov_velocity_map {
  struct gkyl_rect_grid grid_vel; // Velocity-space grid.
  struct gkyl_range local_vel; // Velocity-space range the arrays are defined on.
  struct gkyl_basis basis_vel; // Velocity-space basis (b_type and poly_order drive kernel dispatch).
  struct gkyl_basis basis_pgkyl; // Basis of the I/O representation vmap_pgkyl_host (p=3 serendipity); used to build the write metadata.
  enum gkyl_vlasov_vmap_rep rep; // Representation of the stored map.
  bool is_identity; // True if no user mapping was given in any direction.
  bool is_mapped; // Always true: the map is always created and consumed.
                  // Retained for consumers that set use_vmap from it; to be
                  // removed once those call sites are simplified.

  // Solver arrays; device-resident when created with use_gpu=true.
  struct gkyl_array *vmap; // C^1 cubic representation of mapping in each velocity dimension (vdim*4 components: vx, then vy, then vz).
  struct gkyl_array *jacob_vel; // Jacobian (derivative of vmap) in each velocity dimension at 1V Gauss-Legendre quadrature points (vdim*(p+1) components, p the velocity basis poly order).
  struct gkyl_array *jacob_vel_surf; // Jacobian in each velocity dimension at the 1V Gauss-Legendre quadrature points used by surface updates: the count matches the active kernel variant (p+2 nodes, except the tensor p=1 hybrid's low-order variant, which uses p+1 = 3 nodes).
  struct gkyl_array *jacob_vel_gauss; // Total velocity-space Jacobian at the full Gauss-Legendre quadrature points (tensor(vdim,p).num_basis components).

  // Host mirrors of the solver arrays (acquired aliases when !use_gpu).
  struct gkyl_array *vmap_host;
  struct gkyl_array *jacob_vel_host;
  struct gkyl_array *jacob_vel_surf_host;
  struct gkyl_array *jacob_vel_gauss_host;

  // I/O-only arrays; always host-resident.
  struct gkyl_array *vmap_pgkyl_host; // C^1 cubic representation of mapping for I/O (defined in the full 1V, 2V, or 3V).
  struct gkyl_array *vmap_avg_pgkyl_host; // Cell average of the mapping for I/O.

  uint32_t flags;
  struct gkyl_ref_count ref_count;
};

/**
 * Create a new velocity map based on the input mapping functions, using a
 * C^1 cubic representation for tensor velocity bases and a C^0 piecewise
 * linear representation for Serendipity velocity bases. The tensor p=1
 * hybrid phase basis is p=2 in velocity space, so hybrid species pass their
 * p=2 tensor velocity basis here and get the cubic representation.
 * Directions with eval_vmap==NULL get the identity map, so this constructor
 * is valid for every species and model. All arrays (including host mirrors
 * and I/O arrays) are allocated and populated here; when use_gpu=true the
 * solver arrays are copied onto device.
 *
 * @param vgrid Velocity-space grid object.
 * @param vrange Velocity-space range the arrays are defined on.
 * @param vel_basis Velocity-space basis (b_type and poly_order drive kernel dispatch).
 * @param inp_vmap Velocity mapping input (function and context) in each velocity-space dimension.
 * @param use_lo Whether the run uses the low-order surface kernel variant;
 *   with a tensor p=2 velocity basis (the tensor p=1 hybrid) this stores
 *   jacob_vel_surf at 3 surface-quadrature nodes instead of 4. Ignored for
 *   every other basis (only the hybrid has distinct lo/ho variants).
 * @param use_gpu Whether to place the solver arrays on device.
 * @return New velocity map object.
 */
struct gkyl_vlasov_velocity_map* gkyl_vlasov_velocity_map_new(
  const struct gkyl_rect_grid *vgrid, const struct gkyl_range *vrange,
  const struct gkyl_basis *vel_basis,
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM],
  bool use_lo, bool use_gpu);

/**
 * Acquire a pointer to the velocity map object, incrementing its reference
 * count. Release with gkyl_vlasov_velocity_map_release.
 *
 * @param vvm Velocity map object.
 * @return Acquired velocity map object.
 */
struct gkyl_vlasov_velocity_map* gkyl_vlasov_velocity_map_acquire(
  const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Release the velocity map object, decrementing its reference count and
 * freeing all contained arrays when the count reaches zero.
 *
 * @param vvm Velocity map object.
 */
void gkyl_vlasov_velocity_map_release(const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Check whether the solver arrays in this velocity map are device-resident.
 *
 * @param vvm Velocity map object.
 * @return True if the solver arrays live on device.
 */
bool gkyl_vlasov_velocity_map_is_cu_dev(const struct gkyl_vlasov_velocity_map *vvm);

/**
 * Evaluate the (DG) computational-to-physical velocity map at a computational
 * velocity coordinate. Used as the velocity part of the c2p coordinate
 * transform when projecting phase-space initial conditions, so the projection
 * samples the user function at the same physical velocities the solver kernels
 * assume (the stored map is the cubic/linear approximation the operators use,
 * not the exact analytic eval_vmap). Host-side only.
 *
 * @param vvm Velocity map object.
 * @param vc Computational velocity coordinate (vdim components).
 * @param vp On output, the physical velocity coordinate (vdim components).
 */
void gkyl_vlasov_velocity_map_eval_c2p(const struct gkyl_vlasov_velocity_map *vvm,
  const double *vc, double *vp);

/**
 * Write the velocity map to %s-%s_vmap.gkyl and its cell average to
 * %s-%s_vmap_avg.gkyl, on rank 0 of comm. The map is static in time so this is
 * intended to be called once (e.g. at frame 0); it is always written, uniform
 * (identity) grids included, and both files are always written so cell-average
 * post-processing has a matching map. The metadata is constructed internally
 * from the stored I/O basis — the files hold the p=3 serendipity pgkyl
 * representation (p=0 for the cell average), not the species' velocity basis,
 * so caller-supplied app metadata would mislabel the data.
 *
 * @param vvm Velocity map object.
 * @param comm Communicator of the species this map belongs to.
 * @param app_name Name of the app.
 * @param species_name Name of the species.
 */
void gkyl_vlasov_velocity_map_write(const struct gkyl_vlasov_velocity_map *vvm,
  struct gkyl_comm *comm, const char *app_name, const char *species_name);

/**
 * Divide out the velocity-space Jacobian from Jf to obtain f for use in the
 * DG update. The division is performed at Gauss-Legendre quadrature points
 * because, by construction, we know Jv at these points and Jv is static in
 * time, so we obtain the weakly equivalent f from Jf in our DG expansion
 * (since we can always transform between the Gauss-Legendre nodal basis and
 * our modal basis). For the C^0 linear representation (Serendipity bases) the
 * Jacobian is constant in the cell so the division is exact, modally and
 * nodally. Runs on device when the map was created with use_gpu=true.
 *
 * @param vvm Velocity map object.
 * @param conf_basis Configuration-space basis.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param Jf Input array Jf.
 * @param f_no_J Output array f with velocity-space Jacobian divided out.
 */
void gkyl_vlasov_velocity_map_divide_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J);

/**
 * Multiply f by the velocity-space Jacobian to obtain Jf for use in the DG
 * update. Performed at Gauss-Legendre quadrature points; the inverse of
 * gkyl_vlasov_velocity_map_divide_jacobvel.
 *
 * @param vvm Velocity map object.
 * @param conf_basis Configuration-space basis.
 * @param phase_basis Phase-space basis.
 * @param phase_range Phase-space range to update.
 * @param f_no_J Input array f with velocity-space Jacobian divided out.
 * @param Jf Output array Jf.
 */
void gkyl_vlasov_velocity_map_rescale_jacobvel(const struct gkyl_vlasov_velocity_map *vvm,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf);
