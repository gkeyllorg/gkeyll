// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_flux_vlasov_kernels.h>
#include <gkyl_vlasov_surf_node_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>
#include <gkyl_vlasov_hyb_build.h>

// Force producer in sum-factorized form. The force at surface node (outer
// transverse-configuration node i, inner velocity node j) is
//   alpha(i, j) = sum_t O[t*NO + i] * I[t*NI + j]
// over the producer's terms t. _shared is called once per surface-node thread
// tid: it fills the outer factors O[(off + t)*NO + tid] if tid < NO and the
// inner factors I[(off + t)*NI + tid] if tid < NI (max(NO, NI) <= NO*NI, so one
// thread per surface node covers it), off being where its terms start, and
// returns its number of terms; called with
// O == NULL it only returns the term count (used to size the buffers).
// hamil_pt_edge selects the face (-1: this cell's lower face, +1: upper).
typedef int (*hamil_alpha_shared_conf_t)(
  int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
);

// Nodal Lax-Friedrichs flux in three stages. Stages 1 and 3 are cooperative:
// node thread tid does item tid (the launcher gives one thread per surface
// node, the CPU dispatch loops over them; items beyond a stage's count are
// no-ops, and there are never more items than surface nodes). Stage 1 of the
// sum-factorized nodal evaluation of f: item = (outer shape, inner node) pair,
// fills G_l[j*NA + a],
// G_r[j*NA + a] from f_l, f_r.
typedef void (*lax_g_t)(
  int tid, const double *f_l, const double *f_r, double *GKYL_RESTRICT G_l,
  double *GKYL_RESTRICT G_r
);

// Stage 2: the Lax flux at node (i, j) from the filled G arrays, written into
// the nodal buffer (num_nodes entries); returns |alpha| for the CFL reduction.
typedef double (*lax_flux_nodal_t)(
  int i, int j, const double *jacob_pos_l, const double *jacob_pos_r, double alpha,
  const double *G_l, const double *G_r, double *GKYL_RESTRICT Fhat_nodal
);

// Stage 3: projection of the nodal buffer onto the surface modal basis, stored
// at the direction's modal offset of the flux array; item = one inner surface
// mode (or one surface mode, by the kernel's choice).
typedef void (*lax_prj_t)(int tid, const double *Fhat_nodal, double *GKYL_RESTRICT flux);

// CFL estimate of the surface from the reduced alpha_max.
typedef double (*lax_cfl_t)(
  const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max
);

typedef double (*conf_flux_surf_t)(
  struct gkyl_dg_vlasov_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int hamil_pt_edge, const double *vmap, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  const double *f_l, const double *f_r, double *GKYL_RESTRICT conf_flux_surf
);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D static struct {
  int vdim[4];
} cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1, 0, 1, 2}, // 1x kernel indices
  {-1, 3, 4, 5}, // 2x kernel indices
  {-1, -1, -1, 6} // 3x kernel indices
};

// for use in kernel tables
typedef struct {
  hamil_alpha_shared_conf_t kernels[4];
} gkyl_hamil_alpha_quad_kern_list;
typedef struct {
  lax_g_t kernels[4];
} gkyl_lax_g_kern_list;
typedef struct {
  lax_flux_nodal_t kernels[4];
} gkyl_lax_flux_nodal_kern_list;
typedef struct {
  lax_prj_t kernels[4];
} gkyl_lax_prj_kern_list;
typedef struct {
  lax_cfl_t kernels[4];
} gkyl_lax_cfl_kern_list;

// Largest force-factor buffer, in doubles, over all kernels: alpha_nterms_max
// * (num_nodes_conf + num_nodes_vel) (3x3v tensor p=1 phase: 12 terms x 68 nodes).
#define GKYL_VLASOV_CONF_FLUX_SURF_MAX_ALPHA_FACTORS 816
// Largest number of surface nodes of the per-node dispatch (2x3v p2 higher-order: 4 x 64).
#define GKYL_VLASOV_FLUX_SURF_MAX_NODES 256
// Node threads per cell that do the first level of the CUDA kernels' alpha_max
// (CFL) reduction; the last node thread folds their partials.
#define GKYL_FLUX_SURF_CFL_REDUCERS 8

struct gkyl_dg_vlasov_conf_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double
    skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  int hamil_dim; // Dimensionality of Hamiltonian.
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index.
  struct gkyl_range
    hamil_range; // Range for indexing Hamiltonian (either Velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in Velocity-space Jacobian.
  const struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array
    *vmap; // Velocity map (4-slot cubic rep per direction; borrowed from vel_map).
  const struct gkyl_array *
    jacob_vel_surf; // Velocity-space Jacobian at surface quadrature points (borrowed from vel_map).
  const struct gkyl_array *
    jacob_pos; // Configuration-space (position-map) Jacobian, per-conf-cell constant (borrowed from pos_map; defined on the extended conf range).
  hamil_alpha_shared_conf_t hamil_alpha_shared[3]; // Hamiltonian force Pi . dH/dv.
  lax_g_t lax_g[3]; // Stage 1 of the nodal f evaluation (fills G_l, G_r).
  lax_flux_nodal_t lax_flux_nodal[3]; // Stage 2: Lax-Friedrichs flux at one surface node.
  lax_prj_t lax_prj[3]; // Stage 3: nodal -> surface-modal projection.
  lax_cfl_t lax_cfl[3]; // CFL estimate from the reduced alpha_max.
  int num_nodes_conf; // Remaining configuration-space surface nodes per configuration surface.
  int num_nodes_vel; // Velocity-space volume nodes per configuration surface.
  int num_surf_basis; // Surface modal basis functions per direction in the stored flux.
  int alpha_nterms_max; // Largest force term count over the directions (buffer sizing).
  conf_flux_surf_t
    conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes.

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_vlasov_conf_flux_surf *on_dev; // pointer to itself or device data.
};

static inline int
conf_flux_surf_num_surf_basis_ipow(int b, int e)
{
  int r = 1;
  for (int i = 0; i < e; ++i) {
    r *= b;
  }
  return r;
}

// Number of surface modal basis functions on a configuration-direction surface:
// the (pdim-1)-dimensional basis of the phase-space kernel type, i.e. the
// serendipity or tensor basis of the same order, or for the tensor p=1
// hybrid (p=1 configuration x p=2 velocity) 2^(cdim-1) * 3^vdim. Matches the
// .nk of the kernels' vst_*_prj_* tables.
static inline int
conf_flux_surf_num_surf_basis(enum gkyl_basis_type kernel_b_type, int cdim, int vdim, int poly_order)
{
  struct gkyl_basis sb;
  if (kernel_b_type == GKYL_BASIS_MODAL_TENSOR) {
    if (poly_order == 1) {
      return conf_flux_surf_num_surf_basis_ipow(2, cdim - 1) *
             conf_flux_surf_num_surf_basis_ipow(3, vdim);
    }
    gkyl_cart_modal_tensor(&sb, cdim + vdim - 1, poly_order);
    return sb.num_basis;
  }
  gkyl_cart_modal_serendip(&sb, cdim + vdim - 1, poly_order);
  return sb.num_basis;
}

// Surface node counts of the per-node dispatch: (p+1) points per direction,
// p+2 for the higher-order (anti-aliasing) and tensor (cubic-map) kernels.
// The tensor p=1 hybrid (p=1 conf x p=2 vel) is anisotropic: 2 nodes per
// configuration direction, 3 (lo) or 4 (ho) per velocity direction.
static inline void
conf_flux_surf_num_nodes(
  enum gkyl_basis_type kernel_b_type, int cdim, int vdim, int poly_order, bool use_lo,
  int *num_nodes_conf, int *num_nodes_vel
)
{
  int nq_conf = poly_order + 1, nq_vel = poly_order + 1;
  if ((poly_order > 1) && !use_lo) {
    nq_conf = poly_order + 2;
    nq_vel = poly_order + 2;
  }
  if (kernel_b_type == GKYL_BASIS_MODAL_TENSOR) {
    if (poly_order == 1) {
      nq_conf = 2;
      nq_vel = use_lo ? 3 : 4;
    } else {
      nq_conf = poly_order + 2;
      nq_vel = poly_order + 2;
    }
  }
  *num_nodes_conf = 1;
  for (int d = 0; d < cdim - 1; ++d) {
    *num_nodes_conf *= nq_conf;
  }
  *num_nodes_vel = 1;
  for (int d = 0; d < vdim; ++d) {
    *num_nodes_vel *= nq_vel;
  }
}

// Total number of force terms in direction dir (size query of the producer).
GKYL_CU_DH static inline int
conf_flux_surf_alpha_nterms(const struct gkyl_dg_vlasov_conf_flux_surf *up, int dir)
{
  return up->hamil_alpha_shared[dir](0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

// Empty function pointer for cases where this force does not exist (no terms).
GKYL_CU_DH static int
no_hamil_alpha_shared(
  int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
)
{
  return 0;
}

// Configuration-space flux of one cell surface (CPU dispatch). This is the
// serial form of the algorithm the CUDA kernel runs with one thread per surface
// node: the shared work of the cell (the force factors and stage 1 of the
// nodal f evaluation, one item per node thread), then the per-node work
// (force by the factor dot product, nodal Lax flux), then the projection onto
// the surface modal basis (one item per node thread). Returns the surface CFL
// estimate.
static double
conf_flux_surf_arrays(
  struct gkyl_dg_vlasov_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int hamil_pt_edge, const double *vmap, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  const double *f_l, const double *f_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  const int NO = up->num_nodes_conf, NI = up->num_nodes_vel, num_nodes = NO * NI;
  double alpha_factors[GKYL_VLASOV_CONF_FLUX_SURF_MAX_ALPHA_FACTORS];
  double *O = alpha_factors, *I = alpha_factors + up->alpha_nterms_max * NO;
  // Sized for the largest surface node count (2x3v p2: 4 x 64 nodes); GKYL_DEFAULT_NUM_THREADS
  // is 1 on non-CUDA builds and must not size these buffers.
  double G_l[GKYL_VLASOV_FLUX_SURF_MAX_NODES], G_r[GKYL_VLASOV_FLUX_SURF_MAX_NODES];
  double Fhat_nodal[GKYL_VLASOV_FLUX_SURF_MAX_NODES];

  // Shared work of the cell, as each surface-node thread of the CUDA kernel
  // does it: its row of the outer/inner factors of the force producer and its
  // item of stage 1 of the nodal f evaluation.
  int nt = 0;
  for (int tid = 0; tid < num_nodes; ++tid) {
    nt = up->hamil_alpha_shared[dir](
      tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos_r, jacob_vel_surf, poisson_tensor_conf, hamil,
      O, I
    );
    up->lax_g[dir](tid, f_l, f_r, G_l, G_r);
  }

  // Per-node work: the force at node (i, j) is the dot product of the outer
  // and inner factors over all terms; nodal Lax flux into the nodal buffer.
  double alpha_max = 0.0;
  for (int n = 0; n < num_nodes; ++n) {
    int i = n / NI, j = n % NI;
    double alpha = 0.0;
    for (int t = 0; t < nt; ++t) {
      alpha += O[t * NO + i] * I[t * NI + j];
    }
    alpha_max = fmax(
      alpha_max,
      up->lax_flux_nodal[dir](i, j, jacob_pos_l, jacob_pos_r, alpha, G_l, G_r, Fhat_nodal)
    );
  }

  // Projection of the nodal flux onto the surface modal basis.
  for (int tid = 0; tid < num_nodes; ++tid) {
    up->lax_prj[dir](tid, Fhat_nodal, conf_flux_surf);
  }
  double cflrate = up->lax_cfl[dir](dxv, jacob_pos_l, jacob_pos_r, alpha_max);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
  if (fabs(f_l[0]) < up->skip_cell_thresh && fabs(f_r[0]) < up->skip_cell_thresh) {
    return 0.0;
  }
  return cflrate;
}

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Serendipity basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_node, lax_flux_nodal_x_1x1v_ser_p2_node, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_node, lax_flux_nodal_x_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_node, lax_flux_nodal_x_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_node, lax_flux_nodal_x_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_node, lax_flux_nodal_x_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_x_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_prj, lax_flux_nodal_x_1x1v_ser_p2_prj, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_prj, lax_flux_nodal_x_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_prj, lax_flux_nodal_x_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_prj, lax_flux_nodal_x_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_prj, lax_flux_nodal_x_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_x_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_g, lax_flux_nodal_x_1x1v_ser_p2_g, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_g, lax_flux_nodal_x_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_g, lax_flux_nodal_x_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_g, lax_flux_nodal_x_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_g, lax_flux_nodal_x_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_cfl, lax_flux_nodal_x_1x1v_ser_p2_cfl, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_cfl, lax_flux_nodal_x_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_cfl, lax_flux_nodal_x_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_cfl, lax_flux_nodal_x_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_cfl, lax_flux_nodal_x_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_node, lax_flux_nodal_y_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_node, lax_flux_nodal_y_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_y_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_prj, lax_flux_nodal_y_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_prj, lax_flux_nodal_y_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_y_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_g, lax_flux_nodal_y_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_g, lax_flux_nodal_y_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_cfl, lax_flux_nodal_y_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_cfl, lax_flux_nodal_y_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_z_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_z_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Tensor basis).
// Phase-space Hamiltonian configuration-flux producers (Tensor basis):
// only the p=1 tensor hybrid has a phase-space Hamiltonian representation.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_x_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
  {NULL, hamil_phase_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, hamil_phase_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Tensor p=1 hybrid triad (vel_sparse) conf-flux producers; debug coverage 1x3v/2x2v.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_ho_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_ho_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_sparse_ho_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Tensor p=1 hybrid, dense velocity-space Hamiltonian conf-direction producers.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_ho_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_ho_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_vel_dense_ho_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_tensor_p1_node, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_tensor_p1_node, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_tensor_p1_node, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_tensor_p1_node, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_x_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_tensor_p1_prj, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_tensor_p1_prj, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_tensor_p1_prj, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_tensor_p1_prj, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_x_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_tensor_p1_g, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_tensor_p1_g, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_tensor_p1_g, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_tensor_p1_g, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_tensor_p1_cfl, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_tensor_p1_cfl, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_tensor_p1_cfl, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_tensor_p1_cfl, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_tensor_p1_node, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_y_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_tensor_p1_prj, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_y_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_tensor_p1_g, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_tensor_p1_cfl, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_z_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_z_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Serendipity basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_node, ho_lax_flux_nodal_x_1x1v_ser_p2_node, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_node, ho_lax_flux_nodal_x_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_node, ho_lax_flux_nodal_x_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_node, ho_lax_flux_nodal_x_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_node, ho_lax_flux_nodal_x_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_x_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_prj, ho_lax_flux_nodal_x_1x1v_ser_p2_prj, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_prj, ho_lax_flux_nodal_x_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_prj, ho_lax_flux_nodal_x_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_prj, ho_lax_flux_nodal_x_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_prj, ho_lax_flux_nodal_x_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_x_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_g, ho_lax_flux_nodal_x_1x1v_ser_p2_g, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_g, ho_lax_flux_nodal_x_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_g, ho_lax_flux_nodal_x_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_g, ho_lax_flux_nodal_x_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_g, ho_lax_flux_nodal_x_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_x_1x1v_ser_p1_cfl, ho_lax_flux_nodal_x_1x1v_ser_p2_cfl, NULL}, // 0
  {NULL, lax_flux_nodal_x_1x2v_ser_p1_cfl, ho_lax_flux_nodal_x_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_x_1x3v_ser_p1_cfl, ho_lax_flux_nodal_x_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_x_2x2v_ser_p1_cfl, ho_lax_flux_nodal_x_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_x_2x3v_ser_p1_cfl, ho_lax_flux_nodal_x_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_x_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_node, ho_lax_flux_nodal_y_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_node, ho_lax_flux_nodal_y_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_y_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_prj, ho_lax_flux_nodal_y_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_prj, ho_lax_flux_nodal_y_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_y_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_g, ho_lax_flux_nodal_y_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_g, ho_lax_flux_nodal_y_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_y_2x2v_ser_p1_cfl, ho_lax_flux_nodal_y_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_y_2x3v_ser_p1_cfl, ho_lax_flux_nodal_y_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_y_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_z_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_z_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_z_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Tensor basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_node, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_node, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_node, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_node, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_x_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_x_prj_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_prj, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_prj, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_prj, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_prj, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_x_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_x_g_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_g, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_g, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_g, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_g, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_x_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_cfl, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_cfl, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_cfl, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_cfl, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_x_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_node, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_y_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_y_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_prj, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_y_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_y_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_g, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_y_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_cfl, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_y_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_z_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_z_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_x_kernels[] = {
  // 1x kernels
  {NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_shared, NULL}, // 0
  {NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_x_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_x_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_x_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_x_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_x_kernels[] =
  {
    // 1x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared,
     hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_x_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_x_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_x_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_y_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1_shared,
   hamil_vel_dense_alpha_quad_y_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_y_kernels[] =
  {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_y_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1_shared,
     hamil_vel_sparse_alpha_quad_y_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_z_kernels[] =
  {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_dense_ho_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_sparse_ho_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_x_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_x_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_x_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_x_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_dense_ho_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_y_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1_shared,
     hamil_vel_dense_ho_alpha_quad_y_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_sparse_ho_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_y_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1_shared,
     hamil_vel_sparse_ho_alpha_quad_y_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_dense_ho_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_hamil_vel_sparse_ho_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, NULL, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, NULL, NULL, NULL}, // 5
    // 3x kernels
    {NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_x_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1_shared, hamil_phase_alpha_quad_x_1x1v_ser_p2_shared,
   NULL}, // 0
  {NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1_shared, hamil_phase_alpha_quad_x_1x2v_ser_p2_shared,
   NULL}, // 1
  {NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1_shared, hamil_phase_alpha_quad_x_1x3v_ser_p2_shared,
   NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1_shared, hamil_phase_alpha_quad_x_2x2v_ser_p2_shared,
   NULL}, // 4
  {NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1_shared, hamil_phase_alpha_quad_y_2x2v_ser_p2_shared,
   NULL}, // 4
  {NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_x_1x1v_ser_p2_shared, NULL}, // 0
  {NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_x_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_x_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_x_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_y_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Whole-surface (array-ABI) wrapper tables for the CPU dispatch: same
// coverage as the per-node tables above, entries are the original
// full-surface wrappers (table name gains _arr, entries drop _shared).
