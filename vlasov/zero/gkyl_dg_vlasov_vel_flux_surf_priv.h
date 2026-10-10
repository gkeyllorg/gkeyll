// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_flux_vlasov_kernels.h>
#include <gkyl_vlasov_surf_node_kernels.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>
#include <gkyl_vlasov_hyb_build.h>

// Force producers in sum-factorized form. The force at surface node (outer
// configuration node i, inner velocity node j) is
//   alpha(i, j) = sum_t O[t*NO + i] * I[t*NI + j]
// over the producer's terms t. A producer's _shared is called once per
// surface-node thread tid: it fills the outer factors O[(off + t)*NO + tid] if
// tid < NO and the inner factors I[(off + t)*NI + tid] if tid < NI (max(NO, NI)
// <= NO*NI, so one thread per surface node covers every producer), off being
// where its terms start, and returns its number of terms; the dispatch chains
// the producers, nt += producer(tid, nt, ..., O, I), and forms alpha(i, j) once
// over all terms. Called
// with O == NULL a producer only returns its term count (used to size the buffers).
typedef int (*hamil_alpha_shared_t)(
  int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos,
  const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
);

typedef int (*E_alpha_shared_t)(
  int tid, int off, const double *dxv, const double *qmem, double *GKYL_RESTRICT O,
  double *GKYL_RESTRICT I
);

typedef int (*phi_alpha_shared_t)(
  int tid, int off, const double *dxv, const double *jacob_pos, const double *phi_tot,
  double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
);

typedef int (*B_alpha_shared_t)(
  int tid, int off, const double *dxv, const double *jacob_vel, const double *hamil,
  const double *qmem, double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
);

typedef int (*rad_alpha_shared_t)(
  int tid, int off, const double *dxv, const double *rad, double *GKYL_RESTRICT O,
  double *GKYL_RESTRICT I
);

// Nodal Lax-Friedrichs flux in three stages. Stages 1 and 3 are cooperative:
// node thread tid does item tid (the launcher gives one thread per surface
// node, the CPU dispatch loops over them; items beyond a stage's count are
// no-ops, and there are never more items than surface nodes). Stage 1 of the
// sum-factorized nodal evaluation of f: item = (outer shape, inner node) pair,
// fills G_l[j*NA + a],
// G_r[j*NA + a] from f_l, f_c.
typedef void (*lax_g_t)(
  int tid, const double *f_l, const double *f_c, double *GKYL_RESTRICT G_l,
  double *GKYL_RESTRICT G_r
);

// Stage 2: the Lax flux at node (i, j) from the filled G arrays, written into
// the nodal buffer (num_nodes entries); returns |alpha| for the CFL reduction.
typedef double (*lax_flux_nodal_t)(
  int i, int j, const double *jacob_vel, double alpha, const double *G_l, const double *G_r,
  double *GKYL_RESTRICT Fhat_nodal
);

// Stage 3: projection of the nodal buffer onto the surface modal basis, stored
// at the direction's modal offset of the flux array; item = one inner surface
// mode (or one surface mode, by the kernel's choice).
typedef void (*lax_prj_t)(int tid, const double *Fhat_nodal, double *GKYL_RESTRICT flux);

// CFL estimate of the surface from the reduced alpha_max.
typedef double (*lax_cfl_t)(
  const double *dxv, const double *jacob_vel_l, const double *jacob_vel, double alpha_max
);

typedef double (*vel_flux_surf_t)(
  struct gkyl_dg_vlasov_vel_flux_surf *up, int dir, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_l, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil, const double *qmem, const double *phi_tot,
  const double *rad, const double *f_l, const double *f_c, double *GKYL_RESTRICT vel_flux_surf
);

typedef double (*vel_flux_surf_edge_t)(
  struct gkyl_dg_vlasov_vel_flux_surf *up, int dir, const double *w, const double *dxv,
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *qmem, const double *phi_tot, const double *rad, const double *f_c,
  double *GKYL_RESTRICT vel_flux_surf
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
  hamil_alpha_shared_t kernels[4];
} gkyl_hamil_alpha_quad_kern_list;
typedef struct {
  E_alpha_shared_t kernels[4];
} gkyl_E_alpha_quad_kern_list;
typedef struct {
  phi_alpha_shared_t kernels[4];
} gkyl_phi_alpha_quad_kern_list;
typedef struct {
  B_alpha_shared_t kernels[4];
} gkyl_B_alpha_quad_kern_list;
typedef struct {
  rad_alpha_shared_t kernels[4];
} gkyl_rad_alpha_quad_kern_list;
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
// * (num_nodes_conf + num_nodes_vel) (3x3v tensor p=1 triad: 56 terms x 24 nodes).
#define GKYL_VLASOV_VEL_FLUX_SURF_MAX_ALPHA_FACTORS 1344
// Largest number of surface nodes of the per-node dispatch (2x3v p2 higher-order: 4 x 64).
#define GKYL_VLASOV_FLUX_SURF_MAX_NODES 256
// Node threads per cell that do the first level of the CUDA kernels' alpha_max
// (CFL) reduction; the last node thread folds their partials.
#define GKYL_FLUX_SURF_CFL_REDUCERS 8

struct gkyl_dg_vlasov_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double
    skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  int hamil_dim; // Dimensionality of Hamiltonian.
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index.
  struct gkyl_range
    hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in velocity-space Jacobian.
  const struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array *
    jacob_vel_surf; // Velocity-space Jacobian at surface quadrature points (borrowed from vel_map; host pointer).
  const struct gkyl_array
    *vmap; // Velocity map (4-slot cubic rep per direction; borrowed from vel_map).
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array *
    jacob_pos; // Configuration-space (position-map) Jacobian (borrowed from pos_map; per-conf-cell constant).
  hamil_alpha_shared_t hamil_alpha_shared[3]; // Hamiltonian force -dH/dx (phase-space Hamiltonians).
  E_alpha_shared_t E_alpha_shared[3]; // Lorentz force from the electric field.
  phi_alpha_shared_t phi_alpha_shared[3]; // Scalar potential force -grad(phi).
  B_alpha_shared_t B_alpha_shared[3]; // Lorentz force from the magnetic field.
  rad_alpha_shared_t rad_alpha_shared[3]; // Radiation drag force.
  lax_g_t lax_g[3]; // Stage 1 of the nodal f evaluation (fills G_l, G_r).
  lax_flux_nodal_t lax_flux_nodal[3]; // Stage 2: Lax-Friedrichs flux at one surface node.
  lax_prj_t lax_prj[3]; // Stage 3: nodal -> surface-modal projection.
  lax_cfl_t lax_cfl[3]; // CFL estimate from the reduced alpha_max.
  int num_nodes_conf; // Configuration-space surface nodes per velocity surface.
  int num_nodes_vel; // Transverse velocity-space surface nodes per velocity surface.
  int num_surf_basis; // Surface modal basis functions per direction in the stored flux.
  int alpha_nterms_max; // Largest total force term count over the directions (buffer sizing).
  vel_flux_surf_t
    vel_flux_surf; // Assembly function for computing modal surface expansion of velocity-space fluxes.
  vel_flux_surf_edge_t
    vel_flux_surf_edge; // Edge-of-velocity-space Assembly function for computing velocity-space fluxes.

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_vlasov_vel_flux_surf *on_dev; // pointer to itself or device data.
};

static inline int
vel_flux_surf_num_surf_basis_ipow(int b, int e)
{
  int r = 1;
  for (int i = 0; i < e; ++i) {
    r *= b;
  }
  return r;
}

// Number of surface modal basis functions on a velocity-direction surface:
// the (pdim-1)-dimensional basis of the phase-space kernel type, i.e. the
// serendipity or tensor basis of the same order, or for the tensor p=1
// hybrid (p=1 configuration x p=2 velocity) 2^cdim * 3^(vdim-1). Matches the
// .nk of the kernels' vst_*_prj_* tables.
static inline int
vel_flux_surf_num_surf_basis(enum gkyl_basis_type kernel_b_type, int cdim, int vdim, int poly_order)
{
  struct gkyl_basis sb;
  if (kernel_b_type == GKYL_BASIS_MODAL_TENSOR) {
    if (poly_order == 1) {
      return vel_flux_surf_num_surf_basis_ipow(2, cdim) *
             vel_flux_surf_num_surf_basis_ipow(3, vdim - 1);
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
vel_flux_surf_num_nodes(
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
  for (int d = 0; d < cdim; ++d) {
    *num_nodes_conf *= nq_conf;
  }
  *num_nodes_vel = 1;
  for (int d = 0; d < vdim - 1; ++d) {
    *num_nodes_vel *= nq_vel;
  }
}

// Total number of force terms in direction dir (size query of the producers).
GKYL_CU_DH static inline int
vel_flux_surf_alpha_nterms(const struct gkyl_dg_vlasov_vel_flux_surf *up, int dir)
{
  return up->hamil_alpha_shared[dir](0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) +
         up->E_alpha_shared[dir](0, 0, 0, 0, 0, 0) +
         up->phi_alpha_shared[dir](0, 0, 0, 0, 0, 0, 0) +
         up->B_alpha_shared[dir](0, 0, 0, 0, 0, 0, 0, 0) +
         up->rad_alpha_shared[dir](0, 0, 0, 0, 0, 0);
}

// Empty function pointers for cases where these forces do not exist (no terms).
GKYL_CU_DH static int
no_hamil_alpha_shared(
  int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos,
  const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
)
{
  return 0;
}
GKYL_CU_DH static int
no_E_alpha_shared(
  int tid, int off, const double *dxv, const double *qmem, double *GKYL_RESTRICT O,
  double *GKYL_RESTRICT I
)
{
  return 0;
}
GKYL_CU_DH static int
no_phi_alpha_shared(
  int tid, int off, const double *dxv, const double *jacob_pos, const double *phi,
  double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
)
{
  return 0;
}
GKYL_CU_DH static int
no_B_alpha_shared(
  int tid, int off, const double *dxv, const double *jacob_vel, const double *hamil,
  const double *qmem, double *GKYL_RESTRICT O, double *GKYL_RESTRICT I
)
{
  return 0;
}
GKYL_CU_DH static int
no_rad_alpha_shared(
  int tid, int off, const double *dxv, const double *rad, double *GKYL_RESTRICT O,
  double *GKYL_RESTRICT I
)
{
  return 0;
}
GKYL_CU_DH static double
no_vel_flux_surf_edge(
  struct gkyl_dg_vlasov_vel_flux_surf *up, int dir, const double *w, const double *dxv,
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *qmem, const double *phi_tot, const double *rad, const double *f_c,
  double *GKYL_RESTRICT vel_flux_surf
)
{
  // Zero-flux boundary conditions, so immediately return.
  return 0.0;
}

// Velocity-space flux of one cell surface (CPU dispatch). This is the serial
// form of the algorithm the CUDA kernel runs with one thread per surface node:
// shared work of the cell first (the force factors and stage 1 of the nodal f
// evaluation), then the per-node work (force by the factor dot product, nodal
// Lax flux), then the projection onto the surface modal basis. Returns the
// surface CFL estimate.
static double
vel_flux_surf_arrays(
  struct gkyl_dg_vlasov_vel_flux_surf *up, int dir, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_l, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil, const double *qmem, const double *pot_tot,
  const double *rad, const double *f_l, const double *f_c, double *GKYL_RESTRICT vel_flux_surf
)
{
  const int NO = up->num_nodes_conf, NI = up->num_nodes_vel, num_nodes = NO * NI;
  double alpha_factors[GKYL_VLASOV_VEL_FLUX_SURF_MAX_ALPHA_FACTORS];
  double *O = alpha_factors, *I = alpha_factors + up->alpha_nterms_max * NO;
  // Sized for the largest surface node count (2x3v p2: 4 x 64 nodes); GKYL_DEFAULT_NUM_THREADS
  // is 1 on non-CUDA builds and must not size these buffers.
  double G_l[GKYL_VLASOV_FLUX_SURF_MAX_NODES], G_c[GKYL_VLASOV_FLUX_SURF_MAX_NODES];
  double Fhat_nodal[GKYL_VLASOV_FLUX_SURF_MAX_NODES];

  // Shared work of the cell, as each surface-node thread of the CUDA kernel
  // does it: its row of the outer/inner factors of every force producer (each
  // producer's terms follow the previous one's) and its item of stage 1 of the
  // nodal f evaluation.
  int nt = 0;
  for (int tid = 0; tid < num_nodes; ++tid) {
    nt = 0;
    nt += up->hamil_alpha_shared[dir](
      tid, nt, w, dxv, vmap, jacob_pos, jacob_vel, poisson_tensor_conf, hamil, O, I
    );
    nt += up->E_alpha_shared[dir](tid, nt, dxv, qmem, O, I);
    nt += up->phi_alpha_shared[dir](tid, nt, dxv, jacob_pos, pot_tot, O, I);
    nt += up->B_alpha_shared[dir](tid, nt, dxv, jacob_vel, hamil, qmem, O, I);
    nt += up->rad_alpha_shared[dir](tid, nt, dxv, rad, O, I);
    up->lax_g[dir](tid, f_l, f_c, G_l, G_c);
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
    alpha_max =
      fmax(alpha_max, up->lax_flux_nodal[dir](i, j, jacob_vel, alpha, G_l, G_c, Fhat_nodal));
  }

  // Projection of the nodal flux onto the surface modal basis.
  for (int tid = 0; tid < num_nodes; ++tid) {
    up->lax_prj[dir](tid, Fhat_nodal, vel_flux_surf);
  }
  double cflrate = up->lax_cfl[dir](dxv, jacob_vel_l, jacob_vel, alpha_max);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate.
  if (fabs(f_l[0]) < up->skip_cell_thresh && fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0;
  }
  return cflrate;
}

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_node, lax_flux_nodal_vx_1x1v_ser_p2_node,
   lax_flux_nodal_vx_1x1v_ser_p3_node}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_node, lax_flux_nodal_vx_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_node, lax_flux_nodal_vx_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_node, lax_flux_nodal_vx_2x1v_ser_p2_node,
   lax_flux_nodal_vx_2x1v_ser_p3_node}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_node, lax_flux_nodal_vx_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_node, lax_flux_nodal_vx_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_vx_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_prj, lax_flux_nodal_vx_1x1v_ser_p2_prj,
   lax_flux_nodal_vx_1x1v_ser_p3_prj}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_prj, lax_flux_nodal_vx_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_prj, lax_flux_nodal_vx_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_prj, lax_flux_nodal_vx_2x1v_ser_p2_prj,
   lax_flux_nodal_vx_2x1v_ser_p3_prj}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_prj, lax_flux_nodal_vx_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_prj, lax_flux_nodal_vx_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_vx_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_g, lax_flux_nodal_vx_1x1v_ser_p2_g,
   lax_flux_nodal_vx_1x1v_ser_p3_g}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_g, lax_flux_nodal_vx_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_g, lax_flux_nodal_vx_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_g, lax_flux_nodal_vx_2x1v_ser_p2_g,
   lax_flux_nodal_vx_2x1v_ser_p3_g}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_g, lax_flux_nodal_vx_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_g, lax_flux_nodal_vx_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_cfl, lax_flux_nodal_vx_1x1v_ser_p2_cfl,
   lax_flux_nodal_vx_1x1v_ser_p3_cfl}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_cfl, lax_flux_nodal_vx_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_cfl, lax_flux_nodal_vx_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_cfl, lax_flux_nodal_vx_2x1v_ser_p2_cfl,
   lax_flux_nodal_vx_2x1v_ser_p3_cfl}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_cfl, lax_flux_nodal_vx_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_cfl, lax_flux_nodal_vx_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_node, lax_flux_nodal_vy_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_node, lax_flux_nodal_vy_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_node, lax_flux_nodal_vy_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_node, lax_flux_nodal_vy_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_vy_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_prj, lax_flux_nodal_vy_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_prj, lax_flux_nodal_vy_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_prj, lax_flux_nodal_vy_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_prj, lax_flux_nodal_vy_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_vy_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_g, lax_flux_nodal_vy_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_g, lax_flux_nodal_vy_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_g, lax_flux_nodal_vy_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_g, lax_flux_nodal_vy_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_cfl, lax_flux_nodal_vy_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_cfl, lax_flux_nodal_vy_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_cfl, lax_flux_nodal_vy_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_cfl, lax_flux_nodal_vy_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_node, lax_flux_nodal_vz_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_node, lax_flux_nodal_vz_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_lax_flux_nodal_vz_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_prj, lax_flux_nodal_vz_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_prj, lax_flux_nodal_vz_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_lax_flux_nodal_vz_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_g, lax_flux_nodal_vz_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_g, lax_flux_nodal_vz_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_cfl, lax_flux_nodal_vz_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_cfl, lax_flux_nodal_vz_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Tensor basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_tensor_p1_node, lax_flux_nodal_vx_1x1v_tensor_p2_node,
   lax_flux_nodal_vx_1x1v_tensor_p3_node}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_tensor_p1_node, lax_flux_nodal_vx_1x2v_tensor_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_tensor_p1_node, lax_flux_nodal_vx_1x3v_tensor_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_tensor_p1_node, lax_flux_nodal_vx_2x1v_tensor_p2_node,
   lax_flux_nodal_vx_2x1v_tensor_p3_node}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_tensor_p1_node, lax_flux_nodal_vx_2x2v_tensor_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vx_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_vx_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_tensor_p1_prj, lax_flux_nodal_vx_1x1v_tensor_p2_prj,
   lax_flux_nodal_vx_1x1v_tensor_p3_prj}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_tensor_p1_prj, lax_flux_nodal_vx_1x2v_tensor_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_tensor_p1_prj, lax_flux_nodal_vx_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_tensor_p1_prj, lax_flux_nodal_vx_2x1v_tensor_p2_prj,
   lax_flux_nodal_vx_2x1v_tensor_p3_prj}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_tensor_p1_prj, lax_flux_nodal_vx_2x2v_tensor_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vx_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_vx_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_tensor_p1_g, lax_flux_nodal_vx_1x1v_tensor_p2_g,
   lax_flux_nodal_vx_1x1v_tensor_p3_g}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_tensor_p1_g, lax_flux_nodal_vx_1x2v_tensor_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_tensor_p1_g, lax_flux_nodal_vx_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_tensor_p1_g, lax_flux_nodal_vx_2x1v_tensor_p2_g,
   lax_flux_nodal_vx_2x1v_tensor_p3_g}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_tensor_p1_g, lax_flux_nodal_vx_2x2v_tensor_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vx_3x3v_tensor_p1_g), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_node, lax_flux_nodal_vx_1x1v_tensor_p2_node,
   lax_flux_nodal_vx_1x1v_tensor_p3_node}, // 0
  {NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_node, lax_flux_nodal_vx_1x2v_tensor_p2_node, NULL
  }, // 1
  {NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_node, lax_flux_nodal_vx_1x3v_tensor_p2_node, NULL
  }, // 2
  // 2x kernels
  {NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_node, lax_flux_nodal_vx_2x1v_tensor_p2_node,
   lax_flux_nodal_vx_2x1v_tensor_p3_node}, // 3
  {NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_node, lax_flux_nodal_vx_2x2v_tensor_p2_node, NULL
  }, // 4
  {NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vx_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_vx_prj_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_prj, lax_flux_nodal_vx_1x1v_tensor_p2_prj,
   lax_flux_nodal_vx_1x1v_tensor_p3_prj}, // 0
  {NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_prj, lax_flux_nodal_vx_1x2v_tensor_p2_prj, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_prj, lax_flux_nodal_vx_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_prj, lax_flux_nodal_vx_2x1v_tensor_p2_prj,
   lax_flux_nodal_vx_2x1v_tensor_p3_prj}, // 3
  {NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_prj, lax_flux_nodal_vx_2x2v_tensor_p2_prj, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vx_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_vx_g_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_g, lax_flux_nodal_vx_1x1v_tensor_p2_g,
   lax_flux_nodal_vx_1x1v_tensor_p3_g}, // 0
  {NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_g, lax_flux_nodal_vx_1x2v_tensor_p2_g, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_g, lax_flux_nodal_vx_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_g, lax_flux_nodal_vx_2x1v_tensor_p2_g,
   lax_flux_nodal_vx_2x1v_tensor_p3_g}, // 3
  {NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_g, lax_flux_nodal_vx_2x2v_tensor_p2_g, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vx_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_tensor_p1_cfl, lax_flux_nodal_vx_1x1v_tensor_p2_cfl,
   lax_flux_nodal_vx_1x1v_tensor_p3_cfl}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_tensor_p1_cfl, lax_flux_nodal_vx_1x2v_tensor_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_tensor_p1_cfl, lax_flux_nodal_vx_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_tensor_p1_cfl, lax_flux_nodal_vx_2x1v_tensor_p2_cfl,
   lax_flux_nodal_vx_2x1v_tensor_p3_cfl}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_tensor_p1_cfl, lax_flux_nodal_vx_2x2v_tensor_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vx_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  {NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_cfl, lax_flux_nodal_vx_1x1v_tensor_p2_cfl,
   lax_flux_nodal_vx_1x1v_tensor_p3_cfl}, // 0
  {NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_cfl, lax_flux_nodal_vx_1x2v_tensor_p2_cfl, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_cfl, lax_flux_nodal_vx_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_cfl, lax_flux_nodal_vx_2x1v_tensor_p2_cfl,
   lax_flux_nodal_vx_2x1v_tensor_p3_cfl}, // 3
  {NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_cfl, lax_flux_nodal_vx_2x2v_tensor_p2_cfl, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vx_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_tensor_p1_node, lax_flux_nodal_vy_1x2v_tensor_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_tensor_p1_node, lax_flux_nodal_vy_1x3v_tensor_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_tensor_p1_node, lax_flux_nodal_vy_2x2v_tensor_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vy_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_vy_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_tensor_p1_prj, lax_flux_nodal_vy_1x2v_tensor_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_tensor_p1_prj, lax_flux_nodal_vy_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_tensor_p1_prj, lax_flux_nodal_vy_2x2v_tensor_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vy_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_vy_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_tensor_p1_g, lax_flux_nodal_vy_1x2v_tensor_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_tensor_p1_g, lax_flux_nodal_vy_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_tensor_p1_g, lax_flux_nodal_vy_2x2v_tensor_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vy_3x3v_tensor_p1_g), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_node, lax_flux_nodal_vy_1x2v_tensor_p2_node, NULL
  }, // 1
  {NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_node, lax_flux_nodal_vy_1x3v_tensor_p2_node, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_node, lax_flux_nodal_vy_2x2v_tensor_p2_node, NULL
  }, // 4
  {NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vy_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_vy_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_prj, lax_flux_nodal_vy_1x2v_tensor_p2_prj, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_prj, lax_flux_nodal_vy_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_prj, lax_flux_nodal_vy_2x2v_tensor_p2_prj, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vy_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_vy_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_g, lax_flux_nodal_vy_1x2v_tensor_p2_g, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_g, lax_flux_nodal_vy_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_g, lax_flux_nodal_vy_2x2v_tensor_p2_g, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vy_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_tensor_p1_cfl, lax_flux_nodal_vy_1x2v_tensor_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_tensor_p1_cfl, lax_flux_nodal_vy_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_tensor_p1_cfl, lax_flux_nodal_vy_2x2v_tensor_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vy_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_cfl, lax_flux_nodal_vy_1x2v_tensor_p2_cfl, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_cfl, lax_flux_nodal_vy_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_cfl, lax_flux_nodal_vy_2x2v_tensor_p2_cfl, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vy_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_tensor_p1_node, lax_flux_nodal_vz_1x3v_tensor_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vz_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_lax_flux_nodal_vz_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_tensor_p1_prj, lax_flux_nodal_vz_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vz_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_lax_flux_nodal_vz_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_tensor_p1_g, lax_flux_nodal_vz_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vz_3x3v_tensor_p1_g), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_node, lax_flux_nodal_vz_1x3v_tensor_p2_node, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_node, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vz_3x3v_tensor_p1_node), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list tensor_ho_lax_flux_nodal_vz_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_prj, lax_flux_nodal_vz_1x3v_tensor_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_prj, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vz_3x3v_tensor_p1_prj), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list tensor_ho_lax_flux_nodal_vz_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_g, lax_flux_nodal_vz_1x3v_tensor_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_g, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vz_3x3v_tensor_p1_g), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_tensor_p1_cfl, lax_flux_nodal_vz_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(lax_flux_nodal_vz_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_cfl, lax_flux_nodal_vz_1x3v_tensor_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_cfl, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_vz_3x3v_tensor_p1_cfl), NULL, NULL} // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis).
GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_node, ho_lax_flux_nodal_vx_1x1v_ser_p2_node,
   ho_lax_flux_nodal_vx_1x1v_ser_p3_node}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_node, ho_lax_flux_nodal_vx_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_node, ho_lax_flux_nodal_vx_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_node, ho_lax_flux_nodal_vx_2x1v_ser_p2_node,
   ho_lax_flux_nodal_vx_2x1v_ser_p3_node}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_node, ho_lax_flux_nodal_vx_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_node, ho_lax_flux_nodal_vx_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_vx_prj_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_prj, ho_lax_flux_nodal_vx_1x1v_ser_p2_prj,
   ho_lax_flux_nodal_vx_1x1v_ser_p3_prj}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_prj, ho_lax_flux_nodal_vx_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_prj, ho_lax_flux_nodal_vx_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_prj, ho_lax_flux_nodal_vx_2x1v_ser_p2_prj,
   ho_lax_flux_nodal_vx_2x1v_ser_p3_prj}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_prj, ho_lax_flux_nodal_vx_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_prj, ho_lax_flux_nodal_vx_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_vx_g_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_g, ho_lax_flux_nodal_vx_1x1v_ser_p2_g,
   ho_lax_flux_nodal_vx_1x1v_ser_p3_g}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_g, ho_lax_flux_nodal_vx_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_g, ho_lax_flux_nodal_vx_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_g, ho_lax_flux_nodal_vx_2x1v_ser_p2_g,
   ho_lax_flux_nodal_vx_2x1v_ser_p3_g}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_g, ho_lax_flux_nodal_vx_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_g, ho_lax_flux_nodal_vx_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_nodal_vx_1x1v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x1v_ser_p2_cfl,
   ho_lax_flux_nodal_vx_1x1v_ser_p3_cfl}, // 0
  {NULL, lax_flux_nodal_vx_1x2v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vx_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, lax_flux_nodal_vx_2x1v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x1v_ser_p2_cfl,
   ho_lax_flux_nodal_vx_2x1v_ser_p3_cfl}, // 3
  {NULL, lax_flux_nodal_vx_2x2v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vx_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vx_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_node, ho_lax_flux_nodal_vy_1x2v_ser_p2_node, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_node, ho_lax_flux_nodal_vy_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_node, ho_lax_flux_nodal_vy_2x2v_ser_p2_node, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_node, ho_lax_flux_nodal_vy_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_vy_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_prj, ho_lax_flux_nodal_vy_1x2v_ser_p2_prj, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_prj, ho_lax_flux_nodal_vy_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_prj, ho_lax_flux_nodal_vy_2x2v_ser_p2_prj, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_prj, ho_lax_flux_nodal_vy_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_vy_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_g, ho_lax_flux_nodal_vy_1x2v_ser_p2_g, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_g, ho_lax_flux_nodal_vy_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_g, ho_lax_flux_nodal_vy_2x2v_ser_p2_g, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_g, ho_lax_flux_nodal_vy_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, lax_flux_nodal_vy_1x2v_ser_p1_cfl, ho_lax_flux_nodal_vy_1x2v_ser_p2_cfl, NULL}, // 1
  {NULL, lax_flux_nodal_vy_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vy_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, lax_flux_nodal_vy_2x2v_ser_p1_cfl, ho_lax_flux_nodal_vy_2x2v_ser_p2_cfl, NULL}, // 4
  {NULL, lax_flux_nodal_vy_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vy_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vy_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_node, ho_lax_flux_nodal_vz_1x3v_ser_p2_node, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_node, ho_lax_flux_nodal_vz_2x3v_ser_p2_node, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_node, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_prj_kern_list ser_ho_lax_flux_nodal_vz_prj_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_prj, ho_lax_flux_nodal_vz_1x3v_ser_p2_prj, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_prj, ho_lax_flux_nodal_vz_2x3v_ser_p2_prj, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_prj, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_g_kern_list ser_ho_lax_flux_nodal_vz_g_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_g, ho_lax_flux_nodal_vz_1x3v_ser_p2_g, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_g, ho_lax_flux_nodal_vz_2x3v_ser_p2_g, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_g, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, lax_flux_nodal_vz_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vz_1x3v_ser_p2_cfl, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, lax_flux_nodal_vz_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vz_2x3v_ser_p2_cfl, NULL}, // 5
  // 3x kernels
  {NULL, lax_flux_nodal_vz_3x3v_ser_p1_cfl, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1_shared, hamil_phase_alpha_quad_vx_1x1v_ser_p2_shared,
   hamil_phase_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared, hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared,
   NULL}, // 1
  {NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared, hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared,
   NULL}, // 2
  // 2x kernels
  {NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1_shared, hamil_phase_alpha_quad_vx_2x1v_ser_p2_shared,
   hamil_phase_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared, hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared,
   NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 1
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared, hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared,
   NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2_shared,
   hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2_shared,
   hamil_phase_ho_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 1
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared,
   hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_hamil_alpha_shared, no_hamil_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for Hamil (NC) - vel dependance only - Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for Hamil (NC) - vel dependance only - Hamiltonian forces (Serendipity basis).
// high-order intergation used (_ho)
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_shared,
     nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_shared,
     nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
  {NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1_shared,
   nc_hamil_phase_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis).
// high-order intergation used (_ho)
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_phase_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2_shared, NULL}, // 0
    {NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_phase_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  ser_nc_hamil_phase_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1_shared,
     nc_hamil_phase_ho_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Serendipity basis).
// Phase-space Hamiltonian velocity-flux producers (Tensor basis): only the
// p=1 tensor hybrid has a phase-space Hamiltonian representation; plain
// lists are the lo variant, *_ho_* lists the ho variant.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, hamil_phase_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
  {NULL, hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, hamil_phase_alpha_quad_vx_2x1v_tensor_p1_shared, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_hamil_alpha_shared, NULL, NULL}, // 1
  {NULL, no_hamil_alpha_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, hamil_phase_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_hamil_alpha_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_hamil_alpha_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, hamil_phase_ho_alpha_quad_vx_2x1v_tensor_p1_shared, NULL, NULL}, // 3
    {NULL, hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, no_hamil_alpha_shared, NULL, NULL}, // 1
    {NULL, no_hamil_alpha_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_hamil_phase_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, no_hamil_alpha_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, no_hamil_alpha_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, NULL, NULL}, // 0
  {NULL, B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, NULL, NULL}, // 3
  {NULL, B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, NULL, NULL}, // 0
  {NULL, B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, B_ho_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, NULL, NULL}, // 3
  {NULL, B_ho_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, B_ho_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_ho_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
  {NULL, B_ho_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_ho_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
  {NULL, B_ho_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V_PHASE(B_ho_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, E_alpha_quad_vx_1x1v_ser_p1_shared, E_alpha_quad_vx_1x1v_ser_p2_shared,
   E_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, E_alpha_quad_vx_1x2v_ser_p1_shared, E_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vx_1x3v_ser_p1_shared, E_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, E_alpha_quad_vx_2x1v_ser_p1_shared, E_alpha_quad_vx_2x1v_ser_p2_shared,
   E_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, E_alpha_quad_vx_2x2v_ser_p1_shared, E_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vx_2x3v_ser_p1_shared, E_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, E_alpha_quad_vy_1x2v_ser_p1_shared, E_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vy_1x3v_ser_p1_shared, E_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, E_alpha_quad_vy_2x2v_ser_p1_shared, E_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vy_2x3v_ser_p1_shared, E_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, E_alpha_quad_vz_1x3v_ser_p1_shared, E_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, E_alpha_quad_vz_2x3v_ser_p1_shared, E_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Serendipity basis).
GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, E_alpha_quad_vx_1x1v_ser_p1_shared, E_ho_alpha_quad_vx_1x1v_ser_p2_shared,
   E_ho_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, E_alpha_quad_vx_1x2v_ser_p1_shared, E_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vx_1x3v_ser_p1_shared, E_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, E_alpha_quad_vx_2x1v_ser_p1_shared, E_ho_alpha_quad_vx_2x1v_ser_p2_shared,
   E_ho_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, E_alpha_quad_vx_2x2v_ser_p1_shared, E_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vx_2x3v_ser_p1_shared, E_ho_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, E_alpha_quad_vy_1x2v_ser_p1_shared, E_ho_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vy_1x3v_ser_p1_shared, E_ho_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, E_alpha_quad_vy_2x2v_ser_p1_shared, E_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vy_2x3v_ser_p1_shared, E_ho_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, E_alpha_quad_vz_1x3v_ser_p1_shared, E_ho_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, E_alpha_quad_vz_2x3v_ser_p1_shared, E_ho_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, E_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Tensor basis).
GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, E_alpha_quad_vx_1x1v_tensor_p1_shared, E_alpha_quad_vx_1x1v_tensor_p2_shared,
   E_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, E_alpha_quad_vx_1x2v_tensor_p1_shared, E_alpha_quad_vx_1x2v_tensor_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vx_1x3v_tensor_p1_shared, E_alpha_quad_vx_1x3v_tensor_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, E_alpha_quad_vx_2x1v_tensor_p1_shared, E_alpha_quad_vx_2x1v_tensor_p2_shared,
   E_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, E_alpha_quad_vx_2x2v_tensor_p1_shared, E_alpha_quad_vx_2x2v_tensor_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, E_ho_alpha_quad_vx_1x1v_tensor_p1_shared, E_alpha_quad_vx_1x1v_tensor_p2_shared,
   E_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, E_ho_alpha_quad_vx_1x2v_tensor_p1_shared, E_alpha_quad_vx_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, E_ho_alpha_quad_vx_1x3v_tensor_p1_shared, E_alpha_quad_vx_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, E_ho_alpha_quad_vx_2x1v_tensor_p1_shared, E_alpha_quad_vx_2x1v_tensor_p2_shared,
   E_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, E_ho_alpha_quad_vx_2x2v_tensor_p1_shared, E_alpha_quad_vx_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, E_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, E_alpha_quad_vy_1x2v_tensor_p1_shared, E_alpha_quad_vy_1x2v_tensor_p2_shared, NULL}, // 1
  {NULL, E_alpha_quad_vy_1x3v_tensor_p1_shared, E_alpha_quad_vy_1x3v_tensor_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, E_alpha_quad_vy_2x2v_tensor_p1_shared, E_alpha_quad_vy_2x2v_tensor_p2_shared, NULL}, // 4
  {NULL, E_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, E_ho_alpha_quad_vy_1x2v_tensor_p1_shared, E_alpha_quad_vy_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, E_ho_alpha_quad_vy_1x3v_tensor_p1_shared, E_alpha_quad_vy_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, E_ho_alpha_quad_vy_2x2v_tensor_p1_shared, E_alpha_quad_vy_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, E_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, E_alpha_quad_vz_1x3v_tensor_p1_shared, E_alpha_quad_vz_1x3v_tensor_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, E_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, E_ho_alpha_quad_vz_1x3v_tensor_p1_shared, E_alpha_quad_vz_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, E_ho_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(E_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Serendipity basis).
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, phi_alpha_quad_vx_1x1v_ser_p1_shared, phi_alpha_quad_vx_1x1v_ser_p2_shared,
   phi_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, phi_alpha_quad_vx_1x2v_ser_p1_shared, phi_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, phi_alpha_quad_vx_1x3v_ser_p1_shared, phi_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, phi_alpha_quad_vx_2x1v_ser_p1_shared, phi_alpha_quad_vx_2x1v_ser_p2_shared,
   phi_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, phi_alpha_quad_vx_2x2v_ser_p1_shared, phi_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, phi_alpha_quad_vx_2x3v_ser_p1_shared, phi_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, phi_alpha_quad_vy_2x2v_ser_p1_shared, phi_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, phi_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Serendipity basis).
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, phi_alpha_quad_vx_1x1v_ser_p1_shared, phi_ho_alpha_quad_vx_1x1v_ser_p2_shared,
   phi_ho_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, phi_alpha_quad_vx_1x2v_ser_p1_shared, phi_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, phi_alpha_quad_vx_1x3v_ser_p1_shared, phi_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, phi_alpha_quad_vx_2x1v_ser_p1_shared, phi_ho_alpha_quad_vx_2x1v_ser_p2_shared,
   phi_ho_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, phi_alpha_quad_vx_2x2v_ser_p1_shared, phi_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, phi_alpha_quad_vx_2x3v_ser_p1_shared, phi_ho_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, phi_alpha_quad_vy_2x2v_ser_p1_shared, phi_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, phi_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, phi_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Tensor basis).
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, phi_alpha_quad_vx_1x1v_tensor_p1_shared, phi_alpha_quad_vx_1x1v_tensor_p2_shared,
   phi_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, phi_alpha_quad_vx_1x2v_tensor_p1_shared, phi_alpha_quad_vx_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, phi_alpha_quad_vx_1x3v_tensor_p1_shared, phi_alpha_quad_vx_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, phi_alpha_quad_vx_2x1v_tensor_p1_shared, phi_alpha_quad_vx_2x1v_tensor_p2_shared,
   phi_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, phi_alpha_quad_vx_2x2v_tensor_p1_shared, phi_alpha_quad_vx_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, phi_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, phi_ho_alpha_quad_vx_1x1v_tensor_p1_shared, phi_alpha_quad_vx_1x1v_tensor_p2_shared,
   phi_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, phi_ho_alpha_quad_vx_1x2v_tensor_p1_shared, phi_alpha_quad_vx_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, phi_ho_alpha_quad_vx_1x3v_tensor_p1_shared, phi_alpha_quad_vx_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, phi_ho_alpha_quad_vx_2x1v_tensor_p1_shared, phi_alpha_quad_vx_2x1v_tensor_p2_shared,
   phi_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, phi_ho_alpha_quad_vx_2x2v_tensor_p1_shared, phi_alpha_quad_vx_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, phi_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, phi_alpha_quad_vy_2x2v_tensor_p1_shared, phi_alpha_quad_vy_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, phi_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, phi_ho_alpha_quad_vy_2x2v_tensor_p1_shared, phi_alpha_quad_vy_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, phi_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, no_phi_alpha_shared, no_phi_alpha_shared, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(phi_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Serendipity basis).
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
  {NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
  {NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
  {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
  {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_shared,
   B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_shared,
   B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Serendipity basis).
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_shared,
     B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  ser_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_shared,
     B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
    // 3x kernels
    {NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis).
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
  {NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
  {NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_shared,
   B_hamil_phase_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis).
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
  {NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
  {NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_shared,
   B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Tensor basis).
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_dense_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 0
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, no_B_alpha_shared, no_B_alpha_shared, no_B_alpha_shared}, // 3
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_dense_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_shared, NULL}, // 1
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_shared, NULL}, // 4
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_dense_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_shared,
     B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Sparse-Hamiltonian variant.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_B_alpha_quad_kern_list
  tensor_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_shared,
     B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2_shared, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(B_ho_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Serendipity basis).
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, rad_alpha_quad_vx_1x1v_ser_p1_shared, rad_alpha_quad_vx_1x1v_ser_p2_shared,
   rad_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, rad_alpha_quad_vx_1x2v_ser_p1_shared, rad_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, rad_alpha_quad_vx_1x3v_ser_p1_shared, rad_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, rad_alpha_quad_vx_2x1v_ser_p1_shared, rad_alpha_quad_vx_2x1v_ser_p2_shared,
   rad_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, rad_alpha_quad_vx_2x2v_ser_p1_shared, rad_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, rad_alpha_quad_vx_2x3v_ser_p1_shared, rad_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, rad_alpha_quad_vy_1x2v_ser_p1_shared, rad_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, rad_alpha_quad_vy_1x3v_ser_p1_shared, rad_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, rad_alpha_quad_vy_2x2v_ser_p1_shared, rad_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, rad_alpha_quad_vy_2x3v_ser_p1_shared, rad_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, rad_alpha_quad_vz_1x3v_ser_p1_shared, rad_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, rad_alpha_quad_vz_2x3v_ser_p1_shared, rad_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Serendipity basis).
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, rad_alpha_quad_vx_1x1v_ser_p1_shared, rad_ho_alpha_quad_vx_1x1v_ser_p2_shared,
   rad_ho_alpha_quad_vx_1x1v_ser_p3_shared}, // 0
  {NULL, rad_alpha_quad_vx_1x2v_ser_p1_shared, rad_ho_alpha_quad_vx_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, rad_alpha_quad_vx_1x3v_ser_p1_shared, rad_ho_alpha_quad_vx_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, rad_alpha_quad_vx_2x1v_ser_p1_shared, rad_ho_alpha_quad_vx_2x1v_ser_p2_shared,
   rad_ho_alpha_quad_vx_2x1v_ser_p3_shared}, // 3
  {NULL, rad_alpha_quad_vx_2x2v_ser_p1_shared, rad_ho_alpha_quad_vx_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, rad_alpha_quad_vx_2x3v_ser_p1_shared, rad_ho_alpha_quad_vx_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vx_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, rad_alpha_quad_vy_1x2v_ser_p1_shared, rad_ho_alpha_quad_vy_1x2v_ser_p2_shared, NULL}, // 1
  {NULL, rad_alpha_quad_vy_1x3v_ser_p1_shared, rad_ho_alpha_quad_vy_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, rad_alpha_quad_vy_2x2v_ser_p1_shared, rad_ho_alpha_quad_vy_2x2v_ser_p2_shared, NULL}, // 4
  {NULL, rad_alpha_quad_vy_2x3v_ser_p1_shared, rad_ho_alpha_quad_vy_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vy_3x3v_ser_p1_shared, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, rad_alpha_quad_vz_1x3v_ser_p1_shared, rad_ho_alpha_quad_vz_1x3v_ser_p2_shared, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, rad_alpha_quad_vz_2x3v_ser_p1_shared, rad_ho_alpha_quad_vz_2x3v_ser_p2_shared, NULL}, // 5
  // 3x kernels
  {NULL, rad_alpha_quad_vz_3x3v_ser_p1_shared, NULL, NULL} // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Tensor basis).
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, rad_alpha_quad_vx_1x1v_tensor_p1_shared, rad_alpha_quad_vx_1x1v_tensor_p2_shared,
   rad_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, rad_alpha_quad_vx_1x2v_tensor_p1_shared, rad_alpha_quad_vx_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, rad_alpha_quad_vx_1x3v_tensor_p1_shared, rad_alpha_quad_vx_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, rad_alpha_quad_vx_2x1v_tensor_p1_shared, rad_alpha_quad_vx_2x1v_tensor_p2_shared,
   rad_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, rad_alpha_quad_vx_2x2v_tensor_p1_shared, rad_alpha_quad_vx_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, rad_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  {NULL, rad_ho_alpha_quad_vx_1x1v_tensor_p1_shared, rad_alpha_quad_vx_1x1v_tensor_p2_shared,
   rad_alpha_quad_vx_1x1v_tensor_p3_shared}, // 0
  {NULL, rad_ho_alpha_quad_vx_1x2v_tensor_p1_shared, rad_alpha_quad_vx_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, rad_ho_alpha_quad_vx_1x3v_tensor_p1_shared, rad_alpha_quad_vx_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, rad_ho_alpha_quad_vx_2x1v_tensor_p1_shared, rad_alpha_quad_vx_2x1v_tensor_p2_shared,
   rad_alpha_quad_vx_2x1v_tensor_p3_shared}, // 3
  {NULL, rad_ho_alpha_quad_vx_2x2v_tensor_p1_shared, rad_alpha_quad_vx_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, rad_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, rad_alpha_quad_vy_1x2v_tensor_p1_shared, rad_alpha_quad_vy_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, rad_alpha_quad_vy_1x3v_tensor_p1_shared, rad_alpha_quad_vy_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, rad_alpha_quad_vy_2x2v_tensor_p1_shared, rad_alpha_quad_vy_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, rad_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, rad_ho_alpha_quad_vy_1x2v_tensor_p1_shared, rad_alpha_quad_vy_1x2v_tensor_p2_shared, NULL
  }, // 1
  {NULL, rad_ho_alpha_quad_vy_1x3v_tensor_p1_shared, rad_alpha_quad_vy_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, rad_ho_alpha_quad_vy_2x2v_tensor_p1_shared, rad_alpha_quad_vy_2x2v_tensor_p2_shared, NULL
  }, // 4
  {NULL, rad_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, rad_alpha_quad_vz_1x3v_tensor_p1_shared, rad_alpha_quad_vz_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, rad_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, rad_ho_alpha_quad_vz_1x3v_tensor_p1_shared, rad_alpha_quad_vz_1x3v_tensor_p2_shared, NULL
  }, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, rad_ho_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
  // 3x kernels
  {NULL, GKYL_HYB_3X3V(rad_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

// Whole-surface (array-ABI) wrapper tables for the CPU dispatch: same
// coverage as the per-node tables above, entries are the original
// full-surface wrappers (table name gains _arr, entries drop _shared).
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

// Tensor p=1 hybrid triad velocity-direction producers: dense velocity-space
// Hamiltonian and phase-space Hamiltonian (GR triads).
GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL
    }, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V(nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_phase_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_ho_alpha_quad_vx_kernels[] = {
    // 1x kernels
    {NULL, nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_shared, NULL, NULL}, // 0
    {NULL, nc_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_ho_alpha_quad_vy_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, nc_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1_shared, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, nc_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_shared, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

GKYL_CU_D static const gkyl_hamil_alpha_quad_kern_list
  tensor_nc_hamil_phase_ho_alpha_quad_vz_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    {NULL, NULL, NULL, NULL}, // 1
    {NULL, nc_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1_shared, NULL, NULL}, // 2
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 3
    {NULL, NULL, NULL, NULL}, // 4
    {NULL, nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1_shared, NULL, NULL}, // 5
    // 3x kernels
    {NULL, GKYL_HYB_3X3V_PHASE(nc_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1_shared), NULL, NULL
    } // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// Whole-surface (arr) variants: alpha_v at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis).
// Whole-surface (arr) variants: alpha_v at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis).
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
