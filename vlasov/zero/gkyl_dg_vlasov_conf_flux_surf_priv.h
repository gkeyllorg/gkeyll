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

typedef double (*hamil_alpha_quad_conf_t)(int i, int m, int hamil_pt_edge, const double *w, const double *dxv, 
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil); 

typedef double (*lax_flux_nodal_t)(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf); 

typedef double (*lax_cfl_t)(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max); 


// Whole-surface (array-ABI) kernel types: the original wrappers with
// hard-coded loop bounds INSIDE the kernel body, composed by the CPU
// dispatch (the per-node types above are the GPU 2D-launcher decomposition).
typedef void (*hamil_alpha_quad_conf_arr_t)(const double *w, const double *dxv, const int hamil_pt_edge,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad);

typedef double (*lax_flux_arr_t)(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal);

typedef double (*conf_flux_surf_t)(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  int dir, const double *w, const double *dxv, const int hamil_pt_edge,
  const double *vmap, const double *jacob_pos_l, const double *jacob_pos_r, const double *jacob_vel_surf,
  const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf); 

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { hamil_alpha_quad_conf_t kernels[4]; } gkyl_hamil_alpha_quad_kern_list;     
typedef struct { lax_flux_nodal_t kernels[4]; } gkyl_lax_flux_nodal_kern_list;  
typedef struct { lax_cfl_t kernels[4]; } gkyl_lax_cfl_kern_list;  
typedef struct { hamil_alpha_quad_conf_arr_t kernels[4]; } gkyl_hamil_alpha_quad_arr_kern_list;  
typedef struct { lax_flux_arr_t kernels[4]; } gkyl_lax_flux_arr_kern_list;  

struct gkyl_dg_vlasov_conf_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid. 
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index. 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either Velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in Velocity-space Jacobian. 
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array *vmap; // Velocity map (4-slot cubic rep per direction; borrowed from vel_map).
  const struct gkyl_array *jacob_vel_surf; // Velocity-space Jacobian at surface quadrature points (borrowed from vel_map).
  const struct gkyl_array *jacob_pos; // Configuration-space (position-map) Jacobian, per-conf-cell constant (borrowed from pos_map; defined on the extended conf range).
  hamil_alpha_quad_conf_t hamil_alpha_quad[3]; // Hamiltonian contribution to alpha_c at quadrature points. 
  lax_flux_nodal_t lax_flux_nodal[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion. 
  lax_cfl_t lax_cfl[3]; // CFL finalizer for the node-loop's alpha_max reduction. 
  hamil_alpha_quad_conf_arr_t hamil_alpha_quad_arr[3]; // Whole-surface (hard-coded-loop) variants for the CPU dispatch. 
  lax_flux_arr_t lax_flux_arr[3]; // Whole-surface Lax flux; returns the surface CFL estimate. 
  int num_nodes_conf; // Remaining configuration-space surface nodes per configuration surface. 
  int num_nodes_vel; // Velocity-space volume nodes per configuration surface. 
  conf_flux_surf_t conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes. 

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_vlasov_conf_flux_surf *on_dev; // pointer to itself or device data.  
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static double 
no_hamil_alpha_quad(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil)
{
  return 0.0;
}
GKYL_CU_DH
static void 
no_hamil_alpha_quad_arr(const double *w, const double *dxv, const int hamil_pt_edge,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad)
{
}

// Configuration-space flux assembly (Design A per-node dispatch): the
// Hamiltonian producer returns alpha at each surface node (register), the
// per-node Lax flux consumes it, and the CFL is finalized from the
// node-loop's alpha_max reduction.
GKYL_CU_DH
static double
conf_flux_surf_nodes(struct gkyl_dg_vlasov_conf_flux_surf *up,
  int dir, const double *w, const double *dxv, const int hamil_pt_edge,
  const double *vmap, const double *jacob_pos_l, const double *jacob_pos_r, const double *jacob_vel_surf,
  const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf)
{
  double alpha_max = 0.0;
  for (int i = 0; i < up->num_nodes_conf; ++i) {
    for (int m = 0; m < up->num_nodes_vel; ++m) {
      double alpha = up->hamil_alpha_quad[dir](i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos_r, jacob_vel_surf, poisson_tensor_conf, hamil);
      alpha_max = fmax(alpha_max, up->lax_flux_nodal[dir](i, m, jacob_pos_l, jacob_pos_r, alpha, f_l, f_r, conf_flux_surf));
    }
  }
  double cflrate = up->lax_cfl[dir](dxv, jacob_pos_l, jacob_pos_r, alpha_max);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_r[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  return cflrate; 
}


// Configuration-space flux assembly (whole-surface CPU dispatch): the
// Hamiltonian producer accumulates alpha at all surface nodes via the
// original wrapper (hard-coded loop bounds inside the kernel), then the
// nodal Lax flux wrapper consumes the buffer and returns the surface CFL.
GKYL_CU_DH
static double
conf_flux_surf_arrays(struct gkyl_dg_vlasov_conf_flux_surf *up,
  int dir, const double *w, const double *dxv, const int hamil_pt_edge,
  const double *vmap, const double *jacob_pos_l, const double *jacob_pos_r, const double *jacob_vel_surf,
  const double *poisson_tensor_conf,
  const double *hamil, const double *f_l, const double *f_r, double* GKYL_RESTRICT conf_flux_surf)
{
  double alpha_quad[256];
  const int nn = up->num_nodes_conf*up->num_nodes_vel;
  for (int n = 0; n < nn; ++n) alpha_quad[n] = 0.0;
  up->hamil_alpha_quad_arr[dir](w, dxv, hamil_pt_edge, vmap, jacob_pos_r, jacob_vel_surf, poisson_tensor_conf, hamil, alpha_quad);
  double cflrate = up->lax_flux_arr[dir](dxv, jacob_pos_l, jacob_pos_r, alpha_quad, f_l, f_r, conf_flux_surf);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_r[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  return cflrate; 
}

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1_node, lax_flux_nodal_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1_node, lax_flux_nodal_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1_node, lax_flux_nodal_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1_node, lax_flux_nodal_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1_node, lax_flux_nodal_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1_cfl, lax_flux_nodal_x_1x1v_ser_p2_cfl, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1_cfl, lax_flux_nodal_x_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1_cfl, lax_flux_nodal_x_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1_cfl, lax_flux_nodal_x_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1_cfl, lax_flux_nodal_x_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1_node, lax_flux_nodal_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1_node, lax_flux_nodal_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1_cfl, lax_flux_nodal_y_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1_cfl, lax_flux_nodal_y_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Tensor basis). 
// Phase-space Hamiltonian configuration-flux producers (Tensor basis):
// only the p=1 tensor hybrid has a phase-space Hamiltonian representation.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_ho_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_phase_ho_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_phase_ho_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_phase_ho_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(hamil_phase_ho_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Tensor p=1 hybrid triad (vel_sparse) conf-flux producers; debug coverage 1x3v/2x2v.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_ho_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_ho_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_sparse_ho_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_sparse_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_sparse_ho_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Tensor p=1 hybrid, dense velocity-space Hamiltonian conf-direction producers.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_vel_dense_ho_alpha_quad_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_ho_alpha_quad_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_ho_alpha_quad_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_ho_alpha_quad_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(hamil_vel_dense_ho_alpha_quad_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_vel_dense_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_vel_dense_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(hamil_vel_dense_ho_alpha_quad_z_3x3v_tensor_p1), NULL, NULL }, // 6
};


GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_tensor_p1_cfl, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_tensor_p1_cfl, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_tensor_p1_cfl, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_tensor_p1_cfl, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_x_2x3v_tensor_p1_cfl), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_tensor_p1_cfl, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_y_2x3v_tensor_p1_cfl), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1_node, ho_lax_flux_nodal_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1_node, ho_lax_flux_nodal_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1_node, ho_lax_flux_nodal_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1_node, ho_lax_flux_nodal_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1_node, ho_lax_flux_nodal_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1_cfl, ho_lax_flux_nodal_x_1x1v_ser_p2_cfl, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1_cfl, ho_lax_flux_nodal_x_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1_cfl, ho_lax_flux_nodal_x_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1_cfl, ho_lax_flux_nodal_x_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1_cfl, ho_lax_flux_nodal_x_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1_node, ho_lax_flux_nodal_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1_node, ho_lax_flux_nodal_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1_cfl, ho_lax_flux_nodal_y_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1_cfl, ho_lax_flux_nodal_y_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal Configuration-space flux conversion (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_x_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_x_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_x_cfl_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1_cfl, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1_cfl, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1_cfl, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1_cfl, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_x_2x3v_tensor_p1_cfl), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_y_2x3v_tensor_p1_node), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_y_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1_cfl, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_y_2x3v_tensor_p1_cfl), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_node), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_z_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1_cfl), NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node, hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1_node, hamil_vel_dense_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1_node, hamil_vel_dense_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1_node, hamil_vel_dense_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1_node, hamil_vel_dense_alpha_quad_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node, hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_node, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1_node, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1_node, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1_node, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1_node, hamil_vel_dense_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1_node, hamil_vel_dense_alpha_quad_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1_node, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1_node, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_x_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1_node, hamil_vel_dense_ho_alpha_quad_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1_node, hamil_vel_sparse_ho_alpha_quad_y_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_dense_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_vel_sparse_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1_node, hamil_phase_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1_node, hamil_phase_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1_node, hamil_phase_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1_node, hamil_phase_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1_node, hamil_phase_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_c evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_x_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1_node, hamil_phase_ho_alpha_quad_x_1x1v_ser_p2_node, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1_node, hamil_phase_ho_alpha_quad_x_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1_node, hamil_phase_ho_alpha_quad_x_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1_node, hamil_phase_ho_alpha_quad_x_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1_node, hamil_phase_ho_alpha_quad_y_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// Whole-surface (array-ABI) wrapper tables for the CPU dispatch: same
// coverage as the per-node tables above, entries are the original
// full-surface wrappers (table name gains _arr, entries drop _node).
GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_x_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1, lax_flux_nodal_x_1x1v_ser_p2, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1, lax_flux_nodal_x_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1, lax_flux_nodal_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1, lax_flux_nodal_x_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1, lax_flux_nodal_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1, lax_flux_nodal_y_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1, lax_flux_nodal_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_x_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(lax_flux_nodal_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(lax_flux_nodal_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_x_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_x_1x1v_ser_p1, ho_lax_flux_nodal_x_1x1v_ser_p2, NULL }, // 0
  { NULL, lax_flux_nodal_x_1x2v_ser_p1, ho_lax_flux_nodal_x_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_x_1x3v_ser_p1, ho_lax_flux_nodal_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_x_2x2v_ser_p1, ho_lax_flux_nodal_x_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_x_2x3v_ser_p1, ho_lax_flux_nodal_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_y_2x2v_ser_p1, ho_lax_flux_nodal_y_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_y_2x3v_ser_p1, ho_lax_flux_nodal_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_x_arr_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_x_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_x_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_x_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_x_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_x_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_x_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_y_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(ho_lax_flux_nodal_y_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_y_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(ho_lax_flux_nodal_z_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1, hamil_vel_dense_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1, hamil_vel_dense_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1, hamil_vel_dense_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1, hamil_vel_dense_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1, hamil_vel_dense_alpha_quad_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1, hamil_vel_dense_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1, hamil_vel_dense_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1, hamil_vel_dense_alpha_quad_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_vel_dense_alpha_quad_x_1x2v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_vel_dense_alpha_quad_x_1x3v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_x_2x2v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_x_2x3v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_vel_dense_alpha_quad_x_1x1v_ser_p1, hamil_vel_dense_ho_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1, hamil_vel_sparse_ho_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_vel_sparse_alpha_quad_x_1x3v_ser_p1, hamil_vel_sparse_ho_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_x_2x2v_ser_p1, hamil_vel_sparse_ho_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_x_2x3v_ser_p1, hamil_vel_sparse_ho_alpha_quad_x_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_dense_alpha_quad_y_2x2v_ser_p1, hamil_vel_dense_ho_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_dense_alpha_quad_y_2x3v_ser_p1, hamil_vel_dense_ho_alpha_quad_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_vel_sparse_alpha_quad_y_2x2v_ser_p1, hamil_vel_sparse_ho_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_vel_sparse_alpha_quad_y_2x3v_ser_p1, hamil_vel_sparse_ho_alpha_quad_y_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_dense_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_dense_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_vel_sparse_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_vel_sparse_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1, hamil_phase_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1, hamil_phase_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1, hamil_phase_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1, hamil_phase_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1, hamil_phase_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_x_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_x_1x1v_ser_p1, hamil_phase_ho_alpha_quad_x_1x1v_ser_p2, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_x_1x2v_ser_p1, hamil_phase_ho_alpha_quad_x_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_x_1x3v_ser_p1, hamil_phase_ho_alpha_quad_x_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_x_2x2v_ser_p1, hamil_phase_ho_alpha_quad_x_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_x_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_x_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_y_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_y_2x2v_ser_p1, hamil_phase_ho_alpha_quad_y_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_y_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_y_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_z_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_z_3x3v_ser_p1, NULL, NULL }, // 6
};
