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

typedef double (*hamil_alpha_quad_t)(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil); 

typedef double (*E_alpha_quad_t)(int i, int j, const double *dxv, const double *qmem);  

typedef double (*phi_alpha_quad_t)(int i, int j, const double *dxv, const double *jacob_pos,
  const double *phi_tot);

typedef double (*B_alpha_quad_t)(int i, int j, const double *dxv,
  const double *jacob_vel, const double *hamil, const double *qmem); 

typedef double (*rad_alpha_quad_t)(int i, int j, const double *dxv, const double *rad); 

typedef double (*lax_flux_nodal_t)(int i, int j, const double *jacob_vel, double alpha,
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf); 

typedef double (*lax_cfl_t)(const double *dxv,
  const double *jacob_vel_l, const double *jacob_vel, double alpha_max); 


// Whole-surface (array-ABI) kernel types: the original wrappers with
// hard-coded loop bounds INSIDE the kernel body. The CPU dispatch composes
// these so the compiler can inline, unroll, and hoist loop-invariant work;
// the per-node types above are the GPU 2D-launcher decomposition.
typedef void (*hamil_alpha_quad_arr_t)(const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad);

typedef void (*E_alpha_quad_arr_t)(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad);

typedef void (*phi_alpha_quad_arr_t)(const double *dxv, const double *jacob_pos,
  const double *phi, double* GKYL_RESTRICT alpha_quad);

typedef void (*B_alpha_quad_arr_t)(const double *dxv, const double *jacob_vel,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad);

typedef void (*rad_alpha_quad_arr_t)(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad);

typedef double (*lax_flux_arr_t)(const double *dxv, const double *jacob_vel_l, const double *jacob_vel,
  const double *alpha_quad, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal);

typedef double (*vel_flux_surf_t)(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_l, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf); 

typedef double (*vel_flux_surf_edge_t)(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *w, const double *dxv, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_c, double* GKYL_RESTRICT vel_flux_surf); 

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
typedef struct { hamil_alpha_quad_t kernels[4]; } gkyl_hamil_alpha_quad_kern_list;  
typedef struct { E_alpha_quad_t kernels[4]; } gkyl_E_alpha_quad_kern_list;  
typedef struct { phi_alpha_quad_t kernels[4]; } gkyl_phi_alpha_quad_kern_list;  
typedef struct { B_alpha_quad_t kernels[4]; } gkyl_B_alpha_quad_kern_list;  
typedef struct { rad_alpha_quad_t kernels[4]; } gkyl_rad_alpha_quad_kern_list;    
typedef struct { lax_flux_nodal_t kernels[4]; } gkyl_lax_flux_nodal_kern_list;  
typedef struct { lax_cfl_t kernels[4]; } gkyl_lax_cfl_kern_list;  
typedef struct { hamil_alpha_quad_arr_t kernels[4]; } gkyl_hamil_alpha_quad_arr_kern_list;  
typedef struct { E_alpha_quad_arr_t kernels[4]; } gkyl_E_alpha_quad_arr_kern_list;  
typedef struct { phi_alpha_quad_arr_t kernels[4]; } gkyl_phi_alpha_quad_arr_kern_list;  
typedef struct { B_alpha_quad_arr_t kernels[4]; } gkyl_B_alpha_quad_arr_kern_list;  
typedef struct { rad_alpha_quad_arr_t kernels[4]; } gkyl_rad_alpha_quad_arr_kern_list;  
typedef struct { lax_flux_arr_t kernels[4]; } gkyl_lax_flux_arr_kern_list;  

struct gkyl_dg_vlasov_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid. 
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index. 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  struct gkyl_range vel_range; // Velocity-space range for use in velocity-space Jacobian.
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array *jacob_vel_surf; // Velocity-space Jacobian at surface quadrature points (borrowed from vel_map; host pointer).
  const struct gkyl_array *vmap; // Velocity map (4-slot cubic rep per direction; borrowed from vel_map).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array *jacob_pos; // Configuration-space (position-map) Jacobian (borrowed from pos_map; per-conf-cell constant).
  hamil_alpha_quad_t hamil_alpha_quad[3]; // Hamiltonian contribution to alpha_v at quadrature points.
  E_alpha_quad_t E_alpha_quad[3]; // Lorentz force contribution from electric field to alpha_v at quadrature points. 
  phi_alpha_quad_t phi_alpha_quad[3]; // Scalar potential, -grad(phi), force contribution to alpha_v at quadrature points. 
  B_alpha_quad_t B_alpha_quad[3]; // Lorentz force contribution from magnetic field to alpha_v at quadrature points. 
  rad_alpha_quad_t rad_alpha_quad[3]; // Radiation drag force contribution to alpha_v at quadrature points. 
  lax_flux_nodal_t lax_flux_nodal[3]; // Per-node Lax-Friedrichs flux at surface quadrature points. 
  lax_cfl_t lax_cfl[3]; // CFL finalizer for the node-loop's alpha_max reduction. 
  hamil_alpha_quad_arr_t hamil_alpha_quad_arr[3]; // Whole-surface (hard-coded-loop) variants for the CPU dispatch. 
  E_alpha_quad_arr_t E_alpha_quad_arr[3]; 
  phi_alpha_quad_arr_t phi_alpha_quad_arr[3]; 
  B_alpha_quad_arr_t B_alpha_quad_arr[3]; 
  rad_alpha_quad_arr_t rad_alpha_quad_arr[3]; 
  lax_flux_arr_t lax_flux_arr[3]; // Whole-surface Lax flux; returns the surface CFL estimate. 
  int num_nodes_conf; // Configuration-space surface nodes per velocity surface. 
  int num_nodes_vel; // Transverse velocity-space surface nodes per velocity surface. 
  vel_flux_surf_t vel_flux_surf; // Assembly function for computing modal surface expansion of velocity-space fluxes. 
  vel_flux_surf_edge_t vel_flux_surf_edge; // Edge-of-velocity-space Assembly function for computing velocity-space fluxes. 

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_vlasov_vel_flux_surf *on_dev; // pointer to itself or device data.  
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static double 
no_hamil_alpha_quad(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil)
{
  return 0.0;
}
GKYL_CU_DH
static double 
no_E_alpha_quad(int i, int j, const double *dxv, const double *qmem)
{
  return 0.0;
}
GKYL_CU_DH
static double 
no_phi_alpha_quad(int i, int j, const double *dxv, const double *jacob_pos,
  const double *phi)
{
  return 0.0;
}
GKYL_CU_DH
static double 
no_B_alpha_quad(int i, int j, const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem)
{
  return 0.0;
}
GKYL_CU_DH
static double 
no_rad_alpha_quad(int i, int j, const double *dxv, const double *rad)
{
  return 0.0;
}
GKYL_CU_DH
static void 
no_hamil_alpha_quad_arr(const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_E_alpha_quad_arr(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_phi_alpha_quad_arr(const double *dxv, const double *jacob_pos,
  const double *phi, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_B_alpha_quad_arr(const double *dxv, const double *jacob_vel, 
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static void 
no_rad_alpha_quad_arr(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad)
{
}
GKYL_CU_DH
static double
no_vel_flux_surf_edge(struct gkyl_dg_vlasov_vel_flux_surf *up, 
  int dir, const double *w, const double *dxv, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const double *qmem, const double *phi_tot, const double *rad, 
  const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  // Zero-flux boundary conditions, so immediately return. 
  return 0.0; 
}

// Velocity-space flux assembly (Design A per-node dispatch): compose the
// force producers into a register alpha at each surface node, apply the
// per-node Lax flux, reduce alpha_max, and finalize the CFL estimate.
GKYL_CU_DH
static double
vel_flux_surf_nodes(struct gkyl_dg_vlasov_vel_flux_surf *up,
  int dir, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_l, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad,
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  double alpha_max = 0.0;
  for (int i = 0; i < up->num_nodes_conf; ++i) {
    for (int j = 0; j < up->num_nodes_vel; ++j) {
      double alpha = up->hamil_alpha_quad[dir](i, j, w, dxv, vmap, jacob_pos, jacob_vel, poisson_tensor_conf, hamil)
        + up->E_alpha_quad[dir](i, j, dxv, qmem)
        + up->phi_alpha_quad[dir](i, j, dxv, jacob_pos, pot_tot)
        + up->B_alpha_quad[dir](i, j, dxv, jacob_vel, hamil, qmem)
        + up->rad_alpha_quad[dir](i, j, dxv, rad);
      alpha_max = fmax(alpha_max, up->lax_flux_nodal[dir](i, j, jacob_vel, alpha, f_l, f_c, vel_flux_surf));
    }
  }
  double cflrate = up->lax_cfl[dir](dxv, jacob_vel_l, jacob_vel, alpha_max);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  return cflrate; 
}


// Velocity-space flux assembly (whole-surface CPU dispatch): accumulate the
// force producers into the nodal alpha_quad buffer via the original wrappers
// (hard-coded loop bounds inside each kernel), then apply the nodal Lax flux
// wrapper, which returns the surface CFL estimate.
GKYL_CU_DH
static double
vel_flux_surf_arrays(struct gkyl_dg_vlasov_vel_flux_surf *up,
  int dir, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_l, const double *jacob_vel, const double *poisson_tensor_conf,
  const double *hamil, const double *qmem, const double *pot_tot, const double *rad,
  const double *f_l, const double *f_c, double* GKYL_RESTRICT vel_flux_surf)
{
  double alpha_quad[256];
  const int nn = up->num_nodes_conf*up->num_nodes_vel;
  for (int n = 0; n < nn; ++n) alpha_quad[n] = 0.0;
  up->hamil_alpha_quad_arr[dir](w, dxv, vmap, jacob_pos, jacob_vel, poisson_tensor_conf, hamil, alpha_quad);
  up->E_alpha_quad_arr[dir](dxv, qmem, alpha_quad);
  up->phi_alpha_quad_arr[dir](dxv, jacob_pos, pot_tot, alpha_quad);
  up->B_alpha_quad_arr[dir](dxv, jacob_vel, hamil, qmem, alpha_quad);
  up->rad_alpha_quad_arr[dir](dxv, rad, alpha_quad);
  double cflrate = up->lax_flux_arr[dir](dxv, jacob_vel_l, jacob_vel, alpha_quad, f_l, f_c, vel_flux_surf);

  // Always compute the flux, but if we are below threshold, ignore the stable time step estimate. 
  if (fabs(f_l[0]) < up->skip_cell_thresh && 
      fabs(f_c[0]) < up->skip_cell_thresh) {
    return 0.0; 
  }
  return cflrate; 
}

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1_node, lax_flux_nodal_vx_1x1v_ser_p2_node, lax_flux_nodal_vx_1x1v_ser_p3_node }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1_node, lax_flux_nodal_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1_node, lax_flux_nodal_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1_node, lax_flux_nodal_vx_2x1v_ser_p2_node, lax_flux_nodal_vx_2x1v_ser_p3_node }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1_node, lax_flux_nodal_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1_node, lax_flux_nodal_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1_cfl, lax_flux_nodal_vx_1x1v_ser_p2_cfl, lax_flux_nodal_vx_1x1v_ser_p3_cfl }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1_cfl, lax_flux_nodal_vx_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1_cfl, lax_flux_nodal_vx_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1_cfl, lax_flux_nodal_vx_2x1v_ser_p2_cfl, lax_flux_nodal_vx_2x1v_ser_p3_cfl }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1_cfl, lax_flux_nodal_vx_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1_cfl, lax_flux_nodal_vx_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1_node, lax_flux_nodal_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1_node, lax_flux_nodal_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1_node, lax_flux_nodal_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1_node, lax_flux_nodal_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1_cfl, lax_flux_nodal_vy_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1_cfl, lax_flux_nodal_vy_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1_cfl, lax_flux_nodal_vy_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1_cfl, lax_flux_nodal_vy_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1_node, lax_flux_nodal_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1_node, lax_flux_nodal_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1_cfl, lax_flux_nodal_vz_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1_cfl, lax_flux_nodal_vz_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Tensor basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_tensor_p1_node, lax_flux_nodal_vx_1x1v_tensor_p2_node, lax_flux_nodal_vx_1x1v_tensor_p3_node }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_tensor_p1_node, lax_flux_nodal_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_tensor_p1_node, lax_flux_nodal_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_tensor_p1_node, lax_flux_nodal_vx_2x1v_tensor_p2_node, lax_flux_nodal_vx_2x1v_tensor_p3_node }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_tensor_p1_node, lax_flux_nodal_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_node, lax_flux_nodal_vx_1x1v_tensor_p2_node, lax_flux_nodal_vx_1x1v_tensor_p3_node }, // 0
  { NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_node, lax_flux_nodal_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_node, lax_flux_nodal_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_node, lax_flux_nodal_vx_2x1v_tensor_p2_node, lax_flux_nodal_vx_2x1v_tensor_p3_node }, // 3
  { NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_node, lax_flux_nodal_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_tensor_p1_cfl, lax_flux_nodal_vx_1x1v_tensor_p2_cfl, lax_flux_nodal_vx_1x1v_tensor_p3_cfl }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_tensor_p1_cfl, lax_flux_nodal_vx_1x2v_tensor_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_tensor_p1_cfl, lax_flux_nodal_vx_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_tensor_p1_cfl, lax_flux_nodal_vx_2x1v_tensor_p2_cfl, lax_flux_nodal_vx_2x1v_tensor_p3_cfl }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_tensor_p1_cfl, lax_flux_nodal_vx_2x2v_tensor_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1_cfl, lax_flux_nodal_vx_1x1v_tensor_p2_cfl, lax_flux_nodal_vx_1x1v_tensor_p3_cfl }, // 0
  { NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1_cfl, lax_flux_nodal_vx_1x2v_tensor_p2_cfl, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1_cfl, lax_flux_nodal_vx_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1_cfl, lax_flux_nodal_vx_2x1v_tensor_p2_cfl, lax_flux_nodal_vx_2x1v_tensor_p3_cfl }, // 3
  { NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1_cfl, lax_flux_nodal_vx_2x2v_tensor_p2_cfl, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vx_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_tensor_p1_node, lax_flux_nodal_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_tensor_p1_node, lax_flux_nodal_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_tensor_p1_node, lax_flux_nodal_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_node, lax_flux_nodal_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_node, lax_flux_nodal_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_node, lax_flux_nodal_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_tensor_p1_cfl, lax_flux_nodal_vy_1x2v_tensor_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_tensor_p1_cfl, lax_flux_nodal_vy_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_tensor_p1_cfl, lax_flux_nodal_vy_2x2v_tensor_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1_cfl, lax_flux_nodal_vy_1x2v_tensor_p2_cfl, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1_cfl, lax_flux_nodal_vy_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1_cfl, lax_flux_nodal_vy_2x2v_tensor_p2_cfl, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vy_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_tensor_p1_node, lax_flux_nodal_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list tensor_ho_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_node, lax_flux_nodal_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_tensor_p1_cfl, lax_flux_nodal_vz_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_cfl_kern_list tensor_ho_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1_cfl, lax_flux_nodal_vz_1x3v_tensor_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1_cfl, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vz_3x3v_tensor_p1_cfl, NULL, NULL }, // 6
};

// Nodal Lax-Friedrichs to modal velocity-space flux conversion (Serendipity basis). 
GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vx_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1_node, ho_lax_flux_nodal_vx_1x1v_ser_p2_node, ho_lax_flux_nodal_vx_1x1v_ser_p3_node }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1_node, ho_lax_flux_nodal_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1_node, ho_lax_flux_nodal_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1_node, ho_lax_flux_nodal_vx_2x1v_ser_p2_node, ho_lax_flux_nodal_vx_2x1v_ser_p3_node }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1_node, ho_lax_flux_nodal_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1_node, ho_lax_flux_nodal_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vx_cfl_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x1v_ser_p2_cfl, ho_lax_flux_nodal_vx_1x1v_ser_p3_cfl }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vx_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x1v_ser_p2_cfl, ho_lax_flux_nodal_vx_2x1v_ser_p3_cfl }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vx_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1_node, ho_lax_flux_nodal_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1_node, ho_lax_flux_nodal_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1_node, ho_lax_flux_nodal_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1_node, ho_lax_flux_nodal_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vy_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1_cfl, ho_lax_flux_nodal_vy_1x2v_ser_p2_cfl, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vy_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1_cfl, ho_lax_flux_nodal_vy_2x2v_ser_p2_cfl, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vy_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_nodal_kern_list ser_ho_lax_flux_nodal_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1_node, ho_lax_flux_nodal_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1_node, ho_lax_flux_nodal_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_cfl_kern_list ser_ho_lax_flux_nodal_vz_cfl_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1_cfl, ho_lax_flux_nodal_vz_1x3v_ser_p2_cfl, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1_cfl, ho_lax_flux_nodal_vz_2x3v_ser_p2_cfl, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1_cfl, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1_node, hamil_phase_alpha_quad_vx_1x1v_ser_p2_node, hamil_phase_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, hamil_phase_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, hamil_phase_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1_node, hamil_phase_alpha_quad_vx_2x1v_ser_p2_node, hamil_phase_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, hamil_phase_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, hamil_phase_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for general Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1_node, hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2_node, hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1_node, hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2_node, hamil_phase_ho_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_hamil_phase_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, no_hamil_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};


// alpha_v evaluated at quadrature points for Hamil (NC) - vel dependance only - Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_node, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for Hamil (NC) - vel dependance only - Hamiltonian forces (Serendipity basis).
// high-order intergation used (_ho)
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node, nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_node, nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};





// alpha_v evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for general (NC) Hamiltonian forces (Serendipity basis). 
// high-order intergation used (_ho)
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2_node, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list ser_nc_hamil_phase_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node, nc_hamil_phase_ho_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Serendipity basis). 
// Phase-space Hamiltonian velocity-flux producers (Tensor basis): only the
// p=1 tensor hybrid has a phase-space Hamiltonian representation; plain
// lists are the lo variant, *_ho_* lists the ho variant.
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_tensor_p1_node, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_2x1v_tensor_p1_node, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_hamil_phase_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, NULL, NULL }, // 0
  { NULL, B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, NULL, NULL }, // 3
  { NULL, B_ho_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_ser_p1_node, E_alpha_quad_vx_1x1v_ser_p2_node, E_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, E_alpha_quad_vx_1x2v_ser_p1_node, E_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_ser_p1_node, E_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_ser_p1_node, E_alpha_quad_vx_2x1v_ser_p2_node, E_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, E_alpha_quad_vx_2x2v_ser_p1_node, E_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_ser_p1_node, E_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_ser_p1_node, E_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_ser_p1_node, E_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_ser_p1_node, E_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_ser_p1_node, E_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_ser_p1_node, E_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_ser_p1_node, E_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Serendipity basis). 
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_ser_p1_node, E_ho_alpha_quad_vx_1x1v_ser_p2_node, E_ho_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, E_alpha_quad_vx_1x2v_ser_p1_node, E_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_ser_p1_node, E_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_ser_p1_node, E_ho_alpha_quad_vx_2x1v_ser_p2_node, E_ho_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, E_alpha_quad_vx_2x2v_ser_p1_node, E_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_ser_p1_node, E_ho_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_ser_p1_node, E_ho_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_ser_p1_node, E_ho_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_ser_p1_node, E_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_ser_p1_node, E_ho_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list ser_E_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_ser_p1_node, E_ho_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_ser_p1_node, E_ho_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the electric field Lorentz force (Tensor basis). 
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_tensor_p1_node, E_alpha_quad_vx_1x1v_tensor_p2_node, E_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, E_alpha_quad_vx_1x2v_tensor_p1_node, E_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_tensor_p1_node, E_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_tensor_p1_node, E_alpha_quad_vx_2x1v_tensor_p2_node, E_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, E_alpha_quad_vx_2x2v_tensor_p1_node, E_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, E_ho_alpha_quad_vx_1x1v_tensor_p1_node, E_alpha_quad_vx_1x1v_tensor_p2_node, E_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, E_ho_alpha_quad_vx_1x2v_tensor_p1_node, E_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, E_ho_alpha_quad_vx_1x3v_tensor_p1_node, E_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, E_ho_alpha_quad_vx_2x1v_tensor_p1_node, E_alpha_quad_vx_2x1v_tensor_p2_node, E_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, E_ho_alpha_quad_vx_2x2v_tensor_p1_node, E_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, E_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_tensor_p1_node, E_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_tensor_p1_node, E_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_tensor_p1_node, E_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_ho_alpha_quad_vy_1x2v_tensor_p1_node, E_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, E_ho_alpha_quad_vy_1x3v_tensor_p1_node, E_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_ho_alpha_quad_vy_2x2v_tensor_p1_node, E_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, E_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_tensor_p1_node, E_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_kern_list tensor_E_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_ho_alpha_quad_vz_1x3v_tensor_p1_node, E_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_ho_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_ser_p1_node, phi_alpha_quad_vx_1x1v_ser_p2_node, phi_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_ser_p1_node, phi_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_ser_p1_node, phi_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_ser_p1_node, phi_alpha_quad_vx_2x1v_ser_p2_node, phi_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_ser_p1_node, phi_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_ser_p1_node, phi_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_ser_p1_node, phi_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Serendipity basis). 
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_ser_p1_node, phi_ho_alpha_quad_vx_1x1v_ser_p2_node, phi_ho_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_ser_p1_node, phi_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_ser_p1_node, phi_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_ser_p1_node, phi_ho_alpha_quad_vx_2x1v_ser_p2_node, phi_ho_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_ser_p1_node, phi_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_ser_p1_node, phi_ho_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_ser_p1_node, phi_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list ser_phi_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for scalar potential forces (Tensor basis). 
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_tensor_p1_node, phi_alpha_quad_vx_1x1v_tensor_p2_node, phi_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_tensor_p1_node, phi_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_tensor_p1_node, phi_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_tensor_p1_node, phi_alpha_quad_vx_2x1v_tensor_p2_node, phi_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_tensor_p1_node, phi_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, phi_ho_alpha_quad_vx_1x1v_tensor_p1_node, phi_alpha_quad_vx_1x1v_tensor_p2_node, phi_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, phi_ho_alpha_quad_vx_1x2v_tensor_p1_node, phi_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, phi_ho_alpha_quad_vx_1x3v_tensor_p1_node, phi_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, phi_ho_alpha_quad_vx_2x1v_tensor_p1_node, phi_alpha_quad_vx_2x1v_tensor_p2_node, phi_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, phi_ho_alpha_quad_vx_2x2v_tensor_p1_node, phi_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, phi_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_tensor_p1_node, phi_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_ho_alpha_quad_vy_2x2v_tensor_p1_node, phi_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, phi_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_kern_list tensor_phi_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad, no_phi_alpha_quad, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_node, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1_node, B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1_node, B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, B_hamil_phase_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, B_hamil_phase_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, B_hamil_phase_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1_node, B_hamil_phase_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1_node, B_hamil_phase_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, B_hamil_phase_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node, B_hamil_phase_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list ser_B_ho_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node, B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};


// alpha_v evaluated at quadrature points for the magnetic field Lorentz force (Tensor basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 0
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad, no_B_alpha_quad, no_B_alpha_quad }, // 3
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_node, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_node, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};


// alpha_v evaluated at quadrature points for the radiation drag force (Serendipity basis). 
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_ser_p1_node, rad_alpha_quad_vx_1x1v_ser_p2_node, rad_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_ser_p1_node, rad_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_ser_p1_node, rad_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_ser_p1_node, rad_alpha_quad_vx_2x1v_ser_p2_node, rad_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_ser_p1_node, rad_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_ser_p1_node, rad_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_ser_p1_node, rad_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_ser_p1_node, rad_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_ser_p1_node, rad_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_ser_p1_node, rad_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_ser_p1_node, rad_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_ser_p1_node, rad_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Serendipity basis). 
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_ser_p1_node, rad_ho_alpha_quad_vx_1x1v_ser_p2_node, rad_ho_alpha_quad_vx_1x1v_ser_p3_node }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_ser_p1_node, rad_ho_alpha_quad_vx_1x2v_ser_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_ser_p1_node, rad_ho_alpha_quad_vx_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_ser_p1_node, rad_ho_alpha_quad_vx_2x1v_ser_p2_node, rad_ho_alpha_quad_vx_2x1v_ser_p3_node }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_ser_p1_node, rad_ho_alpha_quad_vx_2x2v_ser_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_ser_p1_node, rad_ho_alpha_quad_vx_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_ser_p1_node, rad_ho_alpha_quad_vy_1x2v_ser_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_ser_p1_node, rad_ho_alpha_quad_vy_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_ser_p1_node, rad_ho_alpha_quad_vy_2x2v_ser_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_ser_p1_node, rad_ho_alpha_quad_vy_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_ser_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list ser_rad_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_ser_p1_node, rad_ho_alpha_quad_vz_1x3v_ser_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_ser_p1_node, rad_ho_alpha_quad_vz_2x3v_ser_p2_node, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_ser_p1_node, NULL, NULL }, // 6
};

// alpha_v evaluated at quadrature points for the radiation drag force (Tensor basis). 
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_tensor_p1_node, rad_alpha_quad_vx_1x1v_tensor_p2_node, rad_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_tensor_p1_node, rad_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_tensor_p1_node, rad_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_tensor_p1_node, rad_alpha_quad_vx_2x1v_tensor_p2_node, rad_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_tensor_p1_node, rad_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, rad_ho_alpha_quad_vx_1x1v_tensor_p1_node, rad_alpha_quad_vx_1x1v_tensor_p2_node, rad_alpha_quad_vx_1x1v_tensor_p3_node }, // 0
  { NULL, rad_ho_alpha_quad_vx_1x2v_tensor_p1_node, rad_alpha_quad_vx_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vx_1x3v_tensor_p1_node, rad_alpha_quad_vx_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, rad_ho_alpha_quad_vx_2x1v_tensor_p1_node, rad_alpha_quad_vx_2x1v_tensor_p2_node, rad_alpha_quad_vx_2x1v_tensor_p3_node }, // 3
  { NULL, rad_ho_alpha_quad_vx_2x2v_tensor_p1_node, rad_alpha_quad_vx_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_tensor_p1_node, rad_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_tensor_p1_node, rad_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_tensor_p1_node, rad_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_ho_alpha_quad_vy_1x2v_tensor_p1_node, rad_alpha_quad_vy_1x2v_tensor_p2_node, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vy_1x3v_tensor_p1_node, rad_alpha_quad_vy_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_ho_alpha_quad_vy_2x2v_tensor_p1_node, rad_alpha_quad_vy_2x2v_tensor_p2_node, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_tensor_p1_node, rad_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_kern_list tensor_rad_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vz_1x3v_tensor_p1_node, rad_alpha_quad_vz_1x3v_tensor_p2_node, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

// Whole-surface (array-ABI) wrapper tables for the CPU dispatch: same
// coverage as the per-node tables above, entries are the original
// full-surface wrappers (table name gains _arr, entries drop _node).
GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1, lax_flux_nodal_vx_1x1v_ser_p2, lax_flux_nodal_vx_1x1v_ser_p3 }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1, lax_flux_nodal_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1, lax_flux_nodal_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1, lax_flux_nodal_vx_2x1v_ser_p2, lax_flux_nodal_vx_2x1v_ser_p3 }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1, lax_flux_nodal_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1, lax_flux_nodal_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1, lax_flux_nodal_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1, lax_flux_nodal_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1, lax_flux_nodal_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1, lax_flux_nodal_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_lax_flux_nodal_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1, lax_flux_nodal_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1, lax_flux_nodal_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_tensor_p1, lax_flux_nodal_vx_1x1v_tensor_p2, lax_flux_nodal_vx_1x1v_tensor_p3 }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_tensor_p1, lax_flux_nodal_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_tensor_p1, lax_flux_nodal_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_tensor_p1, lax_flux_nodal_vx_2x1v_tensor_p2, lax_flux_nodal_vx_2x1v_tensor_p3 }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_tensor_p1, lax_flux_nodal_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, ho_lax_flux_nodal_vx_1x1v_tensor_p1, lax_flux_nodal_vx_1x1v_tensor_p2, lax_flux_nodal_vx_1x1v_tensor_p3 }, // 0
  { NULL, ho_lax_flux_nodal_vx_1x2v_tensor_p1, lax_flux_nodal_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vx_1x3v_tensor_p1, lax_flux_nodal_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, ho_lax_flux_nodal_vx_2x1v_tensor_p1, lax_flux_nodal_vx_2x1v_tensor_p2, lax_flux_nodal_vx_2x1v_tensor_p3 }, // 3
  { NULL, ho_lax_flux_nodal_vx_2x2v_tensor_p1, lax_flux_nodal_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_tensor_p1, lax_flux_nodal_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_tensor_p1, lax_flux_nodal_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_tensor_p1, lax_flux_nodal_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, ho_lax_flux_nodal_vy_1x2v_tensor_p1, lax_flux_nodal_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vy_1x3v_tensor_p1, lax_flux_nodal_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, ho_lax_flux_nodal_vy_2x2v_tensor_p1, lax_flux_nodal_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_lax_flux_nodal_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_tensor_p1, lax_flux_nodal_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list tensor_ho_lax_flux_nodal_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, ho_lax_flux_nodal_vz_1x3v_tensor_p1, lax_flux_nodal_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, ho_lax_flux_nodal_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, ho_lax_flux_nodal_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, lax_flux_nodal_vx_1x1v_ser_p1, ho_lax_flux_nodal_vx_1x1v_ser_p2, ho_lax_flux_nodal_vx_1x1v_ser_p3 }, // 0
  { NULL, lax_flux_nodal_vx_1x2v_ser_p1, ho_lax_flux_nodal_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vx_1x3v_ser_p1, ho_lax_flux_nodal_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, lax_flux_nodal_vx_2x1v_ser_p1, ho_lax_flux_nodal_vx_2x1v_ser_p2, ho_lax_flux_nodal_vx_2x1v_ser_p3 }, // 3
  { NULL, lax_flux_nodal_vx_2x2v_ser_p1, ho_lax_flux_nodal_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vx_2x3v_ser_p1, ho_lax_flux_nodal_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, lax_flux_nodal_vy_1x2v_ser_p1, ho_lax_flux_nodal_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, lax_flux_nodal_vy_1x3v_ser_p1, ho_lax_flux_nodal_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, lax_flux_nodal_vy_2x2v_ser_p1, ho_lax_flux_nodal_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, lax_flux_nodal_vy_2x3v_ser_p1, ho_lax_flux_nodal_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_lax_flux_arr_kern_list ser_ho_lax_flux_nodal_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, lax_flux_nodal_vz_1x3v_ser_p1, ho_lax_flux_nodal_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, lax_flux_nodal_vz_2x3v_ser_p1, ho_lax_flux_nodal_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, lax_flux_nodal_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_tensor_p1, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_2x1v_tensor_p1, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_hamil_phase_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad_arr, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_sparse_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// Tensor p=1 hybrid triad velocity-direction producers: dense velocity-space
// Hamiltonian and phase-space Hamiltonian (GR triads).
GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_vel_dense_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vx_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_node, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1_node, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_node, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1_node, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1_node, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1_node, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list tensor_nc_hamil_phase_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};


GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, NULL, NULL }, // 0
  { NULL, B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, NULL, NULL }, // 3
  { NULL, B_ho_hamil_phase_alpha_quad_vx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vy_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_phase_alpha_quad_vy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_phase_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1, hamil_phase_alpha_quad_vx_1x1v_ser_p2, hamil_phase_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1, hamil_phase_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1, hamil_phase_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1, hamil_phase_alpha_quad_vx_2x1v_ser_p2, hamil_phase_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1, hamil_phase_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1, hamil_phase_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, hamil_phase_alpha_quad_vx_1x1v_ser_p1, hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2, hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, hamil_phase_alpha_quad_vx_1x2v_ser_p1, hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, hamil_phase_alpha_quad_vx_1x3v_ser_p1, hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, hamil_phase_alpha_quad_vx_2x1v_ser_p1, hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2, hamil_phase_ho_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, hamil_phase_alpha_quad_vx_2x2v_ser_p1, hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, hamil_phase_alpha_quad_vy_2x2v_ser_p1, hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_hamil_phase_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_hamil_alpha_quad_arr, no_hamil_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vx_1x1v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_dense_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1, nc_hamil_vel_dense_ho_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_vel_sparse_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1, nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_1x1v_ser_p1, nc_hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vx_1x2v_ser_p1, nc_hamil_phase_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vx_1x3v_ser_p1, nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1, nc_hamil_phase_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, nc_hamil_phase_alpha_quad_vy_1x2v_ser_p1, nc_hamil_phase_ho_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vy_1x3v_ser_p1, nc_hamil_phase_ho_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, nc_hamil_phase_alpha_quad_vy_2x2v_ser_p1, nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_hamil_alpha_quad_arr_kern_list ser_nc_hamil_phase_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, nc_hamil_phase_alpha_quad_vz_1x3v_ser_p1, nc_hamil_phase_ho_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, nc_hamil_phase_alpha_quad_vz_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, nc_hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_ser_p1, E_alpha_quad_vx_1x1v_ser_p2, E_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, E_alpha_quad_vx_1x2v_ser_p1, E_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_ser_p1, E_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_ser_p1, E_alpha_quad_vx_2x1v_ser_p2, E_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, E_alpha_quad_vx_2x2v_ser_p1, E_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_ser_p1, E_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_ser_p1, E_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_ser_p1, E_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_ser_p1, E_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_ser_p1, E_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_ser_p1, E_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_ser_p1, E_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_ser_p1, E_ho_alpha_quad_vx_1x1v_ser_p2, E_ho_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, E_alpha_quad_vx_1x2v_ser_p1, E_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_ser_p1, E_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_ser_p1, E_ho_alpha_quad_vx_2x1v_ser_p2, E_ho_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, E_alpha_quad_vx_2x2v_ser_p1, E_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_ser_p1, E_ho_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_ser_p1, E_ho_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_ser_p1, E_ho_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_ser_p1, E_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_ser_p1, E_ho_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list ser_E_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_ser_p1, E_ho_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_ser_p1, E_ho_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, E_alpha_quad_vx_1x1v_tensor_p1, E_alpha_quad_vx_1x1v_tensor_p2, E_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, E_alpha_quad_vx_1x2v_tensor_p1, E_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, E_alpha_quad_vx_1x3v_tensor_p1, E_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, E_alpha_quad_vx_2x1v_tensor_p1, E_alpha_quad_vx_2x1v_tensor_p2, E_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, E_alpha_quad_vx_2x2v_tensor_p1, E_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, E_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, E_ho_alpha_quad_vx_1x1v_tensor_p1, E_alpha_quad_vx_1x1v_tensor_p2, E_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, E_ho_alpha_quad_vx_1x2v_tensor_p1, E_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, E_ho_alpha_quad_vx_1x3v_tensor_p1, E_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, E_ho_alpha_quad_vx_2x1v_tensor_p1, E_alpha_quad_vx_2x1v_tensor_p2, E_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, E_ho_alpha_quad_vx_2x2v_tensor_p1, E_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, E_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_alpha_quad_vy_1x2v_tensor_p1, E_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, E_alpha_quad_vy_1x3v_tensor_p1, E_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_alpha_quad_vy_2x2v_tensor_p1, E_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, E_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, E_ho_alpha_quad_vy_1x2v_tensor_p1, E_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, E_ho_alpha_quad_vy_1x3v_tensor_p1, E_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, E_ho_alpha_quad_vy_2x2v_tensor_p1, E_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, E_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_alpha_quad_vz_1x3v_tensor_p1, E_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_E_alpha_quad_arr_kern_list tensor_E_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, E_ho_alpha_quad_vz_1x3v_tensor_p1, E_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, E_ho_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, E_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_ser_p1, phi_alpha_quad_vx_1x1v_ser_p2, phi_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_ser_p1, phi_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_ser_p1, phi_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_ser_p1, phi_alpha_quad_vx_2x1v_ser_p2, phi_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_ser_p1, phi_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_ser_p1, phi_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_ser_p1, phi_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_ser_p1, phi_ho_alpha_quad_vx_1x1v_ser_p2, phi_ho_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_ser_p1, phi_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_ser_p1, phi_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_ser_p1, phi_ho_alpha_quad_vx_2x1v_ser_p2, phi_ho_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_ser_p1, phi_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_ser_p1, phi_ho_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_ser_p1, phi_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list ser_phi_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, phi_alpha_quad_vx_1x1v_tensor_p1, phi_alpha_quad_vx_1x1v_tensor_p2, phi_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, phi_alpha_quad_vx_1x2v_tensor_p1, phi_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, phi_alpha_quad_vx_1x3v_tensor_p1, phi_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, phi_alpha_quad_vx_2x1v_tensor_p1, phi_alpha_quad_vx_2x1v_tensor_p2, phi_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, phi_alpha_quad_vx_2x2v_tensor_p1, phi_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, phi_ho_alpha_quad_vx_1x1v_tensor_p1, phi_alpha_quad_vx_1x1v_tensor_p2, phi_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, phi_ho_alpha_quad_vx_1x2v_tensor_p1, phi_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, phi_ho_alpha_quad_vx_1x3v_tensor_p1, phi_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, phi_ho_alpha_quad_vx_2x1v_tensor_p1, phi_alpha_quad_vx_2x1v_tensor_p2, phi_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, phi_ho_alpha_quad_vx_2x2v_tensor_p1, phi_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, phi_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_alpha_quad_vy_2x2v_tensor_p1, phi_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, phi_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, phi_ho_alpha_quad_vy_2x2v_tensor_p1, phi_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, phi_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, phi_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_phi_alpha_quad_arr_kern_list tensor_phi_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, no_phi_alpha_quad_arr, no_phi_alpha_quad_arr, NULL }, // 5
  // 3x kernels
  { NULL, phi_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p1, B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p1, B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Whole-surface (arr) variants: alpha_v at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1, B_hamil_phase_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1, B_hamil_phase_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1, B_hamil_phase_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1, B_hamil_phase_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1, B_hamil_phase_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1, B_hamil_phase_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1, B_hamil_phase_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Whole-surface (arr) variants: alpha_v at quadrature points for the magnetic field Lorentz force with a phase-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_phase_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_phase_alpha_quad_vx_1x2v_ser_p1, B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vx_1x3v_ser_p1, B_ho_hamil_phase_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_phase_alpha_quad_vx_2x2v_ser_p1, B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vx_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_phase_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_phase_alpha_quad_vy_1x2v_ser_p1, B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vy_1x3v_ser_p1, B_ho_hamil_phase_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_phase_alpha_quad_vy_2x2v_ser_p1, B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vy_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list ser_B_ho_hamil_phase_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_phase_alpha_quad_vz_1x3v_ser_p1, B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_phase_alpha_quad_vz_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_phase_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 0
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_alpha_quad_arr, no_B_alpha_quad_arr, no_B_alpha_quad_arr }, // 3
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1, B_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_dense_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p1, B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_B_alpha_quad_arr_kern_list tensor_B_ho_hamil_vel_sparse_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p1, B_hamil_vel_sparse_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, B_ho_hamil_vel_sparse_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_ser_p1, rad_alpha_quad_vx_1x1v_ser_p2, rad_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_ser_p1, rad_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_ser_p1, rad_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_ser_p1, rad_alpha_quad_vx_2x1v_ser_p2, rad_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_ser_p1, rad_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_ser_p1, rad_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_ser_p1, rad_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_ser_p1, rad_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_ser_p1, rad_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_ser_p1, rad_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_ser_p1, rad_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_ser_p1, rad_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_ser_p1, rad_ho_alpha_quad_vx_1x1v_ser_p2, rad_ho_alpha_quad_vx_1x1v_ser_p3 }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_ser_p1, rad_ho_alpha_quad_vx_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_ser_p1, rad_ho_alpha_quad_vx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_ser_p1, rad_ho_alpha_quad_vx_2x1v_ser_p2, rad_ho_alpha_quad_vx_2x1v_ser_p3 }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_ser_p1, rad_ho_alpha_quad_vx_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_ser_p1, rad_ho_alpha_quad_vx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_ser_p1, rad_ho_alpha_quad_vy_1x2v_ser_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_ser_p1, rad_ho_alpha_quad_vy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_ser_p1, rad_ho_alpha_quad_vy_2x2v_ser_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_ser_p1, rad_ho_alpha_quad_vy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list ser_rad_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_ser_p1, rad_ho_alpha_quad_vz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_ser_p1, rad_ho_alpha_quad_vz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_ser_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, rad_alpha_quad_vx_1x1v_tensor_p1, rad_alpha_quad_vx_1x1v_tensor_p2, rad_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, rad_alpha_quad_vx_1x2v_tensor_p1, rad_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vx_1x3v_tensor_p1, rad_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, rad_alpha_quad_vx_2x1v_tensor_p1, rad_alpha_quad_vx_2x1v_tensor_p2, rad_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, rad_alpha_quad_vx_2x2v_tensor_p1, rad_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_ho_alpha_quad_vx_arr_kernels[] = {
  // 1x kernels
  { NULL, rad_ho_alpha_quad_vx_1x1v_tensor_p1, rad_alpha_quad_vx_1x1v_tensor_p2, rad_alpha_quad_vx_1x1v_tensor_p3 }, // 0
  { NULL, rad_ho_alpha_quad_vx_1x2v_tensor_p1, rad_alpha_quad_vx_1x2v_tensor_p2, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vx_1x3v_tensor_p1, rad_alpha_quad_vx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, rad_ho_alpha_quad_vx_2x1v_tensor_p1, rad_alpha_quad_vx_2x1v_tensor_p2, rad_alpha_quad_vx_2x1v_tensor_p3 }, // 3
  { NULL, rad_ho_alpha_quad_vx_2x2v_tensor_p1, rad_alpha_quad_vx_2x2v_tensor_p2, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vx_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vx_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_alpha_quad_vy_1x2v_tensor_p1, rad_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, rad_alpha_quad_vy_1x3v_tensor_p1, rad_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_alpha_quad_vy_2x2v_tensor_p1, rad_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, rad_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_ho_alpha_quad_vy_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, rad_ho_alpha_quad_vy_1x2v_tensor_p1, rad_alpha_quad_vy_1x2v_tensor_p2, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vy_1x3v_tensor_p1, rad_alpha_quad_vy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, rad_ho_alpha_quad_vy_2x2v_tensor_p1, rad_alpha_quad_vy_2x2v_tensor_p2, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vy_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vy_3x3v_tensor_p1, NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_alpha_quad_vz_1x3v_tensor_p1, rad_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};

// High-order variant: at p=1 (the tensor hybrid, the only tensor basis
// with distinct lo/ho variants) this holds the ho_ kernels; tensor p>1
// entries are the single (high-order by design) variant shared with the
// plain tensor list.
GKYL_CU_D
static const gkyl_rad_alpha_quad_arr_kern_list tensor_rad_ho_alpha_quad_vz_arr_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, rad_ho_alpha_quad_vz_1x3v_tensor_p1, rad_alpha_quad_vz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, rad_ho_alpha_quad_vz_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, rad_ho_alpha_quad_vz_3x3v_tensor_p1, NULL, NULL }, // 6
};
