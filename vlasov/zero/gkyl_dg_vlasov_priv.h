#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_hyb_build.h>

// Types for various kernels
typedef void (*hamil_vol_t)(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out);

typedef void (*E_vol_t)(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out);

typedef void (*phi_vol_t)(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out);

typedef void (*B_vol_t)(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out);

typedef void (*rad_vol_t)(const double *w, const double *dxv, 
  const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out);

typedef double (*vlasov_stream_surf_t)(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_stream_boundary_surf_t)(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

typedef double (*vlasov_stream_surf_from_flux_t)(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out);

typedef double (*vlasov_stream_boundary_surf_from_flux_t)(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out);

typedef double (*vlasov_accel_surf_t)(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out);

typedef double (*vlasov_accel_boundary_surf_t)(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static const struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

// Volume forces tables. 
typedef struct { hamil_vol_t kernels[4]; } gkyl_dg_vlasov_hamil_vol_kern_list;
typedef struct { E_vol_t kernels[4]; } gkyl_dg_vlasov_E_vol_kern_list;
typedef struct { phi_vol_t kernels[4]; } gkyl_dg_vlasov_phi_vol_kern_list;
typedef struct { B_vol_t kernels[4]; } gkyl_dg_vlasov_B_vol_kern_list;
typedef struct { rad_vol_t kernels[4]; } gkyl_dg_vlasov_rad_vol_kern_list;

// Surface kernel type tables. 
typedef struct { vlasov_stream_surf_t kernels[4]; } gkyl_dg_vlasov_stream_surf_kern_list;
typedef struct { vlasov_stream_boundary_surf_t kernels[4]; } gkyl_dg_vlasov_stream_boundary_surf_kern_list;
typedef struct { vlasov_stream_surf_from_flux_t kernels[4]; } gkyl_dg_vlasov_stream_surf_from_flux_kern_list;
typedef struct { vlasov_stream_boundary_surf_from_flux_t kernels[4]; } gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list;
typedef struct { vlasov_accel_surf_t kernels[4]; } gkyl_dg_vlasov_accel_surf_kern_list;
typedef struct { vlasov_accel_boundary_surf_t kernels[4]; } gkyl_dg_vlasov_accel_boundary_surf_kern_list;

struct dg_vlasov {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped. 
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index.
  bool use_vmap; // If Using vmap. 
  bool use_conf_flux_surf; // If Using configuration-space fluxes for streaming terms (Triads). 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  struct gkyl_range conf_range; // Configuration-space range for use in indexing forces (EM fields, potentials, etc.)
  struct gkyl_range vel_range; // Velocity-space range for use in velocity-space Jacobian. 
  struct gkyl_range phase_range; // Range for indexing velocity-space flux.
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object (acquired host-side for lifetime safety; 0 if not given).
  const struct gkyl_array *jacob_vel; // Velocity-space Jacobian (borrowed from vel_map; device pointer on GPUs).
  const struct gkyl_array *vmap; // Velocity map (4-slot cubic rep per direction; borrowed from vel_map; device pointer on GPUs).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety; 0 if not given).
  const struct gkyl_array *jacob_pos; // Configuration-space (position-map) Jacobian, per-conf-cell constant (borrowed from pos_map; device pointer on GPUs; defined on the extended conf range).
  const struct gkyl_array *poisson_tensor_conf; // Hamiltonian utilized to compute advection in configuration and velocity space.
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute advection in configuration and velocity space. 
  const struct gkyl_array *qmem; // q/m*(E,B) electromagnetic fields (including external electromagnetic fields and forces).
  const struct gkyl_array *pot_tot; // (q/m*(phi + phi_ext) + m*phi_g, q/m*A_ext) total potentials. 
  const struct gkyl_array *rad; // Radiation force.
  const struct gkyl_array *conf_flux_surf; // Nodal expansion of fluxes at configuration space surfaces. 
  const struct gkyl_array *vel_flux_surf; // Nodal expansion of fluxes at velocity space surfaces.  
  const struct gkyl_array *f_no_J; // Distribution function without velocity-space Jacobian. Used by magnetic field updates.  
  hamil_vol_t hamil_vol; // Volume term for Hamiltonian contribution to update from canonical bracket. 
  E_vol_t E_vol; // Volume term for Lorentz forces due to electric fields. 
  phi_vol_t phi_vol; // Volume term for scalar potential, -grad(phi), forces. 
  B_vol_t Bx_vol; // Volume term for Lorentz forces due to magnetic fields in the x direction. 
  B_vol_t By_vol; // Volume term for Lorentz forces due to magnetic fields in the y direction. 
  B_vol_t Bz_vol; // Volume term for Lorentz forces due to magnetic fields in the z direction. 
  rad_vol_t rad_vol; // Volume term for radiation drag forces. 
  vlasov_stream_surf_t stream_surf[3]; // Surface terms for advection in configuration space.
  vlasov_stream_boundary_surf_t stream_boundary_surf[3]; // Boundary surface terms for advection in configuration space.
  vlasov_stream_surf_from_flux_t stream_surf_from_flux[3]; // Surface terms for streaming using fluxes.
  vlasov_stream_boundary_surf_from_flux_t stream_boundary_surf_from_flux[3]; // Surface terms for streaming using fluxes.
  vlasov_accel_surf_t accel_surf[3]; // Surface terms for acceleration.
  vlasov_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
};

// Empty function pointers for cases where these forces do not exist. 
GKYL_CU_DH
static void 
no_E_vol(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out)
{
}
GKYL_CU_DH
static void 
no_phi_vol(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out)
{
}
GKYL_CU_DH
static void 
no_B_vol(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil, const double *qmem, const double *f, double* GKYL_RESTRICT out)
{
}
GKYL_CU_DH
static void 
no_rad_vol(const double *w, const double *dxv, 
  const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out)
{
}

GKYL_CU_DH
static double
no_stream_surf(const double *w, const double *dxv,
  const double *jacob_pos_l, const double *jacob_pos_c, const double *jacob_pos_r, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out)
{
  return 0.0;
}
GKYL_CU_DH
static double
no_stream_boundary_surf(const double *w, const double *dxv,
  const double *jacob_pos_edge, const double *jacob_pos_skin, const double *jacob_vel,
  const double *poisson_tensor_conf, const double *hamil,
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out)
{
  return 0.0;
}

GKYL_CU_DH
static double
vlasov_vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (fabs(qIn[0]) < vlasov->skip_cell_thresh) {
    return 0.0; // Immediately return if phase-space density is below threshold
  }  

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov->pdim-vlasov->cdim; ++i) {
    idx_vel[i] = idx[vlasov->cdim+i];
  }

  int idx_hamil[GKYL_MAX_DIM];
  for (int i=0; i<vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[vlasov->hamil_offset+i];
  }

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov->vel_range, idx_vel);
  long hidx = gkyl_range_idx(&vlasov->hamil_range, idx_hamil);
  long pidx = gkyl_range_idx(&vlasov->phase_range, idx);

  // f_no_J (velocity Jacobians divided out at Gauss-Legendre points) lets the
  // cubic-map kernels express every velocity-Jacobian factor as an EXACT
  // polynomial multiplication by the complementary J's (division-free), which
  // preserves the discrete bracket antisymmetry and hence exact conservation.
  vlasov->hamil_vol(xc, dx,
    vlasov->vmap ? (const double*) gkyl_array_cfetch(vlasov->vmap, vidx) : 0,
    vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidx) : 0,
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    (const double*) gkyl_array_cfetch(vlasov->poisson_tensor_conf, cidx),
    (const double*) gkyl_array_cfetch(vlasov->hamil, hidx),
    qIn, (const double*) gkyl_array_cfetch(vlasov->f_no_J, pidx), qRhsOut);
  vlasov->E_vol(xc, dx, 
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    vlasov->qmem ? (const double*) gkyl_array_cfetch(vlasov->qmem, cidx) : 0, 
    qIn, qRhsOut); 
  vlasov->phi_vol(xc, dx,
    vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidx) : 0,
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    vlasov->pot_tot ? (const double*) gkyl_array_cfetch(vlasov->pot_tot, cidx) : 0,
    qIn, qRhsOut);
  // Nonuniform mesh kernels utilize f without the velocity-space Jacobian to handle
  // the transverse derivatives in the 1/Jvi grad_vi(H) x B cross product. 
  vlasov->Bx_vol(xc, dx, 
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    (const double*) gkyl_array_cfetch(vlasov->hamil, hidx), 
    vlasov->qmem ? (const double*) gkyl_array_cfetch(vlasov->qmem, cidx) : 0, 
    (const double*) gkyl_array_cfetch(vlasov->f_no_J, pidx), 
    qRhsOut); 
  vlasov->By_vol(xc, dx, 
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    (const double*) gkyl_array_cfetch(vlasov->hamil, hidx), 
    vlasov->qmem ? (const double*) gkyl_array_cfetch(vlasov->qmem, cidx) : 0, 
    (const double*) gkyl_array_cfetch(vlasov->f_no_J, pidx), 
    qRhsOut); 
  vlasov->Bz_vol(xc, dx, 
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    (const double*) gkyl_array_cfetch(vlasov->hamil, hidx), 
    vlasov->qmem ? (const double*) gkyl_array_cfetch(vlasov->qmem, cidx) : 0, 
    (const double*) gkyl_array_cfetch(vlasov->f_no_J, pidx), 
    qRhsOut); 
  vlasov->rad_vol(xc, dx, 
    vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
    vlasov->rad ? (const double*) gkyl_array_cfetch(vlasov->rad, vidx) : 0,
    qIn, qRhsOut); 

  return 0.0;
}

// Hamiltonian, velocity space only, volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_vol_1x1v_ser_p1, vlasov_hamil_vel_dense_vol_1x1v_ser_p2, vlasov_hamil_vel_dense_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_vol_1x2v_ser_p1, vlasov_hamil_vel_dense_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_vol_1x3v_ser_p1, vlasov_hamil_vel_dense_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_vol_2x1v_ser_p1, vlasov_hamil_vel_dense_vol_2x1v_ser_p2, vlasov_hamil_vel_dense_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_vol_2x2v_ser_p1, vlasov_hamil_vel_dense_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_dense_vol_2x3v_ser_p1, vlasov_hamil_vel_dense_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_vol_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_vol_1x1v_ser_p1, vlasov_hamil_vel_dense_vol_1x1v_ser_p2, vlasov_hamil_vel_dense_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_vol_1x2v_ser_p1, vlasov_hamil_vel_sparse_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_vol_1x3v_ser_p1, vlasov_hamil_vel_sparse_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_vol_2x1v_ser_p1, vlasov_hamil_vel_dense_vol_2x1v_ser_p2, vlasov_hamil_vel_dense_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_vol_2x2v_ser_p1, vlasov_hamil_vel_sparse_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_sparse_vol_2x3v_ser_p1, vlasov_hamil_vel_sparse_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Hamiltonian, general phase space, volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_phase_vol_1x1v_ser_p1, vlasov_hamil_phase_vol_1x1v_ser_p2, vlasov_hamil_phase_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_phase_vol_1x2v_ser_p1, vlasov_hamil_phase_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_phase_vol_1x3v_ser_p1, vlasov_hamil_phase_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_phase_vol_2x1v_ser_p1, vlasov_hamil_phase_vol_2x1v_ser_p2, vlasov_hamil_phase_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_phase_vol_2x2v_ser_p1, vlasov_hamil_phase_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_phase_vol_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_phase_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Phase-space Hamiltonian volume kernels (Tensor basis): only the p=1
// tensor hybrid has a phase-space Hamiltonian representation.
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_phase_vol_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_hamil_phase_vol_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_hamil_phase_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_phase_vol_2x1v_tensor_p1, NULL, NULL }, // 3
  { NULL, vlasov_hamil_phase_vol_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_phase_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(vlasov_hamil_phase_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Non-canonical Poisson Bracket, hamil vel space dependance only, volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_nc_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p1, vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x2v_ser_p1, vlasov_nc_hamil_vel_dense_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x3v_ser_p1, vlasov_nc_hamil_vel_dense_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_vel_dense_vol_2x2v_ser_p1, vlasov_nc_hamil_vel_dense_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_nc_hamil_vel_dense_vol_2x3v_ser_p1, vlasov_nc_hamil_vel_dense_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_nc_hamil_vel_dense_vol_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_nc_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p1, vlasov_nc_hamil_vel_dense_vol_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_nc_hamil_vel_sparse_vol_1x2v_ser_p1, vlasov_nc_hamil_vel_sparse_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_nc_hamil_vel_sparse_vol_1x3v_ser_p1, vlasov_nc_hamil_vel_sparse_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_vel_sparse_vol_2x2v_ser_p1, vlasov_nc_hamil_vel_sparse_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_nc_hamil_vel_sparse_vol_2x3v_ser_p1, vlasov_nc_hamil_vel_sparse_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_nc_hamil_vel_sparse_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Tensor p=1 hybrid triad (non-canonical bracket) volume kernels: the C^1
// cubic velocity map's quadratic Jacobian is divided nodally per derivative
// direction in-kernel (1V sparse aliases dense, as for Serendipity).
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_nc_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_nc_hamil_vel_sparse_vol_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_nc_hamil_vel_sparse_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_vel_sparse_vol_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_nc_hamil_vel_sparse_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_nc_hamil_vel_sparse_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Tensor p=1 hybrid triad volumes, dense velocity-space Hamiltonian.
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_nc_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_nc_hamil_vel_dense_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_vel_dense_vol_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_nc_hamil_vel_dense_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_nc_hamil_vel_dense_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Tensor p=1 hybrid triad volumes, phase-space Hamiltonian (GR triads).
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_nc_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_phase_vol_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_nc_hamil_phase_vol_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_nc_hamil_phase_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_phase_vol_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_nc_hamil_phase_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V_PHASE(vlasov_nc_hamil_phase_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Non-canonical Poisson Bracket, hamil phase space dependance, volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list ser_nc_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_nc_hamil_phase_vol_1x1v_ser_p1, vlasov_nc_hamil_phase_vol_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_nc_hamil_phase_vol_1x2v_ser_p1, vlasov_nc_hamil_phase_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_nc_hamil_phase_vol_1x3v_ser_p1, vlasov_nc_hamil_phase_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_nc_hamil_phase_vol_2x2v_ser_p1, vlasov_nc_hamil_phase_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_nc_hamil_phase_vol_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_nc_hamil_phase_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Hamiltonian, velocity space only, volume kernels (Tensor basis).  
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_vol_1x1v_tensor_p1, vlasov_hamil_vel_dense_vol_1x1v_tensor_p2, vlasov_hamil_vel_dense_vol_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_vol_1x2v_tensor_p1, vlasov_hamil_vel_dense_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_vol_1x3v_tensor_p1, vlasov_hamil_vel_dense_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_vol_2x1v_tensor_p1, vlasov_hamil_vel_dense_vol_2x1v_tensor_p2, vlasov_hamil_vel_dense_vol_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_vol_2x2v_tensor_p1, vlasov_hamil_vel_dense_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_dense_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_hamil_vol_kern_list tensor_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_vol_1x1v_tensor_p1, vlasov_hamil_vel_dense_vol_1x1v_tensor_p2, vlasov_hamil_vel_dense_vol_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_vol_1x2v_tensor_p1, vlasov_hamil_vel_sparse_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_vol_1x3v_tensor_p1, vlasov_hamil_vel_sparse_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_vol_2x1v_tensor_p1, vlasov_hamil_vel_dense_vol_2x1v_tensor_p2, vlasov_hamil_vel_dense_vol_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_vol_2x2v_tensor_p1, vlasov_hamil_vel_sparse_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_sparse_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Electric field Lorentz force volume kernels with velocity-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_E_vol_kern_list ser_E_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_E_vol_1x1v_ser_p1, vlasov_E_vol_1x1v_ser_p2, vlasov_E_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_E_vol_1x2v_ser_p1, vlasov_E_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_E_vol_1x3v_ser_p1, vlasov_E_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_E_vol_2x1v_ser_p1, vlasov_E_vol_2x1v_ser_p2, vlasov_E_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_E_vol_2x2v_ser_p1, vlasov_E_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_E_vol_2x3v_ser_p1, vlasov_E_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_E_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Electric field Lorentz force volume kernels with velocity-space Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_E_vol_kern_list tensor_E_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_E_vol_1x1v_tensor_p1, vlasov_E_vol_1x1v_tensor_p2, vlasov_E_vol_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_E_vol_1x2v_tensor_p1, vlasov_E_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_E_vol_1x3v_tensor_p1, vlasov_E_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_E_vol_2x1v_tensor_p1, vlasov_E_vol_2x1v_tensor_p2, vlasov_E_vol_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_E_vol_2x2v_tensor_p1, vlasov_E_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_E_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_E_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Scalar potential, -grad(phi), force volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_phi_vol_kern_list ser_phi_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_phi_vol_1x1v_ser_p1, vlasov_phi_vol_1x1v_ser_p2, vlasov_phi_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_phi_vol_1x2v_ser_p1, vlasov_phi_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_phi_vol_1x3v_ser_p1, vlasov_phi_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_phi_vol_2x1v_ser_p1, vlasov_phi_vol_2x1v_ser_p2, vlasov_phi_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_phi_vol_2x2v_ser_p1, vlasov_phi_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_phi_vol_2x3v_ser_p1, vlasov_phi_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_phi_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Scalar potential, -grad(phi), force volume kernels (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_phi_vol_kern_list tensor_phi_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_phi_vol_1x1v_tensor_p1, vlasov_phi_vol_1x1v_tensor_p2, vlasov_phi_vol_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_phi_vol_1x2v_tensor_p1, vlasov_phi_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_phi_vol_1x3v_tensor_p1, vlasov_phi_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_phi_vol_2x1v_tensor_p1, vlasov_phi_vol_2x1v_tensor_p2, vlasov_phi_vol_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_phi_vol_2x2v_tensor_p1, vlasov_phi_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_phi_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_phi_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Magnetic field in x-direction Lorentz force volume kernels with velocity-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bx_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_Bx_hamil_vel_dense_vol_1x3v_ser_p1, vlasov_Bx_hamil_vel_dense_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_Bx_hamil_vel_dense_vol_2x3v_ser_p1, vlasov_Bx_hamil_vel_dense_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_Bx_hamil_vel_dense_vol_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bx_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_Bx_hamil_vel_sparse_vol_1x3v_ser_p1, vlasov_Bx_hamil_vel_sparse_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_Bx_hamil_vel_sparse_vol_2x3v_ser_p1, vlasov_Bx_hamil_vel_sparse_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_Bx_hamil_vel_sparse_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Magnetic field in x-direction Lorentz force volume kernels with velocity-space Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bx_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_Bx_hamil_vel_dense_vol_1x3v_tensor_p1, vlasov_Bx_hamil_vel_dense_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bx_hamil_vel_dense_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_Bx_hamil_vel_dense_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bx_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_Bx_hamil_vel_sparse_vol_1x3v_tensor_p1, vlasov_Bx_hamil_vel_sparse_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bx_hamil_vel_sparse_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_Bx_hamil_vel_sparse_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Magnetic field in x-direction Lorentz force volume kernels with general Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bx_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_Bx_hamil_phase_vol_1x3v_ser_p1, vlasov_Bx_hamil_phase_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_Bx_hamil_phase_vol_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Magnetic field in x-direction Lorentz force volume kernels with general Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bx_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_Bx_hamil_phase_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bx_hamil_phase_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Magnetic field in y-direction Lorentz force volume kernels with velocity-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_By_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_By_hamil_vel_dense_vol_1x3v_ser_p1, vlasov_By_hamil_vel_dense_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_By_hamil_vel_dense_vol_2x3v_ser_p1, vlasov_By_hamil_vel_dense_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_By_hamil_vel_dense_vol_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_By_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_By_hamil_vel_sparse_vol_1x3v_ser_p1, vlasov_By_hamil_vel_sparse_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_By_hamil_vel_sparse_vol_2x3v_ser_p1, vlasov_By_hamil_vel_sparse_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_By_hamil_vel_sparse_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Magnetic field in y-direction Lorentz force volume kernels with velocity-space Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_By_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_By_hamil_vel_dense_vol_1x3v_tensor_p1, vlasov_By_hamil_vel_dense_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_By_hamil_vel_dense_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_By_hamil_vel_dense_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_By_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_By_hamil_vel_sparse_vol_1x3v_tensor_p1, vlasov_By_hamil_vel_sparse_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_By_hamil_vel_sparse_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_By_hamil_vel_sparse_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Magnetic field in y-direction Lorentz force volume kernels with general Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_By_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, no_B_vol, no_B_vol, NULL }, // 1
  { NULL, vlasov_By_hamil_phase_vol_1x3v_ser_p1, vlasov_By_hamil_phase_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, no_B_vol, no_B_vol, NULL }, // 4
  { NULL, vlasov_By_hamil_phase_vol_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Magnetic field in y-direction Lorentz force volume kernels with general Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_By_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_By_hamil_phase_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_By_hamil_phase_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Magnetic field in z-direction Lorentz force volume kernels with velocity-space Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bz_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, vlasov_Bz_hamil_vel_dense_vol_1x2v_ser_p1, vlasov_Bz_hamil_vel_dense_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_Bz_hamil_vel_dense_vol_1x3v_ser_p1, vlasov_Bz_hamil_vel_dense_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, vlasov_Bz_hamil_vel_dense_vol_2x2v_ser_p1, vlasov_Bz_hamil_vel_dense_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_Bz_hamil_vel_dense_vol_2x3v_ser_p1, vlasov_Bz_hamil_vel_dense_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_Bz_hamil_vel_dense_vol_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bz_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_1x2v_ser_p1, vlasov_Bz_hamil_vel_sparse_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_1x3v_ser_p1, vlasov_Bz_hamil_vel_sparse_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_2x2v_ser_p1, vlasov_Bz_hamil_vel_sparse_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_2x3v_ser_p1, vlasov_Bz_hamil_vel_sparse_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Magnetic field in z-direction Lorentz force volume kernels with velocity-space Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bz_hamil_vel_dense_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, vlasov_Bz_hamil_vel_dense_vol_1x2v_tensor_p1, vlasov_Bz_hamil_vel_dense_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_Bz_hamil_vel_dense_vol_1x3v_tensor_p1, vlasov_Bz_hamil_vel_dense_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, vlasov_Bz_hamil_vel_dense_vol_2x2v_tensor_p1, vlasov_Bz_hamil_vel_dense_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bz_hamil_vel_dense_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_Bz_hamil_vel_dense_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bz_hamil_vel_sparse_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_1x2v_tensor_p1, vlasov_Bz_hamil_vel_sparse_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_1x3v_tensor_p1, vlasov_Bz_hamil_vel_sparse_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, vlasov_Bz_hamil_vel_sparse_vol_2x2v_tensor_p1, vlasov_Bz_hamil_vel_sparse_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bz_hamil_vel_sparse_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_Bz_hamil_vel_sparse_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Magnetic field in z-direction Lorentz force volume kernels with general Hamiltonian (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list ser_Bz_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 0
  { NULL, vlasov_Bz_hamil_phase_vol_1x2v_ser_p1, vlasov_Bz_hamil_phase_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_Bz_hamil_phase_vol_1x3v_ser_p1, vlasov_Bz_hamil_phase_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, no_B_vol, no_B_vol, no_B_vol }, // 3
  { NULL, vlasov_Bz_hamil_phase_vol_2x2v_ser_p1, vlasov_Bz_hamil_phase_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_Bz_hamil_phase_vol_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Magnetic field in z-direction Lorentz force volume kernels with general Hamiltonian (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_B_vol_kern_list tensor_Bz_hamil_phase_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_Bz_hamil_phase_vol_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_Bz_hamil_phase_vol_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_Bz_hamil_phase_vol_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_Bz_hamil_phase_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, NULL, NULL, NULL }, // 6
};

// Radiation force volume kernels (Serendipity basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_rad_vol_kern_list ser_rad_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_rad_vol_1x1v_ser_p1, vlasov_rad_vol_1x1v_ser_p2, vlasov_rad_vol_1x1v_ser_p3 }, // 0
  { NULL, vlasov_rad_vol_1x2v_ser_p1, vlasov_rad_vol_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_rad_vol_1x3v_ser_p1, vlasov_rad_vol_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_rad_vol_2x1v_ser_p1, vlasov_rad_vol_2x1v_ser_p2, vlasov_rad_vol_2x1v_ser_p3 }, // 3
  { NULL, vlasov_rad_vol_2x2v_ser_p1, vlasov_rad_vol_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_rad_vol_2x3v_ser_p1, vlasov_rad_vol_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_rad_vol_3x3v_ser_p1, NULL, NULL }, // 6
};

// Radiation force volume kernels (Tensor basis). 
GKYL_CU_D
static const gkyl_dg_vlasov_rad_vol_kern_list tensor_rad_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_rad_vol_1x1v_tensor_p1, vlasov_rad_vol_1x1v_tensor_p2, vlasov_rad_vol_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_rad_vol_1x2v_tensor_p1, vlasov_rad_vol_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_rad_vol_1x3v_tensor_p1, vlasov_rad_vol_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_rad_vol_2x1v_tensor_p1, vlasov_rad_vol_2x1v_tensor_p2, vlasov_rad_vol_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_rad_vol_2x2v_tensor_p1, vlasov_rad_vol_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_rad_vol_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_rad_vol_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_dense_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_1x1v_ser_p1, vlasov_hamil_vel_dense_surfx_1x1v_ser_p2, vlasov_hamil_vel_dense_surfx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_surfx_1x2v_ser_p1, vlasov_hamil_vel_dense_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_surfx_1x3v_ser_p1, vlasov_hamil_vel_dense_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_2x1v_ser_p1, vlasov_hamil_vel_dense_surfx_2x1v_ser_p2, vlasov_hamil_vel_dense_surfx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_surfx_2x2v_ser_p1, vlasov_hamil_vel_dense_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_dense_surfx_2x3v_ser_p1, vlasov_hamil_vel_dense_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_sparse_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_1x1v_ser_p1, vlasov_hamil_vel_dense_surfx_1x1v_ser_p2, vlasov_hamil_vel_dense_surfx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_surfx_1x2v_ser_p1, vlasov_hamil_vel_sparse_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_surfx_1x3v_ser_p1, vlasov_hamil_vel_sparse_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_2x1v_ser_p1, vlasov_hamil_vel_dense_surfx_2x1v_ser_p2, vlasov_hamil_vel_dense_surfx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_surfx_2x2v_ser_p1, vlasov_hamil_vel_sparse_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_sparse_surfx_2x3v_ser_p1, vlasov_hamil_vel_sparse_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_dense_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p1, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p2, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_surfx_1x2v_tensor_p1, vlasov_hamil_vel_dense_surfx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_surfx_1x3v_tensor_p1, vlasov_hamil_vel_dense_surfx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p1, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p2, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_surfx_2x2v_tensor_p1, vlasov_hamil_vel_dense_surfx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_dense_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_surfx_3x3v_tensor_p1), NULL, NULL}, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_sparse_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p1, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p2, vlasov_hamil_vel_dense_surfx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_surfx_1x2v_tensor_p1, vlasov_hamil_vel_sparse_surfx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_surfx_1x3v_tensor_p1, vlasov_hamil_vel_sparse_surfx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p1, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p2, vlasov_hamil_vel_dense_surfx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_surfx_2x2v_tensor_p1, vlasov_hamil_vel_sparse_surfx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_sparse_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_surfx_3x3v_tensor_p1), NULL, NULL}, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_dense_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_surf, no_stream_surf, no_stream_surf }, // 3  
  { NULL, vlasov_hamil_vel_dense_surfy_2x2v_ser_p1, vlasov_hamil_vel_dense_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_dense_surfy_2x3v_ser_p1, vlasov_hamil_vel_dense_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_sparse_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_surf, no_stream_surf, no_stream_surf }, // 3  
  { NULL, vlasov_hamil_vel_sparse_surfy_2x2v_ser_p1, vlasov_hamil_vel_sparse_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_sparse_surfy_2x3v_ser_p1, vlasov_hamil_vel_sparse_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_dense_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_surf, no_stream_surf, no_stream_surf }, // 3  
  { NULL, vlasov_hamil_vel_dense_surfy_2x2v_tensor_p1, vlasov_hamil_vel_dense_surfy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_dense_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_sparse_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_surf, no_stream_surf, no_stream_surf }, // 3  
  { NULL, vlasov_hamil_vel_sparse_surfy_2x2v_tensor_p1, vlasov_hamil_vel_sparse_surfy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_sparse_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_dense_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_hamil_vel_sparse_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_dense_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_hamil_vel_sparse_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};


// Streaming general (NC) Hamiltonian surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_surfx_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfx_2x3v_ser_p1, vlasov_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general Hamiltonian surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfy_2x3v_ser_p1, vlasov_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general Hamiltonian surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_ho_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_ho_surfx_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_ho_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_ho_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_ho_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfx_2x3v_ser_p1, vlasov_ho_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general Hamiltonian surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_ho_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_ho_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfy_2x3v_ser_p1, vlasov_ho_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general Hamiltonian surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list ser_stream_ho_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
// 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};


// Streaming from-flux surface kernels (Tensor basis): only the p=1 tensor
// hybrid uses the from-flux streaming path (canonical/phase-Hamiltonian
// models); lo list = plain p1 kernels, ho list = ho_ variants.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_surfx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_surfx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfx_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_ho_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_ho_surfx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_ho_surfx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_ho_surfx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_ho_surfx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_ho_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfx_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_ho_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_ho_surfy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_ho_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_from_flux_kern_list tensor_stream_ho_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfx_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_ho_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_ho_surfx_1x1v_tensor_p1, NULL, NULL }, // 0
  { NULL, vlasov_boundary_ho_surfx_1x2v_tensor_p1, NULL, NULL }, // 1
  { NULL, vlasov_boundary_ho_surfx_1x3v_tensor_p1, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_ho_surfx_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_ho_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfx_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_ho_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_ho_surfy_2x2v_tensor_p1, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_ho_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};

GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list tensor_stream_boundary_ho_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_ser_p1, vlasov_surfvx_1x1v_ser_p2, vlasov_surfvx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_surfvx_1x2v_ser_p1, vlasov_surfvx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfvx_1x3v_ser_p1, vlasov_surfvx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x1v_ser_p1, vlasov_surfvx_2x1v_ser_p2, vlasov_surfvx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_surfvx_2x2v_ser_p1, vlasov_surfvx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfvx_2x3v_ser_p1, vlasov_surfvx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list: vx-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_tensor_p1, vlasov_surfvx_1x1v_tensor_p2, vlasov_surfvx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_surfvx_1x2v_tensor_p1, vlasov_surfvx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_surfvx_1x3v_tensor_p1, vlasov_surfvx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x1v_tensor_p1, vlasov_surfvx_2x1v_tensor_p2, vlasov_surfvx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_surfvx_2x2v_tensor_p1, vlasov_surfvx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_surfvx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfvx_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_ser_p1, vlasov_surfvy_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfvy_1x3v_ser_p1, vlasov_surfvy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvy_2x2v_ser_p1, vlasov_surfvy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfvy_2x3v_ser_p1, vlasov_surfvy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list: vy-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_tensor_p1, vlasov_surfvy_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_surfvy_1x3v_tensor_p1, vlasov_surfvy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvy_2x2v_tensor_p1, vlasov_surfvy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_surfvy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfvy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_ser_p1, vlasov_surfvz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, vlasov_surfvz_2x3v_ser_p1, vlasov_surfvz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list: vz-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_tensor_p1, vlasov_surfvz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_surfvz_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_surfvz_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_ho_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_ser_p1, vlasov_ho_surfvx_1x1v_ser_p2, vlasov_ho_surfvx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_surfvx_1x2v_ser_p1, vlasov_ho_surfvx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfvx_1x3v_ser_p1, vlasov_ho_surfvx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x1v_ser_p1, vlasov_ho_surfvx_2x1v_ser_p2, vlasov_ho_surfvx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_surfvx_2x2v_ser_p1, vlasov_ho_surfvx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfvx_2x3v_ser_p1, vlasov_ho_surfvx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_ho_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_ser_p1, vlasov_ho_surfvy_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_surfvy_1x3v_ser_p1, vlasov_ho_surfvy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvy_2x2v_ser_p1, vlasov_ho_surfvy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_surfvy_2x3v_ser_p1, vlasov_ho_surfvy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_ho_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_ser_p1, vlasov_ho_surfvz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, vlasov_surfvz_2x3v_ser_p1, vlasov_ho_surfvz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration surface kernel list, high-order variant: vx-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_ho_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_ho_surfvx_1x1v_tensor_p1, vlasov_surfvx_1x1v_tensor_p2, vlasov_surfvx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_ho_surfvx_1x2v_tensor_p1, vlasov_surfvx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_ho_surfvx_1x3v_tensor_p1, vlasov_surfvx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_ho_surfvx_2x1v_tensor_p1, vlasov_surfvx_2x1v_tensor_p2, vlasov_surfvx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_ho_surfvx_2x2v_tensor_p1, vlasov_surfvx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_ho_surfvx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfvx_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list, high-order variant: vy-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_ho_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_ho_surfvy_1x2v_tensor_p1, vlasov_surfvy_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_ho_surfvy_1x3v_tensor_p1, vlasov_surfvy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_ho_surfvy_2x2v_tensor_p1, vlasov_surfvy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_ho_surfvy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfvy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration surface kernel list, high-order variant: vz-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_ho_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_ho_surfvz_1x3v_tensor_p1, vlasov_surfvz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_ho_surfvz_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_ho_surfvz_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_dense_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p2, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x2v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x3v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p2, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x2v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x3v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_sparse_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p2, vlasov_hamil_vel_dense_boundary_surfx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_1x2v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_1x3v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p1, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p2, vlasov_hamil_vel_dense_boundary_surfx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_2x2v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_2x3v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_dense_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p2, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x2v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x3v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p2, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x2v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_dense_boundary_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_boundary_surfx_3x3v_tensor_p1), NULL, NULL}, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_sparse_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p2, vlasov_hamil_vel_dense_boundary_surfx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_1x2v_tensor_p1, vlasov_hamil_vel_sparse_boundary_surfx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_1x3v_tensor_p1, vlasov_hamil_vel_sparse_boundary_surfx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p2, vlasov_hamil_vel_dense_boundary_surfx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_hamil_vel_sparse_boundary_surfx_2x2v_tensor_p1, vlasov_hamil_vel_sparse_boundary_surfx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_sparse_boundary_surfx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_boundary_surfx_3x3v_tensor_p1), NULL, NULL}, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_dense_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_boundary_surf, no_stream_boundary_surf, no_stream_boundary_surf }, // 3  
  { NULL, vlasov_hamil_vel_dense_boundary_surfy_2x2v_ser_p1, vlasov_hamil_vel_dense_boundary_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_dense_boundary_surfy_2x3v_ser_p1, vlasov_hamil_vel_dense_boundary_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_sparse_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_boundary_surf, no_stream_boundary_surf, no_stream_boundary_surf }, // 3  
  { NULL, vlasov_hamil_vel_sparse_boundary_surfy_2x2v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_hamil_vel_sparse_boundary_surfy_2x3v_ser_p1, vlasov_hamil_vel_sparse_boundary_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_boundary_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_dense_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_boundary_surf, no_stream_boundary_surf, no_stream_boundary_surf }, // 3  
  { NULL, vlasov_hamil_vel_dense_boundary_surfy_2x2v_tensor_p1, vlasov_hamil_vel_dense_boundary_surfy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_dense_boundary_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_boundary_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_sparse_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, no_stream_boundary_surf, no_stream_boundary_surf, no_stream_boundary_surf }, // 3  
  { NULL, vlasov_hamil_vel_sparse_boundary_surfy_2x2v_tensor_p1, vlasov_hamil_vel_sparse_boundary_surfy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_hamil_vel_sparse_boundary_surfy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_boundary_surfy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_dense_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_dense_boundary_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_hamil_vel_sparse_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_hamil_vel_sparse_boundary_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming velocity-space Hamiltonian boundary surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_dense_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_dense_boundary_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_hamil_vel_sparse_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_hamil_vel_sparse_boundary_surfz_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfx_1x1v_ser_p1, vlasov_boundary_surfx_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_boundary_surfx_1x2v_ser_p1, vlasov_boundary_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfx_1x3v_ser_p1, vlasov_boundary_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfx_2x2v_ser_p1, vlasov_boundary_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfx_2x3v_ser_p1, vlasov_boundary_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, vlasov_boundary_surfy_2x2v_ser_p1, vlasov_boundary_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfy_2x3v_ser_p1, vlasov_boundary_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_ho_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfx_1x1v_ser_p1, vlasov_boundary_ho_surfx_1x1v_ser_p2, NULL }, // 0
  { NULL, vlasov_boundary_surfx_1x2v_ser_p1, vlasov_boundary_ho_surfx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfx_1x3v_ser_p1, vlasov_boundary_ho_surfx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfx_2x2v_ser_p1, vlasov_boundary_ho_surfx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfx_2x3v_ser_p1, vlasov_boundary_ho_surfx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_ho_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, vlasov_boundary_surfy_2x2v_ser_p1, vlasov_boundary_ho_surfy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfy_2x3v_ser_p1, vlasov_boundary_ho_surfy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Streaming general (NC) Hamiltonian boundary surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_from_flux_kern_list ser_stream_boundary_ho_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3  
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, NULL, NULL, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_ser_p1, vlasov_boundary_surfvx_1x1v_ser_p2, vlasov_boundary_surfvx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_ser_p1, vlasov_boundary_surfvx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_ser_p1, vlasov_boundary_surfvx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x1v_ser_p1, vlasov_boundary_surfvx_2x1v_ser_p2, vlasov_boundary_surfvx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_boundary_surfvx_2x2v_ser_p1, vlasov_boundary_surfvx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfvx_2x3v_ser_p1, vlasov_boundary_surfvx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_tensor_p1, vlasov_boundary_surfvx_1x1v_tensor_p2, vlasov_boundary_surfvx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_tensor_p1, vlasov_boundary_surfvx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_tensor_p1, vlasov_boundary_surfvx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x1v_tensor_p1, vlasov_boundary_surfvx_2x1v_tensor_p2, vlasov_boundary_surfvx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_boundary_surfvx_2x2v_tensor_p1, vlasov_boundary_surfvx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_surfvx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfvx_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_ser_p1, vlasov_boundary_surfvy_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_ser_p1, vlasov_boundary_surfvy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvy_2x2v_ser_p1, vlasov_boundary_surfvy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfvy_2x3v_ser_p1, vlasov_boundary_surfvy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_tensor_p1, vlasov_boundary_surfvy_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_tensor_p1, vlasov_boundary_surfvy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvy_2x2v_tensor_p1, vlasov_boundary_surfvy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_surfvy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfvy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_ser_p1, vlasov_boundary_surfvz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, vlasov_boundary_surfvz_2x3v_ser_p1, vlasov_boundary_surfvz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_tensor_p1, vlasov_boundary_surfvz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_surfvz_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_surfvz_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_ho_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_ser_p1, vlasov_boundary_ho_surfvx_1x1v_ser_p2, vlasov_boundary_ho_surfvx_1x1v_ser_p3 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_ser_p1, vlasov_boundary_ho_surfvx_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_ser_p1, vlasov_boundary_ho_surfvx_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x1v_ser_p1, vlasov_boundary_ho_surfvx_2x1v_ser_p2, vlasov_boundary_ho_surfvx_2x1v_ser_p3 }, // 3
  { NULL, vlasov_boundary_surfvx_2x2v_ser_p1, vlasov_boundary_ho_surfvx_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfvx_2x3v_ser_p1, vlasov_boundary_ho_surfvx_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_ho_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_ser_p1, vlasov_boundary_ho_surfvy_1x2v_ser_p2, NULL }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_ser_p1, vlasov_boundary_ho_surfvy_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvy_2x2v_ser_p1, vlasov_boundary_ho_surfvy_2x2v_ser_p2, NULL }, // 4
  { NULL, vlasov_boundary_surfvy_2x3v_ser_p1, vlasov_boundary_ho_surfvy_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_ho_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_ser_p1, vlasov_boundary_ho_surfvz_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, vlasov_boundary_surfvz_2x3v_ser_p1, vlasov_boundary_ho_surfvz_2x3v_ser_p2, NULL }, // 5
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_ser_p1, NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list, high-order variant: vx-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_ho_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_ho_surfvx_1x1v_tensor_p1, vlasov_boundary_surfvx_1x1v_tensor_p2, vlasov_boundary_surfvx_1x1v_tensor_p3 }, // 0
  { NULL, vlasov_boundary_ho_surfvx_1x2v_tensor_p1, vlasov_boundary_surfvx_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_boundary_ho_surfvx_1x3v_tensor_p1, vlasov_boundary_surfvx_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_ho_surfvx_2x1v_tensor_p1, vlasov_boundary_surfvx_2x1v_tensor_p2, vlasov_boundary_surfvx_2x1v_tensor_p3 }, // 3
  { NULL, vlasov_boundary_ho_surfvx_2x2v_tensor_p1, vlasov_boundary_surfvx_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_ho_surfvx_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfvx_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list, high-order variant: vy-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_ho_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_ho_surfvy_1x2v_tensor_p1, vlasov_boundary_surfvy_1x2v_tensor_p2, NULL }, // 1
  { NULL, vlasov_boundary_ho_surfvy_1x3v_tensor_p1, vlasov_boundary_surfvy_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_ho_surfvy_2x2v_tensor_p1, vlasov_boundary_surfvy_2x2v_tensor_p2, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_ho_surfvy_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfvy_3x3v_tensor_p1), NULL, NULL }, // 6
};

// Acceleration boundary surface kernel (zero-flux BCs) list, high-order variant: vz-direction (Tensor basis).
// At p=1 (the tensor hybrid, the only tensor basis with distinct lo/ho
// variants) this list holds the ho_ kernels; tensor p>1 has a single
// (high-order by design) variant shared with the plain tensor list.
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_ho_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_ho_surfvz_1x3v_tensor_p1, vlasov_boundary_surfvz_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL, NULL }, // 4
  { NULL, GKYL_HYB_2X3V(vlasov_boundary_ho_surfvz_2x3v_tensor_p1), NULL, NULL }, // 5
  // 3x kernels
  { NULL, GKYL_HYB_3X3V(vlasov_boundary_ho_surfvz_3x3v_tensor_p1), NULL, NULL }, // 6
};

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_vlasov_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (fabs(qInL[0]) < vlasov->skip_cell_thresh && 
      fabs(qInC[0]) < vlasov->skip_cell_thresh && 
      fabs(qInR[0]) < vlasov->skip_cell_thresh) {
    return 0.0; // Immediately return if phase-space density is below threshold.
  }

  if (dir < vlasov->cdim) {
    int idx_conf[GKYL_MAX_DIM], idx_confL[GKYL_MAX_DIM], idx_confR[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->cdim; ++i) {
      idx_conf[i] = idxC[i];
      idx_confL[i] = idxL[i];
      idx_confR[i] = idxR[i];
    }
    long cidx = gkyl_range_idx(&vlasov->conf_range, idx_conf);
    // Conf indices of the left/center/right cells: a conf-space surface
    // separates different conf cells, so the streaming kernel needs each cell's
    // (per-cell constant) position-map Jacobian.
    long cidxL = gkyl_range_idx(&vlasov->conf_range, idx_confL);
    long cidxR = gkyl_range_idx(&vlasov->conf_range, idx_confR);

    int idx_vel[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->pdim-vlasov->cdim; ++i) {
      idx_vel[i] = idxC[vlasov->cdim+i];
    }
    long vidx = gkyl_range_idx(&vlasov->vel_range, idx_vel);

    int idx_hamil[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->hamil_dim; ++i) {
      idx_hamil[i] = idxC[vlasov->hamil_offset+i];
    }
    long hidx = gkyl_range_idx(&vlasov->hamil_range, idx_hamil);

    if (vlasov->use_conf_flux_surf) {
      // Each cell owns the *lower* edge surface flux
      long pidxC = gkyl_range_idx(&vlasov->phase_range, idxC);
      long pidxR = gkyl_range_idx(&vlasov->phase_range, idxR);
      const double* conf_flux_surf_l = (const double*) gkyl_array_cfetch(vlasov->conf_flux_surf, pidxC);
      const double* conf_flux_surf_r = (const double*) gkyl_array_cfetch(vlasov->conf_flux_surf, pidxR);

      return vlasov->stream_surf_from_flux[dir](xcC, dxC,
        conf_flux_surf_l, conf_flux_surf_r, qRhsOut);
    }
    else {
      return vlasov->stream_surf[dir](xcC, dxC,
        vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidxL) : 0,
        vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidx) : 0,
        vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidxR) : 0,
        vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
        (const double*) gkyl_array_cfetch(vlasov->poisson_tensor_conf, cidx),
        (const double*) gkyl_array_cfetch(vlasov->hamil, hidx),
        qInL, qInC, qInR, qRhsOut);
    }
  }
  else {
    // Each cell owns the *lower* edge surface flux
    long pidxC = gkyl_range_idx(&vlasov->phase_range, idxC);
    long pidxR = gkyl_range_idx(&vlasov->phase_range, idxR);
    const double* vel_flux_surf_l = (const double*) gkyl_array_cfetch(vlasov->vel_flux_surf, pidxC);
    const double* vel_flux_surf_r = (const double*) gkyl_array_cfetch(vlasov->vel_flux_surf, pidxR);

    return vlasov->accel_surf[dir-vlasov->cdim](xcC, dxC,
      vel_flux_surf_l, vel_flux_surf_r, qRhsOut);
  }
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (fabs(qInEdge[0]) < vlasov->skip_cell_thresh && 
      fabs(qInSkin[0]) < vlasov->skip_cell_thresh) {
    return 0.0; // Immediately return if phase-space density is below threshold.
  }

  if (dir < vlasov->cdim) {
    int idx_conf[GKYL_MAX_DIM], idx_conf_edge[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->cdim; ++i) {
      idx_conf[i] = idxSkin[i];
      idx_conf_edge[i] = idxEdge[i];
    }
    long cidx = gkyl_range_idx(&vlasov->conf_range, idx_conf);
    // Conf index of the edge (ghost) cell for its position-map Jacobian.
    long cidx_edge = gkyl_range_idx(&vlasov->conf_range, idx_conf_edge);

    int idx_vel[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->pdim-vlasov->cdim; ++i) {
      idx_vel[i] = idxSkin[vlasov->cdim+i];
    }
    long vidx = gkyl_range_idx(&vlasov->vel_range, idx_vel);

    int idx_hamil[GKYL_MAX_DIM];
    for (int i=0; i<vlasov->hamil_dim; ++i) {
      idx_hamil[i] = idxSkin[vlasov->hamil_offset+i];
    }
    long hidx = gkyl_range_idx(&vlasov->hamil_range, idx_hamil);

    if (vlasov->use_conf_flux_surf) {
      const double* conf_flux_surf;
      if (edge == -1) {
        // Lower boundary needs to index edge value of flux
        long pidx = gkyl_range_idx(&vlasov->phase_range, idxEdge);
        conf_flux_surf = (const double*) gkyl_array_cfetch(vlasov->conf_flux_surf, pidx);
      }
      else {
        // Upper boundary needs to index skin value of flux
        long pidx = gkyl_range_idx(&vlasov->phase_range, idxSkin);
        conf_flux_surf = (const double*) gkyl_array_cfetch(vlasov->conf_flux_surf, pidx);      
      }

      return vlasov->stream_boundary_surf_from_flux[dir](xcSkin, dxSkin,
        edge, conf_flux_surf, qRhsOut);
    } 
    else {
      return vlasov->stream_boundary_surf[dir](xcSkin, dxSkin,
        vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidx_edge) : 0,
        vlasov->jacob_pos ? (const double*) gkyl_array_cfetch(vlasov->jacob_pos, cidx) : 0,
        vlasov->jacob_vel ? (const double*) gkyl_array_cfetch(vlasov->jacob_vel, vidx) : 0,
        (const double*) gkyl_array_cfetch(vlasov->poisson_tensor_conf, cidx),
        (const double*) gkyl_array_cfetch(vlasov->hamil, hidx), edge, qInEdge, qInSkin, qRhsOut);
    } 
  } 
  else if (dir >= vlasov->cdim) {
    const double* vel_flux_surf;
    if (edge == -1) {
      // Lower boundary needs to index edge value of flux
      long pidx = gkyl_range_idx(&vlasov->phase_range, idxEdge);
      vel_flux_surf = (const double*) gkyl_array_cfetch(vlasov->vel_flux_surf, pidx);
    }
    else {
      // Upper boundary needs to index skin value of flux
      long pidx = gkyl_range_idx(&vlasov->phase_range, idxSkin);
      vel_flux_surf = (const double*) gkyl_array_cfetch(vlasov->vel_flux_surf, pidx);      
    }

    return vlasov->accel_boundary_surf[dir-vlasov->cdim](xcSkin, dxSkin,
      edge, vel_flux_surf, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_diag(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return boundary_surf(eqn, dir, xcEdge, xcSkin, dxEdge, dxSkin, idxEdge, idxSkin, edge, qInEdge, qInSkin, qRhsOut);
}
