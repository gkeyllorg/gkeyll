#pragma once

#include <gkyl_amr_block_priv.h>
#include <gkyl_amr_patch_priv.h>

// Definitions of private structs and APIs attached to these objects, for use in the coupled patch AMR subsystem.

// Patch-structured data (including grid, range and boundary condition information) for the coupled five-moment equations.
struct five_moment_patch_data {
  struct gkyl_rect_grid grid;

  gkyl_fv_proj *fv_proj_elc;
  gkyl_fv_proj *fv_proj_ion;
  gkyl_fv_proj *fv_proj_maxwell;

  struct gkyl_range ext_range;
  struct gkyl_range range;

  struct gkyl_array *fdup_elc;
  struct gkyl_array *fdup_ion;
  struct gkyl_array *fdup_maxwell;

  struct gkyl_array *f_elc[2];
  struct gkyl_array *f_ion[2];
  struct gkyl_array *f_maxwell[2];

  struct gkyl_array *app_accel_elc;
  struct gkyl_array *app_accel_ion;

  struct gkyl_array *rhs_source_elc;
  struct gkyl_array *rhs_source_ion;

  struct gkyl_array *nT_source_elc;
  struct gkyl_array *nT_source_ion;

  struct gkyl_array *app_current;
  struct gkyl_array *ext_em;

  struct gkyl_wave_geom *geom;

  struct gkyl_wv_eqn *euler_elc;
  struct gkyl_wv_eqn *euler_ion;
  struct gkyl_wv_eqn *maxwell;

  gkyl_wave_prop *slvr_elc[1];
  gkyl_wave_prop *slvr_ion[1];
  gkyl_wave_prop *slvr_maxwell[1];

  struct skin_ghost_ranges_patch skin_ghost;

  struct gkyl_array *bc_buffer_elc;
  struct gkyl_array *bc_buffer_ion;
  struct gkyl_array *bc_buffer_maxwell;

  struct gkyl_wv_apply_bc *lower_bc_elc[1];
  struct gkyl_wv_apply_bc *lower_bc_ion[1];
  struct gkyl_wv_apply_bc *lower_bc_maxwell[1];

  struct gkyl_wv_apply_bc *upper_bc_elc[1];
  struct gkyl_wv_apply_bc *upper_bc_ion[1];
  struct gkyl_wv_apply_bc *upper_bc_maxwell[1];

  gkyl_moment_em_coupling *src_slvr;
};

// Job pool information context for updating patch-structured data for the coupled five-moment equations using threads.
struct five_moment_update_patch_ctx {
  const struct five_moment_patch_data *pdata;
  int pidx;
  int dir;
  double t_curr;
  double dt;

  struct gkyl_wave_prop_status stat_elc;
  struct gkyl_wave_prop_status stat_ion;
  struct gkyl_wave_prop_status stat_maxwell;

  int nstrang;
};

/**
* Initialize patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the coupled five-moment equations.
*
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void five_moment_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the coupled five-moment equations.
*
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void five_moment_nested_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the coupled ten-moment equations.
*
* @param pdata Patch-structured data for the coupled ten-moment equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void ten_moment_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the coupled ten-moment equations.
*
* @param pdata Patch-structured data for the coupled ten-moment equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void ten_moment_nested_patch_bc_updaters_init(struct five_moment_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Release patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the coupled five-moment equations.
*
* @param pdata Patch-structured data for the coupled five-moment equations.
*/
void five_moment_patch_bc_updaters_release(struct five_moment_patch_data* pdata);

/**
* Apply both physical (outer-patch) and non-physical (inter-patch) patch AMR boundary conditions for the coupled five-moment equations.
*
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param tm Simulation time at which the boundary conditions are applied.
* @param fld_elc Output array (electrons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void five_moment_patch_bc_updaters_apply(const struct five_moment_patch_data* pdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array *fld_ion, struct gkyl_array* fld_maxwell);

/**
* Coarse-to-fine projection operator for coupled, patch-structured AMR, assuming a lower coarse patch and a lower fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, patch-structured AMR, assuming a lower fine patch and a lower coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, patch-structured AMR, assuming a lower coarse patch and an upper fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, patch-structured AMR, assuming a lower fine patch and an upper coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
 const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
 struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, patch-structured AMR, assuming an upper coarse patch and a lower fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, patch-structured AMR, assuming an upper fine patch and a lower coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, patch-structured AMR, assuming an upper coarse patch and an upper fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, patch-structured AMR, assuming an upper fine patch and an upper coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void patch_coupled_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_patch_data pdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Synchronize all patches in the patch AMR hierarchy by applying all appropriate physical (outer-patch) and non-physical (inter-patch)
* boundary conditions for the coupled five-moment equations.
*
* @param ptopo Topology/connectivity information for the patch hierarchy.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param fld_elc Output array (electrons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void five_moment_sync_patches(const struct gkyl_block_topo* ptopo, const struct five_moment_patch_data pdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Write patch-structured AMR simulation data for the coupled five-moment equations onto disk.
*
* @param file_nm_elc File name schema to use for the electron simulation output.
* @param file_nm_ion File name schema to use for the ion simulation output.
* @param file_nm_maxwell File name schema to use for the Maxwell field simulation output.
* @param pdata Patch-structured data for the coupled five-moment equations.
*/
void five_moment_patch_data_write(const char* file_nm_elc, const char* file_nm_ion, const char* file_nm_maxwell, const struct five_moment_patch_data* pdata);

/**
* Calculate the maximum stable time-step for the patch-structured, coupled five-moment equations.
*
* @param pdata Patch-structured data for the coupled five-moment equations.
* @return Maximum stable time-step.
*/
double five_moment_patch_data_max_dt(const struct five_moment_patch_data* pdata);

/**
* Update the patch-structured AMR simulation data for the coupled five-moment equations using the thread-based job pool.
*
* @param ctx Context to pass to the function.
*/
void five_moment_update_patch_job_func(void* ctx);

/**
* Update the source terms of the patch-structured AMR simulation data for the coupled five-moment equations using the thread-based job pool.
*
* @param ctx Context to pass to the function.
*/
void five_moment_update_patch_job_func_source(void* ctx);

/**
* Update all patches in the patch AMR hierarchy by using the thread-based job pool for the coupled five-moment equations.
*
* @param job_pool Job pool for updating patch-structured data for the coupled five-moment equations using threads.
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step for the simulation.
* @return Status of the update (success and suggested time-step).
*/
struct gkyl_update_status five_moment_update_all_patches(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt);

/**
* Update the source terms for all patches in the patch AMR hierarchy by using the thread-based job pool for the coupled five-moment equations.
*
* @param job_pool Job pool for updating patch-structured data for the coupled five-moment equations using threads.
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step for the simulation.
* @param nstrang Iteration number in the Strang splitting.
*/
void five_moment_update_all_patches_source(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt, int nstrang);

/**
* Initialize a new job in the thread-based job pool for updating the patch-structured AMR simulation data for the coupled five-moment equations.
*
* @param ctx Context to pass to the function.
*/
void five_moment_init_job_func_patch(void* ctx);

/**
* Take a single time-step across the entire patch AMR hierarchy for the coupled five-moment equations.
*
* @param job_pool Job pool for updating patch-structured data for the coupled five-moment equations using threads.
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt0 Initial guess for the maximum stable time-step.
* @param stats Simulation statistics (allowing for tracking of the number of failed time-steps).
* @return Status of the update (success, suggested time-step and actual time-step).
*/
struct gkyl_update_status five_moment_update_patch(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct five_moment_patch_data pdata[], double t_curr, double dt0, struct sim_stats* stats);

/**
* Write the complete simulation output for the entire patch AMR hierarchy for the coupled five-moment equations onto disk.
*
* @param fbase Base file name schema to use for the simulation output.
* @param num_patches Number of patches in the patch hierarchy.
* @param pdata Array of patch-structured data for the coupled five-moment equations.
*/
void five_moment_write_sol_patch(const char* fbase, int num_patches, const struct five_moment_patch_data pdata[]);

/**
* Calculate the maximum stable time-step across all patches in the patch AMR hierarchy for the coupled five-moment equations.
*
* @param num_patches Number of patches in the patch hierarchy.
* @param pdata Array of patch-structured data for the coupled five-moment equations.
* @return Maximum stable time-step
*/
double five_moment_max_dt_patch(int num_patch, const struct five_moment_patch_data pdata[]);