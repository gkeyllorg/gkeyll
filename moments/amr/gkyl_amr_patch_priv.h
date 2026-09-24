#pragma once

#include <gkyl_amr_block_priv.h>

// Definitions of private structs and APIs attached to these objects, for use in the patch AMR subsystem.

// Ranges indicating the extents of the skin and ghost regions for patch AMR, for the purpose of applying inter-patch boundary conditions.
struct skin_ghost_ranges_patch {
  struct gkyl_range lower_skin[1];
  struct gkyl_range lower_ghost[1];

  struct gkyl_range upper_skin[1];
  struct gkyl_range upper_ghost[1];
};

// Patch-structured data (including grid, range and boundary condition information) for the Euler equations.
struct euler_patch_data {
  struct gkyl_rect_grid grid;
  gkyl_fv_proj *fv_proj;

  struct gkyl_range ext_range;
  struct gkyl_range range;
  struct gkyl_array *fdup;
  struct gkyl_array *f[2];

  struct gkyl_wave_geom *geom;
  
  struct gkyl_wv_eqn *euler;
  gkyl_wave_prop *slvr[1];

  struct skin_ghost_ranges_patch skin_ghost;
  struct gkyl_array *bc_buffer;

  struct gkyl_wv_apply_bc *lower_bc[1];
  struct gkyl_wv_apply_bc *upper_bc[1];
};

// Job pool information context for updating patch-structured data for the Euler equations using threads.
struct euler_update_patch_ctx {
  const struct euler_patch_data *pdata;
  int pidx;
  int dir;
  double t_curr;
  double dt;
  struct gkyl_wave_prop_status stat;
};

/**
* Initialize the ranges indicating the extents of the skin and ghost regions for patch AMR, for the purposes of applying inter-patch boundary conditions.
*
* @param sgr Ranges for the skin and ghost regions.
* @param parent Ranges for the parent regions (of which the skin and ghost regions are subregions).
* @param ghost Number of ghost (and therefore skin) cells.
*/
void skin_ghost_ranges_init_patch(struct skin_ghost_ranges_patch* sgr, const struct gkyl_range* parent, const int* ghost);

/**
* Initialize patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the Euler equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the Euler equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void euler_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the Euler equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the Euler equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void euler_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the general relativistic Euler equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the general relativistic Euler equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void gr_euler_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the general relativistic Euler equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the general relativistic Euler equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void gr_euler_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the Euler mixture equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the Euler mixture equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void euler_mixture_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the Euler mixture equations.
*
* @param eqn Base equation object.
* @param pdata Patch-structured data for the Euler mixture equations.
* @param conn Topology/connectivity data for the patch hierarchy.
*/
void euler_mixture_nested_patch_bc_updaters_init(const struct gkyl_wv_eqn* eqn, struct euler_patch_data* pdata, const struct gkyl_block_connections* conn);

/**
* Release patch AMR updaters for both physical (outer-patch) and non-physical (inter-patch) boundary conditions for the Euler equations.
*
* @param pdata Patch-structured data for the Euler equations.
*/
void euler_patch_bc_updaters_release(struct euler_patch_data* pdata);

/**
* Apply both physical (outer-patch) and non-physical (inter-patch) patch AMR boundary conditions for the Euler equations.
*
* @param pdata Patch-structured data for the Euler equations.
* @param tm Simulation time at which the boundary conditions are applied.
* @param fld Output array.
*/
void euler_patch_bc_updaters_apply(const struct euler_patch_data* pdata, double tm, struct gkyl_array* fld);

/**
* Coarse-to-fine projection operator for patch-structured AMR, assuming a lower coarse patch and a lower fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Fine-to-coarse restriction operator for patch-structured AMR, assuming a lower fine patch and a lower coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Coarse-to-fine projection operator for patch-structured AMR, assuming a lower coarse patch and an upper fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Fine-to-coarse restriction operator for patch-structured AMR, assuming a lower fine patch and an upper coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

  /**
* Coarse-to-fine projection operator for patch-structured AMR, assuming an upper coarse patch and a lower fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Fine-to-coarse restriction operator for patch-structured AMR, assuming an upper fine patch and a lower coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Coarse-to-fine projection operator for patch-structured AMR, assuming an upper coarse patch and an upper fine patch.
*
* @param tbid Target (fine) patch ID.
* @param tdir Target (fine) patch direction.
* @param i Reference (coarse) patch ID.
* @param d Reference (coarse) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Fine-to-coarse restriction operator for patch-structured AMR, assuming an upper fine patch and an upper coarse patch.
*
* @param tbid Target (coarse) patch ID.
* @param tdir Target (coarse) patch direction.
* @param i Reference (fine) patch ID.
* @param d Reference (fine) patch direction.
* @param pdata Patch-structured data for the Euler equations.
* @param bc_buffer Buffer for applying boundary conditions.
* @param fld Output array.
*/
void patch_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct euler_patch_data pdata[],
  const struct gkyl_array* bc_buffer, struct gkyl_array* fld[]);

/**
* Synchronize all patches in the patch AMR hierarchy by applying all appropriate physical (outer-patch) and non-physical (inter-patch)
* boundary conditions for the Euler equations.
*
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch=structured data for the Euler equations.
* @param fld Output array.
*/
void euler_sync_patches(const struct gkyl_block_topo* ptopo, const struct euler_patch_data pdata[], struct gkyl_array* fld[]);

/**
* Write patch-structured AMR simulation data for the Euler equations onto disk.
*
* @param file_nm File name schema to use for the simulation output.
* @param pdata Patch-structured data for the Euler equations.
*/
void euler_patch_data_write(const char* file_nm, const struct euler_patch_data* pdata);

/**
* Calculate the maximum stable time-step for the patch-structured Euler equations.
*
* @param pdata Patch-structured data for the Euler equations.
* @return Maximum stable time-step.
*/
double euler_patch_data_max_dt(const struct euler_patch_data* pdata);

/**
* Update the patch-structured AMR simulation data for the Euler equations using the thread-based job pool.
*
* @param ctx Context to pass to the function.
*/
void euler_update_patch_job_func(void* ctx);

/**
* Update all patches in the patch AMR hierarchy by using the thread-based job pool for the Euler equations.
*
* @param job_pool Job pool for updating patch-structured data for the Euler equations using threads.
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch-structured data for the Euler equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step for the simulation.
* @return Status of the update (success and suggested time-step).
*/
struct gkyl_update_status euler_update_all_patches(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct euler_patch_data pdata[], double t_curr, double dt);

/**
* Initialize a new job in the thread-based job pool for updating the patch-structured AMR simulation data for the Euler equations.
*
* @param ctx Context to pass to the function.
*/
void euler_init_job_func_patch(void* ctx);

/**
* Take a single time-step across the entire patch AMR hierarchy for the Euler equations.
*
* @param job_pool Job pool for updating the patch-structured data for the Euler equations using threads.
* @param ptopo Topology/connectivity information for the entire patch hierarchy.
* @param pdata Patch-structured data for the Euler equations.
* @param t_curr Current simulation time.
* @param dt0 Initial guess for the maximum stable time-step.
* @param stats Simulation statistics (allowing for tracking of the number of failed time-steps).
* @return Status of the update (success, suggested time-step and actual time-step).
*/
struct gkyl_update_status euler_update_patch(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* ptopo,
  const struct euler_patch_data pdata[], double t_curr, double dt0, struct sim_stats* stats);

/**
* Write the complete simulation output for the entire patch AMR hierarchy for the Euler equations onto disk.
*
* @param fbase Base file name schema to use for the simulation output.
* @param num_patches Number of patches in the patch hierarchy.
* @param pdata Array of patch-structured data for the Euler equations.
*/
void euler_write_sol_patch(const char* fbase, int num_patches, const struct euler_patch_data pdata[]);

/**
* Calculate the maximum stable time-step across all patches in the patch AMR hierarchy for the Euler equations.
*
* @param num_patches Number of patches in the patch hierarchy.
* @param pdata Array of patch-structured data for the Euler equations.
* @return Maximum stable time-step.
*/
double euler_max_dt_patch(int num_patches, const struct euler_patch_data pdata[]);

/**
* Set up the topology/connectivity information for the patch AMR hierarchy for a mesh containing a single refinement patch.
*/
struct gkyl_block_topo* create_patch_topo();

/**
* Set up the topology/connectivity information for the patch AMR hierarchy for a mesh containing a doubly-nested refinement patch.
*/
struct gkyl_block_topo* create_nested_patch_topo();