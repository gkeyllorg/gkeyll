#pragma once

#include <gkyl_amr_block_priv.h>

// Definitions of private structs and APIs attached to these objects, for use in the coupled block AMR subsystem.

// Block-structured data (including grid, range and boundary condition information) for the coupled five-moment equations.
struct five_moment_block_data {
  struct gkyl_rect_grid grid;

  gkyl_fv_proj *fv_proj_elc;
  gkyl_fv_proj *fv_proj_ion;
  gkyl_fv_proj *fv_proj_maxwell;

  struct gkyl_range ext_range;
  struct gkyl_range range;

  struct gkyl_array *fdup_elc;
  struct gkyl_array *fdup_ion;
  struct gkyl_array *fdup_maxwell;

  struct gkyl_array *f_elc[3];
  struct gkyl_array *f_ion[3];
  struct gkyl_array *f_maxwell[3];

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

  gkyl_wave_prop *slvr_elc[2];
  gkyl_wave_prop *slvr_ion[2];
  gkyl_wave_prop *slvr_maxwell[2];

  struct skin_ghost_ranges_block skin_ghost;

  struct gkyl_array *bc_buffer_elc;
  struct gkyl_array *bc_buffer_ion;
  struct gkyl_array *bc_buffer_maxwell;

  struct gkyl_wv_apply_bc *lower_bc_elc[2];
  struct gkyl_wv_apply_bc *lower_bc_ion[2];
  struct gkyl_wv_apply_bc *lower_bc_maxwell[2];

  struct gkyl_wv_apply_bc *upper_bc_elc[2];
  struct gkyl_wv_apply_bc *upper_bc_ion[2];
  struct gkyl_wv_apply_bc *upper_bc_maxwell[2];

  gkyl_moment_em_coupling *src_slvr;

  bool copy_x;
  bool copy_y;

  bool wall_x;
  bool wall_y;
};

// Job pool information context for updating block-structured data for the coupled five-moment equations using threads.
struct five_moment_update_block_ctx {
  const struct five_moment_block_data *bdata;
  int bidx;
  int dir;
  double t_curr;
  double dt;

  struct gkyl_wave_prop_status stat_elc;
  struct gkyl_wave_prop_status stat_ion;
  struct gkyl_wave_prop_status stat_maxwell;

  int nstrang;
};

// Context for copying job pool information for the coupled five-moment equations.
struct five_moment_copy_job_ctx {
  int bidx;

  const struct gkyl_array *inp_elc;
  const struct gkyl_array *inp_ion;
  const struct gkyl_array *inp_maxwell;

  struct gkyl_array *out_elc;
  struct gkyl_array *out_ion;
  struct gkyl_array *out_maxwell;
};

/**
* Boundary condition function for applying wall boundary conditions for the coupled five-moment equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void five_moment_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying wall boundary conditions for the coupled ten-moment equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void ten_moment_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying wall boundary conditions for the Maxwell equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void maxwell_wall_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying copy boundary conditions for the coupled five-moment equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply copy boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void five_moment_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying copy boundary conditions for the coupled ten-moment equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply copy boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void ten_moment_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying copy boundary conditions for the Maxwell equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply tranmissive boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
void maxwell_copy_bc(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* GKYL_RESTRICT skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Initialize block AMR updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @param conn Topology/connectivity data for the block hierarchy.
*/
void five_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested block AMR updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @param conn Topology/connectivity data for the block hierarchy.
*/
void five_moment_nested_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn);

/**
* Initialize block AMR updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled ten-moment equations.
*
* @param bdata Block-structured data for the coupled ten-moment equations.
* @param conn Topology/connectivity data for the block hierarchy.
*/
void ten_moment_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn);

/**
* Initialize nested block AMR updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled ten-moment equations.
*
* @param bdata Block-structured data for the coupled ten-moment equations.
* @param conn Topology/connectivity data for the block hierarchy.
*/
void ten_moment_nested_block_bc_updaters_init(struct five_moment_block_data* bdata, const struct gkyl_block_connections* conn);

/**
* Release block AMR updaters for both physical (outer-block) and non-physical (inter-block) boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
*/
void five_moment_block_bc_updaters_release(struct five_moment_block_data* bdata);

/**
* Apply both physical (outer-block) and non-physical (inter-block) block AMR boundary conditions for the coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @param tm Simulation time at which the boundary conditions are applied.
* @param fld_elc Output array (electrons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void five_moment_block_bc_updaters_apply(const struct five_moment_block_data* bdata, double tm,
  struct gkyl_array* fld_elc, struct gkyl_array *fld_ion, struct gkyl_array* fld_maxwell);

/**
* Coarse-to-fine projection operator for coupled, block-structured AMR, assuming a lower coarse block and a lower fine block.
*
* @param tbid Target (fine) block ID.
* @param tdir Target (fine) block direction.
* @param i Reference (coarse) block ID.
* @param d Reference (coarse) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_ll_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, block-structured AMR, assuming a lower fine block and a lower coarse block.
*
* @param tbid Target (coarse) block ID.
* @param tdir Target (coarse) block direction.
* @param i Reference (fine) block ID.
* @param d Reference (fine) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_ll_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, block-structured AMR, assuming a lower coarse block and an upper fine block.
*
* @param tbid Target (fine) block ID.
* @param tdir Target (fine) block direction.
* @param i Reference (coarse) block ID.
* @param d Reference (coarse) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_lu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, block-structured AMR, assuming a lower fine block and an upper coarse block.
*
* @param tbid Target (coarse) block ID.
* @param tdir Target (coarse) block direction.
* @param i Reference (fine) block ID.
* @param d Reference (fine) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_lu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
 const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
 struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, block-structured AMR, assuming an upper coarse block and a lower fine block.
*
* @param tbid Target (fine) block ID.
* @param tdir Target (fine) block direction.
* @param i Reference (coarse) block ID.
* @param d Reference (coarse) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_ul_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, block-structured AMR, assuming an upper fine block and a lower coarse block.
*
* @param tbid Target (coarse) block ID.
* @param tdir Target (coarse) block direction.
* @param i Reference (fine) block ID.
* @param d Reference (fine) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_ul_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Coarse-to-fine projection operator for coupled, block-structured AMR, assuming an upper coarse block and an upper fine block.
*
* @param tbid Target (fine) block ID.
* @param tdir Target (fine) block direction.
* @param i Reference (coarse) block ID.
* @param d Reference (coarse) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_uu_projection_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Fine-to-coarse restriction operator for coupled, block-structured AMR, assuming an upper fine block and an upper coarse block.
*
* @param tbid Target (coarse) block ID.
* @param tdir Target (coarse) block direction.
* @param i Reference (fine) block ID.
* @param d Reference (fine) block direction.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param bc_buffer_elc Buffer for applying electron boundary conditions.
* @param bc_buffer_ion Buffer for applying ion boundary conditions.
* @param bc_buffer_maxwell Buffer for applying Maxwell field boundary conditions.
* @param fld_elc Output array (electons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void block_coupled_uu_restriction_op(const int tbid, const int tdir, const int i, const int d, const struct five_moment_block_data bdata[],
  const struct gkyl_array* bc_buffer_elc, const struct gkyl_array* bc_buffer_ion, const struct gkyl_array* bc_buffer_maxwell,
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Synchronize all blocks in the block AMR hierarchy by applying all appropriate physical (outer-block) and non-physical (inter-block)
* boundary conditions for the coupled five-moment equations.
*
* @param btopo Topology/connectivity information for the block hierarchy.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param fld_elc Output array (electrons).
* @param fld_ion Output array (ions).
* @param fld_maxwell Output array (Maxwell field).
*/
void five_moment_sync_blocks(const struct gkyl_block_topo* btopo, const struct five_moment_block_data bdata[],
  struct gkyl_array* fld_elc[], struct gkyl_array* fld_ion[], struct gkyl_array* fld_maxwell[]);

/**
* Write block-structured AMR simulation data for the coupled five-moment equations onto disk.
*
* @param file_nm_elc File name schema to use for the electron simulation output.
* @param file_nm_ion File name schema to use for the ion simulation output.
* @param file_nm_maxwell File name schema to use for the Maxwell field simulation output.
* @param bdata Block-structured data for the coupled five-moment equations.
*/
void five_moment_block_data_write(const char* file_nm_elc, const char* file_nm_ion, const char* file_nm_maxwell, const struct five_moment_block_data* bdata);

/**
* Calculate the maximum stable time-step for the block-structured, coupled five-moment equations.
*
* @param bdata Block-structured data for the coupled five-moment equations.
* @return Maximum stable time-step.
*/
double five_moment_block_data_max_dt(const struct five_moment_block_data* bdata);

/**
* Update the block-structured AMR simulation data for the coupled five-moment equations using the thread-based job pool.
*
* @param ctx Context to pass to the function.
*/
void five_moment_update_block_job_func(void* ctx);

/**
* Update the source terms of the block-structured AMR simulation data for the coupled five-moment equations using the thread-based job pool.
*
* @param ctx Context to pass to the function.
*/
void five_moment_update_block_job_func_source(void* ctx);

/**
* Update all blocks in the block AMR hierarchy by using the thread-based job pool for the coupled five-moment equations.
*
* @param job_pool Job pool for updating block-structured data for the coupled five-moment equations using threads.
* @param btopo Topology/connectivity information for the entire block hierarchy.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step for the simulation.
* @return Status of the update (success and suggested time-step).
*/
struct gkyl_update_status five_moment_update_all_blocks(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt);

/**
* Update the source terms for all blocks in the block AMR hierarchy by using the thread-based job pool for the coupled five-moment equations.
*
* @param job_pool Job pool for updating block-structured data for the coupled five-moment equations using threads.
* @param btopo Topology/connectivity information for the entire block hierarchy.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step for the simulation.
* @param nstrang Iteration number in the Strang splitting.
*/
void five_moment_update_all_blocks_source(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt, int nstrang);

/**
* Initialize a new job in the thread-based job pool for updating the block-structured AMR simulation data for the coupled five-moment equations.
*
* @param ctx Context to pass to the function.
*/
void five_moment_init_job_func_block(void* ctx);

/**
* Copy an existing job between two thread-based job pools for updating the coupled five-moment equations.
*
* @param ctx Context to pass to the function.
*/
void five_moment_copy_job_func(void* ctx);

/**
* Take a single time-step across the entire block AMR hierarchy for the coupled five-moment equations.
*
* @param job_pool Job pool for updating block-structured data for the coupled five-moment equations using threads.
* @param btopo Topology/connectivity information for the entire block hierarchy.
* @param bdata Block-structured data for the coupled five-moment equations.
* @param t_curr Current simulation time.
* @param dt0 Initial guess for the maximum stable time-step.
* @param stats Simulation statistics (allowing for tracking of the number of failed time-steps).
* @return Status of the update (success, suggested time-step and actual time-step).
*/
struct gkyl_update_status five_moment_update_block(const struct gkyl_job_pool* job_pool, const struct gkyl_block_topo* btopo,
  const struct five_moment_block_data bdata[], double t_curr, double dt0, struct sim_stats* stats);

/**
* Write the complete simulation output for the entire block AMR hierarchy for the coupled five-moment equations onto disk.
*
* @param fbase Base file name schema to use for the simulation output.
* @param num_blocks Number of blocks in the block hierarchy.
* @param bdata Array of block-structured data for the coupled five-moment equations.
*/
void five_moment_write_sol_block(const char* fbase, int num_blocks, const struct five_moment_block_data bdata[]);

/**
* Calculate the maximum stable time-step across all blocks in the block AMR hierarchy for the coupled five-moment equations.
*
* @param num_blocks Number of blocks in the block hierarchy.
* @param bdata Array of block-structured data for the coupled five-moment equations.
* @return Maximum stable time-step
*/
double five_moment_max_dt_block(int num_blocks, const struct five_moment_block_data bdata[]);