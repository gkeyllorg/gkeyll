// Private header for use in MB gyrokinetic app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_rrobin_decomp.h>
#include <gkyl_multib_comm_conn.h>

// A multib_comm_conn send/recv pair used for communicating between blocks.
struct gkyl_mbcc_sr {
  struct gkyl_multib_comm_conn **recv, **send;
};

// top-level internal App
struct gkyl_gyrokinetic_multib_app {
  char name[128]; // name of app
  struct gkyl_comm *comm; // global communicator to use
  bool use_gpu; // Whether to use the GPU.
  
  // geometry and topology of all blocks in simulation
  struct gkyl_gk_block_geom *gk_block_geom;
  struct gkyl_block_topo *block_topo;
  
  double cfl_frac; // CFL fraction to use
  double bmag_ref; // Reference magnetic field
  int num_species; // number of species
  int num_neut_species; // number of neutral species

  bool update_field; // true if there solving Poisson equation
  struct gk_multib_field *field; // Field object.

  char species_name[GKYL_MAX_SPECIES][128]; // name of each species
  char neut_species_name[GKYL_MAX_SPECIES][128]; // name of each neutral species  

  struct gkyl_comm **block_comms; // list of block-communicators

  int num_local_blocks; // total number of blocks on current rank
  int *local_blocks; // local blocks IDs handled by current rank
  struct gkyl_gyrokinetic_app **singleb_apps; // App objects: one per local block

  const struct gkyl_rrobin_decomp *round_robin; // round-robin decomp
  struct gkyl_rect_decomp **decomp; // list of decomps (num_blocks)

  struct gkyl_mbcc_sr *mbcc_sync_conf; // Connections for conf-space sync.
  struct gkyl_mbcc_sr *mbcc_sync_charged; // Connections for charged species phase-space.
  struct gkyl_mbcc_sr *mbcc_sync_neut; // Connections for neut species phase-space.

  double tcurr; // current time
  
  struct gkyl_gyrokinetic_stat stat; // statistics

  gkyl_dynvec dts; // Record time step over time.
  bool is_first_dt_write_call; // flag for integrated moments dynvec written first time
  
  struct gkyl_msgpack_map_elem* io_meta_basic; // Metadata for I/O.
  int io_meta_basic_len; // Number of elements in io_meta.
};

// Meta-data for IO
struct gyrokinetic_multib_output_meta {
  int frame; // frame number
  double stime; // output time
  const char *app_name; // name of App
  const char *topo_file_name; // name of topology file
};

// field data
struct gk_multib_field {
  struct gkyl_gyrokinetic_multib_field info; // data for field
  enum gkyl_gkfield_id gkfield_id; // type of field
  int num_local_blocks; // total number of blocks on current rank
  int cdim; // number of configuration space dimensions
  bool half_domain; // For use in double null
                    // Whether to set BCs for simulation of lower half (Z<0)

  int num_blocks_below[GKYL_MAX_BLOCKS];
  int num_blocks_above[GKYL_MAX_BLOCKS];

  struct gkyl_array **phi_local;
  struct gkyl_array **rho_c_local;

  struct gkyl_array **rho_c_global_dg;

  struct gkyl_array **phi_global_smooth;
  struct gkyl_array **rho_c_global_smooth;

  //
  // Objects for parallel smoothing.
  //
  // Comm conn for sends/recvs in allgather.
  struct gkyl_multib_comm_conn **mbcc_allgatherz_send;
  struct gkyl_multib_comm_conn **mbcc_allgatherz_recv;
  struct gkyl_range **multibz_ranges; // Multib ranges.
  struct gkyl_range **multibz_ranges_ext; // Extended multib ranges.
  struct gkyl_range **block_subrangesz; // Ranges for copying from multib to local range.
  struct gkyl_range **parent_subrangesz; // Ranges for copying from multib to global range.
  // Parallel multib potential (DG and smoothed).
  struct gkyl_array **phi_multibz_dg;
  struct gkyl_array **phi_multibz_smooth;
  // Parallel multib charge density (DG and smoothed).
  struct gkyl_array **rho_c_multibz_dg;
  struct gkyl_array **rho_c_multibz_smooth;
  // Multib weights.
  struct gkyl_array **lhs_weight_multibz;
  struct gkyl_array **rhs_weight_multibz;
  struct gkyl_fem_parproj **fem_parproj; // FEM smoothing operator.
  
  //
  // Objects for perpendicular Poisson solve.
  //
  // Comm conn for sends/recvs in allgather.
  struct gkyl_multib_comm_conn **mbcc_allgather_perp_send;
  struct gkyl_multib_comm_conn **mbcc_allgather_perp_recv;
  struct gkyl_range **multib_perp_ranges; // Multib ranges.
  struct gkyl_range **multib_perp_ranges_ext; // Extended multib ranges.
  struct gkyl_range **block_subranges_perp; // Ranges for copying from multib to local range.
  struct gkyl_range **parent_subranges_perp; // Ranges for copying from multib to global range.
  // Potential and charge density on perpendicular multib range.
  struct gkyl_array **phi_multib_perp;
  struct gkyl_array **rho_c_multib_perp;
  struct gkyl_array **epsilon_multib_perp; // Multib polarization weight.
  struct gkyl_fem_poisson_perp **fem_poisson; // Perpendicular Poisson solver.

  // Functions to advance MB field
  // Pointer to function to calculate the potential.
  void (*rhs_func)(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[]);
};

/** Time stepping API */

/**
 * Compute the gyrokinetic fields.
 *
 * @param app Multiblock gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param fin Array of distribution functions (one for each species) .
 * @param bflux Moments of the boundary fluxes (for all blocks, all species, and all boundaries).
 */
void gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[]);

/**
 * Compute the gyrokinetic fields and apply boundary conditions.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each charged species).
 * @param bflux Moments of the boundary fluxes (for all blocks, all species, and all boundaries).
 * @param distf_neut Array of distribution functions (for each neutral species).
 */
void gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[]);

/**
 * Take time-step using the RK3 method. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_multib_update_ssp_rk3(struct gkyl_gyrokinetic_multib_app* app, double dt0);

/** Field API */

/** Initialize multib field object
 * @param mbinp App inputs.
 * @param mbapp Gyrokinetic multib app.
 * return new multib field object
 */
struct gk_multib_field* gk_multib_field_new(const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *mbapp);


/** Compute the electrostatic potential
 * @param mbapp Gyrokinetic multib app.
 * @param mbf Multib field object.
 * @param fin Distribution function (for all local blocks).
 * @param bflux Moments of the boundary fluxes (for all blocks, all species, and all boundaries).
*/
void gk_multib_field_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[]);

/** Releas the resources for the multib field object
 * @param mbf Multib field object.
*/
void gk_multib_field_release(struct gk_multib_field *mbf);
