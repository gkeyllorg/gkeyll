#pragma once

#include <gkyl_gk_block_geom.h>
#include <gkyl_gyrokinetic.h>

typedef struct gkyl_gyrokinetic_multib_app gkyl_gyrokinetic_multib_app;

// Species input per block
struct gkyl_gyrokinetic_multib_species_pb {
  int block_id; // block ID

  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;

  double polarization_density;

  // Source to include.
  struct gkyl_gyrokinetic_source source;
};

// Species input
struct gkyl_gyrokinetic_multib_species {
  char name[128]; // Species name.

  double charge, mass; // Charge and mass.
  int vdim; // Velocity space dimensions.
  double lower[GKYL_MAX_VDIM], upper[GKYL_MAX_VDIM]; // Lower, upper bounds of velocity-space.
  int cells[GKYL_MAX_VDIM]; // Velocity-space cells.
  struct gkyl_mapc2p_inp mapc2p; // Velocity mapping.

  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.
  int num_integrated_diag_moments; // Number of integrated diagnostic moments.
  enum gkyl_distribution_moments integrated_diag_moments[12]; // List of integrated diagnostic moments.
  bool time_rate_diagnostics; // Whether to ouput df/dt diagnostics.

  struct gkyl_phase_diagnostics_inp boundary_flux_diagnostics;

  // Collisionless terms.
  struct gkyl_gyrokinetic_collisionless collisionless;

  // Collisions to include.
  struct gkyl_gyrokinetic_collisions collisions;

  // Anomalous diffusion.
  struct gkyl_gyrokinetic_anomalous_diffusion anomalous_diffusion;

  // Radiation to include.
  struct gkyl_gyrokinetic_radiation radiation;

  // Reactions between plasma species to include.
  struct gkyl_gyrokinetic_react react;
  // Reactions with neutral species to include.
  struct gkyl_gyrokinetic_react react_neut;

  bool duplicate_across_blocks; // Set to true if all blocks are identical.
  // Species inputs per-block: only one is needed if duplicate_across_blocks = true.
  const struct gkyl_gyrokinetic_multib_species_pb *blocks;

  // Physical boundary conditions.
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_bc *bcs;
};

// Neutral species input per block
struct gkyl_gyrokinetic_multib_neut_species_pb {
  int block_id; // block ID

  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;

  // Source to include.
  struct gkyl_gyrokinetic_source source;
};

// Neutral species input
struct gkyl_gyrokinetic_multib_neut_species {
  char name[128]; // Species name.

  double mass; // Mass.
  int vdim; // Velocity space dimensions.
  double lower[GKYL_MAX_VDIM], upper[GKYL_MAX_VDIM]; // Lower, upper bounds of velocity-space.
  int cells[GKYL_MAX_VDIM]; // Velocity-space cells.

  struct gkyl_mapc2p_inp mapc2p;

  bool is_static; // Set to true if neutral species does not change in time.

  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.

  // Reactions with plasma species to include.
  struct gkyl_gyrokinetic_react react_neut;

  bool duplicate_across_blocks; // set to true if all blocks are identical  
  // species inputs per-block: only one is needed if duplicate_across_blocks = true
  const struct gkyl_gyrokinetic_multib_neut_species_pb *blocks;

  // Physical boundary conditions
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_bc *bcs;
};

// Field input per block 
struct gkyl_gyrokinetic_multib_field_pb {
  int block_id; // block ID

  double polarization_bmag; // B  to use in linearized polarization weight
  double kperpSq; // kperp^2 parameter for 1D field equations

  bool time_rate_diagnostics; // Writes the time rate of change of field energy.

};

// Field input 
struct gkyl_gyrokinetic_multib_field {
  enum gkyl_gkfield_id gkfield_id;
  double kperpSq; // kperp^2 parameter for 1D field equations

  // parameters for adiabatic electrons simulations
  double electron_mass, electron_charge, electron_density, electron_temp;

  //struct gkyl_poisson_bc poisson_bcs;

  bool duplicate_across_blocks; // set to true if all blocks are identical  
  // field inputs per-block: only one is needed if duplicate_across_blocks = true
  bool half_domain; // For use in double null simulations. 
                    // If true core BCs will be set for simulation of lower half of domain (Z < 0)
  const struct gkyl_gyrokinetic_multib_field_pb *blocks;

  // Physical boundary conditions
  int num_physical_bcs;
  const struct gkyl_gyrokinetic_bc *bcs;

  bool time_rate_diagnostics; // Writes the time rate of change of field energy.
  
  struct gkyl_poisson_bias_line_list *bias_line_list; // Biased lines constraining the solution.
};

// Top-level app parameters: this
struct gkyl_gyrokinetic_multib {
  char name[128]; // Name of app.

  int cdim; // Configuration space dimensions.
  int poly_order; // Polynomial order.
  enum gkyl_basis_type basis_type; // Type of basis functions to use.
  bool use_gpu; // Flag to indicate if solver should use GPUs

  // Geometry and topology of all blocks in simulation.
  struct gkyl_gk_block_geom *gk_block_geom;

  double cfl_frac; // CFL fraction to use (default 1.0)
  double cfl_frac_omegaH; // CFL fraction to use for omegaH (default 1.7)
  
  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[3]; // List of periodic directions.

  int num_species; // Number of species.
  // Species inputs.
  struct gkyl_gyrokinetic_multib_species species[GKYL_MAX_SPECIES];

  int num_neut_species; // Number of neutral species.
  // Neutral species inputs.
  struct gkyl_gyrokinetic_multib_neut_species neut_species[GKYL_MAX_SPECIES];

  bool skip_field; // Skip field update -> phi = 0 for all time.
  // Field inputs.
  struct gkyl_gyrokinetic_multib_field field;

  // EIRENE inputs
  struct gkyl_gyrokinetic_eirene eirene;

  // Communicator to use.  
  struct gkyl_comm *comm;  

  struct gkyl_gyrokinetic_metadata_inp metadata; // Optional metadata for output files.
};

/**
 * Construct a new gk multi-block app.
 *
 * @param mbinp App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New multi-block gk app object.
 */
gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp);

/**
 * Construct a new gk multi-block app (geom only).
 *
 * @param mbinp App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New multi-block gk app object.
 */
gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp);

/**
 * Initialize species by projecting initial conditions on
 * basis functions, and the field by solving the field equations.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0);

/**
 * Initialize neutral species by projecting initial conditions on basis
 * functions. Neutral species index (sidx) is the same index used to specify
 * the neutral species in the gkyl_gk object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname);

/**
 * Initialize the gyrokinetic app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...);

/**
 * Write block topology to file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app);

/**
 * Write geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species.
 *
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments of the boundary fluxes for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments of the boundary fluxes for charged
 * species to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 */
void gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species source.
 *
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species source.
 *
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write neutral species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write LBO collisional moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write BGK cross moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_multib_app_write_species_bgk_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write radiation drag coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Write radiation emissivity of each species that species sidx collides with
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments of the radiation model.
 *
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm);

/**
 * Write integrated moments of radiation rhs for radiating species 
 * 
 * @param app App object.
 * @param sidx Index of species from which to write radiation.
 */
void gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_iz_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_iz_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_recomb_react(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_recomb_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write cx react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_cx_react_neut(gkyl_gyrokinetic_multib_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for all species (including sources) to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app *app, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for all species (including sources).
 *
 * @param app App object.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm);

/**
 * Write integrated diagnostic moments for all species (including sources)
 * to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app);

/**
 * Write phase space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write configuration space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write both conf and phase-space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app);

/**
 * Print timing of solver components to iostream.
 *
 * @param app App object.
 * @param iostream Where to write timers to (e.g. stdout, stderr);
 */
void
gkyl_gyrokinetic_multib_app_print_timings(gkyl_gyrokinetic_multib_app* app, FILE *iostream);

/**
 * Record the time step (in private dynvector).
 *
 * @param app App object.
 * @param tm Time stamp.
 * @param dt Time step to record (e.g. provided by app's status object).
 */
void
gkyl_gyrokinetic_multib_app_save_dt(gkyl_gyrokinetic_multib_app* app, double tm, double dt);

/**
 * Write the time step over time.
 *
 * @param app App object.
 */
void
gkyl_gyrokinetic_multib_app_write_dt(gkyl_gyrokinetic_multib_app* app);

/**
 * Read geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_multib_app_read_geometry(gkyl_gyrokinetic_multib_app *app);

/**
 * Write output to console: this is mainly for diagnostic messages the
 * driver code wants to write to console. It accounts for parallel
 * output by not messing up the console with messages from each rank.
 *
 * @param app App object
 * @param fp File pointer for open file for output
 * @param fmt Format string for console output
 * @param argp Objects to write
 */
void gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app);

/**
 * Free gk app.
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* app);

/**
 * Free gk app (geom only).
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_multib_app_release_geom(gkyl_gyrokinetic_multib_app* app);
