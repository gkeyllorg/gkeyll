#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>

#include <stdbool.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

// Parameters for the kinetic block of a PKPM species: the grid of the normalized velocity
// w = v_par/xi, the initial condition of the distributions, collisions, sources and the kinetic
// BCs. Declared as the 'kinetic' block of struct gkyl_pkpm_species. The members it shares with
// struct gkyl_vlasov_kinetic_species have the same names and types.
struct gkyl_pkpm_kinetic_species {
  double lower[3], upper[3]; // Lower, upper bounds of the NORMALIZED velocity w = v_par/xi.
  int cells[3]; // Velocity-space cells.

  // Initial condition using the projection routine: one GKYL_PROJ_FUNC projection of a function
  // returning the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)] of the normalized
  // distributions. F is a probability density in w (M0 = 1, M1 = 0) and G its
  // perpendicular-energy partner (M0(G) = T_perp/T) with M2(F) + 2 M0(G) = 3; e.g. a
  // bi-Maxwellian is F = exp(-w^2/(2 r))/sqrt(2 pi r), G = (T_perp/T) F with r = T_par/T,
  // T = (T_par + 2 T_perp)/3. The first-harmonic blocks [M0_perp, N_perp] must be perpendicular
  // to b, and int M0_perp dw must vanish; a gyrotropic start returns zeros for them.
  int num_init; // Number of initial condition functions (1).
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  // Collisions to include (GKYL_LBO_COLLISIONS or GKYL_BGK_COLLISIONS). Without self_nu the
  // self-collision frequency is the Spitzer one computed from the local moments, and so are
  // the cross-collision frequencies without cross_nu,
  //   nu_sr = nu_frac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2),
  // with the constant part alpha_E of Morse's energy-exchange coefficient at the reference
  // parameters (nu_ss = nu_frac alpha_E n_s / (2 vtsq_s)^(3/2) for the self frequency).
  // is_implicit (BGK): one backward-Euler step in each SSP-RK3 stage, over the time weight of
  // the stage, instead of the explicit update (with its CFL limit); not available on the GPU.
  // write_coll_diagnostics and fixed_temp_relax are not available (source_bgk is the
  // temperature source).
  struct gkyl_vlasov_collisions collisions;

  // Source to include: GKYL_PROJ_SOURCE or GKYL_BFLUX_SOURCE with one GKYL_PROJ_VLASOV_LTE
  // projection, the profile [density, isotropic temperature, drift] injected into the
  // normalized kinetic equations (through the configuration pass, so the closures keep the
  // constraints exact) and into the fluid state [rho, M, E]. With GKYL_PROJ_SOURCE the profile
  // is added as given; with GKYL_BFLUX_SOURCE it is scaled every stage by the particle loss
  // rate of `source_species` through its absorbing walls divided by `source_length` (so a
  // density profile integrating to source_length replaces the losses exactly). The drift
  // enters the parallel Maxwellian through its component along b; a perpendicular drift
  // relative to the local flow is not carried by G. write_source writes the fluid image of the
  // profile at frame 0.
  struct gkyl_vlasov_source source;

  // Temperature source to include, S = nu_Q(x) (f_M(n, u, T_Q(x)) - f) (see struct
  // gkyl_vlasov_source_bgk). In the normalized variables the fluid thermal energy relaxes to
  // (3/2) n T_Q at the rate nu_Q and the distributions relax to the unit Maxwellian at the rate
  // nu_Q T_Q/T, which relaxes every pressure component at the rate nu_Q of the source. Terms of
  // order (T_Q/T - 1) times the departure of the distributions from a Maxwellian are dropped:
  // moments of degree three decay at nu_Q T_Q/T instead of nu_Q (3 T_Q/T - 1)/2. The explicit
  // time step is limited by nu_Q and by nu_Q T_Q/T.
  struct gkyl_vlasov_source_bgk source_bgk;

  // Gyrophase-mixing closure of the first harmonic: rate, in units of the species' local |Omega|, at
  // which agyrotropy beyond what one harmonic can represent relaxes back to that bound (0: the
  // default 10; a negative value switches the closure off). The overshoot above the bound is about
  // (drive rate/|Omega|)/gyro_mix; a drive the default cannot hold calls for more Fourier harmonics.
  double gyro_mix;

  // Gauss-Legendre nodes per transverse direction of a configuration face (0 => default 5, at
  // least 3).
  int num_face_nodes;

  // Boundary conditions of the distributions (see struct gkyl_pkpm_species).
  struct gkyl_vlasov_bcs bcx, bcy, bcz;
};

// Parameters for the fluid block of a PKPM species: the initial condition of the fluid state,
// the positivity limiter, momentum diffusion, the applied acceleration and the fluid BCs.
// Declared as the 'fluid' block of struct gkyl_pkpm_species. The members it shares with
// struct gkyl_vlasov_fluid_species have the same names and types.
struct gkyl_pkpm_fluid_species {
  void *ctx; // Context for initial condition init function.
  // Pointer to the initialization function of the fluid state
  // [rho, rho ux, rho uy, rho uz, E], E = rho |u|^2/2 + 3 p/2 = rho |u|^2/2 + (p_par + 2 p_perp)/2.
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // Vertex positivity limiter of the fluid state [rho, M, E] (Zhang-Shu scaling of the
  // slopes toward the cell averages, always applied, after every forward-Euler stage, the
  // implicit coupling and the initial condition): the vertex density and the vertex internal
  // energy are kept above positivity_eps times the cell values, so the thermal speed is real
  // at the Gauss-Legendre nodes. The cells it acts on are recorded in
  // <name>-<species>_pkpm_positivity_<frame>.gkyl (cumulative counts per cell: [density,
  // internal energy, inadmissible average]) and the counts per integrated-moment interval
  // in the dynvec <name>-<species>-positivity-stat.gkyl ([cells scaled for the density, for
  // the internal energy, cells with an inadmissible average, limiter applications], summed
  // over ranks).
  double positivity_eps; // Relative floor (0 => default 1e-2).

  // Diffusion of the momentum to include (also on a position map): order 2 (grad^2, the
  // default) or 4 (grad^4).
  struct gkyl_vlasov_fluid_diffusion diffusion;

  void *app_accel_ctx; // Context for applied acceleration function.
  // Pointer to applied acceleration function.
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // Set to true if applied acceleration function is time dependent.

  // Boundary conditions of the fluid state (see struct gkyl_pkpm_species).
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameters for one PKPM species: the identity (name, charge, mass) and its two blocks, as in
// struct gkyl_vlasov_species. For example
//   struct gkyl_pkpm_kinetic_species elc_kinetic = { .lower = ..., .cells = ..., ... };
//   struct gkyl_pkpm_fluid_species elc_fluid = { .init = ..., ... };
//   .num_species = 1,
//   .species = {
//     { .name = "elc", .charge = -1.0, .mass = 1.0, .type = GKYL_SPECIES_PKPM,
//       .kinetic = elc_kinetic, .fluid = elc_fluid },
//   },
// Boundary conditions: each block carries its own and each is applied to its own state.
// The kinetic BC fills the ghost cells of the distributions: GKYL_SPECIES_COPY, GKYL_SPECIES_ABSORB
// (empty ghost distributions; the wall losses are recorded and T_perp/T of the ghost cell is the
// mirror image of the skin cell), GKYL_SPECIES_REFLECT (the field-aligned mirror image) or
// GKYL_SPECIES_FIXED_FUNC (the initial condition; x direction only). The fluid BC fills the
// ghost cells of [rho, M, E], which set the frame (rho, u, xi) of the ghost cell: GKYL_SPECIES_COPY,
// GKYL_SPECIES_ABSORB (the mirror image of the skin cell, even in the wall direction, so the
// primitive variables and their gradients are one-sided at the wall), GKYL_SPECIES_REFLECT (the
// mirror image with the momentum reflected) or GKYL_SPECIES_FIXED_FUNC. The same type in both
// blocks is the natural pair: an absorbing wall records its losses, and a reflecting wall
// passes no mass and no energy.
struct gkyl_pkpm_species {
  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  enum gkyl_species_type type; // GKYL_SPECIES_PKPM: the species owns both blocks.

  struct gkyl_pkpm_kinetic_species kinetic; // Kinetic block (velocity grid, projection, ...).
  struct gkyl_pkpm_fluid_species fluid; // Fluid block (init, diffusion, ...).
};

// Parameter for EM field
struct gkyl_pkpm_field {
  enum gkyl_field_id field_id; // type of field
    // (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
  // set to true if field does not change in time (applied and ghost currents are then ignored)
  bool is_static;

  double epsilon0, mu0;
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // set to true if external electromagnetic field function is time dependent

  void *app_current_ctx; // context for applied current function
  // pointer to applied current function
  void (*app_current)(double t, const double *xn, double *app_current_out, void *ctx);
  bool app_current_evolve; // set to true if applied current function is time dependent

  // 1D periodic electrostatics: subtract the domain-mean current from Ampere's law, so that
  // the Maxwell field equals the Poisson field (cancels the current of a drifting
  // initial state, e.g. the Buneman instability).
  bool use_ghost_current;

  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter
  bool limit_em; // Optional input parameter for applying limiters to EM fields

  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];
};

// Per-direction configuration-space position map (non-uniform grid). The map in each
// direction is a function of the (scalar) computational coordinate of that direction; its
// C^0 piecewise linear representation gives a per-cell constant Jacobian J, and the
// evolved quantities carry J (J [F, G], J [rho, M, E], J [E, B]) while all output is physical.
struct pkpm_mapc2p_pos {
  void *mapc2p_pos_ctx; // context for the mapping function
  // pointer to the per-direction mapping: xc is the computational coordinate and on output
  // xp the corresponding physical coordinate. NULL => identity map in that direction.
  void (*mapc2p_pos_func)(double t, const double *xc, double *xp, void *ctx);
};

// Largest number of moments a neural network takes as inputs or as outputs.
#define GKYL_PKPM_NN_MAX_MOMS 16

// Neural-network options of a run, read by the driver of the run (the Lua wrapper): training
// one network per species on the diagnostic moments of the simulation (gkyl_pkpm_app_train,
// gkyl_pkpm_app_write_nn) and testing a trained one against them (gkyl_pkpm_app_test). All
// zeros, the default, means no network is used.
struct gkyl_pkpm_nn {
  bool train_nn; // Train neural network on simulation data?
  bool train_ab_initio; // Train a new network (otherwise the one read from train_nn_file)?
  int nn_width; // Number of neurons to use per layer.
  int nn_depth; // Number of layers to use.
  char train_nn_file[256]; // File path of neural network to train.
  int num_trains; // Number of times to train neural network.
  int num_nn_writes; // Number of times to write out neural network.
  int num_input_moms; // Number of "input" moments to train on.
  int input_moms[GKYL_PKPM_NN_MAX_MOMS]; // "Input" moments to train on.
  int num_output_moms; // Number of "output" moments to train on.
  int output_moms[GKYL_PKPM_NN_MAX_MOMS]; // "Output" moments to train on.
  bool test_nn; // Test neural network on simulation data?
  char test_nn_file[256]; // File path of neural network to test.
  int num_tests; // Number of times to test neural network.
};

// Top-level app parameters
struct gkyl_pkpm {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  double lower[3], upper[3]; // lower, upper bounds of config-space
  int cells[3]; // config-space cells
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  // Per-direction non-uniform configuration-space mapping (diagonal coordinate
  // stretch, independent of the general curvilinear mapc2p above which feeds the
  // wave geometry of the EM limiter). Directions without a function use the identity.
  struct pkpm_mapc2p_pos mapc2p_pos[GKYL_MAX_CDIM];

  double cfl_frac; // CFL fraction to use (default 1.0)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_pkpm_species species[GKYL_MAX_SPECIES]; // species objects

  struct gkyl_pkpm_field field; // field object

  // The default is a staggered (Strang) operator split with implicit fluid-EM coupling: between
  // updates the fluid and E lag by half a source step, and output and diagnostics see them
  // synchronized; frames store the steps that make a restart continue the uninterrupted run.
  bool use_explicit_sources; // Use the fully explicit SSP RK3 scheme instead

  struct gkyl_pkpm_nn nn; // Neural-network options (all zeros: no network).

  struct gkyl_app_parallelism_inp parallelism; // Parallelism-related inputs.
};

// Simulation statistics
struct gkyl_pkpm_stat {
  bool use_gpu; // did this sim use GPU?

  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method

  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure

  double total_tm; // time for simulation (not including ICs)
  double rk3_tm; // time for SSP RK3 step
  double pkpm_em_tm; // time for implicit fluid-EM coupling step
  double init_species_tm; // time to initialize all species
  double init_fluid_species_tm; // time to initialize all fluid species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species collisionless RHS
  double fluid_species_rhs_tm; // time to compute fluid species RHS

  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_coll_tm; // total time for collision updater (excluded moments)

  // time to compute the configuration-space variables of a stage: the moments of the
  // normalized distributions, the flow and thermal speed, the kinetic coupling moments
  // [p_par, p_perp, q_par, q_perp, s, r], the constraints and the gradients for the sources
  double species_pkpm_vars_tm;

  double species_bc_tm; // time to compute species BCs
  double fluid_species_bc_tm; // time to compute fluid species BCs
  double field_bc_tm; // time to compute field

  double field_rhs_tm; // time to compute field RHS
  double field_em_vars_tm; // time to compute EM auxiliary variables (e.g., bvar and E x B)
  double current_tm; // time to compute currents and accumulation

  long n_species_omega_cfl; // number of times CFL-omega all-reduce is called
  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl

  long n_field_omega_cfl; // number of times CFL-omega for field all-reduce is called
  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for field

  long n_mom; // calls to moment calculation
  double mom_tm; // time to compute moments

  long n_diag; // calls to diagnostics
  double diag_tm; // time to compute diagnostics

  long n_io; // number of calls to IO
  double io_tm; // time to perform IO
};

// Object representing pkpm app
typedef struct gkyl_pkpm_app gkyl_pkpm_app;

/**
 * Construct a new PKPM app.
 *
 * @param pkpm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New PKPM app object.
 */
gkyl_pkpm_app *gkyl_pkpm_app_new(struct gkyl_pkpm *pkpm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_pkpm_app_apply_ic(gkyl_pkpm_app *app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_pkpm_app_apply_ic_field(gkyl_pkpm_app *app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_pkpm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_pkpm_app_apply_ic_species(gkyl_pkpm_app *app, int sidx, double t0);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_pkpm_app_calc_integrated_mom(gkyl_pkpm_app *app, double tm);

/**
 * Calculate integrated L2 norm of the distribution function, f^2.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_pkpm_app_calc_integrated_L2_f(gkyl_pkpm_app *app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_pkpm_app_calc_field_energy(gkyl_pkpm_app *app, double tm);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_pkpm_app_write(gkyl_pkpm_app *app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_pkpm_app_write_field(gkyl_pkpm_app *app, double tm, int frame);

/**
 * Write the state of a species to files, both read by a restart:
 * 1. <name>_<frame>: the distributions [F, G, M0_perp (x, y, z), N_perp (x, y, z)]
 * 2. <name>_fluid_<frame>: the fluid state [rho, rho ux, rho uy, rho uz, E]
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_pkpm_app_write_species(gkyl_pkpm_app *app, int sidx, double tm, int frame);

/**
 * Write the PKPM diagnostic data of a species to files:
 * 1. _pkpm_moms_: the moments [rho, M1, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
 * 2. _pkpm_prim_: the pressure tensor and the heat flux of the fluid equations,
 *    [P_xx, P_xy, P_xz, P_yy, P_yz, P_zz, q_x, q_y, q_z]
 * 3. _pkpm_vars_: [ux, uy, uz, T_perp/m, a_P (pressure force), div(b), bb : grad(u),
 *    T_perp/m div(b), xi (thermal speed), c (heating-rate closure)]
 * 4. _pkpm_constraints_: [M0(F), M1(F), (M2(F) + 2 M0(G))/3, int M0_perp dw (x, y, z)]
 * 5. _pkpm_positivity_: the cumulative per-cell record of the positivity limiter
 *
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_pkpm_app_write_mom(gkyl_pkpm_app *app, int sidx, double tm, int frame);

/**
 * Write integrated diagnostic moments for species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_pkpm_app_write_integrated_mom(gkyl_pkpm_app *app);

/**
 * Write integrated L2 norm of the species distribution function to file. Integrated
 * L2 norm is appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_pkpm_app_write_integrated_L2_f(gkyl_pkpm_app *app);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_pkpm_app_write_field_energy(gkyl_pkpm_app *app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_pkpm_app_stat_write(gkyl_pkpm_app *app);

/**
 * Train neural network on PKPM moments data for each fluid species.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 * @param num_input_moms Number of "input" moments to train on.
 * @param input_moms Array of "input" moments to train on.
 * @param num_output_moms Number of "output" moments to train on.
 * @param output_moms Array of "output" moments to train on.
 * @param input_data Array of input data to train on (across the computational domain).
 * @param output_data Array of output data to train on (across the computational domain).
 */
void gkyl_pkpm_app_train(
  gkyl_pkpm_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data,
  struct gkyl_kn_vec *output_data
);

/**
 * Train neural network on PKPM moments data for a particular fluid species.
 * 
 * @param app App object.
 * @param sidx Index of fluid species to train on.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 * @param num_input_moms Number of "input" moments to train on.
 * @param input_moms Array of "input" moments to train on.
 * @param num_output_moms Number of "output" moments to train on.
 * @param output_moms Array of "output" moments to train on.
 * @param input_data Array of input data to train on (across the computational domain).
 * @param output_data Array of output data to train on (across the computational domain).
 */
void gkyl_pkpm_app_train_mom(
  gkyl_pkpm_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data, struct gkyl_kn_vec *output_data
);

/**
 * Write out PKPM moments neural network for each fluid species.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 */
void gkyl_pkpm_app_write_nn(gkyl_pkpm_app *app, double tm, int frame, struct gkyl_kann_net **ann);

/**
 * Write out PKPM moments neural network for a particular fluid species.
 * 
 * @param app App object.
 * @param sidx Index of fluid species to write out.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 */
void gkyl_pkpm_app_write_nn_mom(
  gkyl_pkpm_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann
);

/**
 * Test neural network on PKPM moments data for each fluid species, and write moments to a file.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 * @param num_input_moms Number of "input" moments to test on.
 * @param input_moms Array of "input" moments to test on.
 * @param num_output_moms Number of "output" moments to test on.
 * @param output_moms Array of "output" moments to test on.
 * @param input_data_real Array of real input data to validate against (across the computational domain).
 * @param output_data_real Array of real output data to validate against (across the computational domain).
 * @param output_data_predicted Array of predicted output data to validate (across the computational domain).
 */
void gkyl_pkpm_app_test(
  gkyl_pkpm_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data_real,
  struct gkyl_kn_vec *output_data_real, struct gkyl_kn_vec *output_data_predicted
);

/**
 * Test neural network on PKPM moments data for a particular fluid species, and write moments to a file.
 * 
 * @param app App object.
 * @param sidx Index of fluid species to test on.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 * @param num_input_moms Number of "input" moments to test on.
 * @param input_moms Array of "input" moments to test on.
 * @param num_output_moms Number of "output" moments to test on.
 * @param output_moms Array of "output" moments to test on.
 * @param input_data_real Array of real input data to validate against (across the computational domain).
 * @param output_data_real Array of real output data to validate against (across the computational domain).
 * @param output_data_predicted Array of predicted output data to validate (across the computational domain).
 */
void gkyl_pkpm_app_test_mom(
  gkyl_pkpm_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data_real, struct gkyl_kn_vec *output_data_real,
  struct gkyl_kn_vec *output_data_predicted
);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_file_field(gkyl_pkpm_app *app, const char *fname);

/**
 * Initialize pkpm species from file
 *
 * @param app App object
 * @param sidx Species index
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_file_species(
  gkyl_pkpm_app *app, int sidx, const char *fname
);

/**
 * Initialize pkpm fluid species from file
 *
 * @param app App object
 * @param sidx Species index
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_file_fluid_species(
  gkyl_pkpm_app *app, int sidx, const char *fname
);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_frame_field(gkyl_pkpm_app *app, int frame);

/**
 * Initialize pkpm species and fluid species from frame
 *
 * @param app App object
 * @param sidx Species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_frame_species(
  gkyl_pkpm_app *app, int sidx, int frame
);

/**
 * Initialize the pkpm app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_read_from_frame(gkyl_pkpm_app *app, int frame);

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
void gkyl_pkpm_app_cout(const gkyl_pkpm_app *app, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_pkpm_update(gkyl_pkpm_app *app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_pkpm_stat gkyl_pkpm_app_stat(gkyl_pkpm_app *app);

/**
 * Free PKPM app.
 *
 * @param app App to release.
 */
void gkyl_pkpm_app_release(gkyl_pkpm_app *app);
