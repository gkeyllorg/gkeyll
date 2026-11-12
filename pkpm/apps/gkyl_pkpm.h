#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <stdbool.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

// Parameters for species collisions
struct gkyl_pkpm_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)

  void *ctx; // context for collision function
  // Function for the self-collision frequency. Without it the frequency is the Spitzer one
  // computed from the local moments (see cross_nu).
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);

  double nuFrac; // Rescales collision frequencies (default = 1)
  // Reference parameters of the Coulomb logarithm of the computed collision frequencies. Unset
  // entries default to 1, i.e. normalized units.
  double den_ref; // Reference density
  double temp_ref; // Reference temperature
  double hbar, eps0; // Planck's constant/2 pi, vacuum permittivity

  int num_cross_collisions; // number of species to cross-collide with
  char collide_with[GKYL_MAX_SPECIES][128]; // names of species to cross collide with
  // Functions for the cross-collision frequencies, one per species in collide_with. A species
  // gives all of them or none, and so must its partners. With none, the frequencies are the
  // Spitzer ones computed from the local moments of the pair,
  //   nu_sr = nuFrac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2),
  // with the constant part alpha_E of Morse's energy-exchange coefficient at the reference
  // parameters (nu_ss = nuFrac alpha_E n_s / (2 vtsq_s)^(3/2) for the self frequency).
  void (*cross_nu[GKYL_MAX_SPECIES])(double t, const double *xn, double *fout, void *ctx);
  void *cross_nu_ctx[GKYL_MAX_SPECIES]; // contexts for the cross-collision functions

  // BGK collisions: one backward-Euler step after each time step instead of the explicit update
  // (with its CFL limit). Not available on the GPU.
  bool is_implicit;
};

// Parameters for a BGK-form temperature source of a species,
//   S = nu_Q(x) (f_M(n, u, T_Q(x)) - f),
// the form of the gyrokinetic GKYL_SOURCE_BGK_HEATING with a prescribed temperature T_Q. It
// conserves the density and the flow and relaxes the temperature, d_t T = nu_Q (T_Q - T). In
// the normalized variables the fluid thermal energy relaxes to (3/2) n T_Q at the rate nu_Q
// and the distributions relax to the unit Maxwellian at the rate nu_Q T_Q/T, which relaxes
// every pressure component at the rate nu_Q of the source. Terms of order (T_Q/T - 1) times
// the departure of the distributions from a Maxwellian are dropped: moments of degree three
// decay at nu_Q T_Q/T instead of nu_Q (3 T_Q/T - 1)/2. It is independent of `collisions`.
struct gkyl_pkpm_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // GKYL_NO_SOURCE_BGK (default) or GKYL_SOURCE_BGK_HEATING

  void *rate_profile_ctx; // context for the rate function
  // pointer to the source rate nu_Q(x)
  void (*rate_profile)(double t, const double *xn, double *fout, void *ctx);
  void *temp_shape_ctx; // context for the temperature function
  // pointer to the source temperature T_Q(x)
  void (*temp_shape)(double t, const double *xn, double *fout, void *ctx);

  // One backward-Euler step after each time step instead of the explicit update (whose time
  // step is limited by nu_Q and by nu_Q T_Q/T).
  bool is_implicit;
};

// Parameters for a species source: an LTE profile [density, isotropic temperature, drift]
// injected into the normalized kinetic equations (through the configuration pass, so the
// closures keep the constraints exact) and into the fluid state [rho, M, E]. With
// GKYL_PROJ_SOURCE the profile is added as given; with GKYL_BFLUX_SOURCE it is scaled every
// stage by the particle loss rate of `source_species` through its absorbing walls divided
// by `source_length` (so a density profile integrating to source_length replaces the
// losses exactly). The drift enters the parallel Maxwellian through its component along b;
// a perpendicular drift relative to the local flow is not carried by G.
struct gkyl_pkpm_source {
  enum gkyl_source_id source_id; // GKYL_NO_SOURCE (default), GKYL_PROJ_SOURCE, GKYL_BFLUX_SOURCE
  double source_length; // bflux: the length the density profile integrates to
  char source_species[128]; // bflux: species whose wall losses set the rate
  bool write_source; // write the fluid image of the profile at frame 0

  void *ctx_density, *ctx_temp, *ctx_V_drift; // contexts
  void (*density)(double t, const double *xn, double *fout, void *ctx); // number density
  void (*temp)(double t, const double *xn, double *fout, void *ctx); // isotropic temperature
  void (*V_drift)(
    double t, const double *xn, double *fout, void *ctx
  ); // drift velocity (3 components); optional
};

// Parameters of the momentum diffusion of a species (also on a position map)
struct gkyl_pkpm_fluid_diffusion {
  double D; // constant diffusion coefficient
  int order; // order of the diffusion: 2 (grad^2, the default) or 4 (grad^4)
  void *Dij_ctx; // context for applied diffusion function if using general diffusion tensor
  // pointer to applied diffusion function is using general diffusion tensor
  void (*Dij)(double t, const double *xn, double *Dout, void *ctx);
};

// Parameters for PKPM species
struct gkyl_pkpm_species {
  char name[128]; // species name

  enum gkyl_model_id model_id; // type of model
    // (e.g., SR, general geometry, PKPM, see gkyl_eqn_type.h)

  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of the NORMALIZED velocity w = v_par/xi
  int cells[3]; // velocity-space cells

  void *ctx_dist; // context for initial condition init function
  void *ctx_fluid; // context for initial condition init function
  // pointer to the initialization function of the normalized distributions [F, G]:
  // F is a probability density in w (M0 = 1, M1 = 0) and G its perpendicular-energy
  // partner (M0(G) = T_perp/T) with M2(F) + 2 M0(G) = 3; e.g. a bi-Maxwellian is
  // F = exp(-w^2/(2 r))/sqrt(2 pi r), G = (T_perp/T) F with r = T_par/T, T = (T_par + 2 T_perp)/3.
  void (*init_dist)(double t, const double *xn, double *fout, void *ctx);
  void *ctx_harm; // context for the first-harmonic initial condition
  // pointer to the initialization function of the first-harmonic blocks
  // [M0_perp (x, y, z), N_perp (x, y, z)] (optional; zero if not given). They must be perpendicular
  // to b, and int M0_perp dw must vanish.
  void (*init_harm)(double t, const double *xn, double *fout, void *ctx);
  // pointer to the initialization function of the fluid state
  // [rho, rho ux, rho uy, rho uz, E], E = rho |u|^2/2 + 3 p/2 = rho |u|^2/2 + (p_par + 2 p_perp)/2
  void (*init_fluid)(double t, const double *xn, double *fout, void *ctx);

  // collisions to include
  struct gkyl_pkpm_collisions collisions;

  // source to include
  struct gkyl_pkpm_source source;

  // temperature source to include
  struct gkyl_pkpm_source_bgk source_bgk;

  // diffusion coupling to include for momentum
  struct gkyl_pkpm_fluid_diffusion diffusion;

  // Gyrophase-mixing closure of the first harmonic: rate, in units of the species' local |Omega|, at
  // which agyrotropy beyond what one harmonic can represent relaxes back to that bound (0: the
  // default 10; a negative value switches the closure off). The overshoot above the bound is about
  // (drive rate/|Omega|)/gyro_mix; a drive the default cannot hold calls for more Fourier harmonics.
  double gyro_mix;

  void *app_accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*app_accel)(double t, const double *xn, double *app_accel_out, void *ctx);
  bool app_accel_evolve; // set to true if applied acceleration function is time dependent

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameter for EM field
struct gkyl_pkpm_field {
  enum gkyl_field_id field_id; // type of field
    // (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
  bool is_static; // set to true if field does not change in time

  double epsilon0, mu0;
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // set to true if external electromagnetic field function is time dependent

  void *app_current_ctx; // context for external electromagnetic fields function
  // pointer to external electromagnetic fields function
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

  bool use_explicit_source; // Use fully explicit SSP RK3 scheme
    // Default is a staggered (Strang) operator split with implicit fluid-EM coupling: between
    // updates the fluid and E lag by half a source step, and output and diagnostics see them
    // synchronized; frames store the steps that make a restart continue the uninterrupted run.

  // Vertex positivity limiter of the fluid state [rho, M, E] (Zhang-Shu scaling of the
  // slopes toward the cell averages, always applied, after every forward-Euler stage, the
  // implicit coupling and the initial condition): the vertex density, the vertex internal
  // energy and the vertex thermal speed xi (the p=1 interpolant of its Gauss-Legendre node
  // values) are kept above fluid_positivity_eps times the cell values (default 1e-2). The
  // configuration-space flux reads rho and xi on the faces. The cells it acts on are
  // recorded in <name>-<species>_pkpm_positivity_<frame>.gkyl (cumulative counts per cell:
  // [density, internal energy, thermal speed, inadmissible average]) and the counts per
  // integrated-moment interval in the dynvec <name>-<species>-positivity-stat.gkyl ([cells
  // scaled for the density, for the internal energy, for the thermal speed, cells with an
  // inadmissible average, limiter applications], summed over ranks).
  double fluid_positivity_eps; // Relative floor (0 => default 1e-2).
  int num_face_nodes; // Gauss-Legendre nodes per transverse direction of a configuration face (0 => default 5, at least 3).

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
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_pkpm_vars_tm; // time to compute pkpm vars
    // These are the coupling moments [rho, p_par, p_perp], the self-consistent
    // pressure force (div(p_par b_hat)), and the primitive variables
    // along with the acceleration variables in the kinetic equation
    // and the source distribution functions for Laguerre couplings.

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
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_pkpm_app_write_species(gkyl_pkpm_app *app, int sidx, double tm, int frame);

/**
 * Write pkpm auxiliar data to files. Includes:
 * 1. PKPM moments [rho, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp]
 * 2. PKPM fluid variables [rho, rho ux, rho uy, rho uz, 
 * P_xx + rho ux^2, P_xy + rho ux uy, P_xz + rho ux uz,
 * P_yy + rho uy^2, P_yz + rho uy uz, P_zz + rho uz^2]
 * 3. PKPM variables, including primitive variables (ux, uy, uz, T_perp/m, m/T_perp) and 
 * acceleration variables (div(b), 1/rho div(p_par b), T_perp/m div(b), bb : grad(u), and 
 * vperp configuration space characteristics = bb : grad(u) - div(u) - 2 nu.
 * 
 * @param app App object.
 * @param sidx Index of fluid species to initialize.
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
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_pkpm_app_from_file_species(
  gkyl_pkpm_app *app, int sidx, const char *fname
);

/**
 * Initialize pkpm fluid species from file
 *
 * @param app App object
 * @param sidx gk species index
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
 * @param sidx gk species index
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
