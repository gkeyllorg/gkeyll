#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_velocity_map.h>
#include <gkyl_position_map.h>
#include <gkyl_gyrokinetic_comms.h>
#include <gkyl_mom_type.h>
#include <gkyl_gk_bc_type.h>

#include <stdbool.h>


// Parameters for projection
struct gkyl_gyrokinetic_ic_import {
  // Inputs to initialize the species with the distribution from a file (f_in)
  // and to modify that distribution such that f = alpha(x)*f_in+beta(x,v).
  enum gkyl_ic_import_type type;
  char file_name[128]; // Name of file that contains IC, J*f_in.
  char jacobtot_inv_file_name[128]; // Name of file that contains 1/Jacobian. Used to get f from Jf.
  char jacobvel_file_name[128]; // Name of file that contains the velocity-space Jacobian.
  bool enforce_positivity; // =true sets f to 0 where it was negative.
  void *conf_scale_ctx;
  void (*conf_scale)(double t, const double *xn, double *fout, void *ctx); // alpha(x).
};

struct gkyl_gyrokinetic_projection {
  enum gkyl_projection_id proj_id; // type of projection (see gkyl_eqn_type.h)
  enum gkyl_quad_type quad_type; // quadrature scheme to use: defaults to Gaussian
  double f_floor; // Floor value of the distribution.

  union {
    struct {
      // Distribution function to project and its context. For fluid neutrals
      // this function returns mass density, momentum density, and total
      // energy density.
      void (*func)(double t, const double *xn, double *fout, void *ctx); 
      void *ctx_func; 
    };
    struct {
      // For Maxwellians (or BiMaxwellians), specify density, parallel speed
      // and temperature (or parallel and perpendicular temperature) functions,
      // and their context.
      void (*density)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_density;
      void (*upar)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_upar;
      void (*temp)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_temp;
      void (*temppar)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_temppar;
      void (*tempperp)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_tempperp;

      // Optionally read primitive moments from files.
      struct gkyl_gyrokinetic_ic_import maxwellian_moms_import;
      struct gkyl_gyrokinetic_ic_import bimaxwellian_moms_import;
      struct gkyl_gyrokinetic_ic_import density_import;
      struct gkyl_gyrokinetic_ic_import upar_import;
      struct gkyl_gyrokinetic_ic_import temp_import;
      struct gkyl_gyrokinetic_ic_import temppar_import;
      struct gkyl_gyrokinetic_ic_import tempperp_import;

      // For kinetic neutrals, specify density, drift velocity (udrift) and
      // temperature, and their context, if projecting a Maxwellian.
      void *ctx_udrift;
      void (*udrift)(double t, const double *xn, double *fout, void *ctx);

      // Options for a Maxwellian in vel-space w/ Gaussian shape in conf-space.
      double gaussian_mean[GKYL_MAX_CDIM]; // Center in configuration space.
      double gaussian_std_dev[GKYL_MAX_CDIM]; // Sigma in configuration space, function is constant if sigma is 0.
      double total_num_particles; // Total number of particle (M0 moment).
      double total_kin_energy; // Total kinetic energy (0.5*mass*M2 moment).
      double temp_max; // Maximum temperature of the Gaussian Maxwellian distribution.
      double temp_min; // Minimum temperature of the Gaussian Maxwellian distribution.
      
      // Boolean if we are correcting all the moments or only density.
      bool correct_all_moms; 
      double iter_eps; // Error tolerance for moment fixes (density is always exact).
      int max_iter; // Maximum number of iteration.
      bool use_last_converged; // Use last iteration value regardless of convergence?
    };
  };
};

struct gkyl_phase_diagnostics_inp {
  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.
  int num_integrated_diag_moments; // Number of integrated diagnostic moments.
  enum gkyl_distribution_moments integrated_diag_moments[12]; // List of integrated diagnostic moments.
  bool time_integrated; // Whether to use time integrated diags.
};

// Parameters for collisionless terms.
struct gkyl_gyrokinetic_collisionless {
  enum gkyl_gk_collisionless_type type; // Type of collisionless terms.
  bool write_diagnostics; // Whether to output diagnostics.
  double scale_factor; // Factor multiplying collisionless terms (should be > 0).
};

// Parameters for species collisions
struct gkyl_gyrokinetic_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)
  bool write_diagnostics; // Whether to output diagnostics.
  bool not_in_dfdt; // If true, the collision operator will not be added to df/dt.
    // Used to ignore the collisional updates of this species, while updating cross-species collisions.

  double nu_frac; // Rescales collision frequencies (default = 1).

  // Function for computing self-collision frequency.
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);
  void *self_nu_ctx; // Context for self_nu.

  int num_cross_collisions; // Number of species to collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species to collide with.

  // Functions for computing cross-collision frequencies (one for each num_cross_collisions).
  evalf_t cross_nu[GKYL_MAX_SPECIES];
  void *cross_nu_ctx; // Context for cross_nu.

  // Parameters used to compute the Coulomb Logarithm.
  double den_ref; // Reference density.
  double temp_ref; // Regerence temperature.
  double bmag_ref; // Reference magnetic field magnitude.
  double hbar, eps0, eV; // Planck's constant/2 pi, vacuum permittivity, elementary charge.

  // Boolean for using implicit BGK collisions (replaces rk3).
  bool is_implicit; 
};

// Parameters for species diffusion.
struct gkyl_gyrokinetic_diffusion {
  int num_diff_dir; // Number of diffusion directions.
  int diff_dirs[GKYL_MAX_CDIM]; // List of diffusion directions.
  double D[GKYL_MAX_CDIM]; // Constant diffusion coefficient in each direction.
  int order; // Order of the diffusion (4 for grad^4, 6 for grad^6, default is 2).
};

// Structure to hold parameters for adaptive source.
struct gkyl_gyrokinetic_adapt_source {
  bool adapt_particle; // Whether to adapt the particle source.
  bool adapt_energy; // Whether to adapt the energy source.
  char adapt_to_species[16]; // Species to adapt the particle loss to ensure quasi neutrality.
  int num_boundaries; // Number of boundaries to adapt.
  int dir[GKYL_MAX_CDIM*2]; // Direction to adapt.
  enum gkyl_edge_loc edge[GKYL_MAX_CDIM*2]; // Edge to adapt.
};

// Parameters for species source
struct gkyl_gyrokinetic_source {
  enum gkyl_source_id source_id; // Type of source.
  int num_sources;
  bool evolve; // Whether the source is time dependent.
  int num_adapt_sources;
  struct gkyl_gyrokinetic_adapt_source adapt[GKYL_MAX_SOURCES]; // Adaptive source parameters
  // sources using projection routine
  struct gkyl_gyrokinetic_projection projection[GKYL_MAX_SOURCES];

  struct gkyl_phase_diagnostics_inp diagnostics;
};

// Parameters for the anomalous diffusion term, d/dx D(x) df/dx.
struct gkyl_gyrokinetic_anomalous_diffusion {
  enum gkyl_gk_anomalous_diff_id anomalous_diff_id; // Type of diffusion term.
  void (*D_profile)(double t, const double *xn, double *fout, void *ctx); // D(x).
  void *D_profile_ctx;
  bool write_diagnostics; // Whether to output diagnostics.
};

// Parameters for species bgk source term either:
// 1. nu_Q(x)*(f_M(n,upar,T_Q(t)*s_Q(x)/m) - f) 
//    for default model or
// 2. nu(x)*(f_M(n_S,upar,T) - f)
//    nu(x)*(f_M(n,upar_S,T) - f)
//    nu(x)*(f_M(n,upar,T_S) - f)
//    for external model
struct gkyl_gyrokinetic_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // Type of BGK source  term.
  void (*rate_profile)(double t, const double *xn, double *fout, void *ctx); // nu_Q(x).
  void *rate_profile_ctx;
  void (*temp_shape)(double t, const double *xn, double *fout, void *ctx); // s_Q(x).
  void *temp_shape_ctx;
  double power; // Desired heating power (sets T_Q(t)).
  double injection_time;  // Injection time for external source model
                          // nu(x) = 1/coupling time
  double damping_factor;  // For external source model
                          // n_s = max(n_s, damping factor*n)
                          // to prevent driving n negative
  bool write_diagnostics; // Whether to output diagnostics.
};

// Emitting BCs inputs.
struct gkyl_gyrokinetic_emission_inp {
  int num_species;
  char in_species[GKYL_MAX_SPECIES][128];
  double recycling_frac; // Recycling coefficient.
  double emission_temp; // Temperature of emitted species.
};

// Parameters for boundary conditions
struct gkyl_gyrokinetic_bc {
  int dir;  // Direction in which BC is specified.
  enum gkyl_edge_loc edge; // Which edge this BC is for.
  enum gkyl_gyrokinetic_bc_type type; // BC type flag.
  double value[3]; // Meaning depends on type.
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx); // Auxiliary function (e.g. wall potential).
  void *aux_ctx; // Context for aux_profile.
  struct gkyl_gyrokinetic_projection projection; // Projection object input (e.g. for FIXED_FUNC).
  struct gkyl_gyrokinetic_emission_inp emission; 
  bool write_diagnostics; // Whether to output diagnostics.
  int bidx; // Block index (for multiblock solver).
};

struct gkyl_gyrokinetic_geometry {
  enum gkyl_geometry_id geometry_id;
  char geometry_path[128]; // Path to geometry files

  // Pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);
  void *c2p_ctx; // Context for mapc2p function.

  // Pointer to bfield function.
  void (*bfield_func)(double t, const double *xc, double *xp, void *ctx);
  void *bfield_ctx; // Context for bfield function.

  double world[3]; // Extra computational coordinates for cases with reduced dimensionality.

  bool has_LCFS; // Whether the geometry has a last closed flux surface (LCFS).
  double x_LCFS; // x coordinate of the LCFS.

  struct gkyl_efit_inp efit_info; // Context with RZ data such as efit file for a tokamak or mirror.
  struct gkyl_tok_geo_grid_inp tok_grid_info; // Context for tokamak geometry with computational domain info.
  struct gkyl_mirror_geo_grid_inp mirror_grid_info; // Context for mirror geometry with computational domain info.
  struct gkyl_position_map_inp position_map_info; // Position map object.
};

// Parameters for species radiation.
struct gkyl_gyrokinetic_radiation {
  enum gkyl_radiation_id radiation_id; // Type of radiation.

  int num_cross_collisions; // Number of species to cross-collide with.
  char collide_with[2*GKYL_MAX_SPECIES][128]; // Names of species to cross collide with.

  int atomic_Z[2*GKYL_MAX_SPECIES]; // Z of element for desired fit information.
  int charge_state[2*GKYL_MAX_SPECIES]; // Charge state of element for desired fit information.
  int num_of_densities[GKYL_MAX_RAD_DENSITIES]; // Max number of densities to use per charge state

  // Reference, max and min electron densities to specify range of density fits
  // (when more than 1 density).
  double reference_ne, max_ne, min_ne;

  enum gkyl_te_min_model te_min_model; // How the radiation is turned off (constant or varying Te).
  double Te_min; // Minimum temperature (in J) at which to stop radiating.

  bool write_diagnostics; // Whether to write diagnostics.
};

struct gkyl_gyrokinetic_react_type {
  enum gkyl_react_id react_id; // What type of reaction (ionization, charge exchange, recombination).
  enum gkyl_react_self_type type_self; // What is the role of species in this reaction.
  enum gkyl_ion_type ion_id; // What type of ion is reacting.
  char elc_nm[128]; // Names of electron species in reaction.
  char ion_nm[128]; // Name of ion species in reaction.
  char donor_nm[128]; // Name of donor species in reaction.
  char recvr_nm[128]; // Name of receiver species in reaction.
  char partner_nm[128]; // Name of neut species in cx reaction.
  int charge_state; // Charge state of species in reaction.
  double ion_mass; // Mass of ion in reaction.
  double elc_mass; // Mass of electron in reaction.
  double partner_mass; // Mass of neutral in cx reaction.
  double vt_sq_ion_min;
  double vt_sq_partner_min;
};

// Parameters for species reaction
struct gkyl_gyrokinetic_react {
  int num_react; // number of reactions
  // 3 types of reactions supported currently
  // Ionization, Charge exchange, and Recombination
  // GKYL_MAX_SPECIES number of reactions supported per species (8 different reactions)
  struct gkyl_gyrokinetic_react_type react_type[GKYL_MAX_REACT];
  bool write_diagnostics; // Whether to write diagnostics.
};

// Types of scaling operations for a species.
enum gkyl_gk_species_scaling_type {
  GKYL_GK_SPECIES_SCALING_NONE = 0, // No scaling.
  GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE, // Balance between recycling and ionization.
  GKYL_GK_SPECIES_SCALING_FIXED_FRACTION, // Maintains fixed fraction relative to another species.
  GKYL_GK_SPECIES_SCALING_BOLTZMANN, // n_s = n_{s,sheath}*exp(-q_s*(phi-phi_sheath)/T_s).
};

// Input parameters for scaling a species every time step.
struct gkyl_gyrokinetic_scaling_inp {
  enum gkyl_gk_species_scaling_type type; // Type of scaling operation.

  // Info for GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE.
  int num_boundaries; // Number of boundaries.
  int boundaries_dir[GKYL_MAX_CDIM*2]; // Direction of boundaries.
  enum gkyl_edge_loc boundaries_edge[GKYL_MAX_CDIM*2]; // Edge of boundaries.
  char impacting_ion_name[128]; // Name of impacting species.
  enum gkyl_ion_type impacting_ion_id; // Type of impacting ion.
  char electron_name[128]; // Name of electron species.
  double recycling_coeff; // Recycling coefficient.
  
  // Info for GKYL_GK_SPECIES_SCALING_FIXED_FRACTION.
  char ref_species_name[128]; // Name of reference species.
  double fixed_fraction; // Fraction of reference species density.

  bool write_diagnostics; // Whether to write diagnostics.
};

// Parameters in FLR effects.
struct gkyl_gyrokinetic_flr {
  enum gkyl_gk_flr_type type; 
  double Tperp; // Perp temperature used to evaluate gyroradius. 
  double bmag; // Magnetic field used to evaluate gyroradius. If not provided
               // it'll use B in the center of the domain.
};

struct gkyl_gyrokinetic_correct_inp {
  bool correct_all_moms; // boolean if we are correcting all the moments or only density
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iteration
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge.   
};

enum gkyl_gyrokinetic_positivity_type {
  GKYL_GK_POSITIVITY_NONE = 0, // Do not enforce positivity (default).
  GKYL_GK_POSITIVITY_SHIFT, // Shift f to zero if <0 at Gauss-Legendre nodes.
  GKYL_GK_POSITIVITY_MRS_LIMITER, // Use the More-Rossmanith-Seal limiter, and shift when needed.
};

struct gkyl_gyrokinetic_positivity {
  enum gkyl_gyrokinetic_positivity_type type; // Type of positivity enforcement algorithm.
  bool quasineutrality_rescale; // Whether to rescale this species to enforce quasineutrality in the simulation.
  bool write_diagnostics; // Whether to output diagnostics.
};

enum gkyl_gyrokinetic_damping_type {
  GKYL_GK_DAMPING_NONE = 0,
  GKYL_GK_DAMPING_USER_INPUT,
  GKYL_GK_DAMPING_LOW_PASS_FILTER,
};

struct gkyl_gyrokinetic_damping {
  // Add a damping term to the RHS of the gyrokinetic equation
  //   df/dt = - rate(z) * f
  // with the function rate(z) being:
  //   - a function given by the user (type = GKYL_GK_DAMPING_USER_INPUT).
  //   Alternativly, a low-pass filter can be applied (type = GKYL_GK_DAMPING_LOW_PASS_FILTER), which projects f to a lower-dimensional function fbar(z) and applies the damping as df/dt = - rate(z) * (f - fbar).
  enum gkyl_gyrokinetic_damping_type type;

  double rate_const;
  void (*rate_profile)(double t, const double *xn, double *fout, void *ctx);
  void *rate_profile_ctx; // Context for rate_profile function.
  int num_quad; // Number of quadrature points in each direction to use in projecting the rate.
  bool write_rate; // Whether to write damping-rate diagnostics.
  bool write_fbar; // For low-pass filter damping, write fbar diagnostics each output frame.
  bool cellwise_const; // For low-pass filter damping, whether fbar uses a single value per cell.
  bool do_not_reset_fbar; // On runtime reset, preserve the existing low-pass-filter fbar instead of reinitializing from f.
};

enum gkyl_gyrokinetic_fdot_multiplier_type {
  GKYL_GK_FDOT_MULTIPLIER_NONE = 0,
  GKYL_GK_FDOT_MULTIPLIER_USER_INPUT,
  GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE,
};

struct gkyl_gyrokinetic_fdot_multiplier {
  enum gkyl_gyrokinetic_fdot_multiplier_type type;
  void (*profile)(double t, const double *xn, double *fout, void *ctx); // Profile to multiply df/dt by.
  void *profile_ctx; // Context for profile function.
  bool cellwise_const; // Whether the multiplier has a single value per cell.
  bool write_diagnostics; // Whether to output diagnostics.
};

// Parameters for gk species.
struct gkyl_gyrokinetic_species {
  char name[128]; // Species name.

  double charge, mass; // Charge and mass.
  int vdim; // Velocity-space dimensions.
  double lower[GKYL_MAX_VDIM], upper[GKYL_MAX_VDIM]; // Lower, upper bounds of velocity-space.
  int cells[GKYL_MAX_VDIM]; // Velocity-space cells.

  struct gkyl_mapc2p_inp mapc2p;

  bool is_static; // Set to true if species does not change in time.

  struct gkyl_gyrokinetic_positivity positivity; // Positivity enforcement options.

  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;
  // Initial conditions from a file.
  struct gkyl_gyrokinetic_ic_import init_from_file;

  // Phase-space field multiplying df/dt.
  struct gkyl_gyrokinetic_fdot_multiplier time_rate_multiplier;

  double polarization_density; // Density factor in LHS of quasineutrality eqn.

  // Whether to scale the density using a polarization solve
  bool scale_with_polarization;

  int num_diag_moments; // number of diagnostic moments
  enum gkyl_distribution_moments diag_moments[12]; // list of diagnostic moments
  int num_integrated_diag_moments; // Number of integrated diagnostic moments.
  enum gkyl_distribution_moments integrated_diag_moments[12]; // List of integrated diagnostic moments.
  bool time_rate_diagnostics; // Whether to ouput df/dt diagnostics.
  bool write_omega_cfl; // Whether to ouput dt diagnostic for the CFL constraint.

  struct gkyl_gyrokinetic_collisionless collisionless; // Collisionless terms.

  struct gkyl_gyrokinetic_flr flr; // Options for FLR effects.

  // Diagnostics of the fluxes of f at position-space boundaries.
  struct gkyl_phase_diagnostics_inp boundary_flux_diagnostics;

  // Input quantities used by LTE (local thermodynamic equilibrium, or Maxwellian) projection
  // This projection operator is used by BGK collisions and all reactions.
  struct gkyl_gyrokinetic_correct_inp correct; 

  // Elastic collisions.
  struct gkyl_gyrokinetic_collisions collisions;

  // Source of particles/momentum/energy.
  struct gkyl_gyrokinetic_source source;

  // A damping term -rate*f on RHS.
  struct gkyl_gyrokinetic_damping damping; 

  // Anomalous diffusion.
  struct gkyl_gyrokinetic_anomalous_diffusion anomalous_diffusion;

  // BGK source.
  struct gkyl_gyrokinetic_source_bgk source_bgk;

  // Line radiation.
  struct gkyl_gyrokinetic_radiation radiation;

  // Reactions between plasma species.
  struct gkyl_gyrokinetic_react react;
  // Reactions with neutral species.
  struct gkyl_gyrokinetic_react react_neut;

  // Inputs to operation that scales the species every time step.
  struct gkyl_gyrokinetic_scaling_inp scaling;

  // Boundary conditions.
  struct gkyl_gyrokinetic_bc bcs[2*GKYL_MAX_CDIM];
};

// Parameters for neutral species.
struct gkyl_gyrokinetic_neut_species {
  char name[128]; // Species name.

  double mass; // Mass.
  int vdim; // Velocity-space dimensions.
  double lower[GKYL_MAX_VDIM], upper[GKYL_MAX_VDIM]; // Lower, upper bounds of velocity-space.
  int cells[GKYL_MAX_VDIM]; // Velocity-space cells.

  struct gkyl_mapc2p_inp mapc2p;

  bool is_static; // Set to true if neutral species does not change in time.

  struct gkyl_gyrokinetic_positivity positivity; // Positivity enforcement options.
  
  struct gkyl_gyrokinetic_ic_import init_from_file;
  
  // Initial conditions using projection routine.
  struct gkyl_gyrokinetic_projection projection;

  int num_diag_moments; // Number of diagnostic moments.
  enum gkyl_distribution_moments diag_moments[12]; // List of diagnostic moments.

  struct gkyl_gyrokinetic_collisionless collisionless; // Collisionless terms.

  // Diagnostics of the fluxes of f at position-space boundaries.
  struct gkyl_phase_diagnostics_inp boundary_flux_diagnostics;

  // Input quantities used by LTE (local thermodynamic equilibrium, or Maxwellian) projection
  // This projection operator is used by BGK collisions and all reactions.
  struct gkyl_gyrokinetic_correct_inp correct; 

  // Elastic collisions.
  struct gkyl_gyrokinetic_collisions collisions;

  // Source of particles/momentum/energy.
  struct gkyl_gyrokinetic_source source;

  // Reactions with plasma species.
  struct gkyl_gyrokinetic_react react_neut;

  // Boundary conditions.
  struct gkyl_gyrokinetic_bc bcs[2*GKYL_MAX_CDIM];

  double gas_gamma; // Adiabatic index (fluid neutrals).

  // Inputs to operation that scales the species every time step.
  struct gkyl_gyrokinetic_scaling_inp scaling;
};

// Parameter for gk field.
struct gkyl_gyrokinetic_field {
  enum gkyl_gkfield_id gkfield_id;
  bool is_static; // =true field does not change in time.
  bool zero_init_field; // =true doesn't compute the initial field.

  double polarization_bmag; // B factor in the polarization density.
  double kperpSq; // kperp^2 parameter for 1D field equations

  // parameters for adiabatic electrons simulations
  double electron_mass, electron_charge, electron_density, electron_temp;

  struct gkyl_gyrokinetic_bc poisson_bcs[2*GKYL_MAX_CDIM];

  bool time_rate_diagnostics; // Writes the time rate of change of field energy.

  // Initial potential used to compute the total polarization density.
  void (*polarization_potential)(double t, const double *xn, double *out, void *ctx);
  void *polarization_potential_ctx;
  struct gkyl_gyrokinetic_ic_import polarization_potential_import;

  // Interface to read a potential from file.
  struct gkyl_gyrokinetic_ic_import init_from_file;

  // Profile for the field at t=0.
  void (*init_field_profile)(double t, const double *xn, double *out, void *ctx);
  void *init_field_profile_ctx;

  void *phi_wall_lo_ctx; // context for biased wall potential on lower wall
  // pointer to biased wall potential on lower wall function
  void (*phi_wall_lo)(double t, const double *xn, double *phi_wall_lo_out, void *ctx);
  bool phi_wall_lo_evolve; // set to true if biased wall potential on lower wall function is time dependent  

  void *phi_wall_up_ctx; // context for biased wall potential on upper wall
  // pointer to biased wall potential on upper wall function
  void (*phi_wall_up)(double t, const double *xn, double *phi_wall_up_out, void *ctx);
  bool phi_wall_up_evolve; // set to true if biased wall potential on upper wall function is time dependent  

  struct gkyl_poisson_bias_line_list *bias_line_list; // Biased lines constraining the solution.
};

struct gkyl_gyrokinetic_eirene {
  char input_data_path[128]; // Path to EIRENE data
  char output_data_path[128]; // Path to EIRENE data
  double injection_time[GKYL_MAX_SPECIES]; // Injection time for each species
  double damping_factor[GKYL_MAX_SPECIES]; // Damping Factor
  double core_coll_factor[GKYL_MAX_SPECIES]; // core rate is increased by this factor
  int num_coupling_species; // number of species to couple
  char coupling_species[GKYL_MAX_SPECIES][128]; // Names of species to couple
};

// Top-level app parameters
struct gkyl_gk {
  char name[128]; // Name of app: used as output prefix.

  int cdim; // Configuration-space dimensions.
  double lower[3], upper[3]; // Lower, upper bounds of config-space.
  int cells[3]; // Config-space cells.
  int poly_order; // Polynomial order.
  enum gkyl_basis_type basis_type; // Type of basis functions to use.

  struct gkyl_gyrokinetic_geometry geometry; // Geometry input struct.

  double cfl_frac; // CFL fraction to use (default 1.0).
  double cfl_frac_omegaH; // CFL fraction used for the omega_H dt (default 1.0).

  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[3]; // List of periodic directions.

  int num_species; // Number of species.
  struct gkyl_gyrokinetic_species species[GKYL_MAX_SPECIES]; // Species objects.

  int num_neut_species; // Number of species.
  struct gkyl_gyrokinetic_neut_species neut_species[GKYL_MAX_SPECIES]; // Species objects.
  
  struct gkyl_gyrokinetic_field field; // Field object.

  struct gkyl_gyrokinetic_eirene eirene; // EIRENE input

  struct gkyl_app_parallelism_inp parallelism; // Parallelism-related inputs.
};

// Simulation statistics
struct gkyl_gyrokinetic_stat {
  bool use_gpu; // did this sim use GPU?
  
  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
    
  double init_species_tm; // time to initialize all species
  double init_neut_species_tm; // time to initialize all neutral species

  double time_loop_tm; // time for simulation (not including ICs)
  double fwd_euler_tm; // Time for forward euler
  double fwd_euler_step_f_tm; // Time spent on fwd euler step_f.
  double dfdt_dt_reduce_tm; // Time spent on fwd euler dt reduction.

  double species_collisionless_tm; // Time to compute species collisionless RHS (alpha already computed).
  double species_gyroavg_tm; // Time to compute species collisionless RHS.
  double species_lte_tm; // total time for species LTE (local thermodynamic equilibrium) projection updater
  double species_bflux_calc_tm; // Time for species boundary flux calculation.
  double species_bflux_moms_tm; // Time for species boundary flux moments.
  double species_coll_mom_tm; // time needed to compute various moments needed in collisions
  double species_coll_tm; // total time for collision updater (excluded moments)
  double species_damp_tm; // Time to accumulate species damping onto RHS.
  double species_fdot_mult_tm; // Time to spent on the df/dt multiplier.
  double species_diffusion_tm; // Time to compute species diffusion term.
  double species_rad_mom_tm; // total time to compute various moments needed in radiation operator
  double species_rad_tm; // total time for radiation operator
  double species_react_mom_tm; // total time to compute various moments needed in reactions 
  double species_react_tm; // total time for reactions updaters
  double species_src_tm; // Time to accumulate species source onto RHS.
  double species_source_bgk_tm; // Time to compute bgk source term RHS.
  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl

  double neut_species_collisionless_tm; // Time to compute neutral species collisionless RHS.
  double neut_species_lte_tm; // total time for neutral species LTE (local thermodynamic equilibrium) projection updater
  double neut_species_bflux_calc_tm; // Time for neutral species boundary flux calculation.
  double neut_species_bflux_moms_tm; // Time for neutral species boundary flux moments.
  double neut_species_coll_tm; // total time for neutral self-collisions updater (excluded moments)
  double neut_species_coll_mom_tm; // time needed to compute various moments needed in neutral self-collisions
  double neut_species_react_mom_tm; // total time to compute various moments needed in neutral reactions
  double neut_species_react_tm; // total time for neutral reactions updaters
  double neut_species_src_tm; // Time to accumulate neutral species source onto RHS.
  double neut_species_omega_cfl_tm; // time spent in all-reduce for omega-cfl for neutrals

  double time_rate_diags_tm; // Time spent on time rate of change diagnostics.
  double fdot_tm; // Time spent on computing \dot{f} diagnostics.
  double phidot_tm; // Time spent on computing \dot{phi} diagnostics.

  double field_tm; // Time to compute fields.
  double field_phi_rhs_tm; // Time spent on poisson eqn RHS.
  double field_phi_solve_tm;   // Time spent to solve poisson eqn.

  double bc_tm; // Time to compute BCs.
  double species_bc_tm; // Time to compute species BCs.
  double neut_species_bc_tm; // Time to compute neutral species BCs.

  double time_stepper_arithmetic_tm; // Time spent on arithmetic ops in time stepper.

  double pos_shift_tm; // Time spent on positivity shift.
  double species_pos_shift_tm; // Time spent on species positivity shift.
  double neut_species_pos_shift_tm; // Time spent on neutral species positivity shift.
  double pos_shift_quasineut_tm; // Time spent on positivity shift quasineutrality.

  // Group timers: additions of several of the above timers.
  double fwd_euler_sum_tm;
  double field_sum_tm;
  double bc_sum_tm;
  double time_rate_diags_sum_tm;
  double pos_shift_sum_tm;
  double time_stepper_sum_tm;
  double io_sum_tm;

  double species_io_tm; // Time to write the species distribution.
  double species_diag_calc_tm; // Time to compute species diagnostics.
  double species_diag_io_tm; // Time to write species diagnostics.

  double neut_species_io_tm; // Time to write the neutral species distribution.
  double neut_species_diag_calc_tm; // Time to compute neutral species diagnostics.
  double neut_species_diag_io_tm; // Time to write neutral species diagnostics.

  double field_io_tm; // Time to write the fields.
  double field_diag_calc_tm; // Time to compute field diagnostics.
  double field_diag_io_tm; // Time to write field diagnostics.

  double app_io_tm; // Time to write common diagnostics.
  double io_tm; // Time to write common diagnostics.

  long n_iter_corr[GKYL_MAX_SPECIES]; // total number of iterations used to correct species LTE projection
  long num_corr[GKYL_MAX_SPECIES]; // total number of times correction updater for species LTE projection is called
  long neut_n_iter_corr[GKYL_MAX_SPECIES]; // total number of iterations used to correct neutral species LTE projection
  long neut_num_corr[GKYL_MAX_SPECIES]; // total number of times correction updater for neutral species LTE projection is called

  long n_species_omega_cfl; // number of times CFL-omega all-reduce is called
  long n_mom; // total number of calls to moment updater routines
  long n_diag; // total number of calls to diagnostics
  long n_io; // number of calls to IO
  long n_diag_io; // number of calls to IO for diagnostics

  long n_neut_species_omega_cfl; // number of times CFL-omega all-reduce is called for neutrals
  long n_neut_mom; // total number of calls to neutrals moment updater routines
  long n_neut_diag; // total number of calls to diagnostics for neutral species
  long n_neut_io; // number of calls to IO for neutral species
  long n_neut_diag_io; // number of calls to IO for neutral species diagnostics

  long n_field_diag; // total number of calls to diagnostics for field
  long n_field_io; // number of calls to IO for field
  long n_field_diag_io; // number of calls to IO for field diagnostics
};

// Object representing gk app
typedef struct gkyl_gyrokinetic_app gkyl_gyrokinetic_app;

/**
 * Construct a new gk app.
 *
 * @param gk App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New gk app object.
 */
gkyl_gyrokinetic_app* gkyl_gyrokinetic_app_new(struct gkyl_gk *gk);

/**
 * Construct a new gk app (geometry only).
 *
 * @param gk App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New gk app object.
 */
gkyl_gyrokinetic_app* gkyl_gyrokinetic_app_new_geom(struct gkyl_gk *gk);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_gyrokinetic_app_apply_ic(gkyl_gyrokinetic_app* app, double t0);

/**
 * Write geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_geometry(gkyl_gyrokinetic_app *app, struct gkyl_gk_geometry_inp *geometry_inp);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_app_write_field(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_gyrokinetic_app_calc_field_energy(gkyl_gyrokinetic_app* app, double tm);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_field_energy(gkyl_gyrokinetic_app* app);

/**
 * Write eirene data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_app_write_eirene_diagnostics(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Calculate eirene integrated diagnostics.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 */
void gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(gkyl_gyrokinetic_app* app, double tm);

/**
 * Write eirene integrated diagnostics.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(gkyl_gyrokinetic_app* app);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species data to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose moments to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose moments to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species whose integrated moments to compute.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species.
 *
 * @param app App object.
 * @param sidx Index of neutral species whose integrated moments to compute.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate the L2 norm of a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species whose L2 norm to compute.
 * @param tm Time at which L2 norm is to be computed.
 */
void gkyl_gyrokinetic_app_calc_species_L2norm(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments of the boundary fluxes for a plasma species.
 *
 * @param app App object.
 * @param sidx Index of species whose integrated moments to compute.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments of the boundary fluxes for a neutral species.
 *
 * @param app App object.
 * @param sidx Index of species whose integrated moments to compute.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_neut_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write the L2 norm of a charged species to file. Subsequent calculations of
 * L2 norm are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose L2 norm to write out.
 */
void gkyl_gyrokinetic_app_write_species_L2norm(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write integrated diagnostic moments of the boundary fluxes for charged
 * species to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write integrated diagnostic moments of the boundary fluxes for neutral
 * species to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_neutral_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose source to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose source to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species source to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose source moments to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_source_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for neutral species source to file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose source moments to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species_source_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for a plasma species source.
 *
 * @param app App object.
 * @param sidx Index of species whose integrated moments to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a neutral species source.
 *
 * @param app App object.
 * @param sidx Index of neutral species whose integrated moments to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Calculate integrated diagnostic moments for a plasma species BGK source.
 *
 * @param app App object.
 * @param sidx Index of species whose integrated moments to write.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm);

/**
 * Write integrated diagnostic moments for charged species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose source integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write integrated diagnostic moments for charged species BGK source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose source integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write integrated diagnostic moments for neutral species source to file. Integrated
 * moments are appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose source integrated moments to write.
 */
void gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write LBO collisional moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose LBO moments to write.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_app_write_species_lbo_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write BGK collisional moments for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp.
 * @param frame Frame number.
 */
void gkyl_gyrokinetic_app_write_species_bgk_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of species whose Maxwellian correction status to write.
 */
void gkyl_gyrokinetic_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write neutral species integrated correct Maxwellian status of the to file. 
 * Correct Maxwellian status is appended to the same file.
 * 
 * @param app App object.
 * @param sidx Index of neutral species whose Maxwellian correction status to write.
 */
void gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write radiation drag coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose radiation drag coefficients to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_rad_drag(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Write radiation emissivity of each species that species sidx collides with
 * 
 * @param app App object.
 * @param sidx Index of species whose radiation emissivity to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_rad_emissivity(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame);

/**
 * Calculate integrated diagnostic moments of the radiation model.
 *
 * @param app App object.
 * @param sidx Index of species to whose radiation int moments to compute.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx, double tm);

/**
 * Write integrated moments of radiation rhs for radiating species 
 * 
 * @param app App object.
 * @param sidx Index of species to whose radiation int moments to write.
 */
void gkyl_gyrokinetic_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_iz_react(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write iz react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species to whose ionization diagnostics to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_iz_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose recombination diagnostics to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_recomb_react(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write recomb react rate coefficients for species to file.
 * 
 * @param app App object.
 * @param sidx Index of species whose recombination diagnostics to write.
 * @param ridx Index of reaction to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_recomb_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write the phase-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a charged species.
 * 
 * @param app App object.
 * @param sidx Index of species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write the conf-space diagnostics for a neutral species.
 * 
 * @param app App object.
 * @param sidx Index of neutral species to write.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_neut_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for all species (including sources) to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_mom(gkyl_gyrokinetic_app *app, double tm, int frame);

/**
 * Calculate integrated diagnostic moments for all species (including sources).
 *
 * @param app App object.
 * @param tm Time at which integrated diagnostics are to be computed
 */
void gkyl_gyrokinetic_app_calc_integrated_mom(gkyl_gyrokinetic_app* app, double tm);

/**
 * Calculate the L2 norm for all species.
 *
 * @param app App object.
 * @param tm Time at which L2 norms are to be computed.
 */
void gkyl_gyrokinetic_app_calc_L2norm(gkyl_gyrokinetic_app* app, double tm);

/**
 * Write integrated diagnostic moments for all species (including sources)
 * to file. Integrated moments are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_integrated_mom(gkyl_gyrokinetic_app *app);

/**
 * Write the L2 norm for all species (including sources)
 * to file. L2 norms are appended to the same file.
 * 
 * @param app App object.
 */
void gkyl_gyrokinetic_app_write_L2norm(gkyl_gyrokinetic_app *app);

/**
 * Write phase space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_phase(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Write configuration space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write_conf(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Write both conf and phase-space diagnostics to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_gyrokinetic_app_write(gkyl_gyrokinetic_app* app, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_app_stat_write(gkyl_gyrokinetic_app* app);

/**
 * Print timing of solver components to iostream.
 *
 * @param app App object.
 * @param iostream Where to write timers to (e.g. stdout, stderr);
 */
void
gkyl_gyrokinetic_app_print_timings(gkyl_gyrokinetic_app* app, FILE *iostream);

/**
 * Record the time step (in private dynvector).
 *
 * @param app App object.
 * @param tm Time stamp.
 * @param dt Time step to record (e.g. provided by app's status object).
 */
void
gkyl_gyrokinetic_app_save_dt(gkyl_gyrokinetic_app* app, double tm, double dt);

/**
 * Write the time step over time.
 *
 * @param app App object.
 */
void
gkyl_gyrokinetic_app_write_dt(gkyl_gyrokinetic_app* app);

/**
 * Read geometry file.
 *
 * @param app App object.
 */
void gkyl_gyrokinetic_app_read_geometry(gkyl_gyrokinetic_app *app, struct gkyl_gk_geometry_inp *geometry_inp);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_file_field(gkyl_gyrokinetic_app *app, const char *fname);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param fname file to read
 */
struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_neut_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname);

/**
 * Initialize the gyrokinetic app from a specific frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_app_read_from_frame(gkyl_gyrokinetic_app *app, int frame);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_field(gkyl_gyrokinetic_app *app, int frame);

/**
 * Initialize gyrokinetic species from file
 *
 * @param app App object
 * @param sidx gk species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_species(gkyl_gyrokinetic_app *app, int sidx, int frame);

/**
 * Initialize neutral species from file
 *
 * @param app App object
 * @param sidx neut species index
 * @param frame frame to read
 */
struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_neut_species(gkyl_gyrokinetic_app *app, int sidx, int frame);

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
void gkyl_gyrokinetic_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_gyrokinetic_update(gkyl_gyrokinetic_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics object.
 */
struct gkyl_gyrokinetic_stat gkyl_gyrokinetic_app_stat(gkyl_gyrokinetic_app* app);

/**
 * Free gk app.
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_app_release(gkyl_gyrokinetic_app* app);

/**
 * Free gk app (geom only).
 *
 * @param app App to release.
 */
void gkyl_gyrokinetic_app_release_geom(gkyl_gyrokinetic_app* app);

/**
 * Reset the CFL factor for the omega_H frequency.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param cfl_frac_omegaH New CFL factor to use for the omega_H mode.
 */
void gkyl_gyrokinetic_app_reset_cfl_frac_omegaH(gkyl_gyrokinetic_app* app, double tm,
  double cfl_frac_omegaH);

/**
 * Reset the df/dt multiplier operator for a given species.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param species_name Name of the species to reset.
 * @param fdot_mult_inp Input struct for the fdot_multiplier.
 */
void gkyl_gyrokinetic_app_reset_species_fdot_multiplier(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp);

/**
 * Reset the collisionless multiplier for a given species.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param species_name Name of the species to reset.
 * @param gkcls Input parameters for collisionless terms.
 */
void gkyl_gyrokinetic_app_reset_species_collisionless(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_collisionless gkcls_inp);

/**
 * Reset enforce_positivity for a given species.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param pos_inp Positivity input parameters.
 */
void gkyl_gyrokinetic_app_reset_species_positivity(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_positivity pos_inp);

/**
 * Reset the damping for a given species.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param species_name Name of the species to reset.
 * @param damping_inp Input parameters for damping terms.
 */
void
gkyl_gyrokinetic_app_reset_species_damping(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_damping damping_inp);

/**
 * Reset the field solver.
 *
 * @param app App object.
 * @param tm Time-stamp.
 * @param field_inp Input struct for the field object.
 */
void gkyl_gyrokinetic_app_reset_field(gkyl_gyrokinetic_app* app, double tm,
  struct gkyl_gyrokinetic_field field_inp);
