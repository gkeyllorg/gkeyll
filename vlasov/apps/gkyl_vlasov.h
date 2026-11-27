#pragma once

#include <gkyl_app.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_vlasov_comms.h>
#include <gkyl_mom_type.h>
#include <gkyl_vlasov_triad_geom.h>

#include <stdbool.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

// Parameters for projection
struct gkyl_vlasov_projection {
  enum gkyl_projection_id proj_id; // type of projection (see gkyl_eqn_type.h)
  enum gkyl_quad_type quad_type; // quadrature scheme to use: defaults to Gaussian

  union {
    struct {
      // pointer and context to initialization function
      void *ctx_func;
      void (*func)(double t, const double *xn, double *fout, void *ctx);
    };
    struct {
      // pointers and contexts to initialization functions for LTE distribution projection
      // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
      void *ctx_density;
      void (*density)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_V_drift;
      void (*V_drift)(double t, const double *xn, double *fout, void *ctx);
      void *ctx_temp;
      void (*temp)(double t, const double *xn, double *fout, void *ctx);

      // boolean if we are correcting all the moments or only density
      bool correct_all_moms;
      double iter_eps; // error tolerance for moment fixes (density is always exact)
      int max_iter; // maximum number of iteration
      bool use_last_converged; // use last iteration value regardless of convergence?
    };
  };
};

// Parameters for species collisions
struct gkyl_vlasov_collisions {
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)
  double nu_frac; // Rescales collision frequencies (default = 1).
  bool write_coll_diagnostics; // Whether to output diagnostics.

  // Function for computing self-collision frequency.
  void (*self_nu)(double t, const double *xn, double *fout, void *ctx);
  void *self_nu_ctx; // Context for self_nu.

  int num_cross_collisions; // Number of species to collide with.
  // Names of species to collide with. Kinetic species collide with kinetic species; PKPM species
  // and Euler fluids with gas_gamma = 5/3 collide with each other; a kinetic species with LBO
  // collisions collides with PKPM species (cross_nu given on both sides, explicit collisions of
  // the PKPM species).
  char collide_with[GKYL_MAX_SPECIES][128];
  // Functions for computing cross-collision frequencies (one for each num_cross_collisions).
  // Either all partners specify a function or none do (Spitzer cross frequencies).
  evalf_t cross_nu[GKYL_MAX_SPECIES];
  void *cross_nu_ctx[GKYL_MAX_SPECIES]; // Context for cross_nu.

  // Reference parameters for the Coulomb logarithm of computed (Spitzer) collision
  // frequencies. Unset entries default to 1, i.e. normalized units, as in the Lua
  // interface.
  double den_ref; // Reference density.
  double temp_ref; // Reference temperature.
  double hbar, eps0; // Planck's constant/2 pi, vacuum permittivity.

  // BGK collisions specific inputs
  bool fixed_temp_relax; // Are BGK collisions relaxing to a fixed input temperature?
  bool is_implicit; // Boolean for using implicit BGK collisions (backward Euler in each RK3 stage).
};

// Parameters for a BGK-form temperature source of a species,
//   S = nu_Q(x) (f_M(n, u, T_Q(x)) - f),
// with a prescribed temperature T_Q. The target Maxwellian has the density and the drift of the
// species, so the source conserves both and relaxes the temperature, d_t T = nu_Q (T_Q - T), and
// every other moment of the departure from the target at the rate nu_Q. The explicit form does
// not depend on the collision operator.
struct gkyl_vlasov_source_bgk {
  enum gkyl_source_bgk_id source_bgk_id; // GKYL_NO_SOURCE_BGK (default) or GKYL_SOURCE_BGK_HEATING

  void *rate_profile_ctx; // context for the rate function
  // pointer to the source rate nu_Q(x)
  void (*rate_profile)(double t, const double *xn, double *fout, void *ctx);
  void *temp_shape_ctx; // context for the temperature function
  // pointer to the source temperature T_Q(x)
  void (*temp_shape)(double t, const double *xn, double *fout, void *ctx);

  // One backward-Euler step in each SSP-RK3 stage, over the time weight of the stage, instead of
  // the explicit update (whose time step is limited by nu_Q). A kinetic species with implicit
  // BGK collisions takes the collisions and the source in the same step, with the collisions
  // relaxing to the temperature the source gives the species over the step.
  bool is_implicit;
};

// Parameters for species radiation
struct gkyl_vlasov_radiation {
  enum gkyl_vlasov_radiation_id radiation_id; // type of radiation (see gkyl_eqn_type.h)
  double t_cool; // Cooling time in radiation operator rad_force ~ -1/t_cool*drag
  double p0; // (four-) velocity to cool to.
};

// Parameters for species source
struct gkyl_vlasov_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // Are we writing out the source?
  bool evolve_source; // Are our sources time-dependent?
  int num_sources;

  double source_length; // required for boundary flux source
  // Boundary-flux source: the species whose wall losses set the rate (any species with absorbing
  // walls: kinetic, PKPM, or an Euler or isothermal Euler fluid; a PKPM species must name a PKPM
  // one).
  char source_species[128];

  // Sources using projection routine.
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  // Adaptive source input parameters.
  int num_cross_source; // Number of species that we are sourcing with.
  char source_with[GKYL_MAX_SPECIES][128]; // Names of species that we are using for cross sources.
  double source_with_v_thresh
    [GKYL_MAX_SPECIES]; // Threshold velocity if re-scaling density based on partial moments.
  double source_with_f_thresh
    [GKYL_MAX_SPECIES]; // Threshold f for accumulating partial moments for re-scaling density.
  bool source_with_upper_half
    [GKYL_MAX_SPECIES]; // Are you using the upper-half or lower-half plane for partial moments?
  int source_with_proj
    [GKYL_MAX_SPECIES]; // Which projection function is being used with this adaptive source?
  bool filter; // Are we filtering the rescaled density?
  int num_filters; // Are we filtering repeatedly?
};

// Parameters for boundary conditions
struct gkyl_vlasov_bc {
  enum gkyl_species_bc_type type;
  void *aux_ctx;
  void (*aux_profile)(double t, const double *xn, double *fout, void *ctx);
  double aux_parameter;
};

struct gkyl_vlasov_bcs {
  struct gkyl_vlasov_bc lower, upper;
};

// Source of a fluid species, added to the right-hand side of every stage.
//   GKYL_FUNC_SOURCE: profile gives the rate of every component of the state (num_equations values).
//   GKYL_PROJ_SOURCE (Euler and isothermal Euler): the fluid image of the LTE projections, with
//     density n_s, temperature T_s and drift u_s (three components; V_drift may be NULL):
//     [m n_s, m n_s u_s, m n_s |u_s|^2/2 + n_s T_s/(gas_gamma - 1)] (no energy row for an
//     isothermal fluid), the moments of the Maxwellian source of a kinetic species.
//   GKYL_BFLUX_SOURCE (Euler and isothermal Euler): that image times the particle loss rate of
//     source_species through its absorbing walls over source_length. With an n_s that integrates
//     to source_length over the domain, the source replaces the losses, as for a kinetic species.
// The profiles are sampled at the physical coordinates of a position map.
struct gkyl_vlasov_fluid_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // write the source rate (with the frames)

  void *ctx; // context for source function
  // function for computing source profile (GKYL_FUNC_SOURCE)
  void (*profile)(double t, const double *xn, double *aout, void *ctx);

  double source_length; // length of the boundary-flux source (GKYL_BFLUX_SOURCE)
  // Boundary-flux source: the species whose wall losses set the rate (any species with absorbing
  // walls: kinetic, PKPM, or an Euler or isothermal Euler fluid).
  char source_species[128];

  // LTE projections (proj_id GKYL_PROJ_VLASOV_LTE) of GKYL_PROJ_SOURCE and GKYL_BFLUX_SOURCE,
  // summed.
  int num_sources;
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];
};

// Parameters for fluid species advection
struct gkyl_vlasov_fluid_advection {
  void *velocity_ctx; // context for applied advection function
  // pointer to applied advection velocity function
  void (*velocity)(double t, const double *xn, double *aout, void *ctx);
};

// Parameters for fluid species diffusion
struct gkyl_vlasov_fluid_diffusion {
  double D; // constant diffusion coefficient
  int order; // integer for order of the diffusion (4 for grad^4, 6 for grad^6, default is grad^2)
  void *Dij_ctx; // context for applied diffusion function if using general diffusion tensor
  // pointer to applied diffusion function is using general diffusion tensor. The tensor
  // is evaluated at nodes (assumed continuous), and its cross terms need corner ghost
  // cells: a communicator passed to the app must be created with sync_corners = true.
  void (*Dij)(double t, const double *xn, double *Dout, void *ctx);
};

struct gkyl_vlasov_correct_inp {
  bool correct_all_moms; // Boolean if we are correcting all the moments or only density.
  double iter_eps; // Error tolerance for moment fixes (density is always exact).
  int max_iter; // Maximum number of iteration
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
    // *even if* the scheme fails to converge.
  bool output_f_lte; // Boolean for writing out f_lte (used for calculating transport coeff.).
};

struct vlasov_mapc2p_vel {
  void *mapc2p_vel_ctx; // context for mapc2p function for velocity space
  // pointer to mapc2p function for velocity space:
  // xc are the computational space coordinates and on output
  // xp are the corresponding physical space coordinates.
  void (*mapc2p_vel_func)(double t, const double *xc, double *xp, void *ctx);
};

struct vlasov_mapc2p_pos {
  void *mapc2p_pos_ctx; // context for mapc2p function for configuration space
  // pointer to (per-direction) mapc2p function for configuration space:
  // xc is the computational space coordinate and on output xp is the
  // corresponding physical space coordinate. NULL => identity map.
  void (*mapc2p_pos_func)(double t, const double *xc, double *xp, void *ctx);
};

// Parameters for Vlasov geometry.
struct gkyl_vlasov_geom {
  double spin_bh, mass_bh; // Black hole spin and mass.
  bool use_preset_geom; // Use the preset triad geometry given by triad_preset_geom_type.
  enum gkyl_triad_preset_geom_type
    triad_preset_geom_type; // geom type for preset geometries for triads
};

// Which evolved quantities a species owns: a kinetic distribution function, a
// fluid moment vector, or both.
enum gkyl_species_type {
  GKYL_SPECIES_VLASOV, // kinetic distribution only
  GKYL_SPECIES_FLUID, // fluid moments only
  GKYL_SPECIES_PKPM // both: a parallel-kinetic-perpendicular-moment species
};

// Parameters for the kinetic block of a Vlasov species: the velocity-space
// grid, distribution-function initial conditions, collisions, sources, and
// kinetic BCs. Declared as the 'kinetic' block of struct gkyl_vlasov_species.
struct gkyl_vlasov_kinetic_species {
  enum gkyl_model_id model_id; // Type of model
    // (e.g., SR, general geometry, see gkyl_eqn_type.h).

  double lower[3], upper[3]; // Lower, upper bounds of velocity-space.
  int cells[3]; // Velocity-space cells.

  // Velocity-space mapping in each velocity-space dimension
  struct vlasov_mapc2p_vel mapc2p_vel[GKYL_MAX_CDIM];

  // Initial conditions using projection routine.
  int num_init; // Number of initial condition functions.
  struct gkyl_vlasov_projection projection[GKYL_MAX_PROJ];

  // Held background (First_Fourier_Harmonic_Notes/background_density_gradient.tex): a projection,
  // in the format of those of the initial condition, of the state the species is held about. Once,
  // in gkyl_vlasov_app_apply_ic, the right-hand side of the species is evaluated on its background
  // (the other species at their initial conditions, the evolved field zero): its momentum moment
  // divided by the density is a body force, applied with app_accel in every stage, and what is left
  // is subtracted in every stage. The background is then a steady state of the discrete system, and
  // the force acts on everything that evolves about it. The initial condition may carry a
  // perturbation of the background. The background magnetic field belongs in the external field
  // (ext_em); a background at rest along the field (or unmagnetized) is held by the force alone,
  // while a gradient across the field needs the diamagnetic flow in the background. Non-relativistic
  // species without velocity or position maps. A PKPM species gives the eight normalized
  // distributions here (one GKYL_PROJ_FUNC projection) and its fluid background in the fluid block.
  struct gkyl_vlasov_projection background;

  bool is_static; // Set to true if species does not change in time.
  bool no_collisionless_terms; // Set to true to turn off collisionles terms.
  bool write_omega_cfl; // Whether to output dt diagnostic for the CFL constraint.
  bool write_cell_avg; // Also write the cell average of f.
  bool use_lo; // bool to determine if using low-order kernels for non-canonical Hamiltonian models.

  bool
    use_vierbein; // bool to determine if we are using vierbein input or by default using tangent vectors/ triads
  bool
    use_extended_hamil_def; // bool to determine if we are using the extended hamil definitions which includes potentials

  // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  double skip_cell_thresh;

  int num_diag_moments; // number of diagnostic moments
  enum gkyl_distribution_moments diag_moments[16]; // list of diagnostic moments

  void *app_accel_ctx; // Context for applied acceleration function.
  // Pointer to applied acceleration function.
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // Set to true if applied acceleration function is time dependent.

  void *cov_tangent_basis_ctx; // Context for covariant tangent basis function.
  // Pointer to covariant tangent basis function.
  void (*cov_tangent_basis)(double t, const double *xn, double *aout, void *ctx);

  void *triad_basis_ctx; // Context for triad basis function.
  // Pointer to triad basis function.
  void (*triad_basis)(double t, const double *xn, double *aout, void *ctx);

  void *triad_basis_gradient_ctx; // Context for triad basis gradient function.
  // Pointer to triad basis gradient function.
  void (*triad_basis_gradient)(double t, const double *xn, double *aout, void *ctx);

  void *vierbein_ctx; // Context for vierbein function.
  // Pointer to vierbein function.
  void (*vierbein)(double t, const double *xn, double *aout, void *ctx);

  void *vierbein_gradient_ctx; // Context for vierbein gradient function.
  // Pointer to vierbein gradient function.
  void (*vierbein_gradient)(double t, const double *xn, double *aout, void *ctx);

  void *hamil_ctx; // Context for general (phase space) Hamiltonian function.
  // Pointer to general (phase space) Hamilonian function.
  void (*hamil)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_ctx; // Context for spatial metric function (covariant).
  // Pointer to metric (covariant components) function.
  void (*h_ij)(double t, const double *xn, double *aout, void *ctx);

  void *h_ij_inv_ctx; // Context for spatial metric function (contravariant).
  // Pointer to metric (contravariant components) function.
  void (*h_ij_inv)(double t, const double *xn, double *aout, void *ctx);

  void *det_h_ctx; // Context for determinant of the spatial metric.
  // Pointer to the determinant of the spatial metric.
  void (*det_h)(double t, const double *xn, double *aout, void *ctx);

  void *background_flows_ctx; // Context for background flows in can-pb
  // Pointer to the background flows.
  void (*background_flows)(double t, const double *xn, double *aout, void *ctx);

  void *effective_potential_ctx; // Context for effective potential in can-pb
  // Pointer to the effective potential
  void (*effective_potential)(double t, const double *xn, double *aout, void *ctx);

  // Input quantities used by LTE (local thermodynamic equilibrium) projection
  // This projection operator is used by BGK collisions.
  struct gkyl_vlasov_correct_inp correct;

  // Collisions to include.
  struct gkyl_vlasov_collisions collisions;

  // Radiation to include.
  struct gkyl_vlasov_radiation radiation;

  // Source to include.
  struct gkyl_vlasov_source source;

  // Temperature source to include.
  struct gkyl_vlasov_source_bgk source_bgk;

  // PKPM species only.
  // Gyrophase-mixing closure of the first harmonic: rate, in units of the species' local |Omega|, at
  // which agyrotropy beyond what one harmonic can represent relaxes back to that bound (0: the
  // default 10; a negative value switches the closure off). The overshoot above the bound is about
  // (drive rate/|Omega|)/gyro_mix; a drive the default cannot hold calls for more Fourier harmonics.
  double gyro_mix;

  // Gauss-Legendre nodes per transverse direction of a configuration face (0 => default 5, at
  // least 3).
  int num_face_nodes;

  // Boundary conditions.
  struct gkyl_vlasov_bcs bcx, bcy, bcz;
};

// Parameters for the fluid block of a Vlasov species: the equation object,
// fluid initial condition, advection/diffusion, and fluid BCs. Declared as the
// 'fluid' block of struct gkyl_vlasov_species.
struct gkyl_vlasov_fluid_species {
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
  // Held background of an Euler or isothermal Euler fluid, or of the fluid block of a PKPM species,
  // in the format of init (see the background of the kinetic block).
  void *background_ctx; // context for the background function
  void (*background)(double t, const double *xn, double *fout, void *ctx);
  // Evaluate init and background at the nodes of the basis instead of projecting them: for p = 1
  // the state is continuous, so a smooth profile has no jumps at the faces (no numerical flux
  // dissipation on it, and for a PKPM species the change of frame at a face is the identity).
  bool eval_on_nodes;

  // Hyper-diffusion coefficient
  double nuHyp;

  // Equation object: Euler (gkyl_wv_euler), isothermal Euler (gkyl_wv_iso_euler; p = vt^2 rho,
  // state [rho, rho ux, rho uy, rho uz], not on the GPU), scalar advection, or a canonical
  // Poisson bracket fluid. Its Riemann solver is the numerical flux of an Euler or isothermal
  // Euler fluid.
  struct gkyl_wv_eqn *equation;
  // Limiters of an Euler or isothermal Euler fluid. The vertex positivity limiter (Zhang-Shu,
  // polynomial order 1) is always applied: after every stage and after the implicit fluid-EM source
  // steps of an update, the vertex density and internal energy are kept above positivity_eps times
  // the cell values, so the sound speed is real everywhere in the cell (the density only for an
  // isothermal fluid). It scales the slopes toward the cell averages, keeps the cell averages, and
  // acts only in cells with a vertex value below that floor. With limit_fluid the characteristic
  // (minmod) slope limiter is applied before it in every stage: the slopes of the characteristic
  // variables are bounded by limiter_fac times the differences of the neighboring cell averages,
  // which also flattens smooth extrema. At polynomial order 2 and above there is no positivity
  // limiter and the characteristic limiter is always applied.
  bool limit_fluid; // Apply the characteristic slope limiter at polynomial order 1 (default false).
  double limiter_fac; // Factor of the cell-average differences (0 => default 1/sqrt(3)).
  double positivity_eps; // Relative floor of the positivity limiter (0 => default 1e-2).

  // advection coupling (if scalar advection)
  // gkyl_eqn_type eqn_id = GKYL_EQN_ADVECTION
  struct gkyl_vlasov_fluid_advection advection;

  // source term
  struct gkyl_vlasov_fluid_source source;

  // diffusion coupling to include
  struct gkyl_vlasov_fluid_diffusion diffusion;

  // Collisions with other species: the momentum and energy exchange of the fluid states of an
  // Euler fluid with gas_gamma = 5/3, with PKPM species and with other such fluids, one
  // backward-Euler step in each SSP-RK3 stage (see the collisions of a PKPM species below; not
  // on the GPU). A fluid has no collision operator of its own: read are num_cross_collisions,
  // collide_with, cross_nu and cross_nu_ctx (without them the Spitzer frequencies), nu_frac and
  // the reference parameters. Every partner must list this species.
  struct gkyl_vlasov_collisions collisions;

  void *can_pb_n0_ctx; // context for background density function in canonical PB fluid systems
  // pointer to background density function in canonical PB fluid systems
  void (*can_pb_n0)(double t, const double *xn, double *aout, void *ctx);

  // Temperature source of an Euler fluid (GKYL_SOURCE_BGK_HEATING): the moments of the BGK-form
  // source of a kinetic species (see gkyl_vlasov_source_bgk), which relaxes the pressure toward
  // n T_Q at the rate nu_Q, d_t p = nu_Q (n T_Q - p), with the density and the momentum unchanged.
  // Explicit, in every stage (is_implicit is not available for a fluid).
  struct gkyl_vlasov_source_bgk source_bgk;

  // Applied acceleration [ax, ay, az] (external force per unit mass) of an Euler or isothermal
  // Euler fluid and of the fluid block of a PKPM species: rho a in the momentum and rho u.a in the
  // energy, in the right-hand side of every SSP-RK3 stage (also with the implicit fluid-EM
  // coupling, whose source step does not take it).
  void *app_accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*app_accel)(double t, const double *xn, double *aout, void *ctx);
  bool app_accel_evolve; // set to true if applied acceleration function is time dependent

  // Boundary conditions. An absorbing wall (GKYL_SPECIES_ABSORB) of an Euler or isothermal Euler
  // fluid is a boundary flux: that of the particles of the Maxwellian of the boundary cell (its
  // density, flow and temperature) that move toward the wall, none coming back, the fluid
  // counterpart of the kinetic flux of a kinetic or PKPM species against the empty ghost cells of
  // its absorbing walls. The wall losses are recorded for boundary-flux sources. The ghost cells
  // hold a copy of the boundary cell for the slope limiter.
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// A PKPM species (GKYL_SPECIES_PKPM) owns both blocks: the kinetic block holds the grid of the
// NORMALIZED velocity w = v_par/xi and the normalized distributions, the fluid block the state
// [rho, rho ux, rho uy, rho uz, E], E = rho |u|^2/2 + 3 p/2, whose flux is the moment of the
// kinetic flux (no equation object). It needs the Maxwell field (GKYL_FIELD_E_B), which defines
// b, and polynomial order 1 (either basis type: with kinetic species in the same run the tensor
// type gives those the hybrid phase-space basis, which PKPM species always use). A PKPM species
// has three velocity degrees of freedom, and so must the species it runs together with: kinetic
// species with three velocity dimensions (vdim = 3 describes them; the grid of w comes from the
// kinetic block of the PKPM species) and Euler fluids with gas_gamma = 5/3. For example
//   struct gkyl_vlasov_kinetic_species elc_kinetic = { .lower = ..., .cells = ..., ... };
//   struct gkyl_vlasov_fluid_species elc_fluid = { .init = ..., ... };
//   .species = {
//     { .name = "elc", .charge = -1.0, .mass = 1.0, .type = GKYL_SPECIES_PKPM,
//       .kinetic = elc_kinetic, .fluid = elc_fluid },
//   },
// Kinetic block of a PKPM species:
//   Initial condition using the projection routine: one GKYL_PROJ_FUNC projection of a function
//   returning the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)] of the normalized
//   distributions. F is a probability density in w (M0 = 1, M1 = 0) and G its
//   perpendicular-energy partner (M0(G) = T_perp/T) with M2(F) + 2 M0(G) = 3; e.g. a
//   bi-Maxwellian is F = exp(-w^2/(2 r))/sqrt(2 pi r), G = (T_perp/T) F with r = T_par/T,
//   T = (T_par + 2 T_perp)/3. The first-harmonic blocks [M0_perp, N_perp] must be perpendicular
//   to b, and int M0_perp dw must vanish; a gyrotropic start returns zeros for them.
//   Collisions to include (GKYL_LBO_COLLISIONS or GKYL_BGK_COLLISIONS). Without self_nu the
//   self-collision frequency is the Spitzer one computed from the local moments, and so are
//   the cross-collision frequencies without cross_nu,
//   nu_sr = nu_frac alpha_E (m_s + m_r) / (2 m_s) n_r / (vtsq_s + vtsq_r)^(3/2),
//   with the constant part alpha_E of Morse's energy-exchange coefficient at the reference
//   parameters (nu_ss = nu_frac alpha_E n_s / (2 vtsq_s)^(3/2) for the self frequency).
//   is_implicit (BGK): one backward-Euler step in each SSP-RK3 stage, over the time weight of
//   the stage, instead of the explicit update (with its CFL limit); not available on the GPU.
//   write_coll_diagnostics and fixed_temp_relax are not available (source_bgk is the
//   temperature source).
//   Source to include: GKYL_PROJ_SOURCE or GKYL_BFLUX_SOURCE with one GKYL_PROJ_VLASOV_LTE
//   projection, the profile [density, isotropic temperature, drift] injected into the
//   normalized kinetic equations (through the configuration pass, so the closures keep the
//   constraints exact) and into the fluid state [rho, M, E]. With GKYL_PROJ_SOURCE the profile
//   is added as given; with GKYL_BFLUX_SOURCE it is scaled every stage by the particle loss
//   rate of `source_species` through its absorbing walls divided by `source_length` (so a
//   density profile integrating to source_length replaces the losses exactly). `source_species`
//   must be a PKPM species: the rate is formed before the right-hand sides of the stage, when
//   only the losses of PKPM species are known (a kinetic species may name a PKPM one). The drift
//   enters the parallel Maxwellian through its component along b; a perpendicular drift
//   relative to the local flow is not carried by G. write_source writes the fluid image of the
//   profile at frame 0.
//   Temperature source to include, S = nu_Q(x) (f_M(n, u, T_Q(x)) - f) (see struct
//   gkyl_vlasov_source_bgk). In the normalized variables the fluid thermal energy relaxes to
//   (3/2) n T_Q at the rate nu_Q and the distributions relax to the unit Maxwellian at the rate
//   nu_Q T_Q/T, which relaxes every pressure component at the rate nu_Q of the source. Terms of
//   order (T_Q/T - 1) times the departure of the distributions from a Maxwellian are dropped:
//   moments of degree three decay at nu_Q T_Q/T instead of nu_Q (3 T_Q/T - 1)/2. The explicit
//   time step is limited by nu_Q and by nu_Q T_Q/T.
// Fluid block of a PKPM species:
//   Vertex positivity limiter of the fluid state [rho, M, E] (Zhang-Shu scaling of the
//   slopes toward the cell averages, always applied, after every forward-Euler stage, the
//   implicit coupling and the initial condition): the vertex density and the vertex internal
//   energy are kept above positivity_eps times the cell values, so the thermal speed is real
//   at the Gauss-Legendre nodes. The cells it acts on are recorded in
//   <name>-<species>_pkpm_positivity_<frame>.gkyl (cumulative counts per cell: [density,
//   internal energy, inadmissible average]) and the counts per integrated-moment interval
//   in the dynvec <name>-<species>-positivity-stat.gkyl ([cells scaled for the density, for
//   the internal energy, cells with an inadmissible average, limiter applications], summed
//   over ranks). Both hold the limiter of the initial condition and of the steps taken: a step
//   that is retaken with a smaller time step leaves no count.
//   diffusion: of the momentum (also on a position map), order 2 (grad^2, the default) or 4.
// Boundary conditions: each block carries its own and each is applied to its own state.
// The kinetic BC fills the ghost cells of the distributions: GKYL_SPECIES_COPY, GKYL_SPECIES_ABSORB
// (empty ghost distributions; the wall losses are recorded and T_perp/T of the ghost cell is the
// mirror image of the skin cell), GKYL_SPECIES_REFLECT (the field-aligned mirror image) or
// GKYL_SPECIES_FIXED_FUNC (the initial condition). The fluid BC fills the
// ghost cells of [rho, M, E], which set the frame (rho, u, xi) of the ghost cell: GKYL_SPECIES_COPY,
// GKYL_SPECIES_ABSORB (the mirror image of the skin cell, even in the wall direction, so the
// primitive variables and their gradients are one-sided at the wall), GKYL_SPECIES_REFLECT (the
// mirror image with the momentum reflected) or GKYL_SPECIES_FIXED_FUNC. The same type in both
// blocks is the natural pair: an absorbing wall records its losses, and a reflecting wall
// passes no mass and no energy.

// Parameters for one Vlasov species. The type selects which block is read
// (GKYL_SPECIES_VLASOV: 'kinetic', GKYL_SPECIES_FLUID: 'fluid', GKYL_SPECIES_PKPM:
// both, see above); the identity
// (name, charge, mass) is declared here and not in the blocks. Names must be
// unique. For example, a kinetic species and a fluid species are declared as
//   struct gkyl_vlasov_kinetic_species elc = { .lower = ..., .cells = ..., ... };
//   struct gkyl_vlasov_fluid_species ion = { .equation = eqn, .init = ..., ... };
//   .num_species = 2,
//   .species = {
//     { .name = "elc", .charge = -1.0, .mass = 1.0, .type = GKYL_SPECIES_VLASOV, .kinetic = elc },
//     { .name = "ion", .charge = 1.0, .mass = 1836.0, .type = GKYL_SPECIES_FLUID, .fluid = ion },
//   },
struct gkyl_vlasov_species {
  char name[128]; // Species name.
  double charge, mass; // Charge and mass.

  enum gkyl_species_type type; // Which aspect(s) this species owns.

  struct gkyl_vlasov_kinetic_species
    kinetic; // Kinetic block (velocity grid, projections, collisions, ...).
  struct gkyl_vlasov_fluid_species fluid; // Fluid block (equation, init, diffusion, ...).
};

// Parameter for electromagnetic fields.
struct gkyl_vlasov_field {
  // Type of field, which selects the solver: GKYL_FIELD_E_B (Maxwell, the default),
  // GKYL_FIELD_GR_D_B (GR Maxwell), GKYL_FIELD_PHI (Poisson) or GKYL_FIELD_NULL (no
  // field, as skip_field). See gkyl_eqn_type.h.
  enum gkyl_field_id field_id;
  bool is_static; // Set to true if field does not change in time.

  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  // Gravitational coupling 4 pi G in code units: Vlasov-Poisson solves the
  // self-gravity Poisson equation nabla^2 phi_g = alpha_g * sum_s m_s n_s and
  // accelerates every species by -grad(phi_g) when alpha_g > 0 (and some species
  // has mass). The electrostatic Poisson equation is solved when epsilon0 > 0
  // (and some species is charged); either, both, or neither may be present. In
  // normalized units with epsilon0 = 1 and a reference species of unit charge
  // and mass, alpha_g = 4 pi G epsilon0 m^2/q^2 is the dimensionless ratio of
  // the gravitational to the electrostatic coupling between two reference
  // particles.
  double alpha_g;
  // Screening of the electrostatic potential by a massive (dark) photon: with
  // mu_sq > 0 the Vlasov-Poisson electrostatic solve is the screened Poisson
  // (Helmholtz) equation -nabla^2 phi + mu_sq phi = rho_c/epsilon0, where mu =
  // m_A c/hbar is the inverse screening length (the photon's Compton wave
  // number), and the field energy diagnostic becomes int |grad phi|^2 + mu_sq phi^2.
  double mu_sq;
  // Correction speeds as a fraction of the speed of light for div(E)/div(B) errors. The div(E)
  // correction is for vacuum fields: it is rejected in a run with charged species.
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;
  double K_phi, K_psi; // Geometric source coefficients for electric/magnetic cleaning fields.

  void *ctx; // Context for initial condition init function for Vlasov-Maxwell.
  // Pointer to initialization function for Vlasov-Maxwell fields.
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  double limiter_fac; // Optional input parameter for adjusting diffusion in slope limiter.
  bool limit_em; // Optional input parameter for applying limiters to EM fields.

  // Are we using ghost currents to correct dE/dt = -J in 1x? The average current of every
  // species, kinetic or fluid, is removed from Ampere's law.
  bool use_ghost_current;
  bool use_geom_sources; // Include the geometric source terms of the GR Maxwell (D, B) field.

  // Vlasov-Maxwell boundary conditions.
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];

  // Vlasov-Poisson boundary conditions.
  struct gkyl_poisson_bc poisson_bcs;

  void *sigma_ctx; // context for the resistive layer
  // Resistive layer sigma(x) >= 0: E is damped at this rate (B is not). sigma is sampled at the
  // Gauss-Legendre points of each cell; the damping is a quadrature projection on the stage
  // right-hand side, or exact at each node of the source solve with the fluid-EM coupling.
  void (*sigma)(double t, const double *xn, double *fout, void *ctx);

  void *ext_em_ctx; // Context for external electromagnetic fields function (E,B).
  // Pointer to external electromagnetic fields function.
  void (*ext_em)(double t, const double *xn, double *ext_em_out, void *ctx);
  bool ext_em_evolve; // Set to true if external electromagnetic fields are time dependent.

  void *app_current_ctx; // Context for applied currents function.
  // Pointer to applied currents function.
  void (*app_current)(double t, const double *xn, double *app_current_out, void *ctx);
  bool app_current_evolve; // Set to true if applied currents are time dependent.

  void *external_potentials_ctx; // Context for external (phi,A) potentials.
  // Pointer to function defining external potentials (phi,A).
  void (*external_potentials)(double t, const double *xn, double *ext_pot, void *ctx);
  bool external_potentials_evolve; // Set to true if external potentials are time dependent.
};

// Largest number of moments a neural network takes as inputs or as outputs.
#define GKYL_VLASOV_NN_MAX_MOMS 16

// Neural-network options of a run, read by the driver of the run (the Lua wrapper): training
// one network per species on the diagnostic moments of the simulation (gkyl_vlasov_app_train,
// gkyl_vlasov_app_write_nn) and testing a trained one against them (gkyl_vlasov_app_test). All
// zeros, the default, means no network is used.
struct gkyl_vlasov_nn {
  bool train_nn; // Train neural network on simulation data?
  bool train_ab_initio; // Train a new network (otherwise the one read from train_nn_file)?
  int nn_width; // Number of neurons to use per layer.
  int nn_depth; // Number of layers to use.
  char train_nn_file[256]; // File path of neural network to train.
  int num_trains; // Number of times to train neural network.
  int num_nn_writes; // Number of times to write out neural network.
  int num_input_moms; // Number of "input" moments to train on.
  int input_moms[GKYL_VLASOV_NN_MAX_MOMS]; // "Input" moments to train on.
  int num_output_moms; // Number of "output" moments to train on.
  int output_moms[GKYL_VLASOV_NN_MAX_MOMS]; // "Output" moments to train on.
  bool test_nn; // Test neural network on simulation data?
  char test_nn_file[256]; // File path of neural network to test.
  int num_tests; // Number of times to test neural network.
};

// Top-level app parameters
struct gkyl_vm {
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

  // Per-direction non-uniform configuration-space mapping (C^0 piecewise
  // linear). Each direction with a NULL func is the identity map. This is a
  // diagonal coordinate stretch and is independent of the general curvilinear
  // mapc2p above (which feeds the wave_geom object).
  struct vlasov_mapc2p_pos mapc2p_pos[GKYL_MAX_CDIM];

  double cfl_frac; // CFL fraction to use (default 1.0)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  struct gkyl_vlasov_geom geom; // geom object

  // Species list: kinetic and fluid species in any order, at most
  // GKYL_MAX_SPECIES in total. Declaration order sets the index used by the
  // per-species app functions (files are named by the species name). Kinetic-
  // only diagnostics (distribution moments, L2 norm) skip fluid species.
  int num_species; // total number of species (kinetic + fluid)
  struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // species objects

  bool skip_field; // No field; the field input must then be left unset.
  struct gkyl_vlasov_field field; // Field input; its field_id selects the solver.

  // Momentum-carrying (Euler) fluid species with the Maxwell field: the default is a staggered
  // (Strang) operator split with implicit fluid-EM coupling. Between updates the fluids and E lag
  // by half a source step, output and diagnostics see them synchronized, and the field frames
  // store the steps that make a restart continue the uninterrupted run. Kinetic species may be
  // present (hybrid runs, polynomial order 1): the source step takes their current as known, the
  // transport step moves them with its E, and the step resolves their plasma frequency
  // (cfl_frac / omega_p). The applied accelerations (app_accel) of all species act in the
  // transport step, with the pressure gradients they may balance, and not in the source step.
  // With use_explicit_sources the Lorentz force is in the fluid RHS and the fluid currents on the
  // RHS of Ampere's law in every SSP-RK3 stage, and in each cell the time-step resolves the larger
  // of the rate of the hyperbolic terms and the rate of the source of each coupled fluid (see
  // gkyl_dg_calc_fluid_vars_em_source): omega_p + |Omega_c|, or its acceleration over its
  // thermal speed.
  // With the Poisson field (GKYL_FIELD_PHI) the same split is a kick-drift-kick scheme: the
  // source step kicks the fluid momenta with -q/m grad(phi + phi_ext) and -grad(phi_g), solved
  // from the densities of all species at its start (which it does not change), implicit in the
  // external magnetic field, and the transport step moves the fluids with no force of the
  // potentials (kinetic species solve for them in each stage). The coupling of the densities to
  // the potentials is explicit, so the step resolves the plasma frequency (and the Jeans rate) of
  // every species, cfl_frac / sqrt(sum_s q_s^2 n_s/(epsilon0 m_s) + alpha_g m_s n_s). A PKPM
  // species needs a static external magnetic field (ext_em), which sets b; the fluid aspects
  // take the scalar external potential and not its vector part. On a position map (PKPM and
  // kinetic species) the kick takes the physical gradient, with the per-cell J_d of the map that the
  // metric of the Poisson solve uses. use_explicit_sources is not supported with the Poisson field.
  bool use_explicit_sources; // Explicit fluid-Maxwell coupling in the SSP-RK3 stages.

  struct gkyl_vlasov_nn nn; // Neural-network options of PKPM species (all zeros: no network).

  struct gkyl_app_parallelism_inp parallelism; // Parallelism-related inputs.
};

// Simulation statistics
struct gkyl_vlasov_stat {
  bool use_gpu; // did this sim use GPU?

  long nup; // calls to update
  long nfeuler; // calls to forward-Euler method

  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s
  long nfl_em_retake; // number of repeated implicit fluid-EM coupling steps

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure

  double total_tm; // time for simulation (not including ICs)
  double rk3_tm; // time for SSP RK3 step
  double fl_em_tm; // time for implicit fluid-EM coupling step (with its synchronization for output)
  double init_species_tm; // time to initialize all species
  double init_fluid_species_tm; // time to initialize all fluid species
  double init_field_tm; // time to initialize fields

  // time to compute species collisionless RHS (PKPM species: the streaming and velocity passes;
  // the other terms of its kinetic equations have the timers species_pkpm_*_tm below)
  double species_rhs_tm;
  double fluid_species_rhs_tm; // time to compute fluid species RHS
  double fluid_species_vars_tm; // time to compute fluid variables (flow velocity and pressure)

  long
    n_iter_corr[GKYL_MAX_SPECIES]; // total number of iterations used to correct species LTE projection
  long num_corr
    [GKYL_MAX_SPECIES]; // total number of times correction updater for species LTE projection is called
  double species_coll_mom_tm; // time needed to compute various moments needed in LBO
  double species_lbo_coll_drag_tm[GKYL_MAX_SPECIES]; // time to compute LBO drag terms
  double species_lbo_coll_diff_tm[GKYL_MAX_SPECIES]; // time to compute LBO diffusion terms
  double species_coll_tm; // total time for collision updater (excluded moments)

  double species_rad_tm; // total time for radiation updater

  double species_lte_tm; // time needed to compute the lte equilibrium

  double species_bc_tm; // time to compute species BCs
  double fluid_species_bc_tm; // time to compute fluid species BCs
  double field_bc_tm; // time to compute field BCs

  double field_rhs_tm; // time to compute field RHS
  double current_tm; // time to compute currents and accumulation
  // time to compute the configuration-space variables of a stage of a PKPM species: the
  // moments of the normalized distributions, the flow and thermal speed, the kinetic coupling
  // moments, the constraints and the gradients for the sources
  double species_pkpm_vars_tm;
  // time to compute the configuration-space face fluxes of a stage of a PKPM species (kinetic
  // and fluid, with their rates) and the in-cell gradients of u, rho and xi
  double species_pkpm_conf_flux_tm;
  // time for the first harmonic of a PKPM species: the projection perpendicular to b, the
  // coefficients, w fluxes and sources of the blocks and the integrating factor
  double species_pkpm_harm_tm;
  // time to compute the closures of a PKPM species (moments of the partial RHS and the solve)
  double species_pkpm_closure_tm;
  double species_pkpm_vel_flux_tm; // time to compute the velocity-space face fluxes (PKPM species)
  double field_em_vars_tm; // time to compute b = B/|B| and its gradient (PKPM species)

  double species_omega_cfl_tm; // time spent in all-reduce for omega-cfl for species
  double species_io_tm; // Time to write the species distribution.
  double species_diag_calc_tm; // Time to compute species diagnostics.
  double species_diag_io_tm; // Time to write species diagnostics.

  double field_omega_cfl_tm; // time spent in all-reduce for omega-cfl for fields
  double field_io_tm; // Time to write the fields.
  double field_diag_calc_tm; // Time to compute field diagnostics.
  double field_diag_io_tm; // Time to write field diagnostics.

  long n_species_omega_cfl; // number of times CFL-omega all-reduce is called for species
  long n_field_omega_cfl; // number of times CFL-omega all-reduce is called for fields
  long n_mom; // total number of calls to moment updater routines
  long n_diag; // total number of calls to diagnostics
  long n_io; // number of calls to IO
  long n_diag_io; // number of calls to IO for diagnostics
  long n_field_diag; // total number of calls to diagnostics for field
  long n_field_io; // number of calls to IO for field
  long n_field_diag_io; // number of calls to IO for field diagnostics
};

// Object representing Vlasov app
typedef struct gkyl_vlasov_app gkyl_vlasov_app;

/**
 * Construct a new Vlasov app.
 *
 * @param vm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov app object.
 */
gkyl_vlasov_app *gkyl_vlasov_app_new(struct gkyl_vm *vm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_app_apply_ic(gkyl_vlasov_app *app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app *app, double t0);

/**
 * Initialize a species by projecting initial conditions on basis functions.
 * Species are indexed in the order they appear in gkyl_vm's species[] list
 * (kinetic and fluid alike); see gkyl_vlasov_app_find_species to resolve a
 * species by name.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app *app, int sidx, double t0);

/**
 * Find a species by name.
 *
 * @param app App object.
 * @param nm Species name.
 * @return Index of the species (as used by the per-species app functions), -1 if not found.
 */
int gkyl_vlasov_app_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Initialize field from file
 *
 * @param app App object
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_file_field(
  gkyl_vlasov_app *app, const char *fname
);

/**
 * Initialize Vlasov species from file
 *
 * @param app App object
 * @param sidx Index of species (see gkyl_vlasov_app_apply_ic_species).
 * @param fname file to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_file_species(
  gkyl_vlasov_app *app, int sidx, const char *fname
);

/**
 * Initialize field from frame
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_frame_field(gkyl_vlasov_app *app, int frame);

/**
 * Initialize Vlasov species from frame
 *
 * @param app App object
 * @param sidx Index of species (see gkyl_vlasov_app_apply_ic_species).
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_from_frame_species(
  gkyl_vlasov_app *app, int sidx, int frame
);

/**
 * Initialize the Vlasov app from a specific frame. With the implicit fluid-EM coupling the
 * frame's restart files of the field and the fluids are read as well, so that the run continues
 * as the run that wrote the frame.
 *
 * @param app App object
 * @param frame frame to read
 */
struct gkyl_app_restart_status gkyl_vlasov_app_read_from_frame(gkyl_vlasov_app *app, int frame);

/**
 * Calculate integrated diagnostic moments.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app *app, double tm);

/**
 * Calculate integrated L2 norm of the distribution function, f^2.
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app *app, double tm);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app *app, double tm);

/**
 * Write field and species data to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write field data to file. With the implicit fluid-EM coupling the field the stepper holds
 * between updates is written too (<name>-field_restart_<frame>), for a restart.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_field(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write species data to file. With the implicit fluid-EM coupling the fluid state the stepper
 * holds between updates is written too (<name>-<species>_restart_<frame>), for a restart, and so
 * are the ghost cells of a kinetic species at an emitting wall
 * (<name>-<species>_bc_{lo,up}_restart_<frame>).
 *
 * @param app App object.
 * @param sidx Index of species to write (see gkyl_vlasov_app_apply_ic_species).
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species(gkyl_vlasov_app *app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 *
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_mom(gkyl_vlasov_app *app, double tm, int frame);

/**
 * Write integrated diagnostic moments for all species (kinetic and fluid) to file.
 * Integrated moments are appended to the same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app);

/**
 * Write integrated L2 norm of the species distribution function to file. Integrated
 * L2 norm is appended to the same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_write_integrated_L2_f(gkyl_vlasov_app *app);

/**
 * Write integrated correct lte status of the species distribution function to file. Correct
 * lte status is appended to the same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_write_lte_corr_status(gkyl_vlasov_app *app);

/**
 * Write field energy to file. Field energy data is appended to the
 * same file.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app *app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_stat_write(gkyl_vlasov_app *app);

/**
 * Train neural network on PKPM moments data for each PKPM species.
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
void gkyl_vlasov_app_train(
  gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data,
  struct gkyl_kn_vec *output_data
);

/**
 * Train neural network on PKPM moments data for a particular PKPM species.
 * 
 * @param app App object.
 * @param sidx Index of the PKPM species to train on.
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
void gkyl_vlasov_app_train_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data, struct gkyl_kn_vec *output_data
);

/**
 * Write out PKPM moments neural network for each PKPM species.
 * 
 * @param app App object.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 */
void gkyl_vlasov_app_write_nn(
  gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann
);

/**
 * Write out PKPM moments neural network for a particular PKPM species.
 * 
 * @param app App object.
 * @param sidx Index of the PKPM species to write out.
 * @param tm Time-stamp.
 * @param frame Frame number.
 * @param ann Neural network architecture.
 */
void gkyl_vlasov_app_write_nn_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann
);

/**
 * Test neural network on PKPM moments data for each PKPM species, and write moments to a file.
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
void gkyl_vlasov_app_test(
  gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data_real,
  struct gkyl_kn_vec *output_data_real, struct gkyl_kn_vec *output_data_predicted
);

/**
 * Test neural network on PKPM moments data for a particular PKPM species, and write moments to a file.
 * 
 * @param app App object.
 * @param sidx Index of the PKPM species to test on.
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
void gkyl_vlasov_app_test_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data_real, struct gkyl_kn_vec *output_data_real,
  struct gkyl_kn_vec *output_data_predicted
);

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
void gkyl_vlasov_app_cout(const gkyl_vlasov_app *app, FILE *fp, const char *fmt, ...);

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
struct gkyl_update_status gkyl_vlasov_update(gkyl_vlasov_app *app, double dt);

/**
 * Return simulation statistics.
 *
 * @return Return statistics object.
 */
struct gkyl_vlasov_stat gkyl_vlasov_app_stat(gkyl_vlasov_app *app);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app *app);
