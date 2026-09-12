#pragma once

#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_efit.h>
#include <gkyl_evalf_def.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_position_map.h>


typedef struct gk_geometry gk_geometry;

// Some cumulative statistics
struct gkyl_tok_geo_stat {
  long nquad_cont_calls; // num calls from quadrature
  long nroot_cont_calls; // num calls from root-finder
};  

typedef void (*plate_func)(double s, double* RZ);

// Explicit material target on the unchanged EQDSK limiter. Segment i joins
// vertex i to (i+1)%limiter_n. Supply a connected, nonclosed sequence in either
// file order or its reverse; no device-dependent divertor labels are inferred.
struct gkyl_tok_geo_wall_target {
  int num_segments;
  const int *segments;
};

// Type of flux surface
enum gkyl_tok_geo_type {
  // Full blocks to be used as stand alone simulations
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT, // Full Outboard SOL of double-null (DN) configuration
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN, // Full Inboard SOL of DN configuration
  GKYL_GEOMETRY_TOKAMAK_LSN_SOL, // Full SOL of a lower single-null (LSN) configuration
  GKYL_GEOMETRY_TOKAMAK_USN_UP, // Full SOL of an upper single-null (USN) configuration -- not yet implemented
  GKYL_GEOMETRY_TOKAMAK_CORE, // Full core

  // 6 SOL Block Types for DN multi-block simulations
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO,  // Section of outboard SOL below lower xpt
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID, // Section of outboard SOL between xpts
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP,  // Section of outboard SOL above upper xpt
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO,   // Section of inboard SOL below lower xpt
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID,  // Section of inboard SOL between xpts
  GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP,   // Section of inboard SOL above upper xpt 
  
  // 3 SOL Block Types for LSN multi-block simulations
  GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO, // Outboard divertor leg of LSN
  GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID, // Middle portion of LSN SOL between X-points
  GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP, // Inboard divertor leg of LSN

  // PF Block types that can be used with SN or DN configurations in multi-block simulations
  GKYL_GEOMETRY_TOKAMAK_PF_UP_L, // Left half of Private flux region at top (inboard upper plate to upper xpt)
  GKYL_GEOMETRY_TOKAMAK_PF_UP_R, // Right half of Private flux region at top (upper xpt to outboard upper plate)
  GKYL_GEOMETRY_TOKAMAK_PF_LO_L, // Left half of Private flux region at bottom (lower xpt to inboard lower plate)
  GKYL_GEOMETRY_TOKAMAK_PF_LO_R, // Right half of Private flux region at bottom (outboard lower plate to lower xpt)

  // Core Block types that can be used with SN or DN configurations in multi-block simulations 
  GKYL_GEOMETRY_TOKAMAK_CORE_L, // Left half of core (upper to lower xpt)
  GKYL_GEOMETRY_TOKAMAK_CORE_R, // Right half of core (lower to upper xpt)

  GKYL_GEOMETRY_TOKAMAK_IWL, // Inner Wall Limited
};  

struct gkyl_tok_geo {
  struct gkyl_efit* efit;

  struct gkyl_rect_grid rzgrid; // RZ grid on which psi(R,Z) is defined
  struct gkyl_range rzlocal; // local range over which psiRZ is defined
  struct gkyl_range rzlocal_ext; // extended range
  struct gkyl_rect_grid rzgrid_cubic; // RZ grid on which the cubic rep of psi(R,Z) is defined
  struct gkyl_range rzlocal_cubic; // local range over which the cubic rep of psiRZ is defined
  struct gkyl_range rzlocal_cubic_ext; // extended range
  struct gkyl_basis rzbasis; // basis functions for R,Z grid
  struct gkyl_basis rzbasis_cubic; // cubic basis functions for R,Z grid
  int num_rzbasis; // number of basis functions in RZ
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_array *psiRZ_cubic; // cubic psi(R,Z) DG representation
  // Conservative enclosure of psi over each RZ cell (2 components, lo and
  // hi). Lets the R-root scan skip cells that cannot hold the requested
  // level set. NULL when the enclosure does not apply.
  struct gkyl_array *psi_cell_bounds;
  // Coarse enclosure over runs of psi_block_size neighbouring R cells, so a
  // whole run can be rejected with one test. NULL when unavailable.
  double *psi_block_bounds;
  int psi_block_size;
  int psi_num_blocks;
  struct gkyl_basis_ops_evalf *evf ; // wrapper for cubic evaluation
                   
  struct gkyl_rect_grid fgrid; // flux grid for fpol
  struct gkyl_range frange; // flux range
  struct gkyl_range frange_ext; // extended range
  struct gkyl_basis fbasis; // psi basis for fpol
  const struct gkyl_array *fpoldg; // fpol(psi) dg rep
  const struct gkyl_array *fpolprimedg; // fpol'(psi) dg rep
  const struct gkyl_array *qdg; // q(psi) dg rep                                   

  double sibry; // psi of separatrix as given by EFIT
  double psisep; // psi of separatrix as calculated from the DG psi(R,Z)
  double zmaxis; // z of magnetic axis
  // rleft : If you are in a circular kind of region (like a single null SOL or the core) and
  // theta is greater than theta of the upper turning point (so we have already traced the
  // entire right half of the surface), nodes/roots with R closest to rleft will be chosen.
  double rleft; 
  double rright;
  // rmin : No root with R < rmin will ever be chosen. rmin is interpreted as maybe a machine
  // boundary. So, no node will be placed at r < rmin.
  double rmin;
  double rmax;

  // Flag and functions to specify the plate location/shape in RZ coordinates
  // The functions should specify R(s) and Z(s) on the plate where s is a parameter \in [0,1]
  //
  // RESOLUTION FLOOR. The strike point is found by sampling the plate at 512
  // uniform values of s and bracketing sign changes of psi(plate(s)) - psi0.
  // Two consequences a plate author needs to know, both measured (gate 4):
  //
  //   * Two crossings of the same flux surface separated by less than the
  //     sample spacing, 1/512 of the plate parameter, are not resolved.
  //   * A crossing that is TANGENT -- the plate grazes the surface without
  //     passing through it -- is not found at all, because there is no sign
  //     change to bracket.
  //
  // In both cases the library REFUSES the domain (TOK_GEO_ACTIVE_PLATE_FAILED)
  // rather than guessing a strike point, so the failure is conservative. Shape
  // a plate with features coarser than 1/512 of its length, and avoid grazing
  // incidence, or expect a rejection you have to diagnose.
  // For single null, the "lower" plate is the outboard plate and the "upper plate" is the inboard plate
  // For IWL, when s=0, the plate function must return the coordinates of the corner of the limiter plate
  // at the inboard midplane which lies on the LCFS.
  bool plate_spec;
  // Continue an insufficient plate on the actual EQDSK limiter. Default off.
  // Environment EXTEND_TO_LIMITER=0 or 1 overrides this input.
  bool extend_to_limiter;
  plate_func plate_func_lower;
  plate_func plate_func_upper;
  struct gkyl_tok_geo_wall_target divertor_wall[2]; // owned copies; callback lower/upper slots

  struct { int max_iter; double eps; } root_param;
  struct { int max_level; double eps; } quad_param;

  bool inexact_roots; // If true we will allow approximate roots when no root is found
  bool use_cubics; // If true will use the cubic rep of psi rather than the quadratic representation
  bool use_hyperbolic_numbers; // If true will use the hyperbolic numbers to do cubic root finding (much faster)

  // pointer to root finder (depends on polyorder)
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z,
    double xc[2], double dx[2]);

  double (*calc_grad_psi)(const double *psih, const double eta[2], const double dx[2]);

  struct gkyl_tok_geo_stat stat; 
  struct gkyl_array* mc2p_nodal_fd;
  struct gkyl_range* nrange;
  double* dzc;
};



// Inputs to create a new GK geometry creation object

enum gkyl_tok_geo_xpt_seam_trial_failure {
  GKYL_XPT_SEAM_TRIAL_OK = 0,
  GKYL_XPT_SEAM_TRIAL_INVALID_PARAMETER,
  GKYL_XPT_SEAM_TRIAL_CONTOUR,
  GKYL_XPT_SEAM_TRIAL_BRANCH,
  GKYL_XPT_SEAM_TRIAL_TRACE_ORDERING,
  GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP,
  GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN,
  GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN,
  GKYL_XPT_SEAM_TRIAL_REMOTE_FAILURE,
};

// Status shared by one diagnostic optimizer trial and the mapping/metric
// builders. Applications should leave the corresponding trial pointer unset.
struct gkyl_tok_geo_xpt_seam_trial_status {
  bool contour_valid;
  bool branch_valid;
  bool trace_ordering_valid;
  bool finite_map_valid;
  bool cell_jacobian_valid;
  bool jacobian_valid;
  int jacobian_sign;
  int first_failure_reason;
  double max_realized_displacement;
  double min_cell_jacobian_margin;
};


// Inputs to create geometry for a specific computational grid
struct gkyl_tok_geo_grid_inp {
  // Experimental internal view of an actual radial neighbor. Owned by the
  // multiblock declaration; used only during geometry construction.
  const struct gkyl_tok_geo_grid_inp *shared_theta_peer;
  int shared_theta_radial_edge;
  bool shared_theta_reverse;
  struct gkyl_rect_grid cgrid;
  struct gkyl_basis cbasis;
  enum gkyl_tok_geo_type ftype; // type of geometry
  bool half_domain; // For use in double null simulations
                    // If true, will set the domain to be the lower
                    // half of the tokamak (below Z=0)
  
  double rclose; // closest R to region of interest to discriminate
  double rleft; // closest R to inboard SOL
  double rright; // closest R to outboard SOL
  double rmin, rmax; // Minimum and Maximum R of the machine
  double zmin, zmax; // extents of Z for integration
  double zmin_left, zmin_right; // for lower single null and PF cases diff b/t in and outboard side
  double zmax_left, zmax_right; // for upper single null and PF cases diff b/t in and outboard side

  // Specifications for divertor plate
  bool plate_spec;
  // Continue an insufficient plate on the actual EQDSK limiter. Default off.
  // Environment EXTEND_TO_LIMITER=0 or 1 overrides this input.
  bool extend_to_limiter; // whether a shape function is provided for divertor plates
  // With extension enabled, explicit wall arcs replace the corresponding
  // material plate target. Slot 0 is plate_func_lower, slot 1 plate_func_upper
  // (LSN: outboard/inboard, not geometric lower/upper). Zero entries preserve
  // the original plate policy. Inputs are copied. Native material-face maps
  // follow the wall bends; strike-point chords only visualize that curved face.
  struct gkyl_tok_geo_wall_target divertor_wall[2];

  // Explicit acknowledgement that this equilibrium supplies NO usable vessel
  // outline, so geometry is built without wall enforcement. Required because
  // silence must never disable a safety constraint: an equilibrium whose
  // outline is absent (limiter_status 0) or degenerate (status 2) and which
  // carries no acknowledgement is REJECTED. Every block built under this
  // acknowledgement reports TOK_GEO_WALL_NOT_ENFORCED, so an unenforced wall is
  // always visible in the run log.
  //
  // The acknowledgement is checked, not trusted, and is refused two ways:
  //   * a usable outline exists  -> contradictory declaration, rejected;
  //   * the record is MALFORMED (status -1) -> unreadable or non-finite data is
  //     an input error to fix, never something to declare away.
  // It grants no permission to leave a wall that does exist.
  bool no_vessel_outline;
  plate_func plate_func_lower; // lower plate specification. Gives R,Z in terms of s \in [0,1]
  plate_func plate_func_upper; // upper plate specification. Gives R,Z in terms of s \in [0,1]
                               // In a lower single null "lower" is the outer divertor and
                               // "upper" is the inner divertor

  bool inexact_roots; // If true we will allow approximate roots when no root is found
  bool use_cubics; // If true will use the cubic rep of psi rather than the quadratic representation
  bool use_hyperbolic_numbers; // If true will use the hyperbolic numbers to do cubic root finding (much faster)
  bool straight_xpt_ray; // Align supported half-domain, full-domain double-null,
                         // and lower-single-null block interfaces on straight,
                         // flux-surface-intersecting rays from their X points to
                         // the requested far radial surface.
  bool straight_core_xpt_ray; // Align the CORE_L/CORE_R interface on a straight,
                              // flux-surface-intersecting ray from the lower X-point
                              // to the nearest point on the innermost core surface.
                              // Retained for backward compatibility; straight_xpt_ray
                              // enables the same construction for every supported region.
  bool relaxed_xpt_seam; // Route supported half-domain X-point seams through
                         // the delta-s parameterization. The initial mode
                         // enforces delta_s=0 and is exactly the straight ray.
  bool relaxed_xpt_seam_sweep; // Master enable for a nonzero delta-s: with
                               // this false, relaxed_xpt_seam_delta_s_coeff
                               // is ignored and the seam stays exactly the
                               // straight ray regardless of its value. Used
                               // both for diagnostic-only sweeps/trials and,
                               // once a coefficient has been selected and
                               // guard-checked, for the applied production
                               // choice -- it does not by itself imply
                               // either.
  double relaxed_xpt_seam_delta_s_coeff; // Peak B1 coefficient [m].
  double relaxed_xpt_seam_delta_s_bound; // Hard displacement bound [m].
  bool relaxed_xpt_seam_optimize; // Run the diagnostic bounded optimizer.
  bool relaxed_xpt_seam_optimizer_trial; // Internal trial marker.
  struct gkyl_tok_geo_xpt_seam_trial_status
    *relaxed_xpt_seam_trial_status; // Internal trial status; otherwise NULL.

  // Parameters for root finder: leave unset to use defaults
  struct {
    int max_iter; // typically 20
    double eps; // typically 1e-10
  } root_param;
  // Parameters for nmumerical quadrature: leave unset to use default
  struct {
    int max_levels; // typically 6-7    
    double eps; // typically 1e-10
  } quad_param;
};

// A usable vessel outline: enough finite vertices to bound a region. This is the
// single predicate behind every wall decision; the enforcement sites must not
// re-derive it, or they will disagree the way the reader and its consumers did.
static inline bool
gkyl_tok_wall_usable(const struct gkyl_efit *e)
{
  return e && e->limiter_status == 1 && e->limiter_n >= 3 && e->limiter_R && e->limiter_Z;
}

enum gkyl_tok_wall_policy {
  GKYL_TOK_WALL_ENFORCE = 0,          // usable outline: enforce containment
  GKYL_TOK_WALL_NOT_ENFORCED,         // no usable outline, acknowledged: build unenforced
  GKYL_TOK_WALL_REJECT_UNDECLARED,    // no usable outline and no acknowledgement
  GKYL_TOK_WALL_REJECT_CONTRADICTED,  // acknowledged, yet a usable outline exists
  GKYL_TOK_WALL_REJECT_MALFORMED,     // unreadable/non-finite record: never declarable
};

// Decide the vessel-outline policy for one block. Pure: the caller reports, with
// its own context. Absence of a wall is deliberately NOT treated as being
// outside one -- that conflation is what made wall-less equilibria unbuildable.
static inline enum gkyl_tok_wall_policy
gkyl_tok_wall_policy_for(const struct gkyl_tok_geo_grid_inp *inp, const struct gkyl_efit *e)
{
  bool usable = gkyl_tok_wall_usable(e);
  if (inp && inp->no_vessel_outline) {
    if (usable) return GKYL_TOK_WALL_REJECT_CONTRADICTED;
    if (e && e->limiter_status == -1) return GKYL_TOK_WALL_REJECT_MALFORMED;
    return GKYL_TOK_WALL_NOT_ENFORCED;
  }
  return usable ? GKYL_TOK_WALL_ENFORCE : GKYL_TOK_WALL_REJECT_UNDECLARED;
}

// Name for diagnostics, so every site reports the same reason string.
static inline const char *
gkyl_tok_wall_policy_reason(enum gkyl_tok_wall_policy p)
{
  switch (p) {
    case GKYL_TOK_WALL_ENFORCE: return "enforced";
    case GKYL_TOK_WALL_NOT_ENFORCED: return "acknowledged_absent_outline";
    case GKYL_TOK_WALL_REJECT_UNDECLARED: return "outline_unusable_and_undeclared";
    case GKYL_TOK_WALL_REJECT_CONTRADICTED: return "declared_absent_but_outline_usable";
    case GKYL_TOK_WALL_REJECT_MALFORMED: return "outline_malformed_not_declarable";
  }
  return "unknown";
}


/**
 * Create new updater to compute the geometry needed in GK
 * simulations.
 *
 * @param efit_inp Input parameters related to EFIT data
 * @param grid_inp Input parameters related to computational grid
 */
struct gkyl_tok_geo *gkyl_tok_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_tok_geo_grid_inp *grid_inp);

/**
 * Get R(psi,Z) for a specified psi and Z value. Multiple values may
 * be returned (or none). The R(psi,Z) and dR/dZ are stored in the R
 * and dR arrays which be allocated by the caller.
 *
 * @param geo Geometry object
 * @param psi Psi value
 * @param Z Z value
 * @param nmaxroots Maximum number of roots
 * @param R on output, R(psi,Z)
 * @param dR on output, dR/dZ
 */
int gkyl_tok_geo_R_psiZ(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dRdZ, double *dR, double *dZ);

/**
 * Integrate along a specified psi countour and return its length. The
 * contour must lie completely inside the RZ domain of the psiRZ DG
 * field. The @a rclose parameter is used to select amongst the
 * multiple possible countours with the same psi. Foe example, to
 * select a flux surface on the outboard side of a double-null
 * configuration choose rclose to be Rmax.
 *
 * @param geo Geometry object
 * @param psi Psi value of contour
 * @param zmin Starting z location
 * @param zmax Ending z location
 * @param rclose Value of radial coordinate to discrimate between multiple
 *    contours
 * @return Length of contour
 */
double gkyl_tok_geo_integrate_psi_contour(const struct gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose);

/**
 * Compute physical coordinates (mapc2p)  given computational coordinates
 *
 * @param geo Geometry object
 * @param xn computational coordinates
 * @param ret physical coordinates
 */
void gkyl_tok_geo_mapc2p(const struct gkyl_tok_geo *geo, const struct gkyl_tok_geo_grid_inp *inp,
    const double *xn, double *ret);

/**
 * Compute geometry (mapc2p) on a specified computational grid.
 *
 * @param up gk_geometry object
 * @param nodal range of computational grid
 * @param dzc grid spacing of nodal range
 * @param geo gkyl_tok_geo object with efit dats and root finder specs 
 * @param inp tok_geo_grid_inp Input structure for creating mapc2p
 * @param position_map position map object
 */
void gkyl_tok_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, 
  struct gkyl_tok_geo* geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map);

/**
 * Compute geometry (mapc2p) on a specified computational grid.
 *
 * @param up gk_geometry object
 * @param nodal range of computational grid
 * @param dzc grid spacing of nodal range
 * @param geo gkyl_tok_geo object with efit dats and root finder specs 
 * @param inp tok_geo_grid_inp Input structure for creating mapc2p
 * @param position_map position map object
 */
void gkyl_tok_geo_calc_interior(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  struct gkyl_tok_geo* geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map);

/**
 * Compute geometry (mapc2p) on a specified computational grid.
 *
 * @param up gk_geometry object
 * @param nodal range of computational grid
 * @param dzc grid spacing of nodal range
 * @param geo gkyl_tok_geo object with efit dats and root finder specs 
 * @param inp tok_geo_grid_inp Input structure for creating mapc2p
 * @param position_map position map object
 */
void gkyl_tok_geo_calc_surface(struct gk_geometry* up, int dir, struct gkyl_range *nrange, double dzc[3], 
  struct gkyl_tok_geo* geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map);


/*
 * Get grid extents for a block type based on a global normalization factor
 * and a cut at the arc length of the X-point on the separatrix
 * @param inp grid input
 * @param geo tokamak geometry object
 * @param theta_lo on output the lower grid extent
 * @param theta_up on output the upper grid extent
 * */
void
gkyl_tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up);

/**
 * Return cumulative statistics from geometry computations
 *
 * @param geo Geometry object
 * @return Cumulative statistics
 */
struct gkyl_tok_geo_stat gkyl_tok_geo_get_stat(const struct gkyl_tok_geo *geo);

/**
 * Does this block take the extended, topology-aware construction?
 *
 * Exposed because it is a property of a block that its NEIGHBOURS need to
 * agree with.  The extended path reparameterizes a block's separatrix row (the
 * theta ladder, and the |grad psi| poloidal measure) while the legacy path does
 * not, so two blocks that share that row and disagree about taking it will
 * place different nodes along it -- a seam error that is invisible per block
 * and grows with resolution.  A multiblock consistency check needs to ask this
 * question of both sides of an interface, so the predicate lives here rather
 * than being restated where it is checked.
 *
 * @param inp Grid input for one block
 * @return true if this block uses the extended construction
 */
bool gkyl_tok_geo_uses_extended_construction(const struct gkyl_tok_geo_grid_inp *inp);

/**
 * Delete updater.
 *
 * @param geo Geometry object to delete
 */
void gkyl_tok_geo_release(struct gkyl_tok_geo *geo);

/** Check active material endpoints at a requested flux, without writing a grid.
 * Uses the same native psi evaluator and connected limiter continuation as mapping.
 * False means that this requested surface cannot terminate on its active plate.
 */
bool gkyl_tok_geo_check_plate_coverage(const struct gkyl_tok_geo *geo,
  const struct gkyl_tok_geo_grid_inp *inp, double psi);
