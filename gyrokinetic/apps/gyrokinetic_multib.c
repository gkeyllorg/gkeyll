#include <gkyl_gk_block_geom.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>
#include <gkyl_null_comm.h>
#include <gkyl_rho_wall_adjust_priv.h>
#include <gkyl_tok_geo_wall_trial_priv.h>
#include <gkyl_tok_geo_priv.h>

#include <mpack.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void gyrokinetic_multib_optimize_xpt_seams(
  const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *mbapp);

// Compute total number of ranges specified by cuts.
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// Simple linear search to check if val occurs in lst.
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

// Compute total and maximum number of cuts.
static void
calc_tot_and_max_cuts(const struct gkyl_gk_block_geom *gk_block_geom, int tot_max[2])
{
  int ndim = gkyl_gk_block_geom_ndim(gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(gk_block_geom);

  int max_cuts = 0, tot_cuts = 0;
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(gk_block_geom, i);
    int ncuts = calc_cuts(ndim, bgi->cuts);
    max_cuts = ncuts > max_cuts ? ncuts : max_cuts;
    tot_cuts += ncuts;
  }
  tot_max[0] = tot_cuts;
  tot_max[1] = max_cuts;
}

// Construct single-block App geometry for given block ID.
static struct gkyl_gyrokinetic_app *
singleb_app_new_geom_from_block(const struct gkyl_gyrokinetic_multib *mbinp,
  int bid, const struct gkyl_gyrokinetic_multib_app *mbapp,
  const struct gkyl_gk_block_geom_info *bgi, bool write_geometry)
{
  const char *tok_geo_trace = getenv("GKYL_TOK_GEO_TRACE");
  bool trace_tok_geo = write_geometry && tok_geo_trace &&
    tok_geo_trace[0] != '\0' && tok_geo_trace[0] != '0';
  if (trace_tok_geo) {
    char block_id[32];
    snprintf(block_id, sizeof block_id, "%d", bid);
    setenv("GKYL_TOK_GEO_TRACE_BLOCK", block_id, 1);
    fprintf(stderr, "GKYL_TOK_GEO_TRACE block_start bid=%d\n", bid);
    fflush(stderr);
  }

  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct top-level single-block input on the heap; this struct is large
  // enough to overflow the default stack in geometry-only multiblock tests.
  struct gkyl_gk *app_inp = gkyl_calloc(1, sizeof(*app_inp));

  strcpy(app_inp->name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp->name, app_name.str);
    cstr_drop(&app_name);
  }

  // Set the configuration-space extents, cells, and geometry.
  app_inp->cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp->lower[i] = bgi->lower[i];
    app_inp->upper[i] = bgi->upper[i];
    app_inp->cells[i] = bgi->cells[i];
  }

  // Set z dir grid extents based on tokamak global normalization
  if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK || bgi->geometry.geometry_id == GKYL_GEOMETRY_FROMFILE) {
    gkyl_gk_geometry_tok_set_grid_extents(bgi->geometry.efit_info, bgi->geometry.tok_grid_info, &app_inp->lower[cdim-1], &app_inp->upper[cdim-1]);
    gkyl_gk_block_geom_reset_block_extents(mbapp->gk_block_geom, bid, app_inp->lower, app_inp->upper);
  }

  app_inp->geometry = bgi->geometry;
  // Measure-only prototype: obtain the radial partner from the declaration,
  // keeping the legacy block's off-separatrix cuts and radial domain intact.
  const char *shared_theta = getenv("GKYL_TOK_SHARED_SEP_THETA");
  if (shared_theta && shared_theta[0] && shared_theta[0] != '0' &&
      bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK &&
      !gkyl_tok_geo_uses_extended_construction(&bgi->geometry.tok_grid_info)) {
    for (int e=0; e<2; ++e) {
      const struct gkyl_target_edge *edge = &bgi->connections[0][e];
      if (edge->edge == GKYL_PHYSICAL || edge->bid == bid)
        continue;
      const struct gkyl_gk_block_geom_info *peer =
        gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, edge->bid);
      if (peer->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK ||
          !gkyl_tok_geo_uses_extended_construction(&peer->geometry.tok_grid_info))
        continue;
      const struct gkyl_efit_inp *a = &bgi->geometry.efit_info;
      const struct gkyl_efit_inp *b = &peer->geometry.efit_info;
      if (app_inp->geometry.tok_grid_info.shared_theta_peer || edge->dir != 0 ||
          strcmp(a->filepath, b->filepath) || a->reflect != b->reflect ||
          a->rz_poly_order != b->rz_poly_order || a->flux_poly_order != b->flux_poly_order ||
          a->xpt_bound_n != 0 || b->xpt_bound_n != 0 ||
          bgi->geometry.tok_grid_info.use_cubics != peer->geometry.tok_grid_info.use_cubics) {
        fprintf(stderr, "TOK_SHARED_THETA unsupported interface descriptor block=%d peer=%d\n", bid, edge->bid);
        abort();
      }
      app_inp->geometry.tok_grid_info.shared_theta_peer = &peer->geometry.tok_grid_info;
      app_inp->geometry.tok_grid_info.shared_theta_radial_edge = e;
      app_inp->geometry.tok_grid_info.shared_theta_reverse =
        edge->edge == GKYL_LOWER_NEGATIVE || edge->edge == GKYL_UPPER_NEGATIVE;
      fprintf(stderr, "TOK_SHARED_THETA declared block=%d ftype=%d edge=%d peer=%d peer_ftype=%d\n",
        bid, bgi->geometry.tok_grid_info.ftype, e, edge->bid, peer->geometry.tok_grid_info.ftype);
    }
  }
  // This constructor only builds geometry; species and neutral inputs are not
  // populated in the single-block geometry app.
  app_inp->num_species = 0;
  app_inp->num_neut_species = 0;

  app_inp->poly_order = mbinp->poly_order;
  app_inp->basis_type = mbinp->basis_type;
  app_inp->cfl_frac = mbinp->cfl_frac;
  app_inp->cfl_frac_omegaH = mbinp->cfl_frac_omegaH;

  struct gkyl_comm *comm = mbapp->block_comms[bid];

  struct gkyl_app_parallelism_inp parallel_inp = {};
  parallel_inp.use_gpu = mbinp->use_gpu;
  for (int d=0; d<cdim; ++d) parallel_inp.cuts[d] = bgi->cuts[d];
  parallel_inp.comm = comm;
  // Copy parallelism input into app input.
  memcpy(&app_inp->parallelism, &parallel_inp, sizeof(struct gkyl_app_parallelism_inp));

  app_inp->num_periodic_dir = mbinp->num_periodic_dir;
  for(int i = 0; i < mbinp->cdim; i++)
    app_inp->periodic_dirs[i] = mbinp->periodic_dirs[i];

  app_inp->metadata.num_attributes = mbapp->io_meta_basic_len;
  app_inp->metadata.attributes = mbapp->io_meta_basic;

  struct gkyl_gyrokinetic_app *app = write_geometry
    ? gkyl_gyrokinetic_app_new_geom(app_inp)
    : gkyl_gyrokinetic_app_new_geom_no_write(app_inp);
  gkyl_free(app_inp);
  if (trace_tok_geo) {
    fprintf(stderr, "GKYL_TOK_GEO_TRACE block_done bid=%d\n", bid);
    fflush(stderr);
  }
  return app;
}

static struct gkyl_gyrokinetic_app *
singleb_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp)
{
  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);
  return singleb_app_new_geom_from_block(mbinp, bid, mbapp, bgi, false);
}

// Construct single-block App solver for given block ID.
static void
singleb_app_new_solver(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp, struct gkyl_gyrokinetic_app *app)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);

  // Construct top-level single-block input on the heap; this struct is large
  // enough to overflow the default stack in multiblock solver construction.
  struct gkyl_gk *app_inp = gkyl_calloc(1, sizeof(*app_inp));

  // Set the configuration-space extents, cells.
  app_inp->cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp->lower[i] = bgi->lower[i];
    app_inp->upper[i] = bgi->upper[i];
    app_inp->cells[i] = bgi->cells[i];
  }

  int num_species = app_inp->num_species = mbinp->num_species;
  int num_neut_species = app_inp->num_neut_species = mbinp->num_neut_species; 

  app_inp->poly_order = mbinp->poly_order;
  app_inp->basis_type = mbinp->basis_type;
  app_inp->cfl_frac = mbinp->cfl_frac; 
  app_inp->cfl_frac_omegaH = mbinp->cfl_frac_omegaH; 
  app_inp->eirene = mbinp->eirene;

  for (int i=0; i<num_species; ++i) {
    const struct gkyl_gyrokinetic_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_gyrokinetic_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;

    // Velocity-space information
    int vdim = species_inp.vdim = sp->vdim;
    for (int v=0; v<vdim; ++v) {
      species_inp.lower[v] = sp->lower[v];
      species_inp.upper[v] = sp->upper[v];
      species_inp.cells[v] = sp->cells[v];      
    }
    species_inp.mapc2p = sp->mapc2p;

    // Species physics modules.
    species_inp.collisionless = sp->collisionless;
    species_inp.collisions = sp->collisions;
    species_inp.anomalous_diffusion = sp->anomalous_diffusion;
    species_inp.radiation = sp->radiation;
    species_inp.react = sp->react;
    species_inp.react_neut = sp->react_neut; 

    // Species diagnostics
    species_inp.num_diag_moments = sp->num_diag_moments;
    for (int n=0; n<species_inp.num_diag_moments; ++n) {
      species_inp.diag_moments[n] = sp->diag_moments[n];
    }
    species_inp.num_integrated_diag_moments = sp->num_integrated_diag_moments;
    for (int n=0; n<species_inp.num_integrated_diag_moments; ++n) {
      species_inp.integrated_diag_moments[n] = sp->integrated_diag_moments[n];
    }
    species_inp.time_rate_diagnostics = sp->time_rate_diagnostics;
    species_inp.boundary_flux_diagnostics = sp->boundary_flux_diagnostics;

    // Choose proper block-specific species input.
    const struct gkyl_gyrokinetic_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
      }
    }
    species_inp.projection = sp_pb->projection;
    species_inp.source = sp_pb->source;
    species_inp.polarization_density = sp_pb->polarization_density;  

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    int pardir = cdim-1;
    int num_below = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_BELOW);
    int num_above = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_ABOVE);

    int bc_count_sp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_sp[i] = 0;

    if (cdim ==3) {
      if (num_below > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_LOWER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }

      if (num_above > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_UPPER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }
    }

    // Set species physical BCs.
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {
        species_inp.bcs[bc_count_sp[bid]].dir = sp->bcs[i].dir;
        species_inp.bcs[bc_count_sp[bid]].edge = sp->bcs[i].edge;
        species_inp.bcs[bc_count_sp[bid]].type = sp->bcs[i].type;
        species_inp.bcs[bc_count_sp[bid]].aux_profile = sp->bcs[i].aux_profile;
        species_inp.bcs[bc_count_sp[bid]].aux_ctx = sp->bcs[i].aux_ctx;
        species_inp.bcs[bc_count_sp[bid]].projection = sp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          species_inp.bcs[bc_count_sp[bid]].value[k] = sp->bcs[i].value[k];

        bc_count_sp[bid] += 1;
      }
    }

    // Copy species input into app input.
    memcpy(&app_inp->species[i], &species_inp, sizeof(struct gkyl_gyrokinetic_species));
  }

  for (int i=0; i<num_neut_species; ++i) {
    const struct gkyl_gyrokinetic_multib_neut_species *nsp = &mbinp->neut_species[i];
    
    struct gkyl_gyrokinetic_neut_species neut_species_inp = { };
    strcpy(neut_species_inp.name, nsp->name);

    neut_species_inp.mass = nsp->mass; 
    neut_species_inp.is_static = nsp->is_static; 

    // Velocity space information (neutrals are 3V)
    int vdim = neut_species_inp.vdim = nsp->vdim;
    for (int v=0; v<vdim; ++v) {
      neut_species_inp.lower[v] = nsp->lower[v];
      neut_species_inp.upper[v] = nsp->upper[v];
      neut_species_inp.cells[v] = nsp->cells[v];
    }
    neut_species_inp.mapc2p = nsp->mapc2p;

    // Neutral species physics modules
    neut_species_inp.react_neut = nsp->react_neut;

    // Neutral species diagnostics
    neut_species_inp.num_diag_moments = nsp->num_diag_moments;
    for (int n=0; n<neut_species_inp.num_diag_moments; ++n) {
      neut_species_inp.diag_moments[n] = nsp->diag_moments[n];
    }

    // choose proper block-specific species input
    const struct gkyl_gyrokinetic_multib_neut_species_pb *nsp_pb = &nsp->blocks[0];
    if (!nsp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == nsp->blocks[i].block_id) {
          nsp_pb = &nsp->blocks[i];
          break;
        }
      }
    }
    neut_species_inp.projection = nsp_pb->projection;
    neut_species_inp.source = nsp_pb->source;

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    // Set species physical BCs.
    int bc_count_nsp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_nsp[i] = 0;

    for (int i=0; i<nsp->num_physical_bcs; ++i) {
      if (bid == nsp->bcs[i].bidx) {
        neut_species_inp.bcs[bc_count_nsp[bid]].dir = nsp->bcs[i].dir;
        neut_species_inp.bcs[bc_count_nsp[bid]].edge = nsp->bcs[i].edge;
        neut_species_inp.bcs[bc_count_nsp[bid]].type = nsp->bcs[i].type;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_profile = nsp->bcs[i].aux_profile;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_ctx = nsp->bcs[i].aux_ctx;
        neut_species_inp.bcs[bc_count_nsp[bid]].projection = nsp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          neut_species_inp.bcs[bc_count_nsp[bid]].value[k] = nsp->bcs[i].value[k];

        bc_count_nsp[bid] += 1;
      }
    }

    // Copy neutral species input into app input.
    memcpy(&app_inp->neut_species[i], &neut_species_inp, sizeof(struct gkyl_gyrokinetic_neut_species));
  } 

  // Initialize the single-block field solver (only used for num_blocks=1).
  const struct gkyl_gyrokinetic_multib_field *fld = &mbinp->field;
  struct gkyl_gyrokinetic_field field_inp = { };
  field_inp.gkfield_id = fld->gkfield_id;
  field_inp.kperpSq = fld->kperpSq; 
  field_inp.time_rate_diagnostics = fld->time_rate_diagnostics; 
  
  // Adiabatic electron inputs.
  field_inp.electron_mass = fld->electron_mass;
  field_inp.electron_charge = fld->electron_charge;
  field_inp.electron_density = fld->electron_density; 
  field_inp.electron_temp = fld->electron_temp; 
  
  // BCs.
  for (int d=0; d<cdim; d++) {
    for (int i=0; i<2; i++) {
      // Set it to Dirichlet first, reset below. This avoids problems in
      // creating the single block field solve, which may not be used.
      field_inp.poisson_bcs[2*d+i].dir = d;
      field_inp.poisson_bcs[2*d+i].edge = i==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
      field_inp.poisson_bcs[2*d+i].type = GKYL_BC_GK_FIELD_DIRICHLET;
      for (int k=0; k<3; ++k)
        field_inp.poisson_bcs[i].value[k] = -1.0e3;
    }
  }
  
  for (int i=0; i<fld->num_physical_bcs; i++) { 
    if (bid == fld->bcs[i].bidx) {
      struct gkyl_gyrokinetic_bc *bc_curr = gk_fetch_bc_with_dir_edge(field_inp.poisson_bcs, 2*cdim, fld->bcs[i].dir, fld->bcs[i].edge);
      bc_curr->type = fld->bcs[i].type;
      bc_curr->aux_profile = fld->bcs[i].aux_profile;
      bc_curr->aux_ctx = fld->bcs[i].aux_ctx;
      for (int k=0; k<3; ++k)
        bc_curr->value[k] = fld->bcs[i].value[k];
    }
  }
  
  const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[0];
  // Choose proper block-specific field input.
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  field_inp.polarization_bmag = fld_pb->polarization_bmag ? fld_pb->polarization_bmag : mbapp->bmag_ref;
  field_inp.kperpSq = fld_pb->kperpSq;
  field_inp.time_rate_diagnostics = fld_pb->time_rate_diagnostics; 
  
  field_inp.phi_wall_lo_ctx = fld_pb->phi_wall_lo_ctx; 
  field_inp.phi_wall_lo = fld_pb->phi_wall_lo; 
  field_inp.phi_wall_lo_evolve = fld_pb->phi_wall_lo_evolve; 
  
  field_inp.phi_wall_up_ctx = fld_pb->phi_wall_up_ctx; 
  field_inp.phi_wall_up = fld_pb->phi_wall_up; 
  field_inp.phi_wall_up_evolve = fld_pb->phi_wall_up_evolve;   
  

  // Copy field input into app input.
  memcpy(&app_inp->field, &field_inp, sizeof(struct gkyl_gyrokinetic_field));  

  gkyl_gyrokinetic_app_new_solver(app_inp, app);
  gkyl_free(app_inp);
}

// The grid's psi values are the MAPPED ones: tok_geo applies maps[0] to psi
// before using it. Sampling the raw computational bounds therefore inspects
// surfaces the grid does not have, and goes wrong in both directions -- passing
// a domain that cannot be built, and refusing one that can while naming a psi
// the grid never uses.
//
// Only some maps can be evaluated this early. USER_INPUT carries the caller's
// own function and is usable immediately; the constant-dB and X-point maps are
// installed by gkyl_position_map_optimize, which needs the geometry this check
// runs before. Returns NULL when the map cannot be resolved yet, so callers can
// say the check is advisory rather than assert a verdict they cannot support.
static mc2nu_t
tok_preflight_radial_map(const struct gkyl_gk_block_geom_info *bi, void **ctx)
{
  const struct gkyl_position_map_inp *pmi = &bi->geometry.position_map_info;
  *ctx = pmi->ctxs[0];
  if (!pmi->maps[0]) return 0;   // unset: identity, and raw == mapped
  if (pmi->id == GKYL_PMAP_USER_INPUT ||
      pmi->id == GKYL_PMAP_USER_INPUT_W_DERIVATIVE)
    return pmi->maps[0];
  return 0;
}

// True when a map exists but cannot be resolved at preflight time.
static bool
tok_preflight_map_unresolved(const struct gkyl_gk_block_geom_info *bi)
{
  const struct gkyl_position_map_inp *pmi = &bi->geometry.position_map_info;
  return pmi->id != GKYL_PMAP_USER_INPUT &&
         pmi->id != GKYL_PMAP_USER_INPUT_W_DERIVATIVE;
}

// Cheap material-domain rejection before allocating/writing any block. The
// actual corner/interior/face passes also enforce the material boundary.
static bool
gyrokinetic_multib_material_preflight(const struct gkyl_gyrokinetic_multib *inp)
{
  const struct gkyl_gk_block_geom *bg=inp->gk_block_geom;
  for (int b=0; b<gkyl_gk_block_geom_num_blocks(bg); ++b) {
    const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
    if (bi->geometry.geometry_id!=GKYL_GEOMETRY_TOKAMAK) continue;
    const struct gkyl_tok_geo_grid_inp *ti=&bi->geometry.tok_grid_info;
    struct gkyl_tok_geo *geo=gkyl_tok_geo_new(&bi->geometry.efit_info,ti);
    enum gkyl_tok_wall_policy policy=gkyl_tok_wall_policy_for(ti,geo->efit);
    bool ok=policy!=GKYL_TOK_WALL_REJECT_UNDECLARED &&
            policy!=GKYL_TOK_WALL_REJECT_CONTRADICTED &&
            policy!=GKYL_TOK_WALL_REJECT_MALFORMED;
    if (!ok)
      fprintf(stderr,"TOK_GEO_WALL_UNAVAILABLE block=%d limiter_status=%d vertices=%d reason=%s\n",
        b,geo->efit->limiter_status,geo->efit->limiter_n,gkyl_tok_wall_policy_reason(policy));
    // An acknowledged absent outline builds without wall enforcement. Say so on
    // every such block, so an unenforced wall is never silent in the log.
    if (policy==GKYL_TOK_WALL_NOT_ENFORCED)
      fprintf(stderr,"TOK_GEO_WALL_NOT_ENFORCED block=%d limiter_status=%d vertices=%d\n",
        b,geo->efit->limiter_status,geo->efit->limiter_n);
    // Includes both radial bounds. This early check does not substitute for
    // checks at the actual mapped radial coordinates during construction.
    // Plate coverage is a material-target question, so it is only meaningful
    // where a wall is enforced.
    void *pmctx=0;
    mc2nu_t pmap=tok_preflight_radial_map(bi,&pmctx);
    bool unresolved=tok_preflight_map_unresolved(bi);
    if (unresolved)
      fprintf(stderr,"GKYL_MATERIAL_PREFLIGHT scope=advisory block=%d "
        "reason=position_map_not_resolved_yet id=%d\n",
        b,(int)bi->geometry.position_map_info.id);
    // An advisory check cannot reject based on raw surfaces the resolved map
    // may never use. Actual mapped geometry still enforces plate/wall limits.
    for (int i=0; i<=32 && ok && !unresolved && policy==GKYL_TOK_WALL_ENFORCE; ++i) {
      double psi=bi->lower[0]+(bi->upper[0]-bi->lower[0])*i/32.0;
      if (pmap) pmap(0.0,&psi,&psi,pmctx);
      ok=gkyl_tok_geo_check_plate_coverage(geo,ti,psi);
    }
    gkyl_tok_geo_release(geo);
    if (!ok) {
      fprintf(stderr,"GKYL_MATERIAL_PREFLIGHT status=FAIL block=%d\n",b);
      return false;
    }
  }
  fprintf(stderr,"GKYL_MATERIAL_PREFLIGHT status=PASS\n");
  return true;
}

struct rho_wall_bound {
  int family, edge, group, steps;
  double requested, requested_psi, other, axis, sep;
};

static bool
rho_wall_isfinite(double x)
{
  uint64_t bits;
  memcpy(&bits, &x, sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

static bool
rho_wall_same(double a, double b)
{
  return rho_wall_isfinite(a) && rho_wall_isfinite(b) &&
    fabs(a-b) <= 128.0*DBL_EPSILON*fmax(1.0, fmax(fabs(a),fabs(b)));
}

static int
rho_wall_group(struct rho_wall_bound *bounds, int b)
{
  while (bounds[b].group != b) b = bounds[b].group;
  return b;
}

// Return an owned effective declaration. Every trial uses the actual mapper,
// fresh compression/map state and the requested cell counts. Trials have no
// output and are always released. Only the subsequent ordinary constructor
// can expose a geometry; its wall checks are never in collection mode.
static struct gkyl_gk_block_geom *
gyrokinetic_multib_adjust_wall(const struct gkyl_gyrokinetic_multib *inp)
{
  const int n = gkyl_gk_block_geom_num_blocks(inp->gk_block_geom);
  const int ndim = gkyl_gk_block_geom_ndim(inp->gk_block_geom);
  // Every allocation below is sized by n, and n reaches memset() as
  // n*sizeof(*p). Nothing here bounded n, so the compiler had to assume a
  // negative value was reachable and warned that the length could wrap
  // (-Wstringop-overflow at the two memsets in the iteration loop). State the
  // invariant once, before anything is allocated, rather than silencing it.
  if (n <= 0 || ndim <= 0) {
    fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=empty_block_geometry blocks=%d ndim=%d\n",n,ndim);
    return 0;
  }
  struct gkyl_gk_block_geom *bg = gkyl_gk_block_geom_new(ndim,n);
  struct rho_wall_bound *bounds = gkyl_calloc(n,sizeof(*bounds));
  bool *failed = gkyl_calloc(n,sizeof(*failed));
  bool *bad_block = gkyl_calloc(n,sizeof(*bad_block));
  int *order = gkyl_malloc(n*sizeof(*order));
  for (int b=0; b<n; ++b) order[b]=b;
  struct gkyl_gyrokinetic_multib *probe_inp = gkyl_malloc(sizeof(*probe_inp));
  *probe_inp = *inp;
  probe_inp->gk_block_geom = bg;
  probe_inp->use_gpu = false;
  struct gkyl_gyrokinetic_multib_app *probe = gkyl_calloc(1,sizeof(*probe));
  probe->gk_block_geom = bg;
  probe->block_comms = gkyl_calloc(n,sizeof(*probe->block_comms));
  struct gkyl_comm *serial = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp) {
    .use_gpu = false,
  });
  bool ok = false;

  // Initial opt-in scope: p1, serial application, identity or fresh built-in
  // X-point compression. Other modes need their own trial/transfer contract.
  int ranks = 0;
  gkyl_comm_get_size(inp->comm,&ranks);
  if (ranks != 1 || inp->poly_order != 1) {
    fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=unsupported_configuration ranks=%d poly_order=%d\n",ranks,inp->poly_order);
    goto cleanup;
  }

  for (int b=0; b<n; ++b) {
    const struct gkyl_gk_block_geom_info *bi = gkyl_gk_block_geom_get_block(inp->gk_block_geom,b);
    gkyl_gk_block_geom_set_block(bg,b,bi);
    bounds[b].group = b;
    probe->block_comms[b] = serial;
    if (bi->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK) continue;
    if (bi->geometry.tok_grid_info.relaxed_xpt_seam_optimize) {
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=unsupported_seam_optimizer block=%d\n",b);
      goto cleanup;
    }
    const struct gkyl_position_map_inp *pm = &bi->geometry.position_map_info;
    if ((pm->id != GKYL_PMAP_USER_INPUT && pm->id != GKYL_PMAP_XPT_COMPRESSION) ||
        pm->maps[0] || pm->maps[1] || pm->maps[2]) {
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=unsupported_position_map block=%d\n",b);
      goto cleanup;
    }
    struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&bi->geometry.efit_info,&bi->geometry.tok_grid_info);
    double axis=geo->efit->simag, sep=geo->psisep;
    // Adjustment moves a bound so the domain fits inside the wall, so it is
    // only meaningful where a wall exists. Same predicate as everywhere else.
    bool wall_ok=gkyl_tok_wall_usable(geo->efit);
    gkyl_tok_geo_release(geo);
    if (!wall_ok || !rho_wall_isfinite(axis) || !rho_wall_isfinite(sep) || axis==sep) {
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=wall_or_flux_unavailable block=%d\n",b);
      goto cleanup;
    }
    bounds[b].axis=axis; bounds[b].sep=sep;
    const int family=gkyl_rho_wall_family(bi->geometry.tok_grid_info.ftype);
    if (!family) continue; // Core and other types are checked but never moved.
    double r2lo=(bi->lower[0]-axis)/(sep-axis), r2hi=(bi->upper[0]-axis)/(sep-axis);
    if (!rho_wall_isfinite(r2lo) || !rho_wall_isfinite(r2hi) || r2lo<0 || r2hi<0) {
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=invalid_rho block=%d\n",b);
      goto cleanup;
    }
    double rlo=sqrt(r2lo), rhi=sqrt(r2hi);
    int edge = family==1 ? (rhi>rlo) : (rhi<rlo);
    // Only move a physical radial boundary; never move a declared radial join.
    if (bi->connections[0][edge].edge != GKYL_PHYSICAL) continue;
    bounds[b].family=family; bounds[b].edge=edge;
    bounds[b].requested=edge ? rhi : rlo;
    bounds[b].requested_psi=edge ? bi->upper[0] : bi->lower[0];
    bounds[b].other=edge ? rlo : rhi;
  }

  // Same-family theta/alpha neighbors own the same outer boundary. Keep
  // unrelated PF/core contours distinct even when their numeric flux agrees.
  for (int b=0; b<n; ++b) {
    const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
    if (!bounds[b].family) continue;
    for (int d=1; d<ndim; ++d) for (int e=0; e<2; ++e) {
      const struct gkyl_target_edge *te=&bi->connections[d][e];
      if (te->edge==GKYL_PHYSICAL || te->bid==b) continue;
      int j=te->bid;
      const struct gkyl_gk_block_geom_info *bj=gkyl_gk_block_geom_get_block(bg,j);
      if (bj->geometry.geometry_id==GKYL_GEOMETRY_TOKAMAK &&
          gkyl_rho_wall_family(bj->geometry.tok_grid_info.ftype)==bounds[b].family && !bounds[j].family) {
        fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=mixed_boundary_ownership block=%d peer=%d\n",b,j);
        goto cleanup;
      }
      if (bounds[b].family!=bounds[j].family) continue;
      const struct gkyl_efit_inp *a=&bi->geometry.efit_info, *c=&bj->geometry.efit_info;
      if (!rho_wall_same(bounds[b].requested,bounds[j].requested) ||
          !rho_wall_same(bounds[b].other,bounds[j].other) ||
          !rho_wall_same(bounds[b].axis,bounds[j].axis) || !rho_wall_same(bounds[b].sep,bounds[j].sep) ||
          strcmp(a->filepath,c->filepath) || a->reflect!=c->reflect ||
          a->rz_poly_order!=c->rz_poly_order || a->flux_poly_order!=c->flux_poly_order ||
          bi->geometry.tok_grid_info.use_cubics!=bj->geometry.tok_grid_info.use_cubics) {
        fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=incompatible_boundary_group block=%d peer=%d\n",b,j);
        goto cleanup;
      }
      int rb=rho_wall_group(bounds,b), rj=rho_wall_group(bounds,j);
      if (rb!=rj) bounds[rj].group=rb;
    }
  }
  for (int b=0; b<n; ++b) bounds[b].group=rho_wall_group(bounds,b);

  for (int iteration=0; iteration<=10000; ++iteration) {
    memset(failed,0,n*sizeof(*failed));
    memset(bad_block,0,n*sizeof(*bad_block));
    bool any=false;
    // An explicitly supplied wall target can lose its far-bound root because
    // the requested contour exits elsewhere on the wall (ASDEX wide PF).
    // Permit the same inward rho search before constructing that contour.
    // Original plate failures and fixed-bound target failures still reject.
    for (int b=0;b<n;++b) {
      const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
      if (bi->geometry.geometry_id!=GKYL_GEOMETRY_TOKAMAK) continue;
      const struct gkyl_tok_geo_grid_inp *ti=&bi->geometry.tok_grid_info;
      struct gkyl_tok_geo *geo=gkyl_tok_geo_new(&bi->geometry.efit_info,ti);
      double fixed=bounds[b].edge ? bi->lower[0] : bi->upper[0];
      if (tok_plate_coverage_status(geo,ti,fixed)) {
        gkyl_tok_geo_release(geo);
        fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=fixed_boundary_target_unavailable block=%d\n",b);
        goto cleanup;
      }
      bool fatal=false;
      void *pmctx=0;
      mc2nu_t pmap=tok_preflight_radial_map(bi,&pmctx);
      for (int j=0;j<=32;++j) {
        double psi=bi->lower[0]+(bi->upper[0]-bi->lower[0])*j/32.0;
        if (pmap) pmap(0.0,&psi,&psi,pmctx);
        int status=tok_plate_coverage_status(geo,ti,psi);
        if (!status) continue;
        if (status!=1 || !bounds[b].family) { fatal=true;break; }
        failed[bounds[b].group]=true;bad_block[b]=true;any=true;
        fprintf(stderr,"TOK_RHO_WALL_TARGET_RETRY block=%d family=%s psi=%.17g reason=explicit_target_root_unavailable\n",
          b,bounds[b].family==1 ? "SOL" : "PF",psi);
        // Other active endpoints/fluxes are checked on the next candidate;
        // no geometry is constructed while any target lacks coverage.
        break;
      }
      gkyl_tok_geo_release(geo);
      if (fatal) {
        fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=original_plate_coverage block=%d\n",b);
        goto cleanup;
      }
    }
    if (any) goto advance_wall_bounds;
    for (int k=0; k<n; ++k) {
      int b=order[k];
      // Once a group fails, this candidate cannot be selected. Probe its
      // previous offender first next time, and defer the other members and
      // unchanged core until a candidate can actually pass. Iteration zero
      // visits the original ordering and captures every requested contour.
      if (iteration>0 && (failed[bounds[b].group] || (any && !bounds[b].family))) continue;
      struct gkyl_gk_block_geom_info bi=*gkyl_gk_block_geom_get_block(bg,b);
      if (bi.geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK) continue;
      for (int d=0; d<ndim; ++d) bi.cuts[d]=1;
      fprintf(stderr,"TOK_RHO_WALL_TRIAL iteration=%d block=%d\n",iteration,b);
      tok_wall_trial_begin(bounds[b].family ? bounds[b].edge : -1);
      if (iteration==0 && bounds[b].family)
        tok_wall_trial_capture_requested(b,bounds[b].requested);
      struct gkyl_gyrokinetic_app *app=singleb_app_new_geom_from_block(probe_inp,b,probe,&bi,false);
      long violations=tok_wall_trial_end();
      bool fixed_violation=tok_wall_trial_has_fixed_violation();
      gkyl_gyrokinetic_app_release_geom(app);
      if (violations) {
        fprintf(stderr,"TOK_RHO_WALL_TRIAL_REJECTED iteration=%d block=%d wall_checks_failed=%ld\n",iteration,b,violations);
        if (fixed_violation) {
          fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=fixed_boundary_wall_violation block=%d\n",b);
          goto cleanup;
        }
        if (!bounds[b].family) {
          fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=unadjustable_wall_violation block=%d\n",b);
          goto cleanup;
        }
        failed[bounds[b].group]=true;
        bad_block[b]=true;
        any=true;
      }
    }
    if (!any) {
      for (int b=0; b<n; ++b) if (bounds[b].family) {
        const struct gkyl_gk_block_geom_info *bi=gkyl_gk_block_geom_get_block(bg,b);
        double psi=bounds[b].edge ? bi->upper[0] : bi->lower[0];
        double rho=sqrt((psi-bounds[b].axis)/(bounds[b].sep-bounds[b].axis));
        fprintf(stderr,"TOK_RHO_WALL_BOUND block=%d family=%s edge=%d requested_rho=%.17g effective_rho=%.17g requested_psi=%.17g effective_psi=%.17g steps=%d\n",
          b,bounds[b].family==1 ? "SOL" : "PF",bounds[b].edge,bounds[b].requested,rho,bounds[b].requested_psi,psi,bounds[bounds[b].group].steps);
      }
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_SELECTED iterations=%d step_rho=%g; rebuilding with hard wall guards\n",iteration,GKYL_RHO_WALL_STEP);
      ok=true;
      break;
    }
advance_wall_bounds:
    if (iteration==10000) break;
    int next=0;
    for (int b=0; b<n; ++b) if (bad_block[b]) order[next++]=b;
    for (int b=0; b<n; ++b) if (!bad_block[b] && bounds[b].family) order[next++]=b;
    for (int b=0; b<n; ++b) if (!bounds[b].family) order[next++]=b;
    for (int g=0; g<n; ++g) if (failed[g]) ++bounds[g].steps;
    for (int b=0; b<n; ++b) {
      int g=bounds[b].group;
      if (!failed[g]) continue;
      double rho,psi;
      if (!gkyl_rho_wall_next(bounds[g].requested,bounds[g].other,bounds[g].axis,bounds[g].sep,
          bounds[g].family,bounds[g].steps,&rho,&psi)) {
        fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED reason=no_admissible_increment group=%d requested_rho=%.17g steps=%d\n",
          g,bounds[g].requested,bounds[g].steps);
        goto cleanup;
      }
      struct gkyl_gk_block_geom_info bi=*gkyl_gk_block_geom_get_block(bg,b);
      if (bounds[b].edge) bi.upper[0]=psi; else bi.lower[0]=psi;
      gkyl_gk_block_geom_set_block(bg,b,&bi);
      fprintf(stderr,"TOK_RHO_WALL_ADJUST_STEP block=%d family=%s step=%d rho=%.17g psi=%.17g\n",
        b,bounds[b].family==1 ? "SOL" : "PF",bounds[g].steps,rho,psi);
    }
  }

cleanup:
  gkyl_comm_release(serial);
  gkyl_free(probe->block_comms); gkyl_free(probe); gkyl_free(probe_inp);
  gkyl_free(bounds); gkyl_free(failed);
  gkyl_free(bad_block); gkyl_free(order);
  if (!ok) { gkyl_gk_block_geom_release(bg); return 0; }
  return bg;
}

static gkyl_gyrokinetic_multib_app *gyrokinetic_multib_app_new_geom_impl(const struct gkyl_gyrokinetic_multib *);
static gkyl_gyrokinetic_multib_app *gyrokinetic_multib_app_new_impl(const struct gkyl_gyrokinetic_multib *);

static gkyl_gyrokinetic_multib_app *
gyrokinetic_multib_app_wall_wrapper(const struct gkyl_gyrokinetic_multib *inp, bool geometry_only)
{
  const char *adjust=getenv("ADJUST_IF_EXCEEDING_WALL");
  if (adjust && strcmp(adjust,"0") && strcmp(adjust,"1")) {
    fprintf(stderr,"TOK_RHO_WALL_ADJUST_FAILED ADJUST_IF_EXCEEDING_WALL must be 0 or 1\n");
    return 0;
  }
  if (!adjust || !strcmp(adjust,"0"))
    return geometry_only ? gyrokinetic_multib_app_new_geom_impl(inp) : gyrokinetic_multib_app_new_impl(inp);
  if (!gkyl_gyrokinetic_multib_app_geometry_preflight(inp)) return 0;
  struct gkyl_gk_block_geom *bg=gyrokinetic_multib_adjust_wall(inp);
  if (!bg) return 0;
  struct gkyl_gyrokinetic_multib *effective=gkyl_malloc(sizeof(*effective));
  *effective=*inp; effective->gk_block_geom=bg;
  gkyl_gyrokinetic_multib_app *app=geometry_only ?
    gyrokinetic_multib_app_new_geom_impl(effective) : gyrokinetic_multib_app_new_impl(effective);
  gkyl_free(effective);
  gkyl_gk_block_geom_release(bg);
  return app;
}

gkyl_gyrokinetic_multib_app *
gkyl_gyrokinetic_multib_app_new_geom(const struct gkyl_gyrokinetic_multib *inp)
{
  return gyrokinetic_multib_app_wall_wrapper(inp,true);
}

gkyl_gyrokinetic_multib_app *
gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *inp)
{
  return gyrokinetic_multib_app_wall_wrapper(inp,false);
}

static gkyl_gyrokinetic_multib_app*
gyrokinetic_multib_app_new_geom_impl(const struct gkyl_gyrokinetic_multib *mbinp)
{
  // Reject invalid declarations before communicator access or app allocation.
  if (!gkyl_gyrokinetic_multib_app_geometry_preflight(mbinp))
    return 0;
  if (!gyrokinetic_multib_material_preflight(mbinp))
    return 0;

  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);

  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  // Create multiblock metadata to pass to each block, and combine with user's metadata.
  const char *fmt_btopo = "%s-block_topo.gkyl";
  int sz = gkyl_calc_strlen(fmt_btopo, mbapp->name);
  char fileNm_btopo[sz+1]; // ensures no buffer overflow
  snprintf(fileNm_btopo, sizeof fileNm_btopo, fmt_btopo, mbapp->name);

  // Basic metadata for I/O (including metadata optional from user).
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_default[] = {
    { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
    { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
    { .key = "is_multib", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1 },
    { .key = "topo_file", .elem_type = GKYL_MP_STRING, .cval = fileNm_btopo },
  };
  const struct gkyl_msgpack_map_elem *io_meta_union[] = {io_meta_default, mbinp->metadata.attributes};
  int io_meta_union_len[] = {sizeof(io_meta_default)/sizeof(io_meta_default[0]), mbinp->metadata.num_attributes};

  mbapp->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_union)/sizeof(io_meta_union[0]),
    io_meta_union_len, io_meta_union, &mbapp->io_meta_basic_len);

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).


    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i) {
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);
  }

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  gyrokinetic_multib_optimize_xpt_seams(mbinp, mbapp);

  // NOTE on multiblock interface metric continuity: the raw per-block interior
  // geometry (geo_int.jacobgeo, g_ij, etc.) is never reconciled across a
  // z-direction block interface -- each block's DG geometry is constructed
  // independently, and gkyl_multib_comm_conn_array_transfer below only
  // ghost-copies (gkyl_array_copy_range_to_range in multib_comm_conn_*.c),
  // never averages or blends. Real, sometimes large (tens of percent)
  // discontinuities in the raw interior metric at a shared interface are
  // therefore expected, not a bug -- see
  // ixf-vs-main-comp/IXF_VS_MAIN_CHANGES.md section 6.2 and
  // ixf-vs-main-comp/check_interface_continuity.py for quantified bounds
  // (identical between this branch and main wherever neither side's block
  // uses straight_xpt_ray, confirming this is pre-existing architecture,
  // not something this branch introduced or changed).
  // What IS explicitly reconciled below is the *surface* Jacobian used by
  // flux-relevant kernels (jacobgeo_ratio: each block's own surface
  // jacobgeo is ghost-synced, then divided into the local jacobgeo to give
  // a per-interface correction factor) and the surface arc length
  // (deltats, for tokamak geometry) -- both predate this branch and are the
  // actual mechanism protecting evolved-field flux computations from the
  // raw metric discontinuity above. Their sufficiency (whether the
  // correction fully compensates for the observed discontinuity
  // magnitudes) has not been independently verified.
  //
  // Sync the conf-space volume Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs_vol[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_ghost;
  }
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_vol, jacs_vol);

  // Sync the surface conf-space Jacobian, compute its reciprocal, and store its
  // product with the Jacobian of this block in the ghost cell.
  for (int d = 0; d<cdim; d++) {
    struct gkyl_array *jacs[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      jacs[b] = geo_surf.jacobgeo_ratio;
      gkyl_array_copy_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range_to_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *jacgeo = mkarr(mbapp->use_gpu, geo_surf.jacobgeo_ratio->ncomp, geo_surf.jacobgeo_ratio->size);

      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_lower_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_upper_ghost[d]);
      gkyl_array_copy_range_to_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_lower_ghost[d], &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_upper_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_array_clear_range(geo_surf.jacobgeo_ratio, 0.0, &sbapp->local);
      gkyl_array_shiftc_range(geo_surf.jacobgeo_ratio, pow(sqrt(2.0),sbapp->cdim), 0, &sbapp->local);

      gkyl_array_release(jacgeo);
    }
  }

  const struct gkyl_gk_block_geom_info *bgi0 = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, 0);

  if (cdim > 1 && bgi0->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    for (int d = 0; d<cdim; d++) {
      struct gkyl_array *deltats[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
        deltats[b] = geo_surf.deltats;
        gkyl_array_copy_range_to_range(deltats[b], deltats[b], &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, deltats, deltats);
    }
  }

  // No partial grid: all local blocks and their hard material checks must
  // finish on every rank before any block's geometry is written.
  gkyl_comm_barrier(mbapp->comm);
  gkyl_gyrokinetic_multib_app_write_geometry(mbapp);
  return mbapp;
}

static gkyl_gyrokinetic_multib_app* gyrokinetic_multib_app_new_impl(const struct gkyl_gyrokinetic_multib *mbinp)
{
  // Reject invalid declarations before communicator access or app allocation.
  if (!gkyl_gyrokinetic_multib_app_geometry_preflight(mbinp))
    return 0;
  if (!gyrokinetic_multib_material_preflight(mbinp))
    return 0;

  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);

  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  // Create multiblock metadata to pass to each block, and combine with user's metadata.
  const char *fmt_btopo = "%s-block_topo.gkyl";
  int sz = gkyl_calc_strlen(fmt_btopo, mbapp->name);
  char fileNm_btopo[sz+1]; // ensures no buffer overflow
  snprintf(fileNm_btopo, sizeof fileNm_btopo, fmt_btopo, mbapp->name);

  // Basic metadata for I/O (including metadata optional from user).
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_default[] = {
    { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
    { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
    { .key = "is_multib", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1 },
    { .key = "topo_file", .elem_type = GKYL_MP_STRING, .cval = fileNm_btopo },
  };
  const struct gkyl_msgpack_map_elem *io_meta_union[] = {io_meta_default, mbinp->metadata.attributes};
  int io_meta_union_len[] = {sizeof(io_meta_default)/sizeof(io_meta_default[0]), mbinp->metadata.num_attributes};

  mbapp->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_union)/sizeof(io_meta_union[0]),
    io_meta_union_len, io_meta_union, &mbapp->io_meta_basic_len);

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).

    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);

  // Set bmag_ref.
  double bmag_min_local = DBL_MAX;
  double bmag_min_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_min = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MIN);
    bmag_min_local = GKYL_MIN2(bmag_min_local, bmag_min);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MIN, 1, &bmag_min_local, &bmag_min_global);

  double bmag_max_local = -DBL_MAX;
  double bmag_max_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_max = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MAX);
    bmag_max_local = GKYL_MAX2(bmag_max_local, bmag_max);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);

  mbapp->bmag_ref = (bmag_max_global + bmag_min_global)/2.0;

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, 
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  gyrokinetic_multib_optimize_xpt_seams(mbinp, mbapp);

  // See the metric-continuity NOTE in gkyl_gyrokinetic_multib_app_new_geom's
  // identical sync sequence above (geometry-only path): the raw interior
  // geometry synced here is ghost-copied, not reconciled, across interfaces.
  //
  // Sync the conf-space volume Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs_vol[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_ghost;
  }
  // Sync across blocks.
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_vol, jacs_vol);

  // Sync the surface conf-space Jacobian, compute its reciprocal, and 
  // store its product with the Jacobian of this block (in the ghost cell).
  for (int d = 0; d<cdim; d++) {
    // Sync jacobgeo.
    struct gkyl_array *jacs[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      jacs[b] = geo_surf.jacobgeo_ratio;
      gkyl_array_copy_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range_to_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *jacgeo = mkarr(mbapp->use_gpu, geo_surf.jacobgeo_ratio->ncomp, geo_surf.jacobgeo_ratio->size);

      // Compute 1/jacobgeo in ghost cells.
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_lower_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_upper_ghost[d]);
      // Multiply by the Jacobian of this block.
      gkyl_array_copy_range_to_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_lower_ghost[d], &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_upper_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      // Set the ratio to 1 in the interior (shouldn't be in use).
      gkyl_array_clear_range(geo_surf.jacobgeo_ratio, 0.0, &sbapp->local);
      gkyl_array_shiftc_range(geo_surf.jacobgeo_ratio, pow(sqrt(2.0),sbapp->cdim), 0, &sbapp->local);

      gkyl_array_release(jacgeo);
    }

  }

  const struct gkyl_gk_block_geom_info *bgi0 = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, 0);

  if (cdim > 1 && bgi0->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    // Sync the surface deltats. Need to copy the upper ghost into the upper skin before syncing.
    for (int d = 0; d<cdim; d++) {
      struct gkyl_array *deltats[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
        deltats[b] = geo_surf.deltats;
        gkyl_array_copy_range_to_range(deltats[b], deltats[b], &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, deltats, deltats);
    }
    // Accumulate the appropriate shift
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      int par_dir = cdim-1;
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gkyl_array *delta_ts = sbapp->gk_geom->geo_surf[par_dir].deltats;
      struct gkyl_array *buffer = mkarr(sbapp->use_gpu, delta_ts->ncomp, delta_ts->size);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_upper_skin[par_dir], &sbapp->local_upper_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_upper_skin[par_dir]);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_lower_skin[par_dir], &sbapp->local_lower_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_lower_skin[par_dir]);

      gkyl_array_release(buffer);

      // Deflate delta_ts.
      gyrokinetic_deflate_delta_ts(sbapp, delta_ts);

      // Write the twistshift shift.
      gyrokinetic_app_write_ts_shift(sbapp);
    }
  }


  // Create the rest of the single-block solvers.
  for (int i=0; i<num_local_blocks; ++i)
    singleb_app_new_solver(mbinp, mbapp->local_blocks[i], mbapp, mbapp->singleb_apps[i]);

  // Create the MB field app.
  mbapp->field = gk_multib_field_new(mbinp, mbapp);


  // Create connections needed for syncing charged species phase-space quantities.
  mbapp->mbcc_sync_charged = gkyl_malloc(mbapp->num_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_species; ++i) {
    mbapp->mbcc_sync_charged[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_charged[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_species; ++i) {
      mbapp->mbcc_sync_charged[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_charged[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_charged[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_charged[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_species *gks = &sbapp->species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Create connections needed for syncing neutral species phase-space quantities.
  mbapp->mbcc_sync_neut = gkyl_malloc(mbapp->num_neut_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_neut_species; ++i) {
    mbapp->mbcc_sync_neut[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_neut[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_neut_species; ++i) {
      mbapp->mbcc_sync_neut[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_neut[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_neut[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_neut[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Sync the effective diffusivity of the anomalous diffusion operator.
  // Assume they either all have anomalous diffusion or none of them do.
  bool any_anomalous_diff = false;
  for (int i=0; i<mbapp->num_species; ++i) {
    bool has_anomalous_diff = true;
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      has_anomalous_diff = has_anomalous_diff && sbapp->species[i].anom_diff.anom_diff_id;
    }

    if (has_anomalous_diff) {
      any_anomalous_diff = true;
      // Divide diffD by dz before the transfer
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, sbapp->grid.dx[sbapp->cdim-1], &sbapp->global);
      }

      // Sync
      struct gkyl_array *gkad_nu[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        gkad_nu[b] = sbapp->species[i].anom_diff.diffD;
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, gkad_nu, gkad_nu);

      // Multiply by diffD dz after the transfer to achieve rescaling
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, 1.0/sbapp->grid.dx[sbapp->cdim-1], &sbapp->global_ext);
      }

    }
  }

  if (any_anomalous_diff) {
    // Sync the interior conf-space reciprocal Jacobian.
    struct gkyl_array *jacs_inv_vol[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      jacs_inv_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_inv;
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_inv_vol, jacs_inv_vol);
  }

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {};

  mbapp->dts = gkyl_dynvec_new(GKYL_DOUBLE, 1); // Dynvector to store time steps.
  mbapp->is_first_dt_write_call = true;

  gkyl_free(rank_list);
  gkyl_free(branks);

  // No partial grid: all local blocks and their hard material checks must
  // finish on every rank before any block's geometry is written.
  gkyl_comm_barrier(mbapp->comm);
  gkyl_gyrokinetic_multib_app_write_geometry(mbapp);
  return mbapp;
}

void
gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wtm = gkyl_wall_clock();
  // Compute fields.
  if (app->update_field) {
    // Solve the field equation.
    gk_multib_field_rhs(app, app->field, fin, bflux);
  }
  app->stat.field_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gyrokinetic_multib_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{
  // Apply boundary conditions in each block (including intrablock sync).
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      gk_species_apply_bc(sbapp, &sbapp->species[i], distf[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      gk_neut_species_apply_bc(sbapp, &sbapp->neut_species[i], distf_neut[li_neut+i]);
    }
  }

  // Sync blocks.
  struct timespec wst = gkyl_wall_clock();
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  for (int i=0; i<app->num_species; ++i) {
    // Sync charged species.
    struct gkyl_array *fs[app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      int li_charged = b * app->num_species;
      fs[b] = distf[li_charged+i];
    }
    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_charged[i].send, app->mbcc_sync_charged[i].recv, fs, fs);
  }
  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);

  struct timespec wst_neut = gkyl_wall_clock();
  struct gkyl_gyrokinetic_app *sbapp0 = app->singleb_apps[0];
  for (int i=0; i<app->num_neut_species; ++i) {
    // Sync neutral species.
    if (!sbapp0->neut_species[i].info.is_static) {
      struct gkyl_array *fs[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        int li_neut = b * app->num_neut_species;
        fs[b] = distf_neut[li_neut+i];
      }
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_neut[i].send, app->mbcc_sync_neut[i].recv, fs, fs);
    }
  }
  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst_neut);

}

void
gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[])
{
  // Compute fields and apply BCs.

  // Compute the field.
  // MF 2024/09/27/: Need the cast here for consistency. Fixing
  // this may require removing 'const' from a lot of places.
  gyrokinetic_multib_calc_field(app, tcurr, (const struct gkyl_array **) distf, bflux);

  // Apply boundary conditions.
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_multib_apply_bc(app, tcurr, distf, distf_neut);
  app->stat.bc_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_species(app, i, t0);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, i, t0);
  }  

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
  struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[li_charged+i] = sbapp->species[i].f;
      bflux[li_charged+i] = sbapp->species[i].bflux.f1;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[li_neut+i] = sbapp->neut_species[i].f;
    }
  }
  gyrokinetic_multib_calc_field_and_apply_bc(app, t0, distf, bflux, distf_neut);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_neut_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    for (int i=0; i<app->num_neut_species; i++) {
      int neut_frame = frame;
      if (sbapp->neut_species[i].info.is_static) {
        neut_frame = 0;
      }
      rstat = gkyl_gyrokinetic_app_from_frame_neut_species(sbapp, i, neut_frame);
    }
    for (int i=0; i<app->num_species; i++) {
      rstat = gkyl_gyrokinetic_app_from_frame_species(sbapp, i, frame);
    }
  }
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species * app->num_local_blocks];
    struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
    struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int li_charged = b * app->num_species;
      int li_neut = b * app->num_neut_species;
      for (int i=0; i<app->num_species; ++i) {
        distf[li_charged+i] = sbapp->species[i].f;
        bflux[li_charged+i] = sbapp->species[i].bflux.f;
      }
      for (int i=0; i<app->num_neut_species; ++i) {
        distf_neut[li_neut+i] = sbapp->neut_species[i].f;
      }
    }
//    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
//      for (int i=0; i<app->num_species; ++i) {
//        struct gk_species *s = &app->species[i];
//
//        // Compute advection speeds so we can compute the initial boundary flux.
//        gkyl_gk_collisionless_flux_alpha_surf(s->calc_gk_vars, 
//          &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
//          s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
//
//        // Compute and store (in the ghost cell of of out) the boundary fluxes.
//        // NOTE: this overwrites ghost cells that may be used for sourcing.
//        gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
//      }
//    }
    gyrokinetic_multib_calc_field_and_apply_bc(app, rstat.stime, distf, bflux, distf_neut);
  }

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  sbapp->field->is_first_energy_write_call = false; // Append to existing diagnostic.
                                                    //
  app->is_first_dt_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

// private function to handle variable argument list for printing
static void
v_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gyrokinetic_multib_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-block_topo.gkyl", app->name);
    gkyl_block_topo_write(app->block_topo, file_name.str);
    cstr_drop(&file_name);
  }
}

static struct gkyl_gk_geometry_inp
gyrokinetic_multib_geometry_inp_from_block(const struct gkyl_gk_block_geom_info *bgi,
  struct gkyl_comm *comm)
{
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = bgi->geometry.geometry_id,
    .c2p_ctx = bgi->geometry.c2p_ctx,
    .mapc2p = bgi->geometry.mapc2p,
    .bfield_ctx = bgi->geometry.bfield_ctx,
    .bfield_func = bgi->geometry.bfield_func,
    .efit_info = bgi->geometry.efit_info,
    .tok_grid_info = bgi->geometry.tok_grid_info,
    .mirror_grid_info = bgi->geometry.mirror_grid_info,
    .position_map = 0,
    .comm = comm,
    .world = {
      bgi->geometry.world[0],
      bgi->geometry.world[1],
      bgi->geometry.world[2],
    },
    .has_LCFS = bgi->geometry.has_LCFS,
    .x_LCFS = bgi->geometry.x_LCFS,
  };
  strncpy(geometry_inp.geometry_path, bgi->geometry.geometry_path,
    sizeof(geometry_inp.geometry_path)-1);
  geometry_inp.geometry_path[sizeof(geometry_inp.geometry_path)-1] = '\0';
  return geometry_inp;
}

static bool
jacobgeo_ratio_diag_coeff_is_set(double val)
{
  union { double d; uint64_t u; } bits = { .d = val };
  union { double d; uint64_t u; } huge = { .d = 1.0e300 };
  uint64_t abs_bits = bits.u & UINT64_C(0x7fffffffffffffff);
  return abs_bits < huge.u;
}

static void
gyrokinetic_multib_app_write_jacobgeo_ratio_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-jacobgeo_ratio_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-jacobgeo_ratio_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,num_cells,raw_coeff0_min,raw_coeff0_max,raw_coeff0_mean,ratio_min,ratio_max,ratio_mean\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const double surf_modal_one = pow(sqrt(2.0), cdim-1);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        double raw_min = DBL_MAX, raw_max = -DBL_MAX, raw_sum = 0.0;
        double ratio_min = DBL_MAX, ratio_max = -DBL_MAX, ratio_sum = 0.0;
        long count = 0;

        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, range);
        while (gkyl_range_iter_next(&iter)) {
          long loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
          const double *val = gkyl_array_cfetch(sbapp->gk_geom->geo_surf[d].jacobgeo_ratio, loc);
          double raw = val[0];
          if (!jacobgeo_ratio_diag_coeff_is_set(raw))
            continue;
          double ratio = raw/surf_modal_one;
          raw_min = GKYL_MIN2(raw_min, raw);
          raw_max = GKYL_MAX2(raw_max, raw);
          raw_sum += raw;
          ratio_min = GKYL_MIN2(ratio_min, ratio);
          ratio_max = GKYL_MAX2(ratio_max, ratio);
          ratio_sum += ratio;
          count += 1;
        }

        if (count > 0) {
          fprintf(fp, "%s,%d,%d,%d,%s,%ld,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
            app->name, rank, bid, d, edge, count,
            raw_min, raw_max, raw_sum/count,
            ratio_min, ratio_max, ratio_sum/count);
        }
      }
    }
  }

  fclose(fp);
}

static void
gyrokinetic_multib_app_write_flux_weight_rescale_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-flux_weight_rescale_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-flux_weight_rescale_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,quantity,num_cells,coeff0_raw_min,coeff0_raw_max,coeff0_raw_mean,coeff0_ratio_min,coeff0_ratio_max,coeff0_ratio_mean,coeff0_ratio_scaled_min,coeff0_ratio_scaled_max,coeff0_ratio_scaled_mean\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const double surf_modal_one = pow(sqrt(2.0), cdim-1);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      struct {
        const char *name;
        const struct gkyl_array *array;
      } quantities[] = {
        { "jacobtot_inv", geo_surf->jacobtot_inv },
        { "B3", geo_surf->B3 },
        { "lenr", geo_surf->lenr },
        { "cmag", geo_surf->cmag },
      };

      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
          int surf_ncomp = geo_surf->jacobgeo_ratio->ncomp;
          int side_offset = quantities[q].array->ncomp == 2*surf_ncomp
            ? (e == 0 ? 0 : surf_ncomp) : 0;
          double raw_min = DBL_MAX, raw_max = -DBL_MAX, raw_sum = 0.0;
          double ratio_min = DBL_MAX, ratio_max = -DBL_MAX, ratio_sum = 0.0;
          double scaled_min = DBL_MAX, scaled_max = -DBL_MAX, scaled_sum = 0.0;
          long count = 0;

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;
            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_val = gkyl_array_cfetch(quantities[q].array, skin_loc);
            const double *ratio_val = gkyl_array_cfetch(geo_surf->jacobgeo_ratio, ghost_loc);
            double raw = raw_val[side_offset];
            double ratio_raw = ratio_val[0];
            if (!jacobgeo_ratio_diag_coeff_is_set(raw) || !jacobgeo_ratio_diag_coeff_is_set(ratio_raw))
              continue;

            double ratio = ratio_raw/surf_modal_one;
            double scaled = raw*ratio;
            raw_min = GKYL_MIN2(raw_min, raw);
            raw_max = GKYL_MAX2(raw_max, raw);
            raw_sum += raw;
            ratio_min = GKYL_MIN2(ratio_min, ratio);
            ratio_max = GKYL_MAX2(ratio_max, ratio);
            ratio_sum += ratio;
            scaled_min = GKYL_MIN2(scaled_min, scaled);
            scaled_max = GKYL_MAX2(scaled_max, scaled);
            scaled_sum += scaled;
            count += 1;
          }

          if (count > 0) {
            fprintf(fp, "%s,%d,%d,%d,%s,%s,%ld,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
              app->name, rank, bid, d, edge, quantities[q].name, count,
              raw_min, raw_max, raw_sum/count,
              ratio_min, ratio_max, ratio_sum/count,
              scaled_min, scaled_max, scaled_sum/count);
          }
        }
      }
    }
  }

  fclose(fp);
}

static void
gyrokinetic_multib_app_write_flux_weight_modal_product_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-flux_weight_modal_product_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-flux_weight_modal_product_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,quantity,cell_ord,idx0,idx1,idx2,coeff,side_offset,raw_side_coeff,ratio_coeff,product_coeff\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      int surf_ncomp = geo_surf->jacobgeo_ratio->ncomp;
      struct {
        const char *name;
        const struct gkyl_array *array;
      } quantities[] = {
        { "jacobtot_inv", geo_surf->jacobtot_inv },
        { "B3", geo_surf->B3 },
        { "lenr", geo_surf->lenr },
        { "cmag", geo_surf->cmag },
      };

      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
          int side_offset = quantities[q].array->ncomp == 2*surf_ncomp
            ? (e == 0 ? 0 : surf_ncomp) : 0;
          struct gkyl_array *raw_side = mkarr(app->use_gpu, surf_ncomp, quantities[q].array->size);
          struct gkyl_array *product = mkarr(app->use_gpu, surf_ncomp, quantities[q].array->size);
          gkyl_array_clear(raw_side, 0.0);
          gkyl_array_clear(product, 0.0);

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;
            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_val = gkyl_array_cfetch(quantities[q].array, skin_loc);
            double *raw_side_val = gkyl_array_fetch(raw_side, ghost_loc);
            for (int k=0; k<surf_ncomp; ++k)
              raw_side_val[k] = raw_val[side_offset+k];
          }

          gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, product,
            0, raw_side, 0, geo_surf->jacobgeo_ratio, range);

          long cell_ord = 0;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_side_val = gkyl_array_cfetch(raw_side, ghost_loc);
            const double *ratio_val = gkyl_array_cfetch(geo_surf->jacobgeo_ratio, ghost_loc);
            const double *product_val = gkyl_array_cfetch(product, ghost_loc);
            for (int k=0; k<surf_ncomp; ++k) {
              if (!jacobgeo_ratio_diag_coeff_is_set(raw_side_val[k]) ||
                  !jacobgeo_ratio_diag_coeff_is_set(ratio_val[k]) ||
                  !jacobgeo_ratio_diag_coeff_is_set(product_val[k]))
                continue;
              fprintf(fp, "%s,%d,%d,%d,%s,%s,%ld,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e\n",
                app->name, rank, bid, d, edge, quantities[q].name, cell_ord,
                iter.idx[0], iter.idx[1], iter.idx[2], k, side_offset,
                raw_side_val[k], ratio_val[k], product_val[k]);
            }
            cell_ord += 1;
          }

          gkyl_array_release(raw_side);
          gkyl_array_release(product);
        }
      }
    }
  }

  fclose(fp);
}

static const char*
edge_name_from_index(int edge)
{
  return edge == 0 ? "lower" : "upper";
}

static const char*
oriented_edge_name(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE: return "lower_positive";
  case GKYL_LOWER_NEGATIVE: return "lower_negative";
  case GKYL_UPPER_POSITIVE: return "upper_positive";
  case GKYL_UPPER_NEGATIVE: return "upper_negative";
  case GKYL_PHYSICAL: return "physical";
  }
  return "unknown";
}

static int
oriented_edge_index(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE:
  case GKYL_LOWER_NEGATIVE:
    return 0;
  case GKYL_UPPER_POSITIVE:
  case GKYL_UPPER_NEGATIVE:
    return 1;
  case GKYL_PHYSICAL:
    return -1;
  }
  return -1;
}

static const char*
orientation_name(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE:
  case GKYL_UPPER_POSITIVE:
    return "positive";
  case GKYL_LOWER_NEGATIVE:
  case GKYL_UPPER_NEGATIVE:
    return "negative";
  case GKYL_PHYSICAL:
    return "physical";
  }
  return "unknown";
}

static int
orientation_sign(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE:
  case GKYL_UPPER_POSITIVE:
    return 1;
  case GKYL_LOWER_NEGATIVE:
  case GKYL_UPPER_NEGATIVE:
    return -1;
  case GKYL_PHYSICAL:
    return 0;
  }
  return 0;
}

enum pointwise_geometry_field {
  POINTWISE_VALID,
  POINTWISE_R,
  POINTWISE_Z,
  POINTWISE_PHI,
  POINTWISE_PSI_REQUESTED,
  POINTWISE_PSI_EVALUATED,
  POINTWISE_PSI_RESIDUAL,
  POINTWISE_EPSI_R,
  POINTWISE_EPSI_PHI,
  POINTWISE_EPSI_Z,
  POINTWISE_ETHETA_R,
  POINTWISE_ETHETA_PHI,
  POINTWISE_ETHETA_Z,
  POINTWISE_JACOBIAN_SIGNED,
  POINTWISE_JACOBIAN_ABS,
  POINTWISE_NUM_FIELDS
};

// Canonical partner convention for a same-direction 2D connection:
//
//   local lower -> partner upper     local upper -> partner lower
//   positive: (cell,q) -> (cell,q), orientation sign +1
//   negative: (cell,q) -> (Ncell-1-cell,p-q), orientation sign -1
//
// Lower/upper selects the partner face but does not negate the coordinate
// derivative.  A negative connection reverses the one tangential logical
// coordinate.  Thus a psi-normal interface reverses partner e_theta, while a
// theta-normal interface reverses partner e_psi.  In either case the partner
// signed Jacobian receives the same orientation sign; |J| is unchanged.
static void
pointwise_canonicalize_partner_record(int interface_dir, int orient_sign,
  const double *native, double *canonical)
{
  memcpy(canonical, native, sizeof(double[POINTWISE_NUM_FIELDS]));
  if (orient_sign >= 0)
    return;

  int first = interface_dir == 0 ? POINTWISE_ETHETA_R : POINTWISE_EPSI_R;
  int last = interface_dir == 0 ? POINTWISE_ETHETA_Z : POINTWISE_EPSI_Z;
  for (int f=first; f<=last; ++f)
    canonical[f] = -canonical[f];
  canonical[POINTWISE_JACOBIAN_SIGNED] =
    -canonical[POINTWISE_JACOBIAN_SIGNED];
}

static int
pointwise_int_pow(int base, int exponent)
{
  int result = 1;
  for (int i=0; i<exponent; ++i)
    result *= base;
  return result;
}

static struct gkyl_array*
pointwise_host_copy(const struct gkyl_array *src)
{
  struct gkyl_array *dest = gkyl_array_new(src->type, src->ncomp, src->size);
  gkyl_array_copy(dest, src);
  return dest;
}

// Evaluate the same equilibrium representation selected by tok_geo.c.  In
// particular, do not use the cubic interpolant when the map was constructed
// from the quadratic/DG representation.
static double
pointwise_eval_psi(const struct gkyl_efit *efit, bool use_cubics,
  double R, double Z)
{
  if (!efit || !jacobgeo_ratio_diag_coeff_is_set(R)
      || !jacobgeo_ratio_diag_coeff_is_set(Z))
    return 0.0;

  if (use_cubics) {
    double xn[2] = { R, Z }, psi = 0.0;
    efit->evf->eval_cubic(0.0, xn, &psi, efit->evf->ctx);
    return psi;
  }

  int idx[2];
  idx[0] = GKYL_MIN2(efit->rzlocal.upper[0], GKYL_MAX2(efit->rzlocal.lower[0],
    efit->rzlocal.lower[0]
      +(int) floor((R-efit->rzgrid.lower[0])/efit->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(efit->rzlocal.upper[1], GKYL_MAX2(efit->rzlocal.lower[1],
    efit->rzlocal.lower[1]
      +(int) floor((Z-efit->rzgrid.lower[1])/efit->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&efit->rzlocal, idx);
  const double *coeffs = gkyl_array_cfetch(efit->psizr, loc);
  double xc[2], eta[2];
  gkyl_rect_grid_cell_center(&efit->rzgrid, idx, xc);
  eta[0] = (R-xc[0])/(0.5*efit->rzgrid.dx[0]);
  eta[1] = (Z-xc[1])/(0.5*efit->rzgrid.dx[1]);
  return efit->rzbasis.eval_expand(eta, coeffs);
}

// Recover the radial logical coordinate used by gkyl_tok_geo_calc_surface.
// The current surface-node construction is reliable for p1, which is also
// the polynomial order used by the multiblock X-point cases diagnosed here.
static double
pointwise_requested_psi(const struct gkyl_gyrokinetic_app *sbapp,
  int dir, const int *surf_idx)
{
  int poly_order = sbapp->basis.poly_order;
  if (poly_order != 1)
    return 0.0;

  double dx = sbapp->grid.dx[0];
  double logical_lower = sbapp->grid.lower[0]
    +(sbapp->local.lower[0]-sbapp->global.lower[0])*dx;
  double logical_psi;
  if (dir == 0) {
    int node = surf_idx[0]-sbapp->gk_geom->nrange_surf[dir].lower[0];
    logical_psi = logical_lower+node*dx;
  }
  else {
    int node = surf_idx[0]-sbapp->gk_geom->nrange_surf[dir].lower[0];
    int cell = node/2, quad_node = node%2;
    double eta = quad_node == 0 ? -1.0/sqrt(3.0) : 1.0/sqrt(3.0);
    logical_psi = logical_lower+(cell+0.5*(1.0+eta))*dx;
  }

  double requested_psi = 0.0;
  sbapp->position_map->maps[0](0.0, &logical_psi, &requested_psi,
    sbapp->position_map->ctxs[0]);
  return requested_psi;
}

static void
pointwise_fill_geometry_record(const struct gkyl_gyrokinetic_app *sbapp,
  const struct gkyl_efit *efit, bool use_cubics, int dir,
  const int *surf_idx, const struct gkyl_array *mc2p_nodal_fd,
  const struct gkyl_array *dxdz_nodal,
  const struct gkyl_array *jacobgeo_signed_nodal,
  const struct gkyl_array *jacobgeo_nodal, double *record)
{
  long loc = gkyl_range_idx(&sbapp->gk_geom->nrange_surf[dir], surf_idx);
  const double *position = gkyl_array_cfetch(mc2p_nodal_fd, loc);
  const double *tangents = gkyl_array_cfetch(dxdz_nodal, loc);
  const double *jacobian_signed = gkyl_array_cfetch(jacobgeo_signed_nodal, loc);
  const double *jacobian_abs = gkyl_array_cfetch(jacobgeo_nodal, loc);

  double R = position[0], Z = position[1], phi = position[2];
  double cos_phi = cos(phi), sin_phi = sin(phi);
  double requested_psi = pointwise_requested_psi(sbapp, dir, surf_idx);
  double evaluated_psi = pointwise_eval_psi(efit, use_cubics, R, Z);

  record[POINTWISE_R] = R;
  record[POINTWISE_Z] = Z;
  record[POINTWISE_PHI] = phi;
  record[POINTWISE_PSI_REQUESTED] = requested_psi;
  record[POINTWISE_PSI_EVALUATED] = evaluated_psi;
  record[POINTWISE_PSI_RESIDUAL] = evaluated_psi-requested_psi;

  // dxdz_nodal stores Cartesian e_psi in components 0:2 and Cartesian
  // e_theta in components 6:8.  Rotate those physical vectors to the
  // cylindrical basis at the stored phi before writing them.
  record[POINTWISE_EPSI_R] = cos_phi*tangents[0]+sin_phi*tangents[1];
  record[POINTWISE_EPSI_PHI] = -sin_phi*tangents[0]+cos_phi*tangents[1];
  record[POINTWISE_EPSI_Z] = tangents[2];
  record[POINTWISE_ETHETA_R] = cos_phi*tangents[6]+sin_phi*tangents[7];
  record[POINTWISE_ETHETA_PHI] = -sin_phi*tangents[6]+cos_phi*tangents[7];
  record[POINTWISE_ETHETA_Z] = tangents[8];
  record[POINTWISE_JACOBIAN_SIGNED] = jacobian_signed[0];
  record[POINTWISE_JACOBIAN_ABS] = jacobian_abs[0];

  bool valid = sbapp->basis.poly_order == 1;
  for (int f=POINTWISE_R; f<POINTWISE_NUM_FIELDS; ++f)
    valid = valid && jacobgeo_ratio_diag_coeff_is_set(record[f]);
  record[POINTWISE_VALID] = valid ? 1.0 : 0.0;
}

static void
pointwise_fprint_geometry_record(FILE *fp, const double *record)
{
  for (int i=0; i<POINTWISE_NUM_FIELDS; ++i)
    fprintf(fp, ",%.17e", record[i]);
}

static void
pointwise_write_geometry_header(FILE *fp)
{
  fprintf(fp, "app,rank,block,dir,edge,partner_block,partner_dir,partner_edge,orientation,comparison_status,orientation_transform,interface_cell,node_ordinal,interface_node_index,partner_interface_cell,partner_node_ordinal,partner_interface_node_index,cell_idx0,cell_idx1,cell_idx2,surface_idx0,surface_idx1,surface_idx2,quad_idx0,quad_idx1,quad_idx2");
  const char *prefixes[] = { "local", "partner_native", "partner" };
  const char *fields[] = {
    "valid", "R", "Z", "phi", "psi_requested", "psi_evaluated", "psi_residual",
    "e_psi_R", "e_psi_phi", "e_psi_Z", "e_theta_R", "e_theta_phi",
    "e_theta_Z", "jacobian_signed", "jacobian_abs"
  };
  for (int p=0; p<3; ++p)
    for (int f=0; f<POINTWISE_NUM_FIELDS; ++f)
      fprintf(fp, ",%s_%s", prefixes[p], fields[f]);
  fputc('\n', fp);
}

// Return whether logical-theta edge (0 or 1) of this block is a candidate
// X-point-adjacent seam under the relaxed-seam workflow. Most ftypes have
// exactly one such edge (the other is a physical/plate boundary); a few
// double-null ftypes (CORE_L/CORE_R, the DN SOL "MID" blocks that sit
// between both X-points) and the self-periodic single-null CORE block are
// X-point-adjacent at both edges. At most one of a doubly-adjacent block's
// two candidate edges can actually be claimed by a pair in a given geometry
// build (xpt_optimizer_discover_pairs claims at most one pair per block);
// reporting both here and letting ascending (bid, edge) discovery order
// resolve the conflict is deterministic, not a partial/best-effort guess.
// No longer requires half_domain: the relaxed_xpt_seam displacement is now
// wired into both the half-domain and full-domain (extended) construction
// paths in tok_geo.c, so this recognizes seams under either.
static bool
pointwise_xpt_seam_edge_active(const struct gkyl_gk_block_geom_info *bgi,
  int edge)
{
  if (bgi->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
    return false;
  if (edge != 0 && edge != 1)
    return false;
  const struct gkyl_tok_geo_grid_inp *inp = &bgi->geometry.tok_grid_info;
  if (!inp->straight_xpt_ray || !inp->relaxed_xpt_seam)
    return false;

  switch (inp->ftype) {
    // Double-null lower-X-point-adjacent single-seam ftypes.
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      return edge == 1;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      return edge == 0;
    // Double-null upper-X-point-adjacent single-seam ftypes.
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
      return edge == 0;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      return edge == 1;
    // Double-null blocks adjacent to both X-points (one edge each).
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
      return true;
    // Self-periodic single-null core: both edges are the same
    // X-point-adjacent periodic wrap (see tok_relax_xpt_seam_point's
    // uniform-weight handling of this ftype).
    case GKYL_GEOMETRY_TOKAMAK_CORE:
      return true;
    default:
      return false;
  }
}

static bool
pointwise_xpt_ftype_pair(enum gkyl_tok_geo_type local,
  enum gkyl_tok_geo_type partner)
{
  return
    ((local == GKYL_GEOMETRY_TOKAMAK_CORE_L &&
        partner == GKYL_GEOMETRY_TOKAMAK_CORE_R) ||
      (local == GKYL_GEOMETRY_TOKAMAK_CORE_R &&
        partner == GKYL_GEOMETRY_TOKAMAK_CORE_L)) ||
    (local == GKYL_GEOMETRY_TOKAMAK_CORE &&
      partner == GKYL_GEOMETRY_TOKAMAK_CORE) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_PF_LO_R &&
        partner == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) ||
      (local == GKYL_GEOMETRY_TOKAMAK_PF_LO_L &&
        partner == GKYL_GEOMETRY_TOKAMAK_PF_LO_R)) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_PF_UP_R &&
        partner == GKYL_GEOMETRY_TOKAMAK_PF_UP_L) ||
      (local == GKYL_GEOMETRY_TOKAMAK_PF_UP_L &&
        partner == GKYL_GEOMETRY_TOKAMAK_PF_UP_R)) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) ||
      (local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO)) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) ||
      (local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID)) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) ||
      (local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)) ||
    ((local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) ||
      (local == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID &&
        partner == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP));
}

static bool
pointwise_interface_is_owner(int bid, int dir, int edge, int partner_bid,
  int partner_dir, int partner_edge)
{
  if (bid != partner_bid)
    return bid < partner_bid;
  if (dir != partner_dir)
    return dir < partner_dir;
  return edge < partner_edge;
}

struct pointwise_xpt_mismatch {
  bool local_valid, partner_valid;
  bool local_abs_signed_consistent, partner_abs_signed_consistent;
  bool objective_sample_valid, geometry_constraints_valid;
  double position, position_scale;
  double local_flux_residual, partner_flux_residual;
  double requested_flux_mismatch;
  double local_signed_j, local_abs_j, partner_signed_j, partner_abs_j;
  double epsi_rz, epsi_phi_raw;
  double etheta_direction, etheta_magnitude, etheta_phi_raw;
  double jacobian_magnitude;
};

static struct pointwise_xpt_mismatch
pointwise_measure_xpt_mismatch(const double *local, const double *partner)
{
  struct pointwise_xpt_mismatch mismatch = { 0 };
  mismatch.local_valid = local[POINTWISE_VALID] == 1.0;
  mismatch.partner_valid = partner[POINTWISE_VALID] == 1.0;
  if (!mismatch.local_valid || !mismatch.partner_valid)
    return mismatch;

  double dr = local[POINTWISE_R]-partner[POINTWISE_R];
  double dz = local[POINTWISE_Z]-partner[POINTWISE_Z];
  mismatch.position = hypot(dr, dz);
  mismatch.position_scale = fmax(1.0,
    fmax(hypot(local[POINTWISE_R], local[POINTWISE_Z]),
      hypot(partner[POINTWISE_R], partner[POINTWISE_Z])));
  mismatch.local_flux_residual = local[POINTWISE_PSI_RESIDUAL];
  mismatch.partner_flux_residual = partner[POINTWISE_PSI_RESIDUAL];
  mismatch.requested_flux_mismatch = fabs(local[POINTWISE_PSI_REQUESTED]
    -partner[POINTWISE_PSI_REQUESTED]);

  mismatch.local_signed_j = local[POINTWISE_JACOBIAN_SIGNED];
  mismatch.local_abs_j = local[POINTWISE_JACOBIAN_ABS];
  mismatch.partner_signed_j = partner[POINTWISE_JACOBIAN_SIGNED];
  mismatch.partner_abs_j = partner[POINTWISE_JACOBIAN_ABS];
  mismatch.local_abs_signed_consistent = mismatch.local_abs_j > 0.0 &&
    mismatch.local_signed_j != 0.0 &&
    fabs(mismatch.local_abs_j-fabs(mismatch.local_signed_j))
      <= 1e-12*fmax(1.0, mismatch.local_abs_j);
  mismatch.partner_abs_signed_consistent = mismatch.partner_abs_j > 0.0 &&
    mismatch.partner_signed_j != 0.0 &&
    fabs(mismatch.partner_abs_j-fabs(mismatch.partner_signed_j))
      <= 1e-12*fmax(1.0, mismatch.partner_abs_j);

  double local_epsi = hypot(local[POINTWISE_EPSI_R],
    local[POINTWISE_EPSI_Z]);
  double partner_epsi = hypot(partner[POINTWISE_EPSI_R],
    partner[POINTWISE_EPSI_Z]);
  double epsi_scale = hypot(local_epsi, partner_epsi)/sqrt(2.0);
  double local_etheta = hypot(local[POINTWISE_ETHETA_R],
    local[POINTWISE_ETHETA_Z]);
  double partner_etheta = hypot(partner[POINTWISE_ETHETA_R],
    partner[POINTWISE_ETHETA_Z]);
  double etheta_scale = hypot(local_etheta, partner_etheta)/sqrt(2.0);
  mismatch.objective_sample_valid =
    isfinite(epsi_scale) && epsi_scale > 0.0 &&
    isfinite(etheta_scale) && etheta_scale > 0.0 &&
    isfinite(local_etheta) && local_etheta > 0.0 &&
    isfinite(partner_etheta) && partner_etheta > 0.0 &&
    isfinite(mismatch.local_abs_j) && mismatch.local_abs_j > 0.0 &&
    isfinite(mismatch.partner_abs_j) && mismatch.partner_abs_j > 0.0;
  if (mismatch.objective_sample_valid) {
    double depsi_r = local[POINTWISE_EPSI_R]
      -partner[POINTWISE_EPSI_R];
    double depsi_z = local[POINTWISE_EPSI_Z]
      -partner[POINTWISE_EPSI_Z];
    double depsi_phi = local[POINTWISE_EPSI_PHI]
      -partner[POINTWISE_EPSI_PHI];
    double depsi_r_normalized = depsi_r/epsi_scale;
    double depsi_z_normalized = depsi_z/epsi_scale;
    double depsi_phi_normalized = depsi_phi/epsi_scale;
    mismatch.epsi_rz = depsi_r_normalized*depsi_r_normalized
      +depsi_z_normalized*depsi_z_normalized;
    mismatch.epsi_phi_raw = depsi_phi_normalized*depsi_phi_normalized;

    double cos_angle =
      (local[POINTWISE_ETHETA_R]/local_etheta)
        *(partner[POINTWISE_ETHETA_R]/partner_etheta)
      +(local[POINTWISE_ETHETA_Z]/local_etheta)
        *(partner[POINTWISE_ETHETA_Z]/partner_etheta);
    cos_angle = fmax(-1.0, fmin(1.0, cos_angle));
    mismatch.etheta_direction = 2.0*(1.0-cos_angle);
    double log_etheta_ratio = log(local_etheta)-log(partner_etheta);
    mismatch.etheta_magnitude = log_etheta_ratio*log_etheta_ratio;
    double detheta_phi = local[POINTWISE_ETHETA_PHI]
      -partner[POINTWISE_ETHETA_PHI];
    double detheta_phi_normalized = detheta_phi/etheta_scale;
    mismatch.etheta_phi_raw =
      detheta_phi_normalized*detheta_phi_normalized;
    double log_j_ratio = log(mismatch.local_abs_j)
      -log(mismatch.partner_abs_j);
    mismatch.jacobian_magnitude = log_j_ratio*log_j_ratio;
    mismatch.objective_sample_valid = isfinite(mismatch.epsi_rz) &&
      isfinite(mismatch.epsi_phi_raw) &&
      isfinite(mismatch.etheta_direction) &&
      isfinite(mismatch.etheta_magnitude) &&
      isfinite(mismatch.etheta_phi_raw) &&
      isfinite(mismatch.jacobian_magnitude);
  }

  double flux_scale = fmax(1.0,
    fmax(fabs(local[POINTWISE_PSI_REQUESTED]),
      fabs(partner[POINTWISE_PSI_REQUESTED])));
  mismatch.geometry_constraints_valid =
    mismatch.position <= 1e-10*mismatch.position_scale &&
    fabs(mismatch.local_flux_residual) <= 1e-10*flux_scale &&
    fabs(mismatch.partner_flux_residual) <= 1e-10*flux_scale &&
    mismatch.requested_flux_mismatch <= 1e-12*flux_scale &&
    mismatch.local_abs_signed_consistent &&
    mismatch.partner_abs_signed_consistent;
  return mismatch;
}

static void
pointwise_write_xpt_objective_header(FILE *fp)
{
  fprintf(fp, "app,rank,block,dir,edge,ftype,partner_block,partner_dir,partner_edge,partner_ftype,orientation,comparison_status,interface_node_index,global_interface_node_count,sample_weight,requested_mode,requested_extension,local_relaxed_enabled,partner_relaxed_enabled,local_sweep_enabled,partner_sweep_enabled,local_coefficient_m,partner_coefficient_m,local_bound_m,partner_bound_m,topology_pair_valid,transfer_valid,canonicalization_valid,parameter_valid,local_valid,partner_valid,position_m,position_scale,local_flux_residual,partner_flux_residual,requested_flux_mismatch,local_signed_j,local_abs_j,partner_signed_j,partner_abs_j,local_abs_signed_consistent,partner_abs_signed_consistent,epsi_rz_normalized_sq,epsi_phi_raw_normalized_sq,etheta_direction_normalized_sq,etheta_log_magnitude_sq,etheta_phi_raw_normalized_sq,jacobian_log_magnitude_sq,smoothness_dq_integral_m2,smoothness_d2q_integral_m2,objective_sample_valid,pointwise_constraints_valid,endpoint_constraint_status,branch_constraint_status,radial_ordering_constraint_status,block_jacobian_constraint_status\n");
}

static void
pointwise_write_xpt_objective_node(FILE *fp, const char *app_name,
  int rank, int bid, int edge, int partner_bid, int partner_dir,
  enum gkyl_oriented_edge partner_oriented_edge,
  const char *comparison_status, int transfer_status,
  bool can_canonicalize, long interface_node, long interface_node_count,
  const struct gkyl_gk_block_geom_info *local_bgi,
  const struct gkyl_gk_block_geom_info *partner_bgi,
  const double *local, const double *partner)
{
  const struct gkyl_tok_geo_grid_inp *local_inp =
    &local_bgi->geometry.tok_grid_info;
  const struct gkyl_tok_geo_grid_inp *partner_inp =
    &partner_bgi->geometry.tok_grid_info;
  int partner_edge = oriented_edge_index(partner_oriented_edge);
  bool topology_pair_valid =
    local_bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK &&
    partner_bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK &&
    pointwise_xpt_seam_edge_active(local_bgi, edge) &&
    pointwise_xpt_seam_edge_active(partner_bgi, partner_edge) &&
    partner_dir == 1 &&
    pointwise_xpt_ftype_pair(local_inp->ftype, partner_inp->ftype);

  double coefficient_scale = fmax(1.0,
    fmax(fabs(local_inp->relaxed_xpt_seam_delta_s_coeff),
      fabs(partner_inp->relaxed_xpt_seam_delta_s_coeff)));
  double bound_scale = fmax(1.0,
    fmax(fabs(local_inp->relaxed_xpt_seam_delta_s_bound),
      fabs(partner_inp->relaxed_xpt_seam_delta_s_bound)));
  bool parameters_match =
    fabs(local_inp->relaxed_xpt_seam_delta_s_coeff
      -partner_inp->relaxed_xpt_seam_delta_s_coeff)
      <= 64.0*DBL_EPSILON*coefficient_scale &&
    fabs(local_inp->relaxed_xpt_seam_delta_s_bound
      -partner_inp->relaxed_xpt_seam_delta_s_bound)
      <= 64.0*DBL_EPSILON*bound_scale &&
    local_inp->relaxed_xpt_seam == partner_inp->relaxed_xpt_seam &&
    local_inp->relaxed_xpt_seam_sweep ==
      partner_inp->relaxed_xpt_seam_sweep;
  double coefficient = local_inp->relaxed_xpt_seam_delta_s_coeff;
  double bound = local_inp->relaxed_xpt_seam_delta_s_bound;
  bool requested_nonzero = coefficient != 0.0;
  bool parameter_valid = parameters_match && isfinite(coefficient) &&
    isfinite(bound) && local_inp->relaxed_xpt_seam &&
    ((!requested_nonzero && bound >= 0.0) ||
      (bound > 0.0 && fabs(coefficient) <=
        bound*(1.0+64.0*DBL_EPSILON) &&
        local_inp->relaxed_xpt_seam_sweep));
  const char *requested_mode = requested_nonzero
    ? "requested_candidate" : "requested_zero";

  // Reuse the same per-node mismatch helper the diagnostic optimizer uses
  // (pointwise_measure_xpt_mismatch), so this CSV and the optimizer's
  // in-memory objective cannot silently drift apart.
  struct pointwise_xpt_mismatch mismatch =
    pointwise_measure_xpt_mismatch(local, partner);
  bool local_valid = mismatch.local_valid;
  bool partner_valid = mismatch.partner_valid;
  bool pointwise_constraints_valid = topology_pair_valid &&
    transfer_status == 0 && can_canonicalize && parameter_valid &&
    mismatch.geometry_constraints_valid;

  double sample_weight = interface_node_count > 0
    ? 1.0/interface_node_count : 0.0;
  // B1(q)=4q(1-q): these unweighted analytic integrals are independent of
  // the displacement bound. The bound is a hard constraint, not a scale.
  double smoothness_dq = (16.0/3.0)*coefficient*coefficient;
  double smoothness_d2q = 64.0*coefficient*coefficient;

  fprintf(fp, "%s,%d,%d,1,%s,%d,%d,%d,%s,%d,%s,%s,%ld,%ld,%.17e,%s,requested_fixed_edge_linear_arc_blend",
    app_name, rank, bid, edge_name_from_index(edge), local_inp->ftype,
    partner_bid, partner_dir, oriented_edge_name(partner_oriented_edge),
    partner_inp->ftype, orientation_name(partner_oriented_edge),
    comparison_status, interface_node, interface_node_count, sample_weight,
    requested_mode);
  fprintf(fp, ",%d,%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%d,%d,%d,%d",
    local_inp->relaxed_xpt_seam,
    partner_inp->relaxed_xpt_seam, local_inp->relaxed_xpt_seam_sweep,
    partner_inp->relaxed_xpt_seam_sweep, coefficient,
    partner_inp->relaxed_xpt_seam_delta_s_coeff, bound,
    partner_inp->relaxed_xpt_seam_delta_s_bound, topology_pair_valid,
    transfer_status == 0, can_canonicalize, parameter_valid);
  fprintf(fp, ",%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d",
    local_valid, partner_valid, mismatch.position, mismatch.position_scale,
    mismatch.local_flux_residual, mismatch.partner_flux_residual,
    mismatch.requested_flux_mismatch, mismatch.local_signed_j,
    mismatch.local_abs_j, mismatch.partner_signed_j, mismatch.partner_abs_j,
    mismatch.local_abs_signed_consistent,
    mismatch.partner_abs_signed_consistent);
  fprintf(fp, ",%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d,not_recorded,not_recorded,not_recorded,not_recorded\n",
    mismatch.epsi_rz, mismatch.epsi_phi_raw, mismatch.etheta_direction,
    mismatch.etheta_magnitude, mismatch.etheta_phi_raw,
    mismatch.jacobian_magnitude, smoothness_dq, smoothness_d2q,
    mismatch.objective_sample_valid, pointwise_constraints_valid);
}


static void
gyrokinetic_multib_app_write_interface_pointwise_diag(
  gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-interface_pointwise_geometry_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-interface_pointwise_geometry_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  if (!fp) {
    int saved_errno = errno;
    fprintf(stderr, "Unable to open pointwise interface geometry diagnostic '%s': %s\n",
      file_name.str, strerror(saved_errno));
    cstr_drop(&file_name);
    return;
  }
  cstr_drop(&file_name);
  pointwise_write_geometry_header(fp);

  cstr objective_file_name = rank == 0
    ? cstr_from_fmt("%s-xpt_seam_delta_s_objective.csv", app->name)
    : cstr_from_fmt("%s-xpt_seam_delta_s_objective_rank%d.csv", app->name, rank);
  FILE *objective_fp = fopen(objective_file_name.str, "w");
  if (!objective_fp) {
    int saved_errno = errno;
    fprintf(stderr, "Unable to open X-point seam objective diagnostic '%s': %s\n",
      objective_file_name.str, strerror(saved_errno));
  }
  else {
    pointwise_write_xpt_objective_header(objective_fp);
  }
  cstr_drop(&objective_file_name);

  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  if (cdim < 2 || cdim > 3) {
    fprintf(stderr, "Pointwise interface geometry diagnostic supports cdim=2 or 3; got cdim=%d\n", cdim);
    if (objective_fp)
      fclose(objective_fp);
    fclose(fp);
    return;
  }

  struct gkyl_efit *efit[app->num_local_blocks];
  struct gkyl_efit_inp efit_inp[app->num_local_blocks];
  bool use_cubics[app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    efit[b] = 0;
    use_cubics[b] = false;
    if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      efit_inp[b] = bgi->geometry.efit_info;
      efit_inp[b].use_gpu = false;
      efit[b] = gkyl_efit_new(&efit_inp[b]);
      use_cubics[b] = bgi->geometry.tok_grid_info.use_cubics;
    }
  }

  for (int d=0; d<cdim; ++d) {
    int poly_order = app->singleb_apps[0]->basis.poly_order;
    int num_quad = poly_order+1;
    int nodes_per_face = pointwise_int_pow(num_quad, cdim-1);
    int packed_ncomp = 2*nodes_per_face*POINTWISE_NUM_FIELDS;
    struct gkyl_array *side[app->num_local_blocks];
    struct gkyl_array *side_host[app->num_local_blocks];

    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      side[b] = mkarr(app->use_gpu, packed_ncomp, sbapp->local_ext.volume);
      side_host[b] = app->use_gpu
        ? mkarr(false, packed_ncomp, sbapp->local_ext.volume) : side[b];
      gkyl_array_clear(side_host[b], 0.0);

      if (!efit[b]) {
        if (app->use_gpu)
          gkyl_array_copy(side[b], side_host[b]);
        continue;
      }

      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *mc2p_nodal_fd = pointwise_host_copy(geo_surf->mc2p_nodal_fd);
      struct gkyl_array *dxdz_nodal = pointwise_host_copy(geo_surf->dxdz_nodal);
      struct gkyl_array *jacobgeo_signed_nodal =
        pointwise_host_copy(geo_surf->jacobgeo_signed_nodal);
      struct gkyl_array *jacobgeo_nodal = pointwise_host_copy(geo_surf->jacobgeo_nodal);

      for (int e=0; e<2; ++e) {
        bool on_block_edge = e == 0
          ? sbapp->local.lower[d] == sbapp->global.lower[d]
          : sbapp->local.upper[d] == sbapp->global.upper[d];
        if (!on_block_edge)
          continue;

        const struct gkyl_range *skin_range = e == 0
          ? &sbapp->local_lower_skin[d] : &sbapp->local_upper_skin[d];
        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, skin_range);
        while (gkyl_range_iter_next(&iter)) {
          long cell_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
          double *packed = gkyl_array_fetch(side_host[b], cell_loc);

          for (int n=0; n<nodes_per_face; ++n) {
            int remainder = n;
            int surf_idx[GKYL_MAX_DIM] = { 0 };
            for (int td=0; td<cdim; ++td) {
              if (td == d) {
                surf_idx[td] = e == 0
                  ? sbapp->gk_geom->nrange_surf[d].lower[td]
                  : sbapp->gk_geom->nrange_surf[d].upper[td];
              }
              else {
                int q = remainder%num_quad;
                remainder /= num_quad;
                surf_idx[td] = sbapp->gk_geom->nrange_surf[d].lower[td]
                  +(iter.idx[td]-sbapp->local.lower[td])*num_quad+q;
              }
            }
            pointwise_fill_geometry_record(sbapp, efit[b], use_cubics[b], d,
              surf_idx, mc2p_nodal_fd, dxdz_nodal,
              jacobgeo_signed_nodal, jacobgeo_nodal,
              &packed[(e*nodes_per_face+n)*POINTWISE_NUM_FIELDS]);
          }
        }
      }

      gkyl_array_release(mc2p_nodal_fd);
      gkyl_array_release(dxdz_nodal);
      gkyl_array_release(jacobgeo_signed_nodal);
      gkyl_array_release(jacobgeo_nodal);
      if (app->use_gpu)
        gkyl_array_copy(side[b], side_host[b]);
    }

    int transfer_status = gkyl_multib_comm_conn_array_transfer(app->comm,
      app->num_local_blocks,
      app->local_blocks, app->mbcc_sync_conf->send, app->mbcc_sync_conf->recv,
      side, side);
    if (transfer_status != 0)
      fprintf(stderr, "Pointwise interface geometry transfer failed in direction %d with status %d\n",
        d, transfer_status);
    if (app->use_gpu)
      for (int b=0; b<app->num_local_blocks; ++b)
        gkyl_array_copy(side_host[b], side[b]);

    for (int b=0; b<app->num_local_blocks; ++b) {
      if (!efit[b])
        continue;

      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int bid = app->local_blocks[b];
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *partner =
          &app->block_topo->conn[bid].connections[d][e];
        if (partner->edge == GKYL_PHYSICAL)
          continue;

        bool on_block_edge = e == 0
          ? sbapp->local.lower[d] == sbapp->global.lower[d]
          : sbapp->local.upper[d] == sbapp->global.upper[d];
        if (!on_block_edge)
          continue;

        int partner_edge = oriented_edge_index(partner->edge);
        bool same_direction = partner->dir == d;
        bool complementary_edge = partner_edge == 1-e;
        int orient_sign = orientation_sign(partner->edge);
        const struct gkyl_range *partner_global =
          &app->decomp[partner->bid]->parent_range;
        bool conforming = same_direction;
        bool full_tangential_range = true;
        long interface_cell_count = 1, interface_node_count = 1;
        for (int td=0; td<cdim; ++td) {
          if (td == d)
            continue;
          int num_cells = gkyl_range_shape(&sbapp->global, td);
          conforming = conforming
            && num_cells == gkyl_range_shape(partner_global, td);
          full_tangential_range = full_tangential_range
            && gkyl_range_shape(&sbapp->local, td) == num_cells;
          interface_cell_count *= num_cells;
          interface_node_count *= num_cells*(sbapp->basis.poly_order+1);
        }

        bool can_canonicalize = sbapp->basis.poly_order == 1
          && same_direction && partner_edge >= 0 && orient_sign != 0
          && complementary_edge && conforming
          && (orient_sign > 0 || (cdim == 2 && full_tangential_range));
        const char *comparison_status = "canonicalized";
        if (sbapp->basis.poly_order != 1)
          comparison_status = "unsupported_poly_order";
        else if (!same_direction)
          comparison_status = "unsupported_cross_direction";
        else if (!complementary_edge)
          comparison_status = "unsupported_noncomplementary_edge";
        else if (!conforming)
          comparison_status = "unsupported_nonconforming_interface";
        else if (orient_sign < 0 && cdim != 2)
          comparison_status = "unsupported_negative_orientation_3d";
        else if (orient_sign < 0 && !full_tangential_range)
          comparison_status = "unsupported_negative_orientation_tangential_decomposition";
        else if (!can_canonicalize)
          comparison_status = "unsupported_orientation";
        const char *orientation_transform = can_canonicalize
          ? (orient_sign < 0 ? "reverse_tangential_coordinate" : "identity")
          : "none";
        const struct gkyl_gk_block_geom_info *local_bgi =
          gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
        const struct gkyl_gk_block_geom_info *partner_bgi =
          gkyl_gk_block_geom_get_block(app->gk_block_geom, partner->bid);
        bool unique_xpt_seam = objective_fp && d == 1 &&
          cdim == 2 &&
          (pointwise_xpt_seam_edge_active(local_bgi, e) ||
            pointwise_xpt_seam_edge_active(partner_bgi, partner_edge)) &&
          pointwise_interface_is_owner(bid, d, e, partner->bid,
            partner->dir, partner_edge);
        const struct gkyl_range *ghost_range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, ghost_range);
        while (gkyl_range_iter_next(&iter)) {
          int skin_idx[GKYL_MAX_DIM] = { 0 };
          gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
          skin_idx[d] += e == 0 ? 1 : -1;
          long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
          const double *local_packed = gkyl_array_cfetch(side_host[b], skin_loc);

          int partner_ghost_idx[GKYL_MAX_DIM] = { 0 };
          gkyl_copy_int_arr(cdim, iter.idx, partner_ghost_idx);
          if (can_canonicalize && orient_sign < 0) {
            int tangent_dir = 1-d;
            int cell_offset = skin_idx[tangent_dir]-sbapp->global.lower[tangent_dir];
            partner_ghost_idx[tangent_dir] =
              sbapp->global.upper[tangent_dir]-cell_offset;
          }
          long partner_ghost_loc =
            gkyl_range_idx(&sbapp->local_ext, partner_ghost_idx);
          const double *partner_packed =
            gkyl_array_cfetch(side_host[b], partner_ghost_loc);

          long interface_cell = 0, cell_stride = 1;
          for (int td=0; td<cdim; ++td) {
            if (td == d)
              continue;
            interface_cell += (skin_idx[td]-sbapp->global.lower[td])*cell_stride;
            cell_stride *= gkyl_range_shape(&sbapp->global, td);
          }

          for (int n=0; n<nodes_per_face; ++n) {
            int remainder = n;
            int surf_idx[GKYL_MAX_DIM] = { -1, -1, -1 };
            int quad_idx[GKYL_MAX_DIM] = { -1, -1, -1 };
            long interface_node = 0, node_stride = 1;
            for (int td=0; td<cdim; ++td) {
              if (td == d) {
                surf_idx[td] = e == 0
                  ? sbapp->gk_geom->nrange_surf[d].lower[td]
                  : sbapp->gk_geom->nrange_surf[d].upper[td];
              }
              else {
                quad_idx[td] = remainder%num_quad;
                remainder /= num_quad;
                surf_idx[td] = sbapp->gk_geom->nrange_surf[d].lower[td]
                  +(skin_idx[td]-sbapp->local.lower[td])*num_quad+quad_idx[td];
                long global_node = (skin_idx[td]-sbapp->global.lower[td])*num_quad
                  +quad_idx[td];
                interface_node += global_node*node_stride;
                node_stride *= gkyl_range_shape(&sbapp->global, td)*num_quad;
              }
            }

            const double *local_record =
              &local_packed[(e*nodes_per_face+n)*POINTWISE_NUM_FIELDS];
            double missing_partner[POINTWISE_NUM_FIELDS] = { 0.0 };
            int partner_node_ordinal = -1;
            long partner_interface_cell = -1, partner_interface_node = -1;
            if (can_canonicalize) {
              partner_node_ordinal = orient_sign < 0 ? nodes_per_face-1-n : n;
              partner_interface_cell = orient_sign < 0
                ? interface_cell_count-1-interface_cell : interface_cell;
              partner_interface_node = orient_sign < 0
                ? interface_node_count-1-interface_node : interface_node;
            }
            const double *partner_native = can_canonicalize
              ? &partner_packed[(partner_edge*nodes_per_face+partner_node_ordinal)
                  *POINTWISE_NUM_FIELDS]
              : missing_partner;
            double partner_canonical[POINTWISE_NUM_FIELDS];
            for (int f=0; f<POINTWISE_NUM_FIELDS; ++f)
              partner_canonical[f] = 0.0;
            if (can_canonicalize)
              pointwise_canonicalize_partner_record(d, orient_sign,
                partner_native, partner_canonical);
            if (unique_xpt_seam)
              pointwise_write_xpt_objective_node(objective_fp, app->name,
                rank, bid, e, partner->bid, partner->dir, partner->edge,
                comparison_status, transfer_status, can_canonicalize,
                interface_node, interface_node_count, local_bgi,
                partner_bgi, local_record, partner_canonical);

            fprintf(fp, "%s,%d,%d,%d,%s,%d,%d,%s,%s,%s,%s,%ld,%d,%ld,%ld,%d,%ld,%d,%d,%d,%d,%d,%d,%d,%d,%d",
              app->name, rank, bid, d, edge_name_from_index(e), partner->bid,
              partner->dir, oriented_edge_name(partner->edge),
              orientation_name(partner->edge), comparison_status,
              orientation_transform, interface_cell, n, interface_node,
              partner_interface_cell, partner_node_ordinal,
              partner_interface_node,
              cdim > 0 ? skin_idx[0] : -1,
              cdim > 1 ? skin_idx[1] : -1,
              cdim > 2 ? skin_idx[2] : -1,
              surf_idx[0], surf_idx[1], surf_idx[2],
              quad_idx[0], quad_idx[1], quad_idx[2]);
            pointwise_fprint_geometry_record(fp, local_record);
            pointwise_fprint_geometry_record(fp, partner_native);
            pointwise_fprint_geometry_record(fp, partner_canonical);
            fputc('\n', fp);
          }
        }
      }
    }

    for (int b=0; b<app->num_local_blocks; ++b) {
      if (app->use_gpu)
        gkyl_array_release(side_host[b]);
      gkyl_array_release(side[b]);
    }
  }

  for (int b=0; b<app->num_local_blocks; ++b)
    if (efit[b])
      gkyl_efit_release(efit[b]);
  if (objective_fp)
    fclose(objective_fp);
  fclose(fp);
}

enum xpt_optimizer_reject_reason {
  XPT_OPT_ACCEPTED = 0,
  XPT_OPT_UNSUPPORTED_CONFIGURATION,
  XPT_OPT_INVALID_TOPOLOGY,
  XPT_OPT_INVALID_PARAMETERS,
  XPT_OPT_TRIAL_MAPPING_FAILED,
  XPT_OPT_FIXED_INTERFACE_CHANGED,
  XPT_OPT_RADIAL_ORDERING_FAILED,
  XPT_OPT_JACOBIAN_GUARD_FAILED,
  XPT_OPT_TRANSFER_FAILED,
  XPT_OPT_INCOMPLETE_INTERFACE,
  XPT_OPT_POINTWISE_CONSTRAINT_FAILED,
  XPT_OPT_OBJECTIVE_INVALID,
};

static const char*
xpt_optimizer_reject_name(enum xpt_optimizer_reject_reason reason)
{
  switch (reason) {
    case XPT_OPT_ACCEPTED: return "accepted";
    case XPT_OPT_UNSUPPORTED_CONFIGURATION:
      return "unsupported_configuration";
    case XPT_OPT_INVALID_TOPOLOGY: return "invalid_topology";
    case XPT_OPT_INVALID_PARAMETERS: return "invalid_parameters";
    case XPT_OPT_TRIAL_MAPPING_FAILED: return "trial_mapping_failed";
    case XPT_OPT_FIXED_INTERFACE_CHANGED:
      return "fixed_interface_changed";
    case XPT_OPT_RADIAL_ORDERING_FAILED:
      return "radial_ordering_failed";
    case XPT_OPT_JACOBIAN_GUARD_FAILED:
      return "jacobian_guard_failed";
    case XPT_OPT_TRANSFER_FAILED: return "transfer_failed";
    case XPT_OPT_INCOMPLETE_INTERFACE: return "incomplete_interface";
    case XPT_OPT_POINTWISE_CONSTRAINT_FAILED:
      return "pointwise_constraint_failed";
    case XPT_OPT_OBJECTIVE_INVALID: return "objective_invalid";
  }
  return "unknown";
}

struct xpt_optimizer_pair {
  int bid[2], edge[2];
  enum gkyl_oriented_edge connection;
  double bound;
  bool valid;
  enum xpt_optimizer_reject_reason reject_reason;
};

struct xpt_optimizer_candidate {
  int candidate_index, refinement_round;
  double coefficient, bound;
  bool valid, selected, at_bound;
  enum xpt_optimizer_reject_reason reject_reason;
  int trial_failure_reason;
  long sample_count, expected_sample_count;
  bool zero_map_equivalent, fixed_interfaces_valid;
  bool radial_ordering_valid, jacobian_valid;
  bool pointwise_constraints_valid;
  double max_realized_displacement, min_cell_jacobian_margin;
  double min_jacobian_ratio, max_jacobian_ratio;
  double max_fixed_position_difference, max_anchor_position_difference;
  double epsi_rz, etheta_direction;
  double etheta_log_magnitude, jacobian_log_magnitude;
  double objective;
};

enum { XPT_OPTIMIZER_MAX_PAIRS = 64 };

// Verify that the partner's own connection record points back to (bid,
// dir=1, edge) with the same orientation.  gkyl_gk_block_geom_check_consistency
// checks this redundant data for the whole mesh during app preflight. The
// X-point seam optimizer also checks the specific interface it will use.
static bool
xpt_optimizer_reciprocal_topology(const struct gkyl_gyrokinetic_multib_app *app,
  int bid, int edge, int partner_bid, int partner_dir,
  enum gkyl_oriented_edge partner_oriented_edge)
{
  int partner_edge = oriented_edge_index(partner_oriented_edge);
  if (partner_dir != 1 || partner_edge < 0)
    return false;
  const struct gkyl_target_edge *back =
    &app->block_topo->conn[partner_bid].connections[1][partner_edge];
  return back->bid == bid && back->dir == 1 &&
    oriented_edge_index(back->edge) == edge &&
    orientation_sign(back->edge) == orientation_sign(partner_oriented_edge);
}

// Build and fully validate the seam pair owned by (bid, edge): every static
// topology/parameter requirement in the mission spec is checked here, with
// pair.reject_reason recording the first failure so it can still be
// reported even when the pair is not eligible for search.  Uses only
// globally replicated block-geometry/topology data (never app->local_blocks
// or any per-rank state), so every rank calls this with the same (bid,
// edge) sequence and reaches the identical conclusion.
static struct xpt_optimizer_pair
xpt_optimizer_make_pair(const struct gkyl_gyrokinetic_multib_app *app,
  int bid, int edge)
{
  struct xpt_optimizer_pair pair = {
    .bid = { bid, -1 }, .edge = { edge, -1 }, .connection = GKYL_PHYSICAL,
    .bound = 0.0, .valid = false,
    .reject_reason = XPT_OPT_INVALID_TOPOLOGY,
  };
  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
  const struct gkyl_target_edge *target =
    &app->block_topo->conn[bid].connections[1][edge];
  if (target->edge == GKYL_PHYSICAL)
    return pair;

  int partner_bid = target->bid, partner_dir = target->dir;
  int partner_edge = oriented_edge_index(target->edge);
  int orient_sign = orientation_sign(target->edge);
  int num_blocks = gkyl_gk_block_geom_num_blocks(app->gk_block_geom);
  if (partner_bid < 0 || partner_bid >= num_blocks || partner_edge < 0 ||
      orient_sign == 0)
    return pair;

  bool same_direction = partner_dir == 1;
  bool complementary_edge = partner_edge == 1-edge;
  bool reciprocal = xpt_optimizer_reciprocal_topology(app, bid, edge,
    partner_bid, partner_dir, target->edge);
  if (!same_direction || !complementary_edge || !reciprocal)
    return pair;

  pair.bid[1] = partner_bid;
  pair.edge[1] = partner_edge;
  pair.connection = target->edge;

  const struct gkyl_gk_block_geom_info *partner_bgi =
    gkyl_gk_block_geom_get_block(app->gk_block_geom, partner_bid);
  if (!pointwise_xpt_seam_edge_active(bgi, edge) ||
      !pointwise_xpt_seam_edge_active(partner_bgi, partner_edge))
    return pair;
  if (!pointwise_xpt_ftype_pair(bgi->geometry.tok_grid_info.ftype,
      partner_bgi->geometry.tok_grid_info.ftype))
    return pair;

  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  bool conforming = true;
  for (int td=0; td<cdim; ++td) {
    if (td == 1)
      continue;
    conforming = conforming && bgi->cells[td] == partner_bgi->cells[td];
  }
  // Negative orientation reverses the radial index of the partner cell
  // in-place (see xpt_optimizer_measure_pair); that reversed cell is only
  // guaranteed to be locally addressable when neither block is radially
  // decomposed across ranks, so require cuts[0]==1 on both sides rather
  // than attempt cross-rank addressing here. Theta decomposition
  // (cuts[1] != 1) is excluded too: xpt_optimizer_measure_pair packs each
  // rank's own local skin range directly into a full-block-sized side[]
  // array without an intra-block sync first, so when the seam-adjacent
  // skin cells and this rank's ghost cells are owned by different ranks
  // the samples are silently incomplete -- confirmed by a two-rank
  // regression run with cuts=(1,2), which produced
  // incomplete_interface/fixed_interface_changed on every candidate
  // (including the always-should-succeed zero-coefficient one) even though
  // the same mesh with cuts=(1,1) and blocks split one-per-rank applies
  // correctly. See "regression coverage" commit for the reproduction.
  bool supported_decomp = bgi->cuts[0] == 1 && partner_bgi->cuts[0] == 1 &&
    bgi->cuts[1] == 1 && partner_bgi->cuts[1] == 1;
  if (cdim != 2 || !conforming || !supported_decomp) {
    pair.reject_reason = XPT_OPT_UNSUPPORTED_CONFIGURATION;
    return pair;
  }

  const struct gkyl_tok_geo_grid_inp *local_inp = &bgi->geometry.tok_grid_info;
  const struct gkyl_tok_geo_grid_inp *partner_inp =
    &partner_bgi->geometry.tok_grid_info;
  bool both_request_optimize = local_inp->relaxed_xpt_seam_optimize &&
    partner_inp->relaxed_xpt_seam_optimize;
  // half_domain is deliberately not required here: the relaxed_xpt_seam
  // displacement is wired into both the half-domain and full-domain
  // (extended) construction paths in tok_geo.c (see
  // tok_relax_xpt_seam_point), so this workflow now applies regardless of
  // which domain-splitting mode a block uses.
  bool both_present =
    local_inp->straight_xpt_ray && partner_inp->straight_xpt_ray &&
    local_inp->relaxed_xpt_seam && partner_inp->relaxed_xpt_seam;
  if (!both_request_optimize || !both_present) {
    pair.reject_reason = XPT_OPT_UNSUPPORTED_CONFIGURATION;
    return pair;
  }

  double bound_scale = fmax(1.0,
    fmax(fabs(local_inp->relaxed_xpt_seam_delta_s_bound),
      fabs(partner_inp->relaxed_xpt_seam_delta_s_bound)));
  bool bound_matches =
    isfinite(local_inp->relaxed_xpt_seam_delta_s_bound) &&
    isfinite(partner_inp->relaxed_xpt_seam_delta_s_bound) &&
    local_inp->relaxed_xpt_seam_delta_s_bound > 0.0 &&
    fabs(local_inp->relaxed_xpt_seam_delta_s_bound
      -partner_inp->relaxed_xpt_seam_delta_s_bound)
      <= 64.0*DBL_EPSILON*bound_scale;
  bool baseline_zero =
    local_inp->relaxed_xpt_seam_delta_s_coeff == 0.0 &&
    partner_inp->relaxed_xpt_seam_delta_s_coeff == 0.0 &&
    !local_inp->relaxed_xpt_seam_sweep && !partner_inp->relaxed_xpt_seam_sweep;
  if (!bound_matches || !baseline_zero) {
    pair.reject_reason = XPT_OPT_INVALID_PARAMETERS;
    return pair;
  }

  pair.bound = local_inp->relaxed_xpt_seam_delta_s_bound;
  pair.valid = true;
  pair.reject_reason = XPT_OPT_ACCEPTED;
  return pair;
}

// Discover every X-point seam pair in the (globally replicated) block
// topology, in ascending (bid, edge) order.  Every rank calls this with the
// same gk_block_geom/block_topo regardless of its own local block
// assignment, so every rank produces the identical ordered pair list --
// required because candidate evaluation below is a collective operation
// over app->comm.
static int
xpt_optimizer_discover_pairs(const struct gkyl_gyrokinetic_multib_app *app,
  struct xpt_optimizer_pair pairs[XPT_OPTIMIZER_MAX_PAIRS])
{
  int num_blocks = gkyl_gk_block_geom_num_blocks(app->gk_block_geom);
  int num_pairs = 0;
  bool claimed[num_blocks];
  for (int b=0; b<num_blocks; ++b)
    claimed[b] = false;

  for (int bid=0; bid<num_blocks; ++bid) {
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    // Most ftypes have exactly one active seam edge; a few (CORE_L/CORE_R,
    // the DN SOL "MID" blocks, and self-periodic CORE) are candidates at
    // both edges (see pointwise_xpt_seam_edge_active). Try edge 0 then
    // edge 1: the "claimed" bookkeeping below means at most one actually
    // becomes a pair, and ascending order makes that outcome deterministic
    // rather than a race.
    for (int seam_edge=0; seam_edge<2; ++seam_edge) {
      if (!pointwise_xpt_seam_edge_active(bgi, seam_edge))
        continue;
      const struct gkyl_target_edge *target =
        &app->block_topo->conn[bid].connections[1][seam_edge];
      if (target->edge == GKYL_PHYSICAL)
        continue;
      int partner_edge = oriented_edge_index(target->edge);
      if (!pointwise_interface_is_owner(bid, 1, seam_edge, target->bid,
          target->dir, partner_edge))
        continue;
      if (num_pairs >= XPT_OPTIMIZER_MAX_PAIRS) {
        fprintf(stderr, "X-point seam optimizer: too many candidate pairs "
          "(max %d); skipping remainder\n", XPT_OPTIMIZER_MAX_PAIRS);
        break;
      }

      struct xpt_optimizer_pair pair =
        xpt_optimizer_make_pair(app, bid, seam_edge);
      // One optimizer pair per block: a block already claimed by an earlier
      // pair (either a distinct block claimed on a prior iteration, or this
      // same block's other seam edge claiming it moments ago) makes this
      // pair invalid rather than silently reused.
      bool already_claimed = claimed[bid] ||
        (pair.bid[1] >= 0 && claimed[pair.bid[1]]);
      if (already_claimed) {
        pair.valid = false;
        pair.reject_reason = XPT_OPT_INVALID_TOPOLOGY;
      }
      else {
        claimed[bid] = true;
        if (pair.bid[1] >= 0)
          claimed[pair.bid[1]] = true;
      }
      pairs[num_pairs++] = pair;
    }
  }
  return num_pairs;
}

static struct gkyl_tok_geo_xpt_seam_trial_status
xpt_optimizer_trial_status_init(void)
{
  return (struct gkyl_tok_geo_xpt_seam_trial_status) {
    .contour_valid = true,
    .branch_valid = true,
    .trace_ordering_valid = true,
    .finite_map_valid = true,
    .cell_jacobian_valid = true,
    .jacobian_valid = true,
    .jacobian_sign = 0,
    .first_failure_reason = GKYL_XPT_SEAM_TRIAL_OK,
    .max_realized_displacement = 0.0,
    .min_cell_jacobian_margin = DBL_MAX,
  };
}

static int
xpt_optimizer_local_block_index(const gkyl_gyrokinetic_multib_app *app,
  int bid)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    if (app->local_blocks[b] == bid)
      return b;
  return -1;
}

static bool
xpt_optimizer_array_equal(const struct gkyl_array *left,
  const struct gkyl_array *right)
{
  if (left->size != right->size || left->ncomp != right->ncomp)
    return false;
  for (long i=0; i<left->size; ++i) {
    const double *lv = gkyl_array_cfetch(left, i);
    const double *rv = gkyl_array_cfetch(right, i);
    if (memcmp(lv, rv, left->ncomp*sizeof(double)) != 0)
      return false;
  }
  return true;
}

static bool
xpt_optimizer_same_position(const double *baseline, const double *trial,
  double *max_difference)
{
  double difference = hypot(trial[0]-baseline[0], trial[1]-baseline[1]);
  double scale = fmax(1.0,
    fmax(hypot(baseline[0], baseline[1]), hypot(trial[0], trial[1])));
  *max_difference = fmax(*max_difference, difference);
  return isfinite(difference) && difference <= 64.0*DBL_EPSILON*scale;
}

static bool
xpt_optimizer_array_finite(const struct gkyl_array *array)
{
  for (long i=0; i<array->size; ++i) {
    const double *values = gkyl_array_cfetch(array, i);
    for (int c=0; c<array->ncomp; ++c)
      if (!jacobgeo_ratio_diag_coeff_is_set(values[c]))
        return false;
  }
  return true;
}

// Reject a candidate whose Jacobian magnitude drifts by more than this
// factor, at any single node, relative to the corresponding baseline node.
// This is a pointwise ratio, not a comparison of global minima: a small
// bound displacement should only ever perturb the local Jacobian modestly,
// so a factor of 10 in either direction is already generous headroom while
// still catching an incipient local fold long before jacobian_valid would
// otherwise be tripped by a sign change or nonfinite value.
static const double XPT_OPTIMIZER_MIN_JACOBIAN_RATIO = 0.1;
static const double XPT_OPTIMIZER_MAX_JACOBIAN_RATIO = 10.0;

struct xpt_optimizer_local_guard {
  bool zero_map_equivalent;
  bool fixed_interfaces_valid;
  bool radial_ordering_valid;
  bool jacobian_valid;
  int baseline_jacobian_sign, trial_jacobian_sign;
  double baseline_min_jacobian, trial_min_jacobian;
  double min_jacobian_ratio, max_jacobian_ratio;
  double max_fixed_position_difference;
  double max_anchor_position_difference;
};

static double
xpt_optimizer_eval_corner_component(const struct gkyl_gyrokinetic_app *app,
  const int *cell_idx, const double *eta, int component)
{
  long loc = gkyl_range_idx(&app->local_ext, cell_idx);
  const double *coeff = gkyl_array_cfetch(app->gk_geom->geo_corn.mc2p, loc);
  return app->basis.eval_expand(eta,
    &coeff[component*app->basis.num_basis]);
}

static struct xpt_optimizer_local_guard
xpt_optimizer_check_local_block(
  const struct gkyl_gyrokinetic_app *baseline,
  const struct gkyl_gyrokinetic_app *trial, int seam_edge,
  double coefficient)
{
  struct xpt_optimizer_local_guard guard = {
    .zero_map_equivalent = true,
    .fixed_interfaces_valid = true,
    .radial_ordering_valid = true,
    .jacobian_valid = true,
    .baseline_min_jacobian = DBL_MAX,
    .trial_min_jacobian = DBL_MAX,
    .min_jacobian_ratio = DBL_MAX,
    .max_jacobian_ratio = 0.0,
  };
  if (coefficient == 0.0) {
    guard.zero_map_equivalent =
      xpt_optimizer_array_equal(baseline->gk_geom->geo_corn.mc2p,
        trial->gk_geom->geo_corn.mc2p) &&
      xpt_optimizer_array_equal(baseline->gk_geom->geo_int.mc2p,
        trial->gk_geom->geo_int.mc2p);
    for (int d=0; d<baseline->cdim; ++d)
      guard.zero_map_equivalent = guard.zero_map_equivalent &&
        xpt_optimizer_array_equal(
          baseline->gk_geom->geo_surf[d].mc2p_nodal_fd,
          trial->gk_geom->geo_surf[d].mc2p_nodal_fd);
  }

  // Every radial face and the theta edge opposite the X-point seam are
  // fixed.  Compare central physical R-Z positions at all their surface
  // quadrature nodes.
  for (int d=0; d<baseline->cdim; ++d) {
    const struct gkyl_range *range = &baseline->gk_geom->nrange_surf[d];
    for (int edge=0; edge<2; ++edge) {
      if (d == 1 && edge == seam_edge)
        continue;
      int edge_idx = edge == 0 ? range->lower[d] : range->upper[d];
      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, range);
      while (gkyl_range_iter_next(&iter)) {
        if (iter.idx[d] != edge_idx)
          continue;
        long loc = gkyl_range_idx(range, iter.idx);
        const double *base_position = gkyl_array_cfetch(
          baseline->gk_geom->geo_surf[d].mc2p_nodal_fd, loc);
        const double *trial_position = gkyl_array_cfetch(
          trial->gk_geom->geo_surf[d].mc2p_nodal_fd, loc);
        guard.fixed_interfaces_valid = guard.fixed_interfaces_valid &&
          xpt_optimizer_same_position(base_position, trial_position,
            &guard.max_fixed_position_difference);
      }
    }
  }

  // The two radial endpoints of the seam are cell vertices, not surface
  // quadrature nodes.  Evaluate the p1 corner expansion there explicitly.
  struct gkyl_range_iter cell_iter;
  gkyl_range_iter_init(&cell_iter, &baseline->local);
  while (gkyl_range_iter_next(&cell_iter)) {
    bool on_theta_cell = seam_edge == 0
      ? cell_iter.idx[1] == baseline->global.lower[1]
      : cell_iter.idx[1] == baseline->global.upper[1];
    if (!on_theta_cell)
      continue;
    for (int radial_edge=0; radial_edge<2; ++radial_edge) {
      bool on_radial_cell = radial_edge == 0
        ? cell_iter.idx[0] == baseline->global.lower[0]
        : cell_iter.idx[0] == baseline->global.upper[0];
      if (!on_radial_cell)
        continue;
      double eta[2] = {
        radial_edge == 0 ? -1.0 : 1.0,
        seam_edge == 0 ? -1.0 : 1.0,
      };
      double base_position[2], trial_position[2];
      for (int c=0; c<2; ++c) {
        base_position[c] = xpt_optimizer_eval_corner_component(baseline,
          cell_iter.idx, eta, c);
        trial_position[c] = xpt_optimizer_eval_corner_component(trial,
          cell_iter.idx, eta, c);
      }
      guard.fixed_interfaces_valid = guard.fixed_interfaces_valid &&
        xpt_optimizer_same_position(base_position, trial_position,
          &guard.max_anchor_position_difference);
    }
  }

  // Strict radial progress on the seam, measured against the established
  // straight-map direction at every adjacent p1 surface node.
  const struct gkyl_range *seam_range =
    &baseline->gk_geom->nrange_surf[1];
  int theta_idx = seam_edge == 0
    ? seam_range->lower[1] : seam_range->upper[1];
  for (int ip=seam_range->lower[0]; ip<seam_range->upper[0]; ++ip) {
    int idx0[2] = { ip, theta_idx }, idx1[2] = { ip+1, theta_idx };
    const double *b0 = gkyl_array_cfetch(
      baseline->gk_geom->geo_surf[1].mc2p_nodal_fd,
      gkyl_range_idx(seam_range, idx0));
    const double *b1 = gkyl_array_cfetch(
      baseline->gk_geom->geo_surf[1].mc2p_nodal_fd,
      gkyl_range_idx(seam_range, idx1));
    const double *t0 = gkyl_array_cfetch(
      trial->gk_geom->geo_surf[1].mc2p_nodal_fd,
      gkyl_range_idx(seam_range, idx0));
    const double *t1 = gkyl_array_cfetch(
      trial->gk_geom->geo_surf[1].mc2p_nodal_fd,
      gkyl_range_idx(seam_range, idx1));
    double bdr = b1[0]-b0[0], bdz = b1[1]-b0[1];
    double tdr = t1[0]-t0[0], tdz = t1[1]-t0[1];
    double base_distance = hypot(bdr, bdz);
    double trial_distance = hypot(tdr, tdz);
    double dot = bdr*tdr+bdz*tdz;
    double position_scale = fmax(1.0,
      fmax(hypot(t0[0], t0[1]), hypot(t1[0], t1[1])));
    bool ordered = isfinite(base_distance) && isfinite(trial_distance) &&
      trial_distance > 256.0*DBL_EPSILON*position_scale &&
      dot > 256.0*DBL_EPSILON*base_distance*trial_distance;
    guard.radial_ordering_valid = guard.radial_ordering_valid && ordered;
  }

  // Compare every node against its own baseline counterpart (a pointwise
  // ratio), rather than comparing the global minima of the two maps: a
  // single node collapsing locally would not necessarily move the global
  // minimum, especially when the baseline map already has small Jacobian
  // elsewhere (e.g. near a compressed far surface).
  const struct gkyl_array *base_j =
    baseline->gk_geom->geo_int.jacobgeo_nodal;
  const struct gkyl_array *trial_j = trial->gk_geom->geo_int.jacobgeo_nodal;
  if (base_j->size != trial_j->size)
    guard.jacobian_valid = false;
  else {
    for (long i=0; i<base_j->size; ++i) {
      const double *base_values = gkyl_array_cfetch(base_j, i);
      const double *trial_values = gkyl_array_cfetch(trial_j, i);
      double bj = base_values[0];
      double tj = trial_values[0];
      bool node_valid = isfinite(bj) && bj > 0.0 && isfinite(tj) && tj > 0.0;
      guard.jacobian_valid = guard.jacobian_valid && node_valid;
      guard.baseline_min_jacobian =
        fmin(guard.baseline_min_jacobian, bj);
      guard.trial_min_jacobian = fmin(guard.trial_min_jacobian, tj);
      if (node_valid) {
        double ratio = tj/bj;
        guard.min_jacobian_ratio = fmin(guard.min_jacobian_ratio, ratio);
        guard.max_jacobian_ratio = fmax(guard.max_jacobian_ratio, ratio);
      }
    }
  }
  for (int d=0; d<baseline->cdim; ++d) {
    const struct gkyl_array *base_signed =
      baseline->gk_geom->geo_surf[d].jacobgeo_signed_nodal;
    const struct gkyl_array *trial_signed =
      trial->gk_geom->geo_surf[d].jacobgeo_signed_nodal;
    if (base_signed->size != trial_signed->size) {
      guard.jacobian_valid = false;
      continue;
    }
    for (long i=0; i<base_signed->size; ++i) {
      const double *base_values = gkyl_array_cfetch(base_signed, i);
      const double *trial_values = gkyl_array_cfetch(trial_signed, i);
      double bj = base_values[0];
      double tj = trial_values[0];
      int bsign = bj < 0.0 ? -1 : 1;
      int tsign = tj < 0.0 ? -1 : 1;
      if (guard.baseline_jacobian_sign == 0)
        guard.baseline_jacobian_sign = bsign;
      if (guard.trial_jacobian_sign == 0)
        guard.trial_jacobian_sign = tsign;
      bool node_valid =
        isfinite(bj) && bj != 0.0 && isfinite(tj) && tj != 0.0 &&
        bsign == guard.baseline_jacobian_sign &&
        tsign == guard.trial_jacobian_sign;
      guard.jacobian_valid = guard.jacobian_valid && node_valid;
      guard.baseline_min_jacobian =
        fmin(guard.baseline_min_jacobian, fabs(bj));
      guard.trial_min_jacobian =
        fmin(guard.trial_min_jacobian, fabs(tj));
      if (node_valid) {
        // Same-signed ratio: fabs(J_trial)/fabs(J_baseline) at this node.
        double ratio = tj/bj;
        guard.min_jacobian_ratio = fmin(guard.min_jacobian_ratio, ratio);
        guard.max_jacobian_ratio = fmax(guard.max_jacobian_ratio, ratio);
      }
    }
    guard.jacobian_valid = guard.jacobian_valid &&
      xpt_optimizer_array_finite(
        trial->gk_geom->geo_surf[d].dxdz_nodal);
  }
  guard.jacobian_valid = guard.jacobian_valid &&
    guard.baseline_jacobian_sign != 0 &&
    guard.trial_jacobian_sign == guard.baseline_jacobian_sign &&
    isfinite(guard.min_jacobian_ratio) && isfinite(guard.max_jacobian_ratio) &&
    guard.min_jacobian_ratio >= XPT_OPTIMIZER_MIN_JACOBIAN_RATIO &&
    guard.max_jacobian_ratio <= XPT_OPTIMIZER_MAX_JACOBIAN_RATIO &&
    xpt_optimizer_array_finite(trial->gk_geom->geo_int.dxdz_nodal);
  return guard;
}

static void
xpt_optimizer_measure_pair(gkyl_gyrokinetic_multib_app *app,
  struct gkyl_gyrokinetic_app **trial_apps,
  const struct xpt_optimizer_pair *pair, struct gkyl_efit **efit,
  struct xpt_optimizer_candidate *candidate)
{
  enum { NUM_QUAD = 2, NODES_PER_FACE = 2 };
  int packed_ncomp = 2*NODES_PER_FACE*POINTWISE_NUM_FIELDS;
  struct gkyl_array *side[app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = trial_apps[b];
    side[b] = mkarr(false, packed_ncomp, sbapp->local_ext.volume);
    gkyl_array_clear(side[b], 0.0);
    int bid = app->local_blocks[b];
    if ((bid != pair->bid[0] && bid != pair->bid[1]) || !efit[b])
      continue;
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    const struct gk_geom_surf *surf = &sbapp->gk_geom->geo_surf[1];
    for (int edge=0; edge<2; ++edge) {
      bool on_block_edge = edge == 0
        ? sbapp->local.lower[1] == sbapp->global.lower[1]
        : sbapp->local.upper[1] == sbapp->global.upper[1];
      if (!on_block_edge)
        continue;
      const struct gkyl_range *skin = edge == 0
        ? &sbapp->local_lower_skin[1] : &sbapp->local_upper_skin[1];
      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, skin);
      while (gkyl_range_iter_next(&iter)) {
        long cell_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
        double *packed = gkyl_array_fetch(side[b], cell_loc);
        for (int n=0; n<NODES_PER_FACE; ++n) {
          int surf_idx[2] = {
            sbapp->gk_geom->nrange_surf[1].lower[0]
              +(iter.idx[0]-sbapp->local.lower[0])*NUM_QUAD+n,
            edge == 0 ? sbapp->gk_geom->nrange_surf[1].lower[1]
              : sbapp->gk_geom->nrange_surf[1].upper[1],
          };
          pointwise_fill_geometry_record(sbapp, efit[b],
            bgi->geometry.tok_grid_info.use_cubics, 1, surf_idx,
            surf->mc2p_nodal_fd, surf->dxdz_nodal,
            surf->jacobgeo_signed_nodal, surf->jacobgeo_nodal,
            &packed[(edge*NODES_PER_FACE+n)*POINTWISE_NUM_FIELDS]);
        }
      }
    }
  }

  int transfer_status = gkyl_multib_comm_conn_array_transfer(app->comm,
    app->num_local_blocks, app->local_blocks, app->mbcc_sync_conf->send,
    app->mbcc_sync_conf->recv, side, side);
  int64_t transfer_failed_local = transfer_status != 0;
  int64_t transfer_failed = 0;
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, 1,
    &transfer_failed_local, &transfer_failed);

  double local_sum[4] = { 0.0 }, global_sum[4] = { 0.0 };
  int64_t local_count = 0, global_count = 0;
  int64_t local_invalid = 0, global_invalid = 0;
  int owner_local = xpt_optimizer_local_block_index(app, pair->bid[0]);
  if (owner_local >= 0) {
    struct gkyl_gyrokinetic_app *sbapp = trial_apps[owner_local];
    int edge = pair->edge[0], partner_edge = pair->edge[1];
    int orient_sign = orientation_sign(pair->connection);
    const struct gkyl_range *ghost = edge == 0
      ? &sbapp->local_lower_ghost[1] : &sbapp->local_upper_ghost[1];
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, ghost);
    while (gkyl_range_iter_next(&iter)) {
      int skin_idx[2] = { iter.idx[0], iter.idx[1] };
      skin_idx[1] += edge == 0 ? 1 : -1;
      const double *local_packed = gkyl_array_cfetch(side[owner_local],
        gkyl_range_idx(&sbapp->local_ext, skin_idx));
      int partner_ghost_idx[2] = { iter.idx[0], iter.idx[1] };
      if (orient_sign < 0) {
        int cell_offset = skin_idx[0]-sbapp->global.lower[0];
        partner_ghost_idx[0] = sbapp->global.upper[0]-cell_offset;
      }
      const double *partner_packed = gkyl_array_cfetch(side[owner_local],
        gkyl_range_idx(&sbapp->local_ext, partner_ghost_idx));
      for (int n=0; n<NODES_PER_FACE; ++n) {
        int partner_node = orient_sign < 0 ? NODES_PER_FACE-1-n : n;
        const double *local_record =
          &local_packed[(edge*NODES_PER_FACE+n)*POINTWISE_NUM_FIELDS];
        const double *partner_native =
          &partner_packed[(partner_edge*NODES_PER_FACE+partner_node)
            *POINTWISE_NUM_FIELDS];
        double partner_record[POINTWISE_NUM_FIELDS] = { 0.0 };
        pointwise_canonicalize_partner_record(1, orient_sign,
          partner_native, partner_record);
        struct pointwise_xpt_mismatch mismatch =
          pointwise_measure_xpt_mismatch(local_record, partner_record);
        local_count += 1;
        bool valid = mismatch.objective_sample_valid &&
          mismatch.geometry_constraints_valid;
        local_invalid = local_invalid || !valid;
        if (mismatch.objective_sample_valid) {
          local_sum[0] += mismatch.epsi_rz;
          local_sum[1] += mismatch.etheta_direction;
          local_sum[2] += mismatch.etheta_magnitude;
          local_sum[3] += mismatch.jacobian_magnitude;
        }
      }
    }
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 4,
    local_sum, global_sum);
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_SUM, 1,
    &local_count, &global_count);
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, 1,
    &local_invalid, &global_invalid);
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_array_release(side[b]);

  const struct gkyl_gk_block_geom_info *owner_bgi =
    gkyl_gk_block_geom_get_block(app->gk_block_geom, pair->bid[0]);
  candidate->expected_sample_count =
    2L*owner_bgi->cells[0];
  candidate->sample_count = global_count;
  candidate->pointwise_constraints_valid =
    !transfer_failed && !global_invalid &&
    global_count == candidate->expected_sample_count;
  if (global_count > 0) {
    candidate->epsi_rz = global_sum[0]/global_count;
    candidate->etheta_direction = global_sum[1]/global_count;
    candidate->etheta_log_magnitude = global_sum[2]/global_count;
    candidate->jacobian_log_magnitude = global_sum[3]/global_count;
    candidate->objective = candidate->epsi_rz
      +candidate->etheta_direction
      +0.5*(candidate->etheta_log_magnitude
        +candidate->jacobian_log_magnitude);
  }
  if (transfer_failed)
    candidate->reject_reason = XPT_OPT_TRANSFER_FAILED;
  else if (global_count != candidate->expected_sample_count)
    candidate->reject_reason = XPT_OPT_INCOMPLETE_INTERFACE;
  else if (global_invalid)
    candidate->reject_reason = XPT_OPT_POINTWISE_CONSTRAINT_FAILED;
  else if (!isfinite(candidate->objective))
    candidate->reject_reason = XPT_OPT_OBJECTIVE_INVALID;
}

static bool
xpt_optimizer_status_valid(
  const struct gkyl_tok_geo_xpt_seam_trial_status *status)
{
  return status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK &&
    status->contour_valid && status->branch_valid &&
    status->trace_ordering_valid && status->finite_map_valid &&
    status->cell_jacobian_valid && status->jacobian_valid &&
    status->jacobian_sign != 0 &&
    isfinite(status->min_cell_jacobian_margin) &&
    status->min_cell_jacobian_margin > 256.0*DBL_EPSILON &&
    isfinite(status->max_realized_displacement);
}

static struct xpt_optimizer_candidate
xpt_optimizer_evaluate_candidate(
  const struct gkyl_gyrokinetic_multib *mbinp,
  gkyl_gyrokinetic_multib_app *app,
  const struct xpt_optimizer_pair *pair, double coefficient,
  int candidate_index, int refinement_round, struct gkyl_efit **efit)
{
  struct xpt_optimizer_candidate result = {
    .candidate_index = candidate_index,
    .refinement_round = refinement_round,
    .coefficient = coefficient,
    .bound = pair->bound,
    .zero_map_equivalent = true,
    .fixed_interfaces_valid = true,
    .radial_ordering_valid = true,
    .jacobian_valid = true,
    .pointwise_constraints_valid = false,
    .min_cell_jacobian_margin = DBL_MAX,
    .min_jacobian_ratio = DBL_MAX,
    .max_jacobian_ratio = 0.0,
    .objective = DBL_MAX,
    .reject_reason = XPT_OPT_ACCEPTED,
  };
  result.at_bound = fabs(fabs(coefficient)-pair->bound)
    <= 64.0*DBL_EPSILON*fmax(1.0, pair->bound);
  if (!isfinite(coefficient) || fabs(coefficient) > pair->bound
      +64.0*DBL_EPSILON*fmax(1.0, pair->bound)) {
    result.reject_reason = XPT_OPT_INVALID_PARAMETERS;
    return result;
  }

  fprintf(stderr, "TOK_GEO_OPTIMIZER_SCOPE kind=candidate event=begin\n");
  fflush(stderr);
  struct gkyl_gyrokinetic_app **trial_apps =
    gkyl_malloc(sizeof(*trial_apps)*app->num_local_blocks);
  bool *owned = gkyl_calloc(app->num_local_blocks, sizeof(bool));
  struct gkyl_tok_geo_xpt_seam_trial_status status[2] = {
    xpt_optimizer_trial_status_init(), xpt_optimizer_trial_status_init(),
  };
  int64_t local_mapping_failed = 0, local_zero_failed = 0;
  int64_t local_fixed_failed = 0, local_ordering_failed = 0;
  int64_t local_jacobian_failed = 0, local_failure_reason = 0;
  double local_max_displacement = 0.0;
  double local_min_margin = DBL_MAX;
  double local_min_jacobian_ratio = DBL_MAX, local_max_jacobian_ratio = 0.0;
  double local_max_fixed_position_difference = 0.0;
  double local_max_anchor_position_difference = 0.0;

  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    int pair_side = bid == pair->bid[0] ? 0 : bid == pair->bid[1] ? 1 : -1;
    if (pair_side < 0) {
      trial_apps[b] = app->singleb_apps[b];
      continue;
    }
    const struct gkyl_gk_block_geom_info *base_bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    struct gkyl_gk_block_geom_info trial_bgi = *base_bgi;
    struct gkyl_tok_geo_grid_inp *trial_inp =
      &trial_bgi.geometry.tok_grid_info;
    trial_inp->relaxed_xpt_seam = true;
    trial_inp->relaxed_xpt_seam_sweep = true;
    trial_inp->relaxed_xpt_seam_delta_s_coeff = coefficient;
    trial_inp->relaxed_xpt_seam_delta_s_bound = pair->bound;
    trial_inp->relaxed_xpt_seam_optimize = false;
    trial_inp->relaxed_xpt_seam_optimizer_trial = true;
    trial_inp->relaxed_xpt_seam_trial_status = &status[pair_side];
    trial_apps[b] = singleb_app_new_geom_from_block(mbinp, bid, app,
      &trial_bgi, false);
    owned[b] = true;

    bool status_valid = xpt_optimizer_status_valid(&status[pair_side]);
    local_mapping_failed = local_mapping_failed || !status_valid;
    local_failure_reason = GKYL_MAX2(local_failure_reason,
      status[pair_side].first_failure_reason);
    local_max_displacement = fmax(local_max_displacement,
      status[pair_side].max_realized_displacement);
    local_min_margin = fmin(local_min_margin,
      status[pair_side].min_cell_jacobian_margin);

    struct xpt_optimizer_local_guard guard =
      xpt_optimizer_check_local_block(app->singleb_apps[b], trial_apps[b],
        pair->edge[pair_side], coefficient);
    local_zero_failed = local_zero_failed || !guard.zero_map_equivalent;
    local_fixed_failed = local_fixed_failed ||
      !guard.fixed_interfaces_valid;
    local_ordering_failed = local_ordering_failed ||
      !guard.radial_ordering_valid;
    local_jacobian_failed = local_jacobian_failed ||
      !guard.jacobian_valid || (status_valid &&
        status[pair_side].jacobian_sign != guard.baseline_jacobian_sign);
    local_min_jacobian_ratio =
      fmin(local_min_jacobian_ratio, guard.min_jacobian_ratio);
    local_max_jacobian_ratio =
      fmax(local_max_jacobian_ratio, guard.max_jacobian_ratio);
    local_max_fixed_position_difference = fmax(
      local_max_fixed_position_difference, guard.max_fixed_position_difference);
    local_max_anchor_position_difference = fmax(
      local_max_anchor_position_difference, guard.max_anchor_position_difference);
  }

  int64_t local_flags[5] = {
    local_mapping_failed, local_zero_failed, local_fixed_failed,
    local_ordering_failed, local_jacobian_failed,
  };
  int64_t global_flags[5] = { 0 };
  int64_t global_failure_reason = 0;
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, 5,
    local_flags, global_flags);
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, 1,
    &local_failure_reason, &global_failure_reason);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &local_max_displacement, &result.max_realized_displacement);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1,
    &local_min_margin, &result.min_cell_jacobian_margin);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1,
    &local_min_jacobian_ratio, &result.min_jacobian_ratio);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &local_max_jacobian_ratio, &result.max_jacobian_ratio);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &local_max_fixed_position_difference,
    &result.max_fixed_position_difference);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &local_max_anchor_position_difference,
    &result.max_anchor_position_difference);
  result.trial_failure_reason = global_failure_reason;
  result.zero_map_equivalent = !global_flags[1];
  result.fixed_interfaces_valid = !global_flags[2];
  result.radial_ordering_valid = !global_flags[3];
  result.jacobian_valid = !global_flags[4];
  double displacement_tolerance =
    64.0*DBL_EPSILON*fmax(1.0, pair->bound);
  if (result.max_realized_displacement > pair->bound
      +displacement_tolerance)
    global_flags[0] = 1;

  if (global_flags[0])
    result.reject_reason = XPT_OPT_TRIAL_MAPPING_FAILED;
  else if (global_flags[1] || global_flags[2])
    result.reject_reason = XPT_OPT_FIXED_INTERFACE_CHANGED;
  else if (global_flags[3])
    result.reject_reason = XPT_OPT_RADIAL_ORDERING_FAILED;
  else if (global_flags[4])
    result.reject_reason = XPT_OPT_JACOBIAN_GUARD_FAILED;
  else
    xpt_optimizer_measure_pair(app, trial_apps, pair, efit, &result);

  result.valid = result.reject_reason == XPT_OPT_ACCEPTED &&
    result.pointwise_constraints_valid && isfinite(result.objective);
  for (int b=0; b<app->num_local_blocks; ++b)
    if (owned[b])
      gkyl_gyrokinetic_app_release_geom(trial_apps[b]);
  gkyl_free(owned);
  gkyl_free(trial_apps);
  fprintf(stderr, "TOK_GEO_OPTIMIZER_SCOPE kind=candidate event=end accepted=0\n");
  fflush(stderr);
  return result;
}

enum { XPT_OPTIMIZER_MAX_CANDIDATES = 32 };

struct xpt_optimizer_search {
  int count;
  struct xpt_optimizer_candidate candidates[XPT_OPTIMIZER_MAX_CANDIDATES];
  int zero_index;        // index of the coefficient==0 candidate, or -1
  int best_valid_index;  // argmin objective among valid candidates, or -1
  int recommended_index; // best_valid_index if it improves on zero, else -1
  double zero_objective;
};

static bool
xpt_optimizer_already_tried(const struct xpt_optimizer_search *search,
  double coefficient, double bound)
{
  double tol = 64.0*DBL_EPSILON*fmax(1.0, bound);
  for (int i=0; i<search->count; ++i)
    if (fabs(search->candidates[i].coefficient-coefficient) <= tol)
      return true;
  return false;
}

// Deterministic tie-break for the diagnostic search: strictly smaller
// objective wins; ties are broken first by smaller |coefficient| (the more
// conservative displacement), then by numerical coefficient value.
static bool
xpt_optimizer_candidate_better(const struct xpt_optimizer_candidate *a,
  const struct xpt_optimizer_candidate *b)
{
  double tol = 1e-12*fmax(1.0, fmax(fabs(a->objective), fabs(b->objective)));
  if (fabs(a->objective-b->objective) > tol)
    return a->objective < b->objective;
  if (fabs(fabs(a->coefficient)-fabs(b->coefficient)) > 64.0*DBL_EPSILON)
    return fabs(a->coefficient) < fabs(b->coefficient);
  return a->coefficient < b->coefficient;
}

static void
xpt_optimizer_try_candidate(const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *app,
  const struct xpt_optimizer_pair *pair, struct gkyl_efit **efit,
  double coefficient, int refinement_round,
  struct xpt_optimizer_search *search)
{
  if (search->count >= XPT_OPTIMIZER_MAX_CANDIDATES)
    return;
  double clipped = fmax(-pair->bound, fmin(pair->bound, coefficient));
  if (xpt_optimizer_already_tried(search, clipped, pair->bound))
    return;
  int idx = search->count++;
  search->candidates[idx] = xpt_optimizer_evaluate_candidate(mbinp, app, pair,
    clipped, idx, refinement_round, efit);
  if (clipped == 0.0)
    search->zero_index = idx;
  if (search->candidates[idx].valid &&
      (search->best_valid_index < 0 || xpt_optimizer_candidate_better(
        &search->candidates[idx], &search->candidates[search->best_valid_index])))
    search->best_valid_index = idx;
}

// Deterministic bounded search: nine coarse samples spanning the full
// bound (this always includes the zero baseline and both endpoints),
// followed by three bounded local-refinement rounds around the best valid
// coefficient found so far, halving the step each round.  Every candidate
// is recorded, valid or not.
static struct xpt_optimizer_search
xpt_optimizer_run_search(const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *app,
  const struct xpt_optimizer_pair *pair, struct gkyl_efit **efit)
{
  struct xpt_optimizer_search search = {
    .count = 0, .zero_index = -1, .best_valid_index = -1,
    .recommended_index = -1,
  };

  for (int k=-4; k<=4; ++k)
    xpt_optimizer_try_candidate(mbinp, app, pair, efit,
      k*pair->bound/4.0, 0, &search);

  for (int round=1; round<=3; ++round) {
    double center = search.best_valid_index >= 0
      ? search.candidates[search.best_valid_index].coefficient : 0.0;
    double step = (pair->bound/4.0)/(1 << round);
    xpt_optimizer_try_candidate(mbinp, app, pair, efit,
      center-step, round, &search);
    xpt_optimizer_try_candidate(mbinp, app, pair, efit,
      center+step, round, &search);
  }

  if (search.zero_index >= 0) {
    search.zero_objective = search.candidates[search.zero_index].objective;
    if (search.best_valid_index >= 0 &&
        search.best_valid_index != search.zero_index &&
        search.candidates[search.best_valid_index].objective <
          search.zero_objective) {
      search.recommended_index = search.best_valid_index;
      search.candidates[search.recommended_index].selected = true;
    }
  }
  return search;
}

static void
xpt_optimizer_write_header(FILE *fp)
{
  fprintf(fp,
    "app,rank,pair_index,bid0,bid1,edge0,edge1,ftype0,ftype1,orientation,"
    "pair_valid,pair_reject_reason,bound,"
    "candidate_index,refinement_round,coefficient,at_bound,"
    "trial_failure_reason,zero_map_equivalent,fixed_interfaces_valid,"
    "radial_ordering_valid,jacobian_valid,pointwise_constraints_valid,"
    "sample_count,expected_sample_count,"
    "epsi_rz,etheta_direction,etheta_log_magnitude,jacobian_log_magnitude,"
    "objective,max_realized_displacement,min_cell_jacobian_margin,"
    "min_jacobian_ratio,max_jacobian_ratio,"
    "max_fixed_position_difference,max_anchor_position_difference,"
    "candidate_valid,candidate_reject_reason,is_zero_candidate,selected,"
    "relative_improvement_vs_zero\n");
}

static void
xpt_optimizer_write_candidate_row(FILE *fp, const char *app_name, int rank,
  int pair_index, const struct xpt_optimizer_pair *pair,
  const struct gkyl_gk_block_geom_info *bgi0,
  const struct gkyl_gk_block_geom_info *bgi1,
  const struct xpt_optimizer_candidate *c, const char *candidate_status,
  bool is_zero_candidate, double zero_objective)
{
  int ftype0 = bgi0->geometry.tok_grid_info.ftype;
  int ftype1 = bgi1 ? bgi1->geometry.tok_grid_info.ftype : -1;
  double relative_improvement =
    isfinite(zero_objective) && zero_objective > 0.0 && isfinite(c->objective)
      ? (zero_objective-c->objective)/zero_objective : 0.0;

  fprintf(fp,
    "%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%s,%.17e,"
    "%d,%d,%.17e,%d,"
    "%d,%d,%d,%d,%d,%d,"
    "%ld,%ld,"
    "%.17e,%.17e,%.17e,%.17e,"
    "%.17e,%.17e,%.17e,"
    "%.17e,%.17e,"
    "%.17e,%.17e,"
    "%d,%s,%d,%d,%.17e\n",
    app_name, rank, pair_index, pair->bid[0], pair->bid[1],
    pair->edge[0], pair->edge[1], ftype0, ftype1,
    orientation_name(pair->connection), pair->valid,
    xpt_optimizer_reject_name(pair->reject_reason), pair->bound,
    c->candidate_index, c->refinement_round, c->coefficient, c->at_bound,
    c->trial_failure_reason, c->zero_map_equivalent, c->fixed_interfaces_valid,
    c->radial_ordering_valid, c->jacobian_valid, c->pointwise_constraints_valid,
    c->sample_count, c->expected_sample_count,
    c->epsi_rz, c->etheta_direction, c->etheta_log_magnitude,
    c->jacobian_log_magnitude,
    c->objective, c->max_realized_displacement, c->min_cell_jacobian_margin,
    c->min_jacobian_ratio, c->max_jacobian_ratio,
    c->max_fixed_position_difference, c->max_anchor_position_difference,
    c->valid, candidate_status, is_zero_candidate, c->selected,
    relative_improvement);
}

static void
xpt_optimizer_write_pair_not_evaluated(FILE *fp, const char *app_name,
  int rank, int pair_index, const struct xpt_optimizer_pair *pair,
  const struct gkyl_gk_block_geom_info *bgi0,
  const struct gkyl_gk_block_geom_info *bgi1)
{
  struct xpt_optimizer_candidate placeholder = {
    .candidate_index = -1, .refinement_round = -1,
  };
  xpt_optimizer_write_candidate_row(fp, app_name, rank, pair_index, pair,
    bgi0, bgi1, &placeholder, "not_evaluated", false, 0.0);
}

// Minimum relative reduction in the primary mismatch objective (against the
// zero-coefficient baseline) required before a recommended coefficient is
// actually applied to production geometry. Below this, a displacement is
// reported as diagnostic-only and the straight map is kept: a small or
// possibly noise-level "improvement" is not worth deviating from the
// well-understood straight-ray construction.
static const double XPT_APPLY_MIN_RELATIVE_IMPROVEMENT = 0.05;

enum xpt_apply_status {
  XPT_APPLY_NOT_USEFUL = 0,
  XPT_APPLY_NO_VALID_CANDIDATE,
  XPT_APPLY_APPLIED,
  XPT_APPLY_REVALIDATION_FAILED,
};

static const char*
xpt_apply_status_name(enum xpt_apply_status status)
{
  switch (status) {
    case XPT_APPLY_NOT_USEFUL: return "not_useful";
    case XPT_APPLY_NO_VALID_CANDIDATE: return "no_valid_candidate";
    case XPT_APPLY_APPLIED: return "applied";
    case XPT_APPLY_REVALIDATION_FAILED: return "revalidation_failed";
  }
  return "unknown";
}

// Rebuild write-enabled, production geometry for both (locally-owned) blocks
// in `pair` at `coefficient`, re-running the same hard guards used during
// the search on this fresh construction -- this is the last gate before
// touching mbapp->singleb_apps, so a trial that merely looked good during
// the (no-write) search is not trusted blindly. On success, releases the
// previous (straight) singleb_apps for this pair's local blocks and installs
// the new ones; every other block is untouched. On failure, releases the
// newly-built candidates and leaves mbapp->singleb_apps exactly as it was:
// the guarded fallback to the straight map is simply "do nothing".
static enum xpt_apply_status
xpt_optimizer_apply_selection(const struct gkyl_gyrokinetic_multib *mbinp,
  gkyl_gyrokinetic_multib_app *app, const struct xpt_optimizer_pair *pair,
  double coefficient)
{
  fprintf(stderr, "TOK_GEO_OPTIMIZER_SCOPE kind=revalidation event=begin\n");
  fflush(stderr);
  struct gkyl_gyrokinetic_app **candidate_apps =
    gkyl_malloc(sizeof(*candidate_apps)*app->num_local_blocks);
  bool *owned = gkyl_calloc(app->num_local_blocks, sizeof(bool));
  struct gkyl_tok_geo_xpt_seam_trial_status status[2] = {
    xpt_optimizer_trial_status_init(), xpt_optimizer_trial_status_init(),
  };
  int64_t local_mapping_failed = 0, local_fixed_failed = 0;
  int64_t local_ordering_failed = 0, local_jacobian_failed = 0;
  double local_max_displacement = 0.0;

  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    int pair_side = bid == pair->bid[0] ? 0 : bid == pair->bid[1] ? 1 : -1;
    if (pair_side < 0)
      continue;
    const struct gkyl_gk_block_geom_info *base_bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    struct gkyl_gk_block_geom_info production_bgi = *base_bgi;
    struct gkyl_tok_geo_grid_inp *production_inp =
      &production_bgi.geometry.tok_grid_info;
    production_inp->relaxed_xpt_seam = true;
    production_inp->relaxed_xpt_seam_sweep = true;
    production_inp->relaxed_xpt_seam_delta_s_coeff = coefficient;
    production_inp->relaxed_xpt_seam_delta_s_bound = pair->bound;
    production_inp->relaxed_xpt_seam_optimize = false;
    production_inp->relaxed_xpt_seam_optimizer_trial = true;
    production_inp->relaxed_xpt_seam_trial_status = &status[pair_side];
    candidate_apps[b] = singleb_app_new_geom_from_block(mbinp, bid, app,
      &production_bgi, true);
    owned[b] = true;

    bool status_valid = xpt_optimizer_status_valid(&status[pair_side]);
    local_mapping_failed = local_mapping_failed || !status_valid;
    local_max_displacement = fmax(local_max_displacement,
      status[pair_side].max_realized_displacement);

    struct xpt_optimizer_local_guard guard = xpt_optimizer_check_local_block(
      app->singleb_apps[b], candidate_apps[b], pair->edge[pair_side],
      coefficient);
    local_fixed_failed = local_fixed_failed || !guard.fixed_interfaces_valid;
    local_ordering_failed = local_ordering_failed ||
      !guard.radial_ordering_valid;
    local_jacobian_failed = local_jacobian_failed || !guard.jacobian_valid ||
      (status_valid &&
        status[pair_side].jacobian_sign != guard.baseline_jacobian_sign);
  }

  int64_t local_flags[4] = {
    local_mapping_failed, local_fixed_failed, local_ordering_failed,
    local_jacobian_failed,
  };
  int64_t global_flags[4] = { 0 };
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, 4, local_flags,
    global_flags);
  double max_realized_displacement = 0.0;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1,
    &local_max_displacement, &max_realized_displacement);
  double displacement_tolerance = 64.0*DBL_EPSILON*fmax(1.0, pair->bound);
  bool bound_violated =
    max_realized_displacement > pair->bound+displacement_tolerance;
  bool success = !global_flags[0] && !global_flags[1] && !global_flags[2] &&
    !global_flags[3] && !bound_violated;

  for (int b=0; b<app->num_local_blocks; ++b) {
    if (!owned[b])
      continue;
    if (success) {
      gkyl_gyrokinetic_app_release_geom(app->singleb_apps[b]);
      app->singleb_apps[b] = candidate_apps[b];
    }
    else
      gkyl_gyrokinetic_app_release_geom(candidate_apps[b]);
  }
  gkyl_free(owned);
  gkyl_free(candidate_apps);

  if (success) {
    // gk_block_geom is replicated identically on every rank (it is queried
    // for any block ID regardless of locality elsewhere in this file), so
    // every rank -- not just the ones owning bid[0]/bid[1] -- must record
    // the selection to keep that replica consistent; otherwise diagnostics
    // and any future reconstruction of these blocks would keep reporting
    // the pre-selection zero coefficient even though the live geometry has
    // actually been displaced.
    gkyl_gk_block_geom_apply_xpt_seam_selection(app->gk_block_geom,
      pair->bid[0], coefficient, pair->bound);
    gkyl_gk_block_geom_apply_xpt_seam_selection(app->gk_block_geom,
      pair->bid[1], coefficient, pair->bound);
  }
  fprintf(stderr, "TOK_GEO_OPTIMIZER_SCOPE kind=revalidation event=end accepted=%d\n", (int) success);
  fflush(stderr);
  return success ? XPT_APPLY_APPLIED : XPT_APPLY_REVALIDATION_FAILED;
}

static void
xpt_optimizer_write_applied_header(FILE *fp)
{
  fprintf(fp, "app,pair_index,bid0,bid1,edge0,edge1,bound,"
    "useful_threshold,zero_objective,recommended_coefficient,"
    "recommended_objective,relative_improvement,useful,applied,"
    "apply_status\n");
}

static void
xpt_optimizer_write_applied_row(FILE *fp, const char *app_name,
  int pair_index, const struct xpt_optimizer_pair *pair,
  double zero_objective, double coefficient, double objective,
  double relative_improvement, bool useful, bool applied,
  enum xpt_apply_status status)
{
  fprintf(fp, "%s,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d,%s\n",
    app_name, pair_index, pair->bid[0], pair->bid[1], pair->edge[0],
    pair->edge[1], pair->bound, XPT_APPLY_MIN_RELATIVE_IMPROVEMENT,
    zero_objective, coefficient, objective, relative_improvement, useful,
    applied, xpt_apply_status_name(status));
}

// Evaluates a deterministic, bounded set of delta-s candidates per X-point
// seam pair and records the result (including every rejected candidate and
// reason) to a CSV.  A recommendation that clears
// XPT_APPLY_MIN_RELATIVE_IMPROVEMENT is then re-validated against a fresh,
// write-enabled production build and, only if that also passes every hard
// guard, swapped into mbapp->singleb_apps for that pair's blocks; otherwise
// (not useful, no valid candidate, or the production re-validation itself
// fails) the straight map already in mbapp->singleb_apps is left completely
// untouched -- the guarded fallback requires no explicit action.  The
// coefficient is always identical on both blocks sharing a seam, since a
// single search produces one coefficient per pair.
static void
gyrokinetic_multib_optimize_xpt_seams(
  const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *mbapp)
{
  int rank;
  gkyl_comm_get_rank(mbapp->comm, &rank);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);

  // Any block that asks for optimization but hits an unsupported build
  // configuration must fail visibly, never silently.  The host-side array
  // inspection throughout this optimizer is not GPU safe, and only
  // cdim=2, poly_order=1 tokamak geometry is supported so far.
  bool any_requests_optimize = false;
  for (int bid=0; bid<num_blocks; ++bid) {
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);
    if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK &&
        bgi->geometry.tok_grid_info.relaxed_xpt_seam_optimize)
      any_requests_optimize = true;
  }
  if (!any_requests_optimize)
    return;
  if (mbapp->use_gpu) {
    fprintf(stderr, "X-point seam optimizer: relaxed_xpt_seam_optimize was "
      "requested but GPU builds are not yet supported; skipping diagnostic "
      "optimization entirely.\n");
    return;
  }
  if (cdim != 2 || mbinp->poly_order != 1) {
    fprintf(stderr, "X-point seam optimizer: relaxed_xpt_seam_optimize was "
      "requested but only cdim=2, poly_order=1 is supported (got cdim=%d, "
      "poly_order=%d); skipping diagnostic optimization entirely.\n",
      cdim, mbinp->poly_order);
    return;
  }

  struct xpt_optimizer_pair pairs[XPT_OPTIMIZER_MAX_PAIRS];
  int num_pairs = xpt_optimizer_discover_pairs(mbapp, pairs);
  if (num_pairs == 0)
    return;

  struct gkyl_efit *efit[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    int bid = mbapp->local_blocks[b];
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);
    efit[b] = 0;
    if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      struct gkyl_efit_inp efit_inp = bgi->geometry.efit_info;
      efit_inp.use_gpu = false;
      efit[b] = gkyl_efit_new(&efit_inp);
    }
  }

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-xpt_seam_optimizer_diagnostics.csv", mbapp->name)
    : cstr_from_fmt("%s-xpt_seam_optimizer_diagnostics_rank%d.csv",
        mbapp->name, rank);
  FILE *fp = rank == 0 ? fopen(file_name.str, "w") : NULL;
  if (rank == 0 && !fp) {
    int saved_errno = errno;
    fprintf(stderr, "Unable to open X-point seam optimizer diagnostic "
      "'%s': %s\n", file_name.str, strerror(saved_errno));
  }
  cstr_drop(&file_name);
  if (fp)
    xpt_optimizer_write_header(fp);

  cstr applied_file_name = rank == 0
    ? cstr_from_fmt("%s-xpt_seam_applied_diagnostics.csv", mbapp->name)
    : cstr_from_fmt("%s-xpt_seam_applied_diagnostics_rank%d.csv",
        mbapp->name, rank);
  FILE *applied_fp = rank == 0 ? fopen(applied_file_name.str, "w") : NULL;
  if (rank == 0 && !applied_fp) {
    int saved_errno = errno;
    fprintf(stderr, "Unable to open X-point seam applied diagnostic "
      "'%s': %s\n", applied_file_name.str, strerror(saved_errno));
  }
  cstr_drop(&applied_file_name);
  if (applied_fp)
    xpt_optimizer_write_applied_header(applied_fp);

  for (int p=0; p<num_pairs; ++p) {
    const struct xpt_optimizer_pair *pair = &pairs[p];
    const struct gkyl_gk_block_geom_info *bgi0 =
      gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, pair->bid[0]);
    const struct gkyl_gk_block_geom_info *bgi1 = pair->bid[1] >= 0
      ? gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, pair->bid[1]) : 0;

    if (!pair->valid) {
      if (fp)
        xpt_optimizer_write_pair_not_evaluated(fp, mbapp->name, rank, p,
          pair, bgi0, bgi1);
      continue;
    }

    struct xpt_optimizer_search search =
      xpt_optimizer_run_search(mbinp, mbapp, pair, efit);
    if (fp)
      for (int c=0; c<search.count; ++c) {
        const char *status = search.candidates[c].valid ? "accepted"
          : xpt_optimizer_reject_name(search.candidates[c].reject_reason);
        xpt_optimizer_write_candidate_row(fp, mbapp->name, rank, p, pair,
          bgi0, bgi1, &search.candidates[c], status,
          c == search.zero_index, search.zero_objective);
      }

    bool has_recommendation = search.recommended_index >= 0;
    const struct xpt_optimizer_candidate *best = has_recommendation
      ? &search.candidates[search.recommended_index] : 0;
    double relative_improvement =
      has_recommendation && search.zero_objective > 0.0
        ? (search.zero_objective-best->objective)/search.zero_objective
        : 0.0;
    bool useful = has_recommendation &&
      relative_improvement >= XPT_APPLY_MIN_RELATIVE_IMPROVEMENT;

    enum xpt_apply_status apply_status;
    bool applied = false;
    if (!has_recommendation)
      apply_status = XPT_APPLY_NO_VALID_CANDIDATE;
    else if (!useful)
      apply_status = XPT_APPLY_NOT_USEFUL;
    else {
      apply_status = xpt_optimizer_apply_selection(mbinp, mbapp, pair,
        best->coefficient);
      applied = apply_status == XPT_APPLY_APPLIED;
    }

    if (applied_fp)
      xpt_optimizer_write_applied_row(applied_fp, mbapp->name, p, pair,
        search.zero_objective, has_recommendation ? best->coefficient : 0.0,
        has_recommendation ? best->objective : search.zero_objective,
        relative_improvement, useful, applied, apply_status);

    if (rank == 0) {
      if (applied)
        fprintf(stdout, "X-point seam optimizer: pair (block %d edge %d)-"
          "(block %d edge %d) APPLIED delta_s = %.6e m to production "
          "geometry (objective %.6e -> %.6e, %.2f%% reduction).\n",
          pair->bid[0], pair->edge[0], pair->bid[1], pair->edge[1],
          best->coefficient, search.zero_objective, best->objective,
          100.0*relative_improvement);
      else if (has_recommendation)
        fprintf(stdout, "X-point seam optimizer: pair (block %d edge %d)-"
          "(block %d edge %d) found delta_s = %.6e m (objective %.6e -> "
          "%.6e, %.2f%% reduction) but %s; keeping the straight map.\n",
          pair->bid[0], pair->edge[0], pair->bid[1], pair->edge[1],
          best->coefficient, search.zero_objective, best->objective,
          100.0*relative_improvement,
          apply_status == XPT_APPLY_NOT_USEFUL
            ? "this did not clear the usefulness threshold"
            : "the production re-validation failed");
    }
  }

  if (fp)
    fclose(fp);
  if (applied_fp)
    fclose(applied_fp);
  for (int b=0; b<mbapp->num_local_blocks; ++b)
    if (efit[b])
      gkyl_efit_release(efit[b]);
}

static const struct gkyl_array*
gyrokinetic_multib_surf_quantity(const struct gk_geom_surf *geo_surf, const char *name)
{
  if (strcmp(name, "jacobgeo") == 0) return geo_surf->jacobgeo;
  if (strcmp(name, "jacobgeo_signed") == 0)
    return geo_surf->jacobgeo_signed;
  if (strcmp(name, "bmag") == 0) return geo_surf->bmag;
  if (strcmp(name, "jacobtot_inv") == 0) return geo_surf->jacobtot_inv;
  if (strcmp(name, "B3") == 0) return geo_surf->B3;
  if (strcmp(name, "cmag") == 0) return geo_surf->cmag;
  if (strcmp(name, "lenr") == 0) return geo_surf->lenr;
  return 0;
}

static void
gyrokinetic_multib_copy_surf_quantity_skin(gkyl_gyrokinetic_multib_app *app,
  int dir, const char *quantity, struct gkyl_array **side)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[dir];
    const struct gkyl_array *arr = gyrokinetic_multib_surf_quantity(geo_surf, quantity);

    gkyl_array_clear(side[b], 0.0);
    gkyl_array_copy_range(side[b], arr, &sbapp->local_lower_skin[dir]);
    gkyl_array_copy_range_to_range(side[b], arr,
      &sbapp->local_upper_skin[dir], &sbapp->local_upper_ghost[dir]);
  }
}

static void
gyrokinetic_multib_app_write_interface_partner_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-interface_partner_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-interface_partner_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,partner_block,partner_dir,partner_edge,quantity,num_cells,num_coeff,local_min,local_max,local_mean,partner_min,partner_max,partner_mean,diff_min,diff_max,diff_l1,diff_l2,diff_maxabs,rel_l1,rel_l2\n");

  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const char *quantities[] = { "jacobgeo", "jacobgeo_signed", "bmag", "jacobtot_inv", "B3", "cmag", "lenr" };

  for (int d=0; d<cdim; ++d) {
    for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
      struct gkyl_array *side[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
        const struct gkyl_array *arr =
          gyrokinetic_multib_surf_quantity(&sbapp->gk_geom->geo_surf[d], quantities[q]);
        side[b] = mkarr(app->use_gpu, arr->ncomp, arr->size);
      }

      gyrokinetic_multib_copy_surf_quantity_skin(app, d, quantities[q], side);
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_conf->send, app->mbcc_sync_conf->recv, side, side);

      for (int b=0; b<app->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
        int bid = app->local_blocks[b];
        int ncomp = side[b]->ncomp;

        for (int e=0; e<2; ++e) {
          const struct gkyl_target_edge *partner =
            &app->block_topo->conn[bid].connections[d][e];
          if (partner->edge == GKYL_PHYSICAL)
            continue;

          const struct gkyl_range *ghost_range = e == 0
            ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

          double local_min = DBL_MAX, local_max = -DBL_MAX, local_sum = 0.0;
          double partner_min = DBL_MAX, partner_max = -DBL_MAX, partner_sum = 0.0;
          double diff_min = DBL_MAX, diff_max = -DBL_MAX, diff_l1 = 0.0, diff_l2 = 0.0, diff_maxabs = 0.0;
          double partner_l1 = 0.0, partner_l2 = 0.0;
          long count = 0, cell_count = 0;

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, ghost_range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;

            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *local_val = gkyl_array_cfetch(side[b], skin_loc);
            const double *partner_val = gkyl_array_cfetch(side[b], ghost_loc);

            for (int k=0; k<ncomp; ++k) {
              double lv = local_val[k], pv = partner_val[k];
              if (!jacobgeo_ratio_diag_coeff_is_set(lv) || !jacobgeo_ratio_diag_coeff_is_set(pv))
                continue;

              double diff = lv-pv;
              local_min = GKYL_MIN2(local_min, lv);
              local_max = GKYL_MAX2(local_max, lv);
              local_sum += lv;
              partner_min = GKYL_MIN2(partner_min, pv);
              partner_max = GKYL_MAX2(partner_max, pv);
              partner_sum += pv;
              diff_min = GKYL_MIN2(diff_min, diff);
              diff_max = GKYL_MAX2(diff_max, diff);
              diff_l1 += fabs(diff);
              diff_l2 += diff*diff;
              diff_maxabs = GKYL_MAX2(diff_maxabs, fabs(diff));
              partner_l1 += fabs(pv);
              partner_l2 += pv*pv;
              count += 1;
            }
            cell_count += 1;
          }

          if (count > 0) {
            double rel_l1 = partner_l1 > 0.0 ? diff_l1/partner_l1 : 0.0;
            double rel_l2 = partner_l2 > 0.0 ? sqrt(diff_l2/partner_l2) : 0.0;
            fprintf(fp, "%s,%d,%d,%d,%s,%d,%d,%s,%s,%ld,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
              app->name, rank, bid, d, edge_name_from_index(e),
              partner->bid, partner->dir, oriented_edge_name(partner->edge), quantities[q],
              cell_count, ncomp, local_min, local_max, local_sum/count,
              partner_min, partner_max, partner_sum/count,
              diff_min, diff_max, diff_l1, sqrt(diff_l2), diff_maxabs,
              rel_l1, rel_l2);
          }
        }
      }

      for (int b=0; b<app->num_local_blocks; ++b)
        gkyl_array_release(side[b]);
    }
  }

  fclose(fp);
}

void
gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    struct gkyl_gk_geometry_inp geometry_inp =
      gyrokinetic_multib_geometry_inp_from_block(bgi, app->block_comms[bid]);
    gkyl_gyrokinetic_app_write_geometry(app->singleb_apps[b], &geometry_inp);
  }

  gyrokinetic_multib_app_write_jacobgeo_ratio_diag(app);
  gyrokinetic_multib_app_write_flux_weight_rescale_diag(app);
  gyrokinetic_multib_app_write_flux_weight_modal_product_diag(app);
  gyrokinetic_multib_app_write_interface_partner_diag(app);
  gyrokinetic_multib_app_write_interface_pointwise_diag(app);
}

void
gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_field(app->singleb_apps[b], tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_calc_field_energy(app->singleb_apps[b], tm);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_write_field_energy(app->singleb_apps[b]);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_eirene(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{

  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_diagnostics(app->singleb_apps[b], tm, frame);

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  cstr fileNm = cstr_from_fmt("%snew_data_flag", sbapp->eirene->info.output_data_path);
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    FILE *fp = fopen(fileNm.str, "w");
    if (fp == NULL)
        return;
    fprintf(fp, "%d\n", frame);
    fclose(fp);
  }
  cstr_drop(&fileNm);
}

void
gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(app->singleb_apps[b], tm);
}

void
gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(app->singleb_apps[b]);
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. BGK Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx);
  }
}

//
// ............. LTE outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lbo_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_bgk_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_bgk_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_drag(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Functions that group several outputs for a single species ............... //
//
void
gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_rad_drag(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_bgk_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(app, i, tm);
  }

  gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(app, tm);
}

void
gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(app, i);
    gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(app, i);
  }

  gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(app);
}

void
gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_conf(app, i, tm, frame);
  }

  gkyl_gyrokinetic_multib_app_write_eirene(app, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_multib_app_write_conf(app, tm, frame);
}

//
// ............. End of write functions ............... //
// 

struct gkyl_update_status
gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = gyrokinetic_multib_update_ssp_rk3(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.time_loop_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());

  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    app->stat.n_iter_corr[i] = 0;
    app->stat.num_corr[i] = 0;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    app->stat.neut_n_iter_corr[i] = 0;
    app->stat.neut_num_corr[i] = 0;
  }

  for (int b=0; b<app->num_local_blocks; ++b) {
    // Add time spent on various operations for each local block.
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    struct gkyl_gyrokinetic_stat sb_stat = gkyl_gyrokinetic_app_stat(sbapp);

    app->stat.field_phi_rhs_tm += sb_stat.field_phi_rhs_tm;

    app->stat.species_lte_tm += sb_stat.species_lte_tm;
    app->stat.species_coll_mom_tm += sb_stat.species_coll_mom_tm;
    app->stat.species_coll_tm += sb_stat.species_coll_tm;
    app->stat.species_rad_mom_tm += sb_stat.species_rad_mom_tm;
    app->stat.species_rad_tm += sb_stat.species_rad_tm;
    app->stat.species_react_mom_tm += sb_stat.species_react_mom_tm;
    app->stat.species_react_tm += sb_stat.species_react_tm;
    app->stat.species_omega_cfl_tm += sb_stat.species_omega_cfl_tm;

    app->stat.neut_species_lte_tm += sb_stat.neut_species_lte_tm;
    app->stat.neut_species_coll_mom_tm += sb_stat.neut_species_coll_mom_tm;
    app->stat.neut_species_coll_tm += sb_stat.neut_species_coll_tm;
    app->stat.neut_species_react_mom_tm += sb_stat.neut_species_react_mom_tm;
    app->stat.neut_species_react_tm += sb_stat.neut_species_react_tm;
    app->stat.neut_species_omega_cfl_tm += sb_stat.neut_species_omega_cfl_tm;

    app->stat.species_bc_tm += sb_stat.species_bc_tm;
    app->stat.neut_species_bc_tm += sb_stat.neut_species_bc_tm;

    app->stat.species_io_tm += sb_stat.species_io_tm;
    app->stat.species_diag_calc_tm += sb_stat.species_diag_calc_tm;
    app->stat.species_diag_io_tm += sb_stat.species_diag_io_tm;

    app->stat.neut_species_io_tm += sb_stat.neut_species_io_tm;
    app->stat.neut_species_diag_calc_tm += sb_stat.neut_species_diag_calc_tm;
    app->stat.neut_species_diag_io_tm += sb_stat.neut_species_diag_io_tm;

    app->stat.n_species_omega_cfl += sb_stat.n_species_omega_cfl;
    app->stat.n_mom += sb_stat.n_mom;    
    app->stat.n_diag += sb_stat.n_diag;
    app->stat.n_io += sb_stat.n_io;
    app->stat.n_diag_io += sb_stat.n_diag_io;

    app->stat.n_neut_species_omega_cfl += sb_stat.n_neut_species_omega_cfl;
    app->stat.n_neut_mom += sb_stat.n_neut_mom;  
    app->stat.n_neut_diag += sb_stat.n_neut_diag;
    app->stat.n_neut_io += sb_stat.n_neut_io;
    app->stat.n_neut_diag_io += sb_stat.n_neut_diag_io;

    for (int i=0; i<app->num_species; ++i) {
      app->stat.n_iter_corr[i] += sb_stat.n_iter_corr[i];
      app->stat.num_corr[i] += sb_stat.num_corr[i];
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      app->stat.neut_n_iter_corr[i] += sb_stat.neut_n_iter_corr[i];
      app->stat.neut_num_corr[i] += sb_stat.neut_num_corr[i];
    }

    // Timers not yet computed in app directly.
    app->stat.time_rate_diags_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;
    app->stat.pos_shift_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;
    app->stat.io_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
      + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
      + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;
    
    // Additions of several timers.
    app->stat.fwd_euler_sum_tm += sb_stat.species_coll_mom_tm + sb_stat.species_react_mom_tm + sb_stat.neut_species_coll_mom_tm
      + sb_stat.neut_species_react_mom_tm + sb_stat.species_rad_mom_tm + sb_stat.species_gyroavg_tm + sb_stat.species_collisionless_tm
      + sb_stat.species_coll_tm + sb_stat.species_diffusion_tm + sb_stat.species_rad_tm + sb_stat.species_react_tm
      + sb_stat.species_bflux_calc_tm+sb_stat.species_bflux_moms_tm + sb_stat.species_omega_cfl_tm + sb_stat.species_src_tm
      + sb_stat.neut_species_collisionless_tm + sb_stat.neut_species_coll_tm + sb_stat.neut_species_react_tm
      + sb_stat.neut_species_omega_cfl_tm + sb_stat.neut_species_src_tm + sb_stat.dfdt_dt_reduce_tm;

    app->stat.time_rate_diags_sum_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;

    app->stat.pos_shift_sum_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;

    app->stat.io_sum_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
    + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
    + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;

  }

  app->stat.fwd_euler_sum_tm += app->stat.fwd_euler_step_f_tm + app->stat.dfdt_dt_reduce_tm;
  app->stat.bc_sum_tm = app->stat.species_bc_tm + app->stat.neut_species_bc_tm;
  app->stat.field_sum_tm = app->stat.field_phi_rhs_tm  + app->stat.field_phi_solve_tm;
  app->stat.time_stepper_sum_tm = app->stat.fwd_euler_tm  + app->stat.field_tm + app->stat.bc_tm
    + app->stat.time_rate_diags_tm + app->stat.pos_shift_tm + app->stat.time_stepper_arithmetic_tm;

  return app->stat;
}

static inline
double
ratio_to_percent(double num, double den, double alt)
{
  return den > 1e-12 ? 100.*num/den : alt;
}

void
gkyl_gyrokinetic_multib_app_print_timings(gkyl_gyrokinetic_multib_app* app, FILE *iostream)
{
  struct gkyl_gyrokinetic_stat *stat = &app->stat;

  double bflux_tm = stat->species_bflux_calc_tm+stat->species_bflux_moms_tm;

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "Timing:\n");
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - Time loop:                         %.4e sec.\n", stat->time_loop_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Forward Euler:                   %.4e sec. / %4.2f %%.\n", stat->fwd_euler_tm, ratio_to_percent(stat->fwd_euler_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (charged):   %.4e sec. / %4.2f %%.\n", stat->species_coll_mom_tm          , ratio_to_percent(stat->species_coll_mom_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (charged):    %.4e sec. / %4.2f %%.\n", stat->species_react_mom_tm         , ratio_to_percent(stat->species_react_mom_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (neutral):   %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_mom_tm     , ratio_to_percent(stat->neut_species_coll_mom_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_react_mom_tm    , ratio_to_percent(stat->neut_species_react_mom_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation moments:             %.4e sec. / %4.2f %%.\n", stat->species_rad_mom_tm           , ratio_to_percent(stat->species_rad_mom_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species gyroaverage:           %.4e sec. / %4.2f %%.\n", stat->species_gyroavg_tm           , ratio_to_percent(stat->species_gyroavg_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (charged): %.4e sec. / %4.2f %%.\n", stat->species_collisionless_tm     , ratio_to_percent(stat->species_collisionless_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (charged):     %.4e sec. / %4.2f %%.\n", stat->species_coll_tm              , ratio_to_percent(stat->species_coll_tm              ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Diffusion (charged):           %.4e sec. / %4.2f %%.\n", stat->species_diffusion_tm         , ratio_to_percent(stat->species_diffusion_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation terms:               %.4e sec. / %4.2f %%.\n", stat->species_rad_tm               , ratio_to_percent(stat->species_rad_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (charged):      %.4e sec. / %4.2f %%.\n", stat->species_react_tm             , ratio_to_percent(stat->species_react_tm             ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Boundary fluxes (charged):     %.4e sec. / %4.2f %%.\n", bflux_tm                           , ratio_to_percent(bflux_tm                           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (charged):           %.4e sec. / %4.2f %%.\n", stat->species_omega_cfl_tm         , ratio_to_percent(stat->species_omega_cfl_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (charged):             %.4e sec. / %4.2f %%.\n", stat->species_src_tm               , ratio_to_percent(stat->species_src_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (neutral): %.4e sec. / %4.2f %%.\n", stat->neut_species_collisionless_tm, ratio_to_percent(stat->neut_species_collisionless_tm,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_tm         , ratio_to_percent(stat->neut_species_coll_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (neutral):      %.4e sec. / %4.2f %%.\n", stat->neut_species_react_tm        , ratio_to_percent(stat->neut_species_react_tm        ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (neutral):           %.4e sec. / %4.2f %%.\n", stat->neut_species_omega_cfl_tm    , ratio_to_percent(stat->neut_species_omega_cfl_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_src_tm          , ratio_to_percent(stat->neut_species_src_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Time step reduction:           %.4e sec. / %4.2f %%.\n", stat->dfdt_dt_reduce_tm            , ratio_to_percent(stat->dfdt_dt_reduce_tm            ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Step f:                        %.4e sec. / %4.2f %%.\n", stat->fwd_euler_step_f_tm          , ratio_to_percent(stat->fwd_euler_step_f_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->fwd_euler_sum_tm, stat->fwd_euler_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field solves:                    %.4e sec. / %4.2f %%.\n", stat->field_tm, ratio_to_percent(stat->field_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn RHS:                   %.4e sec. / %4.2f %%.\n", stat->field_phi_rhs_tm  , ratio_to_percent(stat->field_phi_rhs_tm  ,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn solve:                 %.4e sec. / %4.2f %%.\n", stat->field_phi_solve_tm, ratio_to_percent(stat->field_phi_solve_tm,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->field_sum_tm, stat->field_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Boundary conditions::            %.4e sec. / %4.2f %%.\n", stat->bc_tm, ratio_to_percent(stat->bc_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_bc_tm     , ratio_to_percent(stat->species_bc_tm     ,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_bc_tm, ratio_to_percent(stat->neut_species_bc_tm,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->bc_sum_tm, stat->bc_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time rate diagnostics:           %.4e sec. / %4.2f %%.\n", stat->time_rate_diags_tm, ratio_to_percent(stat->time_rate_diags_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Charged species:               %.4e sec. / %4.2f %%.\n", stat->fdot_tm  , ratio_to_percent(stat->fdot_tm  ,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi:                           %.4e sec. / %4.2f %%.\n", stat->phidot_tm, ratio_to_percent(stat->phidot_tm,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->time_rate_diags_sum_tm, stat->time_rate_diags_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Positivity shift:                %.4e sec. / %4.2f %%.\n", stat->pos_shift_tm, ratio_to_percent(stat->pos_shift_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_pos_shift_tm     , ratio_to_percent(stat->species_pos_shift_tm     ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_pos_shift_tm, ratio_to_percent(stat->neut_species_pos_shift_tm,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Quasineutrality:               %.4e sec. / %4.2f %%.\n", stat->pos_shift_quasineut_tm   , ratio_to_percent(stat->pos_shift_quasineut_tm   ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->pos_shift_sum_tm, stat->pos_shift_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time stepper arithmetic:         %.4e sec. / %4.2f %%.\n", stat->time_stepper_arithmetic_tm, ratio_to_percent(stat->time_stepper_arithmetic_tm,stat->time_loop_tm, 0.0));

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->time_stepper_sum_tm, stat->time_loop_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - I/O:                               %.4e sec.\n", stat->io_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (charged):               %.4e sec. / %4.2f %%.\n", stat->species_io_tm            , ratio_to_percent(stat->species_io_tm            , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (charged):     %.4e sec. / %4.2f %%.\n", stat->species_diag_calc_tm     , ratio_to_percent(stat->species_diag_calc_tm     , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (charged):    %.4e sec. / %4.2f %%.\n", stat->species_diag_io_tm       , ratio_to_percent(stat->species_diag_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (neutral):               %.4e sec. / %4.2f %%.\n", stat->neut_species_io_tm       , ratio_to_percent(stat->neut_species_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_calc_tm, ratio_to_percent(stat->neut_species_diag_calc_tm, stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_io_tm  , ratio_to_percent(stat->neut_species_diag_io_tm  , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field write:                     %.4e sec. / %4.2f %%.\n", stat->field_io_tm              , ratio_to_percent(stat->field_io_tm              , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag calc:                 %.4e sec. / %4.2f %%.\n", stat->field_diag_calc_tm       , ratio_to_percent(stat->field_diag_calc_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag write:                %.4e sec. / %4.2f %%.\n", stat->field_diag_io_tm         , ratio_to_percent(stat->field_diag_io_tm         , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Common write:                    %.4e sec. / %4.2f %%.\n", stat->app_io_tm                , ratio_to_percent(stat->app_io_tm                , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->io_sum_tm, stat->io_tm, 100.0));
}

void
gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    gkyl_gyrokinetic_app_stat_write(sbapp);
  }
}

void
gkyl_gyrokinetic_multib_app_write_dt(gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "dt");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, "dt");

    struct timespec wtm = gkyl_wall_clock();
    if (app->is_first_dt_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Time step size." }
      };
      int io_meta_len[] = {app->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(app->dts, fileNm, mt);
      app->is_first_dt_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(app->dts, fileNm);
    }
    app->stat.app_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 1;
  }
  gkyl_dynvec_clear(app->dts);
}

void
gkyl_gyrokinetic_multib_app_save_dt(gkyl_gyrokinetic_multib_app* app, double tm, double dt)
{
  gkyl_dynvec_append(app->dts, tm, &dt);
}

void
gkyl_gyrokinetic_multib_app_release_geom(gkyl_gyrokinetic_multib_app* mbapp)
{
  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release_geom(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_msgpack_map_elem_release(mbapp->io_meta_basic_len, mbapp->io_meta_basic);

  gkyl_free(mbapp);
}

void
gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  for (int i=0; i<mbapp->num_neut_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_neut[i].send);
    gkyl_free(mbapp->mbcc_sync_neut[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_neut);

  for (int i=0; i<mbapp->num_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_charged[i].send);
    gkyl_free(mbapp->mbcc_sync_charged[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_charged);

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  gk_multib_field_release(mbapp->field);

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_dynvec_release(mbapp->dts);

  gkyl_free(mbapp);
}
