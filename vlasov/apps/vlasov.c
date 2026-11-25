#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_vlasov_priv.h>
#include <gkyl_app_priv.h>
#include <gkyl_wv_euler.h>

#include <mpack.h>

// returned gkyl_array_meta must be freed using vlasov_array_meta_release
struct gkyl_msgpack_data *
vlasov_array_meta_new(struct vlasov_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);

  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.stime);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

  // the staggered operator split's last and suggested steps (for a restart)
  if (meta.dt_prev > 0.0 || meta.dt_suggested > 0.0) {
    mpack_write_cstr(&writer, "dt_prev");
    mpack_write_double(&writer, meta.dt_prev);
    mpack_write_cstr(&writer, "dt_suggested");
    mpack_write_double(&writer, meta.dt_suggested);
  }

  // the phase-space split of a hybrid basis
  if (meta.num_vdim > 0) {
    mpack_write_cstr(&writer, "num_cdim");
    mpack_write_i64(&writer, meta.num_cdim);
    mpack_write_cstr(&writer, "num_vdim");
    mpack_write_i64(&writer, meta.num_vdim);
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

void
vlasov_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) {
    return;
  }
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

struct vlasov_output_meta
vlasov_meta_from_mpack(struct gkyl_msgpack_data *mt)
{
  struct vlasov_output_meta meta = {.frame = 0, .stime = 0.0, .dt_prev = 0.0, .dt_suggested = 0.0};

  if (mt->meta_sz > 0) {
    mpack_tree_t tree;
    mpack_tree_init_data(&tree, mt->meta, mt->meta_sz);
    mpack_tree_parse(&tree);
    mpack_node_t root = mpack_tree_root(&tree);

    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    meta.stime = mpack_node_double(tm_node);

    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    meta.frame = mpack_node_i64(fr_node);

    mpack_node_t po_node = mpack_node_map_cstr(root, "polyOrder");
    meta.poly_order = mpack_node_i64(po_node);

    mpack_node_t bt_node = mpack_node_map_cstr(root, "basisType");
    char *basis_type = mpack_node_cstr_alloc(bt_node, 64);
    strcpy(meta.basis_type_nm, basis_type);
    MPACK_FREE(basis_type);

    mpack_node_t dtp_node = mpack_node_map_cstr_optional(root, "dt_prev");
    meta.dt_prev = mpack_node_is_missing(dtp_node) ? 0.0 : mpack_node_double(dtp_node);
    mpack_node_t dts_node = mpack_node_map_cstr_optional(root, "dt_suggested");
    meta.dt_suggested = mpack_node_is_missing(dts_node) ? 0.0 : mpack_node_double(dts_node);

    mpack_tree_destroy(&tree);
  }
  return meta;
}

// Collisions of a PKPM species: LBO or BGK.
static bool
vlasov_pkpm_has_collisions(const struct gkyl_vlasov_collisions *coll)
{
  return coll->collision_id == GKYL_LBO_COLLISIONS || coll->collision_id == GKYL_BGK_COLLISIONS;
}

// An Euler fluid with gas_gamma = 5/3: the fluid with three velocity degrees of freedom.
static bool
vlasov_fluid_is_euler_3dof(const struct gkyl_vlasov_fluid_species *fluid)
{
  return fluid->equation->type == GKYL_EQN_EULER &&
         fabs(gkyl_wv_euler_gas_gamma(fluid->equation) - 5.0 / 3.0) <= 1.0e-12;
}

// The species a species names must exist and be usable by its operators.
// Collision partners: kinetic species collide with kinetic species through their operators;
// PKPM species and fluid species collide with each other through the exchange of their fluid
// states, which has three degrees of freedom, so a fluid in it is an Euler fluid with gas_gamma
// = 5/3. A kinetic and a PKPM species collide through the LBO operator of the kinetic species,
// which carries the exchange of the pair in the stages. A PKPM species in a collision pair needs
// LBO or BGK collisions of its own. The partners of an adaptive source are read through the data
// of the kinetic route.
// The source species of a boundary-flux source needs wall losses (a kinetic or PKPM species); a
// PKPM species forms its source rate before the right-hand sides, when only the losses of PKPM
// species are known, so its source species must be a PKPM species.
static void
vlasov_check_species_partners(const gkyl_vlasov_app *app, const struct gkyl_vm *vm)
{
  for (int i = 0; i < vm->num_species; ++i) {
    const struct gkyl_vlasov_species *inp = &vm->species[i];
    bool is_fluid = inp->type == GKYL_SPECIES_FLUID;
    const struct gkyl_vlasov_collisions *coll = is_fluid ? &inp->fluid.collisions :
                                                           &inp->kinetic.collisions;
    for (int j = 0; j < coll->num_cross_collisions; ++j) {
      const struct vlasov_species *other = vlasov_find_species(app, coll->collide_with[j]);
      bool both_kinetic = other && inp->type == GKYL_SPECIES_VLASOV &&
                          other->type == GKYL_SPECIES_VLASOV;
      bool both_exchange = other && inp->type != GKYL_SPECIES_VLASOV &&
                           other->type != GKYL_SPECIES_VLASOV;
      bool kinetic_pkpm = other && inp->type != GKYL_SPECIES_FLUID &&
                          other->type != GKYL_SPECIES_FLUID && !both_kinetic && !both_exchange;
      if (!both_kinetic && !both_exchange && !kinetic_pkpm) {
        gkyl_exit(
          "gkyl_vlasov_app_new: a kinetic species collides with kinetic and PKPM species, a fluid species with PKPM and fluid species."
        );
      }
      const struct gkyl_vlasov_species *other_inp = &vm->species[other - app->species];
      const struct gkyl_vlasov_collisions *other_coll = other->type == GKYL_SPECIES_FLUID ?
                                                          &other_inp->fluid.collisions :
                                                          &other_inp->kinetic.collisions;
      if ((inp->type == GKYL_SPECIES_PKPM && !vlasov_pkpm_has_collisions(coll)) ||
          (other->type == GKYL_SPECIES_PKPM && !vlasov_pkpm_has_collisions(other_coll))) {
        gkyl_exit(
          "gkyl_vlasov_app_new: a PKPM species in a collision pair needs LBO or BGK collisions."
        );
      }
      if (kinetic_pkpm) {
        const struct gkyl_vlasov_collisions *kin_coll =
          inp->type == GKYL_SPECIES_VLASOV ? coll :
                                             &vm->species[other - app->species].kinetic.collisions;
        const struct gkyl_vlasov_collisions *pkpm_coll =
          inp->type == GKYL_SPECIES_PKPM ? coll :
                                           &vm->species[other - app->species].kinetic.collisions;
        if (kin_coll->collision_id != GKYL_LBO_COLLISIONS || app->use_gpu ||
            (pkpm_coll->collision_id == GKYL_BGK_COLLISIONS && pkpm_coll->is_implicit)) {
          gkyl_exit(
            "gkyl_vlasov_app_new: collisions between a kinetic and a PKPM species need LBO collisions of the kinetic species and explicit collisions of the PKPM species, on the CPU."
          );
        }
        if (!coll->cross_nu[j]) {
          gkyl_exit(
            "gkyl_vlasov_app_new: collisions between a kinetic and a PKPM species need the cross-collision frequencies (cross_nu) of both."
          );
        }
      }
    }
    if (is_fluid && coll->num_cross_collisions > 0 && !vlasov_fluid_is_euler_3dof(&inp->fluid)) {
      gkyl_exit(
        "gkyl_vlasov_app_new: a fluid species with collisions must be an Euler fluid with gas_gamma = 5/3."
      );
    }
    if (is_fluid) {
      continue;
    }
    const struct gkyl_vlasov_source *src = &inp->kinetic.source;
    if (src->source_id == GKYL_BFLUX_SOURCE) {
      const struct vlasov_species *other = vlasov_find_species(app, src->source_species);
      if (!other || other->type == GKYL_SPECIES_FLUID) {
        gkyl_exit(
          "gkyl_vlasov_app_new: the source species of a boundary-flux source must be a kinetic or a PKPM species."
        );
      }
      if (inp->type == GKYL_SPECIES_PKPM && other->type != GKYL_SPECIES_PKPM) {
        gkyl_exit(
          "gkyl_vlasov_app_new: the source species of the boundary-flux source of a PKPM species must be a PKPM species."
        );
      }
    }
    if (src->source_id == GKYL_PROJ_ADAPT_DENSITY_SOURCE) {
      for (int j = 0; j < src->num_cross_source; ++j) {
        const struct vlasov_species *other = vlasov_find_species(app, src->source_with[j]);
        if (!other || other->type != inp->type) {
          gkyl_exit(
            "gkyl_vlasov_app_new: the partners of an adaptive source must be of the type of the species."
          );
        }
      }
    }
  }
}

gkyl_vlasov_app *
gkyl_vlasov_app_new(struct gkyl_vm *vm)
{
  disable_denorm_float();

  assert(vm->num_species <= GKYL_MAX_SPECIES);

  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));

  int ntot = vm->num_species;

  int cdim = app->cdim = vm->cdim;
  int vdim = app->vdim = vm->vdim;
  int pdim = cdim + vdim;
  int poly_order = app->poly_order = vm->poly_order;
  app->num_species = ntot;

  double cfl_frac = vm->cfl_frac == 0 ? 1.0 : vm->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm->parallelism.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = vm->num_periodic_dir;
  for (int d = 0; d < cdim; ++d) {
    app->periodic_dirs[d] = vm->periodic_dirs[d];
  }

  strcpy(app->name, vm->name);
  app->tcurr = 0.0; // reset on init

  // PKPM species need polynomial order 1 in configuration space. The serendipity and tensor
  // bases are the same there; the type selects the phase-space basis of the kinetic species of
  // the run (tensor: the hybrid basis, which PKPM species always use).
  app->has_pkpm = false;
  app->num_pkpm_species = 0; // gathered in vlasov_species_init_pkpm
  for (int i = 0; i < ntot; ++i) {
    app->has_pkpm = app->has_pkpm || vm->species[i].type == GKYL_SPECIES_PKPM;
  }
  if (app->has_pkpm && poly_order != 1) {
    gkyl_exit("gkyl_vlasov_app_new: PKPM species need polyOrder = 1.");
  }
  // A PKPM species has three velocity degrees of freedom, and so must the species it shares a run
  // with: kinetic species with three velocity dimensions, Euler fluids with gas_gamma = 5/3.
  for (int i = 0; app->has_pkpm && i < ntot; ++i) {
    if (vm->species[i].type == GKYL_SPECIES_VLASOV && vm->vdim != 3) {
      gkyl_exit(
        "gkyl_vlasov_app_new: a kinetic species in a run with PKPM species needs three velocity dimensions."
      );
    }
    if (vm->species[i].type == GKYL_SPECIES_FLUID &&
        !vlasov_fluid_is_euler_3dof(&vm->species[i].fluid)) {
      gkyl_exit(
        "gkyl_vlasov_app_new: a fluid species in a run with PKPM species must be an Euler fluid with gas_gamma = 5/3."
      );
    }
  }
  // So does the exchange between colliding species with a fluid aspect, which is at the nodes.
  bool has_fluid_collisions = false;
  for (int i = 0; i < ntot; ++i) {
    has_fluid_collisions = has_fluid_collisions ||
                           (vm->species[i].type == GKYL_SPECIES_FLUID &&
                            vm->species[i].fluid.collisions.num_cross_collisions > 0);
  }
  if (has_fluid_collisions && (poly_order != 1 || app->use_gpu)) {
    gkyl_exit("gkyl_vlasov_app_new: collisions of fluid species need polyOrder = 1 on the CPU.");
  }

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  } else {
    app->basis_on_dev = &app->basis;
  }

  // basis functions
  switch (vm->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->basis, cdim, poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev, cdim, poly_order);
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->basis, cdim, poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_tensor_cu_dev(app->basis_on_dev, cdim, poly_order);
      }
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, vm->lower, vm->upper, vm->cells);

  int ghost[] = {1, 1, 1};
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (vm->parallelism.comm == 0) {
    int cuts[3] = {1, 1, 1};
    app->decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);

    // The general diffusion tensor has cross terms whose stencil reaches corner ghost
    // cells, so these need to be synced when a fluid species uses one.
    bool sync_corners = false;
    for (int i = 0; i < vm->num_species; ++i) {
      if (vm->species[i].type == GKYL_SPECIES_FLUID && vm->species[i].fluid.diffusion.Dij) {
        sync_corners = true;
      }
    }

    app->comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
      .decomp = app->decomp,
      .use_gpu = app->use_gpu,
      .sync_corners = sync_corners,
    });

    // Global and local ranges are same, and so just copy them.
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));
  } else {
    // Create decomp.
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->cdim, vm->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(vm->parallelism.comm, 0, app->decomp);

    // Create local and local_ext.
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_create_ranges(&app->decomp->ranges[rank], ghost, &app->local_ext, &app->local);
  }

  // local skin and ghost ranges for configuration space fields
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost
    );
  }

  // Configuration space geometry initialization
  // The geometry has a p=1 DG representation.
  app->c2p_ctx = app->mapc2p = 0;
  app->has_mapc2p = vm->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = vm->c2p_ctx;
    app->mapc2p = vm->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, cdim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, cdim * basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p =
      gkyl_eval_on_nodes_new(&app->grid, &basis, cdim, vm->mapc2p, vm->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    cstr fileNm = cstr_from_fmt("%s-mapc2p.gkyl", app->name);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, c2p, fileNm.str);
    cstr_drop(&fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom =
    gkyl_wave_geom_new(&app->grid, &app->local_ext, app->mapc2p, app->c2p_ctx, app->use_gpu);

  // Create the configuration-space position map (C^0 piecewise linear,
  // diagonal). Always created; directions without a user map are the identity.
  // Shared by all species, which acquire a reference.
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    inp_pmap[d].eval_pmap = vm->mapc2p_pos[d].mapc2p_pos_func;
    inp_pmap[d].ctx = vm->mapc2p_pos[d].mapc2p_pos_ctx;
  }
  app->pos_map = gkyl_vlasov_position_map_new(
    &app->grid, &app->local, &app->local_ext, &app->basis, inp_pmap, app->use_gpu
  );

  // Finalize the position-map ghost cells. The constructor filled every ghost
  // with the adjacent local-skin value (a copy BC, valid at physical
  // boundaries); now overwrite the periodic and inter-rank ghosts with the
  // wrapped/neighbor values, mirroring the conf-array sync in vm_*_apply_bc.
  struct gkyl_array *pos_map_arrs[] = {
    app->pos_map->pmap, app->pos_map->jacob_pos, app->pos_map->jacob_pos_surf,
    app->pos_map->jacob_pos_gauss
  };
  for (int ia = 0; ia < 4; ++ia) {
    gkyl_comm_array_per_sync(
      app->comm, &app->local, &app->local_ext, app->num_periodic_dir, app->periodic_dirs,
      pos_map_arrs[ia]
    );
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, pos_map_arrs[ia]);
  }

  // The position map is static in time; write it once here (uniform/identity
  // grids included), mirroring the mapc2p write above. Both the p=3 map and its
  // p=0 cell average are written, with metadata built from the map's I/O basis.
  gkyl_vlasov_position_map_write(app->pos_map, app->comm, app->name, "position-map");

  // Conversions to the configuration-nodal, velocity-modal form of the node-based operators of
  // PKPM species and of the collisions of fluid species (CPU runs).
  app->conf_nodal = 0;
  if ((app->has_pkpm || has_fluid_collisions) && !app->use_gpu) {
    // the phase-space basis of the PKPM species (the same for all of them, one velocity dimension)
    struct gkyl_basis phase_basis;
    gkyl_cart_modal_hybrid(&phase_basis, cdim, 1);
    app->conf_nodal = gkyl_dg_conf_nodal_inew(&(struct gkyl_dg_conf_nodal_inp){
      .conf_basis = &app->basis,
      .phase_basis = &phase_basis,
      .use_gpu = false,
    });
  }

  // allocate space to store vlasov-maxwell geometry objects
  app->vm_geom = gkyl_malloc(sizeof(struct vm_geom));
  vm_geom_init(vm, app, app->vm_geom);

  // No field when it is skipped or its type is the null field (vlasov_field_new).
  app->has_field = !(vm->skip_field || vm->field.field_id == GKYL_FIELD_NULL);
  // A field object is always created: the null field when no field is present.
  app->use_explicit_sources = vm->use_explicit_sources;
  app->field = vlasov_field_new(vm, app);

  app->species = ntot > 0 ? gkyl_malloc(sizeof(struct vlasov_species[ntot])) : 0;

  // Construct each species container (allocates its aspect, stores the input,
  // and sets its methods); the aspects are initialized below, once all
  // containers exist, since the initializations look species up by name.
  for (int i = 0; i < ntot; ++i) {
    vlasov_species_init(app, &vm->species[i], &app->species[i]);
  }
  vlasov_check_species_partners(app, vm);
  // Divergence cleaning of the electric field relaxes div(E) to zero, which holds in vacuum
  // only: with a charged species div(E) is its charge density.
  for (int i = 0; app->has_field && i < ntot; ++i) {
    if (app->field->info.elcErrorSpeedFactor != 0.0 && vm->species[i].charge != 0.0) {
      gkyl_exit(
        "gkyl_vlasov_app_new: divergence cleaning of the electric field (elcErrorSpeedFactor) is not supported with charged species."
      );
    }
  }

  // Species initialization runs in three phases: kinetic aspects, then the
  // links between them (emission, cross collisions, sources), then the fluid
  // aspects and their sources. PKPM species build both aspects in one route,
  // then their links and the exchange between them.
  vlasov_species_init_kinetic_aspects(vm, app);
  vlasov_species_link_kinetic_aspects(app);
  vlasov_species_init_fluid_aspects(vm, app);
  if (app->has_pkpm) {
    vlasov_species_init_pkpm(vm, app);
  }
  // Held backgrounds, captured at the end of the initial conditions.
  app->has_background = false;
  for (int i = 0; i < ntot; ++i) {
    app->has_background = app->has_background || vlasov_species_has_background(&app->species[i]);
  }

  // The momentum and energy exchange between the species with a fluid aspect that collide.
  vm_coll_exchange_init(app);
  app->has_coll_stage = app->has_pkpm || app->coll_exchange.num_species > 0;

  // Fluid-EM coupling: supported for the fluid aspects that carry a momentum (PKPM species;
  // fluid species with the Euler or isothermal Euler equations)
  // under the standard Maxwell field. Scalar advection and canonical-PB fluids
  // carry no current and are never coupled; a momentum-carrying fluid under any
  // other field type has no coupling formulation and is rejected. The coupling
  // is implicit (the staggered operator split, see vlasov_update_op_split)
  // unless use_explicit_sources puts the Lorentz force in the fluid RHS and the
  // fluid current on the Ampere RHS of every SSP-RK3 stage.
  int num_coupled_fluids = 0;
  for (int i = 0; i < ntot; ++i) {
    if (vm_fluid_em_coupling_supported(&app->species[i])) {
      num_coupled_fluids += 1;
    }
  }
  app->has_fluid_em_coupling = false;
  app->src_dt_func = vlasov_src_dt_none;
  app->src_omega_p_sq = 0;
  app->op_split_dt_prev = 0.0;
  app->op_split_dt_suggested = 0.0;
  app->op_split_synced = false;
  if (num_coupled_fluids > 0 && app->has_field) {
    if (app->field->field_id != GKYL_FIELD_E_B) {
      gkyl_exit(
        "gkyl_vlasov_app_new: fluid-EM coupling is only supported with the standard Maxwell field (E_B)."
      );
    }
    if (!app->use_explicit_sources) {
      // Kinetic species join the split: the source step takes their current as known and the
      // transport step moves them with its E (see vlasov_update_op_split).
      app->has_fluid_em_coupling = true;
      app->fl_em = vm_fluid_em_coupling_new(app);
      vm_field_set_source_split(app, app->field);
    }
  }
  // Without the split, the currents of all species and the resistive layer are explicit in the
  // stages of a Maxwell field: the time-step resolves the plasma frequency (with the cyclotron
  // frequency in the rate of each coupled fluid) and the rate of the layer. A static field does
  // not respond to them, and with the split the source step resolves them: the plasma frequency
  // the force kernels of the stages read stays zero.
  if (app->has_field && app->field->field_id == GKYL_FIELD_E_B) {
    app->src_omega_p_sq = mkarr(app->use_gpu, 1, app->local_ext.volume);
    gkyl_array_clear(app->src_omega_p_sq, 0.0);
    if (!app->field->info.is_static && !app->has_fluid_em_coupling) {
      app->src_dt_func = vlasov_src_dt_explicit;
      app->src_omega_p_sq_kin = mkarr(app->use_gpu, 1, app->local_ext.volume);
      gkyl_array_clear(app->src_omega_p_sq_kin, 0.0);
      app->src_rho_avg = mkarr(app->use_gpu, 1, app->local_ext.volume);
      app->src_omega_p_sq_max = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) :
                                               gkyl_malloc(sizeof(double));
      app->src_num_kinetic = 0;
      for (int i = 0; i < ntot; ++i) {
        if (app->species[i].type == GKYL_SPECIES_VLASOV && app->species[i].charge != 0.0) {
          app->src_kinetic[app->src_num_kinetic] = i;
          app->src_num_kinetic += 1;
        }
      }
    }
  }

  // Implicit operators of the species (BGK collisions, temperature sources) are
  // applied inside the SSP-RK3 step, see vlasov_stage_implicit.
  app->has_implicit_species_ops = false;
  for (int i = 0; i < ntot; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (app->species[i].type == GKYL_SPECIES_VLASOV &&
        (vms->info.collisions.is_implicit || vms->src_bgk.is_implicit)) {
      app->has_implicit_species_ops = true;
    }
  }

  // Set the appropriate update function for taking a single time step
  // If we have implicit fluid-EM coupling, we perform a staggered (Strang)
  // operator split and treat that coupling implicitly.
  // Otherwise, we default to an SSP-RK3 method.
  if (app->has_fluid_em_coupling) {
    app->update_func = vlasov_update_op_split;
  } else {
    app->update_func = vlasov_update_ssp_rk3;
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat){
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = {DBL_MAX, 0.0},
    .stage_3_dt_diff = {DBL_MAX, 0.0},
  };

  return app;
}

int
vlasov_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (strcmp(nm, app->species[i].name) == 0) {
      return i;
    }
  }
  return -1;
}

struct vlasov_species *
vlasov_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  int i = vlasov_find_species_idx(app, nm);
  return i >= 0 ? &app->species[i] : 0;
}

struct vm_species *
vm_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  int i = vm_find_species_idx(app, nm);
  return i >= 0 ? app->species[i].kinetic : 0;
}

int
vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  int i = vlasov_find_species_idx(app, nm);
  return (i >= 0 && app->species[i].kinetic) ? i : -1;
}

void
vm_pos_c2p(const double *xcomp, double *xphys, void *ctx)
{
  gkyl_vlasov_position_map_eval_mc2p(ctx, xcomp, xphys);
}

void
vm_apply_bc(
  gkyl_vlasov_app *app, double tcurr, struct gkyl_array *distf[], struct gkyl_array *fluid[],
  struct gkyl_array *emfield
)
{
  // distf[] and fluid[] are indexed by species.
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_apply_bc(app, &app->species[i], distf[i], fluid[i], tcurr);
  }
  // No-op for Vlasov-Poisson and the null field (no EM boundary conditions).
  vlasov_field_apply_bc(app, emfield);
}

// Held backgrounds (First_Fourier_Harmonic_Notes/background_density_gradient.tex), captured at the
// end of the initial conditions (no-op without a held species). The held species are put on their
// backgrounds, with the other species at their initial conditions and the evolved Maxwell field zero
// (the background magnetic field is external). The right-hand side of a stage there, without the
// force of the backgrounds, gives the force that holds each held species (its momentum moment divided
// by the density); the right-hand side with that force is the remainder every stage subtracts. The
// initial conditions wait in the fnew buffers (the projection of a PKPM species uses f1), the field in
// em1.
static void
vlasov_capture_backgrounds(gkyl_vlasov_app *app, double t0)
{
  if (!app->has_background) {
    return;
  }
  int num_species = app->num_species;
  struct vm_field *field = app->field;
  bool has_em = field->field_id == GKYL_FIELD_E_B;

  for (int i = 0; i < num_species; ++i) {
    vlasov_species_copy_range(app, &app->species[i], VM_RK_FNEW, VM_RK_F);
  }
  if (has_em) {
    gkyl_array_copy(field->em1, field->em);
    gkyl_array_clear(field->em, 0.0);
    vlasov_field_apply_bc(app, field->em);
  }
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (vlasov_species_has_background(sp)) {
      vlasov_species_background_reset(app, sp);
      vlasov_species_background_project(app, sp, t0);
    }
  }

  struct gkyl_array *fin[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  vlasov_species_gather_rk_state(app, VM_RK_F, VM_RK_F1, fin, fout, fluidin, fluidout);

  // The force, from the right-hand side without it.
  vlasov_stage_rhs(app, 1, t0, fin, fluidin, field->em, fout, fluidout, field->emnew);
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (vlasov_species_has_background(sp)) {
      vlasov_species_background_calc_force(app, sp, fin[i], fluidin[i], fout[i], fluidout[i]);
    }
  }

  // The remainder, the right-hand side with the force. A PKPM species takes the speeds of its
  // velocity-domain edges from the stage before (in the first stage of a step from the last
  // completed step): the stage is repeated until those speeds are the background's own (the
  // dependence on them is through the traces at the edges, so two repeats reach round-off), and the
  // first stage of the run takes them too (complete_step keeps them).
  for (int pass = 0; pass < 3; ++pass) {
    vlasov_stage_rhs(
      app, pass == 0 ? 1 : 2, t0, fin, fluidin, field->em, fout, fluidout, field->emnew
    );
  }
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (vlasov_species_has_background(sp)) {
      vlasov_species_background_store(app, sp, fout[i], fluidout[i]);
    }
  }

  // The initial conditions back; nothing the evaluations recorded is kept, except the speeds of the
  // velocity-domain edges of a held PKPM species.
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    vlasov_species_discard_step(app, sp);
    if (vlasov_species_has_background(sp)) {
      vlasov_species_complete_step(app, sp);
    }
    vlasov_species_copy_range(app, sp, VM_RK_F, VM_RK_FNEW);
  }
  if (has_em) {
    gkyl_array_copy(field->em, field->em1);
  }
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app *app, double t0)
{
  app->tcurr = t0;
  // a new run of the staggered operator split
  app->op_split_dt_prev = 0.0;
  app->op_split_dt_suggested = 0.0;
  app->op_split_synced = false;
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_apply_ic(app, &app->species[i], t0);
  }

  gkyl_vlasov_app_apply_ic_field(app, t0); // no-op for the null field

  // BCs must be done after all species initialize for emission BCs to work.
  int num_species = app->num_species;
  struct gkyl_array *distf[num_species];
  struct gkyl_array *fluid[num_species];
  vlasov_species_gather_state(app, distf, fluid);
  vm_apply_bc(app, t0, distf, fluid, app->field->em);

  // The force and the remainder of the held backgrounds, with the other species at their initial
  // conditions (no-op without a held species).
  vlasov_capture_backgrounds(app, t0);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app *app, double t0)
{
  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();

  // Dispatches to the Maxwell or Poisson IC (no-op for the null field).
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_apply_ic(app, fin, t0);

  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app *app, int sidx, double t0)
{
  assert(sidx < app->num_species);
  vlasov_species_apply_ic(app, &app->species[sidx], t0);
}

int
gkyl_vlasov_app_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  return vlasov_find_species_idx(app, nm);
}

void
gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app *app, double tm)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_calc_integrated_mom(app, &app->species[i], tm);
  }
}

void
gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app *app, double tm)
{
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_calc_integrated_L2_f(app, &app->species[i], tm);
  }
}

void
gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app *app, double tm)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  // The distributions let Vlasov-Poisson solve for the potential at time tm.
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_calc_energy(app, tm, fin); // no-op for the null field
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app *app, double tm, int frame)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  gkyl_vlasov_app_write_field(app, tm, frame); // no-op for the null field
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_vlasov_app_write_species(app, i, tm, frame);
  }
}

void
gkyl_vlasov_app_write_field(gkyl_vlasov_app *app, double tm, int frame)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  // The distributions let Vlasov-Poisson solve for the potential at time tm.
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_write(app, tm, frame, fin);
  vlasov_op_split_write_field(app, tm, frame);
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app *app, int sidx, double tm, int frame)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  vlasov_species_write(app, &app->species[sidx], tm, frame);
  vlasov_op_split_write_species(app, &app->species[sidx], tm, frame);
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app *app, double tm, int frame)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_write_mom(app, &app->species[i], tm, frame);
  }
}

void
gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_write_integrated_mom(app, &app->species[i]);
  }
}

void
gkyl_vlasov_app_write_integrated_L2_f(gkyl_vlasov_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_write_integrated_L2_f(app, &app->species[i]);
  }
}

void
gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app *app)
{
  vlasov_field_write_energy(app); // no-op for the null field
}

void
gkyl_vlasov_app_write_lte_corr_status(gkyl_vlasov_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_write_lte_corr_status(app, &app->species[i]);
  }
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app *app, double dt)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu) {
    checkCuda(cudaGetLastError());
  }
  return status;
}

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app *app)
{
  vm_species_n_iter_corr(app);
  return app->stat;
}

static void
range_stat_write(gkyl_vlasov_app *app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_vlasov_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i = 0; i < r->ndim; ++i) {
    gkyl_vlasov_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  }
  gkyl_vlasov_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(
  const gkyl_vlasov_app *app, const struct gkyl_vlasov_stat *local, struct gkyl_vlasov_stat *global
)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_vlasov_stat));
    return;
  }

  global->use_gpu = local->use_gpu;

  enum { NUP, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, NFL_EM_RETAKE, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail,
    [NFL_EM_RETAKE] = local->nfl_em_retake,
  };

  int64_t l_red_global[L_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

  global->nup = l_red_global[NUP];
  global->nfeuler = l_red_global[NFEULER];
  global->nstage_2_fail = l_red_global[NSTAGE_2_FAIL];
  global->nstage_3_fail = l_red_global[NSTAGE_3_FAIL];
  global->nfl_em_retake = l_red_global[NFL_EM_RETAKE];

  int64_t l_red_n_iter_corr[app->num_species];
  int64_t l_red_num_corr[app->num_species];
  for (int s = 0; s < app->num_species; ++s) {
    l_red_n_iter_corr[s] = local->n_iter_corr[s];
    l_red_num_corr[s] = local->num_corr[s];
  }

  int64_t l_red_global_n_iter_corr[app->num_species];
  int64_t l_red_global_num_corr[app->num_species];
  gkyl_comm_allreduce_host(
    app->comm, GKYL_INT_64, GKYL_MAX, app->num_species, l_red_n_iter_corr, l_red_global_n_iter_corr
  );
  gkyl_comm_allreduce_host(
    app->comm, GKYL_INT_64, GKYL_MAX, app->num_species, l_red_num_corr, l_red_global_num_corr
  );

  for (int s = 0; s < app->num_species; ++s) {
    global->n_iter_corr[s] = l_red_global_n_iter_corr[s];
    global->num_corr[s] = l_red_global_num_corr[s];
  }

  enum {
    TOTAL_TM,
    RK3_TM,
    FL_EM_TM,
    INIT_SPECIES_TM,
    INIT_FLUID_SPECIES_TM,
    INIT_FIELD_TM,
    SPECIES_RHS_TM,
    FLUID_SPECIES_RHS_TM,
    FLUID_SPECIES_VARS_TM,
    SPECIES_PKPM_VARS_TM,
    FIELD_EM_VARS_TM,
    SPECIES_COLL_MOM_TM,
    SPECIES_COL_TM,
    SPECIES_RAD_TM,
    SPECIES_LTE_TM,
    FIELD_RHS_TM,
    CURRENT_TM,
    SPECIES_OMEGA_CFL_TM,
    FIELD_OMEGA_CFL_TM,
    SPECIES_BC_TM,
    FLUID_SPECIES_BC_TM,
    FIELD_BC_TM,
    SPECIES_IO_TM,
    SPECIES_DIAG_CALC_TM,
    SPECIES_DIAG_IO_TM,
    FIELD_IO_TM,
    FIELD_DIAG_CALC_TM,
    FIELD_DIAG_IO_TM,
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [RK3_TM] = local->rk3_tm,
    [FL_EM_TM] = local->fl_em_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FLUID_SPECIES_TM] = local->init_fluid_species_tm,
    [INIT_FIELD_TM] = local->init_field_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FLUID_SPECIES_RHS_TM] = local->fluid_species_rhs_tm,
    [FLUID_SPECIES_VARS_TM] = local->fluid_species_vars_tm,
    [SPECIES_PKPM_VARS_TM] = local->species_pkpm_vars_tm,
    [FIELD_EM_VARS_TM] = local->field_em_vars_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [SPECIES_RAD_TM] = local->species_rad_tm,
    [SPECIES_LTE_TM] = local->species_lte_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [CURRENT_TM] = local->current_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [FIELD_OMEGA_CFL_TM] = local->field_omega_cfl_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FLUID_SPECIES_BC_TM] = local->fluid_species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm,
    [SPECIES_IO_TM] = local->species_io_tm,
    [SPECIES_DIAG_CALC_TM] = local->species_diag_calc_tm,
    [SPECIES_DIAG_IO_TM] = local->species_diag_io_tm,
    [FIELD_IO_TM] = local->field_io_tm,
    [FIELD_DIAG_CALC_TM] = local->field_diag_calc_tm,
    [FIELD_DIAG_IO_TM] = local->field_diag_io_tm,
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);

  global->total_tm = d_red_global[TOTAL_TM];
  global->rk3_tm = d_red_global[RK3_TM];
  global->fl_em_tm = d_red_global[FL_EM_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_fluid_species_tm = d_red_global[INIT_FLUID_SPECIES_TM];
  global->init_field_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->fluid_species_rhs_tm = d_red_global[FLUID_SPECIES_RHS_TM];
  global->fluid_species_vars_tm = d_red_global[FLUID_SPECIES_VARS_TM];
  global->species_pkpm_vars_tm = d_red_global[SPECIES_PKPM_VARS_TM];
  global->field_em_vars_tm = d_red_global[FIELD_EM_VARS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_rad_tm = d_red_global[SPECIES_RAD_TM];
  global->species_lte_tm = d_red_global[SPECIES_LTE_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->current_tm = d_red_global[CURRENT_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->field_omega_cfl_tm = d_red_global[FIELD_OMEGA_CFL_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->fluid_species_bc_tm = d_red_global[FLUID_SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];
  global->species_io_tm = d_red_global[SPECIES_IO_TM];
  global->species_diag_calc_tm = d_red_global[SPECIES_DIAG_CALC_TM];
  global->species_diag_io_tm = d_red_global[SPECIES_DIAG_IO_TM];
  global->field_io_tm = d_red_global[FIELD_IO_TM];
  global->field_diag_calc_tm = d_red_global[FIELD_DIAG_CALC_TM];
  global->field_diag_io_tm = d_red_global[FIELD_DIAG_IO_TM];

  // misc data needing reduction

  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff, global->stage_2_dt_diff
  );
  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff, global->stage_3_dt_diff
  );

  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm
  );
  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
    global->species_lbo_coll_diff_tm
  );
}

void
gkyl_vlasov_app_stat_write(gkyl_vlasov_app *app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  int num_ranks;
  gkyl_comm_get_size(app->comm, &num_ranks);

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  vm_species_n_iter_corr(app);

  struct gkyl_vlasov_stat stat = {};
  comm_reduce_app_stat(app, &app->stat, &stat);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) {
    fp = fopen(fileNm, "a");
  }

  gkyl_vlasov_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm)) {
    gkyl_vlasov_app_cout(app, fp, " date : %s,\n", buff);
  }

  gkyl_vlasov_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_vlasov_app_cout(app, fp, " num_ranks : %d,\n", num_ranks);

  for (int s = 0; s < app->num_species; ++s) {
    if (app->species[s].kinetic) {
      range_stat_write(app, app->species[s].kinetic->name, &app->species[s].kinetic->global, fp);
    }
  }

  gkyl_vlasov_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_vlasov_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_vlasov_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, fp, " nfluid_em_coupling_retake : %ld,\n", stat.nfl_em_retake);

  gkyl_vlasov_app_cout(
    app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n", stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]
  );
  gkyl_vlasov_app_cout(
    app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n", stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]
  );

  gkyl_vlasov_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_vlasov_app_cout(app, fp, " rk3_tm : %lg,\n", stat.rk3_tm);
  gkyl_vlasov_app_cout(app, fp, " fluid_em_coupling_tm : %lg,\n", stat.fl_em_tm);
  gkyl_vlasov_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  if (app->has_field) {
    gkyl_vlasov_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);
  }

  gkyl_vlasov_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);

  for (int s = 0; s < app->num_species; ++s) {
    gkyl_vlasov_app_cout(app, fp, " n_iter_corr[%d] : %ld,\n", s, stat.n_iter_corr[s]);
    gkyl_vlasov_app_cout(app, fp, " num_corr[%d] : %ld,\n", s, stat.num_corr[s]);
    gkyl_vlasov_app_cout(
      app, fp, " species_coll_drag_tm[%d] : %lg,\n", s, stat.species_lbo_coll_drag_tm[s]
    );
    gkyl_vlasov_app_cout(
      app, fp, " species_coll_diff_tm[%d] : %lg,\n", s, stat.species_lbo_coll_diff_tm[s]
    );
  }

  gkyl_vlasov_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, fp, " species_rad_tm : %lg,\n", stat.species_rad_tm);

  gkyl_vlasov_app_cout(app, fp, " species_lte_tm : %lg,\n", stat.species_lte_tm);

  gkyl_vlasov_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);

  gkyl_vlasov_app_cout(app, fp, " fluid_species_rhs_tm : %lg,\n", stat.fluid_species_rhs_tm);

  gkyl_vlasov_app_cout(app, fp, " fluid_species_bc_tm : %lg,\n", stat.fluid_species_bc_tm);

  if (app->has_field) {
    gkyl_vlasov_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
    gkyl_vlasov_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);

    gkyl_vlasov_app_cout(app, fp, " current_tm : %lg,\n", stat.current_tm);
  }

  gkyl_vlasov_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);
  gkyl_vlasov_app_cout(app, fp, " field_omega_cfl_tm : %lg\n", stat.field_omega_cfl_tm);
  gkyl_vlasov_app_cout(app, fp, " species_io_tm : %lg\n", stat.species_io_tm);
  gkyl_vlasov_app_cout(app, fp, " species_diag_tm : %lg\n", stat.species_diag_calc_tm);
  gkyl_vlasov_app_cout(app, fp, " species_diag_io_tm : %lg\n", stat.species_diag_io_tm);
  gkyl_vlasov_app_cout(app, fp, " field_io_tm : %lg\n", stat.field_io_tm);
  gkyl_vlasov_app_cout(app, fp, " field_diag_tm : %lg\n", stat.field_diag_calc_tm);
  gkyl_vlasov_app_cout(app, fp, " field_diag_io_tm : %lg\n", stat.field_diag_io_tm);

  gkyl_vlasov_app_cout(app, fp, " n_species_omega_cfl : %ld,\n", stat.n_species_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " n_field_omega_cfl : %ld,\n", stat.n_field_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " n_mom : %ld,\n", stat.n_mom);
  gkyl_vlasov_app_cout(app, fp, " n_diag : %ld,\n", stat.n_diag);
  gkyl_vlasov_app_cout(app, fp, " n_io : %ld,\n", stat.n_io);
  gkyl_vlasov_app_cout(app, fp, " n_diag_io : %ld,\n", stat.n_diag_io);
  gkyl_vlasov_app_cout(app, fp, " n_field_diag : %ld,\n", stat.n_field_diag);
  gkyl_vlasov_app_cout(app, fp, " n_field_io : %ld,\n", stat.n_field_io);
  gkyl_vlasov_app_cout(app, fp, " n_field_diag_io : %ld,\n", stat.n_field_diag_io);

  if (app->has_pkpm) {
    // PKPM species: the configuration-space variables of the stages and b of the field, and the
    // record of the closure solves of each species.
    gkyl_vlasov_app_cout(app, fp, " species_pkpm_vars_tm : %lg,\n", stat.species_pkpm_vars_tm);
    gkyl_vlasov_app_cout(app, fp, " field_em_vars_tm : %lg,\n", stat.field_em_vars_tm);
    for (int i = 0; i < app->num_species; ++i) {
      if (app->species[i].type == GKYL_SPECIES_PKPM) {
        vlasov_pkpm_species_stat_write(app, &app->species[i], fp);
      }
    }
  }

  gkyl_vlasov_app_cout(app, fp, "}\n");

  if (rank == 0) {
    fclose(fp);
  }
}

struct gkyl_app_restart_status
vlasov_header_from_file(gkyl_vlasov_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = {.io_status = 0};

  FILE *fp = 0;
  with_file(fp, fname, "r")
  {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE) {
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
      }
    }

    struct vlasov_output_meta meta = vlasov_meta_from_mpack(&(struct gkyl_msgpack_data){
      .meta = hdr.meta,
      .meta_sz = hdr.meta_size,
    });

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }

  return rstat;
}

struct vlasov_output_meta
vlasov_meta_from_file(const char *fname)
{
  struct vlasov_output_meta meta = {.frame = 0, .stime = 0.0, .dt_prev = 0.0, .dt_suggested = 0.0};

  FILE *fp = 0;
  with_file(fp, fname, "r")
  {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    if (GKYL_ARRAY_RIO_SUCCESS == gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp)) {
      meta = vlasov_meta_from_mpack(&(struct gkyl_msgpack_data){
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size,
      });
    }
    gkyl_grid_sub_array_header_release(&hdr);
  }

  return meta;
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_file_field(gkyl_vlasov_app *app, const char *fname)
{
  return vlasov_field_from_file(app, fname);
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_file_species(gkyl_vlasov_app *app, int sidx, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_species_from_file(app, &app->species[sidx], fname);
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    vlasov_species_restart_bflux(app, &app->species[sidx]);
    vlasov_species_restart_bc(app, &app->species[sidx]);
  }
  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_field(gkyl_vlasov_app *app, int frame)
{
  return vlasov_field_read_from_frame(app, frame);
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_species(gkyl_vlasov_app *app, int sidx, int frame)
{
  struct gkyl_app_restart_status rstat =
    vlasov_species_read_from_frame(app, &app->species[sidx], frame);
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    vlasov_species_restart_bflux(app, &app->species[sidx]);
    vlasov_species_restart_bc(app, &app->species[sidx]);
  }
  return rstat;
}

struct gkyl_app_restart_status
gkyl_vlasov_app_read_from_frame(gkyl_vlasov_app *app, int frame)
{
  // Species are read first: they set the app clock from their files, and the
  // field then rebuilds its time-dependent inputs at that time (the Poisson
  // field has no file of its own and evaluates its external potentials at the
  // app clock; the potential itself is solved from the distribution whenever it
  // is needed). The first failure is reported; otherwise the status carries the
  // frame and time read from the species files, or from the field file when
  // there are no species.
  // A static kinetic species has its distribution in the file of frame 0 only:
  // it is read after the species and the field that evolve, and the frame and
  // the clock are theirs. When nothing else has a file, the run restarts at the
  // time of frame 0.
  // With PKPM species the Maxwell field (which has a file and a time of its own) is read first:
  // their boundary conditions read b of the restart field.
  // Held backgrounds are not in the frames: they are captured again from the initial conditions
  // and the backgrounds, before the frame is read over them.
  if (app->has_background) {
    gkyl_vlasov_app_apply_ic(app, 0.0);
  }
  struct gkyl_app_restart_status rstat = {.io_status = GKYL_ARRAY_RIO_SUCCESS};
  struct gkyl_app_restart_status fstat = {.io_status = GKYL_ARRAY_RIO_SUCCESS};
  if (app->has_pkpm) {
    fstat = gkyl_vlasov_app_from_frame_field(app, frame);
    if (fstat.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      return fstat;
    }
  }
  int num_evolving = 0;
  for (int i = 0; i < app->num_species; i++) {
    struct vlasov_species *sp = &app->species[i];
    if (sp->kinetic && sp->kinetic->info.is_static) {
      continue;
    }
    rstat = vlasov_species_read_from_frame(app, sp, frame);
    if (rstat.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      return rstat;
    }
    num_evolving += 1;
  }

  if (!app->has_pkpm) {
    fstat = gkyl_vlasov_app_from_frame_field(app, frame);
  }
  if (fstat.io_status != GKYL_ARRAY_RIO_SUCCESS || app->num_species == 0) {
    return fstat;
  }
  bool has_field_file = app->field->field_id == GKYL_FIELD_E_B ||
                        app->field->field_id == GKYL_FIELD_GR_D_B;
  if (num_evolving == 0 && has_field_file) {
    rstat = fstat;
  }

  double tm_frame = app->tcurr;
  for (int i = 0; i < app->num_species; i++) {
    struct vlasov_species *sp = &app->species[i];
    if (sp->kinetic && sp->kinetic->info.is_static) {
      struct gkyl_app_restart_status sstat = vlasov_species_read_from_frame(app, sp, frame);
      if (sstat.io_status != GKYL_ARRAY_RIO_SUCCESS) {
        return sstat;
      }
      if (num_evolving == 0 && !has_field_file) {
        rstat = sstat;
      } else {
        app->tcurr = tm_frame;
      }
    }
  }

  // Rebuild the ghost cells once every species is read, as the initial conditions
  // do: first all the boundary fluxes, since an emitting wall reads those of its
  // impact species (which may come later in the species list), then the BCs.
  for (int i = 0; i < app->num_species; i++) {
    vlasov_species_restart_bflux(app, &app->species[i]);
  }
  for (int i = 0; i < app->num_species; i++) {
    vlasov_species_restart_bc(app, &app->species[i]);
  }

  // The frames hold the synchronized state of the staggered operator split.
  struct gkyl_app_restart_status sstat = vlasov_op_split_restart(app, frame);
  if (sstat.io_status != GKYL_ARRAY_RIO_SUCCESS) {
    return sstat;
  }

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_vlasov_app_cout(const gkyl_vlasov_app *app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_vlasov_app_cout(const gkyl_vlasov_app *app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_vlasov_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app *app)
{
  vm_geom_release(app, app->vm_geom);
  gkyl_free(app->vm_geom);
  for (int i = 0; i < app->num_species; ++i) {
    vlasov_species_release(app, &app->species[i]);
  }
  if (app->species) {
    gkyl_free(app->species);
  }
  // A field object always exists (the null field when no field is present).
  vlasov_field_release(app);
  if (app->has_fluid_em_coupling) {
    vm_fluid_em_coupling_release(app, app->fl_em);
  }
  vm_coll_exchange_release(app);
  if (app->conf_nodal) {
    gkyl_dg_conf_nodal_release(app->conf_nodal);
  }
  if (app->src_omega_p_sq) {
    gkyl_array_release(app->src_omega_p_sq);
  }
  if (app->src_dt_func == vlasov_src_dt_explicit) {
    gkyl_array_release(app->src_omega_p_sq_kin);
    gkyl_array_release(app->src_rho_avg);
    if (app->use_gpu) {
      gkyl_cu_free(app->src_omega_p_sq_max);
    } else {
      gkyl_free(app->src_omega_p_sq_max);
    }
  }

  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);

  gkyl_wave_geom_release(app->geom);
  gkyl_vlasov_position_map_release(app->pos_map);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev);
  }

  gkyl_free(app);
}
