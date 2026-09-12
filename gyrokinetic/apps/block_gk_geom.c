#include <gkyl_gk_block_geom.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_alloc.h>
#include <gkyl_tok_geo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Geometry info for all blocks in simulation
struct gkyl_gk_block_geom {
  int ndim; // dimension
  int num_blocks; // total number of blocks
  struct gkyl_gk_block_geom_info *blocks; // info for each block
  struct gkyl_block_topo *btopo; // topology of blocks
  
  struct gkyl_ref_count ref_count;
};

static void
gk_block_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_gk_block_geom *bgeom = container_of(ref, struct gkyl_gk_block_geom, ref_count);
  gkyl_free(bgeom->blocks);
  gkyl_block_topo_release(bgeom->btopo);
  gkyl_free(bgeom);
}

struct gkyl_gk_block_geom*
gkyl_gk_block_geom_new(int ndim, int nblocks)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_malloc(sizeof(struct gkyl_gk_block_geom));
  bgeom->ndim = ndim;
  bgeom->num_blocks = nblocks;
  bgeom->blocks = gkyl_calloc(sizeof(struct gkyl_gk_block_geom_info), nblocks);

  bgeom->btopo = gkyl_block_topo_new(ndim, nblocks);

  bgeom->ref_count = gkyl_ref_count_init(gk_block_geom_free);

  return bgeom;
}

int
gkyl_gk_block_geom_ndim(const struct gkyl_gk_block_geom *bgeom)
{
  return bgeom->ndim;
}

int
gkyl_gk_block_geom_num_blocks(const struct gkyl_gk_block_geom *bgeom)
{
  return bgeom->num_blocks;
}

void
gkyl_gk_block_geom_set_block(struct gkyl_gk_block_geom *bgeom, int bidx,
  const struct gkyl_gk_block_geom_info *info)
{
  memcpy(&bgeom->blocks[bidx], info, sizeof(struct gkyl_gk_block_geom_info));
  
  for (int d=0; d<bgeom->ndim; ++d)
    bgeom->blocks[bidx].cuts[d] = info->cuts[d] > 0 ? info->cuts[d] : 1;
  
  // set topology information
  for (int i=0; i<bgeom->ndim; ++i)
    for (int e=0; e<2; ++e)
      bgeom->btopo->conn[bidx].connections[i][e] = info->connections[i][e];  
}

void
gkyl_gk_block_geom_reset_block_extents(struct gkyl_gk_block_geom *bgeom, int bidx, double *lower, double *upper)
{
  struct gkyl_gk_block_geom_info *bgi = &bgeom->blocks[bidx];
  for (int i = 0; i < bgeom->ndim; ++i) {
    bgi->lower[i] = lower[i];
    bgi->upper[i] = upper[i];
  }
}

void
gkyl_gk_block_geom_apply_xpt_seam_selection(struct gkyl_gk_block_geom *bgeom,
  int bidx, double coefficient, double bound)
{
  struct gkyl_gk_block_geom_info *bgi = &bgeom->blocks[bidx];
  struct gkyl_tok_geo_grid_inp *inp = &bgi->geometry.tok_grid_info;
  inp->relaxed_xpt_seam = true;
  inp->relaxed_xpt_seam_sweep = true;
  inp->relaxed_xpt_seam_delta_s_coeff = coefficient;
  inp->relaxed_xpt_seam_delta_s_bound = bound;
}

const struct gkyl_gk_block_geom_info*
gkyl_gk_block_geom_get_block(const struct gkyl_gk_block_geom *bgeom, int bidx)
{
  return &bgeom->blocks[bidx];
}

enum gkyl_gk_shared_sep_row_status
gkyl_gk_block_geom_shared_sep_row_status(const struct gkyl_gk_block_geom_info *legacy,
  const struct gkyl_gk_block_geom_info *peer, int src_dir, int tgt_dir)
{
  const char *shared_theta = getenv("GKYL_TOK_SHARED_SEP_THETA");
  if (!(shared_theta && shared_theta[0] && shared_theta[0] != '0'))
    return GKYL_GK_SHARED_SEP_ROW_NONE;

  if (legacy->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK ||
      peer->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
    return GKYL_GK_SHARED_SEP_ROW_NONE;

  // Only a block that does NOT take the extended construction can adopt a
  // peer's row, and only from one that does.  Anything else -- both extended,
  // both legacy, or the pair the wrong way round -- is a different question.
  if (gkyl_tok_geo_uses_extended_construction(&legacy->geometry.tok_grid_info) ||
      !gkyl_tok_geo_uses_extended_construction(&peer->geometry.tok_grid_info))
    return GKYL_GK_SHARED_SEP_ROW_NONE;

  // The shared row is the separatrix row, which two blocks share only across a
  // RADIAL interface.  A theta interface joins two different rows, so nothing
  // there can be taken from one trace builder.
  if (src_dir != 0)
    return GKYL_GK_SHARED_SEP_ROW_NONE;
  if (tgt_dir != 0)
    return GKYL_GK_SHARED_SEP_ROW_UNSUPPORTED;

  // The peer's row is only reusable if it is the SAME curve: same equilibrium,
  // same DG representation of it, and no X-point bounding polygon narrowing
  // either block's critical-point search.
  const struct gkyl_efit_inp *a = &legacy->geometry.efit_info;
  const struct gkyl_efit_inp *b = &peer->geometry.efit_info;
  if (strcmp(a->filepath, b->filepath) || a->reflect != b->reflect ||
      a->rz_poly_order != b->rz_poly_order || a->flux_poly_order != b->flux_poly_order ||
      a->xpt_bound_n != 0 || b->xpt_bound_n != 0 ||
      legacy->geometry.tok_grid_info.use_cubics != peer->geometry.tok_grid_info.use_cubics)
    return GKYL_GK_SHARED_SEP_ROW_UNSUPPORTED;

  return GKYL_GK_SHARED_SEP_ROW_SHARED;
}

// Whether the shared-separatrix-row construction covers a mixed interface,
// asked the way the multiblock app wires it: from the legacy block's own
// radial connections.  A block adopts at most ONE peer row, so if both of its
// radial edges qualify the app refuses; report that as not covered rather than
// guessing which one wins.
static bool
gk_block_geom_seam_row_is_shared(const struct gkyl_gk_block_geom *bgeom,
  int legacy_bid, int peer_bid)
{
  const struct gkyl_gk_block_geom_info *legacy = &bgeom->blocks[legacy_bid];
  int nshared = 0, nedge = 0;
  for (int e=0; e<2; ++e) {
    const struct gkyl_target_edge *te = &bgeom->btopo->conn[legacy_bid].connections[0][e];
    if (te->edge == GKYL_PHYSICAL || te->bid == legacy_bid)
      continue;
    if (te->bid < 0 || te->bid >= bgeom->num_blocks)
      continue; // malformed; the topology check reports this itself
    enum gkyl_gk_shared_sep_row_status st =
      gkyl_gk_block_geom_shared_sep_row_status(legacy,
        &bgeom->blocks[te->bid], 0, te->dir);
    if (st == GKYL_GK_SHARED_SEP_ROW_NONE)
      continue;
    nedge += 1;
    if (te->bid == peer_bid && st == GKYL_GK_SHARED_SEP_ROW_SHARED)
      nshared += 1;
  }
  return nshared == 1 && nedge == 1;
}

// Report interfaces whose two blocks disagree about taking the extended
// construction.
//
// The extended path reparameterizes a block's separatrix row -- the theta
// ladder, and the |grad psi| poloidal measure -- while the legacy path leaves
// it on plain arc length.  Two blocks that SHARE that row and disagree about
// taking it therefore trace the same curve and still place different nodes
// along it.  Measured on asdex/tcv, which declare straight_xpt_ray on 3 of 6
// blocks: the two sides agree on arc length to 7e-05 relative and on both
// endpoints to 1e-14 m, and differ only in the interior, by 35-66 mm -- which
// is 0.6-0.8 cells at x1 and DOUBLES at every refinement, because it is a fixed
// offset rather than a convergence error.
//
// This detects one declaration-level risk without constructing geometry.
// Uniform participation alone does not establish interface alignment or
// interior ordering; those still require checks on the constructed grid.
//
// Note the participation predicate deliberately comes from
// gkyl_tok_geo_uses_extended_construction() rather than being restated: it also
// depends on half_domain and on an environment override, and a second copy
// would drift.
//
// The offset above is what happens when the two sides build the row
// SEPARATELY.  Where the shared-separatrix-row construction applies, the
// legacy block takes the row from its extended peer's trace builder, and then
// no declaration can make the two disagree: measured A/B at x1 with a 0.01
// cell tolerance, asdex goes 0.83 -> 7.7e-05 cells and tcv 2.04 -> 3.1e-03.
// So the interface is reported either way, but only an interface whose row is
// genuinely built twice is a defect -- a deliberate mixed declaration with the
// shared row on is a supported configuration, not a latent one.
//
// This REPORTS by default and does not change the verdict, because existing
// multiblock cases carry mixed interfaces today and would begin failing at
// setup.  Set GKYL_TOK_STRICT_SEAM_PARTICIPATION=1 to make an unshared one an
// error.
static int
gk_block_geom_check_seam_participation(const struct gkyl_gk_block_geom *bgeom,
  int *interfaces_examined, int *mixed, int *unshared)
{
  const char *strict_env = getenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  bool strict = strict_env && strict_env[0] != '\0' && strict_env[0] != '0';
  const char *diag_env = getenv("GKYL_TOK_SEAM_PARTICIPATION_DIAG");
  bool diag = diag_env && diag_env[0] != '\0' && diag_env[0] != '0';
  int nmixed = 0, nunshared = 0, nchecked = 0;

  // "No mixed interface" and "the check never ran" look identical from the
  // outside, and that cost a debugging cycle when this was first wired up.
  // The diagnostic distinguishes them by naming every block and every edge it
  // actually examined.
  if (diag) {
    fprintf(stderr, "TOK_SEAM_PARTICIPATION_DIAG num_blocks=%d ndim=%d\n",
      bgeom->num_blocks, bgeom->ndim);
    for (int i=0; i<bgeom->num_blocks; ++i) {
      const struct gkyl_gk_block_geom_info *bi = &bgeom->blocks[i];
      fprintf(stderr,
        "TOK_SEAM_PARTICIPATION_DIAG block=%d geometry_id=%d ftype=%d "
        "straight_xpt_ray=%d half_domain=%d extended=%d\n",
        i, (int) bi->geometry.geometry_id,
        (int) bi->geometry.tok_grid_info.ftype,
        (int) bi->geometry.tok_grid_info.straight_xpt_ray,
        (int) bi->geometry.tok_grid_info.half_domain,
        (int) gkyl_tok_geo_uses_extended_construction(&bi->geometry.tok_grid_info));
      for (int d=0; d<bgeom->ndim; ++d)
        for (int e=0; e<2; ++e) {
          const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
          fprintf(stderr,
            "TOK_SEAM_PARTICIPATION_DIAG   block=%d dir=%d edge=%s -> bid=%d edgekind=%d\n",
            i, d, e ? "upper" : "lower", te->bid, (int) te->edge);
        }
    }
  }

  for (int i=0; i<bgeom->num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bi = &bgeom->blocks[i];
    if (bi->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
      continue;

    for (int d=0; d<bgeom->ndim; ++d) {
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
        if (te->edge == GKYL_PHYSICAL)
          continue;
        int j = te->bid;
        // Visit each interface once, and skip a block joined to itself: a
        // self-periodic block cannot disagree with itself.
        if (j <= i)
          continue;
        if (j < 0 || j >= bgeom->num_blocks)
          continue; // malformed; the topology check reports this itself

        const struct gkyl_gk_block_geom_info *bj = &bgeom->blocks[j];
        if (bj->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
          continue;

        nchecked += 1;
        bool ext_i = gkyl_tok_geo_uses_extended_construction(&bi->geometry.tok_grid_info);
        bool ext_j = gkyl_tok_geo_uses_extended_construction(&bj->geometry.tok_grid_info);
        if (ext_i == ext_j)
          continue;

        // Mixed participation only misparameterizes the shared row if the two
        // blocks actually build it separately.  When the shared-separatrix-row
        // construction covers this interface the row comes from ONE trace
        // builder, so the disagreement has nothing left to act on -- measured
        // on asdex/tcv, where turning that construction off moves the seam
        // from 7.7e-05/3.1e-03 cells to 0.83/2.04 cells against a 0.01
        // tolerance.  Report the interface either way; fail only on the ones
        // that are genuinely built twice.
        bool shared_row = gk_block_geom_seam_row_is_shared(bgeom,
          ext_i ? j : i, ext_i ? i : j);

        if (nmixed == 0)
          fprintf(stderr,
            "TOK_SEAM_PARTICIPATION mixed extended construction across a declared "
            "interface: the blocks may parameterize their shared row "
            "differently. Uniform participation is a declaration check; "
            "constructed seam alignment and interior ordering still need "
            "validation.\n");

        fprintf(stderr,
          "TOK_SEAM_PARTICIPATION block=%d ftype=%d extended=%d <-> block=%d "
          "ftype=%d extended=%d dir=%d edge=%s shared_row=%d\n",
          i, (int) bi->geometry.tok_grid_info.ftype, (int) ext_i,
          j, (int) bj->geometry.tok_grid_info.ftype, (int) ext_j,
          d, e ? "upper" : "lower", (int) shared_row);
        nmixed += 1;
        if (!shared_row)
          nunshared += 1;
      }
    }
  }

  if (diag)
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION_DIAG interfaces_examined=%d mixed=%d unshared=%d\n",
      nchecked, nmixed, nunshared);

  *interfaces_examined = nchecked;
  *mixed = nmixed;
  *unshared = nunshared;

  if (nmixed > nunshared)
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION %d of %d mixed interface(s) take their shared row "
      "from one trace builder; those cannot be misparameterized by the "
      "declaration and do not fail\n", nmixed - nunshared, nmixed);

  if (nunshared > 0 && strict) {
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION %d mixed interface(s) build their shared row "
      "twice; failing because GKYL_TOK_STRICT_SEAM_PARTICIPATION is set\n",
      nunshared);
    return 0;
  }
  return 1;
}

static int
gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom,
  int *interfaces_examined, int *mixed, int *unshared)
{
  *interfaces_examined = 0;
  *mixed = 0;
  *unshared = 0;
  if (!bgeom || bgeom->ndim < 1 || bgeom->ndim > GKYL_MAX_CDIM ||
      bgeom->num_blocks < 1 || !bgeom->blocks || !bgeom->btopo ||
      !bgeom->btopo->conn || bgeom->btopo->ndim != bgeom->ndim ||
      bgeom->btopo->num_blocks != bgeom->num_blocks) {
    fprintf(stderr, "GKYL_BLOCK_GEOMETRY_INVALID missing or invalid block geometry\n");
    return 0;
  }

  // Validate target indices before following them. The gyrokinetic range
  // transforms still require aligned logical directions.
  int topo_ok = 1;
  for (int i=0; i<bgeom->num_blocks; ++i)
    for (int d=0; d<bgeom->ndim; ++d)
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
        if (te->edge < GKYL_LOWER_POSITIVE || te->edge > GKYL_PHYSICAL ||
            (te->edge != GKYL_PHYSICAL &&
             (te->bid < 0 || te->bid >= bgeom->num_blocks ||
              te->dir != d))) {
          fprintf(stderr,
            "GKYL_BLOCK_GEOMETRY_INVALID block=%d dir=%d edge=%d "
            "target_block=%d target_dir=%d target_edge=%d\n",
            i, d, e, te->bid, te->dir, (int) te->edge);
          topo_ok = 0;
        }
      }
  // Follow the declared target direction/edge when checking reciprocity.
  // Same-side radial edges are valid: TCV CORE and LSN_SOL_MID both use their
  // upper radial edge at the separatrix. Follow the target edge rather than
  // assuming it is the opposite of the source edge.
  if (topo_ok)
    for (int i=0; i<bgeom->num_blocks; ++i)
      for (int d=0; d<bgeom->ndim; ++d)
        for (int e=0; e<2; ++e) {
          const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
          if (te->edge == GKYL_PHYSICAL)
            continue;
          int target_edge = te->edge <= GKYL_LOWER_NEGATIVE ? 0 : 1;
          bool negative = te->edge == GKYL_LOWER_NEGATIVE || te->edge == GKYL_UPPER_NEGATIVE;
          int expected_edge = (e ? GKYL_UPPER_POSITIVE : GKYL_LOWER_POSITIVE) + negative;
          const struct gkyl_target_edge *back =
            &bgeom->btopo->conn[te->bid].connections[te->dir][target_edge];
          if (back->bid != i || back->dir != d || back->edge != expected_edge) {
            fprintf(stderr,
              "GKYL_BLOCK_GEOMETRY_INVALID nonreciprocal block=%d dir=%d edge=%d "
              "target_block=%d target_dir=%d target_edge=%d\n",
              i, d, e, te->bid, te->dir, (int) te->edge);
            topo_ok = 0;
          }
        }
  if (!topo_ok)
    fprintf(stderr, "GKYL_BLOCK_GEOMETRY_INVALID inconsistent topology\n");

  // Participation diagnostics remain useful even when topology is invalid.
  int seam_ok = gk_block_geom_check_seam_participation(bgeom,
    interfaces_examined, mixed, unshared);
  return topo_ok && seam_ok;
}

int
gkyl_gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom)
{
  int interfaces_examined, mixed, unshared;
  return gk_block_geom_check_consistency(bgeom, &interfaces_examined, &mixed, &unshared);
}

int
gkyl_gyrokinetic_multib_app_geometry_preflight(const struct gkyl_gyrokinetic_multib *mbinp)
{
  const struct gkyl_gk_block_geom *bgeom = mbinp ? mbinp->gk_block_geom : 0;
  int interfaces_examined, mixed, unshared;
  int ok = gk_block_geom_check_consistency(bgeom, &interfaces_examined, &mixed, &unshared);
  if (mbinp && bgeom && mbinp->cdim != bgeom->ndim) {
    fprintf(stderr, "GKYL_BLOCK_GEOMETRY_INVALID app_cdim=%d geometry_ndim=%d\n",
      mbinp->cdim, bgeom->ndim);
    ok = 0;
  }
  const char *strict_env = getenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  bool strict = strict_env && strict_env[0] != '\0' && strict_env[0] != '0';
  fprintf(stderr,
    "GKYL_GEOMETRY_PREFLIGHT status=%s scope=declaration num_blocks=%d "
    "strict=%d interfaces_examined=%d mixed=%d unshared=%d\n",
    ok ? "PASS" : "FAIL", bgeom ? bgeom->num_blocks : 0,
    (int) strict, interfaces_examined, mixed, unshared);
  fflush(stderr);
  return ok;
}

struct gkyl_gk_block_geom *
gkyl_gk_block_geom_acquire(const struct gkyl_gk_block_geom* bgeom)
{
  gkyl_ref_count_inc(&bgeom->ref_count);
  return (struct gkyl_gk_block_geom*) bgeom;
}

struct gkyl_block_topo*
gkyl_gk_block_geom_topo(const struct gkyl_gk_block_geom *bgeom)
{
  return gkyl_block_topo_acquire(bgeom->btopo);
}

void
gkyl_gk_block_geom_release(struct gkyl_gk_block_geom* bgeom)
{
  gkyl_ref_count_dec(&bgeom->ref_count);
}
