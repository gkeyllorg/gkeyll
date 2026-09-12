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
// This REPORTS by default and does not change the verdict, because existing
// multiblock cases carry mixed interfaces today and would begin failing at
// setup.  Set GKYL_TOK_STRICT_SEAM_PARTICIPATION=1 to make it an error.
static int
gk_block_geom_check_seam_participation(const struct gkyl_gk_block_geom *bgeom,
  int *interfaces_examined, int *mixed)
{
  const char *strict_env = getenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  bool strict = strict_env && strict_env[0] != '\0' && strict_env[0] != '0';
  const char *diag_env = getenv("GKYL_TOK_SEAM_PARTICIPATION_DIAG");
  bool diag = diag_env && diag_env[0] != '\0' && diag_env[0] != '0';
  int nmixed = 0, nchecked = 0;

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

        if (nmixed == 0)
          fprintf(stderr,
            "TOK_SEAM_PARTICIPATION mixed extended construction across a declared "
            "interface: the blocks may parameterize their shared row "
            "differently. Uniform participation is a declaration check; "
            "constructed seam alignment and interior ordering still need "
            "validation.\n");

        fprintf(stderr,
          "TOK_SEAM_PARTICIPATION block=%d ftype=%d extended=%d <-> block=%d "
          "ftype=%d extended=%d dir=%d edge=%s\n",
          i, (int) bi->geometry.tok_grid_info.ftype, (int) ext_i,
          j, (int) bj->geometry.tok_grid_info.ftype, (int) ext_j,
          d, e ? "upper" : "lower");
        nmixed += 1;
      }
    }
  }

  if (diag)
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION_DIAG interfaces_examined=%d mixed=%d\n",
      nchecked, nmixed);

  *interfaces_examined = nchecked;
  *mixed = nmixed;

  if (nmixed > 0 && strict) {
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION %d mixed interface(s); failing because "
      "GKYL_TOK_STRICT_SEAM_PARTICIPATION is set\n", nmixed);
    return 0;
  }
  return 1;
}

static int
gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom,
  int *interfaces_examined, int *mixed)
{
  *interfaces_examined = 0;
  *mixed = 0;
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
    interfaces_examined, mixed);
  return topo_ok && seam_ok;
}

int
gkyl_gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom)
{
  int interfaces_examined, mixed;
  return gk_block_geom_check_consistency(bgeom, &interfaces_examined, &mixed);
}

int
gkyl_gyrokinetic_multib_app_geometry_preflight(const struct gkyl_gyrokinetic_multib *mbinp)
{
  const struct gkyl_gk_block_geom *bgeom = mbinp ? mbinp->gk_block_geom : 0;
  int interfaces_examined, mixed;
  int ok = gk_block_geom_check_consistency(bgeom, &interfaces_examined, &mixed);
  if (mbinp && bgeom && mbinp->cdim != bgeom->ndim) {
    fprintf(stderr, "GKYL_BLOCK_GEOMETRY_INVALID app_cdim=%d geometry_ndim=%d\n",
      mbinp->cdim, bgeom->ndim);
    ok = 0;
  }
  const char *strict_env = getenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  bool strict = strict_env && strict_env[0] != '\0' && strict_env[0] != '0';
  fprintf(stderr,
    "GKYL_GEOMETRY_PREFLIGHT status=%s scope=declaration num_blocks=%d "
    "strict=%d interfaces_examined=%d mixed=%d\n",
    ok ? "PASS" : "FAIL", bgeom ? bgeom->num_blocks : 0,
    (int) strict, interfaces_examined, mixed);
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
