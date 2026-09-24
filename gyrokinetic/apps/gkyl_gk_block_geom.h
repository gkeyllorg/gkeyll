#pragma once

#include <gkyl_block_topo.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic.h>

// Geometry info a single config-space block
struct gkyl_gk_block_geom_info {
  // lower and upper extents of blocks
  double lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  int cells[GKYL_MAX_CDIM]; // cells extents in each direction
  int cuts[GKYL_MAX_CDIM];  // domain split to use

  struct gkyl_gyrokinetic_geometry geometry; // GK geometry

  struct gkyl_target_edge connections[GKYL_MAX_CDIM][2]; // block connections
};

typedef struct gkyl_gk_block_geom gkyl_gk_block_geom;

/**
 * Construct a new empty N-dim block geometry with given total number
 * of blocks. The returned block geometry needs to constructed using
 * the set_block method before it can be used.
 *
 * @param ndim Dimension of space
 * @param nblocks Total number of blocks in geometry
 * @return New block geometry
 */
struct gkyl_gk_block_geom *gkyl_gk_block_geom_new(int ndim, int nblocks);

/**
 * Return geometry dimension
 *
 * @param bgeom Block geometry
 * @return Dimension
 */
int gkyl_gk_block_geom_ndim(const struct gkyl_gk_block_geom *bgeom);

/**
 * Return number of blocks in domain
 *
 * @param bgeom Block geometry
 * @return number of blocks
 */
int gkyl_gk_block_geom_num_blocks(const struct gkyl_gk_block_geom *bgeom);

/**
 * Acquire pointer to block-geometry. The pointer must be released
 * using release method.
 *
 * @param bgeom Block geometry to which reference is required
 * @return Pointer to acquired block-topo
 */
struct gkyl_gk_block_geom* gkyl_gk_block_geom_acquire(const struct gkyl_gk_block_geom *bgeom);

/**
 * Acquire a pointer to the block topology. The caller must release
 * the returned object by calling the gkyl_block_topo_release method.
 *
 * @param bgeom Geometry object from which to fetch topology
 * @return topology object
 */
struct gkyl_block_topo* gkyl_gk_block_geom_topo(const struct gkyl_gk_block_geom *bgeom);

/**
 * Set geometry and connectivity information about a block.
 *
 * @param bgeom Geometry object
 * @param bidx Block index
 * @param info Geometry info for block @a bidx
 *
 */
void gkyl_gk_block_geom_set_block(struct gkyl_gk_block_geom *bgeom, int bidx,
  const struct gkyl_gk_block_geom_info *info);


/**
 * Reset grid extents for block geometry info
 *
 * @param bgeom Geometry object
 * @param bidx Block index
 * @param lower Lower extents
 * @param upper Upper extents
 */
void
gkyl_gk_block_geom_reset_block_extents(struct gkyl_gk_block_geom *bgeom, int bidx,
  double *lower, double *upper);

/**
 * Get geometry and connectivity information about a block.
 *
 * @param bgeom Geometry object
 * @param bidx Block index
 * @return Geometry info for block @a bidx
 *
 */
const struct gkyl_gk_block_geom_info *gkyl_gk_block_geom_get_block(
  const struct gkyl_gk_block_geom *bgeom, int bidx);
    
/**
 * Check consistency of block geometry: the geometry typically has
 * redundant data in it. This method ensures the redundant data is
 * consistent. It also checks if there are any blocks with unspecified
 * connections.
 *
 * @param bgeom Block geometry to check
 * @return 1 if geometry is consistent, 0 otherwise.
 */
int gkyl_gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom);

/**
 * Free block geometry.
 *
 * @return Block geometry to free
 */
void gkyl_gk_block_geom_release(struct gkyl_gk_block_geom* bgeom);
