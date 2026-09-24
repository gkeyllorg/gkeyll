#include <gkyl_block_geom.h>
#include <gkyl_alloc.h>

#include <string.h>

// Geometry info for all blocks in simulation
struct gkyl_block_geom {
  int ndim; // dimension
  int num_blocks; // total number of blocks
  struct gkyl_block_geom_info *blocks; // info for each block
  struct gkyl_block_topo *btopo; // topology of blocks
  
  struct gkyl_ref_count ref_count;
};

static void
block_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_block_geom *bgeom = container_of(ref, struct gkyl_block_geom, ref_count);
  gkyl_free(bgeom->blocks);
  gkyl_block_topo_release(bgeom->btopo);
  gkyl_free(bgeom);
}

struct gkyl_block_geom*
gkyl_block_geom_new(int ndim, int nblocks)
{
  struct gkyl_block_geom *bgeom = gkyl_malloc(sizeof(struct gkyl_block_geom));
  bgeom->ndim = ndim;
  bgeom->num_blocks = nblocks;
  bgeom->blocks = gkyl_calloc(sizeof(struct gkyl_block_geom_info), nblocks);

  bgeom->btopo = gkyl_block_topo_new(ndim, nblocks);

  bgeom->ref_count = gkyl_ref_count_init(block_geom_free);

  return bgeom;
}

int
gkyl_block_geom_ndim(const struct gkyl_block_geom *bgeom)
{
  return bgeom->ndim;
}

int
gkyl_block_geom_num_blocks(const struct gkyl_block_geom *bgeom)
{
  return bgeom->num_blocks;
}

void
gkyl_block_geom_set_block(struct gkyl_block_geom *bgeom, int bidx,
  const struct gkyl_block_geom_info *info)
{
  memcpy(&bgeom->blocks[bidx], info, sizeof(struct gkyl_block_geom_info));
  
  for (int d=0; d<bgeom->ndim; ++d)
    bgeom->blocks[bidx].cuts[d] = info->cuts[d] > 0 ? info->cuts[d] : 1;
  
  // set topology information
  for (int i=0; i<bgeom->ndim; ++i)
    for (int e=0; e<2; ++e)
      bgeom->btopo->conn[bidx].connections[i][e] = info->connections[i][e];  
}

void
gkyl_block_geom_reset_block_extents(struct gkyl_block_geom *bgeom, int bidx, double *lower, double *upper)
{
  struct gkyl_block_geom_info *bgi = &bgeom->blocks[bidx];
  for (int i = 0; i < bgeom->ndim; ++i) {
    bgi->lower[i] = lower[i];
    bgi->upper[i] = upper[i];
  }
}

const struct gkyl_block_geom_info*
gkyl_block_geom_get_block(const struct gkyl_block_geom *bgeom, int bidx)
{
  return &bgeom->blocks[bidx];
}

int
gkyl_block_geom_check_consistency(const struct gkyl_block_geom *bgeom)
{
  // MORE TESTS ARE NEEDED HERE
  return gkyl_block_topo_check_consistency(bgeom->btopo);
}

struct gkyl_block_geom *
gkyl_block_geom_acquire(const struct gkyl_block_geom* bgeom)
{
  gkyl_ref_count_inc(&bgeom->ref_count);
  return (struct gkyl_block_geom*) bgeom;
}

struct gkyl_block_topo*
gkyl_block_geom_topo(const struct gkyl_block_geom *bgeom)
{
  return gkyl_block_topo_acquire(bgeom->btopo);
}

void
gkyl_block_geom_release(struct gkyl_block_geom* bgeom)
{
  gkyl_ref_count_dec(&bgeom->ref_count);
}
