#include <gkyl_alloc.h>
#include <gkyl_rrobin_decomp.h>

// internal struct to represent the round-robin algo data
struct rrobin_decomp {
  struct gkyl_rrobin_decomp rrobin;
  int *branks; // ranks per block
};

const struct gkyl_rrobin_decomp*
gkyl_rrobin_decomp_new(int total_ranks, int nblocks, const int *branks)
{
  struct rrobin_decomp *rr = gkyl_malloc(sizeof(*rr));
  rr->rrobin.total_ranks = total_ranks;  
  rr->rrobin.nblocks = nblocks;

  rr->branks = gkyl_malloc(sizeof(int[nblocks]));
  for (int i=0; i<nblocks; ++i)
    rr->branks[i] = branks[i];
  
  int tot_branks = 0;
  for (int i=0; i<nblocks; ++i)
    tot_branks += branks[i];

  return &rr->rrobin;
}

int
gkyl_rrobin_decomp_nranks(const struct gkyl_rrobin_decomp *rr, int bn)
{
  struct rrobin_decomp *rrd = container_of(rr, struct rrobin_decomp, rrobin);  
  return rrd->branks[bn];
}

void
gkyl_rrobin_decomp_getranks(const struct gkyl_rrobin_decomp *rr, int bn, int ranks[])
{
  struct rrobin_decomp *rrd = container_of(rr, struct rrobin_decomp, rrobin);

  // ranks are distributed in a round-robin way amongst the blocks
  int loc = 0;
  for (int i=0; i<bn; ++i)
    loc += rrd->branks[i];
  int start = loc % rr->total_ranks;
  for (int i=0; i<rrd->branks[bn]; ++i)
    ranks[i] = (start+i) % rr->total_ranks;
}

void
gkyl_rrobin_decomp_release(const struct gkyl_rrobin_decomp *rr)
{
  struct rrobin_decomp *rrd = container_of(rr, struct rrobin_decomp, rrobin);

  gkyl_free(rrd->branks);
  gkyl_free(rrd);
}
