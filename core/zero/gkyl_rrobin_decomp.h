#pragma once

// Round-robin decomposition
struct gkyl_rrobin_decomp {
  int total_ranks; // total number of ranks in decomp
  int nblocks; // number of blocks in decomposition
};  

/**
 * Create a new round-robin decomposition with given @a total_ranks,
 * number of "blocks" and the ranks needed for each block.
 *
 * @param total_ranks Total number of ranks available
 * @param nblocks Number of "blocks" in system
 * @param branks An array of @a nblocks size, each with ranks in each block
 * @return New round-robin decomposition
 */
const struct gkyl_rrobin_decomp* gkyl_rrobin_decomp_new(int total_ranks, int nblocks,
  const int *branks);

/**
 * Get number of ranks in block @a bn.
 *
 * @param rr Round-robin decomp
 * @param bn Block nummber for which to fetch rank list
 * @return Number of ranks in block @a bn
 */
int gkyl_rrobin_decomp_nranks(const struct gkyl_rrobin_decomp *rr, int bn);

/**
 * Get the list of ranks in block @a bn. The ranks are returned in @a
 * ranks. This must be big enough to store the list.
 *
 * @param rr Round-robin decomp
 * @param bn Block nummber for which to fetch rank list
 * @param ranks List of ranks in @a bn block
 */
void gkyl_rrobin_decomp_getranks(const struct gkyl_rrobin_decomp *rr, int bn, int ranks[]);

/**
 * Release round-robin decomposition.
 *
 * @param rr Round-robin decomp to release
 */
void gkyl_rrobin_decomp_release(const struct gkyl_rrobin_decomp *rr);
