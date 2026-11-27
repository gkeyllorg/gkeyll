#include <gkyl_multib_comm_conn_priv.h>
#include <gkyl_comm_priv.h>
#include <gkyl_null_comm_priv.h>
#include <assert.h>

int
gkyl_multib_comm_conn_array_transfer_null(
  struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv
)
{
  struct null_comm *null = container_of(comm, struct null_comm, priv_comm.pub_comm);

  // post nonblocking recv to get data into ghost-cells
  int nridx = 0;
  for (int bI = 0; bI < num_blocks_local; ++bI) {
    int bid = local_blocks[bI];
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];

    for (int nr = 0; nr < mbcc_r->num_comm_conn; ++nr) {
      int bid_src = mbcc_r->comm_conn[nr].block_id;
      struct gkyl_range *range_dest = &mbcc_r->comm_conn[nr].range;
      int e = mbcc_r->comm_conn[nr].tar_edge;

      // The send connections may not be ordered in the same way as the recvs.
      // Loop over the send connections and find the one from the same block.
      int bid_src_idx = -1;
      for (int cI = 0; cI < num_blocks_local; ++cI) {
        if (local_blocks[cI] == bid_src) {
          bid_src_idx = cI;
          break;
        }
      }
      assert(
        (bid_src_idx >= 0) &&
        "gkyl_multib_comm_conn_array_transfer_null: recv connection references a "
        "source block that is not in local_blocks (block topology is inconsistent "
        "with the local block list)"
      );
      struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bid_src_idx];
      int conn_src_idx = -1;
      for (int ns = 0; ns < mbcc_s->num_comm_conn; ++ns) {
        if (mbcc_s->comm_conn[ns].block_id == bid && mbcc_s->comm_conn[ns].src_edge == e) {
          conn_src_idx = ns;
          break;
        }
      }
      assert(
        (conn_src_idx >= 0) &&
        "gkyl_multib_comm_conn_array_transfer_null: no send connection matches "
        "this recv connection's (block_id, tar_edge) -- the block_geom_info "
        "connections[][] entries for this pair of blocks are not mutually "
        "consistent (e.g. one side's declared target edge doesn't match the "
        "local edge the other side actually uses for this connection)"
      );

      struct gkyl_range *range_src = &mbcc_s->comm_conn[conn_src_idx].range;

      gkyl_array_copy_range_to_range(arr_recv[bI], arr_send[bid_src_idx], range_dest, range_src);
    }
  }

  return 0;
}
