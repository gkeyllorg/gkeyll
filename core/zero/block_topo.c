#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_block_topo.h>
#include <gkyl_elem_type_priv.h>

#include <mpack.h>

// for use in consistency checking
static const enum gkyl_oriented_edge complimentary_edges[] = {
  [0] = 0, // can't happen for fully-specified edges
  [GKYL_LOWER_POSITIVE] = GKYL_UPPER_POSITIVE,
  [GKYL_LOWER_NEGATIVE] = GKYL_UPPER_NEGATIVE,
  [GKYL_UPPER_POSITIVE] = GKYL_LOWER_POSITIVE,
  [GKYL_UPPER_NEGATIVE] = GKYL_LOWER_NEGATIVE,
  [GKYL_PHYSICAL] = GKYL_PHYSICAL,
};

static const char *block_edge_names[] = { "lower", "upper" };

static void
block_topo_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_block_topo *btopo = container_of(ref, struct gkyl_block_topo, ref_count);
  gkyl_free(btopo->conn);
  gkyl_free(btopo);
}

static struct gkyl_msgpack_data *
btopo_create_mpack(const struct gkyl_block_topo *btopo)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));
  mt->meta_sz = 0;
  mt->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);

  mpack_write_cstr(&writer, "ndim");
  mpack_write_i64(&writer, btopo->ndim);
  
  mpack_write_cstr(&writer, "num_blocks");
  mpack_write_i64(&writer, btopo->num_blocks);

  // 3 values written (block ID, direction, edge) for each of the 2 edges:
  size_t num_arr_elems = btopo->num_blocks*btopo->ndim*3*2;
  // write each block connectivity into an array
  mpack_write_cstr(&writer, "connections");
  mpack_start_array(&writer, num_arr_elems);

  // Store connectivity data as a flat array
  for (int i=0; i<btopo->num_blocks; ++i) {
    for (int d=0; d<btopo->ndim; ++d) {
      for (int e=0; e<2; ++e) {
        mpack_write_i64(&writer, btopo->conn[i].connections[d][e].bid);
        mpack_write_i64(&writer, btopo->conn[i].connections[d][e].dir);
        mpack_write_i64(&writer, btopo->conn[i].connections[d][e].edge);
      }
    }
  }

  mpack_finish_array(&writer);
  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }  

  return mt;
}

static void
btopo_array_meta_release(struct gkyl_msgpack_data *amet)
{
  if (!amet) return;
  MPACK_FREE(amet->meta);
  gkyl_free(amet);
}

struct gkyl_block_topo*
gkyl_block_topo_new(int ndim, int nblocks)
{
  struct gkyl_block_topo *btopo = gkyl_malloc(sizeof(struct gkyl_block_topo));
  btopo->ndim = ndim;
  btopo->num_blocks = nblocks;
  btopo->conn = gkyl_calloc(sizeof(struct gkyl_block_connections), nblocks);

  btopo->ref_count = gkyl_ref_count_init(block_topo_free);

  return btopo;
}

int
gkyl_block_topo_check_consistency(const struct gkyl_block_topo *btopo)
{
  for (int i=0; i<btopo->num_blocks; ++i) {
    for (int d=0; d<btopo->ndim; ++d) {
      
      const struct gkyl_target_edge *te = btopo->conn[i].connections[d];
      
      for (int e=0; e<2; ++e) { // 0: lower, 1: upper
        if (te[e].edge < 1) // unspecified edges are defaulted to 0
          return 0;

        // check consistency
        if (te[e].edge != GKYL_PHYSICAL) {
          if (te[e].bid < 0 || te[e].bid >= btopo->num_blocks) // improperly numbered block
            return 0;

          // fetch edge which should point back to ith-block
          const struct gkyl_target_edge *te_back =
            &btopo->conn[te[e].bid].connections[d][(e+1)%2];

          // check if the edge belongs to block 'i'
          if (te_back->bid != i)
            return 0;

          // check if edge orientation is complimentary
          if (te_back->edge != complimentary_edges[te[e].edge])
            return 0;
        }
        
      }
    }
  }

  return 1;
}

int
gkyl_block_topo_write(const struct gkyl_block_topo *btopo, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    struct gkyl_msgpack_data *amet = btopo_create_mpack(btopo);
    if (amet) {
      status = gkyl_header_meta_write_fp( &(struct gkyl_array_header_info) {
          .file_type = gkyl_file_type_int[GKYL_BLOCK_TOPO_DATA_FILE],
          .meta_size = amet->meta_sz,
          .meta = amet->meta
        },
        fp
      );
      btopo_array_meta_release(amet);
    }
    else {
      status = GKYL_ARRAY_RIO_META_FAILED;
    }
  }
  return status;
}

struct gkyl_block_topo *
gkyl_block_topo_read(const char *fname, int *status)
{
  struct gkyl_block_topo *btopo = 0;
  
  enum gkyl_array_rio_status read_status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "r") {
    struct gkyl_array_header_info hdr = { };
    read_status = gkyl_header_meta_read_fp(&hdr, fp);
    if (hdr.file_type != 4) {
      read_status = GKYL_ARRAY_RIO_BAD_VERSION;
      break;
    }
    if (hdr.meta_size > 0) {
      // parse header and construct a new block-topo
      mpack_tree_t tree;
      mpack_tree_init_data(&tree, hdr.meta, hdr.meta_size);
      mpack_tree_parse(&tree);
      mpack_node_t root = mpack_tree_root(&tree);

      mpack_node_t ndim_node = mpack_node_map_cstr(root, "ndim");
      int ndim = mpack_node_i64(ndim_node);

      mpack_node_t num_blocks_node = mpack_node_map_cstr(root, "num_blocks");
      int num_blocks = mpack_node_i64(num_blocks_node);

      btopo = gkyl_block_topo_new(ndim, num_blocks);
      mpack_node_t conn_array_node = mpack_node_map_cstr(root, "connections");

      int array_idx = 0;
      for (int i=0; i<num_blocks; ++i) {

        struct gkyl_block_connections conn = { };
        
        for (int d=0; d<ndim; ++d) {
          for (int e=0; e<2; ++e) {
            mpack_node_t bid_node = mpack_node_array_at(conn_array_node, array_idx++);
            conn.connections[d][e].bid = mpack_node_i64(bid_node);
            
            mpack_node_t dir_node = mpack_node_array_at(conn_array_node, array_idx++);
            conn.connections[d][e].dir = mpack_node_i64(dir_node);

            mpack_node_t edge_node = mpack_node_array_at(conn_array_node, array_idx++);
            conn.connections[d][e].edge = mpack_node_i64(edge_node);
          }
        }
        btopo->conn[i] = conn;
      }
      
      mpack_tree_destroy(&tree);
      gkyl_free(hdr.meta);
    }
  }
  
  *status = read_status;
  return btopo;
}

struct gkyl_block_topo *
gkyl_block_topo_acquire(const struct gkyl_block_topo* btopo)
{
  gkyl_ref_count_inc(&btopo->ref_count);
  return (struct gkyl_block_topo*) btopo;
}     

void
gkyl_block_topo_release(struct gkyl_block_topo* btopo)
{
  gkyl_ref_count_dec(&btopo->ref_count);
}
