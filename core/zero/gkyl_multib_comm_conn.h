#pragma once

#include <gkyl_block_topo.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_ref_count.h>
#include <gkyl_comm.h>
#include <gkyl_array.h>

enum multib_send_recv { GKYL_COMM_CONN_SEND, GKYL_COMM_CONN_RECV };

// Information for block send/recv
struct gkyl_comm_conn {
  enum multib_send_recv sr; // flag indicating send/recv
  int block_id; // send/recv block ID
  int rank; // send/recv range ID in block
  struct gkyl_range range; // send/recv range
  enum gkyl_oriented_edge src_edge;
  enum gkyl_oriented_edge tar_edge;
};

// List of send/recv for a given rank
struct gkyl_multib_comm_conn {
  int num_comm_conn; // number of send/recv
  struct gkyl_comm_conn *comm_conn; // communication connections (size num_comm_conn)
  struct gkyl_ref_count ref_count;  
};

/**
 * Create new communication connection object. The list @a comm_conn
 * is managed by the caller.
 *
 * @param num Number of communication connections
 * @param comm_conn List of individual communication connections
 * @return New communication connection object
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new(int num,
  const struct gkyl_comm_conn *comm_conn);

/**
 * Construct the send communication connections for a rank from its
 * local block rank and topology connection.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param block_conn Topological connections for block
 * @param nghost Number of ghost cells in direction d is nghost[d]
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object for sends
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_send(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp);

/**
 * Construct the received communication connections for a rank from
 * its local block rank and topology connection.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param block_conn Topological connections for block
 * @param nghost Number of ghost cells in direction d is nghost[d]
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object receives
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_recv(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp);

/**
 * Construct the send communication connections for a rank from its
 * local block rank and a list of blocks connected to it along a direction.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param nconnected Number of blocks including self connected along direction
 * @param block list Ordered (based on topology) list of connected block ids (including self)
 * @param dir direction in which blocks are connected
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object for sends
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_send_from_connections(
  int block_id, int block_rank, const int *nghost,
  int nconnected, int* block_list, int dir,
  struct gkyl_rect_decomp **decomp);

/**
 * Construct the recv communication connections for a rank from its
 * local block rank and a list of blocks connected to it along a direction.
 *
 * @param block_id ID of block
 * @param block_rank Local rank in block
 * @param nconnected Number of blocks including self connected along direction
 * @param block list Ordered (based on topology) list of connected block ids (including self)
 * @param dir direction in which blocks are connected
 * @param decomp List of decomposition objects for each block
 * @return New communication connection object for sends
 */
struct gkyl_multib_comm_conn *gkyl_multib_comm_conn_new_recv_from_connections(
  int block_id, int block_rank, const int *nghost,
  int nconnected, int* block_list, int dir,
  struct gkyl_rect_decomp **decomp);

/**
 * Transfer data from 'ain' and to 'aout' according to connections in
 * 'mbcc_send' and 'mbcc_recv'.
 *
 * @param comm Multiblock comm object.
 * @param num_blocks_local Number of blocks in this MPI process.
 * @param local_blocks Block IDs of the local blocks.
 * @param mbcc_send Sending connections between ranks (for each local block).
 * @param mbcc_recv Receiving connections between ranks (for each local block).
 * @param arr_send Array to send from (for each local block). 
 * @param arr_recv Array to receive into (for each local block).
 */
int gkyl_multib_comm_conn_array_transfer(struct gkyl_comm *comm,
  int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv);

/**
 * Create a multib range and extended range that spans
 * a set of blocks in a direction
 * @param multib_range_ext on output extended range spanning all connected blocks in direction
 * @param multib_range on output the range spanning all connected blocks in direction
 * @param nghost Number of ghost cells in direction d is nghost[d]
 * @param nconnected Number of blocks including self connected along direction
 * @param block list Ordered (based on topology) list of connected block ids (including self)
 * @param dir direction in which blocks are connected
 * @param decomp List of decomposition objects for each block
 */
void gkyl_multib_comm_conn_create_multib_ranges_in_dir(struct gkyl_range *multib_range_ext,
  struct gkyl_range *multib_range, const int *nghost, int nconnected, 
  int* block_list, int dir, struct gkyl_rect_decomp **decomp);

/**
 * Sort the connections in ascending order according to rank, and block id.
 *
 * @param comm_conn List of individual communication connections.
 */
void gkyl_multib_comm_conn_sort(struct gkyl_multib_comm_conn *mbcc);

/**
 * Release comm-conn object.
 *
 * @param cconn Object to release
 */
void gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn);
