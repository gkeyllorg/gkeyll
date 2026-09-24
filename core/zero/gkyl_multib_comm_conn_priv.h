#pragma once

#include <gkyl_multib_comm_conn.h>

// Private header file for multib_comm_conn. Do not include in user-facing
// header files.

//
// Functions for a null_comm.
//
int gkyl_multib_comm_conn_array_transfer_null(struct gkyl_comm *comm,
  int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv);

//
// Functions for a mpi_comm.
//
int gkyl_multib_comm_conn_array_transfer_mpi(struct gkyl_comm *comm,
  int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv);

//
// Functions for a nccl_comm.
//
int gkyl_multib_comm_conn_array_transfer_nccl(struct gkyl_comm *comm,
  int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv);
