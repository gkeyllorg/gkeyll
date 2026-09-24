#include <acutest.h>

#ifdef GKYL_HAVE_NCCL

#include <math.h>
#include <stc/cstr.h>
#include <gkyl_util.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_nccl_comm.h>
#include <gkyl_rrobin_decomp.h>

void
nccl_allreduce()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int n_rank, n_sz;
  gkyl_comm_get_rank(comm_dev, &n_rank);
  gkyl_comm_get_size(comm_dev, &n_sz);
  TEST_CHECK( n_rank == m_rank );
  TEST_CHECK( n_sz == m_sz );

  double vals_ho[2];
  if (m_rank == 0) {
    vals_ho[0] = 1.0;
    vals_ho[1] = 2.0;
  }
  if (m_rank == 1) {
    vals_ho[0] = 3.0;
    vals_ho[1] = -1.0;
  }
  double *vals = gkyl_cu_malloc(2*sizeof(double));
  gkyl_cu_memcpy(vals, vals_ho, 2*sizeof(double), GKYL_CU_MEMCPY_H2D);

  double v_max_ho[2], v_min_ho[2], v_sum_ho[2];
  double *v_max = gkyl_cu_malloc(2*sizeof(double));
  double *v_min = gkyl_cu_malloc(2*sizeof(double));
  double *v_sum = gkyl_cu_malloc(2*sizeof(double));

  gkyl_comm_allreduce(comm_dev, GKYL_DOUBLE, GKYL_MAX, 2, vals, v_max);
  gkyl_cu_memcpy(v_max_ho, v_max, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_max_ho[0] == 3.0 );
  TEST_CHECK( v_max_ho[1] == 2.0 );

  gkyl_comm_allreduce(comm_dev, GKYL_DOUBLE, GKYL_MIN, 2, vals, v_min);
  gkyl_cu_memcpy(v_min_ho, v_min, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_min_ho[0] == 1.0 );
  TEST_CHECK( v_min_ho[1] == -1.0 );

  gkyl_comm_allreduce(comm_dev, GKYL_DOUBLE, GKYL_SUM, 2, vals, v_sum);
  gkyl_cu_memcpy(v_sum_ho, v_sum, 2*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( v_sum_ho[0] == 4.0 );
  TEST_CHECK( v_sum_ho[1] == 1.0 );

  gkyl_cu_free(vals);
  gkyl_cu_free(v_max);
  gkyl_cu_free(v_min);
  gkyl_cu_free(v_sum);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm_ho);
  gkyl_comm_release(comm_dev);
}

void
nccl_n2_allgather_1d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range global;
  gkyl_range_init(&global, 1, (int[]) { 1 }, (int[]) { 10 });

  int cuts[] = { 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(global.ndim, cuts, &global);
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_local = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, global.volume);
  struct gkyl_array *arr_local_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global_ho = gkyl_array_new(GKYL_DOUBLE, 1, global.volume);
  gkyl_array_clear(arr_global, 0.0);

  gkyl_array_clear(arr_local_ho, 200005.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_local_ho, idx);
    f[0] = idx+10.0*rank;
  }
  gkyl_array_copy(arr_local, arr_local_ho);

  gkyl_comm_array_allgather(comm, &local, &global, arr_local, arr_global);

  gkyl_array_copy(arr_global_ho, arr_global);
  struct gkyl_range_iter iter_global;
  gkyl_range_iter_init(&iter_global, &global);
  while (gkyl_range_iter_next(&iter_global)) {
    long idx = gkyl_range_idx(&global, iter_global.idx);
    double *f = gkyl_array_fetch(arr_global_ho, idx);
    // first 5 entries are 1-5, second 5 entries are 11-15
    if (idx < local.volume)
      TEST_CHECK( idx+1.0 == f[0] );
    else 
      TEST_CHECK( idx+6.0 == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_local);      
  gkyl_array_release(arr_global); 
  gkyl_array_release(arr_local_ho);      
  gkyl_array_release(arr_global_ho); 
}

void
nccl_n4_allgather_2d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create global range
  int cells[] = { 10, 10 };
  struct gkyl_range global;
  gkyl_create_global_range(2, cells, &global);

  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &global);  
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_local = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, global.volume);
  struct gkyl_array *arr_local_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global_ho = gkyl_array_new(GKYL_DOUBLE, 1, global.volume);
  gkyl_array_clear(arr_global, 0.0);

  gkyl_array_clear(arr_local_ho, 200005.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_local_ho, idx);
    f[0] = iter.idx[0] + iter.idx[1]*(rank+1.0) + 10.0*rank;
  } 
  gkyl_array_copy(arr_local, arr_local_ho);

  gkyl_comm_array_allgather(comm, &local, &global, arr_local, arr_global);

  gkyl_array_copy(arr_global_ho, arr_global);
  struct gkyl_range_iter iter_global;
  gkyl_range_iter_init(&iter_global, &global);
  while (gkyl_range_iter_next(&iter_global)) {
    long idx = gkyl_range_idx(&global, iter_global.idx);
    double *f = gkyl_array_fetch(arr_global_ho, idx);
    // check value of {2, 2} decomp organized as 
    // rank 0 owns {1, 1} to {5, 5}
    // rank 1 owns {1, 6} to {5, 10} 
    // rank 2 owns {6, 1} to {10, 5} 
    // rank 3 owns {6, 6} to {10, 10}
    double val;
    if (iter_global.idx[0] <= cells[0]/cuts[0] && iter_global.idx[1] <= cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1];
    else if (iter_global.idx[0] <= cells[0]/cuts[0] && iter_global.idx[1] > cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1]*2.0 + 10.0;
    else if (iter_global.idx[0] > cells[0]/cuts[0] && iter_global.idx[1] <= cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1]*3.0 + 20.0;
    else 
      val = iter_global.idx[0] + iter_global.idx[1]*4.0 + 30.0;
    TEST_CHECK( val == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_local);
  gkyl_array_release(arr_global);
  gkyl_array_release(arr_local_ho);
  gkyl_array_release(arr_global_ho);
}


void
nccl_n2_allgather_1d_host()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range global;
  gkyl_range_init(&global, 1, (int[]) { 1 }, (int[]) { 10 });

  int cuts[] = { 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(global.ndim, cuts, &global);
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_local_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global_ho = gkyl_array_new(GKYL_DOUBLE, 1, global.volume);

  gkyl_array_clear(arr_local_ho, 200005.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_local_ho, idx);
    f[0] = idx+10.0*rank;
  }

  gkyl_comm_array_allgather_host(comm, &local, &global, arr_local_ho, arr_global_ho);

  struct gkyl_range_iter iter_global;
  gkyl_range_iter_init(&iter_global, &global);
  while (gkyl_range_iter_next(&iter_global)) {
    long idx = gkyl_range_idx(&global, iter_global.idx);
    double *f = gkyl_array_fetch(arr_global_ho, idx);
    // first 5 entries are 1-5, second 5 entries are 11-15
    if (idx < local.volume)
      TEST_CHECK( idx+1.0 == f[0] );
    else 
      TEST_CHECK( idx+6.0 == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_local_ho);      
  gkyl_array_release(arr_global_ho); 
}

void
nccl_n4_allgather_2d_host()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create global range
  int cells[] = { 10, 10 };
  struct gkyl_range global;
  gkyl_create_global_range(2, cells, &global);

  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &global);  
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_local_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr_global_ho = gkyl_array_new(GKYL_DOUBLE, 1, global.volume);

  gkyl_array_clear(arr_local_ho, 200005.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_local_ho, idx);
    f[0] = iter.idx[0] + iter.idx[1]*(rank+1.0) + 10.0*rank;
  } 

  gkyl_comm_array_allgather_host(comm, &local, &global, arr_local_ho, arr_global_ho);

  struct gkyl_range_iter iter_global;
  gkyl_range_iter_init(&iter_global, &global);
  while (gkyl_range_iter_next(&iter_global)) {
    long idx = gkyl_range_idx(&global, iter_global.idx);
    double *f = gkyl_array_fetch(arr_global_ho, idx);
    // check value of {2, 2} decomp organized as 
    // rank 0 owns {1, 1} to {5, 5}
    // rank 1 owns {1, 6} to {5, 10} 
    // rank 2 owns {6, 1} to {10, 5} 
    // rank 3 owns {6, 6} to {10, 10}
    double val;
    if (iter_global.idx[0] <= cells[0]/cuts[0] && iter_global.idx[1] <= cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1];
    else if (iter_global.idx[0] <= cells[0]/cuts[0] && iter_global.idx[1] > cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1]*2.0 + 10.0;
    else if (iter_global.idx[0] > cells[0]/cuts[0] && iter_global.idx[1] <= cells[1]/cuts[1])
      val = iter_global.idx[0] + iter_global.idx[1]*3.0 + 20.0;
    else 
      val = iter_global.idx[0] + iter_global.idx[1]*4.0 + 30.0;
    TEST_CHECK( val == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_local_ho);
  gkyl_array_release(arr_global_ho);
}

// MF 2024/09/12: disable these for now per 498b7d1569eaa9285ae59581bd22dab124672f7b.
// void
// nccl_n2_array_send_irecv_2d()
// {
//   // Test array_send and array_recv with a nonblocking comm.
//   struct gkyl_range range;
//   gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
// 
//   int comm_size;
//   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
//   if (comm_size != 2) return;
//   
//   int cuts[] = { comm_size, 1 };
//   struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
//   
//   struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
//       .mpi_comm = MPI_COMM_WORLD,
//       .decomp = decomp
//     }
//   );
// 
//   int rank, sz;
//   gkyl_comm_get_rank(comm_ho, &rank);
//   gkyl_comm_get_size(comm_ho, &sz);
// 
//   struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
//       .mpi_comm = MPI_COMM_WORLD,
//       .decomp = decomp,
//     }
//   );
// 
//   // Allocate send/recv buffers.
//   double sendval = rank==0? 20005. : 30005.;
//   double recvval = rank==0? 30005. : 20005.;
// 
//   // Assume the range is not decomposed.
//   struct gkyl_array *arrA = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
//   struct gkyl_array *arrB = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
//   struct gkyl_array *recvbuff_ho = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
//   gkyl_array_clear(arrA, sendval*(1-rank));
//   gkyl_array_clear(arrB, sendval*rank);
// 
//   struct gkyl_array *recvbuff = rank==0? arrB : arrA;
//   struct gkyl_array *sendbuff = rank==0? arrA : arrB;
// 
//   struct gkyl_comm_state *cstate = gkyl_comm_state_new(comm_dev);
//   int tag = 13;
//   // Communicate data from rank 0 to rank 1.
//   if (rank == 1)
//     gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
//   if (rank == 0)
//     gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);
// 
//   if (rank == 1) {
//     gkyl_comm_state_wait(comm_dev, cstate);
// 
//     gkyl_array_copy(recvbuff_ho, recvbuff);
//     struct gkyl_range_iter iter;
//     gkyl_range_iter_init(&iter, &range);
//     while (gkyl_range_iter_next(&iter)) {
//       long idx = gkyl_range_idx(&range, iter.idx);
//       const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//       TEST_CHECK( f[0] == recvval );
//     }
//   }
// 
//   // Communicate data from rank 1 to rank 0.
//   if (rank == 0)
//     gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
//   if (rank == 1)
//     gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);
// 
//   if (rank == 0) {
//     gkyl_comm_state_wait(comm_dev, cstate);
// 
//     gkyl_array_copy(recvbuff_ho, recvbuff);
//     struct gkyl_range_iter iter;
//     gkyl_range_iter_init(&iter, &range);
//     while (gkyl_range_iter_next(&iter)) {
//       long idx = gkyl_range_idx(&range, iter.idx);
//       const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//       TEST_CHECK( f[0] == recvval );
//     }
//   }
// 
//   // Communicate data between rank 0 and 1 at the same time.
//   gkyl_array_clear(recvbuff, 0.);
//   gkyl_comm_group_call_start(comm_dev);
//   gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate);
//   gkyl_comm_array_send(comm_dev, sendbuff, (rank+1) % 2, tag);
//   gkyl_comm_group_call_end(comm_dev);
//   gkyl_comm_state_wait(comm_dev, cstate);
// 
//   gkyl_array_copy(recvbuff_ho, recvbuff);
//   struct gkyl_range_iter iter;
//   gkyl_range_iter_init(&iter, &range);
//   while (gkyl_range_iter_next(&iter)) {
//     long idx = gkyl_range_idx(&range, iter.idx);
//     const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//     TEST_CHECK( f[0] == recvval );
//   }
// 
//   gkyl_comm_barrier(comm_dev);
// 
//   gkyl_comm_state_release(comm_dev, cstate);
//   gkyl_array_release(recvbuff_ho);
//   gkyl_array_release(arrA);
//   gkyl_array_release(arrB);
// 
//   gkyl_rect_decomp_release(decomp);
//   gkyl_comm_release(comm_dev);
//   gkyl_comm_release(comm_ho);
// }
// 
// void
// nccl_n2_array_isend_irecv_2d()
// {
//   // Test array_send and array_recv with a nonblocking comm.
//   struct gkyl_range range;
//   gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
// 
//   int comm_size;
//   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
//   if (comm_size != 2) return;
//   
//   int cuts[] = { comm_size, 1 };
//   struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
//   
//   struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
//       .mpi_comm = MPI_COMM_WORLD,
//       .decomp = decomp,
//     }
//   );
// 
//   int rank, sz;
//   gkyl_comm_get_rank(comm_ho, &rank);
//   gkyl_comm_get_size(comm_ho, &sz);
// 
//   struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
//       .mpi_comm = MPI_COMM_WORLD,
//       .decomp = decomp,
//     }
//   );
// 
//   // Allocate send/recv buffers.
//   double sendval = rank==0? 20005. : 30005.;
//   double recvval = rank==0? 30005. : 20005.;
// 
//   // Assume the range is not decomposed.
//   struct gkyl_array *arrA = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
//   struct gkyl_array *arrB = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, range.volume);
//   struct gkyl_array *recvbuff_ho = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
//   gkyl_array_clear(arrA, sendval*(1-rank));
//   gkyl_array_clear(arrB, sendval*rank);
// 
//   struct gkyl_array *recvbuff = rank==0? arrB : arrA;
//   struct gkyl_array *sendbuff = rank==0? arrA : arrB;
// 
//   struct gkyl_comm_state *cstate_r = gkyl_comm_state_new(comm_dev);
//   struct gkyl_comm_state *cstate_s = gkyl_comm_state_new(comm_dev);
//   int tag = 13;
//   // Communicate data from rank 0 to rank 1.
//   if (rank == 1)
//     gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate_r);
//   if (rank == 0)
//     gkyl_comm_array_isend(comm_dev, sendbuff, (rank+1) % 2, tag, cstate_s);
// 
//   gkyl_comm_state_wait(comm_dev, cstate_r);
//   gkyl_comm_state_wait(comm_dev, cstate_s);
// 
//   if (rank == 1) {
//     gkyl_array_copy(recvbuff_ho, recvbuff);
//     struct gkyl_range_iter iter;
//     gkyl_range_iter_init(&iter, &range);
//     while (gkyl_range_iter_next(&iter)) {
//       long idx = gkyl_range_idx(&range, iter.idx);
//       const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//       TEST_CHECK( f[0] == recvval );
//     }
//   }
// 
//   // Communicate data from rank 1 to rank 0.
//   if (rank == 0)
//     gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate_r);
//   if (rank == 1)
//     gkyl_comm_array_isend(comm_dev, sendbuff, (rank+1) % 2, tag, cstate_s);
// 
//   gkyl_comm_state_wait(comm_dev, cstate_r);
//   gkyl_comm_state_wait(comm_dev, cstate_s);
// 
//   if (rank == 0) {
//     gkyl_array_copy(recvbuff_ho, recvbuff);
//     struct gkyl_range_iter iter;
//     gkyl_range_iter_init(&iter, &range);
//     while (gkyl_range_iter_next(&iter)) {
//       long idx = gkyl_range_idx(&range, iter.idx);
//       const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//       TEST_CHECK( f[0] == recvval );
//     }
//   }
// 
//   // Communicate data between rank 0 and 1 at the same time.
//   gkyl_array_clear(recvbuff, 0.);
//   gkyl_comm_group_call_start(comm_dev);
//   gkyl_comm_array_irecv(comm_dev, recvbuff, (rank+1) % 2, tag, cstate_r);
//   gkyl_comm_array_isend(comm_dev, sendbuff, (rank+1) % 2, tag, cstate_s);
//   gkyl_comm_group_call_end(comm_dev);
//   gkyl_comm_state_wait(comm_dev, cstate_r);
//   gkyl_comm_state_wait(comm_dev, cstate_s);
// 
//   gkyl_array_copy(recvbuff_ho, recvbuff);
//   struct gkyl_range_iter iter;
//   gkyl_range_iter_init(&iter, &range);
//   while (gkyl_range_iter_next(&iter)) {
//     long idx = gkyl_range_idx(&range, iter.idx);
//     const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//     TEST_CHECK( f[0] == recvval );
//   }
// 
//   gkyl_comm_barrier(comm_dev);
// 
//   gkyl_comm_state_release(comm_dev, cstate_r);
//   gkyl_comm_state_release(comm_dev, cstate_s);
//   gkyl_array_release(recvbuff_ho);
//   gkyl_array_release(arrA);
//   gkyl_array_release(arrB);
// 
//   gkyl_rect_decomp_release(decomp);
//   gkyl_comm_release(comm_dev);
//   gkyl_comm_release(comm_ho);
// }

void
nccl_n2_sync_1d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 10 });

  int cuts[] = { 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);
  
  struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = false,
    }
  );

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = false,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_ho, idx);
    f[0] = iter.idx[0];
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_sync(comm_dev, &local, &local_ext, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);

  gkyl_array_copy(arr_ho, arr);
  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr_ho, idx);
    
    TEST_CHECK( iter.idx[0] == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_comm_release(comm_dev);
  gkyl_comm_release(comm_ho);
}

void
nccl_n4_sync_2d(bool use_corners)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 10 });

  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm_ho = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = use_corners,
    }
  );

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = use_corners,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x, local_ext_x, local_y, local_ext_y;
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) {1, 0},
    &local_ext_x, &local_x);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, 1 },
    &local_ext_y, &local_y);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005);

  gkyl_comm_barrier(comm_dev);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_ho, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_sync(comm_dev, &local, &local_ext, arr);

  gkyl_array_copy(arr_ho, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr_ho, idx);

    if (use_corners) {
      TEST_CHECK( iter.idx[0] == f[0] );
      TEST_CHECK( iter.idx[1] == f[1] );
    }
    else {
      // excludes corners
      if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
        TEST_CHECK( iter.idx[0] == f[0] );
        TEST_CHECK( iter.idx[1] == f[1] );
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_comm_release(comm_dev);
  gkyl_comm_release(comm_ho);
}

void nccl_n4_sync_2d_no_corner() { nccl_n4_sync_2d(false); }
void nccl_n4_sync_2d_use_corner() { nccl_n4_sync_2d(true); }

void
nccl_n4_sync_1x1v()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 512 });

  int cuts[] = { m_sz };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(1, cuts, &range);

  struct gkyl_comm *comm_dev = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  struct gkyl_range vrange;
  gkyl_range_init(&vrange, 1, (int[]) { 1 }, (int[]) { 64 } );

  struct gkyl_rect_decomp *ext_decomp =
    gkyl_rect_decomp_extended_new(&vrange, decomp);

  int nghost[] = { 1, 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&ext_decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_comm *ext_comm = gkyl_comm_extend_comm(comm_dev, &vrange);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005);
  gkyl_comm_barrier(comm_dev);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr_ho, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_sync(ext_comm, &local, &local_ext, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &ext_decomp->parent_range);

  gkyl_array_copy(arr_ho, arr);
  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr_ho, idx);

    TEST_CHECK( iter.idx[0] == f[0] );
    TEST_CHECK( iter.idx[1] == f[1] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_rect_decomp_release(ext_decomp);
  gkyl_comm_release(comm_dev);
  gkyl_comm_release(ext_comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
}

void
nccl_n1_per_sync_2d_tests(int num_per_dirs, int *per_dirs)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 1) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, idx);

    for (int d=0; d<local.ndim; ++d)
      f[d] = iter.idx[d];
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_per_sync(comm, &local, &local_ext, num_per_dirs, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  int count = 0;
  
  gkyl_array_copy(arr_ho, arr);
  for (int id=0; id<num_per_dirs; ++id) {
    int d = per_dirs[id];
    int ncell = gkyl_range_shape(&local, d);

    gkyl_range_iter_init(&iter, &local_ext_x[d]);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n) idx[n] = iter.idx[n];

        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr_ho, lidx);
        for (int n=0; n<local.ndim; ++n) {
          TEST_CHECK( idx[n] == f[n] );
          TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %d | Produced: %.13e", rank, iter.idx[0], iter.idx[1], idx[n], f[n] );
	}
      }
    }
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

void
nccl_n1_per_sync_2d()
{
  int per_dirs_0[] = {0};
  int per_dirs_1[] = {1};
  int per_dirs_01[] = {0,1};

  nccl_n1_per_sync_2d_tests(1, per_dirs_0);
  nccl_n1_per_sync_2d_tests(1, per_dirs_1);
  nccl_n1_per_sync_2d_tests(2, per_dirs_01);

  nccl_n1_per_sync_2d_tests(1, per_dirs_0);
  nccl_n1_per_sync_2d_tests(1, per_dirs_1);
  nccl_n1_per_sync_2d_tests(2, per_dirs_01);
}

void
nccl_n2_per_sync_2d_tests(int *cuts, int num_per_dirs, int *per_dirs)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  // Redefine local_ext_x so it's local shifted in the right direction
  // so it covers the ghost cells of interest.
  int decomp_dir = cuts[0]>1? 0 : 1;
  int delta[] = {0, 0};
  delta[decomp_dir] = 2*rank-1;
  struct gkyl_range local_ext_x_shifted;
  gkyl_range_shift(&local_ext_x_shifted, &local_ext_x[decomp_dir], delta);
  gkyl_sub_range_init(&local_ext_x[decomp_dir], &local_ext, local_ext_x_shifted.lower, local_ext_x_shifted.upper);

  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, idx);

    for (int d=0; d<local.ndim; ++d)
      f[d] = iter.idx[d];
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_per_sync(comm, &local, &local_ext, num_per_dirs, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  int count = 0;
  
  gkyl_array_copy(arr_ho, arr);
  for (int id=0; id<num_per_dirs; ++id) {
    int d = per_dirs[id];
    int ncell = gkyl_range_shape(&range, d);

    gkyl_range_iter_init(&iter, &local_ext_x[d]);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n) idx[n] = iter.idx[n];

        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else if (idx[d] < local.lower[d])
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr_ho, lidx);
        for (int n=0; n<local.ndim; ++n) {
          TEST_CHECK( idx[n] == f[n] );
          TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %d | Produced: %.13e", rank, iter.idx[0], iter.idx[1], idx[n], f[n] );
	}
      }
    }
  }

  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

void
nccl_n2_per_sync_2d()
{
  int cuts_21[] = {2,1};
  int cuts_12[] = {1,2};
  int per_dirs_0[] = {0};
  int per_dirs_1[] = {1};
  int per_dirs_01[] = {0,1};

  nccl_n2_per_sync_2d_tests(cuts_21, 1, per_dirs_0);
  nccl_n2_per_sync_2d_tests(cuts_21, 1, per_dirs_1);
  nccl_n2_per_sync_2d_tests(cuts_21, 2, per_dirs_01);

  nccl_n2_per_sync_2d_tests(cuts_12, 1, per_dirs_0);
  nccl_n2_per_sync_2d_tests(cuts_12, 1, per_dirs_1);
  nccl_n2_per_sync_2d_tests(cuts_12, 2, per_dirs_01);
}

void
nccl_n4_multicomm_2d()
{
  // Test the use of two gkyl_comm objects simultaneously, mimicing the case
  // where one is used to decompose space and the other species.
  // We sync across the conf communicator, and send/recv across the species
  // comm.
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  struct gkyl_comm *worldcomm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = 0,
    }
  );

  int worldrank;
  gkyl_comm_get_rank(worldcomm, &worldrank);

  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 20 });

  int confcuts[] = { 2, 1 };
  struct gkyl_rect_decomp *confdecomp = gkyl_rect_decomp_new_from_cuts(2, confcuts, &range);  
  
  int confcolor = floor(worldrank/confdecomp->ndecomp);
  struct gkyl_comm *confcomm = gkyl_comm_split_comm(worldcomm, confcolor, confdecomp);
  int confrank;
  gkyl_comm_get_rank(confcomm, &confrank);

  int speciescolor = worldrank % confdecomp->ndecomp;
  struct gkyl_rect_decomp *speciesdecomp = 0;
  struct gkyl_comm *speciescomm = gkyl_comm_split_comm(worldcomm, speciescolor, speciesdecomp);
  int speciesrank;
  gkyl_comm_get_rank(speciescomm, &speciesrank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&confdecomp->ranges[confrank], nghost, &local_ext, &local);

  struct gkyl_array *arrA_ho = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arrB_ho = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arrA = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arrB = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arrA_ho, 0.);
  gkyl_array_clear(arrB_ho, 0.);

  // Sync across the conf-space communicator.
  struct gkyl_array *recvbuff = speciesrank==0? arrB : arrA;
  struct gkyl_array *sendbuff = speciesrank==0? arrA : arrB;
  struct gkyl_array *recvbuff_ho = speciesrank==0? arrB_ho : arrA_ho;
  struct gkyl_array *sendbuff_ho = speciesrank==0? arrA_ho : arrB_ho;

  gkyl_comm_barrier(worldcomm);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(sendbuff_ho, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }
  gkyl_array_copy(sendbuff, sendbuff_ho);

  // Sync sendbuff array and check results.
  gkyl_comm_array_sync(confcomm, &local, &local_ext, sendbuff);

  gkyl_array_copy(sendbuff_ho, sendbuff);
  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);
  struct gkyl_range local_x, local_ext_x, local_y, local_ext_y;
  gkyl_create_ranges(&confdecomp->ranges[confrank], (int[]) {1, 0}, &local_ext_x, &local_x);
  gkyl_create_ranges(&confdecomp->ranges[confrank], (int[]) { 0, 1 }, &local_ext_y, &local_y);
  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double *f = gkyl_array_cfetch(sendbuff_ho, idx);
    // exclude corners
    if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
      TEST_CHECK( iter.idx[0] == f[0] );
      TEST_CHECK( iter.idx[1] == f[1] );
    }
  }

// MF 2024/09/12: disable these for now per 498b7d1569eaa9285ae59581bd22dab124672f7b.
//  // Now send/recv across species communicator and check results.
//  struct gkyl_comm_state *cstate = gkyl_comm_state_new(speciescomm);
//  int tag = 13;
//  // Post irecv before send.
//  gkyl_comm_group_call_start(speciescomm);
//  gkyl_comm_array_irecv(speciescomm, recvbuff, (speciesrank+1) % 2, tag, cstate);
//  gkyl_comm_array_send(speciescomm, sendbuff, (speciesrank+1) % 2, tag);
//  gkyl_comm_group_call_end(speciescomm);
//  gkyl_comm_state_wait(speciescomm, cstate);
//
//  gkyl_array_copy(recvbuff_ho, recvbuff);
//  gkyl_range_iter_init(&iter, &in_range);
//  while (gkyl_range_iter_next(&iter)) {
//    long idx = gkyl_range_idx(&in_range, iter.idx);
//    const double *f = gkyl_array_cfetch(recvbuff_ho, idx);
//    // exclude corners
//    if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
//      TEST_CHECK( iter.idx[0] == f[0] );
//      TEST_CHECK( iter.idx[1] == f[1] );
//    }
//  }
//
//  gkyl_comm_state_release(speciescomm, cstate);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);
  gkyl_array_release(arrA_ho);
  gkyl_array_release(arrB_ho);
  gkyl_comm_release(speciescomm);
  gkyl_comm_release(confcomm);
  gkyl_comm_release(worldcomm);
  gkyl_rect_decomp_release(confdecomp);
}
  
static void
nccl_n4_create_comm_from_ranks_1()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
    }
  );

  int branks[2] =  { 2, 2 };
  bool status = false;
  
  const struct gkyl_rrobin_decomp *rrd =
    gkyl_rrobin_decomp_new(m_sz, 2, branks);

  int rb1[4];
  gkyl_rrobin_decomp_getranks(rrd, 0, rb1);

  struct gkyl_comm *comm_b1 =
    gkyl_comm_create_comm_from_ranks(comm, branks[0], rb1, 0, &status);

  if (rank == rb1[0])
    TEST_CHECK( status );
  if (rank == rb1[1])
    TEST_CHECK( status );

  if (comm_b1) {
    int sz_b1;
    gkyl_comm_get_size(comm_b1, &sz_b1);
    TEST_CHECK( branks[0] == sz_b1);
  }

  int rb2[4];
  gkyl_rrobin_decomp_getranks(rrd, 1, rb2);

  struct gkyl_comm *comm_b2 =
    gkyl_comm_create_comm_from_ranks(comm, branks[1], rb2, 0, &status);

  if (rank == rb2[0])
    TEST_CHECK( status );
  if (rank == rb2[1])
    TEST_CHECK( status );

  if (comm_b2) {
    int sz_b2;
    gkyl_comm_get_size(comm_b2, &sz_b2);
    TEST_CHECK( branks[1] == sz_b2);
  }  

  gkyl_rrobin_decomp_release(rrd);
  gkyl_comm_release(comm);
  gkyl_comm_release(comm_b1);
  gkyl_comm_release(comm_b2);
}

static void
nccl_n4_create_comm_from_ranks_2()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
    }
  );

  int branks[2] =  { 4, 2 };
  bool status = false;
  
  const struct gkyl_rrobin_decomp *rrd =
    gkyl_rrobin_decomp_new(m_sz, 2, branks);

  int rb1[4];
  gkyl_rrobin_decomp_getranks(rrd, 0, rb1);

  struct gkyl_comm *comm_b1 =
    gkyl_comm_create_comm_from_ranks(comm, branks[0], rb1, 0, &status);

  if (rank == rb1[0])
    TEST_CHECK( status );
  if (rank == rb1[1])
    TEST_CHECK( status );
  if (rank == rb1[2])
    TEST_CHECK( status );
  if (rank == rb1[3])
    TEST_CHECK( status );

  if (comm_b1) {
    int sz_b1;
    gkyl_comm_get_size(comm_b1, &sz_b1);
    TEST_CHECK( branks[0] == sz_b1);
  }

  int rb2[4];
  gkyl_rrobin_decomp_getranks(rrd, 1, rb2);

  struct gkyl_comm *comm_b2 =
    gkyl_comm_create_comm_from_ranks(comm, branks[1], rb2, 0, &status);

  if (rank == rb2[0])
    TEST_CHECK( status );
  if (rank == rb2[1])
    TEST_CHECK( status );

  if (comm_b2) {
    int sz_b2;
    gkyl_comm_get_size(comm_b2, &sz_b2);
    TEST_CHECK( branks[1] == sz_b2);
  }

  gkyl_rrobin_decomp_release(rrd);
  gkyl_comm_release(comm);
  gkyl_comm_release(comm_b1);
  gkyl_comm_release(comm_b2);
}

void
nccl_bcast_1d()
{
  int m_sz, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bcast_rank = m_sz > 1? 1 : 0;

  struct gkyl_range global;
  gkyl_range_init(&global, 1, (int[]) { 1 }, (int[]) { 8*27*125 });

  int cuts[] = { m_sz };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(global.ndim, cuts, &global);
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    f[0] = linidx+10.0*rank;
  }
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_bcast(comm, arr, arr, bcast_rank);

  gkyl_array_copy(arr_ho, arr);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    TEST_CHECK( linidx+10.0*bcast_rank == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr); 
  gkyl_array_release(arr_ho); 
}

void
nccl_bcast_2d_test(int *cuts)
{
  int m_sz, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bcast_rank = m_sz > 1? 1 : 0;

  // create global range
  int cells[] = { 4*9*25, 4*9*25 };
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_range global;
  gkyl_create_global_range(ndim, cells, &global);

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(ndim, cuts, &global);  
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  struct gkyl_array *arr = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    f[0] = iter.idx[0] + iter.idx[1]*(rank+1.0) + 10.0*rank;
  } 
  gkyl_array_copy(arr, arr_ho);

  gkyl_comm_array_bcast(comm, arr, arr, bcast_rank);

  gkyl_array_copy(arr_ho, arr);
  struct gkyl_range bcast_rank_local, bcast_rank_local_ext;
  gkyl_create_ranges(&decomp->ranges[bcast_rank], nghost, &bcast_rank_local_ext, &bcast_rank_local);
  gkyl_range_iter_init(&iter, &bcast_rank_local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&bcast_rank_local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    double val = iter.idx[0] + iter.idx[1]*(bcast_rank+1.0) + 10.0*bcast_rank;
    TEST_CHECK( val == f[0] );
    TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %.13e | Produced: %.13e", rank, iter.idx[0], iter.idx[1], val, f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_ho);
}

void
nccl_bcast_2d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);

  if (m_sz == 2) {
    int cuts12[] = {1, 2};
    nccl_bcast_2d_test(cuts12);
  
    int cuts21[] = {2, 1};
    nccl_bcast_2d_test(cuts21);
  
  } else if (m_sz == 3) {
    int cuts13[] = {1, 3};
    nccl_bcast_2d_test(cuts13);
  
    int cuts31[] = {3, 1};
    nccl_bcast_2d_test(cuts31);

  } else if (m_sz == 4) {
    int cuts22[] = {2, 2};
    nccl_bcast_2d_test(cuts22);
  
  }
}

void
nccl_bcast_1d_host()
{
  int m_sz, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bcast_rank = m_sz > 1? 1 : 0;

  struct gkyl_range global;
  gkyl_range_init(&global, 1, (int[]) { 1 }, (int[]) { 8*27*125 });

  int cuts[] = { m_sz };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(global.ndim, cuts, &global);
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    f[0] = linidx+10.0*rank;
  }

  gkyl_comm_array_bcast_host(comm, arr_ho, arr_ho, bcast_rank);

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    TEST_CHECK( linidx+10.0*bcast_rank == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_ho); 
}

void
nccl_bcast_2d_host_test(int *cuts)
{
  int m_sz, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int bcast_rank = m_sz > 1? 1 : 0;

  // create global range
  int cells[] = { 4*9*25, 4*9*25 };
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_range global;
  gkyl_create_global_range(ndim, cells, &global);

  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(ndim, cuts, &global);  
  
  struct gkyl_comm *comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, 1, local_ext.volume);
  gkyl_array_clear(arr_ho, 200005.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    f[0] = iter.idx[0] + iter.idx[1]*(rank+1.0) + 10.0*rank;
  } 

  gkyl_comm_array_bcast_host(comm, arr_ho, arr_ho, bcast_rank);

  struct gkyl_range bcast_rank_local, bcast_rank_local_ext;
  gkyl_create_ranges(&decomp->ranges[bcast_rank], nghost, &bcast_rank_local_ext, &bcast_rank_local);
  gkyl_range_iter_init(&iter, &bcast_rank_local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&bcast_rank_local, iter.idx);
    double *f = gkyl_array_fetch(arr_ho, linidx);
    double val = iter.idx[0] + iter.idx[1]*(bcast_rank+1.0) + 10.0*bcast_rank;
    TEST_CHECK( val == f[0] );
    TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %.13e | Produced: %.13e", rank, iter.idx[0], iter.idx[1], val, f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr_ho);
}

void
nccl_bcast_2d_host()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);

  if (m_sz == 2) {
    int cuts12[] = {1, 2};
    nccl_bcast_2d_host_test(cuts12);
  
    int cuts21[] = {2, 1};
    nccl_bcast_2d_host_test(cuts21);
  
  } else if (m_sz == 3) {
    int cuts13[] = {1, 3};
    nccl_bcast_2d_host_test(cuts13);
  
    int cuts31[] = {3, 1};
    nccl_bcast_2d_host_test(cuts31);

  } else if (m_sz == 4) {
    int cuts22[] = {2, 2};
    nccl_bcast_2d_host_test(cuts22);
  
  }
}

  
TEST_LIST = {
  {"nccl_allreduce", nccl_allreduce},
  {"nccl_n2_allgather_1d", nccl_n2_allgather_1d},
  {"nccl_n4_allgather_2d", nccl_n4_allgather_2d},
  {"nccl_n2_allgather_1d_host", nccl_n2_allgather_1d_host},
  {"nccl_n4_allgather_2d_host", nccl_n4_allgather_2d_host},
//  {"nccl_n2_array_send_irecv_2d", nccl_n2_array_send_irecv_2d},
//  {"nccl_n2_array_isend_irecv_2d", nccl_n2_array_isend_irecv_2d},
  {"nccl_n2_sync_1d", nccl_n2_sync_1d},
  {"nccl_n4_sync_2d_no_corner", nccl_n4_sync_2d_no_corner },
  {"nccl_n4_sync_2d_use_corner", nccl_n4_sync_2d_use_corner},
  {"nccl_n4_sync_1x1v", nccl_n4_sync_1x1v },
  {"nccl_n1_per_sync_2d", nccl_n1_per_sync_2d },
  {"nccl_n2_per_sync_2d", nccl_n2_per_sync_2d },
  {"nccl_n4_multicomm_2d", nccl_n4_multicomm_2d},
  {"nccl_n4_create_comm_from_ranks_1", nccl_n4_create_comm_from_ranks_1 },
  {"nccl_n4_create_comm_from_ranks_2", nccl_n4_create_comm_from_ranks_2 },
  {"nccl_bcast_1d", nccl_bcast_1d},
  {"nccl_bcast_2d", nccl_bcast_2d},
  {"nccl_bcast_1d_host", nccl_bcast_1d_host},
  {"nccl_bcast_2d_host", nccl_bcast_2d_host},
  {NULL, NULL},
};

#else

// nothing to test if not building with NCCL
TEST_LIST = {
  {NULL, NULL},
};

#endif
