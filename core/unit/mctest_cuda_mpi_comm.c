#include <acutest.h>
#include <gkyl_cuda_mpi_comm.h>

#if defined(GKYL_HAVE_MPI) && defined(GKYL_HAVE_CUDA)
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_multib_comm_conn.h>
#include <string.h>

static struct gkyl_comm *new_comm(const struct gkyl_rect_decomp *decomp, bool corners)
{
  return gkyl_cuda_mpi_comm_new(&(
    struct gkyl_cuda_mpi_comm_inp
  ){.mpi_comm = MPI_COMM_WORLD, .decomp = decomp, .sync_corners = corners});
}

static void test_collectives(void)
{
  struct gkyl_comm *comm = new_comm(0, false);
  int rank, size;
  gkyl_comm_get_rank(comm, &rank);
  gkyl_comm_get_size(comm, &size);
  TEST_CHECK(strcmp(comm->id, "cuda_mpi_comm") == 0);
  TEST_CHECK(!comm->has_decomp);
  enum gkyl_elem_type types[] = {GKYL_INT, GKYL_INT_64, GKYL_FLOAT, GKYL_DOUBLE, GKYL_LONG};
  for (int idx = 0; idx < 5; ++idx) {
    // Raw reduction buffers avoid the separate array allocator's lack of an
    // element size for GKYL_INT_64.
    union {
      int integer;
      int64_t integer64;
      float single;
      double real;
      long long_integer;
    } host = {0};
    if (idx == 0) {
      host.integer = rank + 1;
    }
    if (idx == 1) {
      host.integer64 = rank + 1;
    }
    if (idx == 2) {
      host.single = rank + 1;
    }
    if (idx == 3) {
      host.real = rank + 1;
    }
    if (idx == 4) {
      host.long_integer = rank + 1;
    }
    void *gpu = gkyl_cu_malloc(sizeof host);
    gkyl_cu_memcpy(gpu, &host, sizeof host, GKYL_CU_MEMCPY_H2D);
    TEST_CHECK(gkyl_comm_allreduce(comm, types[idx], GKYL_SUM, 1, gpu, gpu) == 0);
    gkyl_cu_memcpy(&host, gpu, sizeof host, GKYL_CU_MEMCPY_D2H);
    double result = idx == 0 ? host.integer :
                    idx == 1 ? host.integer64 :
                    idx == 2 ? host.single :
                    idx == 3 ? host.real :
                               host.long_integer;
    TEST_CHECK(result == size * (size + 1.0) / 2);
    gkyl_cu_free(gpu);
  }
  struct gkyl_array *send = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 128);
  struct gkyl_array *recv = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 128);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 1, 128);
  gkyl_array_clear(send, rank + 1);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_MIN, 128, send->data, recv->data);
  gkyl_array_copy(host, recv);
  for (int idx = 0; idx < 128; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == 1);
  }
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_MAX, 128, send->data, recv->data);
  gkyl_array_copy(host, recv);
  for (int idx = 0; idx < 128; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == size);
  }

  // Work submitted on a nonblocking stream must finish before MPI reads it.
  cudaStream_t producer;
  checkCuda(cudaStreamCreateWithFlags(&producer, cudaStreamNonBlocking));
  checkCuda(cudaMemsetAsync(send->data, 0, send->esznc * send->size, producer));
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, 128, send->data, recv->data);
  gkyl_array_copy(host, recv);
  for (int idx = 0; idx < 128; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == 0);
  }
  checkCuda(cudaStreamDestroy(producer));

  gkyl_array_clear(send, 17);
  gkyl_array_clear(recv, -1);
  TEST_CHECK(gkyl_comm_array_bcast(comm, rank == size - 1 ? send : 0, recv, size - 1) == 0);
  gkyl_array_copy(host, recv);
  for (int idx = 0; idx < 128; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == 17);
  }
  struct gkyl_array *host_send = gkyl_array_new(GKYL_DOUBLE, 1, 128);
  gkyl_array_clear(host_send, 23);
  gkyl_comm_array_bcast_host(comm, rank == 0 ? host_send : 0, host, 0);
  for (int idx = 0; idx < 128; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == 23);
  }
  double value = rank + 1;
  TEST_CHECK(gkyl_comm_allreduce_host(comm, GKYL_USER, GKYL_SUM, 1, &value, &value) == MPI_ERR_TYPE);
  gkyl_comm_allreduce_host(comm, GKYL_DOUBLE, GKYL_SUM, 1, &value, &value);
  TEST_CHECK(value == size * (size + 1.0) / 2);
  gkyl_comm_group_call_start(comm);
  gkyl_comm_group_call_end(comm);
  gkyl_array_release(send);
  gkyl_array_release(recv);
  gkyl_array_release(host);
  gkyl_array_release(host_send);
  gkyl_comm_release(comm);
}

static double coordinate_value(int x, int y, int comp)
{
  return x + 1000.0 * y + 1000000.0 * comp;
}

static void fill_coordinates(struct gkyl_array *host, const struct gkyl_range *local)
{
  gkyl_array_clear(host, -1);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    double *value = gkyl_array_fetch(host, gkyl_range_idx(local, iter.idx));
    for (int comp = 0; comp < host->ncomp; ++comp) {
      value[comp] = coordinate_value(iter.idx[0], iter.idx[1], comp);
    }
  }
}

static void test_gather(void)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts_and_cells(2, (int[]){size, 1}, (int[]){3 * size + 1, 5});
  struct gkyl_comm *comm = new_comm(decomp, false);
  struct gkyl_range local, ext, global, global_ext;
  gkyl_create_ranges(&decomp->ranges[rank], (int[]){1, 1}, &ext, &local);
  gkyl_create_ranges(&decomp->parent_range, (int[]){1, 1}, &global_ext, &global);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 3, ext.volume);
  struct gkyl_array *src = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, ext.volume);
  struct gkyl_array *dst = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, global_ext.volume);
  struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 3, global_ext.volume);
  fill_coordinates(host, &local);
  gkyl_array_copy(src, host);
  for (int device = 0; device < 2; ++device) {
    gkyl_array_clear(out, -1);
    gkyl_array_copy(dst, out);
    int status = device ? gkyl_comm_array_allgather(comm, &local, &global, src, dst) :
                          gkyl_comm_array_allgather_host(comm, &local, &global, host, out);
    TEST_CHECK(status == 0);
    if (device) {
      gkyl_array_copy(out, dst);
    }
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &global_ext);
    while (gkyl_range_iter_next(&iter)) {
      double *value = gkyl_array_fetch(out, gkyl_range_idx(&global_ext, iter.idx));
      bool inside = gkyl_range_contains_idx(&global, iter.idx);
      for (int comp = 0; comp < 3; ++comp) {
        TEST_CHECK(value[comp] == (inside ? coordinate_value(iter.idx[0], iter.idx[1], comp) : -1));
      }
    }
  }
  gkyl_array_release(host);
  gkyl_array_release(src);
  gkyl_array_release(dst);
  gkyl_array_release(out);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

static void test_halo(void)
{
  int rank, size, cuts[2] = {0, 0};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Dims_create(size, 2, cuts);
  int cells[] = {5 * cuts[0] + 1, 5 * cuts[1] + 1};
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts_and_cells(2, cuts, cells);
  struct gkyl_comm *comm = new_comm(decomp, true);
  struct gkyl_range local, ext;
  gkyl_create_ranges(&decomp->ranges[rank], (int[]){2, 2}, &ext, &local);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 3, ext.volume);
  struct gkyl_array *array = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, ext.volume);
  fill_coordinates(host, &local);
  gkyl_array_copy(array, host);
  TEST_CHECK(gkyl_comm_array_sync(comm, &local, &ext, array) == 0);
  gkyl_array_copy(host, array);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ext);
  while (gkyl_range_iter_next(&iter)) {
    if (!gkyl_range_contains_idx(&decomp->parent_range, iter.idx)) {
      continue;
    }
    double *value = gkyl_array_fetch(host, gkyl_range_idx(&ext, iter.idx));
    for (int comp = 0; comp < 3; ++comp) {
      TEST_CHECK(value[comp] == coordinate_value(iter.idx[0], iter.idx[1], comp));
      TEST_MSG(
        "rank %d index %d,%d component %d: %.0f", rank, iter.idx[0], iter.idx[1], comp, value[comp]
      );
    }
  }
  TEST_CHECK(gkyl_comm_array_per_sync(comm, &local, &ext, 2, (int[]){0, 1}, array) == 0);
  gkyl_array_copy(host, array);
  gkyl_range_iter_init(&iter, &ext);
  while (gkyl_range_iter_next(&iter)) {
    int x = ((iter.idx[0] - 1) % cells[0] + cells[0]) % cells[0] + 1;
    int y = ((iter.idx[1] - 1) % cells[1] + cells[1]) % cells[1] + 1;
    double *value = gkyl_array_fetch(host, gkyl_range_idx(&ext, iter.idx));
    for (int comp = 0; comp < 3; ++comp) {
      TEST_CHECK(value[comp] == coordinate_value(x, y, comp));
      TEST_MSG(
        "periodic rank %d index %d,%d component %d: %.0f expected %.0f", rank, iter.idx[0],
        iter.idx[1], comp, value[comp], coordinate_value(x, y, comp)
      );
    }
  }
  gkyl_array_release(host);
  gkyl_array_release(array);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

static void test_zero_halo(void)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts_and_cells(1, (int[]){size}, (int[]){4 * size});
  struct gkyl_comm *comm = new_comm(decomp, false);
  struct gkyl_range local, ext;
  gkyl_create_ranges(&decomp->ranges[rank], (int[]){0}, &ext, &local);
  struct gkyl_array *array = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, ext.volume);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 1, ext.volume);
  gkyl_array_clear(array, rank + 1);
  TEST_CHECK(gkyl_comm_array_sync(comm, &local, &ext, array) == 0);
  TEST_CHECK(gkyl_comm_array_per_sync(comm, &local, &ext, 1, (int[]){0}, array) == 0);
  gkyl_array_copy(host, array);
  for (int idx = 0; idx < host->size; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == rank + 1);
  }
  gkyl_array_release(array);
  gkyl_array_release(host);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

static void test_lifetime(void)
{
  struct gkyl_comm *comm = new_comm(0, true);
  int rank, size;
  gkyl_comm_get_rank(comm, &rank);
  gkyl_comm_get_size(comm, &size);
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]){1}, (int[]){3});
  struct gkyl_comm *extended = gkyl_comm_extend_comm(comm, &range);
  struct gkyl_comm *split = gkyl_comm_split_comm(comm, rank % 2 ? MPI_UNDEFINED : 0, 0);
  TEST_CHECK((split != 0) == (rank % 2 == 0));
  int ranks[] = {size - 1};
  bool valid;
  struct gkyl_comm *subset = gkyl_comm_create_comm_from_ranks(comm, 1, ranks, 0, &valid);
  TEST_CHECK(valid == (rank == size - 1));
  gkyl_comm_release(comm);
  double value = 1, result;
  gkyl_comm_allreduce_host(extended, GKYL_DOUBLE, GKYL_SUM, 1, &value, &result);
  TEST_CHECK(result == size);
  if (split) {
    gkyl_comm_allreduce_host(split, GKYL_DOUBLE, GKYL_SUM, 1, &value, &result);
    TEST_CHECK(result == (size + 1) / 2);
  }
  if (subset) {
    gkyl_comm_get_rank(subset, &rank);
    TEST_CHECK(rank == 0);
    gkyl_comm_get_size(subset, &size);
    TEST_CHECK(size == 1);
  }
  gkyl_comm_release(extended);
  gkyl_comm_release(split);
  gkyl_comm_release(subset);
}

static void test_multiblock(void)
{
  struct gkyl_comm *comm = new_comm(0, false);
  int rank, size;
  gkyl_comm_get_rank(comm, &rank);
  gkyl_comm_get_size(comm, &size);
  int next = (rank + 1) % size, prev = (rank + size - 1) % size;
  struct gkyl_array *src = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 20);
  struct gkyl_array *dst = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 20);
  struct gkyl_array *host = gkyl_array_new(GKYL_DOUBLE, 1, 20);
  for (int idx = 0; idx < 20; ++idx) {
    ((double *)host->data)[idx] = 100 * rank + idx;
  }
  gkyl_array_copy(src, host);
  gkyl_array_clear(dst, -1);
  struct gkyl_range all;
  gkyl_range_init(&all, 1, (int[]){0}, (int[]){19});
  struct gkyl_comm_conn send[5] = {0}, recv[5] = {0};
  for (int idx = 0; idx < 5; ++idx) {
    gkyl_sub_range_init(&send[idx].range, &all, (int[]){4 * idx}, (int[]){4 * idx + 3});
    send[idx].rank = next;
    send[idx].block_id = rank;
    recv[idx] = send[idx];
    recv[idx].rank = prev;
    recv[idx].block_id = prev;
  }
  struct gkyl_multib_comm_conn *sconn = gkyl_multib_comm_conn_new(5, send);
  struct gkyl_multib_comm_conn *rconn = gkyl_multib_comm_conn_new(5, recv);
  TEST_CHECK(gkyl_multib_comm_conn_array_transfer(comm, 1, &rank, &sconn, &rconn, &src, &dst) == 0);
  gkyl_array_copy(host, dst);
  for (int idx = 0; idx < 20; ++idx) {
    TEST_CHECK(((double *)host->data)[idx] == 100 * prev + idx);
  }
  gkyl_multib_comm_conn_release(sconn);
  gkyl_multib_comm_conn_release(rconn);
  gkyl_array_release(src);
  gkyl_array_release(dst);
  gkyl_array_release(host);
  gkyl_comm_release(comm);
}

TEST_LIST = {
  {"cuda_mpi_collectives", test_collectives},
  {"cuda_mpi_uneven_gather", test_gather},
  {"cuda_mpi_halo_corners_periodic", test_halo},
  {"cuda_mpi_zero_width_halo", test_zero_halo},
  {"cuda_mpi_lifetime", test_lifetime},
  {"cuda_mpi_multiblock", test_multiblock},
  {NULL, NULL}
};
#else
TEST_LIST = {{NULL, NULL}};
#endif
