// Completed-operation benchmark. See benchmarks/cuda_mpi/README.md for protocol.
#define _GNU_SOURCE
#include <gkyl_cuda_mpi_comm.h>
#include <gkyl_nccl_comm.h>

#if defined(GKYL_HAVE_NCCL) && defined(GKYL_HAVE_CUDA)
#include <gkyl_array_ops.h>
#include <gkyl_alloc.h>
#include <sched.h>
#include <string.h>
#include <unistd.h>

struct bench {
  int rank, size, partner, samples, launch, only_backend;
  const char *placement;
  ncclComm_t raw_nccl;
  cudaStream_t stream;
  struct gkyl_comm *comm[2];
  struct gkyl_range local, ext;
  struct gkyl_array *src, *dst, *host;
};

static void bench_nccl(ncclResult_t status)
{
  if (status != ncclSuccess) {
    fprintf(stderr, "NCCL failure: %s\n", ncclGetErrorString(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

static void bench_mpi(int status)
{
  if (status != MPI_SUCCESS) {
    MPI_Abort(MPI_COMM_WORLD, status);
  }
}

static void bench_call(struct bench *bench, int backend, int operation, int count)
{
  if (operation == 0) {
    if (backend == 0) {
      bench_mpi(MPI_Sendrecv(
        bench->src->data, count, MPI_DOUBLE, bench->partner, 31, bench->dst->data, count,
        MPI_DOUBLE, bench->partner, 31, MPI_COMM_WORLD, MPI_STATUS_IGNORE
      ));
    } else {
      bench_nccl(ncclGroupStart());
      bench_nccl(ncclRecv(
        bench->dst->data, count, ncclDouble, bench->partner, bench->raw_nccl, bench->stream
      ));
      bench_nccl(ncclSend(
        bench->src->data, count, ncclDouble, bench->partner, bench->raw_nccl, bench->stream
      ));
      bench_nccl(ncclGroupEnd());
    }
  } else if (operation == 1) {
    bench_mpi(gkyl_comm_allreduce(
      bench->comm[backend], GKYL_DOUBLE, GKYL_SUM, count, bench->src->data, bench->dst->data
    ));
  } else {
    bench_mpi(gkyl_comm_array_sync(bench->comm[backend], &bench->local, &bench->ext, bench->src));
  }
  // Include completion, including halo unpack kernels, in BOTH backends.
  checkCuda(cudaDeviceSynchronize());
}

static void bench_verify(struct bench *bench, int operation)
{
  gkyl_array_copy(bench->host, operation == 2 ? bench->src : bench->dst);
  int failed = 0;
  if (operation != 2) {
    double expected = operation == 0 ? bench->partner + 1.0 : bench->size * (bench->size + 1.0) / 2;
    const double *data = bench->host->data;
    for (size_t idx = 0; idx < bench->host->size; ++idx) {
      failed |= data[idx] != expected;
    }
  } else {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &bench->ext);
    while (gkyl_range_iter_next(&iter)) {
      int y = iter.idx[1];
      double expected = y < 1 || y > 4 * bench->size ? -1 : (y - 1) / 4 + 1;
      const double *data = gkyl_array_cfetch(bench->host, gkyl_range_idx(&bench->ext, iter.idx));
      for (int comp = 0; comp < 8; ++comp) {
        failed |= data[comp] != expected;
      }
    }
  }
  int any_failed;
  bench_mpi(MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
  if (any_failed) {
    fprintf(stderr, "Validation failed: rank %d operation %d\n", bench->rank, operation);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
}

static void bench_reset(struct bench *bench, int operation)
{
  // Prevent a correct previous backend/sample from masking a missing transfer.
  if (operation == 2) {
    gkyl_array_clear(bench->host, -1);
    gkyl_array_clear_range(bench->host, bench->rank + 1, &bench->local);
    gkyl_array_copy(bench->src, bench->host);
  } else {
    gkyl_array_clear(bench->dst, -1);
  }
  checkCuda(cudaDeviceSynchronize());
}

static void bench_size(struct bench *bench, int operation, size_t bytes)
{
  int count = bytes / sizeof(double);
  struct gkyl_rect_decomp *decomp = 0;
  if (operation == 2) {
    // y is the fastest varying coordinate: its faces require strided packing.
    decomp = gkyl_rect_decomp_new_from_cuts_and_cells(
      2, (int[]){1, bench->size}, (int[]){bytes / 64, 4 * bench->size}
    );
    gkyl_create_ranges(&decomp->ranges[bench->rank], (int[]){0, 1}, &bench->ext, &bench->local);
    bench->src = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8, bench->ext.volume);
    bench->host = gkyl_array_new(GKYL_DOUBLE, 8, bench->ext.volume);
    gkyl_array_clear(bench->host, -1);
    gkyl_array_clear_range(bench->host, bench->rank + 1, &bench->local);
    gkyl_array_copy(bench->src, bench->host);
    bench->dst = 0;
  } else {
    bench->src = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, count);
    bench->dst = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, count);
    bench->host = gkyl_array_new(GKYL_DOUBLE, 1, count);
    gkyl_array_clear(bench->src, bench->rank + 1);
    gkyl_array_clear(bench->dst, -1);
  }
  bench->comm[0] = bench->only_backend == 1 ?
                     0 :
                     gkyl_cuda_mpi_comm_new(&(
                       struct gkyl_cuda_mpi_comm_inp
                     ){.mpi_comm = MPI_COMM_WORLD, .decomp = decomp, .device_set = true});
  bench->comm[1] =
    bench->only_backend == 0 ?
      0 :
      gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){
        .mpi_comm = MPI_COMM_WORLD, .decomp = decomp, .device_set = true, .custream = bench->stream
      });
  int iterations = bytes <= 8192 ? 500 : bytes <= 1048576 ? 200 : bytes <= 8388608 ? 60 : 20;
  int warmup = 20;
  const char *operations[] = {"duplex", "allreduce", "halo"};
  const char *backends[] = {"cuda_mpi", "nccl"};
  for (int sample = 0; sample < bench->samples; ++sample) {
    for (int position = 0; position < 2; ++position) {
      int backend = (sample + bench->launch + position) % 2;
      if (bench->only_backend >= 0 && backend != bench->only_backend) {
        continue;
      }
      bench_reset(bench, operation);
      for (int iter = 0; iter < warmup; ++iter) {
        bench_call(bench, backend, operation, count);
      }
      bench_verify(bench, operation);
      bench_reset(bench, operation);
      bench_mpi(MPI_Barrier(MPI_COMM_WORLD));
      double start = MPI_Wtime();
      for (int iter = 0; iter < iterations; ++iter) {
        bench_call(bench, backend, operation, count);
      }
      double elapsed = (MPI_Wtime() - start) / iterations;
      double maximum, minimum, sum;
      bench_mpi(MPI_Reduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
      bench_mpi(MPI_Reduce(&elapsed, &minimum, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD));
      bench_mpi(MPI_Reduce(&elapsed, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
      bench_verify(bench, operation);
      if (bench->rank == 0) {
        printf(
          "%s,%d,%d,%s,%s,%zu,%d,%d,%d,%d,%.12g,%.12g,%.12g\n", bench->placement, bench->size,
          bench->launch, operations[operation], backends[backend], bytes, sample, position,
          iterations, warmup, maximum, minimum, sum / bench->size
        );
        fflush(stdout);
      }
    }
  }
  gkyl_comm_release(bench->comm[0]);
  gkyl_comm_release(bench->comm[1]);
  if (decomp) {
    gkyl_rect_decomp_release(decomp);
  }
  gkyl_array_release(bench->src);
  if (bench->dst) {
    gkyl_array_release(bench->dst);
  }
  gkyl_array_release(bench->host);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  struct bench bench = {.samples = 9, .launch = 0, .only_backend = -1, .placement = "unspecified"};
  MPI_Comm_rank(MPI_COMM_WORLD, &bench.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &bench.size);
  int distance = 1, max_power = 26;
  for (int idx = 1; idx < argc; ++idx) {
    if (!strcmp(argv[idx], "--backend") && idx + 1 < argc) {
      const char *backend = argv[++idx];
      if (!strcmp(backend, "cuda_mpi")) {
        bench.only_backend = 0;
      } else if (!strcmp(backend, "nccl")) {
        bench.only_backend = 1;
      } else {
        MPI_Abort(MPI_COMM_WORLD, 3);
      }
    } else if (!strcmp(argv[idx], "--placement") && idx + 1 < argc) {
      bench.placement = argv[++idx];
    } else if (!strcmp(argv[idx], "--launch") && idx + 1 < argc) {
      bench.launch = atoi(argv[++idx]);
    } else if (!strcmp(argv[idx], "--samples") && idx + 1 < argc) {
      bench.samples = atoi(argv[++idx]);
    } else if (!strcmp(argv[idx], "--distance") && idx + 1 < argc) {
      distance = atoi(argv[++idx]);
    } else if (!strcmp(argv[idx], "--max-power") && idx + 1 < argc) {
      max_power = atoi(argv[++idx]);
    } else {
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
  }
  bench.partner = bench.rank ^ distance;
  if (
    bench.size < 2 || bench.partner >= bench.size || bench.samples < 1 || max_power > 26 ||
    max_power < 6
  ) {
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  MPI_Comm shared;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, bench.rank, MPI_INFO_NULL, &shared);
  int local_rank, visible;
  MPI_Comm_rank(shared, &local_rank);
  MPI_Comm_free(&shared);
  checkCuda(cudaGetDeviceCount(&visible));
  checkCuda(cudaSetDevice(local_rank % visible));
  int device;
  checkCuda(cudaGetDevice(&device));
  char hostname[256], bus[32];
  gethostname(hostname, sizeof hostname);
  checkCuda(cudaDeviceGetPCIBusId(bus, sizeof bus, device));
  fprintf(
    stderr, "rank=%d host=%s local_rank=%d device=%d bus=%s cpu=%d visible=%s\n", bench.rank,
    hostname, local_rank, device, bus, sched_getcpu(), getenv("CUDA_VISIBLE_DEVICES")
  );
  if (bench.rank == 0) {
    int runtime, driver, nccl_version, len;
    cudaRuntimeGetVersion(&runtime);
    cudaDriverGetVersion(&driver);
    ncclGetVersion(&nccl_version);
    char library[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Get_library_version(library, &len);
    fprintf(
      stderr, "CUDA runtime=%d driver=%d NCCL=%d MPI=%s\n", runtime, driver, nccl_version, library
    );
    fprintf(
      stderr, "MPICH_GPU_SUPPORT_ENABLED=%s NCCL_NET=%s NCCL_CROSS_NIC=%s\n",
      getenv("MPICH_GPU_SUPPORT_ENABLED"), getenv("NCCL_NET"), getenv("NCCL_CROSS_NIC")
    );
    printf(
      "placement,ranks,launch,operation,backend,bytes,sample,order,iterations,warmup,max_seconds,min_seconds,mean_seconds\n"
    );
  }
  checkCuda(cudaStreamCreate(&bench.stream));
  ncclUniqueId id;
  if (bench.only_backend != 0) {
    if (bench.rank == 0) {
      bench_nccl(ncclGetUniqueId(&id));
    }
    MPI_Bcast(&id, sizeof id, MPI_BYTE, 0, MPI_COMM_WORLD);
    bench_nccl(ncclCommInitRank(&bench.raw_nccl, bench.size, id, bench.rank));
  }
  for (int operation = 0; operation < 3; ++operation) {
    for (int power = operation == 2 ? 6 : 3; power <= max_power;
         power = power == 3 ? 6 : power + 2) {
      bench_size(&bench, operation, (size_t)1 << power);
    }
  }
  if (bench.only_backend != 0) {
    bench_nccl(ncclCommDestroy(bench.raw_nccl));
  }
  checkCuda(cudaStreamDestroy(bench.stream));
  MPI_Finalize();
  return 0;
}
#else
int main(void)
{
  return 77;
}
#endif
