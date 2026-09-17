// Benchmark GPU vs CPU training and inference for RNN architectures
// Compares MLP, GRU, and GRU+Norm (stdnorm) across different sizes,
// including sequential RNN inference with pre-recurrence.
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

// 2D function to fit: sin(2*pi*x)*sin(2*t)*exp(-t)
static inline float
ufunc(float t, float x)
{
  return sinf(2.0f*(float)M_PI*x)*sinf(2.0f*t)*expf(-t);
}

enum arch_type {
  ARCH_MLP,
  ARCH_GRU,
  ARCH_GRU_NORM
};

static const char *arch_name[] = { "MLP", "GRU", "GRU+Norm" };

static kad_node_t*
build_net(enum arch_type arch, int nwidth, int ndepth)
{
  kad_node_t *t = kann_layer_input(2);
  for (int i = 0; i < ndepth; ++i) {
    switch (arch) {
      case ARCH_MLP:
        t = kann_layer_dense(t, nwidth);
        break;
      case ARCH_GRU:
        t = kann_layer_gru(t, nwidth, 0);
        break;
      case ARCH_GRU_NORM:
        t = kann_layer_gru(t, nwidth, KANN_RNN_NORM);
        break;
    }
    t = kad_tanh(t);
  }
  t = kann_layer_cost(t, 1, KANN_C_MSE);
  return t;
}

static void
fill_2d_data(struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out,
  int Nt, int Nx)
{
  float dt = 3.0f / (Nt - 1);
  float dx = 1.0f / (Nx - 1);
  for (int i = 0; i < Nt; ++i)
    for (int j = 0; j < Nx; ++j) {
      int idx = i * Nx + j;
      float t = dt * i;
      float x = dx * j;
      inp->vals[idx][0] = t;
      inp->vals[idx][1] = x;
      out->vals[idx][0] = ufunc(t, x);
    }
}

static double
bench_train(enum arch_type arch, int ntrain_t, int ntrain_x,
  int nwidth, int ndepth, bool use_gpu, const char *save_file)
{
  kad_node_t *cost = build_net(arch, nwidth, ndepth);
  struct gkyl_kann_net *net = gkyl_kann_net_new(cost, use_gpu);

  int N = ntrain_t * ntrain_x;
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 2);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);
  fill_2d_data(inp, out, ntrain_t, ntrain_x);

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(N, 2);
    out_cu = gkyl_kn_vec_cu_dev_new(N, 1);
    gkyl_kn_vec_copy(inp_cu, inp);
    gkyl_kn_vec_copy(out_cu, out);
    inp_t = inp_cu;
    out_t = out_cu;
  }

  struct gkyl_kann_train_params params = {
    .learning_rate = 1e-3f,
    .mini_size = 64,
    .max_epoch = 50,
    .max_drop_streak = 10,
    .frac_val = 0.1f,
  };

  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  gkyl_kann_net_train_fnn1(net, &params, inp_t, out_t);
  clock_gettime(CLOCK_MONOTONIC, &t1);

  double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;

  if (save_file)
    gkyl_kann_net_save(net, save_file);

  if (use_gpu) {
    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  }
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  gkyl_kann_net_release(net);

  return elapsed;
}

static double
bench_infer_batch(const char *model_file, int nvec, bool use_gpu)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(model_file, use_gpu);

  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 2);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

  float dt = 3.0f / (nvec - 1);
  for (int i = 0; i < nvec; ++i) {
    inp->vals[i][0] = dt * i;
    inp->vals[i][1] = 0.35f;
  }

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(nvec, 2);
    out_cu = gkyl_kn_vec_cu_dev_new(nvec, 1);
    gkyl_kn_vec_copy(inp_cu, inp);
    inp_t = inp_cu;
    out_t = out_cu;
  }

  // Warmup
  gkyl_kann_net_apply(net, inp_t, out_t);

  int nreps = 100;
  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  for (int r = 0; r < nreps; ++r)
    gkyl_kann_net_apply(net, inp_t, out_t);
  clock_gettime(CLOCK_MONOTONIC, &t1);

  double elapsed = ((t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9) / nreps;

  if (use_gpu) {
    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  }
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  gkyl_kann_net_release(net);

  return elapsed;
}

static double
bench_infer_rnn(const char *model_file, int nvec, bool use_gpu)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(model_file, use_gpu);

  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 2);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

  float dt = 3.0f / (nvec - 1);
  for (int i = 0; i < nvec; ++i) {
    inp->vals[i][0] = dt * i;
    inp->vals[i][1] = 0.35f;
  }

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(nvec, 2);
    out_cu = gkyl_kn_vec_cu_dev_new(nvec, 1);
    gkyl_kn_vec_copy(inp_cu, inp);
    inp_t = inp_cu;
    out_t = out_cu;
  }

  // Warmup
  gkyl_kann_net_apply_rnn(net, inp_t, out_t);

  int nreps = 100;
  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  for (int r = 0; r < nreps; ++r)
    gkyl_kann_net_apply_rnn(net, inp_t, out_t);
  clock_gettime(CLOCK_MONOTONIC, &t1);

  double elapsed = ((t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9) / nreps;

  if (use_gpu) {
    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  }
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  gkyl_kann_net_release(net);

  return elapsed;
}

int
main(int argc, char *argv[])
{
  int ntrain_t = 51, ntrain_x = 51;
  int ninfer = 101;

  int widths[] = { 32, 64, 128 };
  int depths[] = { 2, 4 };
  int nw = sizeof(widths) / sizeof(widths[0]);
  int nd = sizeof(depths) / sizeof(depths[0]);

  enum arch_type archs[] = { ARCH_MLP, ARCH_GRU, ARCH_GRU_NORM };
  int na = sizeof(archs) / sizeof(archs[0]);

  // ---- Training benchmark ----
  fprintf(stdout, "=== Training Benchmark (ntrain=%dx%d=%d, 50 epochs, mini_size=64) ===\n",
    ntrain_t, ntrain_x, ntrain_t * ntrain_x);
  fprintf(stdout, "%10s %6s %6s %10s %10s %10s\n",
    "arch", "width", "depth", "CPU (s)", "GPU (s)", "speedup");

  for (int ai = 0; ai < na; ++ai) {
    for (int di = 0; di < nd; ++di) {
      for (int wi = 0; wi < nw; ++wi) {
        int w = widths[wi], d = depths[di];

        double t_cpu = bench_train(archs[ai], ntrain_t, ntrain_x,
          w, d, false, NULL);
        double t_gpu = bench_train(archs[ai], ntrain_t, ntrain_x,
          w, d, true, NULL);

        fprintf(stdout, "%10s %6d %6d %10.4f %10.4f %10.2fx\n",
          arch_name[archs[ai]], w, d, t_cpu, t_gpu, t_cpu / t_gpu);
      }
    }
  }

  // ---- Batch inference benchmark ----
  fprintf(stdout, "\n=== Batch Inference Benchmark (nvec=%d, avg of 100 reps) ===\n", ninfer);
  fprintf(stdout, "%10s %6s %6s %10s %10s %10s\n",
    "arch", "width", "depth", "CPU (ms)", "GPU (ms)", "speedup");

  for (int ai = 0; ai < na; ++ai) {
    for (int di = 0; di < nd; ++di) {
      for (int wi = 0; wi < nw; ++wi) {
        int w = widths[wi], d = depths[di];

        // Quick train to get a model file
        bench_train(archs[ai], ntrain_t, ntrain_x, w, d, false,
          "bench_cmp_tmp.kann");

        double t_cpu = bench_infer_batch("bench_cmp_tmp.kann", ninfer, false);
        double t_gpu = bench_infer_batch("bench_cmp_tmp.kann", ninfer, true);

        fprintf(stdout, "%10s %6d %6d %10.4f %10.4f %10.2fx\n",
          arch_name[archs[ai]], w, d, t_cpu * 1000, t_gpu * 1000,
          t_cpu / t_gpu);

        remove("bench_cmp_tmp.kann");
      }
    }
  }

  // ---- Sequential RNN inference benchmark (GRU and GRU+Norm only) ----
  fprintf(stdout, "\n=== Sequential RNN Inference Benchmark (nvec=%d, avg of 100 reps) ===\n", ninfer);
  fprintf(stdout, "%10s %6s %6s %10s %10s %10s\n",
    "arch", "width", "depth", "CPU (ms)", "GPU (ms)", "speedup");

  enum arch_type rnn_archs[] = { ARCH_GRU, ARCH_GRU_NORM };
  int nra = sizeof(rnn_archs) / sizeof(rnn_archs[0]);

  for (int ai = 0; ai < nra; ++ai) {
    for (int di = 0; di < nd; ++di) {
      for (int wi = 0; wi < nw; ++wi) {
        int w = widths[wi], d = depths[di];

        bench_train(rnn_archs[ai], ntrain_t, ntrain_x, w, d, false,
          "bench_cmp_tmp.kann");

        double t_cpu = bench_infer_rnn("bench_cmp_tmp.kann", ninfer, false);
        double t_gpu = bench_infer_rnn("bench_cmp_tmp.kann", ninfer, true);

        fprintf(stdout, "%10s %6d %6d %10.4f %10.4f %10.2fx\n",
          arch_name[rnn_archs[ai]], w, d, t_cpu * 1000, t_gpu * 1000,
          t_cpu / t_gpu);

        remove("bench_cmp_tmp.kann");
      }
    }
  }

  return 0;
}
