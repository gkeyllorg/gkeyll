// Benchmark GPU vs CPU training and inference for different network sizes
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

static inline float
ufunc(float x)
{
  return 1.0f/(1.0f+100.0f*x*x);
}

static double
bench_train(int ntrain, int nwidth, int ndepth, bool use_gpu)
{
  kad_node_t *t_net = kann_layer_input(1);
  for (int i = 0; i < ndepth; ++i) {
    t_net = kann_layer_dense(t_net, nwidth);
    t_net = kad_tanh(t_net);
  }
  t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
  struct gkyl_kann_net *net = gkyl_kann_net_new(t_net, use_gpu);

  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(ntrain, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(ntrain, 1);

  float dx = 2.0f / (ntrain - 1);
  for (int i = 0; i < ntrain; ++i) {
    float x = -1.0f + dx * i;
    inp->vals[i][0] = x;
    out->vals[i][0] = ufunc(x);
  }

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(ntrain, 1);
    out_cu = gkyl_kn_vec_cu_dev_new(ntrain, 1);
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
bench_infer(int nvec, int nwidth, int ndepth, bool use_gpu,
  const char *model_file)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(model_file, use_gpu);

  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

  float dx = 2.0f / (nvec - 1);
  for (int i = 0; i < nvec; ++i)
    inp->vals[i][0] = -1.0f + dx * i;

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(nvec, 1);
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

int
main(int argc, char *argv[])
{
  int ntrain = 1001;
  int ninfer = 1001;

  int widths[] = { 64, 128, 256, 512 };
  int depths[] = { 2, 4 };
  int nw = sizeof(widths) / sizeof(widths[0]);
  int nd = sizeof(depths) / sizeof(depths[0]);

  fprintf(stdout, "=== Training Benchmark (ntrain=%d, 50 epochs, mini_size=64) ===\n", ntrain);
  fprintf(stdout, "%6s %6s %10s %10s %10s\n",
    "width", "depth", "CPU (s)", "GPU (s)", "speedup");

  for (int di = 0; di < nd; ++di) {
    for (int wi = 0; wi < nw; ++wi) {
      int w = widths[wi], d = depths[di];

      double t_cpu = bench_train(ntrain, w, d, false);
      double t_gpu = bench_train(ntrain, w, d, true);

      fprintf(stdout, "%6d %6d %10.4f %10.4f %10.2fx\n",
        w, d, t_cpu, t_gpu, t_cpu / t_gpu);
    }
  }

  fprintf(stdout, "\n=== Inference Benchmark (nvec=%d, avg of 100 reps) ===\n", ninfer);
  fprintf(stdout, "%6s %6s %10s %10s %10s\n",
    "width", "depth", "CPU (ms)", "GPU (ms)", "speedup");

  for (int di = 0; di < nd; ++di) {
    for (int wi = 0; wi < nw; ++wi) {
      int w = widths[wi], d = depths[di];

      // Train a model to benchmark inference with
      kad_node_t *t_net = kann_layer_input(1);
      for (int i = 0; i < d; ++i) {
        t_net = kann_layer_dense(t_net, w);
        t_net = kad_tanh(t_net);
      }
      t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
      struct gkyl_kann_net *net = gkyl_kann_net_new(t_net, false);

      struct gkyl_kn_vec *inp = gkyl_kn_vec_new(ntrain, 1);
      struct gkyl_kn_vec *out = gkyl_kn_vec_new(ntrain, 1);
      float dx = 2.0f / (ntrain - 1);
      for (int i = 0; i < ntrain; ++i) {
        inp->vals[i][0] = -1.0f + dx * i;
        out->vals[i][0] = ufunc(inp->vals[i][0]);
      }
      struct gkyl_kann_train_params params = {
        .learning_rate = 1e-3f, .mini_size = 64,
        .max_epoch = 5, .max_drop_streak = 5, .frac_val = 0.1f,
      };
      gkyl_kann_net_train_fnn1(net, &params, inp, out);
      gkyl_kann_net_save(net, "bench_tmp.kann");
      gkyl_kn_vec_release(inp);
      gkyl_kn_vec_release(out);
      gkyl_kann_net_release(net);

      double t_cpu = bench_infer(ninfer, w, d, false, "bench_tmp.kann");
      double t_gpu = bench_infer(ninfer, w, d, true, "bench_tmp.kann");

      fprintf(stdout, "%6d %6d %10.4f %10.4f %10.2fx\n",
        w, d, t_cpu * 1000, t_gpu * 1000, t_cpu / t_gpu);
    }
  }

  remove("bench_tmp.kann");
  return 0;
}
