// Verify CPU vs GPU inference consistency across network sizes.
// Trains ONE model on CPU, then runs inference on both CPU and GPU
// and compares the outputs. Also trains on GPU and verifies the
// GPU-trained model gives reasonable results on CPU inference.
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <stdlib.h>
#include <math.h>

static inline float
ufunc(float x)
{
  return 1.0f/(1.0f+100.0f*x*x);
}

static void
train_and_save(int ntrain, int nwidth, int ndepth, bool use_gpu,
  const char *filename)
{
  kad_node_t *t_net = kann_layer_input(1);
  for (int i = 0; i < ndepth; ++i) {
    t_net = kann_layer_dense(t_net, nwidth);
    t_net = kad_tanh(t_net);
  }
  t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
  struct gkyl_kann_net *net = gkyl_kann_net_new(t_net, use_gpu);

  int N = ntrain;
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);

  float dx = 2.0f / (N - 1);
  for (int i = 0; i < N; ++i) {
    inp->vals[i][0] = -1.0f + dx * i;
    out->vals[i][0] = ufunc(inp->vals[i][0]);
  }

  struct gkyl_kn_vec *inp_t = inp, *out_t = out;
  struct gkyl_kn_vec *inp_cu = 0, *out_cu = 0;
  if (use_gpu) {
    inp_cu = gkyl_kn_vec_cu_dev_new(N, 1);
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

  gkyl_kann_net_train_fnn1(net, &params, inp_t, out_t);
  gkyl_kann_net_save(net, filename);

  if (use_gpu) {
    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  }
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  gkyl_kann_net_release(net);
}

// Run inference and return outputs in out (always host kn_vec)
static void
infer(const char *filename, bool use_gpu,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(filename, use_gpu);

  if (use_gpu) {
    int nvec = inp->nvec;
    struct gkyl_kn_vec *inp_cu = gkyl_kn_vec_cu_dev_new(nvec, inp->N);
    struct gkyl_kn_vec *out_cu = gkyl_kn_vec_cu_dev_new(nvec, out->N);
    gkyl_kn_vec_copy(inp_cu, inp);
    gkyl_kann_net_apply(net, inp_cu, out_cu);
    gkyl_kn_vec_copy(out, out_cu);
    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  } else {
    gkyl_kann_net_apply(net, inp, out);
  }

  gkyl_kann_net_release(net);
}

int
main(int argc, char *argv[])
{
  int ntrain = 1001;
  int ninfer = 21;

  int widths[] = { 64, 128, 256, 512 };
  int depths[] = { 2, 4 };
  int nw = sizeof(widths) / sizeof(widths[0]);
  int nd = sizeof(depths) / sizeof(depths[0]);

  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(ninfer, 1);
  struct gkyl_kn_vec *out_cpu = gkyl_kn_vec_new(ninfer, 1);
  struct gkyl_kn_vec *out_gpu = gkyl_kn_vec_new(ninfer, 1);

  float dx = 2.0f / (ninfer - 1);
  for (int i = 0; i < ninfer; ++i)
    inp->vals[i][0] = -1.0f + dx * i;

  fprintf(stdout, "=== Test 1: Train on CPU, compare CPU vs GPU inference ===\n");
  fprintf(stdout, "  (Same model, same weights — should match to float precision)\n\n");

  for (int di = 0; di < nd; ++di) {
    for (int wi = 0; wi < nw; ++wi) {
      int w = widths[wi], d = depths[di];

      train_and_save(ntrain, w, d, false, "verify.kann");

      infer("verify.kann", false, inp, out_cpu);
      infer("verify.kann", true, inp, out_gpu);

      float max_diff = 0;
      for (int i = 0; i < ninfer; ++i) {
        float diff = fabsf(out_cpu->vals[i][0] - out_gpu->vals[i][0]);
        if (diff > max_diff) max_diff = diff;
      }

      const char *status = max_diff < 1e-4f ? "PASS" : "FAIL";
      fprintf(stdout, "  width=%3d depth=%d: max|cpu-gpu|=%.3e  [%s]\n",
        w, d, max_diff, status);
    }
  }

  fprintf(stdout, "\n=== Test 2: Train on GPU, verify model quality via CPU inference ===\n");
  fprintf(stdout, "  (GPU-trained model loaded on CPU — should approximate f(x) well)\n\n");

  for (int di = 0; di < nd; ++di) {
    for (int wi = 0; wi < nw; ++wi) {
      int w = widths[wi], d = depths[di];

      train_and_save(ntrain, w, d, true, "verify.kann");

      infer("verify.kann", false, inp, out_cpu);

      float mse = 0;
      for (int i = 0; i < ninfer; ++i) {
        float exact = ufunc(inp->vals[i][0]);
        float diff = out_cpu->vals[i][0] - exact;
        mse += diff * diff;
      }
      mse /= ninfer;

      float val_at_0 = out_cpu->vals[ninfer/2][0];

      const char *status = (mse < 0.05f && fabsf(val_at_0 - 1.0f) < 0.2f)
        ? "PASS" : "FAIL";
      fprintf(stdout, "  width=%3d depth=%d: mse=%.3e  f(0)=%.4f (exact=1.0)  [%s]\n",
        w, d, mse, val_at_0, status);
    }
  }

  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out_cpu);
  gkyl_kn_vec_release(out_gpu);

  remove("verify.kann");

  return 0;
}
