// Compare MLP, GRU and BEACONS in 2D
// Uses the gkyl_kann_net wrapper (identical results to rt_kann_cmp_arch)
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <stdlib.h>
#include <unistd.h>
#include <math.h>

struct xrange {
  float xleft, xright;
  int N;
};

static inline float
xrange_n(struct xrange xr, int n)
{
  float dx = (xr.xright-xr.xleft)/(xr.N-1);
  return xr.xleft + dx*n;
}

// function to fit
static inline float
ufunc(float t, float x)
{
  return sinf(2.0f*M_PI*x)*sinf(2.0f*t)*expf(-t);
}

enum ann_layer_type {
  ANN_DENSE,
  ANN_GRU,
  ANN_GRU_NORM
};

struct train_inp {
  int ntrain[2];
  int nwidth;
  int ndepth;
  float learning_rate;
  enum ann_layer_type layer_type;
  bool use_gpu;
};

void
train_ann(struct train_inp *nn_inp, const char *nn_name)
{
  kad_node_t *t_net;
  t_net = kann_layer_input(2);

  for (int i=0; i<nn_inp->ndepth; ++i) {
    switch(nn_inp->layer_type) {
      case ANN_DENSE:
        t_net = kann_layer_dense(t_net, nn_inp->nwidth);
        break;
      case ANN_GRU:
        t_net = kann_layer_gru(t_net, nn_inp->nwidth, 0);
        break;
      case ANN_GRU_NORM:
        t_net = kann_layer_gru(t_net, nn_inp->nwidth, KANN_RNN_NORM);
        break;
    }
    t_net = kad_tanh(t_net);
  }

  t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
  struct gkyl_kann_net *net = gkyl_kann_net_new(t_net, nn_inp->use_gpu);

  // allocate memory for input/output vectors
  int Nt = nn_inp->ntrain[0], Nx = nn_inp->ntrain[1];
  int N = Nt*Nx;
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 2);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);

  struct xrange tr = {
    .xleft = 0.0f,
    .xright = 3.0f,
    .N = Nt
  };

  struct xrange xr = {
    .xleft = 0.0f,
    .xright = 1.0f,
    .N = Nx
  };

  // initialize input/output mapping
  for (int i=0; i<Nt; ++i)
    for (int j=0; j<Nx; ++j) {
      long idx = i*Nx+j;

      double t = inp->vals[idx][0] = xrange_n(tr, i);
      double x = inp->vals[idx][1] = xrange_n(xr, j);

      out->vals[idx][0] = ufunc(t, x);
    }

  struct gkyl_kann_train_params params = {
    .learning_rate = nn_inp->learning_rate,
    .mini_size = 64,
    .max_epoch = 50,
    .max_drop_streak = 10,
    .frac_val = 0.1f,
  };

  if (nn_inp->use_gpu) {
    struct gkyl_kn_vec *inp_cu = gkyl_kn_vec_cu_dev_new(N, 2);
    struct gkyl_kn_vec *out_cu = gkyl_kn_vec_cu_dev_new(N, 1);
    gkyl_kn_vec_copy(inp_cu, inp);
    gkyl_kn_vec_copy(out_cu, out);

    gkyl_kann_net_train_fnn1(net, &params, inp_cu, out_cu);

    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  } else {
    gkyl_kann_net_train_fnn1(net, &params, inp, out);
  }

  gkyl_kann_net_save(net, nn_name);

  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  gkyl_kann_net_release(net);
}

// run inference on N input values (batch mode)
void
infer_ann(const char *nn_name, bool use_gpu,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(nn_name, use_gpu);

  if (use_gpu) {
    struct gkyl_kn_vec *inp_cu = gkyl_kn_vec_cu_dev_new(inp->nvec, inp->N);
    struct gkyl_kn_vec *out_cu = gkyl_kn_vec_cu_dev_new(out->nvec, out->N);

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

// run sequential RNN inference (one timestep at a time with recurrence)
void
infer_ann_rnn(const char *nn_name, bool use_gpu,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  struct gkyl_kann_net *net = gkyl_kann_net_load(nn_name, use_gpu);

  if (use_gpu) {
    struct gkyl_kn_vec *inp_cu = gkyl_kn_vec_cu_dev_new(inp->nvec, inp->N);
    struct gkyl_kn_vec *out_cu = gkyl_kn_vec_cu_dev_new(out->nvec, out->N);

    gkyl_kn_vec_copy(inp_cu, inp);
    gkyl_kann_net_apply_rnn(net, inp_cu, out_cu);
    gkyl_kn_vec_copy(out, out_cu);

    gkyl_kn_vec_release(inp_cu);
    gkyl_kn_vec_release(out_cu);
  } else {
    gkyl_kann_net_apply_rnn(net, inp, out);
  }

  gkyl_kann_net_release(net);
}

void
write_to_gplot(void)
{
  const char *gp_code =
    "set macros\n"
    "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue\n"
    "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red\n"
    "BLUE = \"1\"\n"
    "RED = \"2\"\n"
    "set grid\n"
    "set title \"Comparing Different Architectures\"\n"
    "plot [0:3] sin(0.35*2*pi)*sin(2*x)*exp(-x) with lines ls @BLUE title \"Exact\" ";

  FILE *fp = 0;
  with_file(fp, "rt_kann_cmp_arch_gkw.gp", "w") {
    fprintf(fp, "%s", gp_code);
    fprintf(fp, ", \"rt_kann_cmp_arch_gkw_mlp.txt\" using 1:2 with points pt 9 ps 3 title \"MLP\" ");
    fprintf(fp, ", \"rt_kann_cmp_arch_gkw_gru.txt\" using 1:2 with points pt 5 ps 2 title \"GRU\" ");
    fprintf(fp, "\n");
  }
}

void
write_infer_data(const char *fname, const struct gkyl_kn_vec *inp, const struct gkyl_kn_vec *out)
{
  FILE *fp = 0;
  with_file(fp, fname, "w") {
    for (int i=0; i<inp->nvec; ++i)
      fprintf(fp, "%.5g %.5g\n", inp->vals[i][0], out->vals[i][0]);
  }
}

int
main(int argc, char *argv[])
{
  int p_train = 0, p_infer = 0, p_verbose = 0, c;
  bool use_gpu = false;
  while ((c = getopt(argc, argv, "+htivg")) != -1) {
    switch (c)
    {
      case 'h':
        fprintf(stdout, "rt_kann_cmp_arch_gkw -i -t -v -g\n");
        fprintf(stdout, "  -t Run Training\n");
        fprintf(stdout, "  -i Run Inference\n");
        fprintf(stdout, "  -v Verbose mode\n");
        fprintf(stdout, "  -g Use GPU\n");
        exit(0);
        break;

      case 't':
        p_train = 1;
        break;

      case 'i':
        p_infer = 1;
        break;

      case 'v':
        p_verbose = 3;
        break;

      case 'g':
        use_gpu = true;
        break;

      case '?':
        break;
    }
  }

  gkyl_kann_net_set_verbose(p_verbose);

  if (p_train) {
    fprintf(stdout, "*** Training MLP%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    train_ann( &(struct train_inp) {
        .ntrain = { 101, 101 },
        .ndepth = 2,
        .nwidth = 64,
        .learning_rate = 1e-3f,
        .layer_type = ANN_DENSE,
        .use_gpu = use_gpu
      },
      "rt_kann_cmp_arch_gkw_mlp.kann"
    );

    fprintf(stdout, "*** Training GRU%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    train_ann( &(struct train_inp) {
        .ntrain = { 101, 101 },
        .ndepth = 2,
        .nwidth = 32,
        .learning_rate = 1e-3f,
        .layer_type = ANN_GRU,
        .use_gpu = use_gpu
      },
      "rt_kann_cmp_arch_gkw_gru.kann"
    );

    fprintf(stdout, "*** Training GRU+Norm%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    train_ann( &(struct train_inp) {
        .ntrain = { 101, 101 },
        .ndepth = 2,
        .nwidth = 32,
        .learning_rate = 1e-3f,
        .layer_type = ANN_GRU_NORM,
        .use_gpu = use_gpu
      },
      "rt_kann_cmp_arch_gkw_gru_norm.kann"
    );
  }

  if (p_infer) {
    int nvec = 11;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 2);
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

    struct xrange tr = { .xleft = 0.0f, .xright = 3.0f, .N = inp->nvec };
    for (int i=0; i<inp->nvec; ++i) {
      inp->vals[i][0] = xrange_n(tr, i);
      inp->vals[i][1] = 0.35f;
    }

    fprintf(stdout, "*** MLP Inference%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    infer_ann("rt_kann_cmp_arch_gkw_mlp.kann", use_gpu, inp, out);
    write_infer_data("rt_kann_cmp_arch_gkw_mlp.txt", inp, out);

    fprintf(stdout, "*** GRU Inference%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    infer_ann("rt_kann_cmp_arch_gkw_gru.kann", use_gpu, inp, out);
    write_infer_data("rt_kann_cmp_arch_gkw_gru.txt", inp, out);

    fprintf(stdout, "*** GRU+Norm Inference%s (gkyl_kann_net wrapper)\n",
      use_gpu ? " (GPU)" : "");
    infer_ann("rt_kann_cmp_arch_gkw_gru_norm.kann", use_gpu, inp, out);
    write_infer_data("rt_kann_cmp_arch_gkw_gru_norm.txt", inp, out);

    // Sequential RNN inference: exercises pre-recurrence (hidden state
    // carries forward across timesteps)
    fprintf(stdout, "*** GRU Sequential RNN Inference%s\n",
      use_gpu ? " (GPU)" : "");
    infer_ann_rnn("rt_kann_cmp_arch_gkw_gru.kann", use_gpu, inp, out);
    write_infer_data("rt_kann_cmp_arch_gkw_gru_rnn.txt", inp, out);

    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);
  }

  write_to_gplot();

  return 0;
}
