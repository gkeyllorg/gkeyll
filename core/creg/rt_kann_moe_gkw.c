#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <stdlib.h>
#include <unistd.h>
#include <gkyl_alloc.h>
#include <math.h>

struct xrange {
  float xleft, xright;
  int N;
};

static inline float
xrange_n(struct xrange xr, int n)
{
  float dx = (xr.xright - xr.xleft) / (xr.N - 1);
  return xr.xleft + dx * n;
}

// Family of pretraining functions.
static inline float
tfunc(int n, float x)
{
  return sinf((2.0f * M_PI * (n + 1) * x) + n);
}

// Function to fit.
static inline float
ufunc(float x)
{
  return 1.0f / (1.0f + 100.0f * x * x);
}

struct train_inp {
  int ntrain;
  int nwidth;
  int ndepth;
  int nexperts;
  float learning_rate;
  bool use_gpu;
};

// Construct an MLP (with tanh activation) to use as a single "expert".
static inline kad_node_t*
single_expert(kad_node_t *input, int n_layers, int n_hidden, int n_output)
{
  kad_node_t *t_net;

  for (int i = 0; i < n_layers; i++) {
    t_net = kann_layer_dense(input, n_hidden);
    t_net = kad_tanh(t_net);
  }
  t_net = kann_layer_dense(t_net, n_output);

  return t_net;
}

// Construct a "weighted expert", with an initial (trainable) scalar weight and a trainable scalar bias.
static inline kad_node_t*
weighted_expert(kad_node_t *expert, int n_output, float init_weight)
{
  kad_node_t *weight, *bias;
  kad_node_t *t_net;

  weight = kann_new_scalar(KAD_VAR, init_weight);
  bias = kann_new_bias(n_output);

  t_net = kad_mul(expert, weight);
  t_net = kad_add(t_net, bias);

  return t_net;
}

// Construct a "mixture of experts" architecture consisting of multiple single "experts" linked together (with tanh activation and MSE cost).
static inline kad_node_t*
mixture_of_experts(int n_input, int n_layers, int n_hidden, int n_experts, int n_output)
{
  kad_node_t *input;
  kad_node_t **experts;
  kad_node_t **weighted_experts;

  kad_node_t *opinion_sum;
  kad_node_t *activation, *truth, *cost;

  input = kann_layer_input(n_input);

  experts = gkyl_malloc(n_experts * sizeof(kad_node_t*));
  for (int i = 0; i < n_experts; i++) {
    experts[i] = single_expert(input, n_layers, n_hidden, n_output);
  }

  weighted_experts = gkyl_malloc(n_experts * sizeof(kad_node_t*));
  for (int i = 0; i < n_experts; i++) {
    weighted_experts[i] = weighted_expert(experts[i], n_output, 1.0f / n_experts);
  }

  opinion_sum = weighted_experts[0];
  for (int i = 1; i < n_experts; i++) {
    opinion_sum = kad_add(opinion_sum, weighted_experts[i]);
  }

  activation = kad_tanh(opinion_sum);
  activation->ext_flag |= KANN_F_OUT;

  truth = kad_feed(2, 1, n_output);
  truth->ext_flag |= KANN_F_TRUTH;

  cost = kad_mse(activation, truth);
  cost->ext_flag |= KANN_F_COST;
  
  gkyl_free(experts);
  gkyl_free(weighted_experts);

  return cost;
}

void
train_mixture(struct train_inp *nn_inp, const char *nn_name)
{
  kad_node_t *t_net = mixture_of_experts(1, nn_inp->ndepth, nn_inp->nwidth, nn_inp->nexperts, 1);
  struct gkyl_kann_net *net = gkyl_kann_net_new(t_net, nn_inp->use_gpu);

  // Hyperparameters for training.
  struct gkyl_kann_train_params params = {
    .learning_rate = nn_inp->learning_rate,
    .mini_size = 64,
    .max_epoch = 50,
    .max_drop_streak = 10,
    .frac_val = 0.1f,
  };

  // Run individual expert training (i.e. pretraining).
  for (int i = 0; i < nn_inp->nexperts; i++) {
    int N_expert = nn_inp->ntrain;

    // Host kn_vecs for filling data for individual experts.
    struct gkyl_kn_vec *inp_expert = gkyl_kn_vec_new(N_expert, 1);
    struct gkyl_kn_vec *out_expert = gkyl_kn_vec_new(N_expert, 1);

    struct xrange xr_expert = {
      .xleft = -1.0,
      .xright = 1.0,
      .N = N_expert
    };

    for (int j = 0; j < N_expert; j++) {
      inp_expert->vals[j][0] = xrange_n(xr_expert, j);
      out_expert->vals[j][0] = tfunc(i, inp_expert->vals[j][0]);
    }

    if (nn_inp->use_gpu) {
      // Create device kn_vecs and copy data H2D for individual experts.
      struct gkyl_kn_vec *inp_cu_expert = gkyl_kn_vec_cu_dev_new(N_expert, 1);
      struct gkyl_kn_vec *out_cu_expert = gkyl_kn_vec_cu_dev_new(N_expert, 1);
      gkyl_kn_vec_copy(inp_cu_expert, inp_expert);
      gkyl_kn_vec_copy(out_cu_expert, out_expert);

      gkyl_kann_net_train_fnn1(net, &params, inp_cu_expert, out_cu_expert);

      gkyl_kn_vec_release(inp_cu_expert);
      gkyl_kn_vec_release(out_cu_expert);
    } else {
      gkyl_kann_net_train_fnn1(net, &params, inp_expert, out_expert);
    }

    gkyl_kn_vec_release(inp_expert);
    gkyl_kn_vec_release(out_expert);
  }

  int N = nn_inp->ntrain;

  // Host kn_vecs for filling data.
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);

  struct xrange xr = {
    .xleft = -1.0,
    .xright = 1.0,
    .N = N
  };

  for (int i = 0; i < N; i++) {
    inp->vals[i][0] = xrange_n(xr, i);
    out->vals[i][0] = ufunc(inp->vals[i][0]);
  }

  if (nn_inp->use_gpu) {
    // Create device kn_vecs and copy data H2D.
    struct gkyl_kn_vec *inp_cu = gkyl_kn_vec_cu_dev_new(N, 1);
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

// Run inference on N input values.
void
infer_ann(const char *nn_name, bool use_gpu,
  struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
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

void
write_to_gplot(const struct gkyl_kn_vec *inp, const struct gkyl_kn_vec *out)
{
  const char *gpcode =
    "set macros\n"
    "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue\n"
    "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red\n"
    "BLUE = \"1\"\n"
    "RED = \"2\"\n"
    "set grid\n"
    "plot \"rt_kann_moe_gkw_data.txt\" using 1:2 with points pt 9 ps 3 title \"NN\", [-1:1] 1/(1+100*x**2) with lines ls @BLUE title \"Exact\"";

  FILE *fp = 0;
  with_file(fp, "rt_kann_moe_gkw.gp", "w") {
    fprintf(fp, "%s", gpcode);
  }

  fp = 0;
  with_file(fp, "rt_kann_moe_gkw_data.txt", "w") {
    for (int i = 0; i < inp->nvec; i++) {
      fprintf(fp, "%.5g %.5g\n", inp->vals[i][0], out->vals[i][0]);
    }
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
        fprintf(stdout, "rt_kann_moe_gkw -i -t -v -g\n");
        fprintf(stdout, "  -t Run Training\n");
        fprintf(stdout, "  -i Run Inference\n");
        fprintf(stdout, "  -v Verbose mode\n");
        fprintf(stdout, "  -g Run on GPU\n");
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
    fprintf(stdout, "*** Training%s\n", use_gpu ? " (GPU)" : "");
    train_mixture( &(struct train_inp) {
        .ntrain = 1001,
        .ndepth = 2,
        .nwidth = 256,
        .nexperts = 3,
        .learning_rate = 1e-3f,
        .use_gpu = use_gpu
      },
      "rt_kann_moe_gkw.kann"
    );
  }

  if (p_infer) {
    fprintf(stdout, "*** Inference%s\n", use_gpu ? " (GPU)" : "");
    int nvec = 11;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 1);
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

    struct xrange xr = { .xleft = -1.0, .xright = 1.0, .N = inp->nvec };
    for (int i = 0; i < inp->nvec; i++) {
      inp->vals[i][0] = xrange_n(xr, i);
    }

    infer_ann("rt_kann_moe_gkw.kann", use_gpu, inp, out);
    write_to_gplot(inp, out);

    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);
  }

  return 0;
}
