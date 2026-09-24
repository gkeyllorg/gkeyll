#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <kann.h>
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
  kann_t *ann = kann_new(t_net, 0);

  // Hyperparameters for training.
  float lr = nn_inp->learning_rate; // Learning rate.
  int mini_size = 64;
  int max_epoch = 50;
  int max_drop_streak = 10;
  float frac_val = 0.1f; // Fraction of samples to use for validation.

  // Run individual expert training (i.e. pretraining).
  for (int i = 0; i < nn_inp->nexperts; i++) {
    // Allocate memory for input/output vectors for individual experts.
    int N_expert = nn_inp->ntrain; // Training samples
    struct gkyl_kn_vec *inp_expert = gkyl_kn_vec_new(N_expert, 1);
    struct gkyl_kn_vec *out_expert = gkyl_kn_vec_new(N_expert, 1);

    struct xrange xr_expert = {
      .xleft = -1.0,
      .xright = 1.0,
      .N = N_expert
    };

    // Initialize input/output mapping for individual experts.
    for (int j = 0; j < N_expert; j++) {
      inp_expert->vals[j][0] = xrange_n(xr_expert, j);
      out_expert->vals[j][0] = tfunc(i, inp_expert->vals[j][0]);
    }

    kann_train_fnn1(ann, lr, mini_size, max_epoch, max_drop_streak, frac_val, N_expert, inp_expert->vals, out_expert->vals);

    gkyl_kn_vec_release(inp_expert);
    gkyl_kn_vec_release(out_expert);
  }

  // Allocate memory for input/output vectors.
  int N = nn_inp->ntrain; // Training samples
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);

  struct xrange xr = {
    .xleft = -1.0,
    .xright = 1.0,
    .N = N
  };

  // Initialize input/output mapping.
  for (int i = 0; i < N; i++) {
    inp->vals[i][0] = xrange_n(xr, i);
    out->vals[i][0] = ufunc(inp->vals[i][0]);
  }
  
  // Run mixture of experts training (i.e. finetuning).
  kann_train_fnn1(ann, lr, mini_size, max_epoch, max_drop_streak, frac_val, N, inp->vals, out->vals);
  kann_save(nn_name, ann); // Save to file.
  
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  kann_delete(ann);  
}

// Run inference on N input values.
void
infer_ann(const char *nn_name, const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  kann_t *ann = kann_load(nn_name);
  const float *ov;
  for (int i = 0; i < inp->nvec; i++) {
    ov = kann_apply1(ann, inp->vals[i]);
    for (int j = 0; j < out->N; j++) {
      out->vals[i][j] = ov[j];
    }
  }
  kann_delete(ann);
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
    "plot \"rt_kann_moe_data.txt\" using 1:2 with points pt 9 ps 3 title \"NN\", [-1:1] 1/(1+100*x**2) with lines ls @BLUE title \"Exact\"";

  FILE *fp = 0;
  with_file(fp, "rt_kann_moe.gp", "w") {
    fprintf(fp, "%s", gpcode);
  }

  fp = 0;
  with_file(fp, "rt_kann_moe_data.txt", "w") {
    for (int i = 0; i < inp->nvec; i++) {
      fprintf(fp, "%.5g %.5g\n", inp->vals[i][0], out->vals[i][0]);
    }
  }
}

int
main(int argc, char *argv[])
{
  int p_train = 0, p_infer = 0, p_verbose = 0, c;
  while ((c = getopt(argc, argv, "+htiv")) != -1) {
    switch (c)
    {
      case 'h':
        fprintf(stdout, "rt_kann_moe -i -t -v\n");
        fprintf(stdout, "  -t Run Training\n");
        fprintf(stdout, "  -i Run Inference\n");
        fprintf(stdout, "  -v Verbose mode\n");
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
      
      case '?':
        break;
    }
  }

  kann_set_verbose_level(p_verbose);

  if (p_train) {
    fprintf(stdout, "*** Training\n");
    train_mixture( &(struct train_inp) {
        .ntrain = 1001,
        .ndepth = 2,
        .nwidth = 256,
        .nexperts = 3,
        .learning_rate = 1e-3f
      },
      "rt_kann_moe.kann"
    );
  }

  if (p_infer) {
    fprintf(stdout, "*** Inference\n");
    // Run inference.
    int nvec = 11;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 1);
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

    struct xrange xr = { .xleft = -1.0, .xright = 1.0, .N = inp->nvec };
    for (int i = 0; i < inp->nvec; i++) {
      inp->vals[i][0] = xrange_n(xr, i);
    }
  
    infer_ann("rt_kann_moe.kann", inp, out);
    write_to_gplot(inp, out);
    
    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);
  }
  
  return 0;
}