// Construct a MLP with KANN and fit a simple 1D function
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <kann.h>
#include <stdlib.h>
#include <unistd.h>

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
ufunc(float x)
{
  return 1.0f/(1.0f+100.0f*x*x);
}

struct train_inp {
  int ntrain;
  int nwidth;
  int ndepth;
  float learning_rate;
};

void
train_ann(struct train_inp *nn_inp, const char *nn_name)
{
  kad_node_t *t_net;
  t_net = kann_layer_input(1);

  for (int i=0; i<nn_inp->ndepth; ++i) {
    t_net = kann_layer_dense(t_net, nn_inp->nwidth);
    t_net = kad_tanh(t_net);
  }
  
  t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
  kann_t *ann = kann_new(t_net, 0);

  // allocate memory for input/output vectors
  int N = nn_inp->ntrain; // training samples
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(N, 1);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(N, 1);

  struct xrange xr = {
    .xleft = -1.0,
    .xright = 1.0,
    .N = N
  };

  // initialize input/output mapping
  for (int i=0; i<N; ++i) {
    inp->vals[i][0] = xrange_n(xr, i);
    out->vals[i][0] = ufunc(inp->vals[i][0]);
  }

  // hyper-parameters for training
  float lr = nn_inp->learning_rate; // learning rate
  int mini_size = 64;
  int max_epoch = 50;
  int max_drop_streak = 10;
  float frac_val = 0.1f; // fraction of samples to use for validation
  
  // run training
  kann_train_fnn1(ann, lr, mini_size, max_epoch, max_drop_streak, frac_val, N, inp->vals, out->vals);
  kann_save(nn_name, ann); // save to file
  
  gkyl_kn_vec_release(inp);
  gkyl_kn_vec_release(out);
  kann_delete(ann);  
}

// run inference on N input values
void
infer_ann(const char *nn_name, const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  kann_t *ann = kann_load(nn_name);
  const float *ov;
  for (int i=0; i<inp->nvec; ++i) {
    ov = kann_apply1(ann, inp->vals[i]);
    for (int j=0; j<out->N; ++j) out->vals[i][j] = ov[j];
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
    "plot \"rt_kann_mlp_data.txt\" using 1:2 with points pt 9 ps 3 title \"NN\", [-1:1] 1/(1+100*x**2) with lines ls @BLUE title \"Exact\"";

  FILE *fp = 0;
  with_file(fp, "rt_kann_mlp.gp", "w") {
    fprintf(fp, "%s", gpcode);
  }

  fp = 0;
  with_file(fp, "rt_kann_mlp_data.txt", "w") {
    for (int i=0; i<inp->nvec; ++i)
      fprintf(fp, "%.5g %.5g\n", inp->vals[i][0], out->vals[i][0]);
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
        fprintf(stdout, "rt_kann_mlp -i -t -v\n");
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
    train_ann( &(struct train_inp) {
        .ntrain = 1001,
        .ndepth = 2,
        .nwidth = 256,
        .learning_rate = 1e-3f
      },
      "rt_kann_mlp.kann"
    );
  }

  if (p_infer) {
    fprintf(stdout, "*** Inference\n");
    // run inference
    int nvec = 11;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 1);
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

    struct xrange xr = { .xleft = -1.0, .xright = 1.0, .N = inp->nvec };
    for (int i=0; i<inp->nvec; ++i)
      inp->vals[i][0] = xrange_n(xr, i);
  
    infer_ann("rt_kann_mlp.kann", inp, out);
    write_to_gplot(inp, out);
    
    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);
  }
  
  return 0;
}
