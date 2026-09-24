// Compare MLP, GRU and BEACONS in 2D
#include <gkyl_knutils.h>
#include <gkyl_util.h>

#include <kann.h>
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
  ANN_GRU
};

struct train_inp {
  int ntrain[2];
  int nwidth;
  int ndepth;
  float learning_rate;
  enum ann_layer_type layer_type;
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
    }
    t_net = kad_tanh(t_net);
  }
  
  t_net = kann_layer_cost(t_net, 1, KANN_C_MSE);
  kann_t *ann = kann_new(t_net, 0);

  // allocate memory for input/output vectors
  int Nt = nn_inp->ntrain[0], Nx = nn_inp->ntrain[1]; // training samples
  struct gkyl_kn_vec *inp = gkyl_kn_vec_new(Nt*Nx, 2);
  struct gkyl_kn_vec *out = gkyl_kn_vec_new(Nt*Nx, 1);

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

  // hyper-parameters for training
  float lr = nn_inp->learning_rate; // learning rate
  int mini_size = 64;
  int max_epoch = 50;
  int max_drop_streak = 10;
  float frac_val = 0.1f; // fraction of samples to use for validation
  
  // run training
  kann_train_fnn1(ann, lr, mini_size, max_epoch, max_drop_streak, frac_val, Nx*Nt, inp->vals, out->vals);
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
  with_file(fp, "rt_kann_cmp_arch.gp", "w") {
    fprintf(fp, "%s", gp_code);
    fprintf(fp, ", \"rt_kann_cmp_arch_mlp.txt\" using 1:2 with points pt 9 ps 3 title \"MLP\" ");
    fprintf(fp, ", \"rt_kann_cmp_arch_gru.txt\" using 1:2 with points pt 5 ps 2 title \"GRU\" ");    
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
  while ((c = getopt(argc, argv, "+htiv")) != -1) {
    switch (c)
    {
      case 'h':
        fprintf(stdout, "rt_kann_cmp_arch -i -t -v\n");
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
    fprintf(stdout, "*** Training MLP\n");
    train_ann( &(struct train_inp) {
        .ntrain = { 101, 101 },
        .ndepth = 2,
        .nwidth = 64,
        .learning_rate = 1e-3f,
        .layer_type = ANN_DENSE
      },
      "rt_kann_cmp_arch_mlp.kann"
    );

    fprintf(stdout, "*** Training GRU\n");
    train_ann( &(struct train_inp) {
        .ntrain = { 101, 101 },
        .ndepth = 2,
        .nwidth = 32,
        .learning_rate = 1e-3f,
        .layer_type = ANN_GRU
      },
      "rt_kann_cmp_arch_gru.kann"
    );    
  }

  if (p_infer) {
    // run inference
    int nvec = 11;
    struct gkyl_kn_vec *inp = gkyl_kn_vec_new(nvec, 2);
    struct gkyl_kn_vec *out = gkyl_kn_vec_new(nvec, 1);

    struct xrange tr = { .xleft = 0.0f, .xright = 3.0f, .N = inp->nvec };
    for (int i=0; i<inp->nvec; ++i) {
      inp->vals[i][0] = xrange_n(tr, i);
      inp->vals[i][1] = 0.35f;
    }

    fprintf(stdout, "*** MLP Inference\n");
    infer_ann("rt_kann_cmp_arch_mlp.kann", inp, out);
    write_infer_data("rt_kann_cmp_arch_mlp.txt", inp, out);

    fprintf(stdout, "*** GRU Inference\n");
    infer_ann("rt_kann_cmp_arch_gru.kann", inp, out);
    write_infer_data("rt_kann_cmp_arch_gru.txt", inp, out);

    gkyl_kn_vec_release(inp);
    gkyl_kn_vec_release(out);
  }
    
  write_to_gplot();
  
  return 0;
}
