#include <gkyl_alloc.h>
#include <gkyl_kann_net.h>
#include <gkyl_kann_net_priv.h>

#ifdef GKYL_HAVE_CUDA
#include <gkyl_kann_net_cu_priv.h>
#endif

#include <assert.h>
#include <float.h>
#include <string.h>

// Mirror of kann_verbose in kann.c (set via gkyl_kann_net_set_verbose)
static int kann_net_verbose = 0;

static void
kann_net_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_kann_net *net = container_of(ref, struct gkyl_kann_net, ref_count);
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(net->flags)) {
    if (net->cg)
      kann_cu_graph_free(net->cg);
    cublasDestroy(net->cublas_h);
  }
#endif
  if (net->ann)
    kann_delete(net->ann);
  gkyl_free(net);
}

// Initialize cached dimensions from the underlying kann_t
static void
kann_net_cache_dims(struct gkyl_kann_net *net)
{
  net->n_in = kann_dim_in(net->ann);
  net->n_out = kann_dim_out(net->ann);
  net->n_var = kann_size_var(net->ann);
  net->n_const = kann_size_const(net->ann);
}

struct gkyl_kann_net*
gkyl_kann_net_new(kad_node_t *cost, bool use_gpu)
{
  struct gkyl_kann_net *net = gkyl_malloc(sizeof(*net));

  net->ann = kann_new(cost, 0);
  if (!net->ann) {
    gkyl_free(net);
    return 0;
  }

  net->flags = 0;
  kann_net_cache_dims(net);
  net->ref_count = gkyl_ref_count_init(kann_net_free);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    GKYL_SET_CU_ALLOC(net->flags);
    cublasCreate_v2(&net->cublas_h);
    net->cg = 0; // lazily created when batch size is known
  }
#else
  assert(use_gpu == false);
#endif

  return net;
}

struct gkyl_kann_net*
gkyl_kann_net_load(const char *filename, bool use_gpu)
{
  kann_t *ann = kann_load(filename);
  if (!ann) return 0;

  struct gkyl_kann_net *net = gkyl_malloc(sizeof(*net));
  net->ann = ann;
  net->flags = 0;
  kann_net_cache_dims(net);
  net->ref_count = gkyl_ref_count_init(kann_net_free);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    GKYL_SET_CU_ALLOC(net->flags);
    cublasCreate_v2(&net->cublas_h);
    net->cg = 0; // lazily created when batch size is known
  }
#else
  assert(use_gpu == false);
#endif

  return net;
}

void
gkyl_kann_net_save(const struct gkyl_kann_net *net, const char *filename)
{
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(net->flags) && net->cg) {
    // Download weights from GPU before saving
    kann_cu_graph_download_vars(net->cg, net->ann);
  }
#endif
  kann_save(filename, net->ann);
}

// GPU training loop: mirrors kann_train_fnn1 but runs forward/backward on device.
// All training data lives on device via kn_vec on_dev; mini-batches are
// assembled on device via gather kernels.
#ifdef GKYL_HAVE_CUDA
static int
kann_net_train_fnn1_cu(struct gkyl_kann_net *net,
  const struct gkyl_kann_train_params *params,
  const struct gkyl_kn_vec *inp, const struct gkyl_kn_vec *out)
{
  kann_t *ann = net->ann;
  int n = inp->nvec;
  int n_in = net->n_in, n_out = net->n_out;
  int mini_size = params->mini_size;
  int n_var = net->n_var, n_const = net->n_const;

  // Lazily create GPU graph sized for this mini-batch
  if (!net->cg || net->cg->max_batch_size < mini_size) {
    if (net->cg) kann_cu_graph_free(net->cg);
    net->cg = kann_cu_graph_new(ann, mini_size);
  }
  struct kann_cu_graph *cg = net->cg;

  // inp and out must be device kn_vecs
  assert(gkyl_kn_vec_is_cu_dev(inp) && gkyl_kn_vec_is_cu_dev(out));
  float *d_inp = inp->data;
  float *d_out = out->data;

  // Shuffle indices on host, upload to device once per epoch
  int *shuf = (int *)malloc(n * sizeof(int));
  kann_shuffle(n, shuf);

  // Build initial train/val split: shuf[0..n_train-1] = train,
  // shuf[n_train..n-1] = val (indices into original data)
  int *idx = (int *)malloc(n * sizeof(int));
  for (int j = 0; j < n; ++j)
    idx[j] = shuf[j];

  int n_val = (int)(n * params->frac_val);
  int n_train = n - n_val;

  // Device buffer for shuffle indices
  int *d_idx = (int *)gkyl_cu_malloc(n * sizeof(int));

  float *min_x = (float *)malloc(n_var * sizeof(float));
  float *min_c = (float *)malloc(n_const * sizeof(float));

  // Upload initial variables and constants to GPU
  kann_cu_graph_upload_vars(cg, ann);
  kann_cu_graph_upload_consts(cg, ann);

  float min_val_cost = FLT_MAX;
  int drop_streak = 0, min_set = 0;
  int epoch;

  for (epoch = 0; epoch < params->max_epoch; ++epoch) {
    int n_proc = 0;
    double train_cost = 0.0;

    // Reshuffle training indices and upload to device
    kann_shuffle(n_train, shuf);
    for (int j = 0; j < n_train; ++j)
      shuf[j] = idx[shuf[j]];
    gkyl_cu_memcpy(d_idx, shuf, n_train * sizeof(int), GKYL_CU_MEMCPY_H2D);

    // Training pass
    while (n_proc < n_train) {
      int ms = n_train - n_proc < mini_size ? n_train - n_proc : mini_size;

      kann_cu_sync_dim(cg, ann, ms);

      // Gather mini-batch on device
      kann_cu_gather_input(cg, d_inp, d_idx, n_proc, ms, n_in);
      kann_cu_gather_truth(cg, d_out, d_idx, n_proc, ms, n_out);

      kann_cu_forward(cg, net->cublas_h, cg->cost_node_idx);

      float cost = kann_cu_get_cost(cg);
      train_cost += cost * ms;

      kann_cu_backward(cg, net->cublas_h, cg->cost_node_idx);
      kann_cu_rmsprop(cg, params->learning_rate, 0.9f);

      n_proc += ms;
    }
    train_cost /= n_train;

    // Validation pass (uses fixed val indices from initial split)
    double val_cost = 0.0;
    if (n_val > 0) {
      // Upload validation indices to device
      gkyl_cu_memcpy(d_idx, &idx[n_train], n_val * sizeof(int),
        GKYL_CU_MEMCPY_H2D);

      n_proc = 0;
      while (n_proc < n_val) {
        int ms = n_val - n_proc < mini_size ? n_val - n_proc : mini_size;

        kann_cu_sync_dim(cg, ann, ms);
        kann_cu_gather_input(cg, d_inp, d_idx, n_proc, ms, n_in);
        kann_cu_gather_truth(cg, d_out, d_idx, n_proc, ms, n_out);
        kann_cu_forward(cg, net->cublas_h, cg->cost_node_idx);

        float cost = kann_cu_get_cost(cg);
        val_cost += cost * ms;
        n_proc += ms;
      }
      val_cost /= n_val;
    }

    if (kann_net_verbose >= 3)
      fprintf(stderr, "epoch: %d; training cost: %g; validation cost: %g\n",
        epoch + 1, train_cost, val_cost);

    if (epoch >= params->max_drop_streak && n_val > 0) {
      if (val_cost < min_val_cost) {
        min_set = 1;
        kann_cu_graph_download_vars(cg, ann);
        memcpy(min_x, ann->x, n_var * sizeof(float));
        memcpy(min_c, ann->c, n_const * sizeof(float));
        drop_streak = 0;
        min_val_cost = (float)val_cost;
      } else if (++drop_streak >= params->max_drop_streak)
        break;
    }
  }

  // Restore best weights
  if (min_set) {
    memcpy(ann->x, min_x, n_var * sizeof(float));
    memcpy(ann->c, min_c, n_const * sizeof(float));
    kann_cu_graph_upload_vars(cg, ann);
    kann_cu_graph_upload_consts(cg, ann);
  } else {
    kann_cu_graph_download_vars(cg, ann);
  }

  gkyl_cu_free(d_idx);
  free(min_c); free(min_x);
  free(idx); free(shuf);
  return epoch;
}
#endif

int
gkyl_kann_net_train_fnn1(struct gkyl_kann_net *net,
  const struct gkyl_kann_train_params *params,
  const struct gkyl_kn_vec *inp, const struct gkyl_kn_vec *out)
{
  assert(inp->nvec == out->nvec);
  assert(inp->N == net->n_in);
  assert(out->N == net->n_out);

#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(net->flags))
    return kann_net_train_fnn1_cu(net, params, inp, out);
#endif

  return kann_train_fnn1(net->ann,
    params->learning_rate,
    params->mini_size,
    params->max_epoch,
    params->max_drop_streak,
    params->frac_val,
    inp->nvec,
    inp->vals,
    out->vals);
}

// GPU inference: feed device input, run forward, copy output to device kn_vec
#ifdef GKYL_HAVE_CUDA
static void
kann_net_apply_cu(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  int n_in = net->n_in, n_out = net->n_out;
  int nvec = inp->nvec;

  // inp and out must be device kn_vecs
  assert(gkyl_kn_vec_is_cu_dev(inp) && gkyl_kn_vec_is_cu_dev(out));

  // Lazily create GPU graph if needed
  if (!net->cg || net->cg->max_batch_size < nvec) {
    if (net->cg) kann_cu_graph_free(net->cg);
    net->cg = kann_cu_graph_new(net->ann, nvec);
  }
  struct kann_cu_graph *cg = net->cg;

  // Sync dimensions and feed directly from device kn_vec data
  kann_cu_sync_dim(cg, net->ann, nvec);
  kann_cu_feed_input_dev(cg, nvec, inp->data, n_in);

  int eval_to = cg->out_node_idx >= 0 ? cg->out_node_idx : cg->cost_node_idx;
  kann_cu_forward(cg, net->cublas_h, eval_to);

  // Copy output from graph directly to device kn_vec
  int idx = cg->out_node_idx >= 0 ? cg->out_node_idx : cg->cost_node_idx;
  struct kann_cu_node *hn = &cg->h_nodes[idx];
  gkyl_cu_memcpy(out->data, cg->x + hn->x_off,
    nvec * n_out * sizeof(float), GKYL_CU_MEMCPY_D2D);
}
#endif

void
gkyl_kann_net_apply(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  assert(inp->N == net->n_in);
  assert(out->N == net->n_out);
  assert(inp->nvec == out->nvec);

#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(net->flags)) {
    kann_net_apply_cu(net, inp, out);
    return;
  }
#endif

  for (int i = 0; i < inp->nvec; ++i) {
    const float *ov = kann_apply1(net->ann, inp->vals[i]);
    memcpy(out->vals[i], ov, net->n_out * sizeof(float));
  }
}

// GPU sequential RNN inference: process one timestep at a time with
// pre-recurrence between steps.
#ifdef GKYL_HAVE_CUDA
static void
kann_net_apply_rnn_cu(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  int n_in = net->n_in, n_out = net->n_out;
  int nvec = inp->nvec;

  assert(gkyl_kn_vec_is_cu_dev(inp) && gkyl_kn_vec_is_cu_dev(out));

  // Graph needs batch_size=1 for sequential processing
  if (!net->cg || net->cg->max_batch_size < 1) {
    if (net->cg) kann_cu_graph_free(net->cg);
    net->cg = kann_cu_graph_new(net->ann, 1);
  }
  struct kann_cu_graph *cg = net->cg;

  kann_cu_sync_dim(cg, net->ann, 1);

  int eval_to = cg->out_node_idx >= 0 ? cg->out_node_idx : cg->cost_node_idx;
  int out_idx = cg->out_node_idx >= 0 ? cg->out_node_idx : cg->cost_node_idx;

  // Zero h0 nodes (initial hidden state)
  for (int i = 0; i < cg->n_pre_pairs; ++i) {
    int h0_idx = cg->h_pre_pairs[2*i+1];
    struct kann_cu_node *h0_n = &cg->h_nodes[h0_idx];
    cudaMemset(cg->x + h0_n->x_off, 0, h0_n->len * sizeof(float));
  }

  for (int t = 0; t < nvec; ++t) {
    // Feed one timestep from device input
    kann_cu_feed_input_dev(cg, 1, inp->data + t * n_in, n_in);

    kann_cu_forward(cg, net->cublas_h, eval_to);

    // Copy output for this timestep
    struct kann_cu_node *hn = &cg->h_nodes[out_idx];
    gkyl_cu_memcpy(out->data + t * n_out, cg->x + hn->x_off,
      n_out * sizeof(float), GKYL_CU_MEMCPY_D2D);

    // Apply pre-recurrence: copy output node x to h0 node x
    kann_cu_apply_pre(cg);
  }
}
#endif

void
gkyl_kann_net_apply_rnn(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out)
{
  assert(inp->N == net->n_in);
  assert(out->N == net->n_out);
  assert(inp->nvec == out->nvec);

#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(net->flags)) {
    kann_net_apply_rnn_cu(net, inp, out);
    return;
  }
#endif

  // CPU path: use kann_rnn_start / kann_apply1 / kann_rnn_end
  kann_rnn_start(net->ann);
  for (int i = 0; i < inp->nvec; ++i) {
    const float *ov = kann_apply1(net->ann, inp->vals[i]);
    memcpy(out->vals[i], ov, net->n_out * sizeof(float));
  }
  kann_rnn_end(net->ann);
}

int
gkyl_kann_net_dim_in(const struct gkyl_kann_net *net)
{
  return net->n_in;
}

int
gkyl_kann_net_dim_out(const struct gkyl_kann_net *net)
{
  return net->n_out;
}

void
gkyl_kann_net_set_verbose(int level)
{
  kann_net_verbose = level;
  kann_set_verbose_level(level);
}

bool
gkyl_kann_net_is_cu_dev(const struct gkyl_kann_net *net)
{
  return GKYL_IS_CU_ALLOC(net->flags);
}

struct gkyl_kann_net*
gkyl_kann_net_acquire(const struct gkyl_kann_net *net)
{
  gkyl_ref_count_inc(&net->ref_count);
  return (struct gkyl_kann_net*) net;
}

void
gkyl_kann_net_release(struct gkyl_kann_net *net)
{
  if (net)
    gkyl_ref_count_dec(&net->ref_count);
}
