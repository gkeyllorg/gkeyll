/* -*- c++ -*- */

// CUDA kernels and device-side graph operations for the gkyl_kann_net GPU path.
// Implements the flattened graph construction, forward/backward passes,
// and training support (RMSprop, gradient zeroing).

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cmath>

#include <cublas_v2.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_kann_net.h>
#include <gkyl_kann_net_priv.h>
#include <gkyl_kann_net_cu_priv.h>
#include <gkyl_util.h>
}

#define KANN_CU_THREADS 256

// ============================================================
// Element-wise CUDA kernels for forward and backward passes
// ============================================================

// --- add (op 1): p = child[0] + child[1] (with broadcasting) ---
__global__ static void
ker_add_fwd(int n0, int n1, const float *x0, const float *x1, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n0;
       i += blockDim.x * gridDim.x)
    px[i] = x0[i] + x1[i % n1];
}

__global__ static void
ker_add_bwd_child0(int n, const float *pg, float *g0)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    g0[i] += pg[i];
}

__global__ static void
ker_add_bwd_child1(int n0, int n1, const float *pg, float *g1)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n1;
       i += blockDim.x * gridDim.x) {
    float s = 0.0f;
    for (int j = i; j < n0; j += n1)
      s += pg[j];
    g1[i] += s;
  }
}

// --- mul (op 2): element-wise multiplication with broadcasting ---
__global__ static void
ker_mul_fwd(int n0, int n1, const float *x0, const float *x1, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n0;
       i += blockDim.x * gridDim.x)
    px[i] = x0[i] * x1[i % n1];
}

__global__ static void
ker_mul_bwd_child0(int n0, int n1, const float *pg, const float *x1, float *g0)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n0;
       i += blockDim.x * gridDim.x)
    g0[i] += pg[i] * x1[i % n1];
}

__global__ static void
ker_mul_bwd_child1(int n0, int n1, const float *pg, const float *x0, float *g1)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n1;
       i += blockDim.x * gridDim.x) {
    float s = 0.0f;
    for (int j = i; j < n0; j += n1)
      s += pg[j] * x0[j];
    g1[i] += s;
  }
}

// --- sigm (op 6): sigmoid ---
__global__ static void
ker_sigm_fwd(int n, const float *qx, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    px[i] = 1.0f / (1.0f + expf(-qx[i]));
}

__global__ static void
ker_sigm_bwd(int n, const float *pg, const float *px, float *qg)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    qg[i] += pg[i] * (px[i] * (1.0f - px[i]));
}

// --- tanh (op 7) ---
__global__ static void
ker_tanh_fwd(int n, const float *qx, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x) {
    if (qx[i] < -20.0f)
      px[i] = -1.0f;
    else {
      float y = expf(-2.0f * qx[i]);
      px[i] = (1.0f - y) / (1.0f + y);
    }
  }
}

__global__ static void
ker_tanh_bwd(int n, const float *pg, const float *px, float *qg)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    qg[i] += pg[i] * (1.0f - px[i] * px[i]);
}

// --- relu (op 8) ---
__global__ static void
ker_relu_fwd(int n, const float *qx, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    px[i] = qx[i] > 0.0f ? qx[i] : 0.0f;
}

__global__ static void
ker_relu_bwd(int n, const float *pg, const float *qx, float *qg)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    if (qx[i] > 0.0f)
      qg[i] += pg[i];
}

// --- 1minus (op 11): p = 1 - child ---
__global__ static void
ker_1minus_fwd(int n, const float *qx, float *px)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    px[i] = 1.0f - qx[i];
}

__global__ static void
ker_1minus_bwd(int n, const float *pg, float *qg)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    qg[i] -= pg[i];
}

// --- softmax (op 14) ---
__global__ static void
ker_softmax_fwd(int d0, int n1, const float *qx, float *px)
{
  int row = blockIdx.x;
  if (row >= d0) return;
  const float *x = qx + row * n1;
  float *y = px + row * n1;

  float mx = -FLT_MAX;
  for (int i = 0; i < n1; ++i)
    mx = fmaxf(mx, x[i]);

  float s = 0.0f;
  for (int i = 0; i < n1; ++i) {
    y[i] = expf(x[i] - mx);
    s += y[i];
  }
  s = 1.0f / s;
  for (int i = 0; i < n1; ++i)
    y[i] *= s;
}

__global__ static void
ker_softmax_bwd(int d0, int n1, const float *pg, const float *px, float *qg)
{
  int row = blockIdx.x;
  if (row >= d0) return;
  const float *g = pg + row * n1;
  const float *y = px + row * n1;
  float *h = qg + row * n1;

  float s = 0.0f;
  for (int i = 0; i < n1; ++i)
    s += g[i] * y[i];
  for (int i = 0; i < n1; ++i)
    h[i] += y[i] * (g[i] - s);
}

// --- mse (op 29): mean square error ---
__global__ static void
ker_mse_fwd(int n, const float *y1, const float *y0, float *cost)
{
  __shared__ float sdata[KANN_CU_THREADS];
  int tid = threadIdx.x;
  float sum = 0.0f;
  for (int i = tid; i < n; i += blockDim.x) {
    float d = y1[i] - y0[i];
    sum += d * d;
  }
  sdata[tid] = sum;
  __syncthreads();

  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s)
      sdata[tid] += sdata[tid + s];
    __syncthreads();
  }
  if (tid == 0)
    cost[0] = sdata[0] / n;
}

__global__ static void
ker_mse_bwd(int n, const float *pg, const float *y1, const float *y0, float *y1g)
{
  float t = 2.0f * pg[0] / n;
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    y1g[i] += t * (y1[i] - y0[i]);
}

// --- stdnorm (op 32): layer normalization ---
// Each row of length n is independently normalized: subtract mean, divide by std.
// si_out stores per-row 1/std for backward pass.
// One block per row, shared-memory reduction.
__global__ static void
ker_stdnorm_fwd(int m, int n, const float *qx, float *px, float *si_out)
{
  int row = blockIdx.x;
  if (row >= m) return;
  const float *in = qx + row * n;
  float *out = px + row * n;

  // Compute mean via shared reduction
  extern __shared__ float smem[];
  float sum = 0.0f;
  for (int i = threadIdx.x; i < n; i += blockDim.x)
    sum += in[i];
  smem[threadIdx.x] = sum;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) smem[threadIdx.x] += smem[threadIdx.x + s];
    __syncthreads();
  }
  float avg = smem[0] / n;

  // Subtract mean
  for (int i = threadIdx.x; i < n; i += blockDim.x)
    out[i] = in[i] - avg;
  __syncthreads();

  // Compute variance
  sum = 0.0f;
  for (int i = threadIdx.x; i < n; i += blockDim.x)
    sum += out[i] * out[i];
  smem[threadIdx.x] = sum;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) smem[threadIdx.x] += smem[threadIdx.x + s];
    __syncthreads();
  }
  float var = smem[0] / n;
  float std_inv = (var == 0.0f) ? 1.0f : rsqrtf(var);

  // Scale by 1/std
  for (int i = threadIdx.x; i < n; i += blockDim.x)
    out[i] *= std_inv;

  if (threadIdx.x == 0)
    si_out[row] = std_inv;
}

// Backward: qg[i] += std_inv * (pg[i] - mean(pg) - px[i] * mean(px * pg))
__global__ static void
ker_stdnorm_bwd(int m, int n, const float *pg, const float *px,
  const float *si, float *qg)
{
  int row = blockIdx.x;
  if (row >= m) return;
  const float *pg_r = pg + row * n;
  const float *px_r = px + row * n;
  float *qg_r = qg + row * n;
  float std_inv = si[row];

  extern __shared__ float smem[];
  // smem[0..blockDim-1] = sum of pg, smem[blockDim..2*blockDim-1] = sum of px*pg
  float *smem_s = smem;
  float *smem_t = smem + blockDim.x;

  float s_sum = 0.0f, t_sum = 0.0f;
  for (int i = threadIdx.x; i < n; i += blockDim.x) {
    s_sum += pg_r[i];
    t_sum += px_r[i] * pg_r[i];
  }
  smem_s[threadIdx.x] = s_sum;
  smem_t[threadIdx.x] = t_sum;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      smem_s[threadIdx.x] += smem_s[threadIdx.x + s];
      smem_t[threadIdx.x] += smem_t[threadIdx.x + s];
    }
    __syncthreads();
  }
  float s_avg = smem_s[0] / n;
  float t_avg = smem_t[0] / n;

  for (int i = threadIdx.x; i < n; i += blockDim.x)
    qg_r[i] += std_inv * (pg_r[i] - s_avg - px_r[i] * t_avg);
}

// ============================================================
// Support kernels
// ============================================================

__global__ static void
ker_zero(int n, float *buf)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    buf[i] = 0.0f;
}

__global__ static void
ker_saxpy(int n, float a, const float *x, float *y)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x)
    y[i] += a * x[i];
}

// RMSprop: r = (1-decay)*g^2 + decay*r; t -= lr/sqrt(r+eps)*g
__global__ static void
ker_rmsprop(int n, float lr, float decay, const float *g, float *t, float *r)
{
  float d1 = 1.0f - decay;
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
       i += blockDim.x * gridDim.x) {
    r[i] = d1 * g[i] * g[i] + decay * r[i];
    t[i] -= lr * rsqrtf(1e-6f + r[i]) * g[i];
  }
}

// ============================================================
// Helper: compute grid size
// ============================================================
static inline int cu_nblocks(int n)
{
  return (n + KANN_CU_THREADS - 1) / KANN_CU_THREADS;
}

// ============================================================
// Graph construction: flatten kad_node_t** into kann_cu_graph
// ============================================================

static int host_kad_len(const kad_node_t *p)
{
  int n = 1;
  for (int i = 0; i < p->n_d; ++i) n *= p->d[i];
  return n;
}

struct kann_cu_graph*
kann_cu_graph_new(kann_t *ann, int max_batch_size)
{
  int n = ann->n;
  kad_node_t **v = ann->v;

  // Sync host graph to max batch size so buffer allocation is large enough
  kad_sync_dim(n, v, max_batch_size);

  struct kann_cu_graph *cg =
    (struct kann_cu_graph *)calloc(1, sizeof(struct kann_cu_graph));
  cg->n_node = n;
  cg->max_batch_size = max_batch_size;
  cg->cur_batch_size = max_batch_size;

  // Build host-side node metadata array (persisted as host mirror)
  cg->h_nodes =
    (struct kann_cu_node *)calloc(n, sizeof(struct kann_cu_node));

  // First pass: assign offsets into flat x/g buffers
  int x_total = 0, g_total = 0;
  for (int i = 0; i < n; ++i) {
    kad_node_t *p = v[i];
    int len = host_kad_len(p);

    cg->h_nodes[i].op = p->op;
    cg->h_nodes[i].n_d = p->n_d;
    cg->h_nodes[i].flag = p->flag;
    cg->h_nodes[i].n_child = p->n_child;
    cg->h_nodes[i].ext_flag = p->ext_flag;
    cg->h_nodes[i].ext_label = p->ext_label;
    cg->h_nodes[i].len = len;
    for (int d = 0; d < KAD_MAX_DIM; ++d)
      cg->h_nodes[i].d[d] = p->d[d];

    // Integer parameter for select op
    cg->h_nodes[i].ptr_i32 = 0;
    if (p->ptr && (p->op == 12))
      cg->h_nodes[i].ptr_i32 = *(int32_t *)p->ptr;

    // Pre-linkage (RNN recurrence)
    cg->h_nodes[i].pre_idx = -1;
    // si_off for stdnorm (set in second pass below)
    cg->h_nodes[i].si_off = -1;

    // x offset (allocated at max batch size)
    cg->h_nodes[i].x_off = x_total;
    x_total += (len > 0 ? len : 1);

    // g offset: only for nodes that participate in backprop
    if (p->flag & KAD_VAR) {
      cg->h_nodes[i].g_off = g_total;
      g_total += (len > 0 ? len : 1);
    } else {
      cg->h_nodes[i].g_off = -1;
    }

    // Resolve child indices
    for (int c = 0; c < p->n_child && c < KANN_CU_MAX_CHILD; ++c) {
      cg->h_nodes[i].child_idx[c] = -1;
      for (int j = 0; j < n; ++j) {
        if (v[j] == p->child[c]) {
          cg->h_nodes[i].child_idx[c] = j;
          break;
        }
      }
      assert(cg->h_nodes[i].child_idx[c] >= 0);
    }
  }

  cg->x_total = x_total;
  cg->g_total = g_total;

  // Find special nodes and build per-variable/per-constant scatter info.
  // Variables are NOT contiguous in the x buffer (internal nodes are
  // interleaved in topological order), so we track each variable node
  // individually for correct upload/download/RMSprop.
  cg->in_node_idx = -1;
  cg->truth_node_idx = -1;
  cg->cost_node_idx = -1;
  cg->out_node_idx = -1;

  // First pass: count variable and constant nodes
  int n_var_nodes = 0, n_const_nodes = 0;
  int n_var = 0, n_const = 0;
  for (int i = 0; i < n; ++i) {
    kad_node_t *p = v[i];
    if (p->ext_flag & KANN_F_IN)   cg->in_node_idx = i;
    if (p->ext_flag & KANN_F_TRUTH) cg->truth_node_idx = i;
    if (p->ext_flag & KANN_F_COST)  cg->cost_node_idx = i;
    if (p->ext_flag & KANN_F_OUT)   cg->out_node_idx = i;

    if (kad_is_var(p)) {
      n_var_nodes++;
      n_var += cg->h_nodes[i].len;
    }
    if (kad_is_const(p)) {
      n_const_nodes++;
      n_const += cg->h_nodes[i].len;
    }
  }

  // Build scatter arrays
  cg->h_vars = (struct kann_cu_var_node *)malloc(
    n_var_nodes * sizeof(struct kann_cu_var_node));
  cg->n_var_nodes = n_var_nodes;
  cg->n_var = n_var;

  cg->h_consts = (struct kann_cu_var_node *)malloc(
    n_const_nodes * sizeof(struct kann_cu_var_node));
  cg->n_const_nodes = n_const_nodes;
  cg->n_const = n_const;

  int vi = 0, ci = 0, r_off = 0;
  for (int i = 0; i < n; ++i) {
    kad_node_t *p = v[i];
    if (kad_is_var(p)) {
      int len = cg->h_nodes[i].len;
      cg->h_vars[vi].x_off = cg->h_nodes[i].x_off;
      cg->h_vars[vi].g_off = cg->h_nodes[i].g_off;
      cg->h_vars[vi].r_off = r_off;
      cg->h_vars[vi].len = len;
      r_off += len;
      vi++;
    }
    if (kad_is_const(p)) {
      int len = cg->h_nodes[i].len;
      cg->h_consts[ci].x_off = cg->h_nodes[i].x_off;
      cg->h_consts[ci].g_off = -1;
      cg->h_consts[ci].r_off = -1;
      cg->h_consts[ci].len = len;
      ci++;
    }
  }

  // Build forward/backward order (internal nodes only)
  int n_internal = 0;
  for (int i = 0; i < n; ++i)
    if (v[i]->n_child > 0) n_internal++;

  cg->h_fwd_order = (int *)malloc(n_internal * sizeof(int));
  cg->h_bwd_order = (int *)malloc(n_internal * sizeof(int));
  int fi = 0;
  for (int i = 0; i < n; ++i)
    if (v[i]->n_child > 0)
      cg->h_fwd_order[fi++] = i;
  for (int i = 0; i < n_internal; ++i)
    cg->h_bwd_order[i] = cg->h_fwd_order[n_internal - 1 - i];

  cg->n_internal = n_internal;

  // Resolve pre-linkage (RNN recurrence) and stdnorm cache offsets
  int n_pre_pairs = 0;
  int si_total = 0;
  for (int i = 0; i < n; ++i) {
    kad_node_t *p = v[i];
    // Resolve pre pointer to node index
    if (p->pre) {
      for (int j = 0; j < n; ++j) {
        if (v[j] == p->pre) {
          cg->h_nodes[i].pre_idx = j;
          n_pre_pairs++;
          break;
        }
      }
    }
    // Assign stdnorm cache offset (m floats per stdnorm node)
    if (p->op == 32 && p->n_child > 0) {
      int child_n = cg->h_nodes[cg->h_nodes[i].child_idx[0]].d[
        cg->h_nodes[cg->h_nodes[i].child_idx[0]].n_d - 1];
      int m = cg->h_nodes[i].len / child_n;
      cg->h_nodes[i].si_off = si_total;
      si_total += m;
    }
  }

  // Build pre-pair list: (output_node_idx, h0_node_idx)
  cg->n_pre_pairs = n_pre_pairs;
  cg->h_pre_pairs = NULL;
  if (n_pre_pairs > 0) {
    cg->h_pre_pairs = (int *)malloc(2 * n_pre_pairs * sizeof(int));
    int pi = 0;
    for (int i = 0; i < n; ++i) {
      if (cg->h_nodes[i].pre_idx >= 0) {
        cg->h_pre_pairs[2*pi] = i;                     // output node
        cg->h_pre_pairs[2*pi+1] = cg->h_nodes[i].pre_idx; // h0 node
        pi++;
      }
    }
  }

  cg->stdnorm_si_total = si_total;

  // Allocate device memory
  cg->nodes = (struct kann_cu_node *)gkyl_cu_malloc(n * sizeof(struct kann_cu_node));
  cg->x = (float *)gkyl_cu_malloc(x_total * sizeof(float));
  cg->g = (float *)gkyl_cu_malloc(g_total * sizeof(float));
  cg->r = (float *)gkyl_cu_malloc(n_var * sizeof(float));
  cg->stdnorm_si = si_total > 0
    ? (float *)gkyl_cu_malloc(si_total * sizeof(float)) : NULL;

  // Upload node metadata to device
  gkyl_cu_memcpy(cg->nodes, cg->h_nodes, n * sizeof(struct kann_cu_node),
    GKYL_CU_MEMCPY_H2D);

  // Zero out device buffers
  cudaMemset(cg->x, 0, x_total * sizeof(float));
  cudaMemset(cg->g, 0, g_total * sizeof(float));
  cudaMemset(cg->r, 0, n_var * sizeof(float));

  // Upload initial values (variables, constants, and any pre-set x)
  for (int i = 0; i < n; ++i) {
    kad_node_t *p = v[i];
    int len = cg->h_nodes[i].len;
    if (p->x && len > 0 && kad_is_ext(p)) {
      gkyl_cu_memcpy(cg->x + cg->h_nodes[i].x_off, p->x,
        len * sizeof(float), GKYL_CU_MEMCPY_H2D);
    }
  }

  // Restore host graph to batch_size=1 so it is not left in a modified state
  kad_sync_dim(n, v, 1);

  return cg;
}

void
kann_cu_graph_free(struct kann_cu_graph *cg)
{
  if (!cg) return;
  gkyl_cu_free(cg->nodes);
  gkyl_cu_free(cg->x);
  gkyl_cu_free(cg->g);
  gkyl_cu_free(cg->r);
  if (cg->stdnorm_si) gkyl_cu_free(cg->stdnorm_si);
  free(cg->h_nodes);
  free(cg->h_vars);
  free(cg->h_consts);
  free(cg->h_fwd_order);
  free(cg->h_bwd_order);
  free(cg->h_pre_pairs);
  free(cg);
}

void
kann_cu_graph_upload_vars(struct kann_cu_graph *cg, const kann_t *ann)
{
  if (!ann->x) return;
  int host_off = 0;
  for (int i = 0; i < cg->n_var_nodes; ++i) {
    struct kann_cu_var_node *vn = &cg->h_vars[i];
    gkyl_cu_memcpy(cg->x + vn->x_off, ann->x + host_off,
      vn->len * sizeof(float), GKYL_CU_MEMCPY_H2D);
    host_off += vn->len;
  }
}

void
kann_cu_graph_download_vars(struct kann_cu_graph *cg, kann_t *ann)
{
  if (!ann->x) return;
  int host_off = 0;
  for (int i = 0; i < cg->n_var_nodes; ++i) {
    struct kann_cu_var_node *vn = &cg->h_vars[i];
    gkyl_cu_memcpy(ann->x + host_off, cg->x + vn->x_off,
      vn->len * sizeof(float), GKYL_CU_MEMCPY_D2H);
    host_off += vn->len;
  }
}

void
kann_cu_graph_upload_consts(struct kann_cu_graph *cg, const kann_t *ann)
{
  if (!ann->c) return;
  int host_off = 0;
  for (int i = 0; i < cg->n_const_nodes; ++i) {
    struct kann_cu_var_node *cn = &cg->h_consts[i];
    gkyl_cu_memcpy(cg->x + cn->x_off, ann->c + host_off,
      cn->len * sizeof(float), GKYL_CU_MEMCPY_H2D);
    host_off += cn->len;
  }
}

// Update node dimensions for a new batch size. Uses the host kann_t
// to run proper kad_sync_dim propagation, then copies updated dims
// to the host mirror and re-uploads to device.
void
kann_cu_sync_dim(struct kann_cu_graph *cg, kann_t *ann, int batch_size)
{
  if (batch_size == cg->cur_batch_size) return;
  assert(batch_size <= cg->max_batch_size);

  // Run proper dimension propagation on host
  kad_sync_dim(ann->n, ann->v, batch_size);

  // Update host mirror dims and lens (offsets stay fixed)
  for (int i = 0; i < cg->n_node; ++i) {
    kad_node_t *p = ann->v[i];
    cg->h_nodes[i].n_d = p->n_d;
    for (int d = 0; d < KAD_MAX_DIM; ++d)
      cg->h_nodes[i].d[d] = p->d[d];
    cg->h_nodes[i].len = host_kad_len(p);
  }

  // Re-upload metadata to device
  gkyl_cu_memcpy(cg->nodes, cg->h_nodes,
    cg->n_node * sizeof(struct kann_cu_node), GKYL_CU_MEMCPY_H2D);

  cg->cur_batch_size = batch_size;
}

// ============================================================
// Forward pass dispatch: host loop, one kernel per node
// ============================================================

static void
dispatch_forward(struct kann_cu_graph *cg, cublasHandle_t cublas_h,
  struct kann_cu_node *hn, int idx)
{
  struct kann_cu_node *p = &hn[idx];
  int len = p->len;
  float *px = cg->x + p->x_off;

  if (p->n_child == 0) return;

  float *cx[KANN_CU_MAX_CHILD];
  int clen[KANN_CU_MAX_CHILD];
  for (int c = 0; c < p->n_child && c < KANN_CU_MAX_CHILD; ++c) {
    int ci = p->child_idx[c];
    cx[c] = cg->x + hn[ci].x_off;
    clen[c] = hn[ci].len;
  }

  int nb;

  switch (p->op) {
  case 1: // add
    nb = cu_nblocks(clen[0]);
    ker_add_fwd<<<nb, KANN_CU_THREADS>>>(clen[0], clen[1], cx[0], cx[1], px);
    break;

  case 2: // mul
    nb = cu_nblocks(clen[0]);
    ker_mul_fwd<<<nb, KANN_CU_THREADS>>>(clen[0], clen[1], cx[0], cx[1], px);
    break;

  case 3: { // cmul: Y = X * W^T via cublasSgemm
    struct kann_cu_node *q0 = &hn[p->child_idx[0]];
    struct kann_cu_node *q1 = &hn[p->child_idx[1]];
    int n_col = q0->d[q0->n_d - 1] > q1->d[q1->n_d - 1]
      ? q0->d[q0->n_d - 1] : q1->d[q1->n_d - 1];
    int n_a_col = 1, n_b_col = 1;
    for (int i = q0->n_d - 1; i >= 0; --i) { if (n_a_col < n_col) n_a_col *= q0->d[i]; }
    for (int i = q1->n_d - 1; i >= 0; --i) { if (n_b_col < n_col) n_b_col *= q1->d[i]; }
    int n_a_row = q0->len / n_a_col;
    int n_b_row = q1->len / n_b_col;

    // Row-major C = A * B^T  <==>  col-major C^T = B * A^T
    float alpha = 1.0f, beta = 0.0f;
    cublasSgemm(cublas_h,
      CUBLAS_OP_T, CUBLAS_OP_N,
      n_b_row, n_a_row, n_col,
      &alpha,
      cx[1], n_b_col,
      cx[0], n_a_col,
      &beta,
      px, n_b_row);
    break;
  }

  case 6: // sigm
    nb = cu_nblocks(len);
    ker_sigm_fwd<<<nb, KANN_CU_THREADS>>>(len, cx[0], px);
    break;

  case 7: // tanh
    nb = cu_nblocks(len);
    ker_tanh_fwd<<<nb, KANN_CU_THREADS>>>(len, cx[0], px);
    break;

  case 8: // relu
    nb = cu_nblocks(len);
    ker_relu_fwd<<<nb, KANN_CU_THREADS>>>(len, cx[0], px);
    break;

  case 10: { // avg: mean of n_child tensors
    float inv = 1.0f / p->n_child;
    nb = cu_nblocks(len);
    ker_zero<<<nb, KANN_CU_THREADS>>>(len, px);
    for (int c = 0; c < p->n_child; ++c) {
      int ci = p->child_idx[c];
      ker_saxpy<<<nb, KANN_CU_THREADS>>>(len, inv, cg->x + hn[ci].x_off, px);
    }
    break;
  }

  case 11: // 1minus
    nb = cu_nblocks(len);
    ker_1minus_fwd<<<nb, KANN_CU_THREADS>>>(len, cx[0], px);
    break;

  case 12: { // select
    int which = p->ptr_i32;
    if (which < 0) which += p->n_child;
    int ci = p->child_idx[which];
    int slen = hn[ci].len;
    cudaMemcpy(px, cg->x + hn[ci].x_off, slen * sizeof(float),
      cudaMemcpyDeviceToDevice);
    break;
  }

  case 14: { // softmax
    int n1 = hn[p->child_idx[0]].d[hn[p->child_idx[0]].n_d - 1];
    int d0 = clen[0] / n1;
    ker_softmax_fwd<<<d0, 1>>>(d0, n1, cx[0], px);
    break;
  }

  case 29: // mse
    ker_mse_fwd<<<1, KANN_CU_THREADS>>>(clen[0], cx[0], cx[1], px);
    break;

  case 32: { // stdnorm (layer normalization)
    struct kann_cu_node *q = &hn[p->child_idx[0]];
    int sn = q->d[q->n_d - 1];
    int sm = len / sn;
    int threads = KANN_CU_THREADS;
    if (threads > sn) threads = sn;
    // Round down to power of 2 for shared-memory reduction
    int t2 = 1;
    while (t2 * 2 <= threads) t2 *= 2;
    ker_stdnorm_fwd<<<sm, t2, t2 * sizeof(float)>>>(
      sm, sn, cx[0], px, cg->stdnorm_si + p->si_off);
    break;
  }

  default:
    fprintf(stderr, "kann_cu_forward: unimplemented op %d\n", p->op);
    assert(0);
  }
}

static void
dispatch_backward(struct kann_cu_graph *cg, cublasHandle_t cublas_h,
  struct kann_cu_node *hn, int idx)
{
  struct kann_cu_node *p = &hn[idx];
  float *pg = cg->g + p->g_off;

  if (p->n_child == 0) return;

  float *cx[KANN_CU_MAX_CHILD], *cg_arr[KANN_CU_MAX_CHILD];
  int clen[KANN_CU_MAX_CHILD];
  uint8_t cflag[KANN_CU_MAX_CHILD];
  for (int c = 0; c < p->n_child && c < KANN_CU_MAX_CHILD; ++c) {
    int ci = p->child_idx[c];
    cx[c] = cg->x + hn[ci].x_off;
    cg_arr[c] = (hn[ci].g_off >= 0) ? (cg->g + hn[ci].g_off) : NULL;
    clen[c] = hn[ci].len;
    cflag[c] = hn[ci].flag;
  }

  int nb;

  switch (p->op) {
  case 1: // add
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      nb = cu_nblocks(clen[0]);
      ker_add_bwd_child0<<<nb, KANN_CU_THREADS>>>(clen[0], pg, cg_arr[0]);
    }
    if ((cflag[1] & KAD_VAR) && cg_arr[1]) {
      nb = cu_nblocks(clen[1]);
      ker_add_bwd_child1<<<nb, KANN_CU_THREADS>>>(clen[0], clen[1], pg, cg_arr[1]);
    }
    break;

  case 2: // mul
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      nb = cu_nblocks(clen[0]);
      ker_mul_bwd_child0<<<nb, KANN_CU_THREADS>>>(clen[0], clen[1], pg, cx[1], cg_arr[0]);
    }
    if ((cflag[1] & KAD_VAR) && cg_arr[1]) {
      nb = cu_nblocks(clen[1]);
      ker_mul_bwd_child1<<<nb, KANN_CU_THREADS>>>(clen[0], clen[1], pg, cx[0], cg_arr[1]);
    }
    break;

  case 3: { // cmul backward
    struct kann_cu_node *q0 = &hn[p->child_idx[0]];
    struct kann_cu_node *q1 = &hn[p->child_idx[1]];
    int n_col = q0->d[q0->n_d - 1] > q1->d[q1->n_d - 1]
      ? q0->d[q0->n_d - 1] : q1->d[q1->n_d - 1];
    int n_a_col = 1, n_b_col = 1;
    for (int i = q0->n_d - 1; i >= 0; --i) { if (n_a_col < n_col) n_a_col *= q0->d[i]; }
    for (int i = q1->n_d - 1; i >= 0; --i) { if (n_b_col < n_col) n_b_col *= q1->d[i]; }
    int n_a_row = q0->len / n_a_col;
    int n_b_row = q1->len / n_b_col;
    float alpha = 1.0f, beta = 1.0f;

    // G_x += G_y * W
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      cublasSgemm(cublas_h,
        CUBLAS_OP_N, CUBLAS_OP_N,
        n_a_col, n_a_row, n_b_row,
        &alpha,
        cx[1], n_b_col,
        pg, n_b_row,
        &beta,
        cg_arr[0], n_a_col);
    }
    // G_w += G_y^T * X
    if ((cflag[1] & KAD_VAR) && cg_arr[1]) {
      cublasSgemm(cublas_h,
        CUBLAS_OP_N, CUBLAS_OP_T,
        n_b_col, n_b_row, n_a_row,
        &alpha,
        cx[0], n_a_col,
        pg, n_b_row,
        &beta,
        cg_arr[1], n_b_col);
    }
    break;
  }

  case 6: // sigm
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      float *px = cg->x + p->x_off;
      nb = cu_nblocks(p->len);
      ker_sigm_bwd<<<nb, KANN_CU_THREADS>>>(p->len, pg, px, cg_arr[0]);
    }
    break;

  case 7: // tanh
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      float *px = cg->x + p->x_off;
      nb = cu_nblocks(p->len);
      ker_tanh_bwd<<<nb, KANN_CU_THREADS>>>(p->len, pg, px, cg_arr[0]);
    }
    break;

  case 8: // relu
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      nb = cu_nblocks(p->len);
      ker_relu_bwd<<<nb, KANN_CU_THREADS>>>(p->len, pg, cx[0], cg_arr[0]);
    }
    break;

  case 10: { // avg backward
    float inv = 1.0f / p->n_child;
    nb = cu_nblocks(p->len);
    for (int c = 0; c < p->n_child; ++c) {
      int ci = p->child_idx[c];
      if ((hn[ci].flag & KAD_VAR) && hn[ci].g_off >= 0)
        ker_saxpy<<<nb, KANN_CU_THREADS>>>(p->len, inv, pg, cg->g + hn[ci].g_off);
    }
    break;
  }

  case 11: // 1minus
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      nb = cu_nblocks(p->len);
      ker_1minus_bwd<<<nb, KANN_CU_THREADS>>>(p->len, pg, cg_arr[0]);
    }
    break;

  case 12: { // select
    int which = p->ptr_i32;
    if (which < 0) which += p->n_child;
    int ci = p->child_idx[which];
    if ((hn[ci].flag & KAD_VAR) && hn[ci].g_off >= 0) {
      int slen = hn[ci].len;
      nb = cu_nblocks(slen);
      ker_saxpy<<<nb, KANN_CU_THREADS>>>(slen, 1.0f, pg, cg->g + hn[ci].g_off);
    }
    break;
  }

  case 14: { // softmax
    struct kann_cu_node *q = &hn[p->child_idx[0]];
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      int n1 = q->d[q->n_d - 1];
      int d0 = q->len / n1;
      float *px = cg->x + p->x_off;
      ker_softmax_bwd<<<d0, 1>>>(d0, n1, pg, px, cg_arr[0]);
    }
    break;
  }

  case 29: // mse
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      nb = cu_nblocks(clen[0]);
      ker_mse_bwd<<<nb, KANN_CU_THREADS>>>(clen[0], pg, cx[0], cx[1], cg_arr[0]);
    }
    break;

  case 32: { // stdnorm backward
    if ((cflag[0] & KAD_VAR) && cg_arr[0]) {
      struct kann_cu_node *q = &hn[p->child_idx[0]];
      int sn = q->d[q->n_d - 1];
      int sm = p->len / sn;
      float *px = cg->x + p->x_off;
      int threads = KANN_CU_THREADS;
      if (threads > sn) threads = sn;
      int t2 = 1;
      while (t2 * 2 <= threads) t2 *= 2;
      ker_stdnorm_bwd<<<sm, t2, 2 * t2 * sizeof(float)>>>(
        sm, sn, pg, px, cg->stdnorm_si + p->si_off, cg_arr[0]);
    }
    break;
  }

  default:
    fprintf(stderr, "kann_cu_backward: unimplemented op %d\n", p->op);
    assert(0);
  }
}

// ============================================================
// Public API implementations
// ============================================================

void
kann_cu_forward(struct kann_cu_graph *cg, void *cublas_h, int to)
{
  struct kann_cu_node *hn = cg->h_nodes; // use persistent host mirror

  int stop = (to >= 0) ? to : cg->n_node - 1;
  for (int fi = 0; fi < cg->n_internal; ++fi) {
    int idx = cg->h_fwd_order[fi];
    if (idx > stop) break;
    dispatch_forward(cg, (cublasHandle_t)cublas_h, hn, idx);
  }
}

void
kann_cu_backward(struct kann_cu_graph *cg, void *cublas_h, int from)
{
  struct kann_cu_node *hn = cg->h_nodes;

  // Zero all gradients
  if (cg->g_total > 0)
    cudaMemset(cg->g, 0, cg->g_total * sizeof(float));

  // Set cost gradient to 1.0
  if (from >= 0 && hn[from].g_off >= 0) {
    float one = 1.0f;
    gkyl_cu_memcpy(cg->g + hn[from].g_off, &one, sizeof(float),
      GKYL_CU_MEMCPY_H2D);
  }

  // Backward pass
  int start = (from >= 0) ? from : cg->n_node - 1;
  for (int bi = 0; bi < cg->n_internal; ++bi) {
    int idx = cg->h_bwd_order[bi];
    if (idx > start) continue;
    dispatch_backward(cg, (cublasHandle_t)cublas_h, hn, idx);
  }
}

void
kann_cu_zero_grad(struct kann_cu_graph *cg)
{
  if (cg->g_total > 0)
    cudaMemset(cg->g, 0, cg->g_total * sizeof(float));
}

void
kann_cu_rmsprop(struct kann_cu_graph *cg, float lr, float decay)
{
  // Variables are scattered in the x buffer, so iterate per variable node
  for (int i = 0; i < cg->n_var_nodes; ++i) {
    struct kann_cu_var_node *vn = &cg->h_vars[i];
    int nb = cu_nblocks(vn->len);
    ker_rmsprop<<<nb, KANN_CU_THREADS>>>(
      vn->len, lr, decay,
      cg->g + vn->g_off,     // this variable's gradients (contiguous in g)
      cg->x + vn->x_off,     // this variable's values (scattered in x)
      cg->r + vn->r_off      // this variable's RMSprop accumulator
    );
  }
}

void
kann_cu_feed_input(struct kann_cu_graph *cg, int batch_size, const float *x_host)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->in_node_idx];
  gkyl_cu_memcpy(cg->x + hn->x_off, x_host,
    hn->len * sizeof(float), GKYL_CU_MEMCPY_H2D);
}

void
kann_cu_feed_truth(struct kann_cu_graph *cg, int batch_size, const float *y_host)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->truth_node_idx];
  gkyl_cu_memcpy(cg->x + hn->x_off, y_host,
    hn->len * sizeof(float), GKYL_CU_MEMCPY_H2D);
}

// Gather kernel: pick rows from src using indices and write contiguous
__global__ static void
ker_gather_rows(int total, int stride, const float *src,
  const int *indices, int offset, float *dst)
{
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < total;
       i += blockDim.x * gridDim.x) {
    int row = i / stride;
    int col = i % stride;
    dst[i] = src[indices[offset + row] * stride + col];
  }
}

void
kann_cu_gather_input(struct kann_cu_graph *cg, const float *data_d,
  const int *indices_d, int offset, int batch_size, int stride)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->in_node_idx];
  int total = batch_size * stride;
  int nb = cu_nblocks(total);
  ker_gather_rows<<<nb, KANN_CU_THREADS>>>(
    total, stride, data_d, indices_d, offset, cg->x + hn->x_off);
}

void
kann_cu_gather_truth(struct kann_cu_graph *cg, const float *data_d,
  const int *indices_d, int offset, int batch_size, int stride)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->truth_node_idx];
  int total = batch_size * stride;
  int nb = cu_nblocks(total);
  ker_gather_rows<<<nb, KANN_CU_THREADS>>>(
    total, stride, data_d, indices_d, offset, cg->x + hn->x_off);
}

void
kann_cu_feed_input_dev(struct kann_cu_graph *cg, int batch_size,
  const float *x_dev, int stride)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->in_node_idx];
  cudaMemcpy(cg->x + hn->x_off, x_dev,
    batch_size * stride * sizeof(float), cudaMemcpyDeviceToDevice);
}

float
kann_cu_get_cost(const struct kann_cu_graph *cg)
{
  struct kann_cu_node *hn = &cg->h_nodes[cg->cost_node_idx];
  float cost;
  gkyl_cu_memcpy(&cost, cg->x + hn->x_off, sizeof(float), GKYL_CU_MEMCPY_D2H);
  return cost;
}

void
kann_cu_get_output(const struct kann_cu_graph *cg, int batch_size, float *out_host)
{
  int idx = cg->out_node_idx >= 0 ? cg->out_node_idx : cg->cost_node_idx;
  struct kann_cu_node *hn = &cg->h_nodes[idx];
  gkyl_cu_memcpy(out_host, cg->x + hn->x_off,
    hn->len * sizeof(float), GKYL_CU_MEMCPY_D2H);
}

void
kann_cu_apply_pre(struct kann_cu_graph *cg)
{
  for (int i = 0; i < cg->n_pre_pairs; ++i) {
    int out_idx = cg->h_pre_pairs[2*i];
    int h0_idx  = cg->h_pre_pairs[2*i+1];
    struct kann_cu_node *out_n = &cg->h_nodes[out_idx];
    struct kann_cu_node *h0_n  = &cg->h_nodes[h0_idx];
    int len = h0_n->len < out_n->len ? h0_n->len : out_n->len;
    cudaMemcpy(cg->x + h0_n->x_off, cg->x + out_n->x_off,
      len * sizeof(float), cudaMemcpyDeviceToDevice);
  }
}
