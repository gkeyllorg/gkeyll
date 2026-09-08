#pragma once

// Flattened graph representation for GPU execution of KANN networks.
// Serializes the pointer-chasing kad_node_t** graph into contiguous arrays
// that can live entirely in device memory.

#include <gkyl_kann_net.h>
#include <gkyl_alloc_flags_priv.h>

#include <stdint.h>

// Maximum children per node in the flattened representation
#define KANN_CU_MAX_CHILD 4

// Per-node metadata in the flattened graph (fixed-size, GPU-friendly)
struct kann_cu_node {
  uint16_t op; // operator index
  uint8_t n_d; // number of dimensions
  uint8_t flag; // KAD_VAR, KAD_CONST, etc.
  int32_t n_child; // number of children
  int32_t child_idx[KANN_CU_MAX_CHILD]; // indices of children in node array
  int32_t d[KAD_MAX_DIM]; // dimensions
  int32_t x_off; // offset into flat x (value) buffer
  int32_t g_off; // offset into flat g (gradient) buffer
  int32_t len; // current number of floats (product of dims)
  int32_t ptr_i32; // integer parameter (e.g. select index)
  uint32_t ext_flag; // external flags (KANN_F_IN, etc.)
  int32_t ext_label; // external label
  int32_t pre_idx; // index of pre-linked node (-1 if none)
  int32_t si_off; // offset into stdnorm_si buffer (-1 if not stdnorm)
};

// Per-variable-node metadata for scattered upload/download/RMSprop
struct kann_cu_var_node {
  int x_off; // offset of this variable's data in the flat x buffer
  int g_off; // offset of this variable's gradient in the flat g buffer
  int r_off; // offset into the RMSprop accumulator buffer
  int len; // number of floats in this variable node
};

// The complete flattened graph living on GPU
struct kann_cu_graph {
  int n_node; // total number of nodes
  struct kann_cu_node *nodes; // device array [n_node] of node metadata
  struct kann_cu_node *h_nodes; // host mirror of node metadata (persistent)

  float *x; // device buffer for all node values
  float *g; // device buffer for all node gradients
  int x_total; // total floats in x buffer (at max batch)
  int g_total; // total floats in g buffer (at max batch)

  // Per-variable-node scatter info (variables are NOT contiguous in x buf)
  struct kann_cu_var_node *h_vars; // host array [n_var_nodes]
  int n_var_nodes; // number of variable nodes
  int n_var; // total number of trainable variable floats

  // Per-constant-node scatter info
  struct kann_cu_var_node *h_consts; // host array [n_const_nodes]
  int n_const_nodes; // number of constant nodes
  int n_const; // total number of constant floats

  // Training support
  float *r; // device buffer for RMSprop running average [n_var]

  // Feed node bookkeeping
  int in_node_idx; // index of the input feed node
  int truth_node_idx; // index of the truth feed node
  int cost_node_idx; // index of the cost node
  int out_node_idx; // index of the output node

  // Topological order for forward/backward
  int *h_fwd_order; // host array: internal node indices in forward order
  int *h_bwd_order; // host array: internal node indices in backward order
  int n_internal; // number of internal (non-leaf) nodes

  // Stdnorm (layer normalization, op 32) support
  float *stdnorm_si; // device buffer for per-row std_inv cache
  int stdnorm_si_total; // total floats in stdnorm_si buffer

  // RNN pre-linkage: pairs of (output_node_idx, h0_node_idx)
  // After forward, copy output's x to h0's x for recurrence
  int *h_pre_pairs; // host array [2 * n_pre_pairs]
  int n_pre_pairs; // number of pre-linked node pairs

  // Batch size tracking
  int max_batch_size; // batch size used for buffer allocation
  int cur_batch_size; // current batch size
};

// Build a flattened GPU graph from a host-side kann_t.
// The host graph should be synced to max_batch_size via kad_sync_dim first.
// Allocates device buffers sized for these dimensions.
struct kann_cu_graph *kann_cu_graph_new(kann_t *ann, int max_batch_size);

// Free all device and host memory associated with the flattened graph.
void kann_cu_graph_free(struct kann_cu_graph *g);

// Copy trainable variables from host kann_t to device graph.
void kann_cu_graph_upload_vars(struct kann_cu_graph *g, const kann_t *ann);

// Copy trainable variables from device graph back to host kann_t.
void kann_cu_graph_download_vars(struct kann_cu_graph *g, kann_t *ann);

// Copy constants from host kann_t to device graph.
void kann_cu_graph_upload_consts(struct kann_cu_graph *g, const kann_t *ann);

// Update node dimensions for current batch size using the host kann_t.
// Calls kad_sync_dim on host, updates host mirror, re-uploads metadata.
void kann_cu_sync_dim(struct kann_cu_graph *g, kann_t *ann, int batch_size);

// Run forward pass on GPU up to and including node 'to'.
void kann_cu_forward(struct kann_cu_graph *g, void *cublas_h, int to);

// Run backward pass on GPU from node 'from'.
void kann_cu_backward(struct kann_cu_graph *g, void *cublas_h, int from);

// Zero all gradients in the graph.
void kann_cu_zero_grad(struct kann_cu_graph *g);

// RMSprop update on device.
void kann_cu_rmsprop(struct kann_cu_graph *g, float lr, float decay);

// Bind input data: copy mini-batch from host to input feed node.
void kann_cu_feed_input(struct kann_cu_graph *g, int batch_size, const float *x_host);

// Bind truth data: copy mini-batch from host to truth feed node.
void kann_cu_feed_truth(struct kann_cu_graph *g, int batch_size, const float *y_host);

// Gather mini-batch from device dataset into input feed node.
// indices_d is a device array of sample indices; offset is where in
// indices_d the current mini-batch starts.
void kann_cu_gather_input(struct kann_cu_graph *g, const float *data_d, const int *indices_d,
  int offset, int batch_size, int stride);

// Gather mini-batch from device dataset into truth feed node.
void kann_cu_gather_truth(struct kann_cu_graph *g, const float *data_d, const int *indices_d,
  int offset, int batch_size, int stride);

// Feed input data directly from a device pointer (contiguous batch).
void kann_cu_feed_input_dev(
  struct kann_cu_graph *g, int batch_size, const float *x_dev, int stride);

// Read back cost scalar from device.
float kann_cu_get_cost(const struct kann_cu_graph *g);

// Read back output values from device.
void kann_cu_get_output(const struct kann_cu_graph *g, int batch_size, float *out_host);

// Apply pre-recurrence: copy output node x to h0 node x for each
// pre-linked pair. Call after forward pass for RNN temporal stepping.
void kann_cu_apply_pre(struct kann_cu_graph *g);
