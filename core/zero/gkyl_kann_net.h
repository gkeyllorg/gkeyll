#pragma once

#include <gkyl_knutils.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <kann.h>

#include <stdbool.h>

/**
 * Gkeyll wrapper around a KANN neural network. Provides a unified
 * interface for host and (future) GPU execution. On host, delegates
 * to the underlying KANN library. On GPU, uses CUDA kernels with
 * cuBLAS for matrix operations.
 */
struct gkyl_kann_net;

/**
 * Hyperparameters for training a feed-forward network.
 */
struct gkyl_kann_train_params {
  float learning_rate;    // learning rate for RMSprop
  int mini_size;          // mini-batch size
  int max_epoch;          // maximum number of training epochs
  int max_drop_streak;    // max epochs without validation improvement
  float frac_val;         // fraction of samples used for validation
};

/**
 * Create a new KANN network wrapper from a compiled computational
 * graph. The cost node must be a scalar (n_d == 0). Takes ownership
 * of the computational graph.
 *
 * @param cost Cost node (root of the computational graph)
 * @param use_gpu If true, allocate network data on GPU
 * @return New network wrapper, or NULL on failure
 */
struct gkyl_kann_net* gkyl_kann_net_new(kad_node_t *cost, bool use_gpu);

/**
 * Load a previously saved KANN network from file.
 *
 * @param filename Path to the saved network file
 * @param use_gpu If true, allocate network data on GPU
 * @return Loaded network wrapper, or NULL on failure
 */
struct gkyl_kann_net* gkyl_kann_net_load(const char *filename, bool use_gpu);

/**
 * Save the network to file. If the network is on GPU, weights are
 * copied to host before saving.
 *
 * @param net Network to save
 * @param filename Path to save to
 */
void gkyl_kann_net_save(const struct gkyl_kann_net *net, const char *filename);

/**
 * Train a feed-forward network with a single input and single output
 * feed node using mini-batch RMSprop.
 *
 * @param net Network to train
 * @param params Training hyperparameters
 * @param inp Training input vectors (nvec samples, each of length N)
 * @param out Training output vectors (nvec samples, each of length N)
 * @return Number of epochs completed
 */
int gkyl_kann_net_train_fnn1(struct gkyl_kann_net *net,
  const struct gkyl_kann_train_params *params,
  const struct gkyl_kn_vec *inp, const struct gkyl_kn_vec *out);

/**
 * Run inference on a batch of input vectors.
 *
 * @param net Network to evaluate
 * @param inp Input vectors (nvec samples, each of length N = n_in)
 * @param out Output vectors (nvec samples, each of length N = n_out).
 *            Filled by this function.
 */
void gkyl_kann_net_apply(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out);

/**
 * Run sequential RNN inference over a sequence of input vectors.
 * Processes one timestep at a time, applying the recurrent pre-linkage
 * between steps so the hidden state carries forward. On CPU this uses
 * kann_rnn_start/kann_apply1/kann_rnn_end; on GPU it uses
 * kann_cu_apply_pre between forward passes.
 *
 * @param net Network to evaluate (must have RNN pre-linkage)
 * @param inp Input vectors (nvec timesteps, each of length N = n_in)
 * @param out Output vectors (nvec timesteps, each of length N = n_out).
 *            Filled by this function.
 */
void gkyl_kann_net_apply_rnn(struct gkyl_kann_net *net,
  const struct gkyl_kn_vec *inp, struct gkyl_kn_vec *out);

/**
 * Get the input dimension of the network (number of floats per input
 * sample, excluding the batch dimension).
 *
 * @param net Network to query
 * @return Input dimension, or -1 if not found
 */
int gkyl_kann_net_dim_in(const struct gkyl_kann_net *net);

/**
 * Get the output dimension of the network.
 *
 * @param net Network to query
 * @return Output dimension, or -1 if not found
 */
int gkyl_kann_net_dim_out(const struct gkyl_kann_net *net);

/**
 * Set the verbose level for training output. Level 3 prints per-epoch
 * training and validation cost.
 *
 * @param level Verbose level (0 = silent, 3 = verbose)
 */
void gkyl_kann_net_set_verbose(int level);

/**
 * Check if the network data lives on GPU.
 *
 * @param net Network to check
 * @return true if on GPU, false otherwise
 */
bool gkyl_kann_net_is_cu_dev(const struct gkyl_kann_net *net);

/**
 * Acquire a reference to the network. Must be released with
 * gkyl_kann_net_release.
 *
 * @param net Network to acquire
 * @return Pointer to acquired network
 */
struct gkyl_kann_net* gkyl_kann_net_acquire(const struct gkyl_kann_net *net);

/**
 * Release a reference to the network. When the reference count
 * reaches zero, the network and all associated memory is freed.
 *
 * @param net Network to release
 */
void gkyl_kann_net_release(struct gkyl_kann_net *net);
