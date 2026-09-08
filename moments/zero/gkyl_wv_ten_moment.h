#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_kann_net.h>
#include <gkyl_knutils.h>

// Input packaged as a struct.
struct gkyl_wv_ten_moment_inp {
  double k0; // Closure parameter.
  bool use_grad_closure; // Should we use gradient-based closure?
  bool use_nn_closure; // Should we use neural network-based closure?
  int poly_order; // Polynomial order of learned DG coefficients.
  struct gkyl_kann_net *ann; // Neural network architecture.
  struct gkyl_wv_embed_geo *embed_geo; // embedded geometry
  bool use_gpu; // Boolean to determine whether wave equation object is on host or device.
};

/**
 * Create a new 10-moment equation object.
 *
 * @param k0 Closure parameter.
 * @param use_grad_closure Should we use gradient-based closure?
 * @param use_nn_closure Should we use neural network-based closure?
 * @param poly_order Polynomial order of learned DG coefficients.
 * @param ann Neural network architecture.
 * @param use_gpu Boolean to determine whether wave equation object is on host or device.
 * @return Pointer to 10-moment equation object.
 */
struct gkyl_wv_eqn *gkyl_wv_ten_moment_new(double k0, bool use_grad_closure, bool use_nn_closure,
  int poly_order, struct gkyl_kann_net *ann, bool use_gpu);

/**
 * Create a new 10-moment equation object.
 *
 * @param inp Input parameters.
 * @return Pointer to 10-moment equation object.
 */
struct gkyl_wv_eqn *gkyl_wv_ten_moment_inew(const struct gkyl_wv_ten_moment_inp *inp);

/**
 * Create a new 10-moment equation object that lives on NV-GPU.
 *
 * @param k0 Closure parameter.
 * @param use_grad_closure Should we use gradient-based closure?
 * @param use_nn_closure Should we use neural network-based closure?
 * @param poly_order Polynomial order of learned DG coefficients.
 * @param ann Neural network architecture.
 * @param use_gpu Boolean to determine whether wave equation object is on host or device.
 * @return Pointer to 10-moment equation object on device.
 */
struct gkyl_wv_eqn *gkyl_wv_ten_moment_cu_dev_new(double k0, bool use_grad_closure,
  bool use_nn_closure, int poly_order, struct gkyl_kann_net *ann, bool use_gpu);

/**
 * Create a new 10-moment equation object that lives on NV-GPU.
 *
 * @param inp Input parameters.
 * @return Pointer to 10-moment equation object on device.
 */
struct gkyl_wv_eqn *gkyl_wv_ten_moment_cu_dev_inew(const struct gkyl_wv_ten_moment_inp *inp);

/**
 * Get closure parameter.
 * 
 * @param wv 10-moment equation object.
 * @return Closure parameter.
 */
double gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn *wv);

/**
 * Should we use gradient-based closure?
 * 
 * @param wv 10-moment equation object.
 * @return True if using gradient-based closure, false otherwise.
 */
bool gkyl_wv_ten_moment_use_grad_closure(const struct gkyl_wv_eqn *wv);

/**
 * Should we use neural netowork-based closure?
 * 
 * @param wv 10-moment equation object.
 * @return True if using neural network-based closure, false otherwise.
 */
bool gkyl_wv_ten_moment_use_nn_closure(const struct gkyl_wv_eqn *wv);

/**
 * Get polynomial order of learned DG coefficients.
 * 
 * @param wv 10-moment equation object.
 * @return Polynomial order of learned DG coefficients.
 */
int gkyl_wv_ten_moment_poly_order(const struct gkyl_wv_eqn *wv);

/**
 * Get neural network architecture.
 * 
 * @param wv 10-moment equation object
 * @return Neural network architecture.
 */
struct gkyl_kann_net *gkyl_wv_ten_moment_ann(const struct gkyl_wv_eqn *wv);
