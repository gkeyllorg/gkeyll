/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_ten_moment.h>    
#include <gkyl_wv_ten_moment_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to ten moment kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_ten_moment_set_cu_dev_ptrs(struct wv_ten_moment *ten_moment)
{
  ten_moment->eqn.waves_func = wave;
  ten_moment->eqn.qfluct_func = qfluct;

  ten_moment->eqn.check_inv_func = check_inv;
  ten_moment->eqn.max_speed_func = max_speed;
  ten_moment->eqn.rotate_to_local_func = rot_to_local;
  ten_moment->eqn.rotate_to_global_func = rot_to_global;

  ten_moment->eqn.cons_to_riem = cons_to_riem;
  ten_moment->eqn.riem_to_cons = riem_to_cons;

  ten_moment->eqn.wall_bc_func = ten_moment_wall;

  ten_moment->eqn.cons_to_diag = gkyl_default_cons_to_diag;
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_cu_dev_inew(const struct gkyl_wv_ten_moment_inp *inp)
{
  double k0 = inp->k0;
  bool use_grad_closure = inp->use_grad_closure;
  bool use_nn_closure = inp->use_nn_closure;
  int poly_order = inp->poly_order;
  struct gkyl_kann_net* ann = inp->ann;
  bool use_gpu = inp->use_gpu;

  struct wv_ten_moment *ten_moment = (struct wv_ten_moment*) gkyl_malloc(sizeof(struct wv_ten_moment));

  ten_moment->k0 = k0;
  ten_moment->use_grad_closure = use_grad_closure;
  ten_moment->use_nn_closure = use_nn_closure;
  ten_moment->poly_order = poly_order;
  ten_moment->ann = ann;

  ten_moment->eqn.type = GKYL_EQN_TEN_MOMENT;
  ten_moment->eqn.num_equations = 10;
  ten_moment->eqn.num_waves = 5;
  ten_moment->eqn.num_diag = 10;

  ten_moment->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(ten_moment->eqn.flags);
  ten_moment->eqn.ref_count = gkyl_ref_count_init(gkyl_ten_moment_free);

  // copy the host struct to device struct
  struct wv_ten_moment *ten_moment_cu = (struct wv_ten_moment*) gkyl_cu_malloc(sizeof(struct wv_ten_moment));
  gkyl_cu_memcpy(ten_moment_cu, ten_moment, sizeof(struct wv_ten_moment), GKYL_CU_MEMCPY_H2D);

  wv_ten_moment_set_cu_dev_ptrs<<<1,1>>>(ten_moment_cu);

  ten_moment->eqn.on_dev = &ten_moment_cu->eqn; // CPU eqn obj points to itself
  return &ten_moment->eqn;
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_cu_dev_new(double k0, bool use_grad_closure, bool use_nn_closure, int poly_order, struct gkyl_kann_net* ann, bool use_gpu)
{
  struct gkyl_wv_ten_moment_inp ten_moment_inp = {
    .k0 = k0,
    .use_grad_closure = use_grad_closure,
    .use_nn_closure = use_nn_closure,
    .poly_order = poly_order,
    .ann = ann,
    .use_gpu = use_gpu,
  };

  return gkyl_wv_ten_moment_cu_dev_inew(&ten_moment_inp);
}