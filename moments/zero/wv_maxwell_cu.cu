/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_maxwell.h>    
#include <gkyl_wv_maxwell_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to maxwell kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_maxwell_set_cu_dev_ptrs(struct wv_maxwell *maxwell)
{
  maxwell->eqn.waves_func = wave;
  maxwell->eqn.qfluct_func = qfluct;

  maxwell->eqn.flux_jump = flux_jump;
  maxwell->eqn.check_inv_func = check_inv;
  maxwell->eqn.max_speed_func = max_speed;
  maxwell->eqn.rotate_to_local_func = rot_to_local;
  maxwell->eqn.rotate_to_global_func = rot_to_global;

  maxwell->eqn.cons_to_riem = cons_to_riem;
  maxwell->eqn.riem_to_cons = riem_to_cons;

  maxwell->eqn.wall_bc_func = maxwell_wall;

  maxwell->eqn.cons_to_diag = maxwell_cons_to_diag;
}

struct gkyl_wv_eqn*
gkyl_wv_maxwell_cu_dev_new(double c, double e_fact, double b_fact)
{
  struct wv_maxwell *maxwell = (struct wv_maxwell*) gkyl_malloc(sizeof(struct wv_maxwell));

  maxwell->eqn.type = GKYL_EQN_MAXWELL;
  maxwell->eqn.num_equations = 8;  
  maxwell->eqn.num_waves = 6;
  maxwell->eqn.num_diag = 6; // Ex^2, Ey^2, Ez^2, Bx^2, By^2, Bz^2
  
  maxwell->c = c;
  maxwell->e_fact = e_fact;
  maxwell->b_fact = b_fact;

  maxwell->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_maxwell_free);

  // copy the host struct to device struct
  struct wv_maxwell *maxwell_cu = (struct wv_maxwell*) gkyl_cu_malloc(sizeof(struct wv_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct wv_maxwell), GKYL_CU_MEMCPY_H2D);

  wv_maxwell_set_cu_dev_ptrs<<<1,1>>>(maxwell_cu);

  maxwell->eqn.on_dev = &maxwell_cu->eqn; // CPU eqn obj points to itself
  return &maxwell->eqn;
}
