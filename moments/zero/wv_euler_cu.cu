/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler.h>    
#include <gkyl_wv_euler_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to euler kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_euler_set_cu_dev_ptrs(enum gkyl_wv_euler_rp rp_type, struct wv_euler *euler)
{
  switch (rp_type) {
    case WV_EULER_RP_ROE:
      euler->eqn.num_waves = 3;  
      euler->eqn.waves_func = wave_roe_l;
      euler->eqn.qfluct_func = qfluct_roe_l;
      break;

    case WV_EULER_RP_HLLC:
      euler->eqn.num_waves = 3;  
      euler->eqn.waves_func = wave_hllc_l;
      euler->eqn.qfluct_func = qfluct_hllc_l;
      break;
      
    case WV_EULER_RP_LAX:
      euler->eqn.num_waves = 2;
      euler->eqn.waves_func = wave_lax_l;
      euler->eqn.qfluct_func = qfluct_lax_l;
      break;   

    case WV_EULER_RP_HLL:
      euler->eqn.num_waves = 2;  
      euler->eqn.waves_func = wave_hll_l;
      euler->eqn.qfluct_func = qfluct_hll_l;
      break;
  }

  euler->eqn.flux_jump = flux_jump;
  euler->eqn.check_inv_func = check_inv;
  euler->eqn.max_speed_func = max_speed;
  euler->eqn.rotate_to_local_func = rot_to_local;
  euler->eqn.rotate_to_global_func = rot_to_global;

  euler->eqn.wall_bc_func = euler_wall;
  euler->eqn.no_slip_bc_func = euler_no_slip;

  euler->eqn.cons_to_riem = cons_to_riem;
  euler->eqn.riem_to_cons = riem_to_cons;

  euler->eqn.cons_to_diag = euler_cons_to_diag;

  euler->eqn.source_func = euler_source;
}

struct gkyl_wv_eqn*
gkyl_wv_euler_cu_dev_inew(const struct gkyl_wv_euler_inp *inp)
{
  struct wv_euler *euler = (struct wv_euler*) gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_diag = 6; // KE and PE stored separate
  
  euler->gas_gamma = inp->gas_gamma;

  euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);

  // copy the host struct to device struct
  struct wv_euler *euler_cu = (struct wv_euler*) gkyl_cu_malloc(sizeof(struct wv_euler));
  gkyl_cu_memcpy(euler_cu, euler, sizeof(struct wv_euler), GKYL_CU_MEMCPY_H2D);

  wv_euler_set_cu_dev_ptrs<<<1,1>>>(inp->rp_type, euler_cu);

  euler->eqn.on_dev = &euler_cu->eqn; // CPU eqn obj points to itself
  return &euler->eqn;
}
