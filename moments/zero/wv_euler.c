#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>

void
gkyl_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  
  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_euler *euler = container_of(base->on_dev, struct wv_euler, eqn);
    gkyl_cu_free(euler);
  }
  
  struct wv_euler *euler = container_of(base, struct wv_euler, eqn);
  gkyl_free(euler);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_inew(const struct gkyl_wv_euler_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if(inp->use_gpu) {
    return gkyl_wv_euler_cu_dev_inew(inp);
  } 
#endif  
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_diag = 6; // KE and PE stored separate
  
  euler->gas_gamma = inp->gas_gamma;

  switch (inp->rp_type) {
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

  euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);
  euler->eqn.on_dev = &euler->eqn; // CPU eqn obj points to itself

  euler->eqn.embed_geo = inp->embed_geo;
  if (euler->eqn.embed_geo) {
    switch (euler->eqn.embed_geo->type) {
      case GKYL_EMBED_ABSORB:
        euler->eqn.embed_geo->embed_func = wave_embed_absorb;
        break;

      case GKYL_EMBED_REFLECT:
        euler->eqn.embed_geo->embed_func = wave_embed_reflect;
        break;

      case GKYL_EMBED_FUNC:
        break; // already set by gkyl_wv_embed_geo_new

      default:
        assert(false);
        break;
    }
  } 

  return &euler->eqn;  
}

struct gkyl_wv_eqn*
gkyl_wv_euler_new(double gas_gamma, bool use_gpu)
{
  return gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = WV_EULER_RP_ROE, 
      .use_gpu = use_gpu
    }
  );
}

double
gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return euler->gas_gamma;
}
