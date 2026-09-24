#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>

void
gkyl_ten_moment_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_ten_moment *ten_moment = container_of(base->on_dev, struct wv_ten_moment, eqn);
    gkyl_cu_free(ten_moment);
  }

  struct wv_ten_moment *ten_moment = container_of(base, struct wv_ten_moment, eqn);
  gkyl_free(ten_moment);
}

static inline void
ten_moment_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 10; i++) {
    sout[i] = 0.0;
  }
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_inew(const struct gkyl_wv_ten_moment_inp *inp)
{
  double k0 = inp->k0;
  bool use_grad_closure = inp->use_grad_closure;
  bool use_nn_closure = inp->use_nn_closure;
  int poly_order = inp->poly_order;
  struct gkyl_kann_net* ann = inp->ann;
  bool use_gpu = inp->use_gpu;

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_wv_ten_moment_cu_dev_new(k0, use_grad_closure, use_nn_closure, poly_order, ann, use_gpu);
  } 
#endif    
  struct wv_ten_moment *ten_moment = gkyl_malloc(sizeof(struct wv_ten_moment));

  ten_moment->k0 = k0;
  ten_moment->use_grad_closure = use_grad_closure;
  ten_moment->use_nn_closure = use_nn_closure;
  ten_moment->poly_order = poly_order;
  ten_moment->ann = ann;

  ten_moment->eqn.type = GKYL_EQN_TEN_MOMENT;
  ten_moment->eqn.num_equations = 10;
  ten_moment->eqn.num_waves = 5;
  ten_moment->eqn.num_diag = 10;
  
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

  ten_moment->eqn.source_func = ten_moment_source;

  ten_moment->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(ten_moment->eqn.flags);
  ten_moment->eqn.ref_count = gkyl_ref_count_init(gkyl_ten_moment_free);
  ten_moment->eqn.on_dev = &ten_moment->eqn; // CPU eqn obj points to itself

  ten_moment->eqn.embed_geo = inp->embed_geo;
  if (ten_moment->eqn.embed_geo) {
    switch (ten_moment->eqn.embed_geo->type) {
      case GKYL_EMBED_ABSORB:
        ten_moment->eqn.embed_geo->embed_func = wave_embed_absorb;
        break;

      case GKYL_EMBED_REFLECT:
        ten_moment->eqn.embed_geo->embed_func = wave_embed_reflect;
        break;

      case GKYL_EMBED_FUNC:
        break; // already set by gkyl_wv_embed_geo_new

      default:
        assert(false);
        break;
    }
  } 

  return &ten_moment->eqn;
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_new(double k0, bool use_grad_closure, bool use_nn_closure, int poly_order, struct gkyl_kann_net* ann, bool use_gpu)
{
  return gkyl_wv_ten_moment_inew( &(struct gkyl_wv_ten_moment_inp) {
      .k0 = k0,
      .use_grad_closure = use_grad_closure,
      .use_nn_closure = use_nn_closure,
      .poly_order = poly_order,
      .ann = ann,
      .use_gpu = use_gpu,
    }
  );
}

double
gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->k0;
}

bool
gkyl_wv_ten_moment_use_grad_closure(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->use_grad_closure;
}

bool
gkyl_wv_ten_moment_use_nn_closure(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->use_nn_closure;
}

int
gkyl_wv_ten_moment_poly_order(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->poly_order;
}

struct gkyl_kann_net*
gkyl_wv_ten_moment_ann(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->ann;
}
