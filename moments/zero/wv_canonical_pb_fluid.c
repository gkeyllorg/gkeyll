#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_wv_canonical_pb_fluid_priv.h>

void
gkyl_wv_can_pb_incompress_euler_free(const struct gkyl_ref_count *ref)
{ 
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_can_pb_incompress_euler *can_pb_incompress_euler = container_of(base, struct wv_can_pb_incompress_euler, eqn);
  gkyl_free(can_pb_incompress_euler);  
}

struct gkyl_wv_eqn*
gkyl_wv_can_pb_incompress_euler_new()
{  
  struct wv_can_pb_incompress_euler *can_pb_incompress_euler = gkyl_malloc(sizeof(struct wv_can_pb_incompress_euler));

  can_pb_incompress_euler->eqn.type = GKYL_EQN_CAN_PB_INCOMPRESS_EULER;
  can_pb_incompress_euler->eqn.num_equations = 1; 

  can_pb_incompress_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(can_pb_incompress_euler->eqn.flags);
  can_pb_incompress_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_can_pb_incompress_euler_free);
  can_pb_incompress_euler->eqn.on_dev = &can_pb_incompress_euler->eqn; // CPU eqn obj points to itself

  return &can_pb_incompress_euler->eqn;
}

void
gkyl_wv_can_pb_hasegawa_mima_free(const struct gkyl_ref_count *ref)
{ 
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_can_pb_hasegawa_mima *can_pb_hasegawa_mima = container_of(base, struct wv_can_pb_hasegawa_mima, eqn);
  gkyl_free(can_pb_hasegawa_mima);  
}

struct gkyl_wv_eqn*
gkyl_wv_can_pb_hasegawa_mima_new()
{  
  struct wv_can_pb_hasegawa_mima *can_pb_hasegawa_mima = gkyl_malloc(sizeof(struct wv_can_pb_hasegawa_mima));

  can_pb_hasegawa_mima->eqn.type = GKYL_EQN_CAN_PB_HASEGAWA_MIMA;
  can_pb_hasegawa_mima->eqn.num_equations = 1; 

  can_pb_hasegawa_mima->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(can_pb_hasegawa_mima->eqn.flags);
  can_pb_hasegawa_mima->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_can_pb_hasegawa_mima_free);
  can_pb_hasegawa_mima->eqn.on_dev = &can_pb_hasegawa_mima->eqn; // CPU eqn obj points to itself

  return &can_pb_hasegawa_mima->eqn;
}

void
gkyl_wv_can_pb_hasegawa_wakatani_free(const struct gkyl_ref_count *ref)
{ 
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_can_pb_hasegawa_wakatani *can_pb_hasegawa_wakatani = container_of(base, struct wv_can_pb_hasegawa_wakatani, eqn);
  gkyl_free(can_pb_hasegawa_wakatani);  
}

struct gkyl_wv_eqn*
gkyl_wv_can_pb_hasegawa_wakatani_new(double alpha, bool is_modified)
{  
  struct wv_can_pb_hasegawa_wakatani *can_pb_hasegawa_wakatani = gkyl_malloc(sizeof(struct wv_can_pb_hasegawa_wakatani));

  can_pb_hasegawa_wakatani->eqn.type = GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI;
  can_pb_hasegawa_wakatani->eqn.num_equations = 2; 
  can_pb_hasegawa_wakatani->alpha = alpha; 
  can_pb_hasegawa_wakatani->is_modified = is_modified; 

  can_pb_hasegawa_wakatani->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(can_pb_hasegawa_wakatani->eqn.flags);
  can_pb_hasegawa_wakatani->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_can_pb_hasegawa_wakatani_free);
  can_pb_hasegawa_wakatani->eqn.on_dev = &can_pb_hasegawa_wakatani->eqn; // CPU eqn obj points to itself

  return &can_pb_hasegawa_wakatani->eqn;
}

double
gkyl_wv_can_pb_hasegawa_wakatani_alpha(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_can_pb_hasegawa_wakatani *can_pb_hasegawa_wakatani = container_of(eqn, struct wv_can_pb_hasegawa_wakatani, eqn);  
  return can_pb_hasegawa_wakatani->alpha;
}

bool
gkyl_wv_can_pb_hasegawa_wakatani_is_modified(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_can_pb_hasegawa_wakatani *can_pb_hasegawa_wakatani = container_of(eqn, struct wv_can_pb_hasegawa_wakatani, eqn);  
  return can_pb_hasegawa_wakatani->is_modified;
}
