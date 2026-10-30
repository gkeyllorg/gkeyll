#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_geom_preset_priv.h>

evalf_t
gkyl_dg_gr_maxwell_preset_lapse(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_lapse_kern(preset_geom_type);
}

evalf_t
gkyl_dg_gr_maxwell_preset_shift(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_shift_kern(preset_geom_type);
}

evalf_t
gkyl_dg_gr_maxwell_preset_geom_factor_con(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_geom_factor_con_kern(preset_geom_type);
}

evalf_t
gkyl_dg_gr_maxwell_preset_h_ij(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_h_ij_kern(preset_geom_type);
}

evalf_t
gkyl_dg_gr_maxwell_preset_h_ij_inv(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_h_ij_inv_kern(preset_geom_type);
}

evalf_t
gkyl_dg_gr_maxwell_preset_det_h(enum gkyl_triad_preset_geom_type preset_geom_type)
{
  return choose_det_h_kern(preset_geom_type);
}
