#pragma once

#include <gkyl_vlasov_position_map.h>

// Device-side launchers for the conf-Jacobian divide/rescale operations
// (implemented in vlasov_position_map_cu.cu). The conf Jacobian is a per-cell
// constant, so these perform an exact scalar divide/multiply on every
// phase-space coefficient.
void gkyl_vlasov_position_map_divide_jacobpos_cu(int cdim,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *jacob_pos_gauss,
  const struct gkyl_array *Jf, struct gkyl_array *f_no_J);

void gkyl_vlasov_position_map_rescale_jacobpos_cu(int cdim,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *jacob_pos_gauss,
  const struct gkyl_array *f_no_J, struct gkyl_array *Jf);
