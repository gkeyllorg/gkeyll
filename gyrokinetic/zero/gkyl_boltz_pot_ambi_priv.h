#pragma once
#include <gkyl_boltz_pot_ambi.h>
#include <gkyl_boltz_potential_priv.h>

struct gkyl_boltz_pot_ambi {
  struct boltz_potential common; // Shared boundary and volume machinery.
};

GKYL_CU_TABLE static const boltz_reference_t boltz_pot_ambi_reference_kernels[3][2] = {
  {boltz_pot_ambi_reference_1x_ser_p1, boltz_pot_ambi_reference_1x_ser_p2},
  {boltz_pot_ambi_reference_2x_ser_p1, boltz_pot_ambi_reference_2x_ser_p2},
  {boltz_pot_ambi_reference_3x_ser_p1, boltz_pot_ambi_reference_3x_ser_p2}
};

#ifdef GKYL_HAVE_CUDA
/**
 * Copy common state to the device and install device kernel pointers.
 *
 * @param up Initialized host updater.
 * @param poly_order Configuration polynomial order.
 */
void boltz_pot_ambi_cu_dev_new(gkyl_boltz_pot_ambi *up, int poly_order);
#endif
