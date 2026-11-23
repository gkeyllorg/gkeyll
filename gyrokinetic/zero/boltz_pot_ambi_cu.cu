extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_boltz_pot_ambi_priv.h>
}

// Kernel addresses must be selected on the device.
__global__ static void
boltz_pot_ambi_set_cu_kernels(struct boltz_potential *up, int poly_order)
{
  boltz_potential_set_kernels(
    up, poly_order, boltz_pot_ambi_reference_kernels[up->cdim - 1][poly_order - 1]
  );
}

void
boltz_pot_ambi_cu_dev_new(gkyl_boltz_pot_ambi *up, int poly_order)
{
  struct boltz_potential *common = &up->common;
  common->on_dev = (struct boltz_potential *)gkyl_cu_malloc(sizeof(*common));
  gkyl_cu_memcpy(common->on_dev, common, sizeof(*common), GKYL_CU_MEMCPY_H2D);
  boltz_pot_ambi_set_cu_kernels<<<1, 1>>>(common->on_dev, poly_order);
}
