/* -*- c++ -*- */

extern "C" {
#include <gkyl_gk_proj_on_basis_c2p_priv.h>
#include <gkyl_position_map.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>
#include <gkyl_velocity_map.h>
}

// Comp. to phys. mapping for phase-space projections (position map in
// configuration space, velocity map in velocity space), callable on the
// device. The context must be a GPU-resident gk_proj_on_basis_c2p_func_ctx
// whose pos_map/vel_map members are the maps' device (on_dev) objects.
GKYL_CU_D static void
gk_proj_on_basis_c2p_phase_func_cu(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = (struct gk_proj_on_basis_c2p_func_ctx *) ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p_dev(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

__global__ static void
gk_proj_on_basis_c2p_phase_func_set_cu_ker(proj_on_basis_c2p_t *fptr_d)
{
  // Assigned in device code so the function pointer is a device address.
  *fptr_d = gk_proj_on_basis_c2p_phase_func_cu;
}

extern "C" proj_on_basis_c2p_t
gk_species_projection_c2p_phase_func_cu_dev_ptr(void)
{
  proj_on_basis_c2p_t *fptr_d, fptr_ho;
  checkCuda(cudaMalloc(&fptr_d, sizeof(proj_on_basis_c2p_t)));
  gk_proj_on_basis_c2p_phase_func_set_cu_ker<<<1,1>>>(fptr_d);
  checkCuda(cudaMemcpy(&fptr_ho, fptr_d, sizeof(proj_on_basis_c2p_t), cudaMemcpyDeviceToHost));
  checkCuda(cudaFree(fptr_d));
  return fptr_ho;
}
