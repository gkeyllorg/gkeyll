#pragma once

#include <gkyl_block_topo.h>
#include <gkyl_comm.h>
#include <gkyl_flash_decomp.h>
#include <gkyl_flash_priv.h>
#include <gkyl_flash_multib.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_rrobin_decomp.h>

// A multib_comm_conn send/recv pair used for communicating between blocks.
struct gkyl_mbcc_sr {
  struct gkyl_multib_comm_conn **recv, **send;
};

struct gkyl_flash_multib_app {
  char name[128];
  int cdim, vdim; // Dimensionality.
  struct gkyl_comm *comm; // Global communicator.
  struct gkyl_comm **block_comms; // Block communicators.

  double cfl_frac; // CFL scaling factor (default 1.0).
  int num_species; // Number of species

  int num_local_blocks; // Number of blocks handled by this app's rank.
  int *local_blocks; // List of block IDs handled by this app's rank.
  struct gkyl_flash_app **block_apps; // List of single-block apps.

  struct gkyl_flash_block_geom *block_geom; // Block geometry object.
  struct gkyl_block_topo *topo; // Block topology and connections.

  enum gkyl_flash_multib_decomp_type decomp_type; // Flag for method of decomposing the global grid into blocks and 
                                                  // assigning ranks per block.

  // Objects to handle decomposition of the grid into blocks and assignment of ranks.
  const struct gkyl_rrobin_decomp *round_robin; // Round-robin decomp.
  const struct gkyl_flash_decomp *lin_decomp; // Row-major ordered linear decomp.

  // Space for grid decomposition pointers.
  struct gkyl_rect_decomp **decomps_kin, **decomps_flu; // Blockwise grid decomposition.

  struct gkyl_mbcc_sr *mbcc_sync_kin; // Connections for phase-space (distribution function) syncing.
  struct gkyl_mbcc_sr *mbcc_sync_conf; // Connections for configuration-space syncing on the kinetic grid.
  struct gkyl_mbcc_sr *mbcc_sync_flu; // Connections for configuration-space syncing on the fluid grid.

  double tcurr;

  struct gkyl_flash_stat stat; // Simulation statistics.

  bool skip_fluid_species; // Flag to not initialize fluid species.
};

/**
 * Advance the hybrid system on all blocks by a single timestep using the
 * SSP-RK3 method. The timestep is the same on every block and rank, is never
 * larger than the suggested timestep dt0, and is returned in the status struct.
 *
 * @param app Multiblock app object.
 * @param dt0 Suggested timestep.
 * @return struct gkyl_update_status
 */
struct gkyl_update_status
gkyl_flash_multib_update_ssp_rk3(struct gkyl_flash_multib_app *app, double dt0);