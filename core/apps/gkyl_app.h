#pragma once

#include <stdbool.h>
#include <gkyl_array_rio.h>

// Update status
struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

// Status of restart
struct gkyl_app_restart_status {
  enum gkyl_array_rio_status io_status; // status of the file read
  int frame; // frame number of file read
  double stime; // simulation time at which data was read
};

// inputs from user for specifying range and communicator to use
struct gkyl_app_comm_low_inp {
  struct gkyl_range local_range; // local range over which App operates
  struct gkyl_comm *comm; // communicator to use
};

// BC for blocks
struct gkyl_block_physical_bcs {
  int bidx; // block index
  int dir; // direction in which BC is specified
  enum gkyl_edge_loc edge; // which edge this BC is for
  int bc_type; // BC code
};

// Parallelization-related inputs.
struct gkyl_app_parallelism_inp {
  bool use_gpu; // Run on the GPU(s).
  int cuts[3]; // Number of subdomain in each dimension.
  struct gkyl_comm *comm; // Communicator to use.
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_SKIP, // Do not apply any BCs to field
  GKYL_SPECIES_REFLECT, // perfect reflector
  GKYL_SPECIES_ABSORB, // Absorbing BCs
  GKYL_SPECIES_NO_SLIP, // no-slip boundary conditions
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_SPECIES_FUNC, // Function boundary conditions
  GKYL_SPECIES_FIXED_FUNC, // Fixed function, time-independent, boundary conditions
  GKYL_SPECIES_EMISSION, // Emission spectrum BCs
  GKYL_SPECIES_ZERO_FLUX, // Zero flux BCs; must be applied on both lower and upper BC
  GKYL_SPECIES_RECYCLE // Recycling BCs
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_SKIP, // Do not apply any BCs to field
  GKYL_FIELD_PEC_WALL, // Maxwell's perfect electrical conductor (zero normal B and zero tangent E)
  GKYL_FIELD_SYM_WALL, // Maxwell's symmetry BC (zero normal E and zero tangent B)
  GKYL_FIELD_RESERVOIR, // Reservoir Maxwell's BCs for heat flux problem
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_FIELD_FUNC, // Function boundary conditions
  GKYL_FIELD_DIRICHLET, // Dirichlet boundary conditions
  GKYL_FIELD_NEUMANN, // Nemann boundary conditions
  GKYL_FIELD_NONE // Do not apply any boundary conditions
};

// Type of file import for initial conditions
enum gkyl_ic_import_type {
  GKYL_IC_IMPORT_NONE = 0,
  GKYL_IC_IMPORT_F, // Import f only.
  GKYL_IC_IMPORT_AF // Import f and scale by alpha(x).
};
