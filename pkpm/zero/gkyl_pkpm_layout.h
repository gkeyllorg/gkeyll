#pragma once

// Layouts of the PKPM arrays shared by the updaters.

// Blocks of the evolved normalized distributions, one phase-space expansion each:
// [F, G, M0_perp (x, y, z), N_perp (x, y, z)]. M0_perp and N_perp are the first gyroangle
// harmonic as Cartesian vectors perpendicular to b.
enum gkyl_pkpm_dist_block {
  PKPM_DIST_F = 0,
  PKPM_DIST_G = 1,
  PKPM_DIST_M0 = 2, // first of the three Cartesian components of M0_perp
  PKPM_DIST_N = 5, // first of the three Cartesian components of N_perp
  PKPM_DIST_NUM_BLOCKS = 8
};

// Components of the fluid state, one configuration-space expansion each:
// [rho, rho u (x, y, z), E].
enum gkyl_pkpm_fluid_comp {
  PKPM_FLUID_RHO = 0,
  PKPM_FLUID_MOM = 1, // first of the three Cartesian components of the momentum density
  PKPM_FLUID_ENERGY = 4,
  PKPM_FLUID_NUM_COMP = 5
};

// Primitive moments at the configuration-space nodes, one value per node each:
// [n, u (x, y, z), vtsq].
enum gkyl_pkpm_prim_comp {
  PKPM_PRIM_N = 0,
  PKPM_PRIM_U = 1, // first of the three Cartesian components of the flow velocity
  PKPM_PRIM_VTSQ = 4,
  PKPM_PRIM_NUM_COMP = 5
};

// Velocity modes per configuration-space mode of the hybrid basis (p = 2 in w): the phase-space
// basis has PKPM_NUM_VEL_MODES times the modes of the configuration-space basis.
#define PKPM_NUM_VEL_MODES 3

// Largest number of modes, and of nodes, of the p = 1 configuration-space basis (3x).
#define PKPM_MAX_CONF_BASIS 8
