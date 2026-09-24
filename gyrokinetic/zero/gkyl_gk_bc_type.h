#pragma once

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_bc_basic.h>
#include <assert.h>

// Boundary conditions on fields and particles in the gyrokinetic solver.
enum gkyl_gyrokinetic_bc_type {
  // Particle BCs.
  GKYL_BC_GK_SKIP = 0, // Do not apply any BCs 
  GKYL_BC_GK_SPECIES_COPY, // Copy skin into ghost.
  GKYL_BC_GK_SPECIES_REFLECT, // Reflect particles.
  GKYL_BC_GK_SPECIES_ABSORB, // Absorbing BCs.
  GKYL_BC_GK_SPECIES_FUNC, // Fill ghost cell using a user-function.
  GKYL_BC_GK_SPECIES_FIXED_FUNC, // Fixed function, time-independent.
  GKYL_BC_GK_SPECIES_ZERO_FLUX, // Zero flux.
  GKYL_BC_GK_SPECIES_SHEATH, // Sheath.
  GKYL_BC_GK_SPECIES_RECYCLE, // Recycling.
  GKYL_BC_GK_SPECIES_PERIODIC, // Periodic.
  GKYL_BC_GK_SPECIES_TWISTSHIFT, // Twist-shift.
  GKYL_BC_GK_SPECIES_BOUNDARY_VALUE, // Skin value at the boundary.
  // Field BCs.
  GKYL_BC_GK_FIELD_PERIODIC, // Periodic.
  GKYL_BC_GK_FIELD_DIRICHLET, // Dirichlet.
  GKYL_BC_GK_FIELD_NEUMANN, // Nemann.
  GKYL_BC_GK_FIELD_DIRICHLET_VARYING, // Spatially varying Dirichlet.
  GKYL_BC_GK_FIELD_BOUNDARY_VALUE, // Skin value at the boundary.
  GKYL_BC_GK_FIELD_TWISTSHIFT, // Twist-shift.
};

// Translates field BC types in gkyl_gyrokinetic_bc_type to
// the BC types used by the FEM Poisson solver.
// This translation needs to match the types in moments/zero/gkyl_fem_poisson_bctype.h
// which at the moment are:
//   enum gkyl_poisson_bc_type {
//     GKYL_POISSON_PERIODIC = 0,
//     GKYL_POISSON_DIRICHLET, // sets the value.
//     GKYL_POISSON_NEUMANN,   // sets the slope normal to the boundary.
//     GKYL_POISSON_ROBIN,  // a combination of dirichlet and neumann.  
//     GKYL_POISSON_DIRICHLET_VARYING, // sets the value, spatially varying.
//   };
static inline enum gkyl_poisson_bc_type
gkyl_gyrokinetic_translate_poisson_bc_type(enum gkyl_gyrokinetic_bc_type bc_type)
{
  enum gkyl_poisson_bc_type poisson_bc_type;
  switch (bc_type) {
    case GKYL_BC_GK_FIELD_PERIODIC:
      poisson_bc_type = GKYL_POISSON_PERIODIC;
      break;
    case GKYL_BC_GK_FIELD_DIRICHLET:
      poisson_bc_type = GKYL_POISSON_DIRICHLET;
      break;
    case GKYL_BC_GK_FIELD_NEUMANN:
      poisson_bc_type = GKYL_POISSON_NEUMANN;
      break;
    case GKYL_BC_GK_FIELD_DIRICHLET_VARYING:
      poisson_bc_type = GKYL_POISSON_DIRICHLET_VARYING;
      break;
    default:
      assert(false);
      break;
  }
  return poisson_bc_type;
}

// Translates species BC types in gkyl_gyrokinetic_bc_type to
// the BC types used by the bc_basic solver for Vlasov neutrals.
// This translation needs to match the types in vlasov/zero/gkyl_bc_basic.h
// which at the moment are:
// BC types in this updater.
//   enum gkyl_bc_basic_type { 
//     GKYL_BC_COPY = 0, 
//     GKYL_BC_ABSORB, 
//     GKYL_BC_DISTF_REFLECT, 
//     GKYL_BC_FIXED_FUNC,
//     GKYL_BC_CONF_BOUNDARY_VALUE, 
//   };
static inline enum gkyl_bc_basic_type
gkyl_gyrokinetic_translate_bc_basic_type(enum gkyl_gyrokinetic_bc_type bc_type)
{
  enum gkyl_bc_basic_type bc_basic_type;
  switch (bc_type) {
    case GKYL_BC_GK_SPECIES_COPY:
      bc_basic_type = GKYL_BC_COPY;
      break;
    case GKYL_BC_GK_SPECIES_ABSORB:
      bc_basic_type = GKYL_BC_ABSORB;
      break;
    case GKYL_BC_GK_SPECIES_REFLECT:
      bc_basic_type = GKYL_BC_DISTF_REFLECT;
      break;
    case GKYL_BC_GK_SPECIES_FIXED_FUNC:
      bc_basic_type = GKYL_BC_FIXED_FUNC;
      break;
    default:
      assert(false);
      break;
  }
  return bc_basic_type;
}
