// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_average.h>
#include <gkyl_basis.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*canonical_pb_fluid_alpha_surf_t)(const double *w, const double *dxv, const double *phi,
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
typedef void (*canonical_pb_fluid_source_t)(const double *dxv, double alpha, 
  const double *phi, const double *n0, 
  const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT rhs); 
typedef void (*canonical_pb_fluid_subtract_zonal_t)(const double *phi_zonal, const double *n_zonal, 
  double* GKYL_RESTRICT adiabatic_coupling_phi_n); 

// for use in kernel tables
typedef struct { canonical_pb_fluid_alpha_surf_t kernels[3]; } gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list;
typedef struct { canonical_pb_fluid_source_t kernels[3]; } gkyl_dg_canonical_pb_fluid_source_kern_list;
typedef struct { canonical_pb_fluid_subtract_zonal_t kernels[3]; } gkyl_dg_canonical_pb_fluid_subtract_zonal_kern_list;

struct gkyl_dg_calc_canonical_pb_fluid_vars {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid for cell spacing and cell center
  struct gkyl_basis conf_basis; // Configuration-space basis
  int cdim; // Configuration-space dimensionality
  double alpha; // Adiabaticity parameter for adiabatic coupling of vorticity and density.
  double kappa; // Constant density gradient scale length (for turbulence drive). 
  bool is_modified; // Boolean parameter for if we are doing the modified Hasegawa-Wakatani. 

  canonical_pb_fluid_alpha_surf_t alpha_surf[3]; // kernel for computing surface expansion of configuration space flux alpha
  canonical_pb_fluid_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of configuration space flux alpha
                                                      // at upper configuration space edge

  canonical_pb_fluid_source_t canonical_pb_fluid_source; // Canonical pb fluid variables source function

  canonical_pb_fluid_subtract_zonal_t subtract_zonal; // Function for subtracting zonal component 
                                                      // of n and phi in modified Hasegawa-Wakatani  

  // Specific memory for Hasegawa-Wakatani to compute adiabatic coupling and if we are 
  // solving modified Hasegawa-Wakatani, subtract the zonal component of n and phi. 
  struct gkyl_array *n; // Density (second component of Hasegawa-Wakatani system)
  struct gkyl_array *adiabatic_coupling_phi_n; // combined (phi, n) array for adiabatic coupling

  struct gkyl_range x_local; // x range for the zonal components after integrating over y
  struct gkyl_range x_local_ext; // extended x range for the zonal components after integrating over y  
  struct gkyl_array *phi_zonal; // 1/Ly int phi dy 
  struct gkyl_array *n_zonal; // 1/Ly int n dy
  struct gkyl_array_average *int_y; // Updater for computing zonal components by averaging in y

  uint32_t flags;
  struct gkyl_dg_calc_canonical_pb_fluid_vars *on_dev; // pointer to itself or device data
};

// Canonical PB fluid alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x_ser_p1, canonical_pb_alpha_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x_ser_p1, canonical_pb_alpha_edge_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x_ser_p1, canonical_pb_alpha_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x_ser_p1, canonical_pb_alpha_edge_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Mima source update (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list ser_canonical_pb_fluid_hasegawa_mima_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Mima source update (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list tensor_canonical_pb_fluid_hasegawa_mima_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani source update (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list ser_canonical_pb_fluid_hasegawa_wakatani_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani source update (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list tensor_canonical_pb_fluid_hasegawa_wakatani_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani subtract zonal component update (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_subtract_zonal_kern_list ser_canonical_pb_fluid_subtract_zonal_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_subtract_zonal_2x_ser_p1, canonical_pb_fluid_subtract_zonal_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani subtract zonal componente update (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_subtract_zonal_kern_list tensor_canonical_pb_fluid_subtract_zonal_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_subtract_zonal_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

GKYL_CU_D
static canonical_pb_fluid_alpha_surf_t
choose_canonical_pb_fluid_alpha_surf_kern(enum gkyl_basis_type b_type, int dir, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_canonical_pb_fluid_alpha_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_fluid_alpha_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_fluid_alpha_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_fluid_alpha_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_alpha_surf_t
choose_canonical_pb_fluid_alpha_edge_surf_kern(enum gkyl_basis_type b_type, int dir, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_canonical_pb_fluid_alpha_edge_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_fluid_alpha_edge_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_fluid_alpha_edge_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_fluid_alpha_edge_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_hasegawa_mima_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_canonical_pb_fluid_hasegawa_mima_source_kernels[cdim-1].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_fluid_hasegawa_mima_source_kernels[cdim-1].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_hasegawa_wakatani_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_canonical_pb_fluid_hasegawa_wakatani_source_kernels[cdim-1].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_fluid_hasegawa_wakatani_source_kernels[cdim-1].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_default_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  return canonical_pb_fluid_default_source;
}

GKYL_CU_D
static canonical_pb_fluid_subtract_zonal_t
choose_canonical_pb_fluid_subtract_zonal_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_canonical_pb_fluid_subtract_zonal_kernels[cdim-1].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_fluid_subtract_zonal_kernels[cdim-1].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}
