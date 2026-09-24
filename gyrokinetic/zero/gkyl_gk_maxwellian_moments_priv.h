// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>

struct gkyl_gk_maxwellian_moments
{
  struct gkyl_basis conf_basis; // Configuration-space basis
  struct gkyl_basis phase_basis; // Phase-space basis
  int num_conf_basis; // Number of configuration-space basis functions
  int vdim_phys; // Number of physical velocity dimensions captured.
  const struct gk_geometry *gk_geom; // Geometry struct
  bool divide_jacobgeo; // Boolean for if we are dividing out the configuration-space Jacobian from density
  double mass; // Species mass
  
  struct gkyl_array *M0; 
  struct gkyl_array *M1;  
  struct gkyl_array *u_par;
  struct gkyl_array *u_par_dot_M1;  
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;
  struct gkyl_array *p_par;
  struct gkyl_array *t_par;
  struct gkyl_array *p_perp;
  struct gkyl_array *t_perp;
  struct gkyl_dg_bin_op_mem *mem;

  struct gkyl_dg_updater_moment *M0_calc; 
  struct gkyl_dg_updater_moment *M1_calc;
  struct gkyl_dg_updater_moment *M2_calc;
  struct gkyl_dg_updater_moment *M2_par_calc;
  struct gkyl_dg_updater_moment *M2_perp_calc;
};
