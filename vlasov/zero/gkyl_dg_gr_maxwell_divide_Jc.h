#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_basis.h>

/**
 * Divide out the configuration-space Jacobian from J*Fields to obtain Fields for use in the DG update. 
 * This update is performed at Gauss-Legendre quadrature points because, by construction,
 * we know Jc at these points and our Jc is static in time, so we can obtain the weakly
 * equivalent Fields from JcFields in our DG expansion (since we can always transform between the
 * Gauss-Legendre nodal basis and our modal basis). 
 * 
 * @param conf_basis Configuration-space basis. 
 * @param conf_range Configuration-space range. 
 * @param det_h Configuration-space Jacobian at Gauss-Legendre quadrature points. 
 * @param field_with_J_con Input array JcFields (contravaraint components) with Jacobian factor. 
 * @param field_no_J_con Output array JcFields (contravaraint components) without Jacobian factor. 
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_gr_maxwell_divide_Jc(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_with_J_con,
  struct gkyl_array *field_no_J_con, bool use_gpu
);

/**
 * Divide out the configuration-space Jacobian from JFields to obtain Fields on
 * NV-GPU. See new() method for documentation.
 */
void gkyl_dg_gr_maxwell_divide_Jc_cu(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_with_J_con,
  struct gkyl_array *field_no_J_con
);

/**
 * Multiply the configuration-space Jacobian by Fields to obtain JFields use in the DG update. 
 * This update is performed at Gauss-Legendre quadrature points because, by construction,
 * we know Jc at these points and our Jc is static in time, so we can obtain the weakly
 * equivalent Fields from JcFields in our DG expansion (since we can always transform between the
 * Gauss-Legendre nodal basis and our modal basis). 
 * 
 * @param conf_basis Configuration-space basis. 
 * @param conf_range Configuration-space range. 
 * @param det_h Configuration-space Jacobian at Gauss-Legendre quadrature points. 
 * @param field_no_J_con Input array JcFields (contravaraint components) without Jacobian factor.
 * @param field_with_J_con Output array JcFields (contravaraint components) with Jacobian factor.  
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_gr_maxwell_rescale_Jc(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_no_J_con,
  struct gkyl_array *field_with_J_con, bool use_gpu
);

/**
 * Multiply the configuration-space Jacobian by f to obtain Jf on
 * NV-GPU. See new() method for documentation.
 */
void gkyl_dg_gr_maxwell_rescale_Jc_cu(
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range,
  const struct gkyl_surf_and_vol_node_arrays *det_h, const struct gkyl_array *field_no_J_con,
  struct gkyl_array *field_with_J_con
);
