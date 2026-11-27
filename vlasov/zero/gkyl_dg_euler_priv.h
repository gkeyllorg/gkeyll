#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_euler_kernels.h>
#include <gkyl_iso_euler_kernels.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wv_eqn.h>

// private header for use in euler DG equation object creation
// functions

// Largest surface expansion of the wall flux: 5 equations times 9 surface modes (3x tensor p2).
#define DG_EULER_MAX_FLUX_WALL 45

// Types for various kernels
typedef double (*euler_vol_t)(
  const double *w, const double *dxv, const double *jacob_pos, double param, const double *u,
  const double *p, const double *fluid, double *GKYL_RESTRICT out
);

typedef double (*euler_surf_t)(
  const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn,
  const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r,
  const double *fluid_l, const double *fluid_c, const double *fluid_r, double *GKYL_RESTRICT out
);

typedef void (*euler_wall_flux_t)(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, const int edge,
  const double *fluid_skin, double *GKYL_RESTRICT flux_wall
);

typedef double (*euler_boundary_surf_t)(
  const double *w, const double *dxv, const double *jacob_pos, const struct gkyl_wv_eqn *wv_eqn,
  const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r,
  const int edge, const double *flux_wall, const double *fluid_edge, const double *fluid_skin,
  double *GKYL_RESTRICT out
);

// for use in kernel tables
typedef struct {
  euler_vol_t kernels[4];
} gkyl_dg_euler_vol_kern_list;
typedef struct {
  euler_surf_t kernels[4];
} gkyl_dg_euler_surf_kern_list;
typedef struct {
  euler_wall_flux_t kernels[4];
} gkyl_dg_euler_wall_flux_kern_list;
typedef struct {
  euler_boundary_surf_t kernels[4];
} gkyl_dg_euler_boundary_surf_kern_list;

struct dg_euler {
  struct gkyl_dg_eqn eqn; // Base object
  euler_vol_t vol; // pointer to the volume kernel
  euler_surf_t surf[3]; // pointers to surface kernels
  euler_wall_flux_t wall_flux[3]; // pointers to the kernels of the absorbing-wall flux
  euler_boundary_surf_t boundary_surf[3]; // pointers to the absorbing-wall surface kernels
  enum gkyl_eqn_type eqn_type; // Equation type
  const struct gkyl_wv_eqn *wv_eqn; // wave equation object (Riemann solver of the interface flux)
  const struct gkyl_wave_geom *geom; // wave geometry
  double gas_gamma; // adiabatic index (Euler)
  double vt; // thermal velocity (isothermal Euler)
  double vol_param; // parameter of the volume kernel: gas_gamma (Euler) or vt (isothermal Euler)
  int cdim; // configuration-space dimensions
  struct gkyl_range conf_range; // configuration space range
  const struct gkyl_array *u; // flow velocity [ux, uy, uz]
  const struct gkyl_array *p; // pressure
  // Position map (acquired; host pointer, 0 on the device copy) and its Jacobian per direction at
  // the 1D quadrature points (borrowed; device pointer on device), constant in each cell: the
  // stride of a direction is poly_order + 1.
  struct gkyl_vlasov_position_map *pos_map;
  const struct gkyl_array *jacob_pos;
  int jacob_stride;
  int num_surf; // modes of the surface expansion of a flux component
};

// Euler: volume kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_vol_kern_list ser_vol_kernels[] = {
  {NULL, euler_vol_1x_ser_p1, euler_vol_1x_ser_p2, euler_vol_1x_ser_p3}, // 0
  {NULL, euler_vol_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_vol_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_surf_x_kernels[] = {
  {NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2, euler_surfx_1x_ser_p3}, // 0
  {NULL, euler_surfx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_surfy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_wall_flux_x_kernels[] = {
  {NULL, euler_wall_fluxx_1x_ser_p1, euler_wall_fluxx_1x_ser_p2, euler_wall_fluxx_1x_ser_p3}, // 0
  {NULL, euler_wall_fluxx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_wall_fluxx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_wall_flux_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_wall_fluxy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_wall_fluxy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_wall_flux_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_wall_fluxz_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_boundary_surf_x_kernels[] = {
  {NULL, euler_boundary_surfx_1x_ser_p1, euler_boundary_surfx_1x_ser_p2,
   euler_boundary_surfx_1x_ser_p3}, // 0
  {NULL, euler_boundary_surfx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_boundary_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_boundary_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_boundary_surfy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, euler_boundary_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_boundary_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_boundary_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Euler: volume kernel list (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_vol_kern_list ten_vol_kernels[] = {
  {NULL, euler_vol_1x_ser_p1, euler_vol_1x_ser_p2, euler_vol_1x_ser_p3}, // 0
  {NULL, euler_vol_2x_ser_p1, euler_vol_2x_tensor_p2, NULL}, // 1
  {NULL, euler_vol_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_surf_x_kernels[] = {
  {NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2, euler_surfx_1x_ser_p3}, // 0
  {NULL, euler_surfx_2x_ser_p1, euler_surfx_2x_tensor_p2, NULL}, // 1
  {NULL, euler_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_surfy_2x_ser_p1, euler_surfy_2x_tensor_p2, NULL}, // 1
  {NULL, euler_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: surf kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_wall_flux_x_kernels[] = {
  {NULL, euler_wall_fluxx_1x_ser_p1, euler_wall_fluxx_1x_ser_p2, euler_wall_fluxx_1x_ser_p3}, // 0
  {NULL, euler_wall_fluxx_2x_ser_p1, euler_wall_fluxx_2x_tensor_p2, NULL}, // 1
  {NULL, euler_wall_fluxx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_wall_flux_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_wall_fluxy_2x_ser_p1, euler_wall_fluxy_2x_tensor_p2, NULL}, // 1
  {NULL, euler_wall_fluxy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: wall flux kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_wall_flux_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_wall_fluxz_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_boundary_surf_x_kernels[] = {
  {NULL, euler_boundary_surfx_1x_ser_p1, euler_boundary_surfx_1x_ser_p2,
   euler_boundary_surfx_1x_ser_p3}, // 0
  {NULL, euler_boundary_surfx_2x_ser_p1, euler_boundary_surfx_2x_tensor_p2, NULL}, // 1
  {NULL, euler_boundary_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_boundary_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, euler_boundary_surfy_2x_ser_p1, euler_boundary_surfy_2x_tensor_p2, NULL}, // 1
  {NULL, euler_boundary_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Euler: boundary surf kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_boundary_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, euler_boundary_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: volume kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_vol_kern_list ser_iso_vol_kernels[] = {
  {NULL, iso_euler_vol_1x_ser_p1, iso_euler_vol_1x_ser_p2, iso_euler_vol_1x_ser_p3}, // 0
  {NULL, iso_euler_vol_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_vol_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_iso_surf_x_kernels[] = {
  {NULL, iso_euler_surfx_1x_ser_p1, iso_euler_surfx_1x_ser_p2, iso_euler_surfx_1x_ser_p3}, // 0
  {NULL, iso_euler_surfx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_iso_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_surfy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ser_iso_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_iso_wall_flux_x_kernels[] = {
  {NULL, iso_euler_wall_fluxx_1x_ser_p1, iso_euler_wall_fluxx_1x_ser_p2,
   iso_euler_wall_fluxx_1x_ser_p3}, // 0
  {NULL, iso_euler_wall_fluxx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_wall_fluxx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_iso_wall_flux_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_wall_fluxy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_wall_fluxy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ser_iso_wall_flux_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_wall_fluxz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: x-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_iso_boundary_surf_x_kernels[] = {
  {NULL, iso_euler_boundary_surfx_1x_ser_p1, iso_euler_boundary_surfx_1x_ser_p2,
   iso_euler_boundary_surfx_1x_ser_p3}, // 0
  {NULL, iso_euler_boundary_surfx_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_boundary_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: y-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_iso_boundary_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_boundary_surfy_2x_ser_p1, NULL, NULL}, // 1
  {NULL, iso_euler_boundary_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: z-direction (Serendipity basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ser_iso_boundary_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_boundary_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: volume kernel list (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_vol_kern_list ten_iso_vol_kernels[] = {
  {NULL, iso_euler_vol_1x_ser_p1, iso_euler_vol_1x_ser_p2, iso_euler_vol_1x_ser_p3}, // 0
  {NULL, iso_euler_vol_2x_ser_p1, iso_euler_vol_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_vol_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_iso_surf_x_kernels[] = {
  {NULL, iso_euler_surfx_1x_ser_p1, iso_euler_surfx_1x_ser_p2, iso_euler_surfx_1x_ser_p3}, // 0
  {NULL, iso_euler_surfx_2x_ser_p1, iso_euler_surfx_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_iso_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_surfy_2x_ser_p1, iso_euler_surfy_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: surf kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_surf_kern_list ten_iso_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_surfz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_iso_wall_flux_x_kernels[] = {
  {NULL, iso_euler_wall_fluxx_1x_ser_p1, iso_euler_wall_fluxx_1x_ser_p2,
   iso_euler_wall_fluxx_1x_ser_p3}, // 0
  {NULL, iso_euler_wall_fluxx_2x_ser_p1, iso_euler_wall_fluxx_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_wall_fluxx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_iso_wall_flux_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_wall_fluxy_2x_ser_p1, iso_euler_wall_fluxy_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_wall_fluxy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: wall flux kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_wall_flux_kern_list ten_iso_wall_flux_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_wall_fluxz_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: x-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_iso_boundary_surf_x_kernels[] = {
  {NULL, iso_euler_boundary_surfx_1x_ser_p1, iso_euler_boundary_surfx_1x_ser_p2,
   iso_euler_boundary_surfx_1x_ser_p3}, // 0
  {NULL, iso_euler_boundary_surfx_2x_ser_p1, iso_euler_boundary_surfx_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_boundary_surfx_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: y-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_iso_boundary_surf_y_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, iso_euler_boundary_surfy_2x_ser_p1, iso_euler_boundary_surfy_2x_tensor_p2, NULL}, // 1
  {NULL, iso_euler_boundary_surfy_3x_ser_p1, NULL, NULL} // 2
};

// Isothermal Euler: boundary surf kernel list: z-direction (Tensor basis)
GKYL_CU_D static const gkyl_dg_euler_boundary_surf_kern_list ten_iso_boundary_surf_z_kernels[] = {
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, iso_euler_boundary_surfz_3x_ser_p1, NULL, NULL} // 2
};

/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_dg_euler_free(const struct gkyl_ref_count *ref);

// Volume term. The kernels take the position-map Jacobian of the cell: its physical size in
// direction d is J_d dx_d.
GKYL_CU_DH static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);
  return euler->vol(
    xc, dx, (const double *)gkyl_array_cfetch(euler->jacob_pos, cidx), euler->vol_param,
    (const double *)gkyl_array_cfetch(euler->u, cidx),
    (const double *)gkyl_array_cfetch(euler->p, cidx), qIn, qRhsOut
  );
}

// Surface term. The wave geometry of a cell is that of its lower face. The flux of a face is the
// same for its two cells (J/J_d is the same on both sides of a face of the separable map), and
// each cell divides it by its physical size.
GKYL_CU_D static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  const struct gkyl_wave_cell_geom *geom_l = gkyl_wave_geom_get(euler->geom, idxC);
  const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(euler->geom, idxR);
  const double *jacob_pos =
    (const double *)gkyl_array_cfetch(euler->jacob_pos, gkyl_range_idx(&euler->conf_range, idxC));
  return euler->surf[dir](
    xcC, dxC, jacob_pos, euler->wv_eqn, geom_l, geom_r, qInL, qInC, qInR, qRhsOut
  );
}

// Absorbing wall: the wall flux of the boundary cell on its wall face and the Riemann flux with
// its interior neighbor on the other face. geom_l and geom_r are the geometries of the lower and
// upper faces of the boundary cell.
GKYL_CU_D static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  int idx_up[GKYL_MAX_DIM];
  for (int d = 0; d < euler->cdim; ++d) {
    idx_up[d] = idxSkin[d];
  }
  idx_up[dir] += 1;
  const struct gkyl_wave_cell_geom *geom_l = gkyl_wave_geom_get(euler->geom, idxSkin);
  const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(euler->geom, idx_up);
  const double *jacob_pos = (const double *)gkyl_array_cfetch(
    euler->jacob_pos, gkyl_range_idx(&euler->conf_range, idxSkin)
  );
  double flux_wall[DG_EULER_MAX_FLUX_WALL];
  euler->wall_flux[dir](euler->wv_eqn, edge == -1 ? geom_l : geom_r, edge, qInSkin, flux_wall);
  return euler->boundary_surf[dir](
    xcSkin, dxSkin, jacob_pos, euler->wv_eqn, geom_l, geom_r, edge, flux_wall, qInEdge, qInSkin,
    qRhsOut
  );
}

// The wall flux of one boundary cell integrated over its wall face (physical area): mode 0 of
// each block of the surface expansion times the integral of the first surface basis function,
// prod_{d != dir} J_d dx_d / sqrt(2).
GKYL_CU_DH static inline void
dg_euler_wall_flux_cell(
  const struct dg_euler *euler, int dir, int edge, const int *idx, const double *dx,
  const double *fluid_skin, double *flux_int
)
{
  int idx_face[GKYL_MAX_DIM];
  for (int d = 0; d < euler->cdim; ++d) {
    idx_face[d] = idx[d];
  }
  if (edge == 1) {
    idx_face[dir] += 1;
  }
  const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(euler->geom, idx_face);
  double flux_wall[DG_EULER_MAX_FLUX_WALL];
  euler->wall_flux[dir](euler->wv_eqn, geom, edge, fluid_skin, flux_wall);
  const double *jac =
    (const double *)gkyl_array_cfetch(euler->jacob_pos, gkyl_range_idx(&euler->conf_range, idx));
  double area = 1.0;
  for (int d = 0; d < euler->cdim; ++d) {
    if (d != dir) {
      area *= jac[d * euler->jacob_stride] * dx[d] / sqrt(2.0);
    }
  }
  for (int k = 0; k < euler->eqn.num_equations; ++k) {
    flux_int[k] = area * flux_wall[k * euler->num_surf];
  }
}
