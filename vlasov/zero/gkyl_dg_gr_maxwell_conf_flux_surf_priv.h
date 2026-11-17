// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <assert.h>

// Largest number of configuration-space basis functions (3x tensor p2).
#define GR_MAXWELL_MAX_CONF_BASIS 27

typedef void (*dg_gr_maxwell_alpha_quad_conf_t)(
  const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, const int geom_edge,
  const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse, const double *shift,
  const double *h_ij, const double *h_ij_inv, const double *det_h, const double *field_con_l,
  const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
  double *GKYL_RESTRICT flux_l, double *GKYL_RESTRICT flux_r, double *GKYL_RESTRICT alpha_quad
);

typedef double (*lax_flux_t)(
  const double *dxv, const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *det_h, const double *flux_l, const double *flux_r, const double *alpha_quad,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
);

typedef double (*conf_flux_surf_t)(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
);

// for use in kernel tables
typedef struct {
  dg_gr_maxwell_alpha_quad_conf_t kernels[4];
} gkyl_dg_gr_maxwell_alpha_quad_kern_list;
typedef struct {
  lax_flux_t kernels[4];
} gkyl_lax_flux_kern_list;
typedef struct {
  conf_flux_surf_t kernels[4];
} gkyl_conf_flux_surf_kern_list;

struct gkyl_dg_gr_maxwell_conf_flux_surf {
  struct gkyl_rect_grid conf_grid; // Conf-space grid.
  gkyl_dg_gr_maxwell_inp gr_maxwell_data; // Parameters needed by GR-Maxwell kernels.
  int cdim; // Configuration-space dimensions.
  dg_gr_maxwell_alpha_quad_conf_t
    dg_gr_maxwell_alpha_quad[3]; // Contribution to alpha_c at quadrature points.
  lax_flux_t lax_flux[3]; // Convert nodal Lax-Friedrichs flux to modal surface expansion.
  conf_flux_surf_t
    conf_flux_surf; // Assembly function for computing modal surface expansion of configuration-space fluxes.
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space mapping object (acquired host-side for lifetime safety).
  const struct gkyl_array
    *jacob_pos; // Position-map Jacobian (borrowed from pos_map; extended conf range).
  int theta_pole_lo
    [GKYL_MAX_CDIM]; // (Lower side) Default zeros, but 1 if any directions use theta-pole BC's.
  int theta_pole_up
    [GKYL_MAX_CDIM]; // (Upper side) Default zeros, but 1 if any directions use theta-pole BC's.
  int num_conf_basis; // Number of configuration-space basis functions.
  // Sign of each configuration-space basis function under a reflection in each direction.
  double mirror_sign[GKYL_MAX_CDIM][GR_MAXWELL_MAX_CONF_BASIS];

  uint32_t flags;
  bool use_gpu;
  struct gkyl_dg_gr_maxwell_conf_flux_surf *on_dev; // pointer to itself or device data.
};

// Set the reflection signs of the configuration-space basis functions (host only: uses the
// function pointers of the basis).
static inline void
conf_flux_surf_set_mirror_sign(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, const struct gkyl_basis *conf_basis
)
{
  assert(conf_basis->num_basis <= GR_MAXWELL_MAX_CONF_BASIS);
  up->num_conf_basis = conf_basis->num_basis;
  double ones[GR_MAXWELL_MAX_CONF_BASIS];
  for (int k = 0; k < conf_basis->num_basis; ++k) {
    ones[k] = 1.0;
  }
  for (int d = 0; d < conf_basis->ndim; ++d) {
    conf_basis->flip_odd_sign(d, ones, up->mirror_sign[d]);
  }
}

// State beyond a theta-pole face: the state of the cell next to the pole reflected across the
// face. E^r, E^phi, B^r, B^phi and the correction potentials are even, E^theta and B^theta odd
// (r-theta-phi ordering of the components, as in the theta-pole ghost-cell boundary condition).
GKYL_CU_DH static inline void
conf_flux_surf_theta_pole_mirror(
  const struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *field,
  double *GKYL_RESTRICT field_mirror
)
{
  const double comp_parity[8] = {1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0};
  int num_basis = up->num_conf_basis;
  for (int c = 0; c < 8; ++c) {
    for (int k = 0; k < num_basis; ++k) {
      field_mirror[c * num_basis + k] =
        comp_parity[c] * up->mirror_sign[dir][k] * field[c * num_basis + k];
    }
  }
}

GKYL_CU_DH static void inline conf_flux_surf_alpha_quad(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT flux_l, double *GKYL_RESTRICT flux_r,
  double *GKYL_RESTRICT alpha_quad
)
{
  up->dg_gr_maxwell_alpha_quad[dir](
    &up->gr_maxwell_data, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv,
    det_h, field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );
}

// Configuration-space flux computation.
// Separated like this to correctly allocate stack-side arrays for quadrature evaluation.
GKYL_CU_DH static double
conf_flux_surf_1x_p1(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 1 quadrature node on a 1x surface, 8 field components.
  double alpha_quad[1] = {0.0};
  double flux_l[8] = {0.0};
  double flux_r[8] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

GKYL_CU_DH static double
conf_flux_surf_1x_p2(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 1 quadrature node on a 1x surface, 8 field components.
  double alpha_quad[1] = {0.0};
  double flux_l[8] = {0.0};
  double flux_r[8] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

GKYL_CU_DH static double
conf_flux_surf_2x_p1(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 2 quadrature nodes on a 2x surface, 8 field components.
  double alpha_quad[2] = {0.0};
  double flux_l[16] = {0.0};
  double flux_r[16] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

GKYL_CU_DH static double
conf_flux_surf_2x_p2(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 3 quadrature nodes on a 2x p2 surface, 8 field components.
  double alpha_quad[3] = {0.0};
  double flux_l[24] = {0.0};
  double flux_r[24] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

GKYL_CU_DH static double
conf_flux_surf_3x_p1(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 4 quadrature nodes on a 3x surface, 8 field components.
  double alpha_quad[4] = {0.0};
  double flux_l[32] = {0.0};
  double flux_r[32] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

GKYL_CU_DH static double
conf_flux_surf_3x_p2(
  struct gkyl_dg_gr_maxwell_conf_flux_surf *up, int dir, const double *w, const double *dxv,
  const int geom_edge, const double *jacob_pos_l, const double *jacob_pos_r, const double *lapse,
  const double *shift, const double *h_ij, const double *h_ij_inv, const double *det_h,
  const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l,
  const double *field_no_J_con_r, double *GKYL_RESTRICT conf_flux_surf
)
{
  // 4 quadrature nodes on a 3x surface, 8 field components.
  double alpha_quad[9] = {0.0};
  double flux_l[72] = {0.0};
  double flux_r[72] = {0.0};

  // Compute the maximum eigenvalue at each quadrature point and retrun the fluxes.
  conf_flux_surf_alpha_quad(
    up, dir, w, dxv, geom_edge, jacob_pos_l, jacob_pos_r, lapse, shift, h_ij, h_ij_inv, det_h,
    field_con_l, field_con_r, field_no_J_con_l, field_no_J_con_r, flux_l, flux_r, alpha_quad
  );

  // Compute nodal Lax-Friedrichs flux
  double cflrate = up->lax_flux[dir](
    dxv, geom_edge, jacob_pos_l, jacob_pos_r, det_h, flux_l, flux_r, alpha_quad, field_con_l,
    field_con_r, field_no_J_con_l, field_no_J_con_r, conf_flux_surf
  );

  // Return the cfl rate
  return cflrate;
}

// Configuration-space flux computation assembly function.
GKYL_CU_D static const gkyl_conf_flux_surf_kern_list conf_flux_surf_kernels[] = {
  // 1x kernels
  {NULL, conf_flux_surf_1x_p1, conf_flux_surf_1x_p2, NULL}, // 0
  // 2x kernels
  {NULL, conf_flux_surf_2x_p1, conf_flux_surf_2x_p2, NULL}, // 1
  // 3x kernels
  {NULL, conf_flux_surf_3x_p1, conf_flux_surf_3x_p2, NULL} // 2
};

// Nodal Lax-Friedrichs flux (Serendipity basis).
GKYL_CU_D static const gkyl_lax_flux_kern_list ser_lax_flux_x_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_x_1x_ser_p1, lax_flux_x_1x_ser_p2, NULL}, // 0
  // 2x kernels
  {NULL, lax_flux_x_2x_ser_p1, lax_flux_x_2x_ser_p2, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_x_3x_ser_p1, lax_flux_x_3x_ser_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_lax_flux_kern_list ser_lax_flux_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  // 2x kernels
  {NULL, lax_flux_y_2x_ser_p1, lax_flux_y_2x_ser_p2, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_y_3x_ser_p1, lax_flux_y_3x_ser_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_lax_flux_kern_list ser_lax_flux_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_z_3x_ser_p1, lax_flux_z_3x_ser_p2, NULL} // 2
};

// Nodal Lax-Friedrichs flux (Tensor basis).
GKYL_CU_D static const gkyl_lax_flux_kern_list ten_lax_flux_x_kernels[] = {
  // 1x kernels
  {NULL, lax_flux_x_1x_tensor_p1, lax_flux_x_1x_tensor_p2, NULL}, // 0
  // 2x kernels
  {NULL, lax_flux_x_2x_tensor_p1, lax_flux_x_2x_tensor_p2, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_x_3x_tensor_p1, lax_flux_x_3x_tensor_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_lax_flux_kern_list ten_lax_flux_y_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  // 2x kernels
  {NULL, lax_flux_y_2x_tensor_p1, lax_flux_y_2x_tensor_p2, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_y_3x_tensor_p1, lax_flux_y_3x_tensor_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_lax_flux_kern_list ten_lax_flux_z_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 1
  // 3x kernels
  {NULL, lax_flux_z_3x_tensor_p1, lax_flux_z_3x_tensor_p2, NULL} // 2
};

// alpha_c evaluated at quadrature points.
GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ser_dg_gr_maxwell_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_1x_ser_p1, dg_gr_maxwell_alpha_quad_x_1x_ser_p2, NULL}, // 0
    // 2x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_2x_ser_p1, dg_gr_maxwell_alpha_quad_x_2x_ser_p2, NULL}, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_3x_ser_p1, dg_gr_maxwell_alpha_quad_x_3x_ser_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ser_dg_gr_maxwell_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    // 2x kernels
    {NULL, dg_gr_maxwell_alpha_quad_y_2x_ser_p1, dg_gr_maxwell_alpha_quad_y_2x_ser_p2, NULL}, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_y_3x_ser_p1, dg_gr_maxwell_alpha_quad_y_3x_ser_p2, NULL} // 2
};

GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ser_dg_gr_maxwell_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_z_3x_ser_p1, dg_gr_maxwell_alpha_quad_z_3x_ser_p2, NULL} // 2
};

// alpha_c evaluated at quadrature points.
GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ten_dg_gr_maxwell_alpha_quad_x_kernels[] = {
    // 1x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_1x_tensor_p1, dg_gr_maxwell_alpha_quad_x_1x_tensor_p2, NULL
    }, // 0
    // 2x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_2x_tensor_p1, dg_gr_maxwell_alpha_quad_x_2x_tensor_p2, NULL
    }, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_x_3x_tensor_p1, dg_gr_maxwell_alpha_quad_x_3x_tensor_p2, NULL
    } // 2
};

GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ten_dg_gr_maxwell_alpha_quad_y_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    // 2x kernels
    {NULL, dg_gr_maxwell_alpha_quad_y_2x_tensor_p1, dg_gr_maxwell_alpha_quad_y_2x_tensor_p2, NULL
    }, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_y_3x_tensor_p1, dg_gr_maxwell_alpha_quad_y_3x_tensor_p2, NULL
    } // 2
};

GKYL_CU_D static const gkyl_dg_gr_maxwell_alpha_quad_kern_list
  ten_dg_gr_maxwell_alpha_quad_z_kernels[] = {
    // 1x kernels
    {NULL, NULL, NULL, NULL}, // 0
    // 2x kernels
    {NULL, NULL, NULL, NULL}, // 1
    // 3x kernels
    {NULL, dg_gr_maxwell_alpha_quad_z_3x_tensor_p1, dg_gr_maxwell_alpha_quad_z_3x_tensor_p2, NULL
    } // 2
};
