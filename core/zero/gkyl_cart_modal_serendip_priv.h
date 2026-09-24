#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_ser_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*ev[4])(const double *z, double *b); } ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_1d_ser_p0, eval_1d_ser_p1, eval_1d_ser_p2, eval_1d_ser_p3 },
  { eval_2d_ser_p0, eval_2d_ser_p1, eval_2d_ser_p2, eval_2d_ser_p3 },
  { eval_3d_ser_p0, eval_3d_ser_p1, eval_3d_ser_p2, eval_3d_ser_p3 },
  { eval_4d_ser_p0, eval_4d_ser_p1, eval_4d_ser_p2, eval_4d_ser_p3 },
  { eval_5d_ser_p0, eval_5d_ser_p1, eval_5d_ser_p2, NULL },
  { eval_6d_ser_p0, eval_6d_ser_p1, NULL, NULL },
};

// Expansion eval for each dimension: eve_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { double (*ev[4])(const double *z, const double *f); } eve_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_expand_1d_ser_p0, eval_expand_1d_ser_p1, eval_expand_1d_ser_p2, eval_expand_1d_ser_p3 },
  { eval_expand_2d_ser_p0, eval_expand_2d_ser_p1, eval_expand_2d_ser_p2, eval_expand_2d_ser_p3 },
  { eval_expand_3d_ser_p0, eval_expand_3d_ser_p1, eval_expand_3d_ser_p2, eval_expand_3d_ser_p3 },
  { eval_expand_4d_ser_p0, eval_expand_4d_ser_p1, eval_expand_4d_ser_p2, eval_expand_4d_ser_p3 },
  { eval_expand_5d_ser_p0, eval_expand_5d_ser_p1, eval_expand_5d_ser_p2, NULL },
  { eval_expand_6d_ser_p0, eval_expand_6d_ser_p1, NULL, NULL },
};

// Expansion eval_grad for each dimension: eveg_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { double (*ev[4])(int dir, const double *z, const double *f); } eveg_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { eval_grad_expand_1d_ser_p0, eval_grad_expand_1d_ser_p1, eval_grad_expand_1d_ser_p2, eval_grad_expand_1d_ser_p3 },
  { eval_grad_expand_2d_ser_p0, eval_grad_expand_2d_ser_p1, eval_grad_expand_2d_ser_p2, eval_grad_expand_2d_ser_p3 },
  { eval_grad_expand_3d_ser_p0, eval_grad_expand_3d_ser_p1, eval_grad_expand_3d_ser_p2, eval_grad_expand_3d_ser_p3 },
  { eval_grad_expand_4d_ser_p0, eval_grad_expand_4d_ser_p1, eval_grad_expand_4d_ser_p2, eval_grad_expand_4d_ser_p3 },
  { eval_grad_expand_5d_ser_p0, eval_grad_expand_5d_ser_p1, eval_grad_expand_5d_ser_p2, NULL },
  { eval_grad_expand_6d_ser_p0, eval_grad_expand_6d_ser_p1, NULL, NULL },
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fos_list[] = {
  {NULL, NULL, NULL, NULL}, // No 0D basis functions
  {flip_odd_sign_1d_ser_p0, flip_odd_sign_1d_ser_p1, flip_odd_sign_1d_ser_p2, flip_odd_sign_1d_ser_p3},
  {flip_odd_sign_2d_ser_p0, flip_odd_sign_2d_ser_p1, flip_odd_sign_2d_ser_p2, flip_odd_sign_2d_ser_p3},
  {flip_odd_sign_3d_ser_p0, flip_odd_sign_3d_ser_p1, flip_odd_sign_3d_ser_p2, flip_odd_sign_3d_ser_p3},
  {flip_odd_sign_4d_ser_p0, flip_odd_sign_4d_ser_p1, flip_odd_sign_4d_ser_p2, flip_odd_sign_4d_ser_p3},
  {flip_odd_sign_5d_ser_p0, flip_odd_sign_5d_ser_p1, flip_odd_sign_5d_ser_p2, NULL},
  {flip_odd_sign_6d_ser_p0, flip_odd_sign_6d_ser_p1, NULL, NULL},
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*fs[4])(int dir, const double *f, double *fout); } fes_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { flip_even_sign_1d_ser_p0, flip_even_sign_1d_ser_p1, flip_even_sign_1d_ser_p2, flip_even_sign_1d_ser_p3 },
  { flip_even_sign_2d_ser_p0, flip_even_sign_2d_ser_p1, flip_even_sign_2d_ser_p2, flip_even_sign_2d_ser_p3 },
  { flip_even_sign_3d_ser_p0, flip_even_sign_3d_ser_p1, flip_even_sign_3d_ser_p2, flip_even_sign_3d_ser_p3 },
  { flip_even_sign_4d_ser_p0, flip_even_sign_4d_ser_p1, flip_even_sign_4d_ser_p2, flip_even_sign_4d_ser_p3 },
  { flip_even_sign_5d_ser_p0, flip_even_sign_5d_ser_p1, flip_even_sign_5d_ser_p2, NULL },
  { flip_even_sign_6d_ser_p0, flip_even_sign_6d_ser_p1, NULL, NULL },
};

// Number of basis functions: num_basis_list[ndim].count[poly_order]
GKYL_CU_D
static struct { int count[4]; } num_basis_list[] = {
  { 1, 1, 1, 1 },
  { 1, 2, 3, 4 },
  { 1, 4, 8, 12 },
  { 1, 8, 20, 32 },
  { 1, 16, 48, 80 },
  { 1, 32, 112, 192 },
  { 1, 64, 256, 0 },
};

// Node list function: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*nl[4])(double * node_list); } nl_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { node_coords_1d_ser_p0, node_coords_1d_ser_p1, node_coords_1d_ser_p2, NULL },
  { node_coords_2d_ser_p0, node_coords_2d_ser_p1, node_coords_2d_ser_p2, NULL },
  { node_coords_3d_ser_p0, node_coords_3d_ser_p1, node_coords_3d_ser_p2, NULL },
  { node_coords_4d_ser_p0, node_coords_4d_ser_p1, node_coords_4d_ser_p2, NULL },
  { node_coords_5d_ser_p0, node_coords_5d_ser_p1, node_coords_5d_ser_p2, NULL },
  { node_coords_6d_ser_p0, node_coords_6d_ser_p1, NULL, NULL },
};

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, nodal_to_modal_1d_ser_p1, nodal_to_modal_1d_ser_p2, NULL },
  { NULL, nodal_to_modal_2d_ser_p1, nodal_to_modal_2d_ser_p2, NULL },
  { NULL, nodal_to_modal_3d_ser_p1, nodal_to_modal_3d_ser_p2, NULL },
  { NULL, nodal_to_modal_4d_ser_p1, nodal_to_modal_4d_ser_p2, NULL },
  { NULL, nodal_to_modal_5d_ser_p1, nodal_to_modal_5d_ser_p2, NULL },
  { NULL, nodal_to_modal_6d_ser_p1, NULL, NULL },
};

// Nodal -> modal conversion functions for nodes on a surface in one direction
// and Gauss-Legendre in the others: ev_list[ndim].ev[poly_order].indir[dir]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_quad_surf_list_x[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, nodal_to_modal_quad_surfx_2d_ser_p1, nodal_to_modal_quad_surfx_2d_ser_p2, NULL },
  { NULL, nodal_to_modal_quad_surfx_3d_ser_p1, nodal_to_modal_quad_surfx_3d_ser_p2, NULL },
};
GKYL_CU_D
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_quad_surf_list_y[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, nodal_to_modal_quad_surfy_2d_ser_p1, nodal_to_modal_quad_surfy_2d_ser_p2, NULL },
  { NULL, nodal_to_modal_quad_surfy_3d_ser_p1, nodal_to_modal_quad_surfy_3d_ser_p2, NULL },
};
GKYL_CU_D
static struct { void (*n2m[4])(const double *fnodal, double *fmodal); } n2m_quad_surf_list_z[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, NULL, NULL, NULL }, // No 2D basis functions
  { NULL, nodal_to_modal_quad_surfz_3d_ser_p1, nodal_to_modal_quad_surfz_3d_ser_p2, NULL },
};

// Node list function: 
// List of nodes on a surface in one direction
// and Gauss-Legendre in the others: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*nl[4])(double * node_list); } nl_quad_surf_list_x[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, node_coords_quad_surfx_2d_ser_p1, node_coords_quad_surfx_2d_ser_p2, NULL },
  { NULL, node_coords_quad_surfx_3d_ser_p1, node_coords_quad_surfx_3d_ser_p2, NULL },
};
GKYL_CU_D
static struct { void (*nl[4])(double * node_list); } nl_quad_surf_list_y[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, node_coords_quad_surfy_2d_ser_p1, node_coords_quad_surfy_2d_ser_p2, NULL },
  { NULL, node_coords_quad_surfy_3d_ser_p1, node_coords_quad_surfy_3d_ser_p2, NULL },
};
GKYL_CU_D
static struct { void (*nl[4])(double * node_list); } nl_quad_surf_list_z[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, NULL, NULL, NULL }, // No 2D basis functions
  { NULL, node_coords_quad_surfz_3d_ser_p1, node_coords_quad_surfz_3d_ser_p2, NULL },
};

// Gauss-Legendre quadrature nodes nodal basis -> modal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fquad, double *fmodal, long linc2); } qn2m_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, quad_to_modal_1d_ser_p1, quad_to_modal_1d_ser_p2, NULL },
  { NULL, quad_to_modal_2d_ser_p1, quad_to_modal_2d_ser_p2, NULL },
  { NULL, quad_to_modal_3d_ser_p1, quad_to_modal_3d_ser_p2, NULL },
  { NULL, quad_to_modal_4d_ser_p1, quad_to_modal_4d_ser_p2, NULL },
  { NULL, quad_to_modal_5d_ser_p1, quad_to_modal_5d_ser_p2, NULL },
  { NULL, quad_to_modal_6d_ser_p1, NULL, NULL },
};

// Number of quadrature nodes: num_quad_list[ndim].count[poly_order]
GKYL_CU_D
static struct { int count[4]; } num_quad_list[] = {
  { 0, 0, 0, 0 },
  { 0, 2, 3, 0 },
  { 0, 4, 9, 0 },
  { 0, 8, 27, 0 },
  { 0, 16, 81, 0 },
  { 0, 32, 243, 0 },
  { 0, 64, 0, 0 },
};

// modal basis -> Gauss-Legendre quadrature nodes nodal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D
static struct { void (*n2m[4])(const double *fmodal, double *fquad, long linc2); } m2qn_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, modal_to_quad_1d_ser_p1, modal_to_quad_1d_ser_p2, NULL },
  { NULL, modal_to_quad_2d_ser_p1, modal_to_quad_2d_ser_p2, NULL },
  { NULL, modal_to_quad_3d_ser_p1, modal_to_quad_3d_ser_p2, NULL },
  { NULL, modal_to_quad_4d_ser_p1, modal_to_quad_4d_ser_p2, NULL },
  { NULL, modal_to_quad_5d_ser_p1, modal_to_quad_5d_ser_p2, NULL },
  { NULL, modal_to_quad_6d_ser_p1, NULL, NULL },
};
