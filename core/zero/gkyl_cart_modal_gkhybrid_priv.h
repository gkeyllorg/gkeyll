#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_gkhyb_kernels.h>

// Basis function eval for each dimension: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*ev[4])(const double *z, double *b);
} ev_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_1x1v_gkhyb_p1, NULL, NULL }, { NULL, eval_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, eval_2x2v_gkhyb_p1, NULL, NULL }, { NULL, eval_3x2v_gkhyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL } // No 6D basis functions
};

// Expansion eval for each dimension: eve_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  double (*ev[4])(const double *z, const double *f);
} eve_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_expand_1x1v_gkhyb_p1, NULL, NULL }, { NULL, eval_expand_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, eval_expand_2x2v_gkhyb_p1, NULL, NULL }, { NULL, eval_expand_3x2v_gkhyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL } // No 6D basis functions
};

// Expansion eval_grad for each dimension: eveg_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  double (*ev[4])(int dir, const double *z, const double *f);
} eveg_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, eval_grad_expand_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, eval_grad_expand_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, eval_grad_expand_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, eval_grad_expand_3x2v_gkhyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL } // No 6D basis functions
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*fs[4])(int dir, const double *f, double *fout);
} fos_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, flip_odd_sign_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, flip_odd_sign_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, flip_odd_sign_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, flip_odd_sign_3x2v_gkhyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL } // No 6D basis functions
};

// Flip-sign functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*fs[4])(int dir, const double *f, double *fout);
} fes_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, flip_even_sign_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, flip_even_sign_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, flip_even_sign_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, flip_even_sign_3x2v_gkhyb_p1, NULL, NULL },
  { NULL, NULL, NULL, NULL } // No 6D basis functions
};

// Number of basis functions: num_basis_list[ndim].count[poly_order]
GKYL_CU_D static struct {
  int count[4];
} num_basis_list[] = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 6, 0, 0 }, { 0, 12, 0, 0 },
  { 0, 24, 0, 0 }, { 0, 48, 0, 0 } };

// Node list function: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*nl[4])(double *node_list);
} nl_list[] = { { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, node_coords_1x1v_gkhyb_p1, NULL, NULL }, { NULL, node_coords_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, node_coords_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, node_coords_3x2v_gkhyb_p1, NULL, NULL } };

// Nodal -> modal conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*n2m[4])(const double *fnodal, double *fmodal);
} n2m_list[] = { { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, nodal_to_modal_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, nodal_to_modal_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, nodal_to_modal_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, nodal_to_modal_3x2v_gkhyb_p1, NULL, NULL } };

// Gauss-Legendre quadrature nodes nodal basis -> modal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*n2m[4])(const double *fquad, double *fmodal, long linc2);
} qn2m_list[] = { { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, quad_to_modal_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, quad_to_modal_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, quad_to_modal_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, quad_to_modal_3x2v_gkhyb_p1, NULL, NULL } };

// Number of quadrature nodes: num_quad_list[ndim].count[poly_order]
GKYL_CU_D static struct {
  int count[4];
} num_quad_list[] = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 6, 0, 0 }, { 0, 12, 0, 0 },
  { 0, 24, 0, 0 }, { 0, 48, 0, 0 } };

// modal basis -> Gauss-Legendre quadrature nodes nodal basis conversion functions: ev_list[ndim].ev[poly_order]
GKYL_CU_D static struct {
  void (*n2m[4])(const double *fmodal, double *fquad, long linc2);
} m2qn_list[] = { { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { NULL, NULL, NULL, NULL }, // No 1D basis functions
  { NULL, modal_to_quad_1x1v_gkhyb_p1, NULL, NULL },
  { NULL, modal_to_quad_1x2v_gkhyb_p1, NULL, NULL },
  { NULL, modal_to_quad_2x2v_gkhyb_p1, NULL, NULL },
  { NULL, modal_to_quad_3x2v_gkhyb_p1, NULL, NULL } };