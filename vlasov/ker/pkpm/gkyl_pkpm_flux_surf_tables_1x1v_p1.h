#pragma once 
// Generated nodal surface-evaluation tables of the PKPM flux-scatter kernels 
// (see gkylcas/pkpm/pkpm/pkpm-flux-tables.mac): five Gauss-Legendre points per 
// transverse configuration direction and three in v_par on configuration faces, five 
// per configuration direction on the v_par face. Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays are 
// elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_nod2mod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_Cw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_Vw[9] = { 
  0.5555555555555556, 
  -1.0758287072798385, 
  -0.8547008547008547, 
  0.8888888888888888, 
  0.0, 
  1.7094017094017093, 
  0.5555555555555556, 
  1.0758287072798385, 
  -0.8547008547008547, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_kamap[3] = { 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_kbmap[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_boff[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_bks[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_out_off[4] = { 
  0, 
  2, 
  4, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_out_mode[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_pm_s[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_pm_a[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_x0_pm_b[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_pm_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_x0_pm_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_pkpm_1x1v_p1_prj_x0_meta = { .nk = 3, .np = 6, .na = 1, .nmv = 3, .no = 1, .ni = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_Cw[10] = { 
  0.23692688505618906, 
  -0.5367459204973625, 
  0.4786286704993664, 
  -0.6443171250014854, 
  0.5688888888888889, 
  0.0, 
  0.4786286704993664, 
  0.6443171250014854, 
  0.23692688505618906, 
  0.5367459204973625, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_kamap[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_kbmap[2] = { 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_boff[2] = { 
  0, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_bks[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_out_off[3] = { 
  0, 
  3, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_out_mode[6] = { 
  0, 
  2, 
  4, 
  1, 
  3, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  1.118033988749895, 
  0.34641016151377546, 
  -0.6, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -1.118033988749895, 
  -0.34641016151377546, 
  -0.6, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_pm_s[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_pm_a[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_pkpm_1x1v_p1_prj_v0_pm_b[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_pm_cl[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  1.118033988749895, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_pkpm_1x1v_p1_prj_v0_pm_cr[6] = { 
  -0.5, 
  -0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  -1.118033988749895, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_pkpm_1x1v_p1_prj_v0_meta = { .nk = 2, .np = 6, .na = 2, .nmv = 1, .no = 5, .ni = 1, .fhat_off = 0 }; 
