#pragma once 
// Generated projection (NODAL->MODAL) tables for the from-flux 
// vlasov_(boundary_)(ho_)surf* consumer kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_nod2mod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_kamap[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_kbmap[2] = { 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_out_off[3] = { 
  0, 
  3, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_out_mode[6] = { 
  0, 
  2, 
  4, 
  1, 
  3, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  1.118033988749895, 
  0.34641016151377546, 
  -0.6, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -1.118033988749895, 
  -0.34641016151377546, 
  -0.6, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_pm_a[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_v0_pm_b[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_pm_cl[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  1.118033988749895, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_v0_pm_cr[6] = { 
  -0.5, 
  -0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  -1.118033988749895, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_tensor_p1_prj_v0_meta = { .nk = 2, .np = 6, .na = 2, .nmv = 1, .no = 2, .ni = 1, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_Cw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_Vw[9] = { 
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
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_kamap[3] = { 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_kbmap[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_out_off[4] = { 
  0, 
  2, 
  4, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_out_mode[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_pm_a[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_prj_x0_pm_b[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_pm_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_prj_x0_pm_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_tensor_p1_prj_x0_meta = { .nk = 3, .np = 6, .na = 1, .nmv = 3, .no = 1, .ni = 3, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_kamap[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_kbmap[2] = { 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_out_off[3] = { 
  0, 
  3, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_out_mode[6] = { 
  0, 
  2, 
  4, 
  1, 
  3, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  1.118033988749895, 
  0.34641016151377546, 
  -0.6, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -1.118033988749895, 
  -0.34641016151377546, 
  -0.6, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_pm_a[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_v0_pm_b[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_pm_cl[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  1.118033988749895, 
  0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_v0_pm_cr[6] = { 
  -0.5, 
  -0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  -1.118033988749895, 
  -0.7745966692414833, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_tensor_p1_ho_prj_v0_meta = { .nk = 2, .np = 6, .na = 2, .nmv = 1, .no = 2, .ni = 1, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_Cw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_Vw[12] = { 
  0.3478548451374539, 
  -0.7488760957794683, 
  -0.8192431415675561, 
  0.6521451548625461, 
  -0.5542924757974403, 
  0.8192431415675563, 
  0.6521451548625461, 
  0.5542924757974403, 
  0.8192431415675563, 
  0.3478548451374539, 
  0.7488760957794683, 
  -0.8192431415675561, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_kamap[3] = { 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_kbmap[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_out_off[4] = { 
  0, 
  2, 
  4, 
  6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_out_mode[6] = { 
  0, 
  1, 
  2, 
  3, 
  4, 
  5, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_out_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_out_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_pm_a[6] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_prj_x0_pm_b[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_pm_cl[6] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
  -0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_prj_x0_pm_cr[6] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
  0.2906888370749726, 
  0.503487835006964, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_tensor_p1_ho_prj_x0_meta = { .nk = 3, .np = 6, .na = 1, .nmv = 3, .no = 1, .ni = 4, .fhat_off = 0 }; 
