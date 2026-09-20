#pragma once 
// Generated projection (NODAL->MODAL) tables for the from-flux 
// vlasov_(boundary_)(ho_)surf* consumer kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_nod2mod_meta.h> 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_kamap[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_kbmap[2] = { 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_out_off[3] = { 
  0, 
  2, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_out_mode[4] = { 
  0, 
  2, 
  1, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_out_cl[4] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_out_cr[4] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_pm_a[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_v0_pm_b[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_pm_cl[4] = { 
  0.5, 
  0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_v0_pm_cr[4] = { 
  -0.5, 
  -0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_ser_p1_prj_v0_meta = { .nk = 2, .np = 4, .na = 2, .nmv = 1, .no = 2, .ni = 1, .fhat_off = 0 }; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_Cw[1] = { 
  1.0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_Vw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_kamap[2] = { 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_kbmap[2] = { 
  0, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_out_off[3] = { 
  0, 
  2, 
  4, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_out_mode[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_out_cl[4] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_out_cr[4] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_pm_a[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const int vst_1x1v_ser_p1_prj_x0_pm_b[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_pm_cl[4] = { 
  0.5, 
  -0.8660254037844386, 
  0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const double vst_1x1v_ser_p1_prj_x0_pm_cr[4] = { 
  -0.5, 
  -0.8660254037844386, 
  -0.34641016151377546, 
  -0.6, 
}; 
GKYL_CU_TABLE __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_ser_p1_prj_x0_meta = { .nk = 2, .np = 4, .na = 1, .nmv = 2, .no = 1, .ni = 2, .fhat_off = 0 }; 
