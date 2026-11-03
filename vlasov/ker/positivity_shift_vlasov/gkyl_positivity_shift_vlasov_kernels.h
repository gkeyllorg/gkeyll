#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH bool positivity_shift_vlasov_conf_pos_check_1x_tensor_p1(const double *fld);

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_1x1v_tensor_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_1x1v_tensor_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_1x2v_tensor_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_1x2v_tensor_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_1x3v_tensor_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_1x3v_tensor_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_vlasov_conf_pos_check_2x_tensor_p1(const double *fld);

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_2x2v_tensor_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_2x2v_tensor_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_2x3v_tensor_p1(double ffloor, double *distf);
GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_2x3v_tensor_p1(double ffloor, double *distf);

GKYL_CU_DH bool positivity_shift_vlasov_conf_pos_check_3x_tensor_p1(const double *fld);



EXTERN_C_END
