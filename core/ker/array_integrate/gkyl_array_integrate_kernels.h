#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void gkyl_array_integrate_op_none_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_none_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_none_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_none_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_none_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_none_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_abs_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);


GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_eps_grad_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);





GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x1v_gkhyb_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x2v_gkhyb_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_2x2v_gkhyb_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_3x2v_gkhyb_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out);

EXTERN_C_END

