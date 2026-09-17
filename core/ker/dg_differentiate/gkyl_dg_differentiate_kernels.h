#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void dg_differentiate_local_ser_1x_p1_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_1x_p1_dir0_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_ser_1x_p2_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_1x_p2_dir0_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_ser_2x_p1_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p1_dir0_order2(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p1_dir1_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p1_dir1_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_ser_2x_p2_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p2_dir0_order2(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p2_dir1_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_2x_p2_dir1_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir0_order2(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir1_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir1_order2(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir2_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_ser_3x_p1_dir2_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_tensor_1x_p2_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_tensor_1x_p2_dir0_order2(double dx, const double *fin, double *out);

GKYL_CU_DH void dg_differentiate_local_tensor_2x_p2_dir0_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_tensor_2x_p2_dir0_order2(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_tensor_2x_p2_dir1_order1(double dx, const double *fin, double *out);
GKYL_CU_DH void dg_differentiate_local_tensor_2x_p2_dir1_order2(double dx, const double *fin, double *out);


EXTERN_C_END
