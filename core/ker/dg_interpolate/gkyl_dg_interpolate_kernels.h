#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void dg_interpolate_1x_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_1x_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_2x_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_2x_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_2x_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_2x_ser_p2_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_3x_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_3x_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_3x_ser_p1_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_3x_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_3x_ser_p2_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_3x_ser_p2_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);


GKYL_CU_DH void dg_interpolate_vlasov_1x1v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x1v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x1v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x1v_ser_p2_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p2_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_ser_p2_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p2_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p2_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_ser_p2_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p2_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p2_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_ser_p2_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p2_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p2_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p2_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_ser_p2_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_ser_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);


GKYL_CU_DH void dg_interpolate_gyrokinetic_1x1v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x1v_ser_p1_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_1x1v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x1v_ser_p2_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p1_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p1_mu(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p2_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_1x2v_ser_p2_mu(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p1_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p1_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p1_mu(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p2_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p2_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p2_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_2x2v_ser_p2_mu(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_gyrokinetic_3x2v_ser_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_3x2v_ser_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_3x2v_ser_p1_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_3x2v_ser_p1_vpar(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_gyrokinetic_3x2v_ser_p1_mu(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);



GKYL_CU_DH void dg_interpolate_vlasov_1x1v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x1v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x2v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x2v_tensor_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_1x3v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_tensor_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_1x3v_tensor_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x2v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_tensor_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x2v_tensor_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_2x3v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_tensor_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_tensor_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_2x3v_tensor_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_x(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_y(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_z(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_vx(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_vy(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
GKYL_CU_DH void dg_interpolate_vlasov_3x3v_tensor_p1_vz(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);




EXTERN_C_END
