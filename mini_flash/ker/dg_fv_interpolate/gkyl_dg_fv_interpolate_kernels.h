#pragma once

#include <gkyl_util.h>

void dg_fv_interpolate_fv_to_modal_1x(const double *fv_data[8], double *dg_out); 
void dg_fv_interpolate_fv_to_modal_2x(const double *fv_data[8], double *dg_out); 
void dg_fv_interpolate_fv_to_modal_3x(const double *fv_data[8], double *dg_out); 

void dg_fv_interpolate_modal_to_fv_1x(const double *dg_data, double *fv_out[8]); 
void dg_fv_interpolate_modal_to_fv_2x(const double *dg_data, double *fv_out[8]); 
void dg_fv_interpolate_modal_to_fv_3x(const double *dg_data, double *fv_out[8]); 
