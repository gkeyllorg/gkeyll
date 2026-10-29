#pragma once

#include <gkyl_dg_fv_interpolate_kernels.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>

// Function type for p=1 modal to finite volume conversion.
typedef void (*dg_fv_interp_modal_to_fv_t)(const double *dg_data, double *fv_out[8]);

// Function type for finite volume to p=1 modal conversion.
typedef void (*dg_fv_interp_fv_to_modal_t)(const double *fv_data[8], double *dg_out);

struct gkyl_dg_fv_interpolate_kernels {
    dg_fv_interp_modal_to_fv_t modal_to_fv;
    dg_fv_interp_fv_to_modal_t fv_to_modal;
};

struct gkyl_dg_fv_interpolate {
    int ndim; // Dimensionality.
    int num_basis; // Number of basis functions.
    struct gkyl_dg_fv_interpolate_kernels kernels;
};

// Kernel lists
static const dg_fv_interp_modal_to_fv_t dg_fv_interp_modal_to_fv_kernels[3] = {
    dg_fv_interpolate_modal_to_fv_1x,
    dg_fv_interpolate_modal_to_fv_2x,
    dg_fv_interpolate_modal_to_fv_3x
};

static const dg_fv_interp_fv_to_modal_t dg_fv_interp_fv_to_modal_kernels[3] = {
    dg_fv_interpolate_fv_to_modal_1x,
    dg_fv_interpolate_fv_to_modal_2x,
    dg_fv_interpolate_fv_to_modal_3x
};

static dg_fv_interp_modal_to_fv_t
choose_dg_fv_interp_modal_to_fv_kernel(int ndim) {
    return dg_fv_interp_modal_to_fv_kernels[ndim-1];
}

static dg_fv_interp_fv_to_modal_t
choose_dg_fv_interp_fv_to_modal_kernel(int ndim) {
    return dg_fv_interp_fv_to_modal_kernels[ndim-1];
}