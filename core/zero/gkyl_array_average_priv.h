#pragma once

// Private header for the array_average updater, not for direct use by user.

#include <gkyl_array_average.h>
#include <gkyl_array_average_kernels.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_util.h>
#include <assert.h>
#include <math.h>

// Function pointer type for array_average kernels.
typedef void (*array_average_t)(
  const double subvol, const double *win, const double *fin, double *GKYL_RESTRICT out
);

// For use in kernel tables.
typedef struct {
  array_average_t kernels[2]; // Reduced-dim and full-dim output kernels.
} array_average_kern_list;
typedef struct {
  array_average_kern_list ops[63]; // Indexed by the averaged dimensions op code minus 1.
} array_average_op_kern_list;
typedef struct {
  array_average_op_kern_list list[2]; // Indexed by poly_order-1.
} array_average_p_kern_list;
typedef struct {
  array_average_op_kern_list vdim[3]; // Indexed by vdim-1.
} array_average_v_kern_list;

// Conf-space serendipity kernels, indexed [ndim-1].list[poly_order-1].ops[op-1].
GKYL_CU_D static const array_average_p_kern_list gkyl_array_average_ker_list_ser[] = {
  // 1x
  {
    .list =
      {{.ops = {{gkyl_array_average_1x_ser_p1_avgx, gkyl_array_average_1x_ser_p1_avgx_fulldim}}},
       {.ops = {{gkyl_array_average_1x_ser_p2_avgx, gkyl_array_average_1x_ser_p2_avgx_fulldim}}}},
  },
  // 2x
  {
    .list =
      {{
         .ops =
           {{gkyl_array_average_2x_ser_p1_avgx, gkyl_array_average_2x_ser_p1_avgx_fulldim},
            {gkyl_array_average_2x_ser_p1_avgy, gkyl_array_average_2x_ser_p1_avgy_fulldim},
            {gkyl_array_average_2x_ser_p1_avgxy, gkyl_array_average_2x_ser_p1_avgxy_fulldim}},
       },
       {
         .ops =
           {{gkyl_array_average_2x_ser_p2_avgx, gkyl_array_average_2x_ser_p2_avgx_fulldim},
            {gkyl_array_average_2x_ser_p2_avgy, gkyl_array_average_2x_ser_p2_avgy_fulldim},
            {gkyl_array_average_2x_ser_p2_avgxy, gkyl_array_average_2x_ser_p2_avgxy_fulldim}},
       }},
  },
  // 3x
  {
    .list =
      {{
         .ops =
           {{gkyl_array_average_3x_ser_p1_avgx, gkyl_array_average_3x_ser_p1_avgx_fulldim},
            {gkyl_array_average_3x_ser_p1_avgy, gkyl_array_average_3x_ser_p1_avgy_fulldim},
            {gkyl_array_average_3x_ser_p1_avgxy, gkyl_array_average_3x_ser_p1_avgxy_fulldim},
            {gkyl_array_average_3x_ser_p1_avgz, gkyl_array_average_3x_ser_p1_avgz_fulldim},
            {gkyl_array_average_3x_ser_p1_avgxz, gkyl_array_average_3x_ser_p1_avgxz_fulldim},
            {gkyl_array_average_3x_ser_p1_avgyz, gkyl_array_average_3x_ser_p1_avgyz_fulldim},
            {gkyl_array_average_3x_ser_p1_avgxyz, gkyl_array_average_3x_ser_p1_avgxyz_fulldim}},
       },
       {
         .ops =
           {{gkyl_array_average_3x_ser_p2_avgx, gkyl_array_average_3x_ser_p2_avgx_fulldim},
            {gkyl_array_average_3x_ser_p2_avgy, gkyl_array_average_3x_ser_p2_avgy_fulldim},
            {gkyl_array_average_3x_ser_p2_avgxy, gkyl_array_average_3x_ser_p2_avgxy_fulldim},
            {gkyl_array_average_3x_ser_p2_avgz, gkyl_array_average_3x_ser_p2_avgz_fulldim},
            {gkyl_array_average_3x_ser_p2_avgxz, gkyl_array_average_3x_ser_p2_avgxz_fulldim},
            {gkyl_array_average_3x_ser_p2_avgyz, gkyl_array_average_3x_ser_p2_avgyz_fulldim},
            {gkyl_array_average_3x_ser_p2_avgxyz, gkyl_array_average_3x_ser_p2_avgxyz_fulldim}},
       }},
  }
};

// Gyrokinetic hybrid (p=1) kernels, indexed [ndim-2].ops[op-1].
GKYL_CU_D static const array_average_op_kern_list gkyl_array_average_ker_list_gkhyb[] = {
  // 1x1v
  {
    .ops =
      {{gkyl_array_average_1x1v_gkhyb_p1_avgx, gkyl_array_average_1x1v_gkhyb_p1_avgx_fulldim},
       {gkyl_array_average_1x1v_gkhyb_p1_avgvpar, gkyl_array_average_1x1v_gkhyb_p1_avgvpar_fulldim},
       {gkyl_array_average_1x1v_gkhyb_p1_avgxvpar, gkyl_array_average_1x1v_gkhyb_p1_avgxvpar_fulldim
       }},
  },
  // 1x2v
  {
    .ops =
      {{gkyl_array_average_1x2v_gkhyb_p1_avgx, gkyl_array_average_1x2v_gkhyb_p1_avgx_fulldim},
       {gkyl_array_average_1x2v_gkhyb_p1_avgvpar, gkyl_array_average_1x2v_gkhyb_p1_avgvpar_fulldim},
       {gkyl_array_average_1x2v_gkhyb_p1_avgxvpar, gkyl_array_average_1x2v_gkhyb_p1_avgxvpar_fulldim
       },
       {gkyl_array_average_1x2v_gkhyb_p1_avgmu, gkyl_array_average_1x2v_gkhyb_p1_avgmu_fulldim},
       {gkyl_array_average_1x2v_gkhyb_p1_avgxmu, gkyl_array_average_1x2v_gkhyb_p1_avgxmu_fulldim},
       {gkyl_array_average_1x2v_gkhyb_p1_avgvparmu,
        gkyl_array_average_1x2v_gkhyb_p1_avgvparmu_fulldim},
       {gkyl_array_average_1x2v_gkhyb_p1_avgxvparmu,
        gkyl_array_average_1x2v_gkhyb_p1_avgxvparmu_fulldim}},
  },
  // 2x2v
  {
    .ops =
      {{gkyl_array_average_2x2v_gkhyb_p1_avgx, gkyl_array_average_2x2v_gkhyb_p1_avgx_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgy, gkyl_array_average_2x2v_gkhyb_p1_avgy_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxy, gkyl_array_average_2x2v_gkhyb_p1_avgxy_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgvpar, gkyl_array_average_2x2v_gkhyb_p1_avgvpar_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxvpar, gkyl_array_average_2x2v_gkhyb_p1_avgxvpar_fulldim
       },
       {gkyl_array_average_2x2v_gkhyb_p1_avgyvpar, gkyl_array_average_2x2v_gkhyb_p1_avgyvpar_fulldim
       },
       {gkyl_array_average_2x2v_gkhyb_p1_avgxyvpar,
        gkyl_array_average_2x2v_gkhyb_p1_avgxyvpar_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgmu, gkyl_array_average_2x2v_gkhyb_p1_avgmu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxmu, gkyl_array_average_2x2v_gkhyb_p1_avgxmu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgymu, gkyl_array_average_2x2v_gkhyb_p1_avgymu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxymu, gkyl_array_average_2x2v_gkhyb_p1_avgxymu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgvparmu,
        gkyl_array_average_2x2v_gkhyb_p1_avgvparmu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxvparmu,
        gkyl_array_average_2x2v_gkhyb_p1_avgxvparmu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgyvparmu,
        gkyl_array_average_2x2v_gkhyb_p1_avgyvparmu_fulldim},
       {gkyl_array_average_2x2v_gkhyb_p1_avgxyvparmu,
        gkyl_array_average_2x2v_gkhyb_p1_avgxyvparmu_fulldim}},
  },
  // 3x2v
  {
    .ops =
      {{gkyl_array_average_3x2v_gkhyb_p1_avgx, gkyl_array_average_3x2v_gkhyb_p1_avgx_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgy, gkyl_array_average_3x2v_gkhyb_p1_avgy_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxy, gkyl_array_average_3x2v_gkhyb_p1_avgxy_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgz, gkyl_array_average_3x2v_gkhyb_p1_avgz_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxz, gkyl_array_average_3x2v_gkhyb_p1_avgxz_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgyz, gkyl_array_average_3x2v_gkhyb_p1_avgyz_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyz, gkyl_array_average_3x2v_gkhyb_p1_avgxyz_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgvpar, gkyl_array_average_3x2v_gkhyb_p1_avgvpar_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxvpar, gkyl_array_average_3x2v_gkhyb_p1_avgxvpar_fulldim
       },
       {gkyl_array_average_3x2v_gkhyb_p1_avgyvpar, gkyl_array_average_3x2v_gkhyb_p1_avgyvpar_fulldim
       },
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyvpar,
        gkyl_array_average_3x2v_gkhyb_p1_avgxyvpar_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgzvpar, gkyl_array_average_3x2v_gkhyb_p1_avgzvpar_fulldim
       },
       {gkyl_array_average_3x2v_gkhyb_p1_avgxzvpar,
        gkyl_array_average_3x2v_gkhyb_p1_avgxzvpar_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgyzvpar,
        gkyl_array_average_3x2v_gkhyb_p1_avgyzvpar_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyzvpar,
        gkyl_array_average_3x2v_gkhyb_p1_avgxyzvpar_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgmu, gkyl_array_average_3x2v_gkhyb_p1_avgmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxmu, gkyl_array_average_3x2v_gkhyb_p1_avgxmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgymu, gkyl_array_average_3x2v_gkhyb_p1_avgymu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxymu, gkyl_array_average_3x2v_gkhyb_p1_avgxymu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgzmu, gkyl_array_average_3x2v_gkhyb_p1_avgzmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxzmu, gkyl_array_average_3x2v_gkhyb_p1_avgxzmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgyzmu, gkyl_array_average_3x2v_gkhyb_p1_avgyzmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyzmu, gkyl_array_average_3x2v_gkhyb_p1_avgxyzmu_fulldim
       },
       {gkyl_array_average_3x2v_gkhyb_p1_avgvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgxvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgyvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgyvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgxyvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgzvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgzvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxzvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgxzvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgyzvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgyzvparmu_fulldim},
       {gkyl_array_average_3x2v_gkhyb_p1_avgxyzvparmu,
        gkyl_array_average_3x2v_gkhyb_p1_avgxyzvparmu_fulldim}},
  }
};

// Vlasov hybrid (p=1) kernels, indexed [cdim-1].vdim[vdim-1].ops[op-1].
GKYL_CU_D static const array_average_v_kern_list gkyl_array_average_ker_list_hyb[] = {
  // 1x
  {
    .vdim =
      {// 1x1v
       {
         .ops =
           {{gkyl_array_average_1x1v_hyb_p1_avgx, gkyl_array_average_1x1v_hyb_p1_avgx_fulldim},
            {gkyl_array_average_1x1v_hyb_p1_avgvx, gkyl_array_average_1x1v_hyb_p1_avgvx_fulldim},
            {gkyl_array_average_1x1v_hyb_p1_avgxvx, gkyl_array_average_1x1v_hyb_p1_avgxvx_fulldim}},
       },
       // 1x2v
       {
         .ops =
           {{gkyl_array_average_1x2v_hyb_p1_avgx, gkyl_array_average_1x2v_hyb_p1_avgx_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgvx, gkyl_array_average_1x2v_hyb_p1_avgvx_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgxvx, gkyl_array_average_1x2v_hyb_p1_avgxvx_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgvy, gkyl_array_average_1x2v_hyb_p1_avgvy_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgxvy, gkyl_array_average_1x2v_hyb_p1_avgxvy_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgvxvy, gkyl_array_average_1x2v_hyb_p1_avgvxvy_fulldim},
            {gkyl_array_average_1x2v_hyb_p1_avgxvxvy,
             gkyl_array_average_1x2v_hyb_p1_avgxvxvy_fulldim}},
       },
       // 1x3v
       {
         .ops =
           {{gkyl_array_average_1x3v_hyb_p1_avgx, gkyl_array_average_1x3v_hyb_p1_avgx_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvx, gkyl_array_average_1x3v_hyb_p1_avgvx_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvx, gkyl_array_average_1x3v_hyb_p1_avgxvx_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvy, gkyl_array_average_1x3v_hyb_p1_avgvy_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvy, gkyl_array_average_1x3v_hyb_p1_avgxvy_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvxvy, gkyl_array_average_1x3v_hyb_p1_avgvxvy_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvxvy,
             gkyl_array_average_1x3v_hyb_p1_avgxvxvy_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvz, gkyl_array_average_1x3v_hyb_p1_avgvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvz, gkyl_array_average_1x3v_hyb_p1_avgxvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvxvz, gkyl_array_average_1x3v_hyb_p1_avgvxvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvxvz,
             gkyl_array_average_1x3v_hyb_p1_avgxvxvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvyvz, gkyl_array_average_1x3v_hyb_p1_avgvyvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvyvz,
             gkyl_array_average_1x3v_hyb_p1_avgxvyvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgvxvyvz,
             gkyl_array_average_1x3v_hyb_p1_avgvxvyvz_fulldim},
            {gkyl_array_average_1x3v_hyb_p1_avgxvxvyvz,
             gkyl_array_average_1x3v_hyb_p1_avgxvxvyvz_fulldim}},
       }
      },
  },
  // 2x
  {
    .vdim =
      {// 2x1v (none)
       {.ops = {{NULL, NULL}}},
       // 2x2v
       {
         .ops =
           {{gkyl_array_average_2x2v_hyb_p1_avgx, gkyl_array_average_2x2v_hyb_p1_avgx_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgy, gkyl_array_average_2x2v_hyb_p1_avgy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxy, gkyl_array_average_2x2v_hyb_p1_avgxy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgvx, gkyl_array_average_2x2v_hyb_p1_avgvx_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxvx, gkyl_array_average_2x2v_hyb_p1_avgxvx_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgyvx, gkyl_array_average_2x2v_hyb_p1_avgyvx_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxyvx, gkyl_array_average_2x2v_hyb_p1_avgxyvx_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgvy, gkyl_array_average_2x2v_hyb_p1_avgvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxvy, gkyl_array_average_2x2v_hyb_p1_avgxvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgyvy, gkyl_array_average_2x2v_hyb_p1_avgyvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxyvy, gkyl_array_average_2x2v_hyb_p1_avgxyvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgvxvy, gkyl_array_average_2x2v_hyb_p1_avgvxvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxvxvy,
             gkyl_array_average_2x2v_hyb_p1_avgxvxvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgyvxvy,
             gkyl_array_average_2x2v_hyb_p1_avgyvxvy_fulldim},
            {gkyl_array_average_2x2v_hyb_p1_avgxyvxvy,
             gkyl_array_average_2x2v_hyb_p1_avgxyvxvy_fulldim}},
       },
       // 2x3v
       {
         .ops =
           {
             {gkyl_array_average_2x3v_hyb_p1_avgx, gkyl_array_average_2x3v_hyb_p1_avgx_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgy, gkyl_array_average_2x3v_hyb_p1_avgy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxy, gkyl_array_average_2x3v_hyb_p1_avgxy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgvx, gkyl_array_average_2x3v_hyb_p1_avgvx_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxvx, gkyl_array_average_2x3v_hyb_p1_avgxvx_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvx, gkyl_array_average_2x3v_hyb_p1_avgyvx_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvx, gkyl_array_average_2x3v_hyb_p1_avgxyvx_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgvy, gkyl_array_average_2x3v_hyb_p1_avgvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxvy, gkyl_array_average_2x3v_hyb_p1_avgxvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvy, gkyl_array_average_2x3v_hyb_p1_avgyvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvy, gkyl_array_average_2x3v_hyb_p1_avgxyvy_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgvxvy, gkyl_array_average_2x3v_hyb_p1_avgvxvy_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgxvxvy,
              gkyl_array_average_2x3v_hyb_p1_avgxvxvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvxvy,
              gkyl_array_average_2x3v_hyb_p1_avgyvxvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvxvy,
              gkyl_array_average_2x3v_hyb_p1_avgxyvxvy_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgvz, gkyl_array_average_2x3v_hyb_p1_avgvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxvz, gkyl_array_average_2x3v_hyb_p1_avgxvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvz, gkyl_array_average_2x3v_hyb_p1_avgyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvz, gkyl_array_average_2x3v_hyb_p1_avgxyvz_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgvxvz, gkyl_array_average_2x3v_hyb_p1_avgvxvz_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgxvxvz,
              gkyl_array_average_2x3v_hyb_p1_avgxvxvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvxvz,
              gkyl_array_average_2x3v_hyb_p1_avgyvxvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvxvz,
              gkyl_array_average_2x3v_hyb_p1_avgxyvxvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgvyvz, gkyl_array_average_2x3v_hyb_p1_avgvyvz_fulldim
             },
             {gkyl_array_average_2x3v_hyb_p1_avgxvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgxvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgyvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgxyvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgvxvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgvxvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxvxvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgxvxvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgyvxvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgyvxvyvz_fulldim},
             {gkyl_array_average_2x3v_hyb_p1_avgxyvxvyvz,
              gkyl_array_average_2x3v_hyb_p1_avgxyvxvyvz_fulldim}
           },
       }
      },
  },
  // 3x
  {
    .vdim =
      {// 3x1v (none)
       {.ops = {{NULL, NULL}}},
       // 3x2v (none)
       {.ops = {{NULL, NULL}}},
       // 3x3v
       {
         .ops =
           {{gkyl_array_average_3x3v_hyb_p1_avgx, gkyl_array_average_3x3v_hyb_p1_avgx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgy, gkyl_array_average_3x3v_hyb_p1_avgy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxy, gkyl_array_average_3x3v_hyb_p1_avgxy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgz, gkyl_array_average_3x3v_hyb_p1_avgz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxz, gkyl_array_average_3x3v_hyb_p1_avgxz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyz, gkyl_array_average_3x3v_hyb_p1_avgyz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyz, gkyl_array_average_3x3v_hyb_p1_avgxyz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvx, gkyl_array_average_3x3v_hyb_p1_avgvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvx, gkyl_array_average_3x3v_hyb_p1_avgxvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvx, gkyl_array_average_3x3v_hyb_p1_avgyvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvx, gkyl_array_average_3x3v_hyb_p1_avgxyvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvx, gkyl_array_average_3x3v_hyb_p1_avgzvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvx, gkyl_array_average_3x3v_hyb_p1_avgxzvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvx, gkyl_array_average_3x3v_hyb_p1_avgyzvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvx,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvx_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvy, gkyl_array_average_3x3v_hyb_p1_avgvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvy, gkyl_array_average_3x3v_hyb_p1_avgxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvy, gkyl_array_average_3x3v_hyb_p1_avgyvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvy, gkyl_array_average_3x3v_hyb_p1_avgxyvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvy, gkyl_array_average_3x3v_hyb_p1_avgzvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvy, gkyl_array_average_3x3v_hyb_p1_avgxzvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvy, gkyl_array_average_3x3v_hyb_p1_avgyzvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvy,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvxvy, gkyl_array_average_3x3v_hyb_p1_avgvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgxvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgyvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgxyvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgzvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgxzvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgyzvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvxvy,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvxvy_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvz, gkyl_array_average_3x3v_hyb_p1_avgvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvz, gkyl_array_average_3x3v_hyb_p1_avgxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvz, gkyl_array_average_3x3v_hyb_p1_avgyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvz, gkyl_array_average_3x3v_hyb_p1_avgxyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvz, gkyl_array_average_3x3v_hyb_p1_avgzvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvz, gkyl_array_average_3x3v_hyb_p1_avgxzvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvz, gkyl_array_average_3x3v_hyb_p1_avgyzvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvxvz, gkyl_array_average_3x3v_hyb_p1_avgvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgxvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgyvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgzvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgxzvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgyzvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvxvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvxvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvyvz, gkyl_array_average_3x3v_hyb_p1_avgvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgyvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgzvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxzvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgyzvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgyvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgzvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgzvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxzvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxzvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgyzvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgyzvxvyvz_fulldim},
            {gkyl_array_average_3x3v_hyb_p1_avgxyzvxvyvz,
             gkyl_array_average_3x3v_hyb_p1_avgxyzvxvyvz_fulldim}},
       }
      },
  }
};
// Phase-space serendipity (p=2) kernels, indexed [cdim-1].vdim[vdim-1].ops[op-1].
GKYL_CU_D static const array_average_v_kern_list gkyl_array_average_ker_list_ser_phase[] = {
  // 1x
  {
    .vdim =
      {// 1x1v
       {
         .ops =
           {{gkyl_array_average_1x1v_ser_p2_avgx, gkyl_array_average_1x1v_ser_p2_avgx_fulldim},
            {gkyl_array_average_1x1v_ser_p2_avgvx, gkyl_array_average_1x1v_ser_p2_avgvx_fulldim},
            {gkyl_array_average_1x1v_ser_p2_avgxvx, gkyl_array_average_1x1v_ser_p2_avgxvx_fulldim}},
       },
       // 1x2v
       {
         .ops =
           {{gkyl_array_average_1x2v_ser_p2_avgx, gkyl_array_average_1x2v_ser_p2_avgx_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgvx, gkyl_array_average_1x2v_ser_p2_avgvx_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgxvx, gkyl_array_average_1x2v_ser_p2_avgxvx_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgvy, gkyl_array_average_1x2v_ser_p2_avgvy_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgxvy, gkyl_array_average_1x2v_ser_p2_avgxvy_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgvxvy, gkyl_array_average_1x2v_ser_p2_avgvxvy_fulldim},
            {gkyl_array_average_1x2v_ser_p2_avgxvxvy,
             gkyl_array_average_1x2v_ser_p2_avgxvxvy_fulldim}},
       },
       // 1x3v
       {
         .ops =
           {{gkyl_array_average_1x3v_ser_p2_avgx, gkyl_array_average_1x3v_ser_p2_avgx_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvx, gkyl_array_average_1x3v_ser_p2_avgvx_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvx, gkyl_array_average_1x3v_ser_p2_avgxvx_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvy, gkyl_array_average_1x3v_ser_p2_avgvy_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvy, gkyl_array_average_1x3v_ser_p2_avgxvy_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvxvy, gkyl_array_average_1x3v_ser_p2_avgvxvy_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvxvy,
             gkyl_array_average_1x3v_ser_p2_avgxvxvy_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvz, gkyl_array_average_1x3v_ser_p2_avgvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvz, gkyl_array_average_1x3v_ser_p2_avgxvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvxvz, gkyl_array_average_1x3v_ser_p2_avgvxvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvxvz,
             gkyl_array_average_1x3v_ser_p2_avgxvxvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvyvz, gkyl_array_average_1x3v_ser_p2_avgvyvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvyvz,
             gkyl_array_average_1x3v_ser_p2_avgxvyvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgvxvyvz,
             gkyl_array_average_1x3v_ser_p2_avgvxvyvz_fulldim},
            {gkyl_array_average_1x3v_ser_p2_avgxvxvyvz,
             gkyl_array_average_1x3v_ser_p2_avgxvxvyvz_fulldim}},
       }
      },
  },
  // 2x
  {
    .vdim =
      {// 2x1v (none)
       {.ops = {{NULL, NULL}}},
       // 2x2v
       {
         .ops =
           {{gkyl_array_average_2x2v_ser_p2_avgx, gkyl_array_average_2x2v_ser_p2_avgx_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgy, gkyl_array_average_2x2v_ser_p2_avgy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxy, gkyl_array_average_2x2v_ser_p2_avgxy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgvx, gkyl_array_average_2x2v_ser_p2_avgvx_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxvx, gkyl_array_average_2x2v_ser_p2_avgxvx_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgyvx, gkyl_array_average_2x2v_ser_p2_avgyvx_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxyvx, gkyl_array_average_2x2v_ser_p2_avgxyvx_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgvy, gkyl_array_average_2x2v_ser_p2_avgvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxvy, gkyl_array_average_2x2v_ser_p2_avgxvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgyvy, gkyl_array_average_2x2v_ser_p2_avgyvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxyvy, gkyl_array_average_2x2v_ser_p2_avgxyvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgvxvy, gkyl_array_average_2x2v_ser_p2_avgvxvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxvxvy,
             gkyl_array_average_2x2v_ser_p2_avgxvxvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgyvxvy,
             gkyl_array_average_2x2v_ser_p2_avgyvxvy_fulldim},
            {gkyl_array_average_2x2v_ser_p2_avgxyvxvy,
             gkyl_array_average_2x2v_ser_p2_avgxyvxvy_fulldim}},
       },
       // 2x3v
       {
         .ops =
           {
             {gkyl_array_average_2x3v_ser_p2_avgx, gkyl_array_average_2x3v_ser_p2_avgx_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgy, gkyl_array_average_2x3v_ser_p2_avgy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxy, gkyl_array_average_2x3v_ser_p2_avgxy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgvx, gkyl_array_average_2x3v_ser_p2_avgvx_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxvx, gkyl_array_average_2x3v_ser_p2_avgxvx_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvx, gkyl_array_average_2x3v_ser_p2_avgyvx_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvx, gkyl_array_average_2x3v_ser_p2_avgxyvx_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgvy, gkyl_array_average_2x3v_ser_p2_avgvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxvy, gkyl_array_average_2x3v_ser_p2_avgxvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvy, gkyl_array_average_2x3v_ser_p2_avgyvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvy, gkyl_array_average_2x3v_ser_p2_avgxyvy_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgvxvy, gkyl_array_average_2x3v_ser_p2_avgvxvy_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgxvxvy,
              gkyl_array_average_2x3v_ser_p2_avgxvxvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvxvy,
              gkyl_array_average_2x3v_ser_p2_avgyvxvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvxvy,
              gkyl_array_average_2x3v_ser_p2_avgxyvxvy_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgvz, gkyl_array_average_2x3v_ser_p2_avgvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxvz, gkyl_array_average_2x3v_ser_p2_avgxvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvz, gkyl_array_average_2x3v_ser_p2_avgyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvz, gkyl_array_average_2x3v_ser_p2_avgxyvz_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgvxvz, gkyl_array_average_2x3v_ser_p2_avgvxvz_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgxvxvz,
              gkyl_array_average_2x3v_ser_p2_avgxvxvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvxvz,
              gkyl_array_average_2x3v_ser_p2_avgyvxvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvxvz,
              gkyl_array_average_2x3v_ser_p2_avgxyvxvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgvyvz, gkyl_array_average_2x3v_ser_p2_avgvyvz_fulldim
             },
             {gkyl_array_average_2x3v_ser_p2_avgxvyvz,
              gkyl_array_average_2x3v_ser_p2_avgxvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvyvz,
              gkyl_array_average_2x3v_ser_p2_avgyvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvyvz,
              gkyl_array_average_2x3v_ser_p2_avgxyvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgvxvyvz,
              gkyl_array_average_2x3v_ser_p2_avgvxvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxvxvyvz,
              gkyl_array_average_2x3v_ser_p2_avgxvxvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgyvxvyvz,
              gkyl_array_average_2x3v_ser_p2_avgyvxvyvz_fulldim},
             {gkyl_array_average_2x3v_ser_p2_avgxyvxvyvz,
              gkyl_array_average_2x3v_ser_p2_avgxyvxvyvz_fulldim}
           },
       }
      },
  }
};

// Primary struct in this updater.
struct gkyl_array_average {
  int ndim; // Dimensionality of the input array.
  int vdim; // Number of velocity-space dimensions (0 for a conf-space array).
  struct gkyl_basis basis; // Basis of the input array.
  struct gkyl_basis basis_avg; // Basis of the output array.
  struct gkyl_range local; // Range of the input array.
  struct gkyl_range local_avg; // Range of the output array.
  bool use_gpu; // Whether to compute on the GPU.
  bool full_dim_out; // Whether the output has the dimensionality of the input.

  int avg_dim[GKYL_MAX_DIM]; // 1 for averaged dimensions, 0 otherwise.
  int dim_remains[GKYL_MAX_DIM]; // 1 for dimensions that remain, 0 otherwise.
  int num_avg_dim; // Number of averaged dimensions.
  int num_dim_remain; // Number of remaining dimensions.

  double subvol; // Cell volume of the averaged dimensions (over the domain
  // volume of these dimensions if the average is unweighted).
  array_average_t kernel; // Single cell average kernel.
  struct gkyl_array_average *on_dev; // Pointer to itself on device.

  bool isweighted; // Whether the average is weighted.
  struct gkyl_array *weight; // Weight (constant if unweighted).
  struct gkyl_array *weight_avg; // Integral of the weight over the averaged dimensions.
  double vol_avg_inv; // Inverse domain volume of the averaged dimensions.
  gkyl_dg_bin_op_mem *div_mem; // Memory for the weak division by weight_avg.
};

GKYL_CU_D static void
array_average_choose_kernel(struct gkyl_array_average *up)
{
  int ndim = up->ndim, vdim = up->vdim, cdim = ndim - vdim;
  int poly_order = up->basis.poly_order;
  int ivar = up->full_dim_out ? 1 : 0;

  // The averaged dimensions are encoded in the bits of op (e.g. op=3 is avgxy, op=5 is avgxz).
  int op = 0;
  for (int d = 0; d < ndim; d++) {
    op += up->avg_dim[d] << d;
  }

  if (vdim == 0) {
    up->kernel =
      gkyl_array_average_ker_list_ser[ndim - 1].list[poly_order - 1].ops[op - 1].kernels[ivar];
  } else {
    switch (up->basis.b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        assert(poly_order == 2);
        up->kernel =
          gkyl_array_average_ker_list_ser_phase[cdim - 1].vdim[vdim - 1].ops[op - 1].kernels[ivar];
        break;
      case GKYL_BASIS_MODAL_GKHYBRID:
        up->kernel = gkyl_array_average_ker_list_gkhyb[ndim - 2].ops[op - 1].kernels[ivar];
        break;
      case GKYL_BASIS_MODAL_HYBRID:
        up->kernel =
          gkyl_array_average_ker_list_hyb[cdim - 1].vdim[vdim - 1].ops[op - 1].kernels[ivar];
        break;
      default:
        assert(false);
        break;
    }
  }
  assert(up->kernel);
}

// Linear index in the output array of the cell receiving the average of the input cell idx.
GKYL_CU_DH static inline long
array_average_out_lidx(const struct gkyl_array_average *up, const int *idx)
{
  int idx_out[GKYL_MAX_DIM];
  if (up->full_dim_out) {
    for (int d = 0; d < up->ndim; d++) {
      idx_out[d] = up->avg_dim[d] ? up->local_avg.lower[d] : idx[d];
    }
  } else if (up->num_dim_remain > 0) {
    int cnt = 0;
    for (int d = 0; d < up->ndim; d++) {
      if (up->dim_remains[d]) {
        idx_out[cnt++] = idx[d];
      }
    }
  } else {
    idx_out[0] = up->local_avg.lower[0];
  }
  return gkyl_range_idx(&up->local_avg, idx_out);
}

#ifdef GKYL_HAVE_CUDA
// Device new function.
struct gkyl_array_average *gkyl_array_average_cu_dev_new(struct gkyl_array_average *up);
// Device advance function.
void gkyl_array_average_advance_cu(
  const struct gkyl_array_average *up, const struct gkyl_array *fin, struct gkyl_array *avgout
);
#endif
