#pragma once

// Private header, not for direct use in user code
#include <gkyl_fpo_vlasov_kernels.h>

// Types for various kernels
typedef double (*fpo_vlasov_diff_surf_t)(const double *dxv,
  const double *diff_coeff_C, const double* diff_coeff_surf_stencil[9], 
  const double* f_stencil[9], double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_fpo_vlasov_diff_vol_kern_list;
typedef struct { fpo_vlasov_diff_surf_t kernels[9]; } fpo_vlasov_diff_surf_kern_list;
typedef struct { fpo_vlasov_diff_surf_kern_list list[3]; } fpo_vlasov_diff_surf_stencil_list;

struct dg_fpo_vlasov_diff {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  int upper_cells[GKYL_MAX_DIM];
  fpo_vlasov_diff_surf_kern_list surf[3][3]; // Generic surface term domain stencil
  struct gkyl_range phase_range; // Configuration space range.
  struct gkyl_dg_fpo_vlasov_diff_auxfields auxfields; // Auxiliary fields.
};

// Volume kernels
GKYL_CU_DH
static double
kernel_fpo_vlasov_diff_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);

  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idx);
  
  return fpo_vlasov_diff_vol_1x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.diff_coeff, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_diff_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);

  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idx);
  
  return fpo_vlasov_diff_vol_1x3v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.diff_coeff, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_diff_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);

  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idx);
  
  return fpo_vlasov_diff_vol_2x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.diff_coeff, pidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_vol_kern_list ser_vol_kernels[] = {
  { kernel_fpo_vlasov_diff_vol_1x3v_ser_p1, kernel_fpo_vlasov_diff_vol_1x3v_ser_p2, NULL },
  { kernel_fpo_vlasov_diff_vol_2x3v_ser_p1, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernels are handled generally using the domain decomp stencil.
// These lists contain both surface and boundary surface
// kernels, and the proper kernel is selected with the index returned from idx_to_inloup_ker 
// in hyper_dg. Diagonal terms only have in/lo/up as options, so keri is either 0/1/2

// Diagonal element kernels
// vx-vx direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vxvx_kernels[] = {
  {{
    { fpo_vlasov_diff_surfvxvx_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_lovx, fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_upvx, NULL, NULL, NULL, NULL, NULL, NULL },
    { fpo_vlasov_diff_surfvxvx_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p2_lovx, fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p2_upvx, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
  }},
  {{
    { fpo_vlasov_diff_surfvxvx_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvx_2x3v_ser_p1_lovx, fpo_vlasov_diff_boundary_surfvxvx_2x3v_ser_p1_upvx, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
   }}
};

// vy-vy direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vyvy_kernels[] = {
  {{
    { fpo_vlasov_diff_surfvyvy_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p1_lovy, fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p1_upvy, NULL, NULL, NULL, NULL, NULL, NULL },
    { fpo_vlasov_diff_surfvyvy_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p2_lovy, fpo_vlasov_diff_boundary_surfvyvy_1x3v_ser_p2_upvy, NULL, NULL, NULL, NULL, NULL, NULL },

    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
  }},
  {{
    { fpo_vlasov_diff_surfvyvy_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvy_2x3v_ser_p1_lovy, fpo_vlasov_diff_boundary_surfvyvy_2x3v_ser_p1_upvy, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
   }}
};

// vz-vz direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vzvz_kernels[] = {
  {{
    { fpo_vlasov_diff_surfvzvz_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p1_lovz, fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p1_upvz, NULL, NULL, NULL, NULL, NULL, NULL },
    { fpo_vlasov_diff_surfvzvz_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p2_lovz, fpo_vlasov_diff_boundary_surfvzvz_1x3v_ser_p2_upvz, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
  }},
  {{
    { fpo_vlasov_diff_surfvzvz_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvz_2x3v_ser_p1_lovz, fpo_vlasov_diff_boundary_surfvzvz_2x3v_ser_p1_upvz, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
   }}
};

// Off-diagonal element kernels
// vx-vy direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vxvy_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvxvy_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_invy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_invy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_invx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_invx_upvy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_lovx_upvy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_upvx_upvy },
    { fpo_vlasov_diff_surfvxvy_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_invy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_invy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_invx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_invx_upvy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_lovx_upvy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_lovy, fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p2_upvx_upvy },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvxvy_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_invy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_invy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_invx_lovy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_invx_upvy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_lovy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_lovx_upvy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_lovy, fpo_vlasov_diff_boundary_surfvxvy_2x3v_ser_p1_upvx_upvy },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

// vx-vz direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vxvz_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvxvz_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_invz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_invz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_invx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_invx_upvz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_lovx_upvz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p1_upvx_upvz },
    { fpo_vlasov_diff_surfvxvz_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_invz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_invz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_invx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_invx_upvz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_lovx_upvz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_lovz, fpo_vlasov_diff_boundary_surfvxvz_1x3v_ser_p2_upvx_upvz },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvxvz_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_invz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_invz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_invx_lovz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_invx_upvz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_lovz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_lovx_upvz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_lovz, fpo_vlasov_diff_boundary_surfvxvz_2x3v_ser_p1_upvx_upvz },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

// vy-vx direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vyvx_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvyvx_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_invy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_invy_upvx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_invx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_invx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_lovy_upvx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p1_upvy_upvx },
    { fpo_vlasov_diff_surfvyvx_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_invy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_invy_upvx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_invx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_invx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_lovx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_lovy_upvx, fpo_vlasov_diff_boundary_surfvyvx_1x3v_ser_p2_upvy_upvx },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvyvx_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_invy_lovx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_invy_upvx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_invx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_invx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_lovx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_lovx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_lovy_upvx, fpo_vlasov_diff_boundary_surfvyvx_2x3v_ser_p1_upvy_upvx },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

// vy-vz direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vyvz_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvyvz_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_invz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_invz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_invy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_invy_upvz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_lovy_upvz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p1_upvy_upvz },
    { fpo_vlasov_diff_surfvyvz_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_invz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_invz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_invy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_invy_upvz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_lovy_upvz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_lovz, fpo_vlasov_diff_boundary_surfvyvz_1x3v_ser_p2_upvy_upvz },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvyvz_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_invz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_invz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_invy_lovz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_invy_upvz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_lovz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_lovy_upvz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_lovz, fpo_vlasov_diff_boundary_surfvyvz_2x3v_ser_p1_upvy_upvz },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

// vz-vx direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vzvx_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvzvx_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_invz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_invz_upvx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_invx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_invx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_lovz_upvx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p1_upvz_upvx },
    { fpo_vlasov_diff_surfvzvx_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_invz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_invz_upvx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_invx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_invx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_lovx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_lovz_upvx, fpo_vlasov_diff_boundary_surfvzvx_1x3v_ser_p2_upvz_upvx },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvzvx_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_invz_lovx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_invz_upvx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_invx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_invx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_lovx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_lovx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_lovz_upvx, fpo_vlasov_diff_boundary_surfvzvx_2x3v_ser_p1_upvz_upvx },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

// vz-vy direction surface kernels
GKYL_CU_D
static const fpo_vlasov_diff_surf_stencil_list ser_surf_vzvy_kernels[] = {
  {{ // 1x3v
    { fpo_vlasov_diff_surfvzvy_1x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_invz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_invz_upvy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_invy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_invy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_lovz_upvy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p1_upvz_upvy },
    { fpo_vlasov_diff_surfvzvy_1x3v_ser_p2, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_invz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_invz_upvy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_invy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_invy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_lovy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_lovz_upvy, fpo_vlasov_diff_boundary_surfvzvy_1x3v_ser_p2_upvz_upvy },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
  }},
  {{ // 2x3v
    { fpo_vlasov_diff_surfvzvy_2x3v_ser_p1, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_invz_lovy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_invz_upvy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_invy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_invy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_lovy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_lovy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_lovz_upvy, fpo_vlasov_diff_boundary_surfvzvy_2x3v_ser_p1_upvz_upvy },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
   }}
};

/**
 * Free fpo_vlasov_diff equation object
 *
 * @param ref Reference counter for constant fpo_vlasov_diff equation
 */
void gkyl_fpo_vlasov_diff_free(const struct gkyl_ref_count* ref);

// Gen surface term called by hyper_dg_gen_stencil_advance
GKYL_CU_D
static double
fpo_diff_gen_surf_term(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  int keri, const int idx[9][GKYL_MAX_DIM], const double* qIn[9],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);

  const double* diff_coeff_surf_stencil[9];
  for (int i=0; i<9; ++i) {
    long lin_offset = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idx[i]);
    diff_coeff_surf_stencil[i] = (const double*)gkyl_array_cfetch(
      fpo_vlasov_diff->auxfields.diff_coeff_surf,
      lin_offset
    );
  }

  long linc = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idxc);
  int cdim = fpo_vlasov_diff->cdim;
  if (dir1 >= cdim && dir2 >= cdim) {
    fpo_vlasov_diff_surf_t *surf_kern_list = fpo_vlasov_diff->surf[dir1-cdim][dir2-cdim].kernels;
    return surf_kern_list[keri](dxc, 
      (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.diff_coeff, linc),
      diff_coeff_surf_stencil,
      qIn, qRhsOut);
  }
  return 0.0;
}
