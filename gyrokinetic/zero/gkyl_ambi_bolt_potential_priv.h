#pragma once

// Private header for ambi_bolt_potential updater, not for direct use in user code.

#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_ambi_bolt_potential_kernels.h>
#include <assert.h>

// Function pointer type for sheath entrance calculations.
typedef void (*sheathker_t)(const double sheathDirDx, double q_e, double m_e, double T_e,
  const double *cmag, const double *jacobtotInv,
  const double *GammaJac_i, const double *m0Ion, const double *m0JacIon, double *out);

// Function pointer type for phi calculation.
typedef void (*phiker_t)(double q_e, double T_e,
  const double *m0JacIon, const double *sheathvals, double *phi);

typedef struct { sheathker_t kernels[2]; } sheath_calc_kern_loc_list;
typedef struct { sheath_calc_kern_loc_list list[2]; } sheath_calc_kern_edge_list;

typedef struct { phiker_t kernels[2]; } phi_calc_kern_list;

// Serendipity sheath_calc kernels.
GKYL_CU_D
static const sheath_calc_kern_edge_list ser_sheath_calc_list[] = {
  // 1x
  { .list = {{ambi_bolt_potential_sheath_calc_lower_1x_ser_p1, ambi_bolt_potential_sheath_calc_upper_1x_ser_p1},
             {ambi_bolt_potential_sheath_calc_lower_1x_ser_p2, ambi_bolt_potential_sheath_calc_upper_1x_ser_p2}}, },
  // 2x
  { .list = {{ambi_bolt_potential_sheath_calc_lower_2x_ser_p1, ambi_bolt_potential_sheath_calc_upper_2x_ser_p1},
             {ambi_bolt_potential_sheath_calc_lower_2x_ser_p2, ambi_bolt_potential_sheath_calc_upper_2x_ser_p2}}, },
//  // 3x
  { .list = {{ambi_bolt_potential_sheath_calc_lower_3x_ser_p1, ambi_bolt_potential_sheath_calc_upper_3x_ser_p1},
             {ambi_bolt_potential_sheath_calc_lower_3x_ser_p2, ambi_bolt_potential_sheath_calc_upper_3x_ser_p2}}, },
};

// Serendipity phi_calc kernels.
GKYL_CU_D
static const phi_calc_kern_list ser_phi_calc_list[] = {
  // 1x kernels
  { ambi_bolt_potential_phi_calc_1x_ser_p1, ambi_bolt_potential_phi_calc_1x_ser_p2 },
  // 2x kernels
  { ambi_bolt_potential_phi_calc_2x_ser_p1, ambi_bolt_potential_phi_calc_2x_ser_p2 },
  // 3x kernels
  { ambi_bolt_potential_phi_calc_3x_ser_p1, ambi_bolt_potential_phi_calc_3x_ser_p2 },
};

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_ambi_bolt_potential_kernels {
  // Pointer to kernel that computes ion density and potential at lower/upper
  // sheath entrance.
  sheathker_t sheath_calc[2];

  // Pointer to kernel that computes the potential everywhere.
  phiker_t phi_calc;
};

// Primary struct in this updater.
struct gkyl_ambi_bolt_potential {
  int cdim;
  double num_basis;
  bool use_gpu;
  double dz;
  double mass_e;  // Electron mass.
  double charge_e;  // Electron charge.
  double temp_e;  // Electron temperature.
  struct gkyl_ambi_bolt_potential_kernels *kernels;  // sheath_calc and phi_calc kernels.
  struct gkyl_ambi_bolt_potential_kernels *kernels_cu;  // device copy.
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CSHEATHK(lst,dim,poly_order,loc) lst[dim-1].list[poly_order-1].kernels[loc]
#define CPHIK(lst,dim,poly_order) lst[dim-1].kernels[poly_order-1]

GKYL_CU_D
static void
ambi_bolt_potential_choose_kernels(const struct gkyl_basis* basis, struct gkyl_ambi_bolt_potential_kernels *kers)
{
  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      for (int k=0; k<2; k++) {
        kers->sheath_calc[k] = CSHEATHK(ser_sheath_calc_list, dim, poly_order, k);
      }
      kers->phi_calc = CPHIK(ser_phi_calc_list, dim, poly_order);
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
      break;
  }
}

#ifdef GKYL_HAVE_CUDA
void
ambi_bolt_potential_choose_kernels_cu(const struct gkyl_basis *basis, struct gkyl_ambi_bolt_potential_kernels *kers);

void
gkyl_ambi_bolt_potential_sheath_calc_cu(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_array *cmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, const struct gkyl_array *Jm0i,
  struct gkyl_array *sheath_vals);

void
gkyl_ambi_bolt_potential_phi_calc_cu(struct gkyl_ambi_bolt_potential *up,
  const struct gkyl_range *local_r, const struct gkyl_range *extlocal_r,
  const struct gkyl_array *m0i, const struct gkyl_array *sheath_vals,
  struct gkyl_array *phi);
#endif
