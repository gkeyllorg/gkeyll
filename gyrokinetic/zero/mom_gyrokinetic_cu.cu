/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>
}

static int
gk_num_mom(int vdim, enum gkyl_distribution_moments mom_type)
{
  int num_mom = 0;
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M2:
    case GKYL_F_MOMENT_M2PAR:
    case GKYL_F_MOMENT_M2PERP:
    case GKYL_F_MOMENT_M3PAR:
    case GKYL_F_MOMENT_M3PERP:
      num_mom = 1;
      break;

    case GKYL_F_MOMENT_M0M1M2:
      num_mom = 3;
      break;      
      
    case GKYL_F_MOMENT_M0M1M2PARM2PERP:
      num_mom = vdim+2;
      break;      
      
    case GKYL_F_MOMENT_HAMILTONIAN:
      num_mom = 3;
      break;

    default: // can't happen
      fprintf(stderr,"Moment option %d not available.\n",mom_type);
      assert(false);
      break;
  }

  return num_mom;
}

__global__
static void
set_cu_ptrs(struct mom_type_gyrokinetic *mom_gk, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  // choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, 
    *m2_par_kernels, *m2_perp_kernels, *m3_par_kernels, *m3_perp_kernels,
    *three_moments_kernels, *four_moments_kernels, *hamiltonian_moments_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1_kernels = ser_m1_kernels;
      m2_kernels = ser_m2_kernels;
      m2_par_kernels = ser_m2_par_kernels;
      m2_perp_kernels = ser_m2_perp_kernels;
      m3_par_kernels = ser_m3_par_kernels;
      m3_perp_kernels = ser_m3_perp_kernels;
      three_moments_kernels = ser_three_moments_kernels;
      four_moments_kernels = ser_four_moments_kernels;
      hamiltonian_moments_kernels = ser_hamiltonian_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      mom_gk->momt.kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M1:
      mom_gk->momt.kernel = m1_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M2:
      mom_gk->momt.kernel = m2_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M2PAR:
      mom_gk->momt.kernel = m2_par_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M2PERP:
      mom_gk->momt.kernel = m2_perp_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M3PAR:
      mom_gk->momt.kernel = m3_par_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M3PERP:
      mom_gk->momt.kernel = m3_perp_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 1;
      break;

    case GKYL_F_MOMENT_M0M1M2:
      mom_gk->momt.kernel = three_moments_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 3;
      break;
      
    case GKYL_F_MOMENT_M0M1M2PARM2PERP:
      mom_gk->momt.kernel = four_moments_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = vdim+2;
      break;
      
    case GKYL_F_MOMENT_HAMILTONIAN:
      mom_gk->momt.kernel = hamiltonian_moments_kernels[tblidx].kernels[poly_order];
      mom_gk->momt.num_mom = 3;
      break;

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_gyrokinetic *mom_gk = (struct mom_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gk->momt.cdim = cdim;
  mom_gk->momt.pdim = pdim;
  mom_gk->momt.poly_order = poly_order;
  mom_gk->momt.num_config = cbasis->num_basis;
  mom_gk->momt.num_phase = pbasis->num_basis;

  mom_gk->momt.num_mom = gk_num_mom(vdim, mom_type); // number of moments

  mom_gk->mass = mass;
  mom_gk->charge = charge;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  mom_gk->gk_geom = geom_ho->on_dev;
  mom_gk->vel_map = vel_map_ho->on_dev;
  mom_gk->phi = 0;
  struct gkyl_array *phi_ho = 0;
  if (phi) {
    phi_ho = gkyl_array_acquire(phi);
    mom_gk->phi = phi_ho->on_dev;
  }

  mom_gk->conf_range = *conf_range;

  mom_gk->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_gk->momt.flags);
  mom_gk->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // Copy struct to device.
  struct mom_type_gyrokinetic *mom_gk_cu = (struct mom_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_gyrokinetic));
  gkyl_cu_memcpy(mom_gk_cu, mom_gk, sizeof(struct mom_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(mom_gk_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_gk->momt.on_dev = &mom_gk_cu->momt;

  // Updater should store host pointers.
  mom_gk->gk_geom = geom_ho; 
  mom_gk->vel_map = vel_map_ho; 
  mom_gk->phi = phi_ho; 
  
  return &mom_gk->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_gyrokinetic* momt, enum gkyl_distribution_moments mom_type,
  enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  // Choose kernel tables based on basis-function type.
  const gkyl_gyrokinetic_mom_kern_list *int_m0_kernels, *int_m1_kernels, *int_m2_par_kernel, *int_m2_perp_kernel, *int_m2_kernels,
    *int_m3_par_kernels, *int_m3_perp_kernels, *int_three_moments_kernels, *int_four_moments_kernels, *int_hamiltonian_moments_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_m0_kernels = ser_int_m0_kernels;
      int_m1_kernels = ser_int_m1_kernels;
      int_m2_par_kernel = ser_int_m2_par_kernels;
      int_m2_perp_kernel = ser_int_m2_perp_kernels;
      int_m2_kernels = ser_int_m2_kernels;
      int_m3_par_kernels = ser_int_m3_par_kernels;
      int_m3_perp_kernels = ser_int_m3_perp_kernels;
      int_three_moments_kernels = ser_int_three_moments_kernels;
      int_four_moments_kernels = ser_int_four_moments_kernels;
      int_hamiltonian_moments_kernels = ser_int_hamiltonian_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_type) {
    case GKYL_F_MOMENT_M0:
      momt->momt.kernel = int_m0_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M1:
      momt->momt.kernel = int_m1_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M2PAR:
      momt->momt.kernel = int_m2_par_kernel[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M2PERP:
      momt->momt.kernel = int_m2_perp_kernel[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M2:
      momt->momt.kernel = int_m2_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M3PAR:
      momt->momt.kernel = int_m3_par_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M3PERP:
      momt->momt.kernel = int_m3_perp_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;
    case GKYL_F_MOMENT_M0M1M2:
      momt->momt.kernel = int_three_moments_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 3;
      break;
      
    case GKYL_F_MOMENT_M0M1M2PARM2PERP:
      momt->momt.kernel = int_four_moments_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim+2;
      break;
      
    case GKYL_F_MOMENT_HAMILTONIAN:
      momt->momt.kernel = int_hamiltonian_moments_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 3;
      break;

    default: // Can't happen.
      assert(false);
      break;
  }
}

struct gkyl_mom_type*
gkyl_int_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, double charge, const struct gkyl_velocity_map* vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_gyrokinetic *momt = (struct mom_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->momt.num_mom = gk_num_mom(vdim, mom_type); // Number of moments.

  momt->mass = mass;
  momt->charge = charge;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  momt->gk_geom = geom_ho->on_dev;
  momt->vel_map = vel_map_ho->on_dev;
  momt->phi = 0;
  struct gkyl_array *phi_ho = 0;
  if (phi) {
    phi_ho = gkyl_array_acquire(phi);
    momt->phi = phi_ho->on_dev;
  }

  momt->conf_range = *conf_range;

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // Copy struct to device.
  struct mom_type_gyrokinetic *momt_cu = (struct mom_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_gyrokinetic));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, mom_type, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;

  // Updater should store host pointers.
  momt->gk_geom = geom_ho; 
  momt->vel_map = vel_map_ho; 
  momt->phi = phi_ho; 
  
  return &momt->momt;
}
