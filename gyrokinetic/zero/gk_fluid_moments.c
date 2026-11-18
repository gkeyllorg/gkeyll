#include <assert.h>
#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_fluid_moments_priv.h>

// Allocate a double array on the requested device.
static inline struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                   gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// The state is J*(rho,Mpar,E,pperp/B), with Mpar=rho*u and u the parallel velocity.
// Moment densities below include J; primitive u and T/m do not.
// Products and quotients of fields are evaluated with weak-DG operations.

// Get the weighted perpendicular pressure: J*pperp = B*state[3].
static void
gk_fluid_moments_calc_Jpperp_from_state(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state
)
{
  gkyl_dg_mul_op_range(&up->basis, 0, up->pperp, 0, up->bmag, 3, state, range);
}

// Construct the weighted perpendicular pressure: J*pperp = (J*rho)*(Tperp/m).
static void
gk_fluid_moments_calc_Jpperp_from_prim(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *prim
)
{
  gkyl_dg_mul_op_range(&up->basis, 0, up->pperp, 0, state, 3, prim, range);
}

// Calculate:
//   upar = (J*Mpar)/(J*rho),
//   J*pperp = B*state[3],
//   J*ppar = 2*J*E-upar*(J*Mpar)-2*J*pperp.
static void
gk_fluid_moments_calc_upar_pressures(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state
)
{
  gkyl_dg_div_op_range(up->mem, &up->basis, 0, up->upar, 1, state, 0, state, range);
  gk_fluid_moments_calc_Jpperp_from_state(up, range, state);
  gkyl_dg_mul_op_range(&up->basis, 0, up->ppar, 0, up->upar, 1, state, range);
  gkyl_array_scale_range(up->ppar, -1.0, range);
  gkyl_array_accumulate_offset_range(up->ppar, 2.0, state, 2 * up->basis.num_basis, range);
  gkyl_array_accumulate_range(up->ppar, -2.0, up->pperp, range);
}

// Compute M0, M1, or M2 by scaling the selected conserved field:
// result = rmass_fac*state[component].
static void
gk_fluid_moments_m0_or_m1_or_m2(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gkyl_array_set_offset_ncoeff_range(
    up->result, up->rmass_fac, state, nb, up->component * nb, 0, range
  );
}

// Compute the perpendicular second moment: M2perp = 2*J*pperp/m.
static void
gk_fluid_moments_m2perp(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gk_fluid_moments_calc_Jpperp_from_state(up, range, state);
  gkyl_array_set_offset_ncoeff_range(up->result, 2.0 / up->mass, up->pperp, nb, 0, 0, range);
}

// Compute the parallel second moment: M2par = 2*(J*E-J*pperp)/m.
static void
gk_fluid_moments_m2par(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gk_fluid_moments_calc_Jpperp_from_state(up, range, state);
  gkyl_array_set_offset_range(up->ppar, 2.0 / up->mass, state, 2 * up->basis.num_basis, range);
  gkyl_array_accumulate_range(up->ppar, -2.0 / up->mass, up->pperp, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, up->ppar, nb, 0, 0, range);
}

// Close the parallel third moment: M3par = upar*(M2par+2*J*ppar/m).
static void
gk_fluid_moments_m3par(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  gk_fluid_moments_calc_upar_pressures(up, range, state);
  gkyl_array_scale_range(up->ppar, 2.0 / up->mass, range);
  gkyl_array_accumulate_offset_range(
    up->ppar, 2.0 / up->mass, state, 2 * up->basis.num_basis, range
  );
  gkyl_array_accumulate_range(up->ppar, -2.0 / up->mass, up->pperp, range);
  gkyl_dg_mul_op_range(&up->basis, 0, up->result, 0, up->upar, 0, up->ppar, range);
}

// Close the perpendicular third moment: M3perp = upar*(2*J*pperp/m).
static void
gk_fluid_moments_m3perp(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  gkyl_dg_div_op_range(up->mem, &up->basis, 0, up->upar, 1, state, 0, state, range);
  gk_fluid_moments_calc_Jpperp_from_state(up, range, state);
  gkyl_array_scale_range(up->pperp, 2.0 / up->mass, range);
  gkyl_dg_mul_op_range(&up->basis, 0, up->result, 0, up->upar, 0, up->pperp, range);
}

// Compute:
//   J*n = J*rho/m,
//   upar = (J*Mpar)/(J*rho),
//   T/m = (J*ppar+2*J*pperp)/(3*J*rho).
static void
gk_fluid_moments_maxwellian(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gk_fluid_moments_calc_upar_pressures(up, range, state);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0 / up->mass, state, nb, 0, 0, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, up->upar, nb, 0, nb, range);
  gkyl_array_accumulate_range(up->ppar, 2.0, up->pperp, range);
  gkyl_array_scale_range(up->ppar, 1.0 / 3.0, range);
  gkyl_dg_div_op_range(up->mem, &up->basis, 2, up->result, 0, up->ppar, 0, state, range);
}

// Compute:
//   J*n = J*rho/m,
//   upar = (J*Mpar)/(J*rho),
//   Tpar/m = (J*ppar)/(J*rho),
//   Tperp/m = (J*pperp)/(J*rho).
static void
gk_fluid_moments_bimaxwellian(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gk_fluid_moments_calc_upar_pressures(up, range, state);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0 / up->mass, state, nb, 0, 0, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, up->upar, nb, 0, nb, range);
  gkyl_dg_div_op_range(up->mem, &up->basis, 2, up->result, 0, up->ppar, 0, state, range);
  gkyl_dg_div_op_range(up->mem, &up->basis, 3, up->result, 0, up->pperp, 0, state, range);
}

// Compute the combined moments: (M0,M1,M2) = (J*rho/m,J*Mpar/m,2*J*E/m).
static void
gk_fluid_moments_m0m1m2(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  for (int c = 0; c < 3; ++c) {
    gkyl_array_set_offset_ncoeff_range(
      up->result, (c == 2 ? 2.0 : 1.0) / up->mass, state, nb, c * nb, c * nb, range
    );
  }
}

// Compute (M0,M1,M2par,M2perp) =
// (J*rho/m,J*Mpar/m,2*(J*E-J*pperp)/m,2*J*pperp/m).
static void
gk_fluid_moments_m0m1m2parperp(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  gk_fluid_moments_m0m1m2(up, range, state, phi);
  gk_fluid_moments_calc_Jpperp_from_state(up, range, state);
  gkyl_array_set_offset_range(up->ppar, 2.0 / up->mass, state, 2 * up->basis.num_basis, range);
  gkyl_array_accumulate_range(up->ppar, -2.0 / up->mass, up->pperp, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, up->ppar, nb, 0, 2 * nb, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 2.0 / up->mass, up->pperp, nb, 0, 3 * nb, range);
}

// Compute density, parallel momentum, and Hamiltonian density:
// (J*n,J*Mpar,J*E+q*phi*(J*n)), with J*n = J*rho/m.
static void
gk_fluid_moments_hamiltonian(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi
)
{
  int nb = up->basis.num_basis;
  assert(phi);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0 / up->mass, state, nb, 0, 0, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, state, nb, nb, nb, range);
  gkyl_dg_mul_op_range(&up->basis, 0, up->ppar, 0, phi, 0, state, range);
  gkyl_array_scale_range(up->ppar, up->charge / up->mass, range);
  gkyl_array_accumulate_offset_range(up->ppar, 1.0, state, 2 * up->basis.num_basis, range);
  gkyl_array_set_offset_ncoeff_range(up->result, 1.0, up->ppar, nb, 0, 2 * nb, range);
}

// Copy the requested DG moment components from the result buffer to out.
static void
gk_fluid_moments_final_copy(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range, struct gkyl_array *out
)
{
  int nb = up->basis.num_basis;
  for (int c = 0; c < up->num_mom_out; ++c) {
    gkyl_array_set_offset_ncoeff_range(out, 1.0, up->result, nb, c * nb, c * nb, range);
  }
}

// Integrate each moment over a cell: integral[c] = product_d(dx[d]/sqrt(2))*q[c,0],
// where q[c,0] is the constant DG coefficient of moment c.
static void
gk_fluid_moments_final_integrate(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range, struct gkyl_array *out
)
{
#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out) == gkyl_array_is_cu_dev(up->result));
  if (gkyl_array_is_cu_dev(out)) {
    gkyl_gk_fluid_moments_final_integrate_cu(up, range, out);
    return;
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *dg = gkyl_array_cfetch(up->result, loc);
    double *integral = gkyl_array_fetch(out, loc);
    for (int c = 0; c < up->num_mom_out; ++c) {
      integral[c] = up->volnorm_fac * dg[c * up->basis.num_basis];
    }
  }
}

struct gkyl_gk_fluid_moments *
gkyl_gk_fluid_moments_inew(const struct gkyl_gk_fluid_moments_inp *inp)
{
  assert(inp->mass > 0 && inp->bmag && inp->num_moments >= 4);

  struct gkyl_gk_fluid_moments *up = gkyl_malloc(sizeof *up);

  *up = (struct gkyl_gk_fluid_moments){
    .basis = *inp->basis,
    .mass = inp->mass,
    .charge = inp->charge,
    .num_mom_out = 1,
    .num_mom_in = inp->num_moments,
    .num_cells = inp->range->volume,
    .rmass_fac = 1.0 / inp->mass,
    .volnorm_fac = 1.0,
    .final_func = inp->integrated ? gk_fluid_moments_final_integrate : gk_fluid_moments_final_copy,
  };

  switch (inp->moment) {
    case GKYL_F_MOMENT_M0:
    case GKYL_F_MOMENT_M1:
    case GKYL_F_MOMENT_M2:
      up->component = inp->moment == GKYL_F_MOMENT_M0 ? 0 : inp->moment == GKYL_F_MOMENT_M1 ? 1 : 2;
      up->rmass_fac *= inp->moment == GKYL_F_MOMENT_M2 ? 2.0 : 1.0;
      up->calc_func = gk_fluid_moments_m0_or_m1_or_m2;
      break;
    case GKYL_F_MOMENT_M2PAR:
      up->calc_func = gk_fluid_moments_m2par;
      break;
    case GKYL_F_MOMENT_M2PERP:
      up->calc_func = gk_fluid_moments_m2perp;
      break;
    case GKYL_F_MOMENT_M3PAR:
      up->calc_func = gk_fluid_moments_m3par;
      break;
    case GKYL_F_MOMENT_M3PERP:
      up->calc_func = gk_fluid_moments_m3perp;
      break;
    case GKYL_F_MOMENT_MAXWELLIAN:
      up->num_mom_out = 3;
      up->calc_func = gk_fluid_moments_maxwellian;
      break;
    case GKYL_F_MOMENT_BIMAXWELLIAN:
      up->num_mom_out = 4;
      up->calc_func = gk_fluid_moments_bimaxwellian;
      break;
    case GKYL_F_MOMENT_M0M1M2:
      up->num_mom_out = 3;
      up->calc_func = gk_fluid_moments_m0m1m2;
      break;
    case GKYL_F_MOMENT_M0M1M2PARM2PERP:
      up->num_mom_out = 4;
      up->calc_func = gk_fluid_moments_m0m1m2parperp;
      break;
    case GKYL_F_MOMENT_HAMILTONIAN:
      up->num_mom_out = 3;
      up->calc_func = gk_fluid_moments_hamiltonian;
      break;
    default:
      gkyl_exit("Unsupported gyrokinetic fluid moment.");
  }

  for (int d = 0; d < inp->grid->ndim; ++d) {
    up->volnorm_fac *= inp->grid->dx[d] / sqrt(2.0);
  }

  int nb = inp->basis->num_basis;

  up->bmag = gkyl_array_acquire(inp->bmag);
  up->upar = mkarr(inp->use_gpu, nb, inp->range_ext->volume);
  up->ppar = mkarr(inp->use_gpu, nb, inp->range_ext->volume);
  up->pperp = mkarr(inp->use_gpu, nb, inp->range_ext->volume);
  up->result = mkarr(inp->use_gpu, 4 * nb, inp->range_ext->volume);
  up->mem = inp->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(inp->range->volume, nb) :
                           gkyl_dg_bin_op_mem_new(inp->range->volume, nb);

  return up;
}

int
gkyl_gk_fluid_moments_num_mom(const struct gkyl_gk_fluid_moments *up)
{
  return up->num_mom_out;
}

void
gkyl_gk_fluid_moments_advance(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *state, const struct gkyl_array *phi, struct gkyl_array *out
)
{
  assert(range->volume == up->num_cells);
  up->calc_func(up, range, state, phi);
  up->final_func(up, range, out);
}

void
gkyl_gk_fluid_moments_from_prim(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *prim, const struct gkyl_array *jacobgeo, struct gkyl_array *state
)
{
  int nb = up->basis.num_basis;
  assert(range->volume == up->num_cells);
  // Weight density first so recovery inverts the same weak products.
  gkyl_dg_mul_op_range(&up->basis, 0, up->upar, 0, jacobgeo, 0, prim, range);
  gkyl_array_set_offset_ncoeff_range(state, up->mass, up->upar, nb, 0, 0, range);
  gkyl_dg_mul_op_range(&up->basis, 1, state, 0, state, 1, prim, range);
  gkyl_dg_mul_op_range(&up->basis, 0, up->ppar, 0, state, 2, prim, range);
  gk_fluid_moments_calc_Jpperp_from_prim(up, range, state, prim);
  gkyl_dg_mul_op_range(&up->basis, 0, up->upar, 1, prim, 1, state, range);
  gkyl_array_accumulate_range(up->upar, 1.0, up->ppar, range);
  gkyl_array_scale_range(up->upar, 0.5, range);
  gkyl_array_accumulate_range(up->upar, 1.0, up->pperp, range);
  gkyl_array_set_offset_ncoeff_range(state, 1.0, up->upar, nb, 0, 2 * nb, range);
  gkyl_dg_div_op_range(up->mem, &up->basis, 3, state, 0, up->pperp, 0, up->bmag, range);
}

void
gkyl_gk_fluid_moments_set_density(
  const struct gkyl_gk_fluid_moments *up, const struct gkyl_range *range,
  const struct gkyl_array *density, struct gkyl_array *state
)
{
  int nb = up->basis.num_basis;
  assert(range->volume == up->num_cells);
  assert(state->ncomp == up->num_mom_in * nb);

  gkyl_array_set_offset_range(up->ppar, up->mass, density, 0, range);
  gkyl_dg_div_op_range(up->mem, &up->basis, 0, up->upar, 0, up->ppar, 0, state, range);
  for (int c = 1; c < up->num_mom_in; ++c) {
    gkyl_dg_mul_op_range(&up->basis, c, state, 0, up->upar, c, state, range);
  }
  gkyl_array_set_offset_ncoeff_range(state, 1.0, up->ppar, nb, 0, 0, range);
}

void
gkyl_gk_fluid_moments_release(struct gkyl_gk_fluid_moments *up)
{
  gkyl_array_release(up->bmag);
  gkyl_array_release(up->upar);
  gkyl_array_release(up->ppar);
  gkyl_array_release(up->pperp);
  gkyl_array_release(up->result);
  gkyl_dg_bin_op_mem_release(up->mem);
  gkyl_free(up);
}
