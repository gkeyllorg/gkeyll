// Test the GK neutral-fluid primitive-variable updater.
//
// The fluid moments are (rho, rho*ux, rho*uy, rho*uz, totalE).  For spatially
// uniform (constant) moments the primitive variables are exact algebraic
// combinations:
//   udrift_i = (rho*u_i)/rho
//   p        = (gas_gamma-1)*(E - 1/2 rho u^2)
//   T        = p/rho * mass            [ = (gas_gamma-1)*(mass*E - 1/2 (rho u)^2)/rho ]
//   flowE    = 1/2 rho u^2
//   thermalE = p/(gas_gamma-1) = E - 1/2 rho u^2
// Using is_integrated=true the updater returns the cell integral of each
// primitive variable, i.e. (value)*(cell volume), which we check exactly.
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_gk_neut_fluid_prim_vars.h>
#include <math.h>

// Exported by the library but (in this build) not declared in the public header.
void gkyl_gk_neut_fluid_prim_vars_flow_energy_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);

// Build a constant DG field of ncomp moment-components on a 1x grid; each
// component is filled with the constant value vals[c].
static struct gkyl_array*
mk_const_moms(struct gkyl_basis *basis, struct gkyl_range *range, int nmom, const double *vals)
{
  int nb = basis->num_basis;
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, nmom*nb, range->volume);
  gkyl_array_clear(arr, 0.0);

  // For a constant function value v, only the 0th DG coefficient is nonzero
  // and equals v*sqrt(2)^cdim.
  double fac = pow(sqrt(2.0), basis->ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    double *d = gkyl_array_fetch(arr, lidx);
    for (int c=0; c<nmom; c++)
      d[c*nb + 0] = vals[c]*fac;
  }
  return arr;
}

struct setup {
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  double cell_vol;
};

static void
make_setup(struct setup *s, int poly_order)
{
  double lower[] = {0.0}, upper[] = {2.0};
  int cells[] = {4};
  int dim = 1;
  gkyl_rect_grid_init(&s->grid, dim, lower, upper, cells);
  int ghost[] = {1};
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  gkyl_cart_modal_serendip(&s->basis, dim, poly_order);
  s->cell_vol = (upper[0]-lower[0])/cells[0];
}

// Read the integrated (cell-integral) scalar from component c of out, at the
// first cell of the local range, and divide out the cell volume to recover the
// pointwise primitive-variable value.
static double
read_val(struct gkyl_array *out, struct gkyl_range *range, int c, double cell_vol)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  gkyl_range_iter_next(&iter);
  long lidx = gkyl_range_idx(range, iter.idx);
  const double *d = gkyl_array_cfetch(out, lidx);
  return d[c]/cell_vol;
}

void
test_prim_vars_udrift_pressure_temp()
{
  struct setup s;
  make_setup(&s, 1);

  double gas_gamma = 5.0/3.0;
  double mass = 2.0;

  // Choose moments.
  double rho = 3.0, ux = 1.5, uy = -0.5, uz = 0.25;
  double rhoux = rho*ux, rhouy = rho*uy, rhouz = rho*uz;
  double usq = ux*ux + uy*uy + uz*uz;
  double p = 4.0; // desired pressure
  // E = p/(gas_gamma-1) + 1/2 rho u^2.
  double E = p/(gas_gamma-1.0) + 0.5*rho*usq;
  double moms[] = {rho, rhoux, rhouy, rhouz, E};

  struct gkyl_array *m = mk_const_moms(&s.basis, &s.local_ext, 5, moms);

  double T = p/rho*mass;             // T = mass*p/rho per kernel definition
  double flowE = 0.5*rho*usq;

  // --- udrift ---
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 3, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_udrift_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), ux, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 1, s.cell_vol), uy, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 2, s.cell_vol), uz, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  // --- pressure ---
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 1, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_pressure_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), p, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  // --- temperature ---
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_TEMP, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 1, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_temp_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), T, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  // --- flow energy ---
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_FLOW_ENERGY, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 1, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_flow_energy_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), flowE, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  gkyl_array_release(m);
}

void
test_prim_vars_combined()
{
  struct setup s;
  make_setup(&s, 1);

  double gas_gamma = 1.4;
  double mass = 1.0;

  double rho = 2.0, ux = 0.5, uy = 1.0, uz = -1.5;
  double usq = ux*ux + uy*uy + uz*uz;
  double p = 6.0;
  double E = p/(gas_gamma-1.0) + 0.5*rho*usq;
  double moms[] = {rho, rho*ux, rho*uy, rho*uz, E};
  double T = p/rho*mass;

  struct gkyl_array *m = mk_const_moms(&s.basis, &s.local_ext, 5, moms);

  // udrift + pressure -> (ux, uy, uz, p)
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_PRESSURE, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 4, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_udrift_pressure_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), ux, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 1, s.cell_vol), uy, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 2, s.cell_vol), uz, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 3, s.cell_vol), p, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  // udrift + temperature -> (ux, uy, uz, T)
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_TEMP, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 4, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_udrift_temp_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), ux, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 1, s.cell_vol), uy, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 2, s.cell_vol), uz, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 3, s.cell_vol), T, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  // LTE -> (n=rho/mass, ux, uy, uz, T/mass)
  {
    struct gkyl_gk_neut_fluid_prim_vars *up = gkyl_gk_neut_fluid_prim_vars_new(
      gas_gamma, mass, &s.basis, &s.grid, &s.local,
      GKYL_GK_NEUT_FLUID_PRIM_VARS_LTE, true, false);
    struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE, 5, s.local_ext.volume);
    gkyl_array_clear(out, 0.0);
    gkyl_gk_neut_fluid_prim_vars_lte_advance(up, m, out, 0);
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 0, s.cell_vol), rho/mass, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 1, s.cell_vol), ux, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 2, s.cell_vol), uy, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 3, s.cell_vol), uz, 1e-12) );
    TEST_CHECK( gkyl_compare(read_val(out, &s.local, 4, s.cell_vol), T/mass, 1e-12) );
    gkyl_array_release(out);
    gkyl_gk_neut_fluid_prim_vars_release(up);
  }

  gkyl_array_release(m);
}

TEST_LIST = {
  { "prim_vars_udrift_pressure_temp", test_prim_vars_udrift_pressure_temp },
  { "prim_vars_combined", test_prim_vars_combined },
  { NULL, NULL },
};
