#include <gkyl_moment_non_ideal_priv.h>

typedef void (*heat_flux_calc_t)(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double dt, double *q);

typedef void (*heat_flux_update_t)(const gkyl_ten_moment_grad_closure *gces,
  const double *q[], double *rhs);

struct gkyl_ten_moment_grad_closure {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  double k0; // damping coefficient
  double cfl; // CFL number to use
  double *cfla;
  bool use_gpu; // Boolean to determine whether wave equation object is on host or device
  
  struct gkyl_comm *comm;
  long offsets_centers[8];
  long offsets_vertices[8];

  heat_flux_calc_t calc_q;
  heat_flux_update_t update_q;

  struct gkyl_ten_moment_grad_closure *on_dev;
};

// Makes indexing cleaner
static const unsigned T11 = 0;
static const unsigned T12 = 1;
static const unsigned T13 = 2;
static const unsigned T22 = 3;
static const unsigned T23 = 4;
static const unsigned T33 = 5;

static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

// 1D stencil locations (L: lower, U: upper)
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper)
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

// 3D stencil locations (L: lower, U: upper)
enum loc_3d {
  LLL_3D, LLU_3D,
  LUL_3D, LUU_3D,
  ULL_3D, ULU_3D,
  UUL_3D, UUU_3D
};

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  int arr1[3] = { -1, -1, -1 }, arr2[3] = { 0, 0, 0 };
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, arr1, arr2);

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  int arr1[3] = { 0, 0, 0 }, arr2[3] = { 1, 1, 1 };
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, arr1, arr2);

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

GKYL_CU_D
static void
var_setup(const gkyl_ten_moment_grad_closure *gces,
  int start, int end,
  const double *fluid_d[],
  double rho[], double p[], double Tij[][6])
{
  for (int j = start; j <= end; ++j) {
    rho[j] = fluid_d[j][RHO];
    p[j] = (fluid_d[j][P11] - fluid_d[j][MX]*fluid_d[j][MX]/fluid_d[j][RHO]
          + fluid_d[j][P22] - fluid_d[j][MY]*fluid_d[j][MY]/fluid_d[j][RHO]
          + fluid_d[j][P33] - fluid_d[j][MZ]*fluid_d[j][MZ]/fluid_d[j][RHO])/3.0;
    Tij[j][T11] = (fluid_d[j][P11] - fluid_d[j][MX]*fluid_d[j][MX]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T12] = (fluid_d[j][P12] - fluid_d[j][MX]*fluid_d[j][MY]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T13] = (fluid_d[j][P13] - fluid_d[j][MX]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T22] = (fluid_d[j][P22] - fluid_d[j][MY]*fluid_d[j][MY]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T23] = (fluid_d[j][P23] - fluid_d[j][MY]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
    Tij[j][T33] = (fluid_d[j][P33] - fluid_d[j][MZ]*fluid_d[j][MZ]/fluid_d[j][RHO])/fluid_d[j][RHO];
  }
}

GKYL_CU_D
static inline double
calc_sym_grad_limiter_2D(double alpha, double a, double b)
{
  double avg = (a + b)/2;
  double min = fmin(alpha*a, a/alpha);
  double max = fmax(alpha*a, a/alpha);
  if (avg <= min) {
    return min;
  } else if (avg >= max) {
    return max;
  } else {
    return avg;
  }
}

GKYL_CU_D
static inline double
calc_sym_grad_limiter_3D(double alpha, double a, double b, double c, double d)
{
  double avg = (a + b + c + d)/4;
  double min = fmin(alpha*a, a/alpha);
  double max = fmax(alpha*a, a/alpha);
  if (avg <= min) {
    return min;
  } else if (avg >= max) {
    return max;
  } else {
    return avg;
  }
}

GKYL_CU_D
static void
calc_unmag_heat_flux_1d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double dt, double *q)
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;

  const double dx = gces->grid.dx[0];
  double dTdx[6] = {0.0};
  double dTdy[6] = {0.0};
  double dTdz[6] = {0.0};
  double Tij[2][6] = {0.0};
  double rho[2] = {0.0};
  double p[2] = {0.0};
  var_setup(gces, L_1D, U_1D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
  p_avg = calc_harmonic_avg_1D(p[L_1D], p[U_1D]);

  dTdx[T11] = calc_sym_grad_1D(dx, Tij[L_1D][T11], Tij[U_1D][T11]);
  dTdx[T12] = calc_sym_grad_1D(dx, Tij[L_1D][T12], Tij[U_1D][T12]);
  dTdx[T13] = calc_sym_grad_1D(dx, Tij[L_1D][T13], Tij[U_1D][T13]);
  dTdx[T22] = calc_sym_grad_1D(dx, Tij[L_1D][T22], Tij[U_1D][T22]);
  dTdx[T23] = calc_sym_grad_1D(dx, Tij[L_1D][T23], Tij[U_1D][T23]);
  dTdx[T33] = calc_sym_grad_1D(dx, Tij[L_1D][T33], Tij[U_1D][T33]);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);
  
  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  q[Q111] = chi*dTdx[T11];
  q[Q112] = chi*2.0*dTdx[T12]/3.0;
  q[Q113] = chi*2.0*dTdx[T13]/3.0;
  q[Q122] = chi*dTdx[T22]/3.0;
  q[Q123] = chi*dTdx[T23]/3.0;
  q[Q133] = chi*dTdx[T33]/3.0;

  double cfla = dt/(dx*dx);
  cflrate[0] = alpha*vth_avg*cfla;
}

GKYL_CU_D
static void
grad_closure_update_1d(const gkyl_ten_moment_grad_closure *gces,
  const double *q[], double *rhs)
{
  double div_qx[6] = {0.0};

  const double dx = gces->grid.dx[0];

  div_qx[0] = calc_sym_grad_1D(dx, q[L_1D][Q111], q[U_1D][Q111]);
  div_qx[1] = calc_sym_grad_1D(dx, q[L_1D][Q112], q[U_1D][Q112]);
  div_qx[2] = calc_sym_grad_1D(dx, q[L_1D][Q113], q[U_1D][Q113]);
  div_qx[3] = calc_sym_grad_1D(dx, q[L_1D][Q122], q[U_1D][Q122]);
  div_qx[4] = calc_sym_grad_1D(dx, q[L_1D][Q123], q[U_1D][Q123]);
  div_qx[5] = calc_sym_grad_1D(dx, q[L_1D][Q133], q[U_1D][Q133]);

  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = div_qx[0];
  rhs[P12] = div_qx[1];
  rhs[P13] = div_qx[2];
  rhs[P22] = div_qx[3];
  rhs[P23] = div_qx[4];
  rhs[P33] = div_qx[5];
}

GKYL_CU_D
static void
calc_unmag_heat_flux_2d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double dt, double *q)
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double limit = 0.75;

  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  double dTx[2][6] = {0.0};
  double dTy[2][6] = {0.0};
  double dTdx[2][6] = {0.0};
  double dTdy[2][6] = {0.0};
  double Tij[4][6] = {0.0};
  double rho[4] = {0.0};
  double p[4] = {0.0};
  var_setup(gces, LL_2D, UU_2D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
  p_avg = calc_harmonic_avg_2D(p[LL_2D], p[LU_2D], p[UL_2D], p[UU_2D]);

  dTx[L_1D][T11] = calc_sym_grad_1D(dx, Tij[LL_2D][T11], Tij[UL_2D][T11]);
  dTx[L_1D][T12] = calc_sym_grad_1D(dx, Tij[LL_2D][T12], Tij[UL_2D][T12]);
  dTx[L_1D][T13] = calc_sym_grad_1D(dx, Tij[LL_2D][T13], Tij[UL_2D][T13]);
  dTx[L_1D][T22] = calc_sym_grad_1D(dx, Tij[LL_2D][T22], Tij[UL_2D][T22]);
  dTx[L_1D][T23] = calc_sym_grad_1D(dx, Tij[LL_2D][T23], Tij[UL_2D][T23]);
  dTx[L_1D][T33] = calc_sym_grad_1D(dx, Tij[LL_2D][T33], Tij[UL_2D][T33]);

  dTx[U_1D][T11] = calc_sym_grad_1D(dx, Tij[LU_2D][T11], Tij[UU_2D][T11]);
  dTx[U_1D][T12] = calc_sym_grad_1D(dx, Tij[LU_2D][T12], Tij[UU_2D][T12]);
  dTx[U_1D][T13] = calc_sym_grad_1D(dx, Tij[LU_2D][T13], Tij[UU_2D][T13]);
  dTx[U_1D][T22] = calc_sym_grad_1D(dx, Tij[LU_2D][T22], Tij[UU_2D][T22]);
  dTx[U_1D][T23] = calc_sym_grad_1D(dx, Tij[LU_2D][T23], Tij[UU_2D][T23]);
  dTx[U_1D][T33] = calc_sym_grad_1D(dx, Tij[LU_2D][T33], Tij[UU_2D][T33]);

  dTdx[L_1D][T11] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T11], dTx[U_1D][T11]);
  dTdx[L_1D][T12] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T12], dTx[U_1D][T12]);
  dTdx[L_1D][T13] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T13], dTx[U_1D][T13]);
  dTdx[L_1D][T22] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T22], dTx[U_1D][T22]);
  dTdx[L_1D][T23] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T23], dTx[U_1D][T23]);
  dTdx[L_1D][T33] = calc_sym_grad_limiter_2D(limit, dTx[L_1D][T33], dTx[U_1D][T33]);

  dTdx[U_1D][T11] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T11], dTx[L_1D][T11]);
  dTdx[U_1D][T12] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T12], dTx[L_1D][T12]);
  dTdx[U_1D][T13] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T13], dTx[L_1D][T13]);
  dTdx[U_1D][T22] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T22], dTx[L_1D][T22]);
  dTdx[U_1D][T23] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T23], dTx[L_1D][T23]);
  dTdx[U_1D][T33] = calc_sym_grad_limiter_2D(limit, dTx[U_1D][T33], dTx[L_1D][T33]);

  dTy[L_1D][T11] = calc_sym_grad_1D(dy, Tij[LL_2D][T11], Tij[LU_2D][T11]);
  dTy[L_1D][T12] = calc_sym_grad_1D(dy, Tij[LL_2D][T12], Tij[LU_2D][T12]);
  dTy[L_1D][T13] = calc_sym_grad_1D(dy, Tij[LL_2D][T13], Tij[LU_2D][T13]);
  dTy[L_1D][T22] = calc_sym_grad_1D(dy, Tij[LL_2D][T22], Tij[LU_2D][T22]);
  dTy[L_1D][T23] = calc_sym_grad_1D(dy, Tij[LL_2D][T23], Tij[LU_2D][T23]);
  dTy[L_1D][T33] = calc_sym_grad_1D(dy, Tij[LL_2D][T33], Tij[LU_2D][T33]);

  dTy[U_1D][T11] = calc_sym_grad_1D(dy, Tij[UL_2D][T11], Tij[UU_2D][T11]);
  dTy[U_1D][T12] = calc_sym_grad_1D(dy, Tij[UL_2D][T12], Tij[UU_2D][T12]);
  dTy[U_1D][T13] = calc_sym_grad_1D(dy, Tij[UL_2D][T13], Tij[UU_2D][T13]);
  dTy[U_1D][T22] = calc_sym_grad_1D(dy, Tij[UL_2D][T22], Tij[UU_2D][T22]);
  dTy[U_1D][T23] = calc_sym_grad_1D(dy, Tij[UL_2D][T23], Tij[UU_2D][T23]);
  dTy[U_1D][T33] = calc_sym_grad_1D(dy, Tij[UL_2D][T33], Tij[UU_2D][T33]);

  dTdy[L_1D][T11] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T11], dTy[U_1D][T11]);
  dTdy[L_1D][T12] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T12], dTy[U_1D][T12]);
  dTdy[L_1D][T13] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T13], dTy[U_1D][T13]);
  dTdy[L_1D][T22] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T22], dTy[U_1D][T22]);
  dTdy[L_1D][T23] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T23], dTy[U_1D][T23]);
  dTdy[L_1D][T33] = calc_sym_grad_limiter_2D(limit, dTy[L_1D][T33], dTy[U_1D][T33]);

  dTdy[U_1D][T11] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T11], dTy[L_1D][T11]);
  dTdy[U_1D][T12] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T12], dTy[L_1D][T12]);
  dTdy[U_1D][T13] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T13], dTy[L_1D][T13]);
  dTdy[U_1D][T22] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T22], dTy[L_1D][T22]);
  dTdy[U_1D][T23] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T23], dTy[L_1D][T23]);
  dTdy[U_1D][T33] = calc_sym_grad_limiter_2D(limit, dTy[U_1D][T33], dTy[L_1D][T33]);

  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  int compx[4] = { L_1D, U_1D, L_1D, U_1D };
  int compy[4] = { L_1D, L_1D, U_1D, U_1D };

  q[LL_2D*10 + Q111] = chi*dTdx[compx[LL_2D]][T11];
  q[LL_2D*10 + Q112] = chi*(2.0*dTdx[compx[LL_2D]][T12]
    + dTdy[compy[LL_2D]][T11])/3.0;
  q[LL_2D*10 + Q113] = chi*2.0*dTdx[compx[LL_2D]][T13]/3.0;
  q[LL_2D*10 + Q122] = chi*(dTdx[compx[LL_2D]][T22]
    + 2.0*dTdy[compy[LL_2D]][T12])/3.0;
  q[LL_2D*10 + Q123] = chi*(dTdx[compx[LL_2D]][T23]
    + dTdy[compy[LL_2D]][T13])/3.0;
  q[LL_2D*10 + Q133] = chi*dTdx[compx[LL_2D]][T33]/3.0;
  q[LL_2D*10 + Q222] = chi*dTdy[compy[LL_2D]][T22];
  q[LL_2D*10 + Q223] = chi*2.0*dTdy[compy[LL_2D]][T23]/3.0;
  q[LL_2D*10 + Q233] = chi*dTdy[compy[LL_2D]][T33]/3.0;

  q[LU_2D*10 + Q111] = chi*dTdx[compx[LU_2D]][T11];
  q[LU_2D*10 + Q112] = chi*(2.0*dTdx[compx[LU_2D]][T12]
    + dTdy[compy[LU_2D]][T11])/3.0;
  q[LU_2D*10 + Q113] = chi*2.0*dTdx[compx[LU_2D]][T13]/3.0;
  q[LU_2D*10 + Q122] = chi*(dTdx[compx[LU_2D]][T22]
    + 2.0*dTdy[compy[LU_2D]][T12])/3.0;
  q[LU_2D*10 + Q123] = chi*(dTdx[compx[LU_2D]][T23]
    + dTdy[compy[LU_2D]][T13])/3.0;
  q[LU_2D*10 + Q133] = chi*dTdx[compx[LU_2D]][T33]/3.0;
  q[LU_2D*10 + Q222] = chi*dTdy[compy[LU_2D]][T22];
  q[LU_2D*10 + Q223] = chi*2.0*dTdy[compy[LU_2D]][T23]/3.0;
  q[LU_2D*10 + Q233] = chi*dTdy[compy[LU_2D]][T33]/3.0;

  q[UL_2D*10 + Q111] = chi*dTdx[compx[UL_2D]][T11];
  q[UL_2D*10 + Q112] = chi*(2.0*dTdx[compx[UL_2D]][T12]
    + dTdy[compy[UL_2D]][T11])/3.0;
  q[UL_2D*10 + Q113] = chi*2.0*dTdx[compx[UL_2D]][T13]/3.0;
  q[UL_2D*10 + Q122] = chi*(dTdx[compx[UL_2D]][T22]
    + 2.0*dTdy[compy[UL_2D]][T12])/3.0;
  q[UL_2D*10 + Q123] = chi*(dTdx[compx[UL_2D]][T23]
    + dTdy[compy[UL_2D]][T13])/3.0;
  q[UL_2D*10 + Q133] = chi*dTdx[compx[UL_2D]][T33]/3.0;
  q[UL_2D*10 + Q222] = chi*dTdy[compy[UL_2D]][T22];
  q[UL_2D*10 + Q223] = chi*2.0*dTdy[compy[UL_2D]][T23]/3.0;
  q[UL_2D*10 + Q233] = chi*dTdy[compy[UL_2D]][T33]/3.0;

  q[UU_2D*10 + Q111] = chi*dTdx[compx[UU_2D]][T11];
  q[UU_2D*10 + Q112] = chi*(2.0*dTdx[compx[UU_2D]][T12]
    + dTdy[compy[UU_2D]][T11])/3.0;
  q[UU_2D*10 + Q113] = chi*2.0*dTdx[compx[UU_2D]][T13]/3.0;
  q[UU_2D*10 + Q122] = chi*(dTdx[compx[UU_2D]][T22]
    + 2.0*dTdy[compy[UU_2D]][T12])/3.0;
  q[UU_2D*10 + Q123] = chi*(dTdx[compx[UU_2D]][T23]
    + dTdy[compy[UU_2D]][T13])/3.0;
  q[UU_2D*10 + Q133] = chi*dTdx[compx[UU_2D]][T33]/3.0;
  q[UU_2D*10 + Q222] = chi*dTdy[compy[UU_2D]][T22];
  q[UU_2D*10 + Q223] = chi*2.0*dTdy[compy[UU_2D]][T23]/3.0;
  q[UU_2D*10 + Q233] = chi*dTdy[compy[UU_2D]][T33]/3.0;

  double da = fmin(dx, dy);
  double cfla = dt/(da*da);
  cflrate[0] = alpha*vth_avg*cfla;
}

GKYL_CU_D
static void
grad_closure_update_2d(const gkyl_ten_moment_grad_closure *gces,
  const double *q[], double *rhs)
{
  double div_qx[6] = {0.0};
  double div_qy[6] = {0.0};

  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];

  div_qx[0] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q111],
    q[LU_2D][UL_2D*10 + Q111], q[UL_2D][LU_2D*10 + Q111], q[UU_2D][LL_2D*10 + Q111]);
  div_qx[1] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q112],
    q[LU_2D][UL_2D*10 + Q112], q[UL_2D][LU_2D*10 + Q112], q[UU_2D][LL_2D*10 + Q112]);
  div_qx[2] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q113],
    q[LU_2D][UL_2D*10 + Q113], q[UL_2D][LU_2D*10 + Q113], q[UU_2D][LL_2D*10 + Q113]);
  div_qx[3] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q122],
    q[LU_2D][UL_2D*10 + Q122], q[UL_2D][LU_2D*10 + Q122], q[UU_2D][LL_2D*10 + Q122]);
  div_qx[4] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q123],
    q[LU_2D][UL_2D*10 + Q123], q[UL_2D][LU_2D*10 + Q123], q[UU_2D][LL_2D*10 + Q123]);
  div_qx[5] = calc_sym_gradx_2D(dx, q[LL_2D][UU_2D*10 + Q133],
    q[LU_2D][UL_2D*10 + Q133], q[UL_2D][LU_2D*10 + Q133], q[UU_2D][LL_2D*10 + Q133]);

  div_qy[0] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q112],
    q[LU_2D][UL_2D*10 + Q112], q[UL_2D][LU_2D*10 + Q112], q[UU_2D][LL_2D*10 + Q112]);
  div_qy[1] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q122],
    q[LU_2D][UL_2D*10 + Q122], q[UL_2D][LU_2D*10 + Q122], q[UU_2D][LL_2D*10 + Q122]);
  div_qy[2] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q123],
    q[LU_2D][UL_2D*10 + Q123], q[UL_2D][LU_2D*10 + Q123], q[UU_2D][LL_2D*10 + Q123]);
  div_qy[3] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q222],
    q[LU_2D][UL_2D*10 + Q222], q[UL_2D][LU_2D*10 + Q222], q[UU_2D][LL_2D*10 + Q222]);
  div_qy[4] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q223],
    q[LU_2D][UL_2D*10 + Q223], q[UL_2D][LU_2D*10 + Q223], q[UU_2D][LL_2D*10 + Q223]);
  div_qy[5] = calc_sym_grady_2D(dy, q[LL_2D][UU_2D*10 + Q233],
    q[LU_2D][UL_2D*10 + Q233], q[UL_2D][LU_2D*10 + Q233], q[UU_2D][LL_2D*10 + Q233]);

  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = div_qx[0] + div_qy[0];
  rhs[P12] = div_qx[1] + div_qy[1];
  rhs[P13] = div_qx[2] + div_qy[2];
  rhs[P22] = div_qx[3] + div_qy[3];
  rhs[P23] = div_qx[4] + div_qy[4];
  rhs[P33] = div_qx[5] + div_qy[5];
}

GKYL_CU_D
static void
calc_unmag_heat_flux_3d(const gkyl_ten_moment_grad_closure *gces,
  const double *fluid_d[], double *cflrate, double dt, double *q)
{
  const int ndim = gces->ndim;
  double rho_avg = 0.0;
  double p_avg = 0.0;
  double limit = 0.75;

  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  const double dz = gces->grid.dx[2];

  double dTx[4][6] = {0.0};
  double dTy[4][6] = {0.0};
  double dTz[4][6] = {0.0};
  double dTdx[4][6] = {0.0};
  double dTdy[4][6] = {0.0};
  double dTdz[4][6] = {0.0};
  double Tij[8][6] = {0.0};
  double rho[8] = {0.0};
  double p[8] = {0.0};
  var_setup(gces, LLL_3D, UUU_3D, fluid_d, rho, p, Tij);

  rho_avg = calc_harmonic_avg_3D(rho[LLL_3D], rho[LLU_3D], rho[LUL_3D], rho[LUU_3D],
                                 rho[ULL_3D], rho[ULU_3D], rho[UUL_3D], rho[UUU_3D]);
  p_avg = calc_harmonic_avg_3D(p[LLL_3D], p[LLU_3D], p[LUL_3D], p[LUU_3D],
                               p[ULL_3D], p[ULU_3D], p[UUL_3D], p[UUU_3D]);

  dTx[LL_2D][T11] = calc_sym_grad_1D(dx, Tij[LLL_3D][T11], Tij[ULL_3D][T11]);
  dTx[LL_2D][T12] = calc_sym_grad_1D(dx, Tij[LLL_3D][T12], Tij[ULL_3D][T12]);
  dTx[LL_2D][T13] = calc_sym_grad_1D(dx, Tij[LLL_3D][T13], Tij[ULL_3D][T13]);
  dTx[LL_2D][T22] = calc_sym_grad_1D(dx, Tij[LLL_3D][T22], Tij[ULL_3D][T22]);
  dTx[LL_2D][T23] = calc_sym_grad_1D(dx, Tij[LLL_3D][T23], Tij[ULL_3D][T23]);
  dTx[LL_2D][T33] = calc_sym_grad_1D(dx, Tij[LLL_3D][T33], Tij[ULL_3D][T33]);

  dTx[LU_2D][T11] = calc_sym_grad_1D(dx, Tij[LUL_3D][T11], Tij[UUL_3D][T11]);
  dTx[LU_2D][T12] = calc_sym_grad_1D(dx, Tij[LUL_3D][T12], Tij[UUL_3D][T12]);
  dTx[LU_2D][T13] = calc_sym_grad_1D(dx, Tij[LUL_3D][T13], Tij[UUL_3D][T13]);
  dTx[LU_2D][T22] = calc_sym_grad_1D(dx, Tij[LUL_3D][T22], Tij[UUL_3D][T22]);
  dTx[LU_2D][T23] = calc_sym_grad_1D(dx, Tij[LUL_3D][T23], Tij[UUL_3D][T23]);
  dTx[LU_2D][T33] = calc_sym_grad_1D(dx, Tij[LUL_3D][T33], Tij[UUL_3D][T33]);
  
  dTx[UL_2D][T11] = calc_sym_grad_1D(dx, Tij[LLU_3D][T11], Tij[ULU_3D][T11]);
  dTx[UL_2D][T12] = calc_sym_grad_1D(dx, Tij[LLU_3D][T12], Tij[ULU_3D][T12]);
  dTx[UL_2D][T13] = calc_sym_grad_1D(dx, Tij[LLU_3D][T13], Tij[ULU_3D][T13]);
  dTx[UL_2D][T22] = calc_sym_grad_1D(dx, Tij[LLU_3D][T22], Tij[ULU_3D][T22]);
  dTx[UL_2D][T23] = calc_sym_grad_1D(dx, Tij[LLU_3D][T23], Tij[ULU_3D][T23]);
  dTx[UL_2D][T33] = calc_sym_grad_1D(dx, Tij[LLU_3D][T33], Tij[ULU_3D][T33]);
  
  dTx[UU_2D][T11] = calc_sym_grad_1D(dx, Tij[LUU_3D][T11], Tij[UUU_3D][T11]);
  dTx[UU_2D][T12] = calc_sym_grad_1D(dx, Tij[LUU_3D][T12], Tij[UUU_3D][T12]);
  dTx[UU_2D][T13] = calc_sym_grad_1D(dx, Tij[LUU_3D][T13], Tij[UUU_3D][T13]);
  dTx[UU_2D][T22] = calc_sym_grad_1D(dx, Tij[LUU_3D][T22], Tij[UUU_3D][T22]);
  dTx[UU_2D][T23] = calc_sym_grad_1D(dx, Tij[LUU_3D][T23], Tij[UUU_3D][T23]);  
  dTx[UU_2D][T33] = calc_sym_grad_1D(dx, Tij[LUU_3D][T33], Tij[UUU_3D][T33]);

  dTdx[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T11], dTx[LU_2D][T11],
    dTx[UL_2D][T11], dTx[UU_2D][T11]);
  dTdx[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T12], dTx[LU_2D][T12],
    dTx[UL_2D][T12], dTx[UU_2D][T12]);
  dTdx[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T13], dTx[LU_2D][T13],
    dTx[UL_2D][T13], dTx[UU_2D][T13]);
  dTdx[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T22], dTx[LU_2D][T22],
    dTx[UL_2D][T22], dTx[UU_2D][T22]);
  dTdx[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T23], dTx[LU_2D][T23],
    dTx[UL_2D][T23], dTx[UU_2D][T23]);
  dTdx[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[LL_2D][T33], dTx[LU_2D][T33],
    dTx[UL_2D][T33], dTx[UU_2D][T33]);

  dTdx[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T11], dTx[LL_2D][T11],
    dTx[UL_2D][T11], dTx[UU_2D][T11]);
  dTdx[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T12], dTx[LL_2D][T12],
    dTx[UL_2D][T12], dTx[UU_2D][T12]);
  dTdx[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T13], dTx[LL_2D][T13],
    dTx[UL_2D][T13], dTx[UU_2D][T13]);
  dTdx[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T22], dTx[LL_2D][T22],
    dTx[UL_2D][T22], dTx[UU_2D][T22]);
  dTdx[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T23], dTx[LL_2D][T23],
    dTx[UL_2D][T23], dTx[UU_2D][T23]);
  dTdx[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[LU_2D][T33], dTx[LL_2D][T33],
    dTx[UL_2D][T33], dTx[UU_2D][T33]);

  dTdx[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T11], dTx[LL_2D][T11],
    dTx[LU_2D][T11], dTx[UU_2D][T11]);
  dTdx[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T12], dTx[LL_2D][T12],
    dTx[LU_2D][T12], dTx[UU_2D][T12]);
  dTdx[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T13], dTx[LL_2D][T13],
    dTx[LU_2D][T13], dTx[UU_2D][T13]);
  dTdx[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T22], dTx[LL_2D][T22],
    dTx[LU_2D][T22], dTx[UU_2D][T22]);
  dTdx[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T23], dTx[LL_2D][T23],
    dTx[LU_2D][T23], dTx[UU_2D][T23]);
  dTdx[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T33], dTx[LL_2D][T33],
    dTx[LU_2D][T33], dTx[UU_2D][T33]);

  dTdx[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T11], dTx[LL_2D][T11],
    dTx[LU_2D][T11], dTx[UU_2D][T11]);
  dTdx[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T12], dTx[LL_2D][T12],
    dTx[LU_2D][T12], dTx[UU_2D][T12]);
  dTdx[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T13], dTx[LL_2D][T13],
    dTx[LU_2D][T13], dTx[UU_2D][T13]);
  dTdx[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T22], dTx[LL_2D][T22],
    dTx[LU_2D][T22], dTx[UU_2D][T22]);
  dTdx[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T23], dTx[LL_2D][T23],
    dTx[LU_2D][T23], dTx[UU_2D][T23]);
  dTdx[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTx[UL_2D][T33], dTx[LL_2D][T33],
    dTx[LU_2D][T33], dTx[UU_2D][T33]);

  dTy[LL_2D][T11] = calc_sym_grad_1D(dy, Tij[LLL_3D][T11], Tij[LUL_3D][T11]);
  dTy[LL_2D][T12] = calc_sym_grad_1D(dy, Tij[LLL_3D][T12], Tij[LUL_3D][T12]);
  dTy[LL_2D][T13] = calc_sym_grad_1D(dy, Tij[LLL_3D][T13], Tij[LUL_3D][T13]);
  dTy[LL_2D][T22] = calc_sym_grad_1D(dy, Tij[LLL_3D][T22], Tij[LUL_3D][T22]);
  dTy[LL_2D][T23] = calc_sym_grad_1D(dy, Tij[LLL_3D][T23], Tij[LUL_3D][T23]);
  dTy[LL_2D][T33] = calc_sym_grad_1D(dy, Tij[LLL_3D][T33], Tij[LUL_3D][T33]);

  dTy[LU_2D][T11] = calc_sym_grad_1D(dy, Tij[ULL_3D][T11], Tij[UUL_3D][T11]);
  dTy[LU_2D][T12] = calc_sym_grad_1D(dy, Tij[ULL_3D][T12], Tij[UUL_3D][T12]);
  dTy[LU_2D][T13] = calc_sym_grad_1D(dy, Tij[ULL_3D][T13], Tij[UUL_3D][T13]);
  dTy[LU_2D][T22] = calc_sym_grad_1D(dy, Tij[ULL_3D][T22], Tij[UUL_3D][T22]);
  dTy[LU_2D][T23] = calc_sym_grad_1D(dy, Tij[ULL_3D][T23], Tij[UUL_3D][T23]);
  dTy[LU_2D][T33] = calc_sym_grad_1D(dy, Tij[ULL_3D][T33], Tij[UUL_3D][T33]);

  dTy[UL_2D][T11] = calc_sym_grad_1D(dy, Tij[LLU_3D][T11], Tij[LUU_3D][T11]);
  dTy[UL_2D][T12] = calc_sym_grad_1D(dy, Tij[LLU_3D][T12], Tij[LUU_3D][T12]);
  dTy[UL_2D][T13] = calc_sym_grad_1D(dy, Tij[LLU_3D][T13], Tij[LUU_3D][T13]);
  dTy[UL_2D][T22] = calc_sym_grad_1D(dy, Tij[LLU_3D][T22], Tij[LUU_3D][T22]);
  dTy[UL_2D][T23] = calc_sym_grad_1D(dy, Tij[LLU_3D][T23], Tij[LUU_3D][T23]);
  dTy[UL_2D][T33] = calc_sym_grad_1D(dy, Tij[LLU_3D][T33], Tij[LUU_3D][T33]);

  dTy[UU_2D][T11] = calc_sym_grad_1D(dy, Tij[ULU_3D][T11], Tij[UUU_3D][T11]);
  dTy[UU_2D][T12] = calc_sym_grad_1D(dy, Tij[ULU_3D][T12], Tij[UUU_3D][T12]);
  dTy[UU_2D][T13] = calc_sym_grad_1D(dy, Tij[ULU_3D][T13], Tij[UUU_3D][T13]);
  dTy[UU_2D][T22] = calc_sym_grad_1D(dy, Tij[ULU_3D][T22], Tij[UUU_3D][T22]);
  dTy[UU_2D][T23] = calc_sym_grad_1D(dy, Tij[ULU_3D][T23], Tij[UUU_3D][T23]);
  dTy[UU_2D][T33] = calc_sym_grad_1D(dy, Tij[ULU_3D][T33], Tij[UUU_3D][T33]);

  dTdy[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T11], dTy[LU_2D][T11],
    dTy[UL_2D][T11], dTy[UU_2D][T11]);
  dTdy[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T12], dTy[LU_2D][T12],
    dTy[UL_2D][T12], dTy[UU_2D][T12]);
  dTdy[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T13], dTy[LU_2D][T13],
    dTy[UL_2D][T13], dTy[UU_2D][T13]);
  dTdy[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T22], dTy[LU_2D][T22],
    dTy[UL_2D][T22], dTy[UU_2D][T22]);
  dTdy[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T23], dTy[LU_2D][T23],
    dTy[UL_2D][T23], dTy[UU_2D][T23]);
  dTdy[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[LL_2D][T33], dTy[LU_2D][T33],
    dTy[UL_2D][T33], dTy[UU_2D][T33]);

  dTdy[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T11], dTy[LL_2D][T11],
    dTy[UL_2D][T11], dTy[UU_2D][T11]);
  dTdy[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T12], dTy[LL_2D][T12],
    dTy[UL_2D][T12], dTy[UU_2D][T12]);
  dTdy[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T13], dTy[LL_2D][T13],
    dTy[UL_2D][T13], dTy[UU_2D][T13]);
  dTdy[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T22], dTy[LL_2D][T22],
    dTy[UL_2D][T22], dTy[UU_2D][T22]);
  dTdy[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T23], dTy[LL_2D][T23],
    dTy[UL_2D][T23], dTy[UU_2D][T23]);
  dTdy[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[LU_2D][T33], dTy[LL_2D][T33],
    dTy[UL_2D][T33], dTy[UU_2D][T33]);

  dTdy[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T11], dTy[LL_2D][T11],
    dTy[LU_2D][T11], dTy[UU_2D][T11]);
  dTdy[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T12], dTy[LL_2D][T12],
    dTy[LU_2D][T12], dTy[UU_2D][T12]);
  dTdy[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T13], dTy[LL_2D][T13],
    dTy[LU_2D][T13], dTy[UU_2D][T13]);
  dTdy[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T22], dTy[LL_2D][T22],
    dTy[LU_2D][T22], dTy[UU_2D][T22]);
  dTdy[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T23], dTy[LL_2D][T23],
    dTy[LU_2D][T23], dTy[UU_2D][T23]);
  dTdy[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T33], dTy[LL_2D][T33],
    dTy[LU_2D][T33], dTy[UU_2D][T33]);

  dTdy[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T11], dTy[LL_2D][T11],
    dTy[LU_2D][T11], dTy[UU_2D][T11]);
  dTdy[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T12], dTy[LL_2D][T12],
    dTy[LU_2D][T12], dTy[UU_2D][T12]);
  dTdy[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T13], dTy[LL_2D][T13],
    dTy[LU_2D][T13], dTy[UU_2D][T13]);
  dTdy[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T22], dTy[LL_2D][T22],
    dTy[LU_2D][T22], dTy[UU_2D][T22]);
  dTdy[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T23], dTy[LL_2D][T23],
    dTy[LU_2D][T23], dTy[UU_2D][T23]);
  dTdy[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTy[UL_2D][T33], dTy[LL_2D][T33],
    dTy[LU_2D][T33], dTy[UU_2D][T33]);

  dTz[LL_2D][T11] = calc_sym_grad_1D(dz, Tij[LLL_3D][T11], Tij[LLU_3D][T11]);
  dTz[LL_2D][T12] = calc_sym_grad_1D(dz, Tij[LLL_3D][T12], Tij[LLU_3D][T12]);
  dTz[LL_2D][T13] = calc_sym_grad_1D(dz, Tij[LLL_3D][T13], Tij[LLU_3D][T13]);
  dTz[LL_2D][T22] = calc_sym_grad_1D(dz, Tij[LLL_3D][T22], Tij[LLU_3D][T22]);
  dTz[LL_2D][T23] = calc_sym_grad_1D(dz, Tij[LLL_3D][T23], Tij[LLU_3D][T23]);
  dTz[LL_2D][T33] = calc_sym_grad_1D(dz, Tij[LLL_3D][T33], Tij[LLU_3D][T33]);

  dTz[LU_2D][T11] = calc_sym_grad_1D(dz, Tij[ULL_3D][T11], Tij[ULU_3D][T11]);
  dTz[LU_2D][T12] = calc_sym_grad_1D(dz, Tij[ULL_3D][T12], Tij[ULU_3D][T12]);
  dTz[LU_2D][T13] = calc_sym_grad_1D(dz, Tij[ULL_3D][T13], Tij[ULU_3D][T13]);
  dTz[LU_2D][T22] = calc_sym_grad_1D(dz, Tij[ULL_3D][T22], Tij[ULU_3D][T22]);
  dTz[LU_2D][T23] = calc_sym_grad_1D(dz, Tij[ULL_3D][T23], Tij[ULU_3D][T23]);
  dTz[LU_2D][T33] = calc_sym_grad_1D(dz, Tij[ULL_3D][T33], Tij[ULU_3D][T33]);

  dTz[UL_2D][T11] = calc_sym_grad_1D(dz, Tij[LUL_3D][T11], Tij[LUU_3D][T11]);
  dTz[UL_2D][T12] = calc_sym_grad_1D(dz, Tij[LUL_3D][T12], Tij[LUU_3D][T12]);
  dTz[UL_2D][T13] = calc_sym_grad_1D(dz, Tij[LUL_3D][T13], Tij[LUU_3D][T13]);
  dTz[UL_2D][T22] = calc_sym_grad_1D(dz, Tij[LUL_3D][T22], Tij[LUU_3D][T22]);
  dTz[UL_2D][T23] = calc_sym_grad_1D(dz, Tij[LUL_3D][T23], Tij[LUU_3D][T23]);
  dTz[UL_2D][T33] = calc_sym_grad_1D(dz, Tij[LUL_3D][T33], Tij[LUU_3D][T33]);

  dTz[UU_2D][T11] = calc_sym_grad_1D(dz, Tij[UUL_3D][T11], Tij[UUU_3D][T11]);
  dTz[UU_2D][T12] = calc_sym_grad_1D(dz, Tij[UUL_3D][T12], Tij[UUU_3D][T12]);
  dTz[UU_2D][T13] = calc_sym_grad_1D(dz, Tij[UUL_3D][T13], Tij[UUU_3D][T13]);
  dTz[UU_2D][T22] = calc_sym_grad_1D(dz, Tij[UUL_3D][T22], Tij[UUU_3D][T22]);
  dTz[UU_2D][T23] = calc_sym_grad_1D(dz, Tij[UUL_3D][T23], Tij[UUU_3D][T23]);
  dTz[UU_2D][T33] = calc_sym_grad_1D(dz, Tij[UUL_3D][T33], Tij[UUU_3D][T33]);

  dTdz[LL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T11], dTz[LU_2D][T11],
    dTz[UL_2D][T11], dTz[UU_2D][T11]);
  dTdz[LL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T12], dTz[LU_2D][T12],
    dTz[UL_2D][T12], dTz[UU_2D][T12]);
  dTdz[LL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T13], dTz[LU_2D][T13],
    dTz[UL_2D][T13], dTz[UU_2D][T13]);
  dTdz[LL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T22], dTz[LU_2D][T22],
    dTz[UL_2D][T22], dTz[UU_2D][T22]);
  dTdz[LL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T23], dTz[LU_2D][T23],
    dTz[UL_2D][T23], dTz[UU_2D][T23]);
  dTdz[LL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[LL_2D][T33], dTz[LU_2D][T33],
    dTz[UL_2D][T33], dTz[UU_2D][T33]);

  dTdz[LU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T11], dTz[LL_2D][T11],
    dTz[UL_2D][T11], dTz[UU_2D][T11]);
  dTdz[LU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T12], dTz[LL_2D][T12],
    dTz[UL_2D][T12], dTz[UU_2D][T12]);
  dTdz[LU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T13], dTz[LL_2D][T13],
    dTz[UL_2D][T13], dTz[UU_2D][T13]);
  dTdz[LU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T22], dTz[LL_2D][T22],
    dTz[UL_2D][T22], dTz[UU_2D][T22]);
  dTdz[LU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T23], dTz[LL_2D][T23],
    dTz[UL_2D][T23], dTz[UU_2D][T23]);
  dTdz[LU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[LU_2D][T33], dTz[LL_2D][T33],
    dTz[UL_2D][T33], dTz[UU_2D][T33]);

  dTdz[UL_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T11], dTz[LL_2D][T11],
    dTz[LU_2D][T11], dTz[UU_2D][T11]);
  dTdz[UL_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T12], dTz[LL_2D][T12],
    dTz[LU_2D][T12], dTz[UU_2D][T12]);
  dTdz[UL_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T13], dTz[LL_2D][T13],
    dTz[LU_2D][T13], dTz[UU_2D][T13]);
  dTdz[UL_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T22], dTz[LL_2D][T22],
    dTz[LU_2D][T22], dTz[UU_2D][T22]);
  dTdz[UL_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T23], dTz[LL_2D][T23],
    dTz[LU_2D][T23], dTz[UU_2D][T23]);
  dTdz[UL_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T33], dTz[LL_2D][T33],
    dTz[LU_2D][T33], dTz[UU_2D][T33]);

  dTdz[UU_2D][T11] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T11], dTz[LL_2D][T11],
    dTz[LU_2D][T11], dTz[UU_2D][T11]);
  dTdz[UU_2D][T12] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T12], dTz[LL_2D][T12],
    dTz[LU_2D][T12], dTz[UU_2D][T12]);
  dTdz[UU_2D][T13] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T13], dTz[LL_2D][T13],
    dTz[LU_2D][T13], dTz[UU_2D][T13]);
  dTdz[UU_2D][T22] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T22], dTz[LL_2D][T22],
    dTz[LU_2D][T22], dTz[UU_2D][T22]);
  dTdz[UU_2D][T23] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T23], dTz[LL_2D][T23],
    dTz[LU_2D][T23], dTz[UU_2D][T23]);
  dTdz[UU_2D][T33] = calc_sym_grad_limiter_3D(limit, dTz[UL_2D][T33], dTz[LL_2D][T33],
    dTz[LU_2D][T33], dTz[UU_2D][T33]);
  
  double alpha = 1.0/gces->k0;
  double vth_avg = sqrt(p_avg/rho_avg);

  // Temperature is actually T/m due to the formulation of var_setup.
  // Thus, the mass density rho is used instead of number density n.
  double chi = alpha*vth_avg*rho_avg;

  int compx[8] = { LL_2D, UL_2D, LU_2D, UU_2D, LL_2D, UL_2D, LU_2D, UU_2D };
  int compy[8] = { LL_2D, UL_2D, LL_2D, UL_2D, LU_2D, UU_2D, LU_2D, UU_2D };
  int compz[8] = { LL_2D, LL_2D, UL_2D, UL_2D, LU_2D, LU_2D, UU_2D, UU_2D };

  q[LLL_3D*10 + Q111] = chi*dTdx[compx[LLL_3D]][T11];
  q[LLL_3D*10 + Q112] = chi*(2.0*dTdx[compx[LLL_3D]][T12]
    + dTdy[compy[LLL_3D]][T11])/3.0;
  q[LLL_3D*10 + Q113] = chi*(2.0*dTdx[compx[LLL_3D]][T13]
    + dTdz[compz[LLL_3D]][T11])/3.0;
  q[LLL_3D*10 + Q122] = chi*(dTdx[compx[LLL_3D]][T22]
    + 2.0*dTdy[compy[LLL_3D]][T12])/3.0;
  q[LLL_3D*10 + Q123] = chi*(dTdx[compx[LLL_3D]][T23]
    + dTdy[compy[LLL_3D]][T13] + dTdz[compz[LLL_3D]][T12])/3.0;
  q[LLL_3D*10 + Q133] = chi*(dTdx[compx[LLL_3D]][T33]
    + 2.0*dTdz[compz[LLL_3D]][T13])/3.0;
  q[LLL_3D*10 + Q222] = chi*dTdy[compy[LLL_3D]][T22];
  q[LLL_3D*10 + Q223] = chi*(2.0*dTdy[compy[LLL_3D]][T23]
    + dTdz[compz[LLL_3D]][T22])/3.0;
  q[LLL_3D*10 + Q233] = chi*(dTdy[compy[LLL_3D]][T33]
    + 2.0*dTdz[compz[LLL_3D]][T23])/3.0;
  q[LLL_3D*10 + Q333] = chi*dTdz[compz[LLL_3D]][T33];

  q[LLU_3D*10 + Q111] = chi*dTdx[compx[LLU_3D]][T11];
  q[LLU_3D*10 + Q112] = chi*(2.0*dTdx[compx[LLU_3D]][T12]
    + dTdy[compy[LLU_3D]][T11])/3.0;
  q[LLU_3D*10 + Q113] = chi*(2.0*dTdx[compx[LLU_3D]][T13]
    + dTdz[compz[LLU_3D]][T11])/3.0;
  q[LLU_3D*10 + Q122] = chi*(dTdx[compx[LLU_3D]][T22]
    + 2.0*dTdy[compy[LLU_3D]][T12])/3.0;
  q[LLU_3D*10 + Q123] = chi*(dTdx[compx[LLU_3D]][T23]
    + dTdy[compy[LLU_3D]][T13] + dTdz[compz[LLU_3D]][T12])/3.0;
  q[LLU_3D*10 + Q133] = chi*(dTdx[compx[LLU_3D]][T33]
    + 2.0*dTdz[compz[LLU_3D]][T13])/3.0;
  q[LLU_3D*10 + Q222] = chi*dTdy[compy[LLU_3D]][T22];
  q[LLU_3D*10 + Q223] = chi*(2.0*dTdy[compy[LLU_3D]][T23]
    + dTdz[compz[LLU_3D]][T22])/3.0;
  q[LLU_3D*10 + Q233] = chi*(dTdy[compy[LLU_3D]][T33]
    + 2.0*dTdz[compz[LLU_3D]][T23])/3.0;
  q[LLU_3D*10 + Q333] = chi*dTdz[compz[LLU_3D]][T33];

  q[LUL_3D*10 + Q111] = chi*dTdx[compx[LUL_3D]][T11];
  q[LUL_3D*10 + Q112] = chi*(2.0*dTdx[compx[LUL_3D]][T12]
    + dTdy[compy[LUL_3D]][T11])/3.0;
  q[LUL_3D*10 + Q113] = chi*(2.0*dTdx[compx[LUL_3D]][T13]
    + dTdz[compz[LUL_3D]][T11])/3.0;
  q[LUL_3D*10 + Q122] = chi*(dTdx[compx[LUL_3D]][T22]
    + 2.0*dTdy[compy[LUL_3D]][T12])/3.0;
  q[LUL_3D*10 + Q123] = chi*(dTdx[compx[LUL_3D]][T23]
    + dTdy[compy[LUL_3D]][T13] + dTdz[compz[LUL_3D]][T12])/3.0;
  q[LUL_3D*10 + Q133] = chi*(dTdx[compx[LUL_3D]][T33]
    + 2.0*dTdz[compz[LUL_3D]][T13])/3.0;
  q[LUL_3D*10 + Q222] = chi*dTdy[compy[LUL_3D]][T22];
  q[LUL_3D*10 + Q223] = chi*(2.0*dTdy[compy[LUL_3D]][T23]
    + dTdz[compz[LUL_3D]][T22])/3.0;
  q[LUL_3D*10 + Q233] = chi*(dTdy[compy[LUL_3D]][T33]
    + 2.0*dTdz[compz[LUL_3D]][T23])/3.0;
  q[LUL_3D*10 + Q333] = chi*dTdz[compz[LUL_3D]][T33];

  q[LUU_3D*10 + Q111] = chi*dTdx[compx[LUU_3D]][T11];
  q[LUU_3D*10 + Q112] = chi*(2.0*dTdx[compx[LUU_3D]][T12]
    + dTdy[compy[LUU_3D]][T11])/3.0;
  q[LUU_3D*10 + Q113] = chi*(2.0*dTdx[compx[LUU_3D]][T13]
    + dTdz[compz[LUU_3D]][T11])/3.0;
  q[LUU_3D*10 + Q122] = chi*(dTdx[compx[LUU_3D]][T22]
    + 2.0*dTdy[compy[LUU_3D]][T12])/3.0;
  q[LUU_3D*10 + Q123] = chi*(dTdx[compx[LUU_3D]][T23]
    + dTdy[compy[LUU_3D]][T13] + dTdz[compz[LUU_3D]][T12])/3.0;
  q[LUU_3D*10 + Q133] = chi*(dTdx[compx[LUU_3D]][T33]
    + 2.0*dTdz[compz[LUU_3D]][T13])/3.0;
  q[LUU_3D*10 + Q222] = chi*dTdy[compy[LUU_3D]][T22];
  q[LUU_3D*10 + Q223] = chi*(2.0*dTdy[compy[LUU_3D]][T23]
    + dTdz[compz[LUU_3D]][T22])/3.0;
  q[LUU_3D*10 + Q233] = chi*(dTdy[compy[LUU_3D]][T33]
    + 2.0*dTdz[compz[LUU_3D]][T23])/3.0;
  q[LUU_3D*10 + Q333] = chi*dTdz[compz[LUU_3D]][T33];

  q[ULL_3D*10 + Q111] = chi*dTdx[compx[ULL_3D]][T11];
  q[ULL_3D*10 + Q112] = chi*(2.0*dTdx[compx[ULL_3D]][T12]
    + dTdy[compy[ULL_3D]][T11])/3.0;
  q[ULL_3D*10 + Q113] = chi*(2.0*dTdx[compx[ULL_3D]][T13]
    + dTdz[compz[ULL_3D]][T11])/3.0;
  q[ULL_3D*10 + Q122] = chi*(dTdx[compx[ULL_3D]][T22]
    + 2.0*dTdy[compy[ULL_3D]][T12])/3.0;
  q[ULL_3D*10 + Q123] = chi*(dTdx[compx[ULL_3D]][T23]
    + dTdy[compy[ULL_3D]][T13] + dTdz[compz[ULL_3D]][T12])/3.0;
  q[ULL_3D*10 + Q133] = chi*(dTdx[compx[ULL_3D]][T33]
    + 2.0*dTdz[compz[ULL_3D]][T13])/3.0;
  q[ULL_3D*10 + Q222] = chi*dTdy[compy[ULL_3D]][T22];
  q[ULL_3D*10 + Q223] = chi*(2.0*dTdy[compy[ULL_3D]][T23]
    + dTdz[compz[ULL_3D]][T22])/3.0;
  q[ULL_3D*10 + Q233] = chi*(dTdy[compy[ULL_3D]][T33]
    + 2.0*dTdz[compz[ULL_3D]][T23])/3.0;
  q[ULL_3D*10 + Q333] = chi*dTdz[compz[ULL_3D]][T33];

  q[ULU_3D*10 + Q111] = chi*dTdx[compx[ULU_3D]][T11];
  q[ULU_3D*10 + Q112] = chi*(2.0*dTdx[compx[ULU_3D]][T12]
    + dTdy[compy[ULU_3D]][T11])/3.0;
  q[ULU_3D*10 + Q113] = chi*(2.0*dTdx[compx[ULU_3D]][T13]
    + dTdz[compz[ULU_3D]][T11])/3.0;
  q[ULU_3D*10 + Q122] = chi*(dTdx[compx[ULU_3D]][T22]
    + 2.0*dTdy[compy[ULU_3D]][T12])/3.0;
  q[ULU_3D*10 + Q123] = chi*(dTdx[compx[ULU_3D]][T23]
    + dTdy[compy[ULU_3D]][T13] + dTdz[compz[ULU_3D]][T12])/3.0;
  q[ULU_3D*10 + Q133] = chi*(dTdx[compx[ULU_3D]][T33]
    + 2.0*dTdz[compz[ULU_3D]][T13])/3.0;
  q[ULU_3D*10 + Q222] = chi*dTdy[compy[ULU_3D]][T22];
  q[ULU_3D*10 + Q223] = chi*(2.0*dTdy[compy[ULU_3D]][T23]
    + dTdz[compz[ULU_3D]][T22])/3.0;
  q[ULU_3D*10 + Q233] = chi*(dTdy[compy[ULU_3D]][T33]
    + 2.0*dTdz[compz[ULU_3D]][T23])/3.0;
  q[ULU_3D*10 + Q333] = chi*dTdz[compz[ULU_3D]][T33];

  q[UUL_3D*10 + Q111] = chi*dTdx[compx[UUL_3D]][T11];
  q[UUL_3D*10 + Q112] = chi*(2.0*dTdx[compx[UUL_3D]][T12]
    + dTdy[compy[UUL_3D]][T11])/3.0;
  q[UUL_3D*10 + Q113] = chi*(2.0*dTdx[compx[UUL_3D]][T13]
    + dTdz[compz[UUL_3D]][T11])/3.0;
  q[UUL_3D*10 + Q122] = chi*(dTdx[compx[UUL_3D]][T22]
    + 2.0*dTdy[compy[UUL_3D]][T12])/3.0;
  q[UUL_3D*10 + Q123] = chi*(dTdx[compx[UUL_3D]][T23]
    + dTdy[compy[UUL_3D]][T13] + dTdz[compz[UUL_3D]][T12])/3.0;
  q[UUL_3D*10 + Q133] = chi*(dTdx[compx[UUL_3D]][T33]
    + 2.0*dTdz[compz[UUL_3D]][T13])/3.0;
  q[UUL_3D*10 + Q222] = chi*dTdy[compy[UUL_3D]][T22];
  q[UUL_3D*10 + Q223] = chi*(2.0*dTdy[compy[UUL_3D]][T23]
    + dTdz[compz[UUL_3D]][T22])/3.0;
  q[UUL_3D*10 + Q233] = chi*(dTdy[compy[UUL_3D]][T33]
    + 2.0*dTdz[compz[UUL_3D]][T23])/3.0;
  q[UUL_3D*10 + Q333] = chi*dTdz[compz[UUL_3D]][T33];

  q[UUU_3D*10 + Q111] = chi*dTdx[compx[UUU_3D]][T11];
  q[UUU_3D*10 + Q112] = chi*(2.0*dTdx[compx[UUU_3D]][T12]
    + dTdy[compy[UUU_3D]][T11])/3.0;
  q[UUU_3D*10 + Q113] = chi*(2.0*dTdx[compx[UUU_3D]][T13]
    + dTdz[compz[UUU_3D]][T11])/3.0;
  q[UUU_3D*10 + Q122] = chi*(dTdx[compx[UUU_3D]][T22]
    + 2.0*dTdy[compy[UUU_3D]][T12])/3.0;
  q[UUU_3D*10 + Q123] = chi*(dTdx[compx[UUU_3D]][T23]
    + dTdy[compy[UUU_3D]][T13] + dTdz[compz[UUU_3D]][T12])/3.0;
  q[UUU_3D*10 + Q133] = chi*(dTdx[compx[UUU_3D]][T33]
    + 2.0*dTdz[compz[UUU_3D]][T13])/3.0;
  q[UUU_3D*10 + Q222] = chi*dTdy[compy[UUU_3D]][T22];
  q[UUU_3D*10 + Q223] = chi*(2.0*dTdy[compy[UUU_3D]][T23]
    + dTdz[compz[UUU_3D]][T22])/3.0;
  q[UUU_3D*10 + Q233] = chi*(dTdy[compy[UUU_3D]][T33]
    + 2.0*dTdz[compz[UUU_3D]][T23])/3.0;
  q[UUU_3D*10 + Q333] = chi*dTdz[compz[UUU_3D]][T33];

  double da = fmin(fmin(dx, dy), dz);
  double cfla = dt/(da*da);
  cflrate[0] = alpha*vth_avg*cfla;
}

GKYL_CU_D
static void
grad_closure_update_3d(const gkyl_ten_moment_grad_closure *gces,
  const double *q[], double *rhs)
{
  double div_qx[6] = {0.0};
  double div_qy[6] = {0.0};
  double div_qz[6] = {0.0};

  const double dx = gces->grid.dx[0];
  const double dy = gces->grid.dx[1];
  const double dz = gces->grid.dx[2];

  div_qx[0] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q111],
    q[LLU_3D][UUL_3D*10 + Q111], q[LUL_3D][ULU_3D*10 + Q111],
    q[LUU_3D][ULL_3D*10 + Q111], q[ULL_3D][LUU_3D*10 + Q111],
    q[ULU_3D][LUL_3D*10 + Q111], q[UUL_3D][LLU_3D*10 + Q111],
    q[UUU_3D][LLL_3D*10 + Q111]);
  div_qx[1] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q112],
    q[LLU_3D][UUL_3D*10 + Q112], q[LUL_3D][ULU_3D*10 + Q112],
    q[LUU_3D][ULL_3D*10 + Q112], q[ULL_3D][LUU_3D*10 + Q112],
    q[ULU_3D][LUL_3D*10 + Q112], q[UUL_3D][LLU_3D*10 + Q112],
    q[UUU_3D][LLL_3D*10 + Q112]);
  div_qx[2] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q113],
    q[LLU_3D][UUL_3D*10 + Q113], q[LUL_3D][ULU_3D*10 + Q113],
    q[LUU_3D][ULL_3D*10 + Q113], q[ULL_3D][LUU_3D*10 + Q113],
    q[ULU_3D][LUL_3D*10 + Q113], q[UUL_3D][LLU_3D*10 + Q113],
    q[UUU_3D][LLL_3D*10 + Q113]);
  div_qx[3] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q122],
    q[LLU_3D][UUL_3D*10 + Q122], q[LUL_3D][ULU_3D*10 + Q122],
    q[LUU_3D][ULL_3D*10 + Q122], q[ULL_3D][LUU_3D*10 + Q122],
    q[ULU_3D][LUL_3D*10 + Q122], q[UUL_3D][LLU_3D*10 + Q122],
    q[UUU_3D][LLL_3D*10 + Q122]);
  div_qx[4] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q123],
    q[LLU_3D][UUL_3D*10 + Q123], q[LUL_3D][ULU_3D*10 + Q123],
    q[LUU_3D][ULL_3D*10 + Q123], q[ULL_3D][LUU_3D*10 + Q123],
    q[ULU_3D][LUL_3D*10 + Q123], q[UUL_3D][LLU_3D*10 + Q123],
    q[UUU_3D][LLL_3D*10 + Q123]);
  div_qx[5] = calc_sym_gradx_3D(dx, q[LLL_3D][UUU_3D*10 + Q133],
    q[LLU_3D][UUL_3D*10 + Q133], q[LUL_3D][ULU_3D*10 + Q133],
    q[LUU_3D][ULL_3D*10 + Q133], q[ULL_3D][LUU_3D*10 + Q133],
    q[ULU_3D][LUL_3D*10 + Q133], q[UUL_3D][LLU_3D*10 + Q133],
    q[UUU_3D][LLL_3D*10 + Q133]);

  div_qy[0] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q112],
    q[LLU_3D][UUL_3D*10 + Q112], q[LUL_3D][ULU_3D*10 + Q112],
    q[LUU_3D][ULL_3D*10 + Q112], q[ULL_3D][LUU_3D*10 + Q112],
    q[ULU_3D][LUL_3D*10 + Q112], q[UUL_3D][LLU_3D*10 + Q112],
    q[UUU_3D][LLL_3D*10 + Q112]);
  div_qy[1] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q122],
    q[LLU_3D][UUL_3D*10 + Q122], q[LUL_3D][ULU_3D*10 + Q122],
    q[LUU_3D][ULL_3D*10 + Q122], q[ULL_3D][LUU_3D*10 + Q122],
    q[ULU_3D][LUL_3D*10 + Q122], q[UUL_3D][LLU_3D*10 + Q122],
    q[UUU_3D][LLL_3D*10 + Q122]);
  div_qy[2] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q123],
    q[LLU_3D][UUL_3D*10 + Q123], q[LUL_3D][ULU_3D*10 + Q123],
    q[LUU_3D][ULL_3D*10 + Q123], q[ULL_3D][LUU_3D*10 + Q123],
    q[ULU_3D][LUL_3D*10 + Q123], q[UUL_3D][LLU_3D*10 + Q123],
    q[UUU_3D][LLL_3D*10 + Q123]);
  div_qy[3] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q222],
    q[LLU_3D][UUL_3D*10 + Q222], q[LUL_3D][ULU_3D*10 + Q222],
    q[LUU_3D][ULL_3D*10 + Q222], q[ULL_3D][LUU_3D*10 + Q222],
    q[ULU_3D][LUL_3D*10 + Q222], q[UUL_3D][LLU_3D*10 + Q222],
    q[UUU_3D][LLL_3D*10 + Q222]);
  div_qy[4] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q223],
    q[LLU_3D][UUL_3D*10 + Q223], q[LUL_3D][ULU_3D*10 + Q223],
    q[LUU_3D][ULL_3D*10 + Q223], q[ULL_3D][LUU_3D*10 + Q223],
    q[ULU_3D][LUL_3D*10 + Q223], q[UUL_3D][LLU_3D*10 + Q223],
    q[UUU_3D][LLL_3D*10 + Q223]);
  div_qy[5] = calc_sym_grady_3D(dy, q[LLL_3D][UUU_3D*10 + Q233],
    q[LLU_3D][UUL_3D*10 + Q233], q[LUL_3D][ULU_3D*10 + Q233],
    q[LUU_3D][ULL_3D*10 + Q233], q[ULL_3D][LUU_3D*10 + Q233],
    q[ULU_3D][LUL_3D*10 + Q233], q[UUL_3D][LLU_3D*10 + Q233],
    q[UUU_3D][LLL_3D*10 + Q233]);

  div_qz[0] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q113],
    q[LLU_3D][UUL_3D*10 + Q113], q[LUL_3D][ULU_3D*10 + Q113],
    q[LUU_3D][ULL_3D*10 + Q113], q[ULL_3D][LUU_3D*10 + Q113],
    q[ULU_3D][LUL_3D*10 + Q113], q[UUL_3D][LLU_3D*10 + Q113],
    q[UUU_3D][LLL_3D*10 + Q113]);
  div_qz[1] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q123],
    q[LLU_3D][UUL_3D*10 + Q123], q[LUL_3D][ULU_3D*10 + Q123],
    q[LUU_3D][ULL_3D*10 + Q123], q[ULL_3D][LUU_3D*10 + Q123],
    q[ULU_3D][LUL_3D*10 + Q123], q[UUL_3D][LLU_3D*10 + Q123],
    q[UUU_3D][LLL_3D*10 + Q123]);
  div_qz[2] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q133],
    q[LLU_3D][UUL_3D*10 + Q133], q[LUL_3D][ULU_3D*10 + Q133],
    q[LUU_3D][ULL_3D*10 + Q133], q[ULL_3D][LUU_3D*10 + Q133],
    q[ULU_3D][LUL_3D*10 + Q133], q[UUL_3D][LLU_3D*10 + Q133],
    q[UUU_3D][LLL_3D*10 + Q133]);
  div_qz[3] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q223],
    q[LLU_3D][UUL_3D*10 + Q223], q[LUL_3D][ULU_3D*10 + Q223],
    q[LUU_3D][ULL_3D*10 + Q223], q[ULL_3D][LUU_3D*10 + Q223],
    q[ULU_3D][LUL_3D*10 + Q223], q[UUL_3D][LLU_3D*10 + Q223],
    q[UUU_3D][LLL_3D*10 + Q223]);
  div_qz[4] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q233],
    q[LLU_3D][UUL_3D*10 + Q233], q[LUL_3D][ULU_3D*10 + Q233],
    q[LUU_3D][ULL_3D*10 + Q233], q[ULL_3D][LUU_3D*10 + Q233],
    q[ULU_3D][LUL_3D*10 + Q233], q[UUL_3D][LLU_3D*10 + Q233],
    q[UUU_3D][LLL_3D*10 + Q233]);
  div_qz[5] = calc_sym_gradz_3D(dz, q[LLL_3D][UUU_3D*10 + Q333],
    q[LLU_3D][UUL_3D*10 + Q333], q[LUL_3D][ULU_3D*10 + Q333],
    q[LUU_3D][ULL_3D*10 + Q333], q[ULL_3D][LUU_3D*10 + Q333],
    q[ULU_3D][LUL_3D*10 + Q333], q[UUL_3D][LLU_3D*10 + Q333],
    q[UUU_3D][LLL_3D*10 + Q333]);

  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = div_qx[0] + div_qy[0] + div_qz[0];
  rhs[P12] = div_qx[1] + div_qy[1] + div_qz[1];
  rhs[P13] = div_qx[2] + div_qy[2] + div_qz[2];
  rhs[P22] = div_qx[3] + div_qy[3] + div_qz[3];
  rhs[P23] = div_qx[4] + div_qy[4] + div_qz[4];
  rhs[P33] = div_qx[5] + div_qy[5] + div_qz[5];
}

GKYL_CU_D
static const heat_flux_calc_t grad_closure_unmag_funcs[3] = { calc_unmag_heat_flux_1d,
  calc_unmag_heat_flux_2d, calc_unmag_heat_flux_3d };

GKYL_CU_D
static const heat_flux_update_t grad_closure_update_funcs[3] = { grad_closure_update_1d,
  grad_closure_update_2d, grad_closure_update_3d };

GKYL_CU_D
static void
grad_closure_calc_q_choose(gkyl_ten_moment_grad_closure *gces)
{
  gces->calc_q = grad_closure_unmag_funcs[gces->ndim - 1];
}

GKYL_CU_D
static void
grad_closure_update_q_choose(gkyl_ten_moment_grad_closure *gces)
{
  gces->update_q = grad_closure_update_funcs[gces->ndim - 1];
}
