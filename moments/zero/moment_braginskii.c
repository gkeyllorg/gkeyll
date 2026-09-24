#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_non_ideal_priv.h>
#include <gkyl_wv_euler_priv.h>

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

struct gkyl_moment_braginskii {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int nfluids; // number of fluids in multi-fluid system
  struct gkyl_moment_braginskii_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters
  double epsilon0; // permittivity of free space
  double coll_fac; // constant multiplicative factor for collision time to increase or decrease collisionality
};

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

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
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

// Fetch input quantities and compute derived quantities for magnetized Braginskii
static void
mag_var_setup(const gkyl_moment_braginskii *bes, int start, int end,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double u[][2][3], double b[][3],
  double T[][2], double tau[][2],
  double eta_par[][2], double eta_perp[][2], double eta_cross[][2],
  double kappa_par[][2], double kappa_perp[][2], double kappa_cross[][2],
  double current_par[][3], double current_cross[][3])
{
  int nfluids = bes->nfluids;

  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };

  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double rho[2] = {};     // Mass density for each species
  double p[2] = {};       // Pressure for each species
  double omega_c[2] = {}; // Cyclotron frequency for each species

  for (int j = start; j <= end; ++j) {
    for (int n = 0; n < nfluids; ++n) {
      // Input density and flow in each cell
      rho[n] = fluid_d[j][n][RHO];
      for (int k = 0; k < 3; ++k)
        u[j][n][k] = fluid_d[j][n][MX + k] / rho[n];

      // Pressure information is different for each equation type
      if (bes->param[n].type_eqn == GKYL_EQN_EULER)
        p[n] = gkyl_euler_pressure(bes->param[n].p_fac, fluid_d[j][n]); // Euler needs to divide out gas_gamma factor to obtain pressure
      else if (bes->param[n].type_eqn == GKYL_EQN_ISO_EULER)
        p[n] = rho[n] * bes->param[n].p_fac * bes->param[n].p_fac; // isothermal Euler input is vth, pressure = rho*vth^2

      T[j][n] = m[n] * p[n] / rho[n];

      omega_c[n] = calc_omega_c(q[n], m[n], em_tot_d[j]);
    }

    // Input magnetic field in each cell
    calc_bhat(em_tot_d[j], b[j]);

    // Derived collision times
    tau[j][ELC] = calc_tau(1.0, bes->coll_fac, bes->epsilon0, q[ELC], q[ION], m[ELC], m[ION], rho[ION], T[j][ELC]);
    tau[j][ION] = calc_tau(1.0, sqrt(2.0)*bes->coll_fac, bes->epsilon0, q[ION], q[ION], m[ION], m[ION], rho[ION], T[j][ION]);

    // Brag-type enum is used to turn coefficients on/off in a branchless fashion
    bool electron_viscosity = (bes->param[ELC].type_brag & GKYL_BRAG_VISC);
    bool electron_heatFlux = (bes->param[ELC].type_brag & GKYL_BRAG_HEATFLUX);
    bool ion_viscosity = (bes->param[ION].type_brag & GKYL_BRAG_VISC);
    bool ion_heatFlux = (bes->param[ION].type_brag & GKYL_BRAG_HEATFLUX);

    eta_par[j][ELC] = electron_viscosity * 1.5*0.73*p[ELC]*tau[j][ELC];
    eta_perp[j][ELC] = electron_viscosity * 0.51*p[ELC]/(tau[j][ELC]*omega_c[ELC]*omega_c[ELC]);
    eta_cross[j][ELC] = electron_viscosity * 0.25*p[ELC]/omega_c[ELC];
    kappa_par[j][ELC] = electron_heatFlux * 3.16*p[ELC]*tau[j][ELC]/m[ELC];
    kappa_perp[j][ELC] = electron_heatFlux * 4.66*p[ELC]/(m[ELC]*tau[j][ELC]*omega_c[ELC]*omega_c[ELC]);
    kappa_cross[j][ELC] = electron_heatFlux * 2.5*p[ELC]/(m[ELC]*omega_c[ELC]);

    eta_par[j][ION] = ion_viscosity * 1.5*0.96*p[ION]*tau[j][ION];
    eta_perp[j][ION] = ion_viscosity * 0.3*p[ION]/(tau[j][ION]*omega_c[ION]*omega_c[ION]);
    eta_cross[j][ION] = ion_viscosity * 0.25*p[ION]/omega_c[ION];
    kappa_par[j][ION] = ion_heatFlux * 3.91*p[ION]*tau[j][ION]/m[ION];
    kappa_perp[j][ION] = ion_heatFlux * 2*p[ION]/(m[ION]*tau[j][ION]*omega_c[ION]*omega_c[ION]);
    kappa_cross[j][ION] = ion_heatFlux * 2.5*p[ION]/(m[ION]*omega_c[ION]);

    double thermal_par = electron_heatFlux * -0.71*p[ELC]; // Parallel thermal force coefficient (same for each species)
    double thermal_perp = electron_heatFlux * 1.5*p[ELC]/(omega_c[ELC]*tau[j][ELC]); // Perpendicular thermal force coefficient (same for each species)
    double b_dot_j = b[j][0]*(u[j][ION][0] - u[j][ELC][0]) + b[j][1]*(u[j][ION][1] - u[j][ELC][1]) + b[j][2]*(u[j][ION][2] - u[j][ELC][2]);
    current_par[j][0] = thermal_par*b[j][0]*b_dot_j;
    current_par[j][1] = thermal_par*b[j][1]*b_dot_j;
    current_par[j][2] = thermal_par*b[j][2]*b_dot_j;
    current_cross[j][0] = thermal_perp*(b[j][1]*(u[j][ION][2] - u[j][ELC][2]) - b[j][2]*(u[j][ION][1] - u[j][ELC][1]));
    current_cross[j][1] = thermal_perp*(b[j][2]*(u[j][ION][0] - u[j][ELC][0]) - b[j][0]*(u[j][ION][2] - u[j][ELC][2]));
    current_cross[j][2] = thermal_perp*(b[j][0]*(u[j][ION][1] - u[j][ELC][1]) - b[j][1]*(u[j][ION][0] - u[j][ELC][0]));
  }
}

static void
mag_brag_calc_vars(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double *cflrate[GKYL_MAX_SPECIES], double *brag_d[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;

  // Allocate some memory on stack regardless of dimensionality or equation type (Euler vs. isothermal Euler)
  // These allocations allow us to make final construction of Braginskii variables dimensionally independent
  // Note: Braginskii implementation currently assumes only two fluids (electrons and an ion species)
  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };

  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double b_avg[3] = {};
  double eta_par_avg[2] = {}, eta_perp_avg[2] = {}, eta_cross_avg[2] = {};
  double w[2][6] = {};

  double pi_par[2][6] = {};
  double pi_perp[2][6] = {};
  double pi_cross[2][6] = {};

  double kappa_par_avg[2] = {}, kappa_perp_avg[2] = {}, kappa_cross_avg[2] = {};
  double current_par_avg[3] = {};   // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
  double current_cross_avg[3] = {}; // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)

  // Temperature gradient, parallel, perp, and cross at cell edges
  double gradxT[2] = {};
  double gradyT[2] = {};
  double gradzT[2] = {};

  double bbgradT[2][3] = {};     // Parallel temperature gradient (b_hat b_hat dot grad T)
  double perp_gradT[2][3] = {};  // Perpendicular temperature gradient (grad T - b_hat b_hat dot grad T)
  double cross_gradT[2][3] = {}; // Cross temperature gradient (b x grad T)

  // Average velocity at cell edges
  double u_avg[2][3] = {0.0};

  if (ndim == 1) {
    double dx = bes->grid.dx[0];

    // Derived quantities
    double u[2][2][3] = {};           // Flow for each species, ux, uy, & uz
    double b[2][3] = {};              // Magnetic field unit vector, Bx/|B|, By/|B|, & Bz/|B|
    double T[2][2] = {};              // Temperature for each species (m*p/rho)
    double tau[2][2] = {};            // Collision times for each species
    double eta_par[2][2] = {};        // Parallel viscosity for each species
    double eta_perp[2][2] = {};       // Perpendicular viscosity for each species
    double eta_cross[2][2] = {};      // Gyro-viscosity for each species
    double kappa_par[2][2] = {};      // Parallel conductivity for each species
    double kappa_perp[2][2] = {};     // Perpendicular conductivity for each species
    double kappa_cross[2][2] = {};    // Gyro-conductivity for each species
    double current_par[2][3] = {};    // b_hat b_hat dot (u_i - u_e) 
    double current_cross[2][3] = {};  // b x (u_i - u_e)

    // Compute derived quantities
    mag_var_setup(bes, L_1D, U_1D,
      fluid_d, em_tot_d,
      u, b, T, tau,
      eta_par, eta_perp, eta_cross,
      kappa_par, kappa_perp, kappa_cross,
      current_par, current_cross);

    // Magnetic field at cell edges (using arithmetic average)
    for (int k = 0; k < 3; ++k)
      b_avg[k] = calc_arithm_avg_1D(b[L_1D][k], b[U_1D][k]);
    
    for (int n = 0; n < nfluids; ++n) {
      // Parallel viscosity, perpendicular viscosity, and gyro-viscosity coefficients at cell edges (using harmonic average)
      eta_par_avg[n] = calc_harmonic_avg_1D(eta_par[L_1D][n], eta_par[U_1D][n]);
      eta_perp_avg[n] = calc_harmonic_avg_1D(eta_perp[L_1D][n], eta_perp[U_1D][n]);
      eta_cross_avg[n] = calc_harmonic_avg_1D(eta_cross[L_1D][n], eta_cross[U_1D][n]);

      // Rate of strain tensor at cell edges
      calc_ros_1D(dx, u[L_1D][n], u[U_1D][n], w[n]);

      // Calculate heat flux and viscous heating if energy variable exists
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        // Parallel conductivity, perpendicular conductivity, and gyro-conductivity coefficients at cell edges (using harmonic average)
        kappa_par_avg[n] = calc_harmonic_avg_1D(kappa_par[L_1D][n], kappa_par[U_1D][n]);
        kappa_perp_avg[n] = calc_harmonic_avg_1D(kappa_perp[L_1D][n], kappa_perp[U_1D][n]);
        kappa_cross_avg[n] = calc_harmonic_avg_1D(kappa_cross[L_1D][n], kappa_cross[U_1D][n]);

        gradxT[n] = calc_sym_grad_1D(dx, T[L_1D][n], T[U_1D][n]);

        for (int k = 0; k < 3; ++k)
          bbgradT[n][k] = b_avg[k] * b_avg[0] * gradxT[n];

        perp_gradT[n][0] = gradxT[n] - bbgradT[n][0];

        cross_gradT[n][1] = b_avg[2] * gradxT[n];
        cross_gradT[n][2] = -b_avg[1] * gradxT[n];

        // Average velocity at cell edges (using arithmetic average)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_1D(u[L_1D][n][k], u[U_1D][n][k]);
      }
    }
    
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) 
      for (int k = 0; k < 3; ++k) {
        // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
        current_par_avg[k] = calc_arithm_avg_1D(current_par[L_1D][k], current_par[U_1D][k]);
        // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)
        current_cross_avg[k] = calc_arithm_avg_1D(current_cross[L_1D][k], current_cross[U_1D][k]);
      }
  }
  else if (ndim == 2) {
    double dx = bes->grid.dx[0];
    double dy = bes->grid.dx[1];

    // Derived quantities
    double u[4][2][3] = {};          // Flow for each species, ux, uy, & uz
    double b[4][3] = {};             // Magnetic field unit vector, Bx/|B|, By/|B|, & Bz/|B|
    double T[4][2] = {};             // Temperature for each species (m*p/rho)
    double tau[4][2] = {};           // Collision times for each species
    double eta_par[4][2] = {};       // Parallel viscosity for each species
    double eta_perp[4][2] = {};      // Perpendicular viscosity for each species
    double eta_cross[4][2] = {};     // Gyro-viscosity for each species
    double kappa_par[4][2] = {};     // Parallel conductivity for each species
    double kappa_perp[4][2] = {};    // Perpendicular conductivity for each species
    double kappa_cross[4][2] = {};   // Gyro-conductivity for each species
    double current_par[4][3] = {};   // b_hat b_hat dot (u_i - u_e) 
    double current_cross[4][3] = {}; // b x (u_i - u_e)

    // Compute derived quantities
    mag_var_setup(bes, LL_2D, UU_2D,
      fluid_d, em_tot_d,
      u, b, T, tau,
      eta_par, eta_perp, eta_cross,
      kappa_par, kappa_perp, kappa_cross,
      current_par, current_cross);

    // Magnetic field at cell vertices (using arithmetic average)
    for (int k = 0; k < 3; ++k)
      b_avg[k] = calc_arithm_avg_2D(b[LL_2D][k], b[LU_2D][k], b[UL_2D][k], b[UU_2D][k]);

    for (int n = 0; n < nfluids; ++n) {
      // Parallel viscosity, perpendicular viscosity, and gyro-viscosity coefficients at cell edges (using harmonic average)
      eta_par_avg[n] = calc_harmonic_avg_2D(eta_par[LL_2D][n], eta_par[LU_2D][n], eta_par[UL_2D][n], eta_par[UU_2D][n]);
      eta_perp_avg[n] = calc_harmonic_avg_2D(eta_perp[LL_2D][n], eta_perp[LU_2D][n], eta_perp[UL_2D][n], eta_perp[UU_2D][n]);
      eta_cross_avg[n] = calc_harmonic_avg_2D(eta_cross[LL_2D][n], eta_cross[LU_2D][n], eta_cross[UL_2D][n], eta_cross[UU_2D][n]);

      // Rate of strain tensor at cell vertices for electrons and ions
      calc_ros_2D(dx, dy, u[LL_2D][n], u[LU_2D][n], u[UL_2D][n], u[UU_2D][n], w[n]);

      // Calculate heat flux and viscous heating if energy variable exists
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        // Parallel conductivity, perpendicular conductivity, and gyro-conductivity coefficients at cell edges (using harmonic average)
        kappa_par_avg[n] = calc_harmonic_avg_2D(kappa_par[LL_2D][n], kappa_par[LU_2D][n], kappa_par[UL_2D][n], kappa_par[UU_2D][n]);
        kappa_perp_avg[n] = calc_harmonic_avg_2D(kappa_perp[LL_2D][n], kappa_perp[LU_2D][n], kappa_perp[UL_2D][n], kappa_perp[UU_2D][n]);
        kappa_cross_avg[n] = calc_harmonic_avg_2D(kappa_cross[LL_2D][n], kappa_cross[LU_2D][n], kappa_cross[UL_2D][n], kappa_cross[UU_2D][n]);

        gradxT[n] = calc_sym_gradx_2D(dx, T[LL_2D][n], T[LU_2D][n], T[UL_2D][n], T[UU_2D][n]);
        gradyT[n] = calc_sym_grady_2D(dy, T[LL_2D][n], T[LU_2D][n], T[UL_2D][n], T[UU_2D][n]);

        for (int k = 0; k < 3; ++k)
          bbgradT[n][k] = b_avg[k] * (b_avg[0]*gradxT[n] + b_avg[1]*gradyT[n]);

        perp_gradT[n][0] = gradxT[n] - bbgradT[n][0];
        perp_gradT[n][1] = gradyT[n] - bbgradT[n][1];

        cross_gradT[n][0] = -b_avg[2]*gradyT[n];
        cross_gradT[n][1] = b_avg[2]*gradxT[n];
        cross_gradT[n][2] = b_avg[0]*gradyT[n] - b_avg[1]*gradxT[n];

        // Average velocity at cell edges (using arithmetic average)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_2D(u[LL_2D][n][k], u[LU_2D][n][k], u[UL_2D][n][k], u[UU_2D][n][k]);
      }
    }
    
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) 
      for (int k = 0; k < 3; ++k) {
        // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
        current_par_avg[k] = calc_arithm_avg_2D(current_par[LL_2D][k], current_par[LU_2D][k], current_par[UL_2D][k], current_par[UU_2D][k]);
        // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)
        current_cross_avg[k] = calc_arithm_avg_2D(current_cross[LL_2D][k], current_cross[LU_2D][k], current_cross[UL_2D][k], current_cross[UU_2D][k]);
      }
  }

  for (int n = 0; n < nfluids; ++n) {
    // Parallel, perpendicular, and gyro-viscous stree tensors at cell edges
    calc_pi_par(eta_par_avg[n], b_avg, w[n], pi_par[n]);
    calc_pi_perp(eta_perp_avg[n], b_avg, w[n], pi_perp[n]);
    calc_pi_cross(eta_cross_avg[n], b_avg, w[n], pi_cross[n]);

    // Total viscous stress tensor
    for (int k = 0; k < 6; ++k)
      brag_d[n][PIXX + k] = pi_par[n][k] + pi_perp[n][k] + pi_cross[n][k];

    // Pi dot u
    double Piu[3] = { brag_d[n][PIXX]*u_avg[n][0] + brag_d[n][PIXY]*u_avg[n][1] + brag_d[n][PIXZ]*u_avg[n][2],
                      brag_d[n][PIXY]*u_avg[n][0] + brag_d[n][PIYY]*u_avg[n][1] + brag_d[n][PIYZ]*u_avg[n][2],
                      brag_d[n][PIXZ]*u_avg[n][0] + brag_d[n][PIYZ]*u_avg[n][1] + brag_d[n][PIZZ]*u_avg[n][2] };

    // Total heat flux + viscous heating
    for (int k = 0; k < 3; ++k)
      brag_d[n][QX + k] = -kappa_par_avg[n]*bbgradT[n][k] - kappa_perp_avg[n]*perp_gradT[n][k] + kappa_cross_avg[n]*cross_gradT[n][k] + Piu[k];
  }

  for (int k = 0; k < 3; ++k)
    brag_d[ELC][QX + k] += current_par_avg[k] + current_cross_avg[k];
}

// Fetch input quantities and compute derived quantities for UNmagnetized Braginskii
static void
unmag_var_setup(const gkyl_moment_braginskii *bes, int start, int end,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double u[][2][3],
  double T[][2], double tau[][2],
  double eta[][2], double kappa[][2],
  double current[][3])
{
  int nfluids = bes->nfluids;

  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };

  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double rho[2] = {};   // Mass density for each species
  double p[2] = {};     // Pressure for each species

  for (int j = start; j <= end; ++j) {
    for (int n = 0; n < nfluids; ++n) {
      // Input density and flow in each cell
      rho[n] = fluid_d[j][n][RHO];
      for (int k = 0; k < 3; ++k)
        u[j][n][k] = fluid_d[j][n][MX + k] / rho[n];

      // Pressure information is different for each equation type
      if (bes->param[n].type_eqn == GKYL_EQN_EULER)
        p[n] = gkyl_euler_pressure(bes->param[n].p_fac, fluid_d[j][n]); // Euler needs to divide out gas_gamma factor to obtain pressure
      else if (bes->param[n].type_eqn == GKYL_EQN_ISO_EULER)
        p[n] = rho[n] * bes->param[n].p_fac * bes->param[n].p_fac; // isothermal Euler input is vth, pressure = rho*vth^2

      T[j][n] = m[n] * p[n] / rho[n];
    }

    tau[j][ELC] = calc_tau(1.0, bes->coll_fac, bes->epsilon0, q[ELC], q[ION], m[ELC], m[ION], rho[ION], T[j][ELC]);
    tau[j][ION] = calc_tau(1.0, sqrt(2.0)*bes->coll_fac, bes->epsilon0, q[ION], q[ION], m[ION], m[ION], rho[ION], T[j][ION]);

    // Brag-type enum is used to turn coefficients on/off in a branchless fashion
    bool electron_viscosity = (bes->param[ELC].type_brag & GKYL_BRAG_VISC);
    bool electron_heatFlux = (bes->param[ELC].type_brag & GKYL_BRAG_HEATFLUX);
    bool ion_viscosity = (bes->param[ION].type_brag & GKYL_BRAG_VISC);
    bool ion_heatFlux = (bes->param[ION].type_brag & GKYL_BRAG_HEATFLUX);

    eta[j][ELC] = electron_viscosity * 0.73*p[ELC]*tau[j][ELC];
    eta[j][ION] = ion_viscosity * 0.96*p[ION]*tau[j][ION];

    kappa[j][ELC] = electron_heatFlux * 3.16*p[ELC]*tau[j][ELC]/m[ELC];
    kappa[j][ION] = ion_heatFlux * 3.91*p[ION]*tau[j][ION]/m[ION];

    double thermal = electron_heatFlux * -0.71*p[ELC]; // thermal force coefficient (same for each species)
    for (int k = 0; k < 3; ++k)
      current[j][k] = thermal*(u[j][ION][k] - u[j][ELC][k]);
  }
}
    
static void
unmag_brag_calc_vars(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double *cflrate[GKYL_MAX_SPECIES], double *brag_d[GKYL_MAX_SPECIES])
{
  const int nfluids = bes->nfluids;
  const int ndim = bes->ndim;

  // Allocate some memory on stack regardless of dimensionality or equation type (Euler vs. isothermal Euler)
  // These allocations allow us to make final construction of Braginskii variables dimensionally independent
  // Note: Braginskii implementation currently assumes only two fluids (electrons and an ion species)
  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };

  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double eta_avg[2] = {};
  double kappa_avg[2] = {};
  double w[2][6] = {};

  // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
  double current_avg[3] = {};

  // Temperature gradient at cell edges
  double gradxT[2] = {};
  double gradyT[2] = {};
  double gradzT[2] = {};

  // Average velocity at cell edges
  double u_avg[2][3] = {};

  if (ndim == 1) {
    double dx = bes->grid.dx[0]; 

    // Derived quantities
    double u[2][2][3] = {};     // Flow for each species, ux, uy, & uz
    double T[2][2] = {};        // Temperature for each species (m*p/rho)
    double tau[2][2] = {};      // Collision times for each species
    double eta[2][2] = {};      // Viscosity for each species
    double kappa[2][2] = {};    // Conductivity for each species
    double current[2][3] = {};  // (u_i - u_e) 

    // Compute derived quantities
    unmag_var_setup(bes, L_1D, U_1D, fluid_d,
      u, T, tau,
      eta, kappa, current);

    for (int n = 0; n < nfluids; ++n) {
      // Viscosity coefficients at cell edges (using harmonic average)
      eta_avg[n] = calc_harmonic_avg_1D(eta[L_1D][n], eta[U_1D][n]);

      // Rate of strain tensor at cell edges for electrons and ions
      calc_ros_1D(dx, u[L_1D][n], u[U_1D][n], w[n]);

      // Calculate heat flux and viscous heating if energy variable exists
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        // Conductivity coefficients at cell edges (using harmonic average)
        kappa_avg[n] = calc_harmonic_avg_1D(kappa[L_1D][n], kappa[U_1D][n]);

        gradxT[n] = calc_sym_grad_1D(dx, T[L_1D][n], T[U_1D][n]);

        // Average velocity at cell edges (using arithmetic average)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_1D(u[L_1D][n][k], u[U_1D][n][k]);
      }
    }

    // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER)
      for (int k = 0; k < 3; ++k)
        current_avg[k] = calc_arithm_avg_1D(current[L_1D][k], current[U_1D][k]);
  }
  else if (ndim == 2) {
    double dx = bes->grid.dx[0];
    double dy = bes->grid.dx[1]; 

    // Derived quantities
    double u[4][2][3] = {};     // Flow for each species, ux, uy, & uz
    double T[4][2] = {};        // Temperature for each species (m*p/rho)
    double tau[4][2] = {};      // Collision times for each species
    double eta[4][2] = {};      // Viscosity for each species
    double kappa[4][2] = {};    // Conductivity for each species
    double current[4][3] = {};  // (u_i - u_e)

    // Compute derived quantities
    unmag_var_setup(bes, LL_2D, UU_2D, fluid_d,
      u, T, tau,
      eta, kappa, current);

    for (int n = 0; n < nfluids; ++n) {
      // Viscosity coefficients at cell edges (using harmonic average)
      eta_avg[n] = calc_harmonic_avg_2D(eta[LL_2D][n], eta[LU_2D][n], eta[UL_2D][n], eta[UU_2D][n]);

      // Rate of strain tensor at cell vertices for electrons and ions
      calc_ros_2D(dx, dy, u[LL_2D][n], u[LU_2D][n], u[UL_2D][n], u[UU_2D][n], w[n]);

      // Calculate heat flux and viscous heating if energy variable exists
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        // Conductivity coefficients at cell edges (using harmonic average)
        kappa_avg[n] = calc_harmonic_avg_2D(kappa[LL_2D][n], kappa[LU_2D][n], kappa[UL_2D][n], kappa[UU_2D][n]);

        gradxT[n] = calc_sym_gradx_2D(dx, T[LL_2D][n], T[LU_2D][n], T[UL_2D][n], T[UU_2D][n]);
        gradyT[n] = calc_sym_grady_2D(dy, T[LL_2D][n], T[LU_2D][n], T[UL_2D][n], T[UU_2D][n]);

        // Average velocity at cell edges (using arithmetic average)
        for (int k = 0; k < 3; ++k)
          u_avg[n][k] = calc_arithm_avg_2D(u[LL_2D][n][k], u[LU_2D][n][k], u[UL_2D][n][k], u[UU_2D][n][k]);
      }
    }
    
    // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER)
      for (int k = 0; k < 3; ++k)
        current_avg[k] = calc_arithm_avg_2D(current[LL_2D][k], current[LU_2D][k], current[UL_2D][k], current[UU_2D][k]);
  }

  for (int n = 0; n < nfluids; ++n) {
    // Total viscous stress tensor
    for (int k = 0; k < 6; ++k)
      brag_d[n][PIXX + k] = -eta_avg[n] * w[n][k];

    double gradT[3] = { gradxT[n], gradyT[n] , gradzT[n] };

    // Pi dot u
    double Piu[3] = { brag_d[n][PIXX]*u_avg[n][0] + brag_d[n][PIXY]*u_avg[n][1] + brag_d[n][PIXZ]*u_avg[n][2],
                      brag_d[n][PIXY]*u_avg[n][0] + brag_d[n][PIYY]*u_avg[n][1] + brag_d[n][PIYZ]*u_avg[n][2],
                      brag_d[n][PIXZ]*u_avg[n][0] + brag_d[n][PIYZ]*u_avg[n][1] + brag_d[n][PIZZ]*u_avg[n][2] };

    // Total heat flux + viscous heating
    for (int k = 0; k < 3; ++k)
      brag_d[n][QX + k] = -kappa_avg[n]*gradT[k] + Piu[k];
  }

  for (int k = 0; k < 3; ++k)
    brag_d[ELC][QX + k] += current_avg[k];
}

static void
brag_calc_update(const gkyl_moment_braginskii *bes,
  const double *brag_d[][GKYL_MAX_SPECIES], double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;
  double div_pi[GKYL_MAX_SPECIES][3] = {};
  double div_q[GKYL_MAX_SPECIES] = {};
  if (ndim == 1) {
    const double dx = bes->grid.dx[0];
    double pi[2][GKYL_MAX_SPECIES][6] = {};
    double q[2][GKYL_MAX_SPECIES][3] = {};
    for (int n=0; n < nfluids; ++n) {
      for (int j = L_1D; j <= U_1D; ++j)
        for (int k = 0; k < 6; ++k)
          pi[j][n][k] = brag_d[j][n][PIXX + k];

      div_pi[n][0] = calc_sym_grad_1D(dx, pi[L_1D][n][0], pi[U_1D][n][0]);
      div_pi[n][1] = calc_sym_grad_1D(dx, pi[L_1D][n][1], pi[U_1D][n][1]);
      div_pi[n][2] = calc_sym_grad_1D(dx, pi[L_1D][n][2], pi[U_1D][n][2]);

      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];

      // If energy variable exists, increment heat flux and viscous heating
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = L_1D; j <= U_1D; ++j)
          for (int k = 0; k < 3; ++k)
            q[j][n][k] = brag_d[j][n][QX + k];

        div_q[n] = calc_sym_grad_1D(dx, q[L_1D][n][0], q[U_1D][n][0]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
  else if (ndim == 2) {
    const double dx = bes->grid.dx[0];
    const double dy = bes->grid.dx[1];
    double pi[4][GKYL_MAX_SPECIES][6] = {};
    double q[4][GKYL_MAX_SPECIES][3] = {};
    for (int n=0; n < nfluids; ++n) {
      for (int j = LL_2D; j <= UU_2D; ++j)
        for (int k = 0; k < 6; ++k)
          pi[j][n][k] = brag_d[j][n][PIXX + k];

      div_pi[n][0] = calc_sym_gradx_2D(dx, pi[LL_2D][n][0], pi[LU_2D][n][0], pi[UL_2D][n][0], pi[UU_2D][n][0])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1]);
      
      div_pi[n][1] = calc_sym_gradx_2D(dx, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][3], pi[LU_2D][n][3], pi[UL_2D][n][3], pi[UU_2D][n][3]);
      
      div_pi[n][2] = calc_sym_gradx_2D(dx, pi[LL_2D][n][2], pi[LU_2D][n][2], pi[UL_2D][n][2], pi[UU_2D][n][2])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][4], pi[LU_2D][n][4], pi[UL_2D][n][4], pi[UU_2D][n][4]);

      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];

      // If energy variable exists, increment heat flux and viscous heating
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = LL_2D; j <= UU_2D; ++j)
          for (int k = 0; k < 3; ++k)
            q[j][n][k] = brag_d[j][n][QX + k];

        div_q[n] = calc_sym_gradx_2D(dx, q[LL_2D][n][0], q[LU_2D][n][0], q[UL_2D][n][0], q[UU_2D][n][0])
                 + calc_sym_grady_2D(dy, q[LL_2D][n][1], q[LU_2D][n][1], q[UL_2D][n][1], q[UU_2D][n][1]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
}

gkyl_moment_braginskii*
gkyl_moment_braginskii_new(struct gkyl_moment_braginskii_inp inp)
{
  gkyl_moment_braginskii *up = gkyl_malloc(sizeof(gkyl_moment_braginskii));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->nfluids = inp.nfluids;
  for (int n=0; n<inp.nfluids; ++n) up->param[n] = inp.param[n];
  up->epsilon0 = inp.epsilon0;
  up->coll_fac = inp.coll_fac;

  return up;
}

static bool
has_mag(const gkyl_moment_braginskii *bes)
{
  bool mag = false;
  for (int n = 0; n < bes->nfluids; ++n)
    if (bes->param[n].type_brag & GKYL_BRAG_MAG)
      mag = true;

  return mag;
}

void
gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes,
  struct gkyl_range brag_vars_range, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate[GKYL_MAX_SPECIES], struct gkyl_array *brag_vars[GKYL_MAX_SPECIES], struct gkyl_array *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  int ndim = update_range.ndim;
  long sz[] = { 2, 4, 8 };

  bool mag = has_mag(bes);

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(&brag_vars_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(&update_range, offsets_centers);
  
  const double* fluid_d[sz[ndim-1]][GKYL_MAX_SPECIES];
  const double* em_tot_d[sz[ndim-1]];
  double *cflrate_d[GKYL_MAX_SPECIES];
  double *brag_vars_d[GKYL_MAX_SPECIES];
  const double* brag_vars_up[sz[ndim-1]][GKYL_MAX_SPECIES];
  double *rhs_d[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, &brag_vars_range);
  while (gkyl_range_iter_next(&iter_vertex)) {
    
    long linc_vertex = gkyl_range_idx(&brag_vars_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_vertex.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]); 
      for (int n=0; n<nfluids; ++n)
        fluid_d[i][n] = gkyl_array_cfetch(fluid[n], linc_center + offsets_vertices[i]);
    }
    for (int n=0; n<nfluids; ++n) {
      brag_vars_d[n] = gkyl_array_fetch(brag_vars[n], linc_vertex);
      cflrate_d[n] = gkyl_array_fetch(cflrate[n], linc_center);
    }

    if (mag)
      mag_brag_calc_vars(bes, fluid_d, em_tot_d, cflrate_d, brag_vars_d);
    else
      unmag_brag_calc_vars(bes, fluid_d, cflrate_d, brag_vars_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, &update_range);
  while (gkyl_range_iter_next(&iter_center)) {
    
    long linc_vertex = gkyl_range_idx(&brag_vars_range, iter_center.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_center.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i)
      for (int n=0; n<nfluids; ++n)
        brag_vars_up[i][n] = gkyl_array_fetch(brag_vars[n], linc_vertex + offsets_centers[i]);
    
    for (int n=0; n<nfluids; ++n)
      rhs_d[n] = gkyl_array_fetch(rhs[n], linc_center);

    brag_calc_update(bes, brag_vars_up, rhs_d);
  }
}

void
gkyl_moment_braginskii_release(gkyl_moment_braginskii* up)
{
  free(up);
}
