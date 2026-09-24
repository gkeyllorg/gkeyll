#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_moment_em_coupling_priv.h>
#include <gkyl_mat.h>

void
explicit_nT_source_update_euler(const double mass, const double dt, double* fluid_old, double* fluid_new, const double* nT_sources)
{
  double rho_old = fluid_old[0];
  double n_old = rho_old / mass;

  double u = fluid_old[1] / rho_old;
  double v = fluid_old[2] / rho_old;
  double w = fluid_old[3] / rho_old;
  double v_sq = (u * u) + (v * v) + (w * w);
  double TT_old = (fluid_old[4] - (0.5 * rho_old * v_sq)) / n_old;

  double n_new = n_old + (dt * nT_sources[0]);
  double TT_new = TT_old + (dt * nT_sources[1]);

  double rho_new = n_new * mass;
  fluid_new[0] = rho_new;
  fluid_new[1] = rho_new * u;
  fluid_new[2] = rho_new * v;
  fluid_new[3] = rho_new * w;
  fluid_new[4] = (n_new * TT_new) + (0.5 * rho_new * v_sq);
}

void
explicit_nT_source_update(const gkyl_moment_em_coupling* mom_em, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* nT_sources_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double *nT_sources = nT_sources_s[i];

    double mass = mom_em->param[i].mass;

    explicit_nT_source_update_euler(mass, dt, f, f, nT_sources);
  }
}

void
explicit_frictional_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double Z, const double T_elc, const double Lambda_ee,
  double t_curr, const double dt, double* f_elc_old, double* f_ion_old, double* f_elc_new, double* f_ion_new)
{
  int nfluids = mom_em->nfluids;
  double pi = M_PI;

  if (nfluids == 2) {
    double mass_elc = mom_em->param[0].mass;
    double epsilon0 = mom_em->epsilon0;
    
    double rho_elc = f_elc_old[0];
    double rho_ion = f_ion_old[0];

    double u_elc = f_elc_old[1], v_elc = f_elc_old[2], w_elc = f_elc_old[3];
    double u_ion = f_ion_old[1], v_ion = f_ion_old[2], w_ion = f_ion_old[3];

    double n_elc = rho_elc / mass_elc;

    double tau_ei = (1.0 / Z) * ((3.0 * sqrt(mass_elc) * ((4.0 * pi * epsilon0) * (4.0 * pi * epsilon0)) * pow(T_elc, 3.0 / 2.0)) /
      (4.0 * sqrt(2.0 * pi) * n_elc * exp(4.0) * log(Lambda_ee)));
    double alpha_par = 1.0 - (pow(Z, 2.0 / 3.0) / ((1.46 * pow(Z, 2.0 / 3.0)) - (0.33 * pow (Z, 1.0 / 3.0)) + 0.888));

    double mom_src_x = -(rho_elc / tau_ei) * (alpha_par * (u_elc - u_ion));
    double mom_src_y = -(rho_elc / tau_ei) * (alpha_par * (v_elc - v_ion));
    double mom_src_z = -(rho_elc / tau_ei) * (alpha_par * (w_elc - w_ion));

    double E_src = (mom_src_x * u_ion) + (mom_src_y * v_ion) + (mom_src_z * w_ion);

    f_elc_new[1] = f_elc_old[1] + (dt * mom_src_x);
    f_elc_new[2] = f_elc_old[2] + (dt * mom_src_y);
    f_elc_new[3] = f_elc_old[3] + (dt * mom_src_z);

    f_ion_new[1] = f_ion_old[1] - (dt * mom_src_x);
    f_ion_new[2] = f_ion_old[2] - (dt * mom_src_y);
    f_ion_new[3] = f_ion_old[3] - (dt * mom_src_z);

    if (mom_em->param[0].type == GKYL_EQN_EULER) {
      f_elc_new[4] = f_elc_old[4] + (dt * E_src);
    }
    if (mom_em->param[1].type == GKYL_EQN_EULER) {
      f_ion_new[4] = f_ion_old[4] - (dt * E_src);
    }

    f_elc_new[0] = f_elc_old[0];
    f_ion_new[0] = f_ion_old[0];
  }
}

void
explicit_frictional_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  if (nfluids == 2) {
    double *f_elc = fluid_s[0];
    double *f_ion = fluid_s[1];

    double Z = mom_em->friction_Z;
    double T_elc = mom_em->friction_T_elc;
    double Lambda_ee = mom_em->friction_Lambda_ee;

    int elc_num_equations = 0;
    int ion_num_equations = 0;

    if (mom_em->param[0].type == GKYL_EQN_EULER) {
      elc_num_equations = 5;
    }
    else if (mom_em->param[0].type == GKYL_EQN_ISO_EULER) {
      elc_num_equations = 4;
    }

    // Right now, we allocate these arrays on the stack with size 5 (current maximum supported number of equations).
    double f_elc_new[5];
    double f_elc_stage1[5];
    double f_elc_stage2[5];
    double f_elc_old[5];

    if (mom_em->param[1].type == GKYL_EQN_EULER) {
      ion_num_equations = 5;
    }
    else if (mom_em->param[1].type == GKYL_EQN_ISO_EULER) {
      ion_num_equations = 4;
    }

    // Right now, we allocate these arrays on the stack with size 5 (current maximum supported number of equations).
    double f_ion_new[5];
    double f_ion_stage1[5];
    double f_ion_stage2[5];
    double f_ion_old[5];

    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_old[i] = f_elc[i];
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_old[i] = f_ion[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr, dt, f_elc_old, f_ion_old, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_stage1[i] = f_elc_new[i];
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_stage1[i] = f_ion_new[i];
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + dt, dt, f_elc_stage1, f_ion_stage1, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc_stage2[i] = (0.75 * f_elc_old[i]) + (0.25 * f_elc_new[i]);
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion_stage2[i] = (0.75 * f_ion_old[i]) + (0.25 * f_ion_new[i]);
    }

    explicit_frictional_source_update_euler(mom_em, Z, T_elc, Lambda_ee, t_curr + (0.5 * dt), dt, f_elc_stage2, f_ion_stage2, f_elc_new, f_ion_new);
    for (int i = 0; i < elc_num_equations; i++) {
      f_elc[i] = ((1.0 / 3.0) * f_elc_old[i]) + ((2.0 / 3.0) * f_elc_new[i]);
    }
    for (int i = 0; i < ion_num_equations; i++) {
      f_ion[i] = ((1.0 / 3.0) * f_ion_old[i]) + ((2.0 / 3.0) * f_ion_new[i]);
    }
  }
}

void
explicit_volume_source_5m_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double U0, const double R0,
  double t_curr, const double dt, double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  double p = (gas_gamma - 1.0) * (fluid_old[4] - (0.5 * rho * ((vx * vx) + (vy * vy) + (vz * vz))));

  double a = 1.0 + ((U0 * t_curr) / R0);

  for (int i = 0; i < 5; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] -= dt * ((2.0 * U0) / (a * R0)) * rho;

  fluid_new[2] -= dt * (U0 / (a * R0)) * rho * vy;
  fluid_new[3] -= dt * (U0 / (a * R0)) * rho * vz;

  double p_new = p;
  p_new -= dt * (gas_gamma * ((2.0 * U0) / (a * R0)) * p);

  fluid_new[4] = (p_new / (gas_gamma - 1.0)) + (0.5 * rho * (vx * vx) + (vy * vy) + (vz * vz));
}

void
explicit_volume_source_10m_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  
  double Pxx = fluid_old[4], Pxy = fluid_old[5], Pxz = fluid_old[6];
  double Pyy = fluid_old[7], Pyz = fluid_old[8], Pzz = fluid_old[9];

  double a = 1.0 + ((U0 * t_curr) / R0);

  for (int i = 0; i < 10; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] -= dt * ((2.0 * U0) / (a * R0)) * rho;

  fluid_new[2] -= dt * (U0 / (a * R0)) * rho * vy;
  fluid_new[3] -= dt * (U0 / (a * R0)) * rho * vz;

  fluid_new[4] -= dt * ((2.0 * U0) / (a * R0)) * Pxx;
  fluid_new[5] -= dt * ((((2.0 * U0) / (a * R0)) * Pxy) + ((U0 / (a * R0)) * Pxy));
  fluid_new[6] -= dt * ((((2.0 * U0) / (a * R0)) * Pxz) + ((U0 / (a * R0)) * Pxz));
  fluid_new[7] -= dt * ((((2.0 * U0) / (a * R0)) * Pyy) + (((2.0 * U0) / (a * R0)) * Pyy));
  fluid_new[8] -= dt * ((((2.0 * U0) / (a * R0)) * Pyz) + (((2.0 * U0) / (a * R0)) * Pyz));
  fluid_new[9] -= dt * ((((2.0 * U0) / (a * R0)) * Pzz) + (((2.0 * U0) / (a * R0)) * Pzz));
}

void
explicit_volume_source_maxwell_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr,
  const double dt, double* em_old, double* em_new, const double* ext_em)
{
  double a = 1.0 + ((U0 * t_curr) / R0);

  double Ex = em_old[0] + ext_em[0];
  double Ey = em_old[1] + ext_em[1];
  double Ez = em_old[2] + ext_em[2];

  double Bx = em_old[3] + ext_em[3];
  double By = em_old[4] + ext_em[4];
  double Bz = em_old[5] + ext_em[5];

  for (int i = 0; i < 8; i++) {
    em_new[i] = em_old[i];
  }

  em_new[0] -= dt * (U0 / (a * R0)) * 2.0 * Ex;
  em_new[1] -= dt * (U0 / (a * R0)) * Ey;
  em_new[2] -= dt * (U0 / (a * R0)) * Ez;

  em_new[3] -= dt * (U0 / (a * R0)) * 2.0 * Bx;
  em_new[4] -= dt * (U0 / (a * R0)) * By;
  em_new[5] -= dt * (U0 / (a * R0)) * Bz;
}

void
explicit_volume_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  double* em, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->volume_gas_gamma;
  double U0 = mom_em->volume_U0;
  double R0 = mom_em->volume_R0;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    if (mom_em->param[i].type == GKYL_EQN_EULER) {
      double f_new[5], f_stage1[5], f_stage2[5], f_old[5];

      for (int j = 0; j < 5; j++) {
        f_old[j] = f[j];
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr, dt, f_old, f_new);
      for (int j = 0; j < 5; j++) {
        f_stage1[j] = f_new[j];
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr + dt, dt, f_stage1, f_new);
      for (int j = 0; j < 5; j++) {
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
      }

      explicit_volume_source_5m_update_euler(mom_em, gas_gamma, U0, R0, t_curr + (0.5 * dt), dt, f_stage2, f_new);
      for (int j = 0; j < 5; j++) {
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
      }
    }
    else if (mom_em->param[i].type == GKYL_EQN_TEN_MOMENT) {
      double f_new[10], f_stage1[10], f_stage2[10], f_old[10];

      for (int j = 0; j < 10; j++) {
        f_old[j] = f[j];
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr, dt, f_old, f_new);
      for (int j = 0; j < 10; j++) {
        f_stage1[j] = f_new[j];
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr + dt, dt, f_stage1, f_new);
      for (int j = 0; j < 10; j++) {
        f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
      }

      explicit_volume_source_10m_update_euler(mom_em, U0, R0, t_curr + (0.5 * dt), dt, f_stage2, f_new);
      for (int j = 0; j < 10; j++) {
        f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
      }
    }
  }

  double em_new[8], em_stage1[8], em_stage2[8], em_old[8];

  for (int i = 0; i < 8; i++) {
    em_old[i] = em[i];
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr, dt, em_old, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em_stage1[i] = em_new[i];
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr + dt, dt, em_stage1, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em_stage2[i] = (0.75 * em_old[i]) + (0.25 * em_new[i]);
  }

  explicit_volume_source_maxwell_update_euler(mom_em, U0, R0, t_curr + (0.5 * dt), dt, em_stage2, em_new, ext_em);
  for (int i = 0; i < 8; i++) {
    em[i] = ((1.0 / 3.0) * em_old[i]) + ((2.0 / 3.0) * em_new[i]);
  }
}

void
explicit_reactive_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double specific_heat_capacity,
  const double energy_of_formation, const double ignition_temperature, const double reaction_rate, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double rho = fluid_old[0];
  double vx = fluid_old[1] / rho;
  double vy = fluid_old[2] / rho;
  double vz = fluid_old[3] / rho;
  double reaction_progress = fluid_old[5] / rho;

  double specific_internal_energy = (fluid_old[4] / rho) - (0.5 * ((vx * vx) + (vy * vy) + (vz * vz))) -
    (energy_of_formation * (reaction_progress - 1.0));
  double temperature = specific_internal_energy / specific_heat_capacity;

  for (int i = 0; i < 6; i++) {
    fluid_new[i] = fluid_old[i];
  }

  if (temperature > ignition_temperature) {
    fluid_new[5] -= dt * (rho * reaction_progress * reaction_rate);
  }
}

void
explicit_reactive_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->reactivity_gas_gamma;
  double specific_heat_capacity = mom_em->reactivity_specific_heat_capacity;
  double energy_of_formation = mom_em->reactivity_energy_of_formation;
  double ignition_temperature = mom_em->reactivity_ignition_temperature;
  double reaction_rate = mom_em->reactivity_reaction_rate;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[6], f_stage1[6], f_stage2[6], f_old[6];

    for (int j = 0; j < 6; j++) {
      f_old[j] = f[j];
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr, dt, f_old, f_new);
    for (int j = 0; j < 6; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 6; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_reactive_source_update_euler(mom_em, gas_gamma, specific_heat_capacity, energy_of_formation, ignition_temperature, reaction_rate,
      t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 6; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_medium_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double kappa, double t_curr,
  const double dt, double* fluid_old, double* fluid_new)
{
  double exp_2a = fluid_old[0];

  double a_dt = fluid_old[1], a_dx = fluid_old[2];
  double b_dt = fluid_old[3], b_dx = fluid_old[4];
  double c_dt = fluid_old[5], c_dx = fluid_old[6];
  
  double a_dt_dx = fluid_old[7], a_dx_dx = fluid_old[8];
  double b_dt_dx = fluid_old[9], b_dx_dx = fluid_old[10];
  double c_dt_dx = fluid_old[11], c_dx_dx = fluid_old[12];

  double Etot = fluid_old[13];
  double mom = fluid_old[14];

  double rho = (1.0 / (gas_gamma - 1.0)) * ((-0.5 * (2.0 - gas_gamma) * Etot) + sqrt((0.25 * (2.0 - gas_gamma) * (2.0 - gas_gamma) * Etot * Etot) +
    ((gas_gamma - 1.0) * ((Etot * Etot) - (mom * mom)))));

  double vel = 0.0;
  if (fabs(mom) > pow(10.0, -8.0)) {
    vel = ((gas_gamma * rho) / (2.0 * mom)) * (sqrt(1.0 + ((4 * mom * mom) / ((gas_gamma * gas_gamma) * (rho * rho)))) - 1.0);
  }

  double p = (gas_gamma - 1.0) * rho;

  double W = 1.0 / sqrt(1.0 - (vel * vel));
  if (vel * vel > 1.0 - pow(1.0, -8.0)) {
    W = 1.0 / sqrt(pow(1.0, -8.0));
  }

  for (int i = 0; i < 15; i++) {
    fluid_new[i] = fluid_old[i];
  }

  fluid_new[0] += dt * (2.0 * a_dt * exp_2a);

  fluid_new[1] += dt * (a_dx_dx + (b_dt * b_dt) - (b_dx * b_dx) - (c_dt * c_dt) + (c_dx * c_dx) - (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))));
  fluid_new[2] += dt * (a_dt_dx);
  fluid_new[3] += dt * (b_dx_dx - (2.0 * (b_dt * b_dt)) + (2.0 * (b_dx * b_dx)) + (0.5 * kappa * exp_2a * (Etot - ((mom * vel) + p))));
  fluid_new[4] += dt * (b_dt_dx);
  fluid_new[5] += dt * (c_dx_dx - (2.0 * ((b_dt * c_dt) - (b_dx * c_dx))));
  fluid_new[6] += dt * (c_dt_dx);

  fluid_new[7] += dt * ((2.0 * (b_dt * b_dt_dx)) - (2.0 * (b_dx * b_dx_dx)) - (2.0 * (c_dt * c_dt_dx)) + (2.0 * (c_dx * c_dx_dx)));
  fluid_new[8] += 0.0;
  fluid_new[9] += dt * (-(4.0 * (b_dt * b_dt_dx)) + (4.0 * (b_dx * b_dx_dx)));
  fluid_new[10] += 0.0;
  fluid_new[11] += dt * (-2.0 * ((b_dt * c_dt_dx) - (b_dx * c_dx_dx) + (b_dt_dx * c_dt) - (b_dx_dx * c_dx)));
  fluid_new[12] += 0.0;

  fluid_new[13] += dt * ((-Etot * (a_dt + (2.0 * b_dt))) - (2.0 * mom * (a_dx + b_dx)) - (((mom * vel) + p) * a_dt) - (2.0 * p * b_dt));
  fluid_new[14] += dt * ((-Etot * a_dx) - (2.0 * mom * (a_dt + b_dt)) - (((mom * vel) + p) * (a_dx + (2.0 * b_dx))) + (2.0 * p * b_dx));
}

void
explicit_medium_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->medium_gas_gamma;
  double kappa = mom_em->medium_kappa;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[15], f_stage1[15], f_stage2[15], f_old[15];

    for (int j = 0; j < 15; j++) {
      f_old[j] = f[j];
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 15; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 15; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_medium_source_update_euler(mom_em, gas_gamma, kappa, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 15; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_gr_ultra_rel_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[4];
  double shift_x = fluid_old[5];
  double shift_y = fluid_old[6];
  double shift_z = fluid_old[7];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[8]; spatial_metric[0][1] = fluid_old[9]; spatial_metric[0][2] = fluid_old[10];
  spatial_metric[1][0] = fluid_old[11]; spatial_metric[1][1] = fluid_old[12]; spatial_metric[1][2] = fluid_old[13];
  spatial_metric[2][0] = fluid_old[14]; spatial_metric[2][1] = fluid_old[15]; spatial_metric[2][2] = fluid_old[16];

  double inv_spatial_metric[3][3];
  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[17]; extrinsic_curvature[0][1] = fluid_old[18]; extrinsic_curvature[0][2] = fluid_old[19];
  extrinsic_curvature[1][0] = fluid_old[20]; extrinsic_curvature[1][1] = fluid_old[21]; extrinsic_curvature[1][2] = fluid_old[22];
  extrinsic_curvature[2][0] = fluid_old[23]; extrinsic_curvature[2][1] = fluid_old[24]; extrinsic_curvature[2][2] = fluid_old[25];

  bool in_excision_region = false;
  if (fluid_old[26] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double Etot = fluid_old[0] / sqrt(spatial_det);
    double momx = fluid_old[1] / sqrt(spatial_det);
    double momy = fluid_old[2] / sqrt(spatial_det);
    double momz = fluid_old[3] / sqrt(spatial_det);

    double mom[3];
    double mom_sq = 0.0;
    mom[0] = momx; mom[1] = momy; mom[2] = momz;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        mom_sq += inv_spatial_metric[i][j] * mom[i] * mom[j];
      }
    }

    double beta = 0.25 * (2.0 - gas_gamma);
    double p = -(2.0 * beta * Etot) + sqrt((4.0 * (beta * beta) * (Etot * Etot)) + ((gas_gamma - 1.0) * ((Etot * Etot) - mom_sq)));
    if (p < pow(10.0, -8.0)) {
      p = pow(10.0, -8.0);
    }
    double rho = p / (gas_gamma - 1.0);

    double cov_vx = momx / (Etot + p);
    double cov_vy = momy / (Etot + p);
    double cov_vz = momz / (Etot + p);

    if (Etot + p < pow(10.0, -8.0)) {
      cov_vx = momx / pow(10.0, -8.0);
      cov_vy = momy / pow(10.0, -8.0);
      cov_vz = momz / pow(10.0, -8.0);
    }

    if (cov_vx > 1.0 - pow(10.0, -8.0)) {
      cov_vx = 1.0 - pow(10.0, -8.0);
    }
    if (cov_vy > 1.0 - pow(10.0, -8.0)) {
      cov_vy = 1.0 - pow(10.0, -8.0);
    }
    if (cov_vz > 1.0 - pow(10.0, -8.0)) {
      cov_vz = 1.0 - pow(10.0, -8.0);
    }

    double cov_vel[3];
    cov_vel[0] = cov_vx; cov_vel[1] = cov_vy; cov_vel[2] = cov_vz;

    double vel[3];
    for (int i = 0; i < 3; i++) {
      vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        vel[i] += inv_spatial_metric[i][j] * cov_vel[j];
      }
    }

    double vx = vel[0];
    double vy = vel[1];
    double vz = vel[2];

    double v_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        v_sq += spatial_metric[i][j] * vel[i] * vel[j];
      }
    }

    double W = 1.0 / sqrt(1.0 - v_sq);
    if (v_sq > 1.0 - pow(10.0, -8.0)) {
      W = 1.0 / sqrt(pow(10.0, -8.0));
    }

    double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift_x * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift_y * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift_z * (W / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double stress_energy[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
      }
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[27];
    lapse_der[1] = fluid_old[28];
    lapse_der[2] = fluid_old[29];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[30]; shift_der[0][1] = fluid_old[31]; shift_der[0][2] = fluid_old[32];
    shift_der[1][0] = fluid_old[33]; shift_der[1][1] = fluid_old[34]; shift_der[1][2] = fluid_old[35];
    shift_der[2][0] = fluid_old[36]; shift_der[2][1] = fluid_old[37]; shift_der[2][2] = fluid_old[38];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[39]; spatial_metric_der[0][0][1] = fluid_old[40]; spatial_metric_der[0][0][2] = fluid_old[41];
    spatial_metric_der[0][1][0] = fluid_old[42]; spatial_metric_der[0][1][1] = fluid_old[43]; spatial_metric_der[0][1][2] = fluid_old[44];
    spatial_metric_der[0][2][0] = fluid_old[45]; spatial_metric_der[0][2][1] = fluid_old[46]; spatial_metric_der[0][2][2] = fluid_old[47];

    spatial_metric_der[1][0][0] = fluid_old[48]; spatial_metric_der[1][0][1] = fluid_old[49]; spatial_metric_der[1][0][2] = fluid_old[50];
    spatial_metric_der[1][1][0] = fluid_old[51]; spatial_metric_der[1][1][1] = fluid_old[52]; spatial_metric_der[1][1][2] = fluid_old[53];
    spatial_metric_der[1][2][0] = fluid_old[54]; spatial_metric_der[1][2][1] = fluid_old[55]; spatial_metric_der[1][2][2] = fluid_old[56];

    spatial_metric_der[2][0][0] = fluid_old[57]; spatial_metric_der[2][0][1] = fluid_old[58]; spatial_metric_der[2][0][2] = fluid_old[59];
    spatial_metric_der[2][1][0] = fluid_old[60]; spatial_metric_der[2][1][1] = fluid_old[61]; spatial_metric_der[2][1][2] = fluid_old[62];
    spatial_metric_der[2][2][0] = fluid_old[63]; spatial_metric_der[2][2][1] = fluid_old[64]; spatial_metric_der[2][2][2] = fluid_old[65];

    mom[0] = (rho + p) * (W * W) * vx;
    mom[1] = (rho + p) * (W * W) * vy;
    mom[2] = (rho + p) * (W * W) * vz;

    for (int i = 0; i < 70; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[0] += dt * (stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[0] += dt * (2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[0] += dt * (stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[0] -= dt * (stress_energy[0][0] * shift[i] * lapse_der[i]);
      fluid_new[0] -= dt * (stress_energy[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[1 + j] -= dt * (stress_energy[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[1 + j] += dt * (0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[1 + j] += dt * (0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[1 + j] += dt * ((mom[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[1 + j] += dt * (stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 70; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_ultra_rel_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->gr_ultra_rel_gas_gamma;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[70], f_stage1[70], f_stage2[70], f_old[70];

    for (int j = 0; j < 70; j++) {
      f_old[j] = f[j];
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 70; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 70; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_gr_ultra_rel_source_update_euler(mom_em, gas_gamma, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 70; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_gr_euler_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[5];
  double shift_x = fluid_old[6];
  double shift_y = fluid_old[7];
  double shift_z = fluid_old[8];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[9]; spatial_metric[0][1] = fluid_old[10]; spatial_metric[0][2] = fluid_old[11];
  spatial_metric[1][0] = fluid_old[12]; spatial_metric[1][1] = fluid_old[13]; spatial_metric[1][2] = fluid_old[14];
  spatial_metric[2][0] = fluid_old[15]; spatial_metric[2][1] = fluid_old[16]; spatial_metric[2][2] = fluid_old[17];

  double inv_spatial_metric[3][3];
  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[18]; extrinsic_curvature[0][1] = fluid_old[19]; extrinsic_curvature[0][2] = fluid_old[20];
  extrinsic_curvature[1][0] = fluid_old[21]; extrinsic_curvature[1][1] = fluid_old[22]; extrinsic_curvature[1][2] = fluid_old[23];
  extrinsic_curvature[2][0] = fluid_old[24]; extrinsic_curvature[2][1] = fluid_old[25]; extrinsic_curvature[2][2] = fluid_old[26];

  bool in_excision_region = false;
  if (fluid_old[27] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double D = fluid_old[0] / sqrt(spatial_det);
    double momx = fluid_old[1] / sqrt(spatial_det);
    double momy = fluid_old[2] / sqrt(spatial_det);
    double momz = fluid_old[3] / sqrt(spatial_det);
    double Etot = fluid_old[4] / sqrt(spatial_det);

    double C = D / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
    double C0 = (D + Etot) / sqrt(((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)));
    if (((Etot + D) * (Etot + D)) - ((momx * momx) + (momy * momy) + (momz * momz)) < pow(10.0, -8.0)) {
      C = D / sqrt(pow(10.0, -8.0));
      C0 = (D + Etot) / sqrt(pow(10.0, -8.0));
    }

    double alpha0 = -1.0 / (gas_gamma * gas_gamma);
    double alpha1 = -2.0 * C * ((gas_gamma - 1.0) / (gas_gamma * gas_gamma));
    double alpha2 = ((gas_gamma - 2.0) / gas_gamma) * ((C0 * C0) - 1.0) + 1.0 - (C * C) * ((gas_gamma - 1.0) / gas_gamma) * ((gas_gamma - 1.0) / gas_gamma);
    double alpha4 = (C0 * C0) - 1.0;
    double eta = 2.0 * C *((gas_gamma - 1.0) / gas_gamma);

    double guess = 1.0;
    int iter = 0;

    while (iter < 100) {
      double poly = (alpha4 * (guess * guess * guess) * (guess - eta)) + (alpha2 * (guess * guess)) + (alpha1 * guess) + alpha0;
      double poly_der = alpha1 + (2.0 * alpha2 * guess) + (4.0 * alpha4 * (guess * guess * guess)) - (3.0 * eta * alpha4 * (guess * guess));

      double guess_new = guess - (poly / poly_der);

      if (fabs(guess - guess_new) < pow(10.0, -8.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        guess = guess_new;
      }
    }

    double W = 0.5 * C0 * guess * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma - 1.0) / gas_gamma) * ((1.0 - (C * guess)) / ((C0 * C0) * (guess * guess))))));
    double h = 1.0 / (C * guess);

    double rho = D / W;
    double vx = momx / (rho * h * (W * W));
    double vy = momy / (rho * h * (W * W));
    double vz = momz / (rho * h * (W * W));
    double p = (rho * h * (W * W)) - D - Etot;

    if (rho < pow(10.0, -8.0)) {
      rho = pow(10.0, -8.0);
    }
    if (p < pow(10.0, -8.0)) {
      p = pow(10.0, -8.0);
    }

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    spacetime_vel[1] = (W * vx) - (shift_x * (W / lapse));
    spacetime_vel[2] = (W * vy) - (shift_y * (W / lapse));
    spacetime_vel[3] = (W * vz) - (shift_z * (W / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double stress_energy[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy[i][j] = (rho * h * spacetime_vel[i] * spacetime_vel[j]) + (p * inv_spacetime_metric[i][j]);
      }
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[28];
    lapse_der[1] = fluid_old[29];
    lapse_der[2] = fluid_old[30];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[31]; shift_der[0][1] = fluid_old[32]; shift_der[0][2] = fluid_old[33];
    shift_der[1][0] = fluid_old[34]; shift_der[1][1] = fluid_old[35]; shift_der[1][2] = fluid_old[36];
    shift_der[2][0] = fluid_old[37]; shift_der[2][1] = fluid_old[38]; shift_der[2][2] = fluid_old[39];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[40]; spatial_metric_der[0][0][1] = fluid_old[41]; spatial_metric_der[0][0][2] = fluid_old[42];
    spatial_metric_der[0][1][0] = fluid_old[43]; spatial_metric_der[0][1][1] = fluid_old[44]; spatial_metric_der[0][1][2] = fluid_old[45];
    spatial_metric_der[0][2][0] = fluid_old[46]; spatial_metric_der[0][2][1] = fluid_old[47]; spatial_metric_der[0][2][2] = fluid_old[48];

    spatial_metric_der[1][0][0] = fluid_old[49]; spatial_metric_der[1][0][1] = fluid_old[50]; spatial_metric_der[1][0][2] = fluid_old[51];
    spatial_metric_der[1][1][0] = fluid_old[52]; spatial_metric_der[1][1][1] = fluid_old[53]; spatial_metric_der[1][1][2] = fluid_old[54];
    spatial_metric_der[1][2][0] = fluid_old[55]; spatial_metric_der[1][2][1] = fluid_old[56]; spatial_metric_der[1][2][2] = fluid_old[57];

    spatial_metric_der[2][0][0] = fluid_old[58]; spatial_metric_der[2][0][1] = fluid_old[59]; spatial_metric_der[2][0][2] = fluid_old[60];
    spatial_metric_der[2][1][0] = fluid_old[61]; spatial_metric_der[2][1][1] = fluid_old[62]; spatial_metric_der[2][1][2] = fluid_old[63];
    spatial_metric_der[2][2][0] = fluid_old[64]; spatial_metric_der[2][2][1] = fluid_old[65]; spatial_metric_der[2][2][2] = fluid_old[66];

    double mom[3];
    mom[0] = (rho + p) * (W * W) * vx;
    mom[1] = (rho + p) * (W * W) * vy;
    mom[2] = (rho + p) * (W * W) * vz;

    for (int i = 0; i < 71; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[4] += dt * (stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[4] -= dt * (stress_energy[0][0] * shift[i] * lapse_der[i]);
      fluid_new[4] -= dt * (stress_energy[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[1 + j] -= dt * (stress_energy[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[1 + j] += dt * (0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[1 + j] += dt * (0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[1 + j] += dt * ((mom[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[1 + j] += dt * (stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 71; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_euler_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->gr_euler_gas_gamma;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[71], f_stage1[71], f_stage2[71], f_old[71];

    for (int j = 0; j < 71; j++) {
      f_old[j] = f[j];
    }

    explicit_gr_euler_source_update_euler(mom_em, gas_gamma, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 71; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_gr_euler_source_update_euler(mom_em, gas_gamma, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 71; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_gr_euler_source_update_euler(mom_em, gas_gamma, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 71; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_gr_twofluid_source_update_elc_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, const double mass_elc,
  const double charge_elc, double t_curr, const double dt, double* fluid_old, double* fluid_new)
{
  for (int i = 0; i < 84; i++) {
    fluid_new[i] = fluid_old[i];
  }

  double Dx = fluid_old[10], Dy = fluid_old[11], Dz = fluid_old[12];
  double Bx = fluid_old[13], By = fluid_old[14], Bz = fluid_old[15];

  double lapse = fluid_old[18];
  double shift_x = fluid_old[19];
  double shift_y = fluid_old[20];
  double shift_z = fluid_old[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[22]; spatial_metric[0][1] = fluid_old[23]; spatial_metric[0][2] = fluid_old[24];
  spatial_metric[1][0] = fluid_old[25]; spatial_metric[1][1] = fluid_old[26]; spatial_metric[1][2] = fluid_old[27];
  spatial_metric[2][0] = fluid_old[28]; spatial_metric[2][1] = fluid_old[29]; spatial_metric[2][2] = fluid_old[30];

  bool in_excision_region = false;
  if (fluid_old[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D_elc = fluid_old[0] / sqrt(spatial_det);
    double momx_elc = fluid_old[1] / sqrt(spatial_det);
    double momy_elc = fluid_old[2] / sqrt(spatial_det);
    double momz_elc = fluid_old[3] / sqrt(spatial_det);
    double Etot_elc = fluid_old[4] / sqrt(spatial_det);

    double C_elc = D_elc / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    double C0_elc = (D_elc + Etot_elc) / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) +
      (momz_elc * momz_elc)));
    if (((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)) < pow(10.0, -8.0)) {
      C_elc = D_elc / sqrt(pow(10.0, -8.0));
      C0_elc = (D_elc + Etot_elc) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_elc = -1.0 / (gas_gamma_elc * gas_gamma_elc);
    double alpha1_elc = -2.0 * C_elc * ((gas_gamma_elc - 1.0) / (gas_gamma_elc * gas_gamma_elc));
    double alpha2_elc = ((gas_gamma_elc - 2.0) / gas_gamma_elc) * ((C0_elc * C0_elc) - 1.0) + 1.0 - (C_elc * C_elc) *
      ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((gas_gamma_elc - 1.0) / gas_gamma_elc);
    double alpha4_elc = (C0_elc * C0_elc) - 1.0;
    double eta_elc = 2.0 * C_elc *((gas_gamma_elc - 1.0) / gas_gamma_elc);

    double guess_elc = 1.0;
    int iter_elc = 0;

    while (iter_elc < 100) {
      double poly_elc = (alpha4_elc * (guess_elc * guess_elc * guess_elc) * (guess_elc - eta_elc)) + (alpha2_elc * (guess_elc * guess_elc)) +
        (alpha1_elc * guess_elc) + alpha0_elc;
      double poly_der_elc = alpha1_elc + (2.0 * alpha2_elc * guess_elc) + (4.0 * alpha4_elc * (guess_elc * guess_elc * guess_elc)) -
        (3.0 * eta_elc * alpha4_elc * (guess_elc * guess_elc));

      double guess_new_elc = guess_elc - (poly_elc / poly_der_elc);

      if (fabs(guess_elc - guess_new_elc) < pow(10.0, -8.0)) {
        iter_elc = 100;
      }
      else {
        iter_elc += 1;
        guess_elc = guess_new_elc;
      }
    }

    double W_elc = 0.5 * C0_elc * guess_elc * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((1.0 - (C_elc * guess_elc)) /
      ((C0_elc * C0_elc) * (guess_elc * guess_elc))))));
    double h_elc = 1.0 / (C_elc * guess_elc);

    double rho_elc = D_elc / W_elc; 
    double vx_elc = momx_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vy_elc = momy_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vz_elc = momz_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double p_elc = (rho_elc * h_elc * (W_elc * W_elc)) - D_elc - Etot_elc;

    if (rho_elc < pow(10.0, -8.0)) {
      rho_elc = pow(10.0, -8.0);
    }
    if (p_elc < pow(10.0, -8.0)) {
      p_elc = pow(10.0, -8.0);
    }

    double D[3], B[3];
    D[0] = Dx; D[1] = Dy; D[2] = Dz;
    B[0] = Bx; B[1] = By; B[2] = Bz;

    double covD[3], covB[3];
    for (int i = 0; i < 3; i++) {
      covD[i] = 0.0;
      covB[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        covD[i] += spatial_metric[i][j] * D[j];
        covB[i] += spatial_metric[i][j] * B[j];
      }
    }

    fluid_new[1] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * covD[0];
    fluid_new[1] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * ((vy_elc * covB[2]) - (vz_elc * covB[1]));
    fluid_new[2] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * covD[1];
    fluid_new[2] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * ((vz_elc * covB[0]) - (vx_elc * covB[2]));
    fluid_new[3] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * covD[2];
    fluid_new[3] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * ((vx_elc * covB[1]) - (vy_elc * covB[0]));

    fluid_new[4] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vx_elc * covD[0];
    fluid_new[4] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vy_elc * covD[1];
    fluid_new[4] += (lapse * sqrt(spatial_det)) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vz_elc * covD[2];
  }
}

void
explicit_gr_twofluid_source_update_ion_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_ion, const double mass_ion,
  const double charge_ion, double t_curr, const double dt, double* fluid_old, double* fluid_new)
{
  for (int i = 0; i < 84; i++) {
    fluid_new[i] = fluid_old[i];
  }

  double Dx = fluid_old[10], Dy = fluid_old[11], Dz = fluid_old[12];
  double Bx = fluid_old[13], By = fluid_old[14], Bz = fluid_old[15];

  double lapse = fluid_old[18];
  double shift_x = fluid_old[19];
  double shift_y = fluid_old[20];
  double shift_z = fluid_old[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[22]; spatial_metric[0][1] = fluid_old[23]; spatial_metric[0][2] = fluid_old[24];
  spatial_metric[1][0] = fluid_old[25]; spatial_metric[1][1] = fluid_old[26]; spatial_metric[1][2] = fluid_old[27];
  spatial_metric[2][0] = fluid_old[28]; spatial_metric[2][1] = fluid_old[29]; spatial_metric[2][2] = fluid_old[30];

  bool in_excision_region = false;
  if (fluid_old[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D_ion = fluid_old[5] / sqrt(spatial_det);
    double momx_ion = fluid_old[6] / sqrt(spatial_det);
    double momy_ion = fluid_old[7] / sqrt(spatial_det);
    double momz_ion = fluid_old[8] / sqrt(spatial_det);
    double Etot_ion = fluid_old[9] / sqrt(spatial_det);

    double C_ion = D_ion / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    double C0_ion = (D_ion + Etot_ion) / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) +
      (momz_ion * momz_ion)));
    if (((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)) < pow(10.0, -8.0)) {
      C_ion = D_ion / sqrt(pow(10.0, -8.0));
      C0_ion = (D_ion + Etot_ion) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_ion = -1.0 / (gas_gamma_ion * gas_gamma_ion);
    double alpha1_ion = -2.0 * C_ion * ((gas_gamma_ion - 1.0) / (gas_gamma_ion * gas_gamma_ion));
    double alpha2_ion = ((gas_gamma_ion - 2.0) / gas_gamma_ion) * ((C0_ion * C0_ion) - 1.0) + 1.0 - (C_ion * C_ion) *
      ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((gas_gamma_ion - 1.0) / gas_gamma_ion);
    double alpha4_ion = (C0_ion * C0_ion) - 1.0;
    double eta_ion = 2.0 * C_ion *((gas_gamma_ion - 1.0) / gas_gamma_ion);

    double guess_ion = 1.0;
    int iter_ion = 0;

    while (iter_ion < 100) {
      double poly_ion = (alpha4_ion * (guess_ion * guess_ion * guess_ion) * (guess_ion - eta_ion)) + (alpha2_ion * (guess_ion * guess_ion)) +
        (alpha1_ion * guess_ion) + alpha0_ion;
      double poly_der_ion = alpha1_ion + (2.0 * alpha2_ion * guess_ion) + (4.0 * alpha4_ion * (guess_ion * guess_ion * guess_ion)) -
        (3.0 * eta_ion * alpha4_ion * (guess_ion * guess_ion));

      double guess_new_ion = guess_ion - (poly_ion / poly_der_ion);

      if (fabs(guess_ion - guess_new_ion) < pow(10.0, -8.0)) {
        iter_ion = 100;
      }
      else {
        iter_ion += 1;
        guess_ion = guess_new_ion;
      }
    }

    double W_ion = 0.5 * C0_ion * guess_ion * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((1.0 - (C_ion * guess_ion)) /
      ((C0_ion * C0_ion) * (guess_ion * guess_ion))))));
    double h_ion = 1.0 / (C_ion * guess_ion);

    double rho_ion = D_ion / W_ion; 
    double vx_ion = momx_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vy_ion = momy_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vz_ion = momz_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double p_ion = (rho_ion * h_ion * (W_ion * W_ion)) - D_ion - Etot_ion;

    if (rho_ion < pow(10.0, -8.0)) {
      rho_ion = pow(10.0, -8.0);
    }
    if (p_ion < pow(10.0, -8.0)) {
      p_ion = pow(10.0, -8.0);
    }
    
    double D[3], B[3];
    D[0] = Dx; D[1] = Dy; D[2] = Dz;
    B[0] = Bx; B[1] = By; B[2] = Bz;

    double covD[3], covB[3];
    for (int i = 0; i < 3; i++) {
      covD[i] = 0.0;
      covB[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        covD[i] += spatial_metric[i][j] * D[j];
        covB[i] += spatial_metric[i][j] * B[j];
      }
    }

    fluid_new[6] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * covD[0];
    fluid_new[6] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * ((vy_ion * covB[2]) - (vz_ion * covB[1]));
    fluid_new[7] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * covD[1];
    fluid_new[7] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * ((vz_ion * covB[0]) - (vx_ion * covB[2]));
    fluid_new[8] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * covD[2];
    fluid_new[8] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * ((vx_ion * covB[1]) - (vy_ion * covB[0]));

    fluid_new[9] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vx_ion * covD[0];
    fluid_new[9] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vy_ion * covD[1];
    fluid_new[9] += (lapse * sqrt(spatial_det)) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vz_ion * covD[2];
  }
}

void
explicit_gr_twofluid_source_update_em_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, const double gas_gamma_ion, const double mass_elc,
  const double charge_elc, const double mass_ion, const double charge_ion, double e_fact, double t_curr, const double dt, double* fluid_old, double* fluid_new)
{
  for (int i = 0; i < 84; i++) {
    fluid_new[i] = fluid_old[i];
  }

  double Dx = fluid_old[10], Dy = fluid_old[11], Dz = fluid_old[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[22]; spatial_metric[0][1] = fluid_old[23]; spatial_metric[0][2] = fluid_old[24];
  spatial_metric[1][0] = fluid_old[25]; spatial_metric[1][1] = fluid_old[26]; spatial_metric[1][2] = fluid_old[27];
  spatial_metric[2][0] = fluid_old[28]; spatial_metric[2][1] = fluid_old[29]; spatial_metric[2][2] = fluid_old[30];

  bool in_excision_region = false;
  if (fluid_old[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

    double D_elc = fluid_old[0] / sqrt(spatial_det);
    double momx_elc = fluid_old[1] / sqrt(spatial_det);
    double momy_elc = fluid_old[2] / sqrt(spatial_det);
    double momz_elc = fluid_old[3] / sqrt(spatial_det);
    double Etot_elc = fluid_old[4] / sqrt(spatial_det);

    double C_elc = D_elc / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    double C0_elc = (D_elc + Etot_elc) / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) +
      (momz_elc * momz_elc)));
    if (((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)) < pow(10.0, -8.0)) {
      C_elc = D_elc / sqrt(pow(10.0, -8.0));
      C0_elc = (D_elc + Etot_elc) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_elc = -1.0 / (gas_gamma_elc * gas_gamma_elc);
    double alpha1_elc = -2.0 * C_elc * ((gas_gamma_elc - 1.0) / (gas_gamma_elc * gas_gamma_elc));
    double alpha2_elc = ((gas_gamma_elc - 2.0) / gas_gamma_elc) * ((C0_elc * C0_elc) - 1.0) + 1.0 - (C_elc * C_elc) *
      ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((gas_gamma_elc - 1.0) / gas_gamma_elc);
    double alpha4_elc = (C0_elc * C0_elc) - 1.0;
    double eta_elc = 2.0 * C_elc *((gas_gamma_elc - 1.0) / gas_gamma_elc);

    double guess_elc = 1.0;
    int iter_elc = 0;

    while (iter_elc < 100) {
      double poly_elc = (alpha4_elc * (guess_elc * guess_elc * guess_elc) * (guess_elc - eta_elc)) + (alpha2_elc * (guess_elc * guess_elc)) +
        (alpha1_elc * guess_elc) + alpha0_elc;
      double poly_der_elc = alpha1_elc + (2.0 * alpha2_elc * guess_elc) + (4.0 * alpha4_elc * (guess_elc * guess_elc * guess_elc)) -
        (3.0 * eta_elc * alpha4_elc * (guess_elc * guess_elc));

      double guess_new_elc = guess_elc - (poly_elc / poly_der_elc);

      if (fabs(guess_elc - guess_new_elc) < pow(10.0, -8.0)) {
        iter_elc = 100;
      }
      else {
        iter_elc += 1;
        guess_elc = guess_new_elc;
      }
    }

    double W_elc = 0.5 * C0_elc * guess_elc * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((1.0 - (C_elc * guess_elc)) /
      ((C0_elc * C0_elc) * (guess_elc * guess_elc))))));
    double h_elc = 1.0 / (C_elc * guess_elc);

    double rho_elc = D_elc / W_elc; 
    double vx_elc = momx_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vy_elc = momy_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vz_elc = momz_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double p_elc = (rho_elc * h_elc * (W_elc * W_elc)) - D_elc - Etot_elc;

    if (rho_elc < pow(10.0, -8.0)) {
      rho_elc = pow(10.0, -8.0);
    }
    if (p_elc < pow(10.0, -8.0)) {
      p_elc = pow(10.0, -8.0);
    }

    double D_ion = fluid_old[5] / sqrt(spatial_det);
    double momx_ion = fluid_old[6] / sqrt(spatial_det);
    double momy_ion = fluid_old[7] / sqrt(spatial_det);
    double momz_ion = fluid_old[8] / sqrt(spatial_det);
    double Etot_ion = fluid_old[9] / sqrt(spatial_det);

    double C_ion = D_ion / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    double C0_ion = (D_ion + Etot_ion) / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) +
      (momz_ion * momz_ion)));
    if (((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)) < pow(10.0, -8.0)) {
      C_ion = D_ion / sqrt(pow(10.0, -8.0));
      C0_ion = (D_ion + Etot_ion) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_ion = -1.0 / (gas_gamma_ion * gas_gamma_ion);
    double alpha1_ion = -2.0 * C_ion * ((gas_gamma_ion - 1.0) / (gas_gamma_ion * gas_gamma_ion));
    double alpha2_ion = ((gas_gamma_ion - 2.0) / gas_gamma_ion) * ((C0_ion * C0_ion) - 1.0) + 1.0 - (C_ion * C_ion) *
      ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((gas_gamma_ion - 1.0) / gas_gamma_ion);
    double alpha4_ion = (C0_ion * C0_ion) - 1.0;
    double eta_ion = 2.0 * C_ion *((gas_gamma_ion - 1.0) / gas_gamma_ion);

    double guess_ion = 1.0;
    int iter_ion = 0;

    while (iter_ion < 100) {
      double poly_ion = (alpha4_ion * (guess_ion * guess_ion * guess_ion) * (guess_ion - eta_ion)) + (alpha2_ion * (guess_ion * guess_ion)) +
        (alpha1_ion * guess_ion) + alpha0_ion;
      double poly_der_ion = alpha1_ion + (2.0 * alpha2_ion * guess_ion) + (4.0 * alpha4_ion * (guess_ion * guess_ion * guess_ion)) -
        (3.0 * eta_ion * alpha4_ion * (guess_ion * guess_ion));

      double guess_new_ion = guess_ion - (poly_ion / poly_der_ion);

      if (fabs(guess_ion - guess_new_ion) < pow(10.0, -8.0)) {
        iter_ion = 100;
      }
      else {
        iter_ion += 1;
        guess_ion = guess_new_ion;
      }
    }

    double W_ion = 0.5 * C0_ion * guess_ion * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((1.0 - (C_ion * guess_ion)) /
      ((C0_ion * C0_ion) * (guess_ion * guess_ion))))));
    double h_ion = 1.0 / (C_ion * guess_ion);

    double rho_ion = D_ion / W_ion; 
    double vx_ion = momx_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vy_ion = momy_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vz_ion = momz_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double p_ion = (rho_ion * h_ion * (W_ion * W_ion)) - D_ion - Etot_ion;

    if (rho_ion < pow(10.0, -8.0)) {
      rho_ion = pow(10.0, -8.0);
    }
    if (p_ion < pow(10.0, -8.0)) {
      p_ion = pow(10.0, -8.0);
    }

    fluid_new[10] -= (4.0 * M_PI) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vx_elc;
    fluid_new[10] -= (4.0 * M_PI) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vx_ion;
    fluid_new[11] -= (4.0 * M_PI) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vy_elc;
    fluid_new[11] -= (4.0 * M_PI) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vy_ion;
    fluid_new[12] -= (4.0 * M_PI) * dt * (charge_elc / mass_elc) * rho_elc * W_elc * vz_elc;
    fluid_new[12] -= (4.0 * M_PI) * dt * (charge_ion / mass_ion) * rho_ion * W_ion * vz_ion;

    fluid_new[16] += dt * e_fact * (charge_elc / mass_elc) * rho_elc;
    fluid_new[16] += dt * e_fact * (charge_ion / mass_ion) * rho_ion;
  }
}

void
explicit_gr_twofluid_source_update_elc_spacetime_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[18];
  double shift_x = fluid_old[19];
  double shift_y = fluid_old[20];
  double shift_z = fluid_old[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[22]; spatial_metric[0][1] = fluid_old[23]; spatial_metric[0][2] = fluid_old[24];
  spatial_metric[1][0] = fluid_old[25]; spatial_metric[1][1] = fluid_old[26]; spatial_metric[1][2] = fluid_old[27];
  spatial_metric[2][0] = fluid_old[28]; spatial_metric[2][1] = fluid_old[29]; spatial_metric[2][2] = fluid_old[30];

  double inv_spatial_metric[3][3];
  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[31]; extrinsic_curvature[0][1] = fluid_old[32]; extrinsic_curvature[0][2] = fluid_old[33];
  extrinsic_curvature[1][0] = fluid_old[34]; extrinsic_curvature[1][1] = fluid_old[35]; extrinsic_curvature[1][2] = fluid_old[36];
  extrinsic_curvature[2][0] = fluid_old[37]; extrinsic_curvature[2][1] = fluid_old[38]; extrinsic_curvature[2][2] = fluid_old[39];

  bool in_excision_region = false;
  if (fluid_old[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double D_elc = fluid_old[0] / sqrt(spatial_det);
    double momx_elc = fluid_old[1] / sqrt(spatial_det);
    double momy_elc = fluid_old[2] / sqrt(spatial_det);
    double momz_elc = fluid_old[3] / sqrt(spatial_det);
    double Etot_elc = fluid_old[4] / sqrt(spatial_det);

    double C_elc = D_elc / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    double C0_elc = (D_elc + Etot_elc) / sqrt(((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)));
    if (((Etot_elc + D_elc) * (Etot_elc + D_elc)) - ((momx_elc * momx_elc) + (momy_elc * momy_elc) + (momz_elc * momz_elc)) < pow(10.0, -8.0)) {
      C_elc = D_elc / sqrt(pow(10.0, -8.0));
      C0_elc = (D_elc + Etot_elc) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_elc = -1.0 / (gas_gamma_elc * gas_gamma_elc);
    double alpha1_elc = -2.0 * C_elc * ((gas_gamma_elc - 1.0) / (gas_gamma_elc * gas_gamma_elc));
    double alpha2_elc = ((gas_gamma_elc - 2.0) / gas_gamma_elc) * ((C0_elc * C0_elc) - 1.0) + 1.0 - (C_elc * C_elc) *
      ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((gas_gamma_elc - 1.0) / gas_gamma_elc);
    double alpha4_elc = (C0_elc * C0_elc) - 1.0;
    double eta_elc = 2.0 * C_elc *((gas_gamma_elc - 1.0) / gas_gamma_elc);

    double guess_elc = 1.0;
    int iter_elc = 0;

    while (iter_elc < 100) {
      double poly_elc = (alpha4_elc * (guess_elc * guess_elc * guess_elc) * (guess_elc - eta_elc)) + (alpha2_elc * (guess_elc * guess_elc)) +
        (alpha1_elc * guess_elc) + alpha0_elc;
      double poly_der_elc = alpha1_elc + (2.0 * alpha2_elc * guess_elc) + (4.0 * alpha4_elc * (guess_elc * guess_elc * guess_elc)) -
        (3.0 * eta_elc * alpha4_elc * (guess_elc * guess_elc));

      double guess_new_elc = guess_elc - (poly_elc / poly_der_elc);

      if (fabs(guess_elc - guess_new_elc) < pow(10.0, -8.0)) {
        iter_elc = 100;
      }
      else {
        iter_elc += 1;
        guess_elc = guess_new_elc;
      }
    }

    double W_elc = 0.5 * C0_elc * guess_elc * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_elc - 1.0) / gas_gamma_elc) * ((1.0 - (C_elc * guess_elc)) /
      ((C0_elc * C0_elc) * (guess_elc * guess_elc))))));
    double h_elc = 1.0 / (C_elc * guess_elc);

    double rho_elc = D_elc / W_elc;
    double vx_elc = momx_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vy_elc = momy_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double vz_elc = momz_elc / (rho_elc * h_elc * (W_elc * W_elc));
    double p_elc = (rho_elc * h_elc * (W_elc * W_elc)) - D_elc - Etot_elc;

    if (rho_elc < pow(10.0, -8.0)) {
      rho_elc = pow(10.0, -8.0);
    }
    if (p_elc < pow(10.0, -8.0)) {
      p_elc = pow(10.0, -8.0);
    }

    double spacetime_vel_elc[4];
    spacetime_vel_elc[0] = W_elc / lapse;
    spacetime_vel_elc[1] = (W_elc * vx_elc) - (shift_x * (W_elc / lapse));
    spacetime_vel_elc[2] = (W_elc * vy_elc) - (shift_y * (W_elc / lapse));
    spacetime_vel_elc[3] = (W_elc * vz_elc) - (shift_z * (W_elc / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double stress_energy_elc[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy_elc[i][j] = (rho_elc * h_elc * spacetime_vel_elc[i] * spacetime_vel_elc[j]) + (p_elc * inv_spacetime_metric[i][j]);
      }
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[41];
    lapse_der[1] = fluid_old[42];
    lapse_der[2] = fluid_old[43];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[44]; shift_der[0][1] = fluid_old[45]; shift_der[0][2] = fluid_old[46];
    shift_der[1][0] = fluid_old[47]; shift_der[1][1] = fluid_old[48]; shift_der[1][2] = fluid_old[49];
    shift_der[2][0] = fluid_old[50]; shift_der[2][1] = fluid_old[51]; shift_der[2][2] = fluid_old[52];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[53]; spatial_metric_der[0][0][1] = fluid_old[54]; spatial_metric_der[0][0][2] = fluid_old[55];
    spatial_metric_der[0][1][0] = fluid_old[56]; spatial_metric_der[0][1][1] = fluid_old[57]; spatial_metric_der[0][1][2] = fluid_old[58];
    spatial_metric_der[0][2][0] = fluid_old[59]; spatial_metric_der[0][2][1] = fluid_old[60]; spatial_metric_der[0][2][2] = fluid_old[61];

    spatial_metric_der[1][0][0] = fluid_old[62]; spatial_metric_der[1][0][1] = fluid_old[63]; spatial_metric_der[1][0][2] = fluid_old[64];
    spatial_metric_der[1][1][0] = fluid_old[65]; spatial_metric_der[1][1][1] = fluid_old[66]; spatial_metric_der[1][1][2] = fluid_old[67];
    spatial_metric_der[1][2][0] = fluid_old[68]; spatial_metric_der[1][2][1] = fluid_old[69]; spatial_metric_der[1][2][2] = fluid_old[70];

    spatial_metric_der[2][0][0] = fluid_old[71]; spatial_metric_der[2][0][1] = fluid_old[72]; spatial_metric_der[2][0][2] = fluid_old[73];
    spatial_metric_der[2][1][0] = fluid_old[74]; spatial_metric_der[2][1][1] = fluid_old[75]; spatial_metric_der[2][1][2] = fluid_old[76];
    spatial_metric_der[2][2][0] = fluid_old[77]; spatial_metric_der[2][2][1] = fluid_old[78]; spatial_metric_der[2][2][2] = fluid_old[79];

    double mom_elc[3];
    mom_elc[0] = (rho_elc + p_elc) * (W_elc * W_elc) * vx_elc;
    mom_elc[1] = (rho_elc + p_elc) * (W_elc * W_elc) * vy_elc;
    mom_elc[2] = (rho_elc + p_elc) * (W_elc * W_elc) * vz_elc;

    for (int i = 0; i < 84; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[4] += dt * (stress_energy_elc[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (2.0 * stress_energy_elc[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (stress_energy_elc[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[4] -= dt * (stress_energy_elc[0][0] * shift[i] * lapse_der[i]);
      fluid_new[4] -= dt * (stress_energy_elc[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[1 + j] -= dt * (stress_energy_elc[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[1 + j] += dt * (0.5 * stress_energy_elc[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[1 + j] += dt * (0.5 * stress_energy_elc[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[1 + j] += dt * ((mom_elc[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[1 + j] += dt * (stress_energy_elc[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_twofluid_source_update_ion_spacetime_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_ion, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[18];
  double shift_x = fluid_old[19];
  double shift_y = fluid_old[20];
  double shift_z = fluid_old[21];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[22]; spatial_metric[0][1] = fluid_old[23]; spatial_metric[0][2] = fluid_old[24];
  spatial_metric[1][0] = fluid_old[25]; spatial_metric[1][1] = fluid_old[26]; spatial_metric[1][2] = fluid_old[27];
  spatial_metric[2][0] = fluid_old[28]; spatial_metric[2][1] = fluid_old[29]; spatial_metric[2][2] = fluid_old[30];

  double inv_spatial_metric[3][3];
  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[31]; extrinsic_curvature[0][1] = fluid_old[32]; extrinsic_curvature[0][2] = fluid_old[33];
  extrinsic_curvature[1][0] = fluid_old[34]; extrinsic_curvature[1][1] = fluid_old[35]; extrinsic_curvature[1][2] = fluid_old[36];
  extrinsic_curvature[2][0] = fluid_old[37]; extrinsic_curvature[2][1] = fluid_old[38]; extrinsic_curvature[2][2] = fluid_old[39];

  bool in_excision_region = false;
  if (fluid_old[40] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double D_ion = fluid_old[5] / sqrt(spatial_det);
    double momx_ion = fluid_old[6] / sqrt(spatial_det);
    double momy_ion = fluid_old[7] / sqrt(spatial_det);
    double momz_ion = fluid_old[8] / sqrt(spatial_det);
    double Etot_ion = fluid_old[9] / sqrt(spatial_det);

    double C_ion = D_ion / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    double C0_ion = (D_ion + Etot_ion) / sqrt(((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)));
    if (((Etot_ion + D_ion) * (Etot_ion + D_ion)) - ((momx_ion * momx_ion) + (momy_ion * momy_ion) + (momz_ion * momz_ion)) < pow(10.0, -8.0)) {
      C_ion = D_ion / sqrt(pow(10.0, -8.0));
      C0_ion = (D_ion + Etot_ion) / sqrt(pow(10.0, -8.0));
    }

    double alpha0_ion = -1.0 / (gas_gamma_ion * gas_gamma_ion);
    double alpha1_ion = -2.0 * C_ion * ((gas_gamma_ion - 1.0) / (gas_gamma_ion * gas_gamma_ion));
    double alpha2_ion = ((gas_gamma_ion - 2.0) / gas_gamma_ion) * ((C0_ion * C0_ion) - 1.0) + 1.0 - (C_ion * C_ion) *
      ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((gas_gamma_ion - 1.0) / gas_gamma_ion);
    double alpha4_ion = (C0_ion * C0_ion) - 1.0;
    double eta_ion = 2.0 * C_ion *((gas_gamma_ion - 1.0) / gas_gamma_ion);

    double guess_ion = 1.0;
    int iter_ion = 0;

    while (iter_ion < 100) {
      double poly_ion = (alpha4_ion * (guess_ion * guess_ion * guess_ion) * (guess_ion - eta_ion)) + (alpha2_ion * (guess_ion * guess_ion)) +
        (alpha1_ion * guess_ion) + alpha0_ion;
      double poly_der_ion = alpha1_ion + (2.0 * alpha2_ion * guess_ion) + (4.0 * alpha4_ion * (guess_ion * guess_ion * guess_ion)) -
        (3.0 * eta_ion * alpha4_ion * (guess_ion * guess_ion));

      double guess_new_ion = guess_ion - (poly_ion / poly_der_ion);

      if (fabs(guess_ion - guess_new_ion) < pow(10.0, -8.0)) {
        iter_ion = 100;
      }
      else {
        iter_ion += 1;
        guess_ion = guess_new_ion;
      }
    }

    double W_ion = 0.5 * C0_ion * guess_ion * (1.0 + sqrt(1.0 + (4.0 * ((gas_gamma_ion - 1.0) / gas_gamma_ion) * ((1.0 - (C_ion * guess_ion)) /
      ((C0_ion * C0_ion) * (guess_ion * guess_ion))))));
    double h_ion = 1.0 / (C_ion * guess_ion);

    double rho_ion = D_ion / W_ion;
    double vx_ion = momx_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vy_ion = momy_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double vz_ion = momz_ion / (rho_ion * h_ion * (W_ion * W_ion));
    double p_ion = (rho_ion * h_ion * (W_ion * W_ion)) - D_ion - Etot_ion;

    if (rho_ion < pow(10.0, -8.0)) {
      rho_ion = pow(10.0, -8.0);
    }
    if (p_ion < pow(10.0, -8.0)) {
      p_ion = pow(10.0, -8.0);
    }

    double spacetime_vel_ion[4];
    spacetime_vel_ion[0] = W_ion / lapse;
    spacetime_vel_ion[1] = (W_ion * vx_ion) - (shift_x * (W_ion / lapse));
    spacetime_vel_ion[2] = (W_ion * vy_ion) - (shift_y * (W_ion / lapse));
    spacetime_vel_ion[3] = (W_ion * vz_ion) - (shift_z * (W_ion / lapse));

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double stress_energy_ion[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy_ion[i][j] = (rho_ion * h_ion * spacetime_vel_ion[i] * spacetime_vel_ion[j]) + (p_ion * inv_spacetime_metric[i][j]);
      }
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[41];
    lapse_der[1] = fluid_old[42];
    lapse_der[2] = fluid_old[43];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[44]; shift_der[0][1] = fluid_old[45]; shift_der[0][2] = fluid_old[46];
    shift_der[1][0] = fluid_old[47]; shift_der[1][1] = fluid_old[48]; shift_der[1][2] = fluid_old[49];
    shift_der[2][0] = fluid_old[50]; shift_der[2][1] = fluid_old[51]; shift_der[2][2] = fluid_old[52];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[53]; spatial_metric_der[0][0][1] = fluid_old[54]; spatial_metric_der[0][0][2] = fluid_old[55];
    spatial_metric_der[0][1][0] = fluid_old[56]; spatial_metric_der[0][1][1] = fluid_old[57]; spatial_metric_der[0][1][2] = fluid_old[58];
    spatial_metric_der[0][2][0] = fluid_old[59]; spatial_metric_der[0][2][1] = fluid_old[60]; spatial_metric_der[0][2][2] = fluid_old[61];

    spatial_metric_der[1][0][0] = fluid_old[62]; spatial_metric_der[1][0][1] = fluid_old[63]; spatial_metric_der[1][0][2] = fluid_old[64];
    spatial_metric_der[1][1][0] = fluid_old[65]; spatial_metric_der[1][1][1] = fluid_old[66]; spatial_metric_der[1][1][2] = fluid_old[67];
    spatial_metric_der[1][2][0] = fluid_old[68]; spatial_metric_der[1][2][1] = fluid_old[69]; spatial_metric_der[1][2][2] = fluid_old[70];

    spatial_metric_der[2][0][0] = fluid_old[71]; spatial_metric_der[2][0][1] = fluid_old[72]; spatial_metric_der[2][0][2] = fluid_old[73];
    spatial_metric_der[2][1][0] = fluid_old[74]; spatial_metric_der[2][1][1] = fluid_old[75]; spatial_metric_der[2][1][2] = fluid_old[76];
    spatial_metric_der[2][2][0] = fluid_old[77]; spatial_metric_der[2][2][1] = fluid_old[78]; spatial_metric_der[2][2][2] = fluid_old[79];

    double mom_ion[3];
    mom_ion[0] = (rho_ion + p_ion) * (W_ion * W_ion) * vx_ion;
    mom_ion[1] = (rho_ion + p_ion) * (W_ion * W_ion) * vy_ion;
    mom_ion[2] = (rho_ion + p_ion) * (W_ion * W_ion) * vz_ion;

    for (int i = 0; i < 84; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[9] += dt * (stress_energy_ion[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[9] += dt * (2.0 * stress_energy_ion[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[9] += dt * (stress_energy_ion[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[9] -= dt * (stress_energy_ion[0][0] * shift[i] * lapse_der[i]);
      fluid_new[9] -= dt * (stress_energy_ion[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[6 + j] -= dt * (stress_energy_ion[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[6 + j] += dt * (0.5 * stress_energy_ion[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[6 + j] += dt * (0.5 * stress_energy_ion[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[6 + j] += dt * ((mom_ion[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[6 + j] += dt * (stress_energy_ion[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 84; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_twofluid_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double mass_elc = mom_em->gr_twofluid_mass_elc;
  double mass_ion = mom_em->gr_twofluid_mass_ion;
  double charge_elc = mom_em->gr_twofluid_charge_elc;
  double charge_ion = mom_em->gr_twofluid_charge_ion;
  double gas_gamma_elc = mom_em->gr_twofluid_gas_gamma_elc;
  double gas_gamma_ion = mom_em->gr_twofluid_gas_gamma_ion;
  double e_fact = mom_em->gr_twofluid_e_fact;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new_elc[84], f_new_elc_spacetime[84], f_new_ion[84], f_new_ion_spacetime[84], f_new[84];
    double f_stage1[84], f_stage2[84], f_old[84];

    for (int j = 0; j < 84; j++) {
      f_old[j] = f[j];
    }

    explicit_gr_twofluid_source_update_elc_euler(mom_em, gas_gamma_elc, mass_elc, charge_elc, t_curr, dt, f_old, f_new_elc);
    explicit_gr_twofluid_source_update_elc_spacetime_euler(mom_em, gas_gamma_elc, t_curr, dt, f_new_elc, f_new_elc_spacetime);
    explicit_gr_twofluid_source_update_ion_euler(mom_em, gas_gamma_ion, mass_ion, charge_ion, t_curr, dt, f_old, f_new_ion);
    explicit_gr_twofluid_source_update_ion_spacetime_euler(mom_em, gas_gamma_ion, t_curr, dt, f_new_ion, f_new_ion_spacetime);
    explicit_gr_twofluid_source_update_em_euler(mom_em, gas_gamma_elc, gas_gamma_ion, mass_elc, charge_elc, mass_ion, charge_ion, e_fact, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 5; j++) {
      f_stage1[j] = f_new_elc_spacetime[j];
    }
    for (int j = 5; j < 10; j++) {
      f_stage1[j] = f_new_ion_spacetime[j];
    }
    for (int j = 10; j < 84; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_gr_twofluid_source_update_elc_euler(mom_em, gas_gamma_elc, mass_elc, charge_elc, t_curr + dt, dt, f_stage1, f_new_elc);
    explicit_gr_twofluid_source_update_elc_spacetime_euler(mom_em, gas_gamma_elc, t_curr + dt, dt, f_new_elc, f_new_elc_spacetime);
    explicit_gr_twofluid_source_update_ion_euler(mom_em, gas_gamma_ion, mass_ion, charge_ion, t_curr + dt, dt, f_stage1, f_new_ion);
    explicit_gr_twofluid_source_update_ion_spacetime_euler(mom_em, gas_gamma_ion, t_curr + dt, dt, f_new_ion, f_new_ion_spacetime);
    explicit_gr_twofluid_source_update_em_euler(mom_em, gas_gamma_elc, gas_gamma_ion, mass_elc, charge_elc, mass_ion, charge_ion, e_fact, t_curr + dt, dt,
      f_stage1, f_new);
    for (int j = 0; j < 5; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new_elc_spacetime[j]);
    }
    for (int j = 5; j < 10; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new_ion_spacetime[j]);
    }
    for (int j = 10; j < 84; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_gr_twofluid_source_update_elc_euler(mom_em, gas_gamma_elc, mass_elc, charge_elc, t_curr + (0.5 * dt), dt, f_stage2, f_new_elc);
    explicit_gr_twofluid_source_update_elc_spacetime_euler(mom_em, gas_gamma_elc, t_curr + (0.5 * dt), dt, f_new_elc, f_new_elc_spacetime);
    explicit_gr_twofluid_source_update_ion_euler(mom_em, gas_gamma_ion, mass_ion, charge_ion, t_curr + (0.5 * dt), dt, f_stage2, f_new_ion);
    explicit_gr_twofluid_source_update_ion_spacetime_euler(mom_em, gas_gamma_ion, t_curr + (0.5 * dt), dt, f_new_ion, f_new_ion_spacetime);
    explicit_gr_twofluid_source_update_em_euler(mom_em, gas_gamma_elc, gas_gamma_ion, mass_elc, charge_elc, mass_ion, charge_ion, e_fact,
      t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 5; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new_elc_spacetime[j]);
    }
    for (int j = 5; j < 10; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new_ion_spacetime[j]);
    }
    for (int j = 10; j < 84; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_vacuum_einstein_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double excision_threshold,
  const enum gkyl_spacetime_slicing spacetime_slicing, const enum gkyl_spacetime_evolution spacetime_evolution, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[0]; spatial_metric[0][1] = fluid_old[1]; spatial_metric[0][2] = fluid_old[2];
  spatial_metric[1][0] = fluid_old[3]; spatial_metric[1][1] = fluid_old[4]; spatial_metric[1][2] = fluid_old[5];
  spatial_metric[2][0] = fluid_old[6]; spatial_metric[2][1] = fluid_old[7]; spatial_metric[2][2] = fluid_old[8];

  double lapse = fluid_old[9];

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[10]; extrinsic_curvature[0][1] = fluid_old[11]; extrinsic_curvature[0][2] = fluid_old[12];
  extrinsic_curvature[1][0] = fluid_old[13]; extrinsic_curvature[1][1] = fluid_old[14]; extrinsic_curvature[1][2] = fluid_old[15];
  extrinsic_curvature[2][0] = fluid_old[16]; extrinsic_curvature[2][1] = fluid_old[17]; extrinsic_curvature[2][2] = fluid_old[18];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = fluid_old[19]; spatial_metric_der[0][0][1] = fluid_old[20]; spatial_metric_der[0][0][2] = fluid_old[21];
  spatial_metric_der[0][1][0] = fluid_old[22]; spatial_metric_der[0][1][1] = fluid_old[23]; spatial_metric_der[0][1][2] = fluid_old[24];
  spatial_metric_der[0][2][0] = fluid_old[25]; spatial_metric_der[0][2][1] = fluid_old[26]; spatial_metric_der[0][2][2] = fluid_old[27];

  spatial_metric_der[1][0][0] = fluid_old[28]; spatial_metric_der[1][0][1] = fluid_old[29]; spatial_metric_der[1][0][2] = fluid_old[30];
  spatial_metric_der[1][1][0] = fluid_old[31]; spatial_metric_der[1][1][1] = fluid_old[32]; spatial_metric_der[1][1][2] = fluid_old[33];
  spatial_metric_der[1][2][0] = fluid_old[34]; spatial_metric_der[1][2][1] = fluid_old[35]; spatial_metric_der[1][2][2] = fluid_old[36];

  spatial_metric_der[2][0][0] = fluid_old[37]; spatial_metric_der[2][0][1] = fluid_old[38]; spatial_metric_der[2][0][2] = fluid_old[39];
  spatial_metric_der[2][1][0] = fluid_old[40]; spatial_metric_der[2][1][1] = fluid_old[41]; spatial_metric_der[2][1][2] = fluid_old[42];
  spatial_metric_der[2][2][0] = fluid_old[43]; spatial_metric_der[2][2][1] = fluid_old[44]; spatial_metric_der[2][2][2] = fluid_old[45];

  double lapse_der[3];
  lapse_der[0] = fluid_old[46];
  lapse_der[1] = fluid_old[47];
  lapse_der[2] = fluid_old[48];

  double aux_vect[3];
  aux_vect[0] = fluid_old[49];
  aux_vect[1] = fluid_old[50];
  aux_vect[2] = fluid_old[51];

  double shift_vect[3];
  shift_vect[0] = fluid_old[52];
  shift_vect[1] = fluid_old[53];
  shift_vect[2] = fluid_old[54];

  double shift_vect_der[3][3];
  shift_vect_der[0][0] = fluid_old[55]; shift_vect_der[0][1] = fluid_old[56]; shift_vect_der[0][2] = fluid_old[57];
  shift_vect_der[1][0] = fluid_old[58]; shift_vect_der[1][1] = fluid_old[59]; shift_vect_der[1][2] = fluid_old[60];
  shift_vect_der[2][0] = fluid_old[61]; shift_vect_der[2][1] = fluid_old[62]; shift_vect_der[2][2] = fluid_old[63];

  bool in_excision_region = false;
  if (lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double inv_spatial_metric[3][3];
    double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
      (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
      (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
    
    double trace = 0.0;
    for (int i = 0; i < 3; i++) {
      trace += spatial_metric[i][i];
    }

    double spatial_metric_sq[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        spatial_metric_sq[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
        }
      }
    }

    double sq_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      sq_trace += spatial_metric_sq[i][i];
    }

    double euclidean_metric[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) {
          euclidean_metric[i][j] = 1.0;
        }
        else {
          euclidean_metric[i][j] = 0.0;
        }
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spatial_metric[i][j] = (1.0 / spatial_det) *
          ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
      }
    }

    double evolution_func = 0.0;
    if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
      evolution_func = 0.0;
    }
    else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
      evolution_func = 1.0;
    }

    double extrinsic_curvature_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];
      }
    }

    double extrinsic_curvature_mixed[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        extrinsic_curvature_mixed[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          extrinsic_curvature_mixed[i][j] += inv_spatial_metric[l][j] * extrinsic_curvature[i][l];
        }
      }
    }

    double extrinsic_curvature_raised[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        extrinsic_curvature_raised[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            extrinsic_curvature_raised[i][j] += inv_spatial_metric[i][l] * inv_spatial_metric[m][j] * extrinsic_curvature[l][m];
          }
        }
      }
    }

    double shift_vect_der_lowered[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        shift_vect_der_lowered[i][j] = 0.0;

        for (int k = 0; k < 3; k++) {
          shift_vect_der_lowered[i][j] += spatial_metric[k][j] * shift_vect_der[i][k];
        }
      }
    }

    double shift_vect_der_switched[3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        shift_vect_der_switched[i][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            shift_vect_der_switched[i][k] += inv_spatial_metric[i][l] * spatial_metric[m][k] * shift_vect_der[l][m];
          }
        }
      }
    }

    double shift_vect_der_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      shift_vect_der_trace += shift_vect_der[i][i];
    }

    double symmetrized_shift[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        symmetrized_shift[i][j] = (1.0 / lapse) * (shift_vect_der_lowered[i][j] + shift_vect_der_lowered[j][i]);
      }
    }

    double slicing_func = 0.0;
    if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
      slicing_func = 0.0;
    }
    else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
      slicing_func = extrinsic_curvature_trace;
    }
    else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
      slicing_func = 2.0 * extrinsic_curvature_trace / lapse;
    }

    double spatial_metric_der_raised1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          spatial_metric_der_raised1[k][i][j] = 0.0;
          
          for (int l = 0; l < 3; l++) {
            spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
          }
        }
      }
    }

    double spatial_metric_der_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          spatial_metric_der_raised3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
          }
        }
      }
    }

    double spatial_metric_der_lowered1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          spatial_metric_der_lowered1[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              spatial_metric_der_lowered1[i][j][k] += inv_spatial_metric[j][l] * inv_spatial_metric[m][k] * spatial_metric_der[i][l][m];
            }
          }
        }
      }
    }

    double spatial_metric_der_lowered3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          spatial_metric_der_lowered3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              spatial_metric_der_lowered3[i][j][k] += inv_spatial_metric[i][l] * inv_spatial_metric[m][j] * spatial_metric_der[l][m][k];
            }
          }
        }
      }
    }

    double aux_vect_raised[3];
    for (int k = 0; k < 3; k++) {
      aux_vect_raised[k] = 0.0;
        
      for (int l = 0; l < 3; l++) {
        aux_vect_raised[k] += inv_spatial_metric[k][l] * aux_vect[l];
      }
    }

    double spatial_christoffel[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          spatial_christoffel[i][k][l] = 0.0;

          for (int m = 0; m < 3; m++) {
            spatial_christoffel[i][k][l] += inv_spatial_metric[i][m] * spatial_metric_der[l][m][k];
            spatial_christoffel[i][k][l] += inv_spatial_metric[i][m] * spatial_metric_der[k][m][l];
            spatial_christoffel[i][k][l] -= inv_spatial_metric[i][m] * spatial_metric_der[m][k][l];
          }
        }
      }
    }

    double spatial_metric_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        spatial_metric_source[i][j] = -2.0 * lapse * (extrinsic_curvature[i][j] - symmetrized_shift[i][j]);

        for (int r = 0; r < 3; r++) {
          spatial_metric_source[i][j] += 2.0 * shift_vect[r] * spatial_metric_der[r][i][j];
        }
      }
    }

    double lapse_source = -(lapse * lapse) * slicing_func;
    for (int r = 0; r < 3; r++) {
      lapse_source += lapse * shift_vect[r] * lapse_der[r];
    }

    double extrinsic_curvature_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        extrinsic_curvature_source[i][j] = 0.0;

        for (int r = 0; r < 3; r++) {
          extrinsic_curvature_source[i][j] += 2.0 * extrinsic_curvature[i][r] * shift_vect_der[j][r];
          extrinsic_curvature_source[i][j] += 2.0 * extrinsic_curvature[j][r] * shift_vect_der[i][r];
          extrinsic_curvature_source[i][j] -= 2.0 * extrinsic_curvature[i][j] * shift_vect_der[r][r];
        }

        for (int k = 0; k < 3; k++) {
          extrinsic_curvature_source[i][j] -= 2.0 * lapse * extrinsic_curvature_mixed[i][k] * extrinsic_curvature[k][j];
          extrinsic_curvature_source[i][j] += lapse * extrinsic_curvature_trace * extrinsic_curvature[i][j];

          for (int r = 0; r < 3; r++) {
            extrinsic_curvature_source[i][j] -= lapse * spatial_christoffel[k][r][i] * spatial_christoffel[r][k][j];

            extrinsic_curvature_source[i][j] += 2.0 * lapse * spatial_metric_der_raised3[i][k][r] * spatial_metric_der_raised3[r][j][k];
            extrinsic_curvature_source[i][j] += 2.0 * lapse * spatial_metric_der_raised3[j][k][r] * spatial_metric_der_raised3[r][i][k];
            extrinsic_curvature_source[i][j] += lapse * spatial_christoffel[k][k][r] * spatial_christoffel[r][i][j];

            extrinsic_curvature_source[i][j] -= lapse * (2.0 * spatial_metric_der_raised3[k][r][k] - lapse_der[r]) * (spatial_metric_der_raised3[i][j][r] + spatial_metric_der_raised3[j][i][r]);
          }

          extrinsic_curvature_source[i][j] += lapse * lapse_der[i] * (aux_vect[j] - (0.5 * spatial_metric_der_raised3[j][k][k]));
          extrinsic_curvature_source[i][j] += lapse * lapse_der[j] * (aux_vect[i] - (0.5 * spatial_metric_der_raised3[i][k][k]));

          extrinsic_curvature_source[i][j] -= lapse * evolution_func * aux_vect_raised[k] * spatial_metric_der[k][i][j];
        }

        for (int k = 0; k < 3; k++) {
          for (int r = 0; r < 3; r++) {
            for (int s = 0; s < 3; s++) {
              extrinsic_curvature_source[i][j] -= (0.25 * evolution_func * lapse * spatial_metric[i][j]) * spatial_metric_der_lowered1[k][r][s] * spatial_christoffel[k][r][s];
              extrinsic_curvature_source[i][j] += (0.25 * evolution_func * lapse * spatial_metric[i][j]) * spatial_metric_der_raised3[k][r][r] * spatial_metric_der_lowered3[k][s][s];
            }
          }

          extrinsic_curvature_source[i][j] -= (0.5 * evolution_func * lapse * spatial_metric[i][j]) * aux_vect_raised[k] * lapse_der[k];
        }

        for (int r = 0; r < 3; r++) {
          for (int s = 0; s < 3; s++) {
            extrinsic_curvature_source[i][j] += (0.25 * evolution_func * lapse * spatial_metric[i][j]) * extrinsic_curvature_raised[r][s] * extrinsic_curvature[r][s];
          }
        }

        extrinsic_curvature_source[i][j] -= (0.25 * evolution_func * lapse * spatial_metric[i][j]) * (extrinsic_curvature_trace * extrinsic_curvature_trace);
      }
    }

    double aux_vect_source[3];
    for (int i = 0; i < 3; i++) {
      aux_vect_source[i] = 0.0;

      for (int r = 0; r < 3; r++) {
        aux_vect_source[i] += lapse * lapse_der[r] * extrinsic_curvature_mixed[i][r];
        
        if (i == r) {
          aux_vect_source[i] -= lapse * lapse_der[r] * extrinsic_curvature_trace;
        }

        for (int s = 0; s < 3; s++) {
          aux_vect_source[i] += lapse * extrinsic_curvature_mixed[s][r] * spatial_metric_der_raised3[i][r][s];
          aux_vect_source[i] -= 2.0 * lapse * extrinsic_curvature_mixed[s][r] * spatial_metric_der_raised3[r][i][s];

          aux_vect_source[i] -= lapse * extrinsic_curvature_mixed[i][r] * spatial_metric_der_raised3[r][s][s];
          aux_vect_source[i] += 2.0 * lapse * extrinsic_curvature_mixed[i][r] * spatial_metric_der_raised3[s][r][s];
        }

        aux_vect_source[i] += 2.0 * shift_vect_der[i][r] * aux_vect[r];

        if (i == r) {
          aux_vect_source[i] -= 2.0 * shift_vect_der_trace * aux_vect[r];
        }

        for (int s = 0; s < 3; s++) {
          aux_vect_source[i] += 2.0 * spatial_metric_der_raised3[r][i][s] * shift_vect_der_switched[r][s];

          if (i == s) {
            for (int j = 0; j < 3; j++) {
              aux_vect_source[i] -= 2.0 * spatial_metric_der_raised1[j][j][r] * shift_vect_der_switched[r][s];
            }
          }
        }
      }
    }

    for (int i = 0; i < 64; i++) {
      fluid_new[i] = fluid_old[i];
    }

    fluid_new[0] += dt * spatial_metric_source[0][0]; fluid_new[1] += dt * spatial_metric_source[0][1]; fluid_new[2] += dt * spatial_metric_source[0][2];
    fluid_new[3] += dt * spatial_metric_source[1][0]; fluid_new[4] += dt * spatial_metric_source[1][1]; fluid_new[5] += dt * spatial_metric_source[1][2];
    fluid_new[6] += dt * spatial_metric_source[2][0]; fluid_new[7] += dt * spatial_metric_source[2][1]; fluid_new[8] += dt * spatial_metric_source[2][2];

    fluid_new[9] += dt * lapse_source;

    fluid_new[10] += dt * extrinsic_curvature_source[0][0]; fluid_new[11] += dt * extrinsic_curvature_source[0][1]; fluid_new[12] += dt * extrinsic_curvature_source[0][2];
    fluid_new[13] += dt * extrinsic_curvature_source[1][0]; fluid_new[14] += dt * extrinsic_curvature_source[1][1]; fluid_new[15] += dt * extrinsic_curvature_source[1][2];
    fluid_new[16] += dt * extrinsic_curvature_source[2][0]; fluid_new[17] += dt * extrinsic_curvature_source[2][1]; fluid_new[18] += dt * extrinsic_curvature_source[2][2];

    fluid_new[49] += dt * aux_vect_source[0];
    fluid_new[50] += dt * aux_vect_source[1];
    fluid_new[51] += dt * aux_vect_source[2];

    if (fluid_new[9] < excision_threshold) {
      for (int i = 0; i < 64; i++) {
        fluid_new[i] = 0.0;
      }
    }
  }
  else {
    for (int i = 0; i < 64; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_vacuum_einstein_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double excision_threshold = mom_em->vacuum_einstein_excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = mom_em->vacuum_einstein_spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = mom_em->vacuum_einstein_spacetime_evolution;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[64], f_stage1[64], f_stage2[64], f_old[64];

    for (int j = 0; j < 64; j++) {
      f_old[j] = f[j];
    }

    explicit_vacuum_einstein_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 64; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_vacuum_einstein_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 64; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_vacuum_einstein_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 64; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_vacuum_einstein_conformal_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double excision_threshold,
  const enum gkyl_spacetime_slicing spacetime_slicing, const enum gkyl_spacetime_evolution spacetime_evolution, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double bssn_conformal_fact = fluid_old[64];
  double conformal_fact = 1.0 / sqrt(bssn_conformal_fact);

  double bssn_conformal_fact_der[3];
  bssn_conformal_fact_der[0] = fluid_old[65]; bssn_conformal_fact_der[1] = fluid_old[66]; bssn_conformal_fact_der[2] = fluid_old[67];
  double conformal_fact_der[3];
  for (int i = 0; i < 3; i++) {
    conformal_fact_der[i] = -(0.5 * bssn_conformal_fact_der[i]) / (pow(bssn_conformal_fact, 1.5) * conformal_fact);
  }

  double bssn_conformal_fact_der2[3][3];
  bssn_conformal_fact_der2[0][0] = fluid_old[68]; bssn_conformal_fact_der2[0][1] = fluid_old[69]; bssn_conformal_fact_der2[0][2] = fluid_old[70];
  bssn_conformal_fact_der2[1][0] = fluid_old[71]; bssn_conformal_fact_der2[1][1] = fluid_old[72]; bssn_conformal_fact_der2[1][2] = fluid_old[73];
  bssn_conformal_fact_der2[2][0] = fluid_old[74]; bssn_conformal_fact_der2[2][1] = fluid_old[75]; bssn_conformal_fact_der2[2][2] = fluid_old[76];
  double conformal_fact_der2[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_fact_der2[i][j] = -(0.5 * bssn_conformal_fact_der2[i][j]) / (pow(bssn_conformal_fact, 1.5) * conformal_fact);
      conformal_fact_der2[i][j] += (0.75 * bssn_conformal_fact_der[i] * bssn_conformal_fact_der[j]) / (pow(bssn_conformal_fact, 2.5) * conformal_fact);
    }
  }

  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = fluid_old[0]; conformal_spatial_metric[0][1] = fluid_old[1]; conformal_spatial_metric[0][2] = fluid_old[2];
  conformal_spatial_metric[1][0] = fluid_old[3]; conformal_spatial_metric[1][1] = fluid_old[4]; conformal_spatial_metric[1][2] = fluid_old[5];
  conformal_spatial_metric[2][0] = fluid_old[6]; conformal_spatial_metric[2][1] = fluid_old[7]; conformal_spatial_metric[2][2] = fluid_old[8];

  double conformal_lapse = fluid_old[9];

  double conformal_extrinsic_curvature[3][3];
  conformal_extrinsic_curvature[0][0] = fluid_old[10]; conformal_extrinsic_curvature[0][1] = fluid_old[11]; conformal_extrinsic_curvature[0][2] = fluid_old[12];
  conformal_extrinsic_curvature[1][0] = fluid_old[13]; conformal_extrinsic_curvature[1][1] = fluid_old[14]; conformal_extrinsic_curvature[1][2] = fluid_old[15];
  conformal_extrinsic_curvature[2][0] = fluid_old[16]; conformal_extrinsic_curvature[2][1] = fluid_old[17]; conformal_extrinsic_curvature[2][2] = fluid_old[18];

  double conformal_spatial_metric_der[3][3][3];
  conformal_spatial_metric_der[0][0][0] = fluid_old[19]; conformal_spatial_metric_der[0][0][1] = fluid_old[20]; conformal_spatial_metric_der[0][0][2] = fluid_old[21];
  conformal_spatial_metric_der[0][1][0] = fluid_old[22]; conformal_spatial_metric_der[0][1][1] = fluid_old[23]; conformal_spatial_metric_der[0][1][2] = fluid_old[24];
  conformal_spatial_metric_der[0][2][0] = fluid_old[25]; conformal_spatial_metric_der[0][2][1] = fluid_old[26]; conformal_spatial_metric_der[0][2][2] = fluid_old[27];

  conformal_spatial_metric_der[1][0][0] = fluid_old[28]; conformal_spatial_metric_der[1][0][1] = fluid_old[29]; conformal_spatial_metric_der[1][0][2] = fluid_old[30];
  conformal_spatial_metric_der[1][1][0] = fluid_old[31]; conformal_spatial_metric_der[1][1][1] = fluid_old[32]; conformal_spatial_metric_der[1][1][2] = fluid_old[33];
  conformal_spatial_metric_der[1][2][0] = fluid_old[34]; conformal_spatial_metric_der[1][2][1] = fluid_old[35]; conformal_spatial_metric_der[1][2][2] = fluid_old[36];

  conformal_spatial_metric_der[2][0][0] = fluid_old[37]; conformal_spatial_metric_der[2][0][1] = fluid_old[38]; conformal_spatial_metric_der[2][0][2] = fluid_old[39];
  conformal_spatial_metric_der[2][1][0] = fluid_old[40]; conformal_spatial_metric_der[2][1][1] = fluid_old[41]; conformal_spatial_metric_der[2][1][2] = fluid_old[42];
  conformal_spatial_metric_der[2][2][0] = fluid_old[43]; conformal_spatial_metric_der[2][2][1] = fluid_old[44]; conformal_spatial_metric_der[2][2][2] = fluid_old[45];

  double conformal_lapse_der[3];
  conformal_lapse_der[0] = fluid_old[46];
  conformal_lapse_der[1] = fluid_old[47];
  conformal_lapse_der[2] = fluid_old[48];

  double conformal_aux_vect[3];
  conformal_aux_vect[0] = fluid_old[49];
  conformal_aux_vect[1] = fluid_old[50];
  conformal_aux_vect[2] = fluid_old[51];

  double conformal_shift_vect[3];
  conformal_shift_vect[0] = fluid_old[52];
  conformal_shift_vect[1] = fluid_old[53];
  conformal_shift_vect[2] = fluid_old[54];

  double conformal_shift_vect_der[3][3];
  conformal_shift_vect_der[0][0] = fluid_old[55]; conformal_shift_vect_der[0][1] = fluid_old[56]; conformal_shift_vect_der[0][2] = fluid_old[57];
  conformal_shift_vect_der[1][0] = fluid_old[58]; conformal_shift_vect_der[1][1] = fluid_old[59]; conformal_shift_vect_der[1][2] = fluid_old[60];
  conformal_shift_vect_der[2][0] = fluid_old[61]; conformal_shift_vect_der[2][1] = fluid_old[62]; conformal_shift_vect_der[2][2] = fluid_old[63];

  bool in_excision_region = false;
  if (conformal_lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double inv_conformal_spatial_metric[3][3];
    double conformal_spatial_det =
      (conformal_spatial_metric[0][0] * ((conformal_spatial_metric[1][1] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[2][1] * conformal_spatial_metric[1][2]))) -
      (conformal_spatial_metric[0][1] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[1][2] * conformal_spatial_metric[2][0]))) +
      (conformal_spatial_metric[0][2] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][1]) - (conformal_spatial_metric[1][1] * conformal_spatial_metric[2][0])));
    
    double trace = 0.0;
    for (int i = 0; i < 3; i++) {
      trace += conformal_spatial_metric[i][i];
    }

    double conformal_spatial_metric_sq[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_spatial_metric_sq[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_sq[i][j] += conformal_spatial_metric[i][k] * conformal_spatial_metric[k][j];
        }
      }
    }

    double sq_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      sq_trace += conformal_spatial_metric_sq[i][i];
    }

    double euclidean_metric[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i == j) {
          euclidean_metric[i][j] = 1.0;
        }
        else {
          euclidean_metric[i][j] = 0.0;
        }
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_conformal_spatial_metric[i][j] = (1.0 / conformal_spatial_det) *
          ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * conformal_spatial_metric[i][j]) + conformal_spatial_metric_sq[i][j]);
      }
    }

    double evolution_func = 0.0;
    if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
      evolution_func = 0.0;
    }
    else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
      evolution_func = 1.0;
    }

    double conformal_extrinsic_curvature_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
      }
    }

    double conformal_extrinsic_curvature_mixed[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_mixed[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          conformal_extrinsic_curvature_mixed[i][j] += inv_conformal_spatial_metric[l][j] * conformal_extrinsic_curvature[i][l];
        }
      }
    }

    double conformal_extrinsic_curvature_raised[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_raised[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_extrinsic_curvature_raised[i][j] += inv_conformal_spatial_metric[i][l] * inv_conformal_spatial_metric[m][j] * conformal_extrinsic_curvature[l][m];
          }
        }
      }
    }

    double conformal_shift_vect_der_lowered[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_shift_vect_der_lowered[i][j] = 0.0;

        for (int k = 0; k < 3; k++) {
          conformal_shift_vect_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_vect_der[i][k];
        }
      }
    }

    double conformal_shift_vect_der_switched[3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        conformal_shift_vect_der_switched[i][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_shift_vect_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_vect_der[l][m];
          }
        }
      }
    }

    double conformal_shift_vect_der_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      conformal_shift_vect_der_trace += conformal_shift_vect_der[i][i];
    }

    double conformal_symmetrized_shift[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_symmetrized_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_vect_der_lowered[i][j] + conformal_shift_vect_der_lowered[j][i]);
      }
    }

    double slicing_func = 0.0;
    if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
      slicing_func = 0.0;
    }
    else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
      slicing_func = conformal_extrinsic_curvature_trace / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
    }
    else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
      slicing_func = 2.0 * conformal_extrinsic_curvature_trace / (conformal_lapse * (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
    }

    double conformal_spatial_metric_der_raised1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised1[k][i][j] = 0.0;
          
          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised1[k][i][j] += inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j];
          }
        }
      }
    }

    double conformal_spatial_metric_der_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised3[i][j][k] += inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l];
          }
        }
      }
    }

    double conformal_spatial_metric_der_lowered1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_lowered1[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_spatial_metric_der_lowered1[i][j][k] += inv_conformal_spatial_metric[j][l] * inv_conformal_spatial_metric[m][k] * conformal_spatial_metric_der[i][l][m];
            }
          }
        }
      }
    }

    double conformal_spatial_metric_der_lowered3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_lowered3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_spatial_metric_der_lowered3[i][j][k] += inv_conformal_spatial_metric[i][l] * inv_conformal_spatial_metric[m][j] * conformal_spatial_metric_der[l][m][k];
            }
          }
        }
      }
    }

    double conformal_aux_vect_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_aux_vect_raised[k] = 0.0;
        
      for (int l = 0; l < 3; l++) {
        conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
      }
    }

    double conformal_spatial_christoffel[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          conformal_spatial_christoffel[i][k][l] = 0.0;

          for (int m = 0; m < 3; m++) {
            conformal_spatial_christoffel[i][k][l] += inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[l][m][k];
            conformal_spatial_christoffel[i][k][l] += inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[k][m][l];
            conformal_spatial_christoffel[i][k][l] -= inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[m][k][l];
          }
        }
      }
    }

    double conformal_lapse_der_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_lapse_der_raised[k] = 0.0;

      for (int l = 0; l < 3; l++) {
        conformal_lapse_der_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_lapse_der[l];
      }
    }

    double conformal_fact_der_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_fact_der_raised[k] = 0.0;

      for (int l = 0; l < 3; l++) {
        conformal_fact_der_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_fact_der[l];
      }
    }

    double Y_tensor[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Y_tensor[i][j] = 2.0 * conformal_fact_der2[i][j];

        for (int k = 0; k < 3; k++) {
          Y_tensor[i][j] += 2.0 * conformal_spatial_metric[i][j] * conformal_fact_der_raised[k] * conformal_fact_der[k];
        }

        for (int r = 0; r < 3; r++) {
          Y_tensor[i][j] -= 2.0 * conformal_fact_der[r] * conformal_spatial_christoffel[r][i][j];
        }
        Y_tensor[i][j] -= 6.0 * conformal_fact_der[i] * conformal_fact_der[j];
      }
    }

    double Y_tensor_mixed[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Y_tensor_mixed[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          Y_tensor_mixed[i][j] += inv_conformal_spatial_metric[l][j] * Y_tensor[i][l];
        }
      }
    }

    double conformal_spatial_metric_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_spatial_metric_source[i][j] = -2.0 * conformal_lapse * (conformal_extrinsic_curvature[i][j] / (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
        conformal_spatial_metric_source[i][j] += 2.0 * conformal_lapse * conformal_symmetrized_shift[i][j];

        for (int r = 0; r < 3; r++) {
          conformal_spatial_metric_source[i][j] += 2.0 * conformal_shift_vect[r] * conformal_spatial_metric_der[r][i][j];
          conformal_spatial_metric_source[i][j] += 4.0 * conformal_shift_vect[r] * conformal_fact_der[r] * conformal_spatial_metric[i][j];
        }
      }
    }

    double conformal_lapse_source = -(conformal_lapse * conformal_lapse) * slicing_func;
    for (int r = 0; r < 3; r++) {
      conformal_lapse_source += conformal_lapse * conformal_shift_vect[r] * conformal_lapse_der[r];
    }

    double conformal_extrinsic_curvature_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_source[i][j] = 0.0;

        for (int r = 0; r < 3; r++) {
          conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_extrinsic_curvature[i][r] * conformal_shift_vect_der[j][r];
          conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_extrinsic_curvature[j][r] * conformal_shift_vect_der[i][r];
          conformal_extrinsic_curvature_source[i][j] -= 2.0 * conformal_extrinsic_curvature[i][j] * conformal_shift_vect_der[r][r];
        }

        for (int k = 0; k < 3; k++) {
          conformal_extrinsic_curvature_source[i][j] -= (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[i][k] * conformal_extrinsic_curvature[k][j]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_extrinsic_curvature_source[i][j] += (conformal_lapse * conformal_extrinsic_curvature_trace * conformal_extrinsic_curvature[i][j]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * conformal_spatial_christoffel[k][r][i] * conformal_spatial_christoffel[r][k][j];

            conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse * conformal_spatial_metric_der_raised3[i][k][r] * conformal_spatial_metric_der_raised3[r][j][k];
            conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse * conformal_spatial_metric_der_raised3[j][k][r] * conformal_spatial_metric_der_raised3[r][i][k];
            conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_spatial_christoffel[k][k][r] * conformal_spatial_christoffel[r][i][j];

            conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * (2.0 * conformal_spatial_metric_der_raised3[k][r][k] - conformal_lapse_der[r]) *
              (conformal_spatial_metric_der_raised3[i][j][r] + conformal_spatial_metric_der_raised3[j][i][r]);
          }

          conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_lapse_der[i] * (conformal_aux_vect[j] - (0.5 * conformal_spatial_metric_der_raised3[j][k][k]));
          conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_lapse_der[j] * (conformal_aux_vect[i] - (0.5 * conformal_spatial_metric_der_raised3[i][k][k]));

          conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * evolution_func * conformal_aux_vect_raised[k] * conformal_spatial_metric_der[k][i][j];
        }

        for (int k = 0; k < 3; k++) {
          for (int r = 0; r < 3; r++) {
            for (int s = 0; s < 3; s++) {
              conformal_extrinsic_curvature_source[i][j] -= (0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
                conformal_spatial_metric_der_lowered1[k][r][s] * conformal_spatial_christoffel[k][r][s];
              conformal_extrinsic_curvature_source[i][j] += (0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
                conformal_spatial_metric_der_raised3[k][r][r] * conformal_spatial_metric_der_lowered3[k][s][s];
            }
          }

          conformal_extrinsic_curvature_source[i][j] -= (0.5 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) * conformal_aux_vect_raised[k] * conformal_lapse_der[k];
        }

        for (int r = 0; r < 3; r++) {
          for (int s = 0; s < 3; s++) {
            conformal_extrinsic_curvature_source[i][j] += ((0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
              conformal_extrinsic_curvature_raised[r][s] * conformal_extrinsic_curvature[r][s]) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          }
        }

        conformal_extrinsic_curvature_source[i][j] -= ((0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
          (conformal_extrinsic_curvature_trace * conformal_extrinsic_curvature_trace)) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

        conformal_extrinsic_curvature_source[i][j] -= Y_tensor[i][j];
        conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse_der[i] * conformal_fact_der[j];
        conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse_der[j] * conformal_fact_der[i];
        
        for (int k = 0; k < 3; k++) {
          conformal_extrinsic_curvature_source[i][j] += conformal_spatial_metric[i][j] * (evolution_func - 1.0) * Y_tensor_mixed[k][k];
          conformal_extrinsic_curvature_source[i][j] -= 2.0 * conformal_spatial_metric[i][j] * conformal_lapse_der_raised[k] * conformal_fact_der[k];
        }
      }
    }

    double conformal_aux_vect_source[3];
    for (int i = 0; i < 3; i++) {
      conformal_aux_vect_source[i] = 0.0;

      for (int r = 0; r < 3; r++) {
        conformal_aux_vect_source[i] += (conformal_lapse * conformal_lapse_der[r] * conformal_extrinsic_curvature_mixed[i][r]) /
          (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        
        if (i == r) {
          conformal_aux_vect_source[i] -= (conformal_lapse * conformal_lapse_der[r] * conformal_extrinsic_curvature_trace) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        }

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect_source[i] += (conformal_lapse * conformal_extrinsic_curvature_mixed[s][r] * conformal_spatial_metric_der_raised3[i][r][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_aux_vect_source[i] -= (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[s][r] * conformal_spatial_metric_der_raised3[r][i][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          conformal_aux_vect_source[i] -= (conformal_lapse * conformal_extrinsic_curvature_mixed[i][r] * conformal_spatial_metric_der_raised3[r][s][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_aux_vect_source[i] += (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[i][r] * conformal_spatial_metric_der_raised3[s][r][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          conformal_aux_vect_source[i] -= (6.0 * conformal_lapse * conformal_fact_der[r] * conformal_extrinsic_curvature_mixed[s][r]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          if (s == r) {
            conformal_aux_vect_source[i] += (2.0 * conformal_lapse * conformal_fact_der[r] * conformal_extrinsic_curvature_trace) /
              (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          }
        }

        conformal_aux_vect_source[i] += 2.0 * conformal_shift_vect_der[i][r] * conformal_aux_vect[r];

        if (i == r) {
          conformal_aux_vect_source[i] -= 2.0 * conformal_shift_vect_der_trace * conformal_aux_vect[r];
        }

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect_source[i] += 2.0 * conformal_spatial_metric_der_raised3[r][i][s] * conformal_shift_vect_der_switched[r][s];

          if (i == s) {
            for (int j = 0; j < 3; j++) {
              conformal_aux_vect_source[i] -= 2.0 * conformal_spatial_metric_der_raised1[j][j][r] * conformal_shift_vect_der_switched[r][s];
            }
          }
        }

        conformal_aux_vect_source[i] += 4.0 * conformal_shift_vect_der[i][r] * conformal_fact_der[r];
        conformal_aux_vect_source[i] -= 4.0 * conformal_shift_vect_der_trace * conformal_fact_der[i];
      }
    }

    for (int i = 0; i < 77; i++) {
      fluid_new[i] = fluid_old[i];
    }

    fluid_new[0] += dt * conformal_spatial_metric_source[0][0]; fluid_new[1] += dt * conformal_spatial_metric_source[0][1]; fluid_new[2] += dt * conformal_spatial_metric_source[0][2];
    fluid_new[3] += dt * conformal_spatial_metric_source[1][0]; fluid_new[4] += dt * conformal_spatial_metric_source[1][1]; fluid_new[5] += dt * conformal_spatial_metric_source[1][2];
    fluid_new[6] += dt * conformal_spatial_metric_source[2][0]; fluid_new[7] += dt * conformal_spatial_metric_source[2][1]; fluid_new[8] += dt * conformal_spatial_metric_source[2][2];

    fluid_new[9] += dt * conformal_lapse_source;

    fluid_new[10] += dt * conformal_extrinsic_curvature_source[0][0]; fluid_new[11] += dt * conformal_extrinsic_curvature_source[0][1]; fluid_new[12] +=
      dt * conformal_extrinsic_curvature_source[0][2];
    fluid_new[13] += dt * conformal_extrinsic_curvature_source[1][0]; fluid_new[14] += dt * conformal_extrinsic_curvature_source[1][1]; fluid_new[15] +=
      dt * conformal_extrinsic_curvature_source[1][2];
    fluid_new[16] += dt * conformal_extrinsic_curvature_source[2][0]; fluid_new[17] += dt * conformal_extrinsic_curvature_source[2][1]; fluid_new[18] +=
      dt * conformal_extrinsic_curvature_source[2][2];

    fluid_new[49] += dt * conformal_aux_vect_source[0];
    fluid_new[50] += dt * conformal_aux_vect_source[1];
    fluid_new[51] += dt * conformal_aux_vect_source[2];

    if (fluid_new[9] < excision_threshold) {
      for (int i = 0; i < 77; i++) {
        fluid_new[i] = 0.0;
      }
    }
  }
  else {
    for (int i = 0; i < 77; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_vacuum_einstein_conformal_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double excision_threshold = mom_em->vacuum_einstein_conformal_excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = mom_em->vacuum_einstein_conformal_spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = mom_em->vacuum_einstein_conformal_spacetime_evolution;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[77], f_stage1[77], f_stage2[77], f_old[77];

    for (int j = 0; j < 77; j++) {
      f_old[j] = f[j];
    }

    explicit_vacuum_einstein_conformal_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 77; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_vacuum_einstein_conformal_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 77; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_vacuum_einstein_conformal_source_update_euler(mom_em, excision_threshold, spacetime_slicing, spacetime_evolution, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 77; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_gr_mhd_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new)
{
  double lapse = fluid_old[9];
  double shift_x = fluid_old[10];
  double shift_y = fluid_old[11];
  double shift_z = fluid_old[12];

  double spatial_metric[3][3];
  spatial_metric[0][0] = fluid_old[13]; spatial_metric[0][1] = fluid_old[14]; spatial_metric[0][2] = fluid_old[15];
  spatial_metric[1][0] = fluid_old[16]; spatial_metric[1][1] = fluid_old[17]; spatial_metric[1][2] = fluid_old[18];
  spatial_metric[2][0] = fluid_old[19]; spatial_metric[2][1] = fluid_old[20]; spatial_metric[2][2] = fluid_old[21];

  double inv_spatial_metric[3][3];
  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inv_spatial_metric[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = fluid_old[22]; extrinsic_curvature[0][1] = fluid_old[23]; extrinsic_curvature[0][2] = fluid_old[24];
  extrinsic_curvature[1][0] = fluid_old[25]; extrinsic_curvature[1][1] = fluid_old[26]; extrinsic_curvature[1][2] = fluid_old[27];
  extrinsic_curvature[2][0] = fluid_old[28]; extrinsic_curvature[2][1] = fluid_old[29]; extrinsic_curvature[2][2] = fluid_old[30];

  bool in_excision_region = false;
  if (fluid_old[31] < pow(10.0, -8.0)) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double D = fluid_old[0] / sqrt(spatial_det);
    double momx = fluid_old[1] / sqrt(spatial_det);
    double momy = fluid_old[2] / sqrt(spatial_det);
    double momz = fluid_old[3] / sqrt(spatial_det);
    double Etot = fluid_old[4] / sqrt(spatial_det);

    double mag_x = fluid_old[5] / sqrt(spatial_det);
    double mag_y = fluid_old[6] / sqrt(spatial_det);
    double mag_z = fluid_old[7] / sqrt(spatial_det);

    double cov_mom[3];
    cov_mom[0] = momx; cov_mom[1] = momy; cov_mom[2] = momz;

    double mom[3];
    for (int i = 0; i < 3; i++) {
      mom[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        mom[i] += inv_spatial_metric[i][j] * cov_mom[j];
      }
    }

    double M_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      M_sq += (mom[i] * cov_mom[i]);
    }

    double mag[3];
    mag[0] = mag_x; mag[1] = mag_y; mag[2] = mag_z;

    double cov_mag[3];
    for (int i = 0; i < 3; i++) {
      cov_mag[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_mag[i] += spatial_metric[i][j] * mag[j];
      }
    }
    
    double mag_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      mag_sq += mag[i] * cov_mag[i];
    }

    double tau_star = 0.0;
    for (int i = 0; i < 3; i++) {
      tau_star += (mag[i] * cov_mom[i]);
    }

    double p_guess = 0.0, p_new = 0.0, rho = 0.0, z = 0.0, W = 0.0;
    int iter = 0;

    while (iter < 100) {
      double a = Etot + D + p_guess + (0.5 * mag_sq);
      double d = 0.5 * ((M_sq * mag_sq) - (tau_star * tau_star));

      double phi = acos((1.0 / a) * sqrt((27.0 * d) / (4.0 * a)));
      double epsilon1 = ((1.0 / 3.0) * a) - ((2.0 / 3.0) * a * cos(((2.0 / 3.0) * phi) + ((2.0 / 3.0) * M_PI)));
      z = epsilon1 - mag_sq;

      double v_sq = ((M_sq * (z * z)) + ((tau_star * tau_star) * (mag_sq + (2.0 * z)))) / ((z * z) * (mag_sq + z) * (mag_sq + z));
      W = 1.0 / sqrt(1.0 - v_sq);
      rho = D / W;
      double h = z / (W * W * rho);

      p_new = rho * (h - 1.0) * ((gas_gamma - 1.0) / gas_gamma);

      if (fabs(p_guess - p_new) < pow(10.0, -15.0)) {
        iter = 100;
      }
      else {
        iter += 1;
        p_guess = p_new;
      }
    }

    if (p_new < pow(10.0, -8.0)) {
      p_new = pow(10.0, -8.0);
    }
    if (rho < pow(10.0, -8.0)) {
      rho = pow(10.0, -8.0);
    }

    double vel[3];
    for (int i = 0; i < 3; i++) {
      vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        vel[i] += (inv_spatial_metric[i][j] * cov_mom[j]) / (z + mag_sq);
        vel[i] += ((mag[j] * cov_mom[j]) * mag[i]) / (z * (z + mag_sq));
      }
    }

    double cov_vel[3];
    for (int i = 0; i < 3; i++) {
      cov_vel[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_vel[i] += spatial_metric[i][j] * vel[j];
      }
    }

    double b0 = 0.0;
    for (int i = 0; i < 3; i++) {
      b0 += W * mag[i] * (cov_vel[i] / lapse);
    }

    double shift[3];
    shift[0] = shift_x; shift[1] = shift_y; shift[2] = shift_z;

    double spacetime_vel[4];
    spacetime_vel[0] = W / lapse;
    for (int i = 0; i < 3; i++) {
      spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse));
    }

    double b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W;
    }

    double b_sq = 0.0;
    for (int i = 0; i < 3; i++) {
      b_sq += (mag[i] * cov_mag[i]) / (W * W);
    }
    b_sq += ((lapse * lapse) * (b0 * b0)) / (W * W);

    double cov_b[3];
    for (int i = 0; i < 3; i++) {
      cov_b[i] = 0.0;

      for (int j = 0; j < 3; j++) {
        cov_b[i] += spatial_metric[i][j] * b[j];
      }
    }

    double h_star = 1.0 + ((p_new / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho);
    double p_star = p_new + (0.5 * b_sq);

    double inv_spacetime_metric[4][4];
    inv_spacetime_metric[0][0] = - (1.0 / (lapse * lapse));
    for (int i = 0; i < 3; i++) {
      inv_spacetime_metric[0][i] = (1.0 / (lapse * lapse)) * shift[i];
      inv_spacetime_metric[i][0] = (1.0 / (lapse * lapse)) * shift[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        inv_spacetime_metric[i][j] = inv_spatial_metric[i][j] - ((1.0 / (lapse * lapse)) * shift[i] * shift[j]);
      }
    }

    double spacetime_b[4];
    spacetime_b[0] = b0;
    for (int i = 0; i < 3; i++) {
      spacetime_b[i + 1] = b[i];
    }

    double stress_energy[4][4];
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        stress_energy[i][j] = (rho * h_star * spacetime_vel[i] * spacetime_vel[j]) + (p_star * inv_spacetime_metric[i][j]) - (spacetime_b[i] * spacetime_b[j]);
      }
    }

    for (int i = 0; i < 3; i++) {
      mom[i] = (rho * h_star * (W * W) * cov_vel[i]) - (lapse * b0 * cov_b[i]);
    }

    double lapse_der[3];
    lapse_der[0] = fluid_old[32];
    lapse_der[1] = fluid_old[33];
    lapse_der[2] = fluid_old[34];

    double shift_der[3][3];
    shift_der[0][0] = fluid_old[35]; shift_der[0][1] = fluid_old[36]; shift_der[0][2] = fluid_old[37];
    shift_der[1][0] = fluid_old[38]; shift_der[1][1] = fluid_old[39]; shift_der[1][2] = fluid_old[40];
    shift_der[2][0] = fluid_old[41]; shift_der[2][1] = fluid_old[42]; shift_der[2][2] = fluid_old[43];

    double spatial_metric_der[3][3][3];
    spatial_metric_der[0][0][0] = fluid_old[44]; spatial_metric_der[0][0][1] = fluid_old[45]; spatial_metric_der[0][0][2] = fluid_old[46];
    spatial_metric_der[0][1][0] = fluid_old[47]; spatial_metric_der[0][1][1] = fluid_old[48]; spatial_metric_der[0][1][2] = fluid_old[49];
    spatial_metric_der[0][2][0] = fluid_old[50]; spatial_metric_der[0][2][1] = fluid_old[51]; spatial_metric_der[0][2][2] = fluid_old[52];

    spatial_metric_der[1][0][0] = fluid_old[53]; spatial_metric_der[1][0][1] = fluid_old[54]; spatial_metric_der[1][0][2] = fluid_old[55];
    spatial_metric_der[1][1][0] = fluid_old[56]; spatial_metric_der[1][1][1] = fluid_old[57]; spatial_metric_der[1][1][2] = fluid_old[58];
    spatial_metric_der[1][2][0] = fluid_old[59]; spatial_metric_der[1][2][1] = fluid_old[60]; spatial_metric_der[1][2][2] = fluid_old[61];

    spatial_metric_der[2][0][0] = fluid_old[62]; spatial_metric_der[2][0][1] = fluid_old[63]; spatial_metric_der[2][0][2] = fluid_old[64];
    spatial_metric_der[2][1][0] = fluid_old[65]; spatial_metric_der[2][1][1] = fluid_old[66]; spatial_metric_der[2][1][2] = fluid_old[67];
    spatial_metric_der[2][2][0] = fluid_old[68]; spatial_metric_der[2][2][1] = fluid_old[69]; spatial_metric_der[2][2][2] = fluid_old[70];

    for (int i = 0; i < 75; i++) {
      fluid_new[i] = fluid_old[i];
    }

    // Energy density source.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fluid_new[4] += dt * (lapse * sqrt(spatial_det)) * (stress_energy[0][0] * shift[i] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (lapse * sqrt(spatial_det)) * (2.0 * stress_energy[0][i + 1] * shift[j] * extrinsic_curvature[i][j]);
        fluid_new[4] += dt * (lapse * sqrt(spatial_det)) * (stress_energy[i + 1][j + 1] * extrinsic_curvature[i][j]);
      }

      fluid_new[4] -= dt * (lapse * sqrt(spatial_det)) * (stress_energy[0][0] * shift[i] * lapse_der[i]);
      fluid_new[4] -= dt * (lapse * sqrt(spatial_det)) * (stress_energy[0][i + 1] * lapse_der[i]);
    }

    // Momentum density sources.
    for (int j = 0; j < 3; j++) {
      fluid_new[1 + j] -= dt * (lapse * sqrt(spatial_det)) * (stress_energy[0][0] * lapse * lapse_der[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          fluid_new[1 + j] += dt * (lapse * sqrt(spatial_det)) * (0.5 * stress_energy[0][0] * shift[k] * shift[l] * spatial_metric_der[j][k][l]);
          fluid_new[1 + j] += dt * (lapse * sqrt(spatial_det)) * (0.5 * stress_energy[k + 1][l + 1] * spatial_metric_der[j][k][l]);
        }

        fluid_new[1 + j] += dt * (lapse * sqrt(spatial_det)) * ((mom[k] / lapse) * shift_der[j][k]);

        for (int i = 0; i < 3; i++) {
          fluid_new[1 + j] += dt * (lapse * sqrt(spatial_det)) * (stress_energy[0][i + 1] * shift[k] * spatial_metric_der[j][i][k]);
        }
      }
    }
  }
  else {
    for (int i = 0; i < 75; i++) {
      fluid_new[i] = fluid_old[i];
    }
  }
}

void
explicit_gr_mhd_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES])
{
  int nfluids = mom_em->nfluids;

  double gas_gamma = mom_em->gr_mhd_gas_gamma;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];

    double f_new[75], f_stage1[75], f_stage2[75], f_old[75];

    for (int j = 0; j < 75; j++) {
      f_old[j] = f[j];
    }

    explicit_gr_mhd_source_update_euler(mom_em, gas_gamma, t_curr, dt, f_old, f_new);
    for (int j = 0; j < 75; j++) {
      f_stage1[j] = f_new[j];
    }

    explicit_gr_mhd_source_update_euler(mom_em, gas_gamma, t_curr + dt, dt, f_stage1, f_new);
    for (int j = 0; j < 75; j++) {
      f_stage2[j] = (0.75 * f_old[j]) + (0.25 * f_new[j]);
    }

    explicit_gr_mhd_source_update_euler(mom_em, gas_gamma, t_curr + (0.5 * dt), dt, f_stage2, f_new);
    for (int j = 0; j < 75; j++) {
      f[j] = ((1.0 / 3.0) * f_old[j]) + ((2.0 / 3.0) * f_new[j]);
    }
  }
}

void
explicit_e_field_source_update_euler(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double e_field_old[3], double* e_field_new,
  double* fluid_s[GKYL_MAX_SPECIES], const double* app_current)
{
  int nfluids = mom_em->nfluids;

  double epsilon0 = mom_em->epsilon0;
  double mu0 = mom_em->mu0;
  double c = 1.0 / sqrt(mu0 * epsilon0);

  e_field_new[0] = e_field_old[0] + (dt * (-(1.0 / epsilon0) * app_current[0]));
  e_field_new[1] = e_field_old[1] + (dt * (-(1.0 / epsilon0) * app_current[1]));
  e_field_new[2] = e_field_old[2] + (dt * (-(1.0 / epsilon0) * app_current[2]));

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double q = mom_em->param[i].charge;
    const double m = mom_em->param[i].mass;

    double rho = f[0];

    if (rho > 0.0) {
      double ux = f[1] / rho;
      double uy = f[2] / rho;
      double uz = f[3] / rho;

      double gamma = sqrt(1.0 + (((ux * ux) + (uy * uy) + (uz * uz)) / (c * c)));

      double vx = ux / gamma;
      double vy = uy / gamma;
      double vz = uz / gamma;

      e_field_new[0] += dt * (-(1.0 / epsilon0) * (q/m) * rho * vx);
      e_field_new[1] += dt * (-(1.0 / epsilon0) * (q/m) * rho * vy);
      e_field_new[2] += dt * (-(1.0 / epsilon0) * (q/m) * rho * vz);
    }
  }
}

void
explicit_e_field_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES], double* em,
  const double* app_current, const double* app_current1, const double* app_current2, const double* ext_em)
{
  double e_field_new[3], e_field_stage1[3], e_field_stage2[3];
  double e_field_old[3];
  e_field_old[0] = em[0]; e_field_old[1] = em[1]; e_field_old[2] = em[2];

  explicit_e_field_source_update_euler(mom_em, t_curr, dt, e_field_old, e_field_new, fluid_s, app_current);
  e_field_stage1[0] = e_field_new[0];
  e_field_stage1[1] = e_field_new[1];
  e_field_stage1[2] = e_field_new[2];

  explicit_e_field_source_update_euler(mom_em, t_curr + dt, dt, e_field_stage1, e_field_new, fluid_s, app_current1);
  e_field_stage2[0] = (0.75 * e_field_old[0]) + (0.25 * e_field_new[0]);
  e_field_stage2[1] = (0.75 * e_field_old[1]) + (0.25 * e_field_new[1]);
  e_field_stage2[2] = (0.75 * e_field_old[2]) + (0.25 * e_field_new[2]);

  explicit_e_field_source_update_euler(mom_em, t_curr + (0.5 * dt), dt, e_field_stage2, e_field_new, fluid_s, app_current2);
  em[0] = ((1.0 / 3.0) * e_field_old[0]) + ((2.0 / 3.0) * e_field_new[0]);
  em[1] = ((1.0 / 3.0) * e_field_old[1]) + ((2.0 / 3.0) * e_field_new[1]);
  em[2] = ((1.0 / 3.0) * e_field_old[2]) + ((2.0 / 3.0) * e_field_new[2]);
}

void
explicit_higuera_cary_push(double* vel, const double q, const double m, const double dt, const double c, const double e_field[3],
  const double b_field[3])
{
  const double q_over_m = (0.5 * dt) * (q / m);
  const double Ex = q_over_m * e_field[0];
  const double Ey = q_over_m * e_field[1];
  const double Ez = q_over_m * e_field[2];
  const double Bx = q_over_m * b_field[0];
  const double By = q_over_m * b_field[1];
  const double Bz = q_over_m * b_field[2];

  const double vel_x_minus = vel[0] + Ex;
  const double vel_y_minus = vel[1] + Ey;
  const double vel_z_minus = vel[2] + Ez;

  const double vel_star = (vel_x_minus * (Bx / c)) + (vel_y_minus * (By / c)) + (vel_z_minus * (Bz / c));
  const double gamma_minus = sqrt(1.0 + (((vel_x_minus * vel_x_minus) + (vel_y_minus * vel_y_minus) + (vel_z_minus * vel_z_minus)) / (c * c)));
  const double dot_tau_tau = (Bx * Bx) + (By * By) + (Bz * Bz);
  const double sigma = (gamma_minus * gamma_minus) - dot_tau_tau;
  const double gamma_new = sqrt(0.5 * (sigma + sqrt((sigma * sigma) + (4.0 * (dot_tau_tau + (vel_star * vel_star))))));

  const double tx = Bx / gamma_new;
  const double ty = By / gamma_new;
  const double tz = Bz / gamma_new;
  const double s =  1.0 / (1.0 + ((tx * tx) + (ty * ty) + (tz * tz)));

  const double t_vel_minus = (vel_x_minus * tx) + (vel_y_minus * ty) + (vel_z_minus * tz);
  const double vel_x_plus = s * (vel_x_minus + (t_vel_minus * tx) + ((vel_y_minus * tz) - (vel_z_minus * ty)));
  const double vel_y_plus = s * (vel_y_minus + (t_vel_minus * ty) + ((vel_z_minus * tx) - (vel_x_minus * tz)));
  const double vel_z_plus = s * (vel_z_minus + (t_vel_minus * tz) + ((vel_x_minus * ty) - (vel_y_minus * tx)));

  vel[0] = vel_x_plus + Ex + ((vel_y_plus * tz) - (vel_z_plus * ty));
  vel[1] = vel_y_plus + Ey + ((vel_z_plus * tx) - (vel_x_plus * tz));
  vel[2] = vel_z_plus + Ez + ((vel_x_plus * ty) - (vel_y_plus * tx));
}

void
explicit_higuera_cary_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* ext_em)
{
  int nfluids = mom_em->nfluids;

  double epsilon0 = mom_em->epsilon0;
  double mu0 = mom_em->mu0;
  double c = 1.0 / sqrt(mu0 * epsilon0);

  double Ex = em[0] + ext_em[0];
  double Ey = em[1] + ext_em[1];
  double Ez = em[2] + ext_em[2];
  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];

  double e_field[3], b_field[3];
  e_field[0] = Ex; e_field[1] = Ey; e_field[2] = Ez;
  b_field[0] = Bx; b_field[1] = By; b_field[2] = Bz;

  for (int i = 0; i < nfluids; i++) {
    double *f = fluid_s[i];
    const double q = mom_em->param[i].charge;
    const double m = mom_em->param[i].mass;

    double rho = f[0];

    if (rho > 0.0) {
      double vel[3];
      vel[0] = f[1] / rho;
      vel[1] = f[2] / rho;
      vel[2] = f[3] / rho;

      explicit_higuera_cary_push(vel, q, m, dt, c, e_field, b_field);

      f[1] = rho * vel[0];
      f[2] = rho * vel[1];
      f[3] = rho * vel[2];
    }
  }
}

void
explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang)
{
  if (nstrang == 0) {
    explicit_e_field_source_update(mom_em, t_curr, dt, fluid_s, em, app_current, app_current1, app_current2, ext_em);
  }
  else if (nstrang == 1) {
    explicit_higuera_cary_update(mom_em, t_curr, dt, fluid_s, app_accel_s, em, ext_em);
  }
}
