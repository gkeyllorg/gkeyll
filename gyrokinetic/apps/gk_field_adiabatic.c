#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <math.h>
#include <time.h>

// Adiabatic species response in the field equation. Methods valid for any adiabatic species
// are named gk_field_adiabatic_*, those specific to the adiabatic electron model
// gk_field_adiab_elc_*.
//
// Adiabatic electron response with a flux-surface average (FSA) in 2x/3x:
//   -div(eps grad_perp phi) + K (phi - <phi>) = rho,  K = (e^2 n0/Te) J,
// with <phi> = int J phi dy dz / int J dy dz. Writing H = L + K, R for the
// FSA reduction and E for the extension of a function of x to a field constant in
// (y,z), the Woodbury identity gives
//   phi = phi1 + H^{-1} K E psi,  H phi1 = rho,  (I - G) psi = R phi1,
// where G = R H^{-1} K E is an m x m matrix (m = Nx*num_basis_x) built
// at init with Nx*num_basis_x 1D Helmholtz solves.
//
// With Pade FLR the solve is for Phi_0 (phi = A Phi_0, A = 1 - rho^2*nabla_perp^2) and the
// response is M (Phi_0 - <Phi_0>), M = K - div(K_0 rho^2 J g^ij grad), so H = L + M with
// L = -div(eps_pol J g^ij grad). Since M = a H + (1-a) K, a = K_0 rho^2/(eps_pol + K_0 rho^2),
// the Woodbury identity with M E in place of K E gives
//   Phi_0 = phi1 + H^{-1} K E psi + (K_0 rho^2/eps_pol) (E psi - phi_rhs0),
// with the same G and psi: the zonal Phi_0 then obeys L Phi_0 = rho, without K_0 rho^2.

// Copy an array w(x,p) into a contiguous vector v(m) with m=1,...,Nx*num_basis_x.
static void
gk_field_adiab_elc_zonal_pack(
  const struct gk_field_adiabatic *ad, const struct gkyl_array *arr_ho, double *vec
)
{
  int nb = ad->basis_x.num_basis;
  int k = 0;
  for (int i = ad->local_x.lower[0]; i <= ad->local_x.upper[0]; i++) {
    const double *c = gkyl_array_cfetch(arr_ho, gkyl_range_idx(&ad->local_x, &i));
    for (int b = 0; b < nb; b++) {
      vec[k++] = c[b];
    }
  }
}

// Copy a v(m) vector into an array w(x,p).
static void
gk_field_adiab_elc_zonal_unpack(
  const struct gk_field_adiabatic *ad, const double *vec, struct gkyl_array *arr_ho
)
{
  int nb = ad->basis_x.num_basis;
  int k = 0;
  for (int i = ad->local_x.lower[0]; i <= ad->local_x.upper[0]; i++) {
    double *c = gkyl_array_fetch(arr_ho, gkyl_range_idx(&ad->local_x, &i));
    for (int b = 0; b < nb; b++) {
      c[b] = vec[k++];
    }
  }
}

// Compute <phi> = int J phi dy dz / int J dy dz into adiab.psi. This is the R operator.
void
gk_field_adiab_elc_fsa(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi
)
{
  const struct gk_field_adiabatic *ad = &gkf->adiab;

  gkyl_dg_mul_op_range(
    &app->basis, 0, ad->jphi, 0, app->gk_geom->geo_int.jacobgeo, 0, phi, &app->local
  );
  gkyl_array_average_advance(ad->fs_avg, ad->jphi, ad->avg_jphi);

  // The average divides by the global (y,z) extent, so the sum over ranks is the global average.
  gkyl_comm_allreduce(
    app->comm, GKYL_DOUBLE, GKYL_SUM, ad->avg_jphi->ncomp * ad->avg_jphi->size, ad->avg_jphi->data,
    ad->avg_jphi_red->data
  );

  gkyl_dg_div_op_range(
    ad->div_mem, &ad->basis_x, 0, ad->psi, 0, ad->avg_jphi_red, 0, ad->avg_j, &ad->local_x
  );
}

// Extend a 1D function of x to a field constant along the other directions. This is the E operator.
void
gk_field_adiab_elc_inflate(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *psi,
  struct gkyl_array *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->adiab;
  if (app->cdim == 2) {
    gkyl_translate_dim_advance(ad->inflate_lo, &ad->local_x, &app->local, psi, 1, out);
  } else {
    gkyl_translate_dim_advance(ad->inflate_lo, &ad->local_x, &ad->local_xy, psi, 1, ad->tmp_xy);
    gkyl_translate_dim_advance(ad->inflate_up, &ad->local_xy, &app->local, ad->tmp_xy, 1, out);
  }
}

// Solve H phi = rhs, smoothing the rhs along z like the charge density.
static void
gk_field_adiabatic_helmholtz_solve(
  gkyl_gyrokinetic_app *app, struct gk_field *gkf, struct gkyl_array *rhs, struct gkyl_array *phi
)
{
  gkf->fem_projection_par_rho_func(app, gkf, rhs, rhs);
  gkyl_fem_poisson_perp_set_rhs(gkf->fem_poisson_perp, rhs);
  gkyl_fem_poisson_perp_solve(gkf->fem_poisson_perp, phi);
}

// Solve H phi = K E psi with zero Dirichlet values and subtract the rhs = 0 solution.
static void
gk_field_adiab_elc_response_solve(
  gkyl_gyrokinetic_app *app, struct gk_field *gkf, const struct gkyl_array *psi,
  struct gkyl_array *phi
)
{
  struct gk_field_adiabatic *ad = &gkf->adiab;
  gk_field_adiab_elc_inflate(app, gkf, psi, phi); // E psi
  gkyl_dg_mul_op_range(&app->basis, 0, ad->rhs2, 0, ad->kJ, 0, phi, &app->local); // K E psi
  gk_field_adiabatic_helmholtz_solve(app, gkf, ad->rhs2, phi); // get phi = H^{-1} K E psi
  gkyl_array_accumulate_range(phi, -1.0, ad->phi_rhs0, &app->local);
}

// Solve (I - G) psi = <phi1> on the host, with <phi1> in psi on input.
static void
gk_field_adiab_elc_zonal_solve(gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->adiab;
  if (app->use_gpu) {
    gkyl_array_copy(ad->psi_ho, ad->psi);
  }
  gk_field_adiab_elc_zonal_pack(ad, ad->psi_ho, gkyl_mat_get_col(ad->rhs_m, 0));
  gkyl_mat_copy(ad->A_lu, ad->A);
  bool status = gkyl_mat_linsolve_lu(ad->A_lu, ad->rhs_m, gkyl_mem_buff_data(ad->ipiv));
  assert(status);
  gk_field_adiab_elc_zonal_unpack(ad, gkyl_mat_get_ccol(ad->rhs_m, 0), ad->psi_ho);
  if (app->use_gpu) {
    gkyl_array_copy(ad->psi, ad->psi_ho);
  }
}

void
gk_field_adiab_elc_rhs_phi_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  struct gk_field_adiabatic *ad = &field->adiab;

  // Solve H phi1 = rho.
  gk_field_adiabatic_helmholtz_solve(app, field, field->rho_c, field->phi_smooth);

  // Zonal correction: (I - G) psi = <phi1>, phi2 = H^{-1} K E psi.
  gk_field_adiab_elc_fsa(app, field, field->phi_smooth); // Compute <phi1> into ad->psi.
  gk_field_adiab_elc_zonal_solve(app, field); // Solve (I - G) psi = <phi1>.
  gk_field_adiab_elc_response_solve(app, field, ad->psi, ad->phi2); // Solve H phi2 = K E psi.
  gkyl_array_accumulate_range(field->phi_smooth, 1.0, ad->phi2, &app->local); // phi = phi1 + phi2.
  if (ad->flr_zonal_fac > 0.0) { // This is not good programming because we introduce a logic in the time loop. We should consider function pointers to avoid this.
    // Pade FLR: add (K_0 rho^2/eps_pol) (E psi - phi_rhs0).
    gk_field_adiab_elc_inflate(app, field, ad->psi, ad->rhs2);
    gkyl_array_accumulate_range(ad->rhs2, -1.0, ad->phi_rhs0, &app->local);
    gkyl_array_accumulate_range(field->phi_smooth, ad->flr_zonal_fac, ad->rhs2, &app->local);
  }

  // Smooth the potential along z.
  field->fem_projection_par_phi_func(app, field, field->phi_smooth, field->phi_smooth);

  field->invert_flr(app, field, field->phi_smooth);
}

// out = phi - <phi>.
void
gk_field_adiab_elc_dphi(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  gk_field_adiab_elc_fsa(app, gkf, phi);
  gk_field_adiab_elc_inflate(app, gkf, gkf->adiab.psi, out);
  gkyl_array_scale_range(out, -1.0, &app->local);
  gkyl_array_accumulate_range(out, 1.0, phi, &app->local);
}

// Add factor*(1/2) int K (phi - <phi>)^2 over the local range to out.
void
gk_field_adiab_elc_energy_accumulate(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, double factor, double *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->adiab;

  gk_field_adiab_elc_dphi(app, gkf, gkf->phi_smooth, ad->phi2);
  gkyl_array_integrate_advance(
    ad->calc_energy, ad->phi2, 0.5 * factor, ad->kJ, &app->local, &app->local, ad->energy_red
  );

  if (app->use_gpu) {
    double en[2];
    gkyl_cu_memcpy(&en[0], out, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(&en[1], ad->energy_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    en[0] += en[1];
    gkyl_cu_memcpy(out, &en[0], sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    out[0] += ad->energy_red[0];
  }
}

// Background charge density q_s*n_s0*J and K = (q_s^2 n_s0/T_s)*J from the t=0 moments
// of the adiabatic species.
void
gk_field_adiabatic_profiles_calc(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->adiab;
  struct gk_species *gks = &app->species[ad->species_idx];
  double q_s = gks->info.charge;

  // Maxwellian moments n_s0*J, upar and T_s/m_s.
  struct gk_species_moment *moms = &gks->lte.moms;
  gk_species_moment_calc(moms, gks->local, app->local, gks->f);

  gkyl_array_set_offset(ad->rho_bg, q_s, moms->marr, 0);
  gkyl_dg_div_op_range(
    moms->mem_geo, &app->basis, 0, ad->kJ, 0, moms->marr, 2, moms->marr, &app->local
  );
  gkyl_array_scale(ad->kJ, q_s * q_s / gks->info.mass);
}

void
gk_field_adiab_elc_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->adiab;

  assert(app->poly_order == 1); // array_average and translate_dim kernels are p=1 only.

  int cdim = app->cdim;
  int nb = app->basis.num_basis;

  ad->phi2 = mkarr(app->use_gpu, nb, app->local_ext.volume);
  ad->rhs2 = mkarr(app->use_gpu, nb, app->local_ext.volume);
  ad->phi_rhs0 = mkarr(app->use_gpu, nb, app->local_ext.volume);

  // 1D (x) basis and ranges.
  gkyl_cart_modal_serendip(&ad->basis_x, 1, app->poly_order);
  gkyl_range_init(&ad->local_x_ext, 1, &app->local_ext.lower[0], &app->local_ext.upper[0]);
  gkyl_sub_range_init(&ad->local_x, &ad->local_x_ext, &app->local.lower[0], &app->local.upper[0]);
  int nb_x = ad->basis_x.num_basis;
  long vol_x = ad->local_x_ext.volume;

  ad->jphi = mkarr(app->use_gpu, nb, app->local_ext.volume);
  ad->avg_jphi = mkarr(app->use_gpu, nb_x, vol_x);
  ad->avg_jphi_red = mkarr(app->use_gpu, nb_x, vol_x);
  ad->avg_j = mkarr(app->use_gpu, nb_x, vol_x);
  ad->psi = mkarr(app->use_gpu, nb_x, vol_x);
  ad->psi_ho = app->use_gpu ? mkarr(false, nb_x, vol_x) : gkyl_array_acquire(ad->psi);
  ad->div_mem = app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(ad->local_x.volume, nb_x) :
                               gkyl_dg_bin_op_mem_new(ad->local_x.volume, nb_x);

  // Flux surface average operator.
  int avg_dim[GKYL_MAX_CDIM] = {0};
  for (int d = 1; d < cdim; d++) {
    avg_dim[d] = 1; // Average all except first dimension (x).
  }
  ad->fs_avg = gkyl_array_average_inew(&(struct gkyl_array_average_inp){
    .grid = &app->grid,
    .basis = app->basis,
    .basis_avg = ad->basis_x,
    .local = &app->local,
    .local_avg = &ad->local_x,
    .local_avg_ext = &ad->local_x_ext,
    .weight = NULL, // We need to divide manually to allow for MPI z decomposition.
    .avg_dim = avg_dim,
    .use_gpu = app->use_gpu,
  });

  // Global (y,z) average of J (time independent).
  gkyl_array_average_advance(ad->fs_avg, app->gk_geom->geo_int.jacobgeo, ad->avg_jphi);
  gkyl_comm_allreduce(
    app->comm, GKYL_DOUBLE, GKYL_SUM, nb_x * vol_x, ad->avg_jphi->data, ad->avg_j->data
  );

  // Updaters extending a 1D function of x to cdim dimensions.
  ad->inflate_up = 0;
  ad->tmp_xy = 0;
  if (cdim == 2) {
    ad->inflate_lo =
      gkyl_translate_dim_new(1, ad->basis_x, 2, app->basis, 0, GKYL_NO_EDGE, app->use_gpu);
  } else {
    struct gkyl_basis basis_xy;
    gkyl_cart_modal_serendip(&basis_xy, 2, app->poly_order);
    gkyl_range_init(&ad->local_xy_ext, 2, app->local_ext.lower, app->local_ext.upper);
    gkyl_sub_range_init(&ad->local_xy, &ad->local_xy_ext, app->local.lower, app->local.upper);
    ad->tmp_xy = mkarr(app->use_gpu, basis_xy.num_basis, ad->local_xy_ext.volume);
    ad->inflate_lo =
      gkyl_translate_dim_new(1, ad->basis_x, 2, basis_xy, 0, GKYL_NO_EDGE, app->use_gpu);
    ad->inflate_up =
      gkyl_translate_dim_new(2, basis_xy, 3, app->basis, 0, GKYL_NO_EDGE, app->use_gpu);
  }

  // Energy diagnostic, (1/2) int K (phi-<phi>)^2.
  ad->calc_energy = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu
  );
  ad->energy_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  // Solution with zero rhs but Dirichlet BCs.
  gkyl_array_clear(ad->rhs2, 0.0);
  gk_field_adiabatic_helmholtz_solve(app, gkf, ad->rhs2, ad->phi_rhs0);

  // Zonal system A = I - G, G = R H^{-1} K E, it requires Nx*num_basis_x 1D Helmholtz solves.
  int m = ad->local_x.volume * nb_x;
  ad->num_zonal = m;
  ad->A = gkyl_mat_new(m, m, 0.0);
  ad->A_lu = gkyl_mat_new(m, m, 0.0);
  ad->rhs_m = gkyl_mat_new(m, 1, 0.0);
  ad->ipiv = gkyl_mem_buff_new(sizeof(long[m]));

  struct timespec wst = gkyl_wall_clock();
  double *ej = gkyl_malloc(sizeof(double[m]));
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < m; i++) {
      ej[i] = i == j ? 1.0 : 0.0;
    }
    gk_field_adiab_elc_zonal_unpack(ad, ej, ad->psi_ho);
    if (app->use_gpu) {
      gkyl_array_copy(ad->psi, ad->psi_ho);
    }
    gk_field_adiab_elc_response_solve(app, gkf, ad->psi, ad->phi2);
    gk_field_adiab_elc_fsa(app, gkf, ad->phi2);
    if (app->use_gpu) {
      gkyl_array_copy(ad->psi_ho, ad->psi);
    }
    gk_field_adiab_elc_zonal_pack(ad, ad->psi_ho, gkyl_mat_get_col(ad->A, j));
  }
  gkyl_free(ej);

  // Spectral radius of G by power iteration (<1 means the fixed-point
  // iteration would converge; the direct solve is exact regardless).
  double *v = gkyl_malloc(sizeof(double[m])), *gv = gkyl_malloc(sizeof(double[m]));
  for (int i = 0; i < m; i++) {
    v[i] = 1.0 + 0.1 * i;
  }
  double grho = 0.0;
  for (int it = 0; it < 200; it++) {
    double vnorm = 0.0;
    for (int i = 0; i < m; i++) {
      vnorm += v[i] * v[i];
    }
    vnorm = sqrt(vnorm);
    for (int i = 0; i < m; i++) {
      v[i] /= vnorm;
    }
    for (int i = 0; i < m; i++) {
      gv[i] = 0.0;
      for (int j = 0; j < m; j++) {
        gv[i] += gkyl_mat_get(ad->A, i, j) * v[j];
      }
    }
    grho = 0.0;
    for (int i = 0; i < m; i++) {
      grho += v[i] * gv[i];
    }
    for (int i = 0; i < m; i++) {
      v[i] = gv[i];
    }
  }
  gkyl_free(v);
  gkyl_free(gv);

  // A = I - G.
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      gkyl_mat_set(ad->A, i, j, (i == j ? 1.0 : 0.0) - gkyl_mat_get(ad->A, i, j));
    }
  }

  gkyl_gyrokinetic_app_cout(
    app, stdout, "Adiabatic species '%s' with FSA: m = %d, rho(G) = %.4e, setup %.3e sec.\n",
    app->species[ad->species_idx].info.name, m, grho, gkyl_time_diff_now_sec(wst)
  );
}

void
gk_field_adiab_elc_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->adiab;

  gkyl_array_release(ad->kSq);
  gkyl_array_release(ad->phi2);
  gkyl_array_release(ad->rhs2);
  gkyl_array_release(ad->phi_rhs0);
  gkyl_array_release(ad->jphi);
  gkyl_array_release(ad->avg_jphi);
  gkyl_array_release(ad->avg_jphi_red);
  gkyl_array_release(ad->avg_j);
  gkyl_array_release(ad->psi);
  gkyl_array_release(ad->psi_ho);
  gkyl_dg_bin_op_mem_release(ad->div_mem);
  gkyl_array_average_release(ad->fs_avg);
  gkyl_translate_dim_release(ad->inflate_lo);
  if (app->cdim == 3) {
    gkyl_translate_dim_release(ad->inflate_up);
    gkyl_array_release(ad->tmp_xy);
  }
  gkyl_array_integrate_release(ad->calc_energy);
  if (app->use_gpu) {
    gkyl_cu_free(ad->energy_red);
  } else {
    gkyl_free(ad->energy_red);
  }
  gkyl_mat_release(ad->A);
  gkyl_mat_release(ad->A_lu);
  gkyl_mat_release(ad->rhs_m);
  gkyl_mem_buff_release(ad->ipiv);
}
