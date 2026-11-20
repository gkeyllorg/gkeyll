#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm_io.h>
#include <gkyl_vlasov_priv.h>

// Neural networks of the moments of PKPM species: training on the diagnostic moments of a run,
// output of the networks, and their test against the moments (a network per PKPM species; the
// arrays of networks are indexed by species).

void
gkyl_vlasov_app_train(
  gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data,
  struct gkyl_kn_vec *output_data
)
{
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_app_train_mom(
        app, i, tm, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms,
        input_data, output_data
      );
    }
  }
}

void
gkyl_vlasov_app_train_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data, struct gkyl_kn_vec *output_data
)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  struct vlasov_species *s = &app->species[sidx];
  // the diagnostic moments need the thermal speed of the current state
  vm_pkpm_species_calc_pkpm_vars(app, s, s->kinetic->f, s->fluid->fluid);
  vm_pkpm_species_moment_diag_calc(
    app, s, s->kinetic->local, app->local, s->kinetic->f, s->fluid->fluid
  );

  if (app->use_gpu) {
    gkyl_array_copy(s->kinetic->pkpm.moms_diag.marr_host, s->kinetic->pkpm.moms_diag.marr);
  }

  int cell_count = 0;
  if (app->cdim == 1) {
    cell_count = app->grid.cells[0];
  } else if (app->cdim == 2) {
    cell_count = app->grid.cells[0] * app->grid.cells[1];
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);

    const double *pkpm_moms_diag_d = gkyl_array_cfetch(s->kinetic->pkpm.moms_diag.marr_host, loc);

    for (int i = 0; i < num_input_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          input_data->vals[count][i * 2] = (float)pkpm_moms_diag_d[input_moms[i] * 2];
          input_data->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 2) + 1];
        } else if (app->cdim == 2) {
          input_data->vals[count][i * 4] = (float)pkpm_moms_diag_d[input_moms[i] * 4];
          input_data->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 1];
          input_data->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 2];
          input_data->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        input_data->vals[count][i * 3] = (float)pkpm_moms_diag_d[input_moms[i] * 3];
        input_data->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 1];
        input_data->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 2];
      }
    }

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          output_data->vals[count][i * 2] = (float)pkpm_moms_diag_d[output_moms[i] * 2];
          output_data->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 2) + 1];
        } else if (app->cdim == 2) {
          output_data->vals[count][i * 4] = (float)pkpm_moms_diag_d[output_moms[i] * 4];
          output_data->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 1];
          output_data->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 2];
          output_data->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        output_data->vals[count][i * 3] = (float)pkpm_moms_diag_d[output_moms[i] * 3];
        output_data->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 1];
        output_data->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 2];
      }
    }

    count += 1;
  }

  struct gkyl_kann_train_params params = {
    .learning_rate = 0.0001f,
    .mini_size = 64,
    .max_epoch = 50,
    .max_drop_streak = 10,
    .frac_val = 0.1f,
  };

  if (gkyl_kann_net_is_cu_dev(ann[sidx])) {
    // Mirror the host training data onto the device and train there.
    struct gkyl_kn_vec *input_data_cu = gkyl_kn_vec_cu_dev_new(input_data->nvec, input_data->N);
    struct gkyl_kn_vec *output_data_cu = gkyl_kn_vec_cu_dev_new(output_data->nvec, output_data->N);
    gkyl_kn_vec_copy(input_data_cu, input_data);
    gkyl_kn_vec_copy(output_data_cu, output_data);

    gkyl_kann_net_train_fnn1(ann[sidx], &params, input_data_cu, output_data_cu);

    gkyl_kn_vec_release(input_data_cu);
    gkyl_kn_vec_release(output_data_cu);
  } else {
    gkyl_kann_net_train_fnn1(ann[sidx], &params, input_data, output_data);
  }
}

void
gkyl_vlasov_app_write_nn(gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann)
{
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_app_write_nn_mom(app, i, tm, frame, ann);
    }
  }
}

void
gkyl_vlasov_app_write_nn_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann
)
{
  struct vlasov_species *s = &app->species[sidx];

  const char *fmt = "%s_moms_nn_%d-%s.dat";
  int sz = gkyl_calc_strlen(fmt, app->name, frame, s->name);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame, s->name);

  gkyl_kann_net_save(ann[sidx], fileNm);
}

void
gkyl_vlasov_app_test(
  gkyl_vlasov_app *app, double tm, int frame, struct gkyl_kann_net **ann, int num_input_moms,
  int *input_moms, int num_output_moms, int *output_moms, struct gkyl_kn_vec *input_data_real,
  struct gkyl_kn_vec *output_data_real, struct gkyl_kn_vec *output_data_predicted
)
{
  for (int i = 0; i < app->num_species; i++) {
    if (app->species[i].type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_app_test_mom(
        app, i, tm, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms,
        input_data_real, output_data_real, output_data_predicted
      );
    }
  }
}

void
gkyl_vlasov_app_test_mom(
  gkyl_vlasov_app *app, int sidx, double tm, int frame, struct gkyl_kann_net **ann,
  int num_input_moms, int *input_moms, int num_output_moms, int *output_moms,
  struct gkyl_kn_vec *input_data_real, struct gkyl_kn_vec *output_data_real,
  struct gkyl_kn_vec *output_data_predicted
)
{
  vlasov_op_split_sync(app); // the staggered split's source-updated state at tcurr
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });

  struct vlasov_species *s = &app->species[sidx];
  // the diagnostic moments need the thermal speed of the current state
  vm_pkpm_species_calc_pkpm_vars(app, s, s->kinetic->f, s->fluid->fluid);
  vm_pkpm_species_moment_diag_calc(
    app, s, s->kinetic->local, app->local, s->kinetic->f, s->fluid->fluid
  );

  if (app->use_gpu) {
    gkyl_array_copy(s->kinetic->pkpm.moms_diag.marr_host, s->kinetic->pkpm.moms_diag.marr);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);

    const double *pkpm_moms_diag_d = gkyl_array_cfetch(s->kinetic->pkpm.moms_diag.marr_host, loc);

    for (int i = 0; i < num_input_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          input_data_real->vals[count][i * 2] = (float)pkpm_moms_diag_d[input_moms[i] * 2];
          input_data_real->vals[count][(i * 2) + 1] =
            (float)pkpm_moms_diag_d[(input_moms[i] * 2) + 1];
        } else if (app->cdim == 2) {
          input_data_real->vals[count][i * 4] = (float)pkpm_moms_diag_d[input_moms[i] * 4];
          input_data_real->vals[count][(i * 4) + 1] =
            (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 1];
          input_data_real->vals[count][(i * 4) + 2] =
            (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 2];
          input_data_real->vals[count][(i * 4) + 3] =
            (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        input_data_real->vals[count][i * 3] = (float)pkpm_moms_diag_d[input_moms[i] * 3];
        input_data_real->vals[count][(i * 3) + 1] =
          (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 1];
        input_data_real->vals[count][(i * 3) + 2] =
          (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 2];
      }
    }

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          output_data_real->vals[count][i * 2] = (float)pkpm_moms_diag_d[output_moms[i] * 2];
          output_data_real->vals[count][(i * 2) + 1] =
            (float)pkpm_moms_diag_d[(output_moms[i] * 2) + 1];
        } else if (app->cdim == 2) {
          output_data_real->vals[count][i * 4] = (float)pkpm_moms_diag_d[output_moms[i] * 4];
          output_data_real->vals[count][(i * 4) + 1] =
            (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 1];
          output_data_real->vals[count][(i * 4) + 2] =
            (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 2];
          output_data_real->vals[count][(i * 4) + 3] =
            (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        output_data_real->vals[count][i * 3] = (float)pkpm_moms_diag_d[output_moms[i] * 3];
        output_data_real->vals[count][(i * 3) + 1] =
          (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 1];
        output_data_real->vals[count][(i * 3) + 2] =
          (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 2];
      }
    }

    count += 1;
  }

  // Batched inference over all sampled cells. On GPU, stage the input and
  // predicted output through device kn_vecs.
  if (gkyl_kann_net_is_cu_dev(ann[sidx])) {
    struct gkyl_kn_vec *input_data_cu =
      gkyl_kn_vec_cu_dev_new(input_data_real->nvec, input_data_real->N);
    struct gkyl_kn_vec *output_data_cu =
      gkyl_kn_vec_cu_dev_new(output_data_predicted->nvec, output_data_predicted->N);
    gkyl_kn_vec_copy(input_data_cu, input_data_real);

    gkyl_kann_net_apply(ann[sidx], input_data_cu, output_data_cu);

    gkyl_kn_vec_copy(output_data_predicted, output_data_cu);
    gkyl_kn_vec_release(input_data_cu);
    gkyl_kn_vec_release(output_data_cu);
  } else {
    gkyl_kann_net_apply(ann[sidx], input_data_real, output_data_predicted);
  }

  struct gkyl_range_iter iter_new;
  gkyl_range_iter_init(&iter_new, &app->local);
  long count_new = 0;
  while (gkyl_range_iter_next(&iter_new)) {
    long loc_new = gkyl_range_idx(&app->local, iter_new.idx);

    double *pkpm_moms_diag_d_new = gkyl_array_fetch(s->kinetic->pkpm.moms_diag.marr_host, loc_new);

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          pkpm_moms_diag_d_new[output_moms[i] * 2] =
            (double)output_data_predicted->vals[count_new][i * 2];
          pkpm_moms_diag_d_new[(output_moms[i] * 2) + 1] =
            (double)output_data_predicted->vals[count_new][(i * 2) + 1];
        } else if (app->cdim == 2) {
          pkpm_moms_diag_d_new[output_moms[i] * 4] =
            (double)output_data_predicted->vals[count_new][i * 4];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 1] =
            (double)output_data_predicted->vals[count_new][(i * 4) + 1];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 2] =
            (double)output_data_predicted->vals[count_new][(i * 4) + 2];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 3] =
            (double)output_data_predicted->vals[count_new][(i * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        pkpm_moms_diag_d_new[output_moms[i] * 3] =
          (double)output_data_predicted->vals[count_new][i * 3];
        pkpm_moms_diag_d_new[(output_moms[i] * 3) + 1] =
          (double)output_data_predicted->vals[count_new][(i * 3) + 1];
        pkpm_moms_diag_d_new[(output_moms[i] * 3) + 2] =
          (double)output_data_predicted->vals[count_new][(i * 3) + 2];
      }
    }

    count_new += 1;
  }

  const char *fmt = "%s-%s_pkpm_moms_predicted_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->name, frame);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->name, frame);

  gkyl_comm_array_write(
    app->comm, &app->grid, &app->local, mt, s->kinetic->pkpm.moms_diag.marr_host, fileNm
  );

  struct gkyl_range_iter iter_old;
  gkyl_range_iter_init(&iter_old, &app->local);
  long count_old = 0;
  while (gkyl_range_iter_next(&iter_old)) {
    long loc_old = gkyl_range_idx(&app->local, iter_old.idx);

    double *pkpm_moms_diag_d_old = gkyl_array_fetch(s->kinetic->pkpm.moms_diag.marr_host, loc_old);

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          pkpm_moms_diag_d_old[output_moms[i] * 2] =
            (double)output_data_real->vals[count_old][i * 2];
          pkpm_moms_diag_d_old[(output_moms[i] * 2) + 1] =
            (double)output_data_real->vals[count_old][(i * 2) + 1];
        } else if (app->cdim == 2) {
          pkpm_moms_diag_d_old[output_moms[i] * 4] =
            (double)output_data_real->vals[count_old][i * 4];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 1] =
            (double)output_data_real->vals[count_old][(i * 4) + 1];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 2] =
            (double)output_data_real->vals[count_old][(i * 4) + 2];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 3] =
            (double)output_data_real->vals[count_old][(i * 4) + 3];
        }
      } else if (app->poly_order == 2) {
        pkpm_moms_diag_d_old[output_moms[i] * 3] = (double)output_data_real->vals[count_old][i * 3];
        pkpm_moms_diag_d_old[(output_moms[i] * 3) + 1] =
          (double)output_data_real->vals[count_old][(i * 3) + 1];
        pkpm_moms_diag_d_old[(output_moms[i] * 3) + 2] =
          (double)output_data_real->vals[count_old][(i * 3) + 2];
      }
    }

    count_old += 1;
  }
}
