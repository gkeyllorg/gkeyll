#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// Allocate the reduction buffers and dynvecs, once the moment calculators exist.
static void
gk_integ_moms_alloc(struct gkyl_gyrokinetic_app *app, struct gk_integ_moms *im)
{
  // The components of the moments are stored one after the other in the
  // reduction buffers.
  im->num_mom_all = 0;
  for (int m = 0; m < im->num; m++) {
    im->offset[m] = im->num_mom_all;
    im->num_mom_all += im->ops[m].num_mom;
  }

  im->red_local = gkyl_calloc_hodev(app->use_gpu, im->num_mom_all, sizeof(double));
  im->red_global = gkyl_calloc_hodev(app->use_gpu, im->num_mom_all, sizeof(double));

  im->diag = gkyl_malloc(sizeof(gkyl_dynvec[im->num]));
  for (int m = 0; m < im->num; m++) {
    im->diag[m] = gkyl_dynvec_new(GKYL_DOUBLE, im->ops[m].num_mom);
  }
  im->is_first_write = true;
}

void
gk_integ_moms_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_integ_moms *im, int num,
  const enum gkyl_distribution_moments *types
)
{
  assert(num <= GKYL_F_MOMENT_COUNT);
  *im = (struct gk_integ_moms){};
  im->num = num;
  if (num == 0) {
    return;
  }

  im->ops = gkyl_malloc(sizeof(struct gk_species_moment[num]));
  for (int m = 0; m < num; m++) {
    gk_species_moment_init(app, gks, &im->ops[m], types[m], GKYL_F_MOMENT_WEIGHT_NONE, true);
  }
  gk_integ_moms_alloc(app, im);
}

void
gk_integ_moms_neut_init(
  struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, struct gk_integ_moms *im, int num,
  const enum gkyl_distribution_moments *types
)
{
  assert(num <= GKYL_F_MOMENT_COUNT);
  *im = (struct gk_integ_moms){};
  im->num = num;
  if (num == 0) {
    return;
  }

  im->ops = gkyl_malloc(sizeof(struct gk_species_moment[num]));
  for (int m = 0; m < num; m++) {
    gk_neut_species_moment_init(app, gkns, &im->ops[m], types[m], true);
  }
  gk_integ_moms_alloc(app, im);
}

int
gk_integ_moms_find(const struct gk_integ_moms *im, enum gkyl_distribution_moments type)
{
  for (int m = 0; m < im->num; m++) {
    if (im->ops[m].type == type) {
      return m;
    }
  }
  return -1;
}

void
gk_integ_moms_calc(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, const struct gk_integ_moms *im,
  const struct gkyl_array *fin
)
{
  for (int m = 0; m < im->num; m++) {
    gk_species_moment_calc(app, &im->ops[m], gks->local, app->local, fin);
  }
  app->stat.n_mom += 1;
}

void
gk_integ_moms_reduce_append(
  gkyl_gyrokinetic_app *app, struct gk_integ_moms *im, struct gkyl_array *const *arrs, double tm
)
{
  int num_mom_all = im->num_mom_all;
  // Each reduction writes the components of its moment, at an offset in red_local.
  for (int m = 0; m < im->num; m++) {
    gkyl_array_reduce_range(im->red_local + im->offset[m], arrs[m], GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom_all, im->red_local, im->red_global);

  double avals_global[num_mom_all];
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, im->red_global, sizeof(double[num_mom_all]), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(avals_global, im->red_global, sizeof(double[num_mom_all]));
  }
  for (int m = 0; m < im->num; m++) {
    gkyl_dynvec_append(im->diag[m], tm, &avals_global[im->offset[m]]);
  }
}

void
gk_integ_moms_calc_append(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_integ_moms *im,
  const struct gkyl_array *fin, double tm
)
{
  gk_integ_moms_calc(app, gks, im, fin);
  struct gkyl_array *arrs[GKYL_MAX2(im->num, 1)];
  for (int m = 0; m < im->num; m++) {
    arrs[m] = im->ops[m].marr;
  }
  gk_integ_moms_reduce_append(app, im, arrs, tm);
}

void
gk_integ_moms_write(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, struct gk_integ_moms *im,
  const char *tag, const char *description
)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int m = 0; m < im->num; m++) {
    if (rank == 0) {
      // File named after the moment, e.g. <name>-<species>_integrated_M0.gkyl.
      const char *fmt = "%s-%s_%sintegrated_%s.gkyl";
      const char *mom_name = gkyl_distribution_moments_strs[im->ops[m].type];
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, tag, mom_name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, tag, mom_name);

      if (im->is_first_write) {
        struct gkyl_msgpack_map_elem io_meta_desc[] = {
          {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = description},
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_desc
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(im->diag[m], fileNm, mt);
        gkyl_msgpack_data_release(mt);
      } else {
        gkyl_dynvec_awrite(im->diag[m], fileNm);
      }
    }
    gkyl_dynvec_clear(im->diag[m]);
    app->stat.n_diag_io += 1;
  }
  im->is_first_write = false;
}

void
gk_integ_moms_release(const struct gkyl_gyrokinetic_app *app, const struct gk_integ_moms *im)
{
  if (im->num == 0) {
    return;
  }

  for (int m = 0; m < im->num; m++) {
    gk_species_moment_release(app, &im->ops[m]);
    gkyl_dynvec_release(im->diag[m]);
  }
  gkyl_free(im->ops);
  gkyl_free(im->diag);
  gkyl_free_hodev(app->use_gpu, im->red_local);
  gkyl_free_hodev(app->use_gpu, im->red_global);
}
