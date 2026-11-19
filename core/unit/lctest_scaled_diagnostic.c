#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_dynvec.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_zero_lw.h>

#include <lauxlib.h>
#include <lualib.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void
write_field(const char *name, double mean, double higher)
{
  double lower[] = {0}, upper[] = {1};
  int cells[] = {1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, cells);
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 2, 1);
  double *v = gkyl_array_fetch(a, 0);
  v[0] = mean;
  v[1] = higher;
  TEST_CHECK(gkyl_grid_sub_array_write(&grid, &range, NULL, a, name) == 0);
  gkyl_array_release(a);
}

static void
write_integrated(const char *name, double m1)
{
  gkyl_dynvec v = gkyl_dynvec_new(GKYL_DOUBLE, 4);
  double values[] = {1e20, m1, 1e30, 1e30};
  gkyl_dynvec_append(v, 0, values);
  gkyl_dynvec_write(v, name);
  TEST_CHECK(access(name, F_OK) == 0);
  gkyl_dynvec_release(v);
}

static void
write_maxwellian(const char *name, double u0, double u1, double temperature)
{
  double lower[] = {0}, upper[] = {1};
  int cells[] = {1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, cells);
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, 6, 1);
  double *v = gkyl_array_fetch(a, 0);
  for (int k = 0; k < 6; ++k) {
    v[k] = 0;
  }
  v[2] = u0;
  v[3] = u1;
  v[4] = temperature;
  TEST_CHECK(gkyl_grid_sub_array_write(&grid, &range, NULL, a, name) == 0);
  gkyl_array_release(a);
}

static void
write_values(const char *name, int ncomp, const double *values)
{
  double lower[] = {0}, upper[] = {1};
  int cells[] = {1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, cells);
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, ncomp, 1);
  memcpy(gkyl_array_fetch(a, 0), values, ncomp * sizeof(double));
  TEST_CHECK(gkyl_grid_sub_array_write(&grid, &range, NULL, a, name) == 0);
  gkyl_array_release(a);
}

static void
write_vlasov_integrated(const char *name, int ncomp, double momentum, double energy)
{
  gkyl_dynvec v = gkyl_dynvec_new(GKYL_DOUBLE, ncomp);
  double values[5] = {1e20, 0, 0, 0, 0};
  for (int k = 1; k < ncomp - 1; ++k) {
    values[k] = momentum;
  }
  values[ncomp - 1] = energy;
  gkyl_dynvec_append(v, 0, values);
  gkyl_dynvec_write(v, name);
  gkyl_dynvec_release(v);
}

static void
write_pair_dynvec(const char *name, const double *values)
{
  gkyl_dynvec v = gkyl_dynvec_new(GKYL_DOUBLE, 2);
  gkyl_dynvec_append(v, 0, values);
  gkyl_dynvec_write(v, name);
  gkyl_dynvec_release(v);
}

static void
test_scale_aware_comparison(void)
{
  char dir[] = "/tmp/gkeyll-scale-diagnostic-XXXXXX";
  TEST_ASSERT(mkdtemp(dir) != NULL);
  char old[4096];
  TEST_ASSERT(getcwd(old, sizeof old) != NULL);
  TEST_ASSERT(chdir(dir) == 0);

  write_field("M0.gkyl", 1e20, 0);
  write_field("M2par.gkyl", 1e30, 0);
  write_field("baseline.gkyl", 0, 0);
  write_field("roundoff.gkyl", 0, 1e17);
  write_field("gradient.gkyl", 0, 1e20);
  write_field("mean.gkyl", 1e20, 0);
  write_field("M2perp.gkyl", 1e30, 0);
  write_field("M3-large.gkyl", 0, 1e29);
  write_maxwellian("maxwell-baseline.gkyl", 0, 0, 1e12);
  write_maxwellian("maxwell-small.gkyl", 0, 0.01, 1e12);
  write_maxwellian("maxwell-large.gkyl", 0, 1, 1e12);
  write_maxwellian("maxwell-temp.gkyl", 0, 0.01, 1e15);
  uint64_t nan_bits = UINT64_C(0x7ff8000000000000);
  double nan_value;
  memcpy(&nan_value, &nan_bits, sizeof nan_value);
  write_field("nonfinite.gkyl", nan_value, 0);
  write_integrated("integrated-baseline.gkyl", 0);
  write_integrated("integrated-small.gkyl", 1e17);
  write_integrated("integrated-large.gkyl", 1e20);
  write_integrated("integrated-nonfinite.gkyl", nan_value);
  write_values("vlasov-M0.gkyl", 2, (double[]){1e20, 0});
  write_values("vlasov-M2.gkyl", 2, (double[]){1e30, 0});
  write_values("vlasov-M1.gkyl", 4, (double[]){0, 0, 0, 0});
  write_values("vlasov-M1-small.gkyl", 4, (double[]){1e17, 1e17, 1e17, 1e17});
  write_values("vlasov-M1-large.gkyl", 4, (double[]){1e20, 0, 0, 0});
  write_values("vlasov-LTE.gkyl", 8, (double[]){1e20, 0, 0, 0, 0, 0, 1e12, 0});
  write_values("vlasov-LTE-small.gkyl", 8, (double[]){1e20, 0, 0.01, 0, 0.01, 0, 1e12, 0});
  write_values("vlasov-LTE-large.gkyl", 8, (double[]){1e20, 0, 1, 0, 0, 0, 1e12, 0});
  write_values("vlasov-LTE-temp.gkyl", 8, (double[]){1e20, 0, 0, 0, 0, 0, 2e12, 0});
  write_vlasov_integrated("vlasov-imom.gkyl", 4, 0, 1e30);
  write_vlasov_integrated("vlasov-imom-small.gkyl", 4, 1e17, 1e30);
  write_vlasov_integrated("vlasov-imom-large.gkyl", 4, 1e20, 1e30);
  write_vlasov_integrated("vlasov-imom-energy.gkyl", 4, 0, 2e30);
  write_vlasov_integrated("vlasov-imom-1v.gkyl", 3, 1e17, 1e30);
  write_vlasov_integrated("vlasov-imom-1v-base.gkyl", 3, 0, 1e30);
  write_vlasov_integrated("vlasov-imom-3v.gkyl", 5, 1e17, 1e30);
  write_vlasov_integrated("vlasov-imom-3v-base.gkyl", 5, 0, 1e30);
  write_vlasov_integrated("vlasov-imom-bad.gkyl", 4, nan_value, 1e30);
  double mixed_base[] = {2.080718738636528e19, -3.5479266818092526e-16};
  double mixed_close[] = {2.0807187386365293e19, 8.2912128896862711e-18};
  double mixed_far[] = {3.080718738636528e19, 8.2912128896862711e-18};
  write_values("mixed-base.gkyl", 2, mixed_base);
  write_values("mixed-close.gkyl", 2, mixed_close);
  write_values("mixed-far.gkyl", 2, mixed_far);
  write_pair_dynvec("mixed-base-dyn.gkyl", mixed_base);
  write_pair_dynvec("mixed-close-dyn.gkyl", mixed_close);
  write_pair_dynvec("mixed-far-dyn.gkyl", mixed_far);

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  gkyl_zero_lw_openlibs(L);
  int status = luaL_dostring(
    L,
    "local z=G0.Zero\n"
    "local g,b=z.arrayNewFromFile('baseline.gkyl')\n"
    "local _,r=z.arrayNewFromFile('roundoff.gkyl')\n"
    "local _,h=z.arrayNewFromFile('gradient.gkyl')\n"
    "local _,v=z.arrayNewFromFile('mean.gkyl')\n"
    "local _,bad=z.arrayNewFromFile('nonfinite.gkyl')\n"
    "local _,n=z.arrayNewFromFile('M0.gkyl')\n"
    "local _,p=z.arrayNewFromFile('M2par.gkyl')\n"
    "local _,q=z.arrayNewFromFile('M2perp.gkyl')\n"
    "local _,m3=z.arrayNewFromFile('M3-large.gkyl')\n"
    "local _,mb=z.arrayNewFromFile('maxwell-baseline.gkyl')\n"
    "local _,ms=z.arrayNewFromFile('maxwell-small.gkyl')\n"
    "local _,ml=z.arrayNewFromFile('maxwell-large.gkyl')\n"
    "local _,mt=z.arrayNewFromFile('maxwell-temp.gkyl')\n"
    "assert(z.scaledMomentDiff(b,r,'M1',1e-7,n,p,nil,g))\n"
    "assert(not z.scaledMomentDiff(b,h,'M1',1e-7,n,p,nil,g))\n"
    "assert(not z.scaledMomentDiff(b,v,'M1',1e-7,n,p,nil,g))\n"
    "assert(z.scaledMomentDiff(b,r,'M1',1e-7,nil,p,nil,g)==nil)\n"
    "assert(z.scaledMomentDiff(b,r,'M3par',1e-7,n,p,nil,g))\n"
    "assert(not z.scaledMomentDiff(b,m3,'M3par',1e-7,n,p,nil,g))\n"
    "assert(z.scaledMomentDiff(b,r,'M3perp',1e-7,n,p,q,g))\n"
    "assert(not z.scaledMomentDiff(b,m3,'M3perp',1e-7,n,p,q,g))\n"
    "assert(z.scaledMomentDiff(mb,ms,'Maxwellian',1e-7,nil,nil,nil,g))\n"
    "assert(not z.scaledMomentDiff(mb,ml,'Maxwellian',1e-7,nil,nil,nil,g))\n"
    "assert(not z.scaledMomentDiff(mb,mt,'Maxwellian',1e-7,nil,nil,nil,g))\n"
    "assert(not z.arrayFinite(bad))\n"
    "assert(not z.scaledMomentDiff(b,bad,'M1',1e-7,n,p,nil,g))\n"
    "assert(z.scaledIntegratedM1Diff('integrated-baseline.gkyl','integrated-small.gkyl',1e-7))\n"
    "assert(not z.scaledIntegratedM1Diff('integrated-baseline.gkyl','integrated-large.gkyl',1e-7))\n"
    "local ok,ratio=z.scaledIntegratedM1Diff('integrated-baseline.gkyl','integrated-nonfinite.gkyl',1e-7)\n"
    "assert(ok==false and ratio==math.huge)\n"
    "local _,vn=z.arrayNewFromFile('vlasov-M0.gkyl')\n"
    "local _,ve=z.arrayNewFromFile('vlasov-M2.gkyl')\n"
    "local _,vb=z.arrayNewFromFile('vlasov-M1.gkyl')\n"
    "local _,vs=z.arrayNewFromFile('vlasov-M1-small.gkyl')\n"
    "local _,vl=z.arrayNewFromFile('vlasov-M1-large.gkyl')\n"
    "local _,vt=z.arrayNewFromFile('vlasov-LTE.gkyl')\n"
    "local _,vts=z.arrayNewFromFile('vlasov-LTE-small.gkyl')\n"
    "local _,vtl=z.arrayNewFromFile('vlasov-LTE-large.gkyl')\n"
    "local _,vtt=z.arrayNewFromFile('vlasov-LTE-temp.gkyl')\n"
    "assert(z.scaledVlasovMomentDiff(vb,vs,'M1',1e-7,1e-7,vn,ve,nil,g))\n"
    "assert(not z.scaledVlasovMomentDiff(vb,vl,'M1',1e-7,1e-7,vn,ve,nil,g))\n"
    "assert(z.scaledVlasovMomentDiff(vb,vs,'M1',1e-7,1e-7,vn,nil,vt,g))\n"
    "assert(z.scaledVlasovMomentDiff(vb,vs,'M1',1e-7,1e-7,vn,nil,nil,g)==nil)\n"
    "assert(z.scaledVlasovMomentDiff(vt,vts,'LTE',1e-7,1e-7,vn,nil,nil,g))\n"
    "assert(not z.scaledVlasovMomentDiff(vt,vtl,'LTE',1e-7,1e-7,vn,nil,nil,g))\n"
    "assert(not z.scaledVlasovMomentDiff(vt,vtt,'LTE',1e-7,1e-7,vn,nil,nil,g))\n"
    "assert(z.scaledVlasovIntegratedDiff('vlasov-imom.gkyl','vlasov-imom-small.gkyl',1e-7,1e-7))\n"
    "assert(not z.scaledVlasovIntegratedDiff('vlasov-imom.gkyl','vlasov-imom-large.gkyl',1e-7,1e-7))\n"
    "assert(not z.scaledVlasovIntegratedDiff('vlasov-imom.gkyl','vlasov-imom-energy.gkyl',1e-7,1e-7))\n"
    "assert(z.scaledVlasovIntegratedDiff('vlasov-imom-1v-base.gkyl','vlasov-imom-1v.gkyl',1e-7,1e-7))\n"
    "assert(z.scaledVlasovIntegratedDiff('vlasov-imom-3v-base.gkyl','vlasov-imom-3v.gkyl',1e-7,1e-7))\n"
    "local vok,vr=z.scaledVlasovIntegratedDiff('vlasov-imom.gkyl','vlasov-imom-bad.gkyl',1e-7,1e-7)\n"
    "assert(vok==false and vr==math.huge)\n"
    "local _,xb=z.arrayNewFromFile('mixed-base.gkyl')\n"
    "local _,xc=z.arrayNewFromFile('mixed-close.gkyl')\n"
    "local _,xf=z.arrayNewFromFile('mixed-far.gkyl')\n"
    "local xr=z.createGridRanges(g,{0,0,0,0,0,0,0})\n"
    "local xd=z.arrayDiff(xb,xc,xr,1e-7,1e-7)\n"
    "assert(math.abs(xd.min_abs_diff)>1e-7 and xd.max_rel_diff==1 and xd.within_combined_tol)\n"
    "assert(not z.arrayDiff(xb,xf,xr,1e-7,1e-7).within_combined_tol)\n"
    "local yd=z.dynvecDiff('mixed-base-dyn.gkyl','mixed-close-dyn.gkyl',1e-7,1e-7)\n"
    "assert(yd.is_compatible and math.abs(yd.min_abs_diff)>1e-7 and yd.max_rel_diff==1 and yd.within_combined_tol)\n"
    "assert(not z.dynvecDiff('mixed-base-dyn.gkyl','mixed-far-dyn.gkyl',1e-7,1e-7).within_combined_tol)\n"
  );
  TEST_CHECK(status == 0);
  TEST_MSG("%s", status ? lua_tostring(L, -1) : "");
  lua_close(L);

  remove("M0.gkyl");
  remove("M2par.gkyl");
  remove("baseline.gkyl");
  remove("roundoff.gkyl");
  remove("gradient.gkyl");
  remove("mean.gkyl");
  remove("M2perp.gkyl");
  remove("M3-large.gkyl");
  remove("maxwell-baseline.gkyl");
  remove("maxwell-small.gkyl");
  remove("maxwell-large.gkyl");
  remove("maxwell-temp.gkyl");
  remove("nonfinite.gkyl");
  remove("integrated-baseline.gkyl");
  remove("integrated-small.gkyl");
  remove("integrated-large.gkyl");
  remove("integrated-nonfinite.gkyl");
  remove("mixed-base.gkyl");
  remove("mixed-close.gkyl");
  remove("mixed-far.gkyl");
  remove("mixed-base-dyn.gkyl");
  remove("mixed-close-dyn.gkyl");
  remove("mixed-far-dyn.gkyl");
  const char *vlasov_files[] = {"vlasov-M0.gkyl",           "vlasov-M2.gkyl",
                                "vlasov-M1.gkyl",           "vlasov-M1-small.gkyl",
                                "vlasov-M1-large.gkyl",     "vlasov-LTE.gkyl",
                                "vlasov-LTE-small.gkyl",    "vlasov-LTE-large.gkyl",
                                "vlasov-LTE-temp.gkyl",     "vlasov-imom.gkyl",
                                "vlasov-imom-small.gkyl",   "vlasov-imom-large.gkyl",
                                "vlasov-imom-energy.gkyl",  "vlasov-imom-1v.gkyl",
                                "vlasov-imom-1v-base.gkyl", "vlasov-imom-3v.gkyl",
                                "vlasov-imom-3v-base.gkyl", "vlasov-imom-bad.gkyl"};
  for (size_t i = 0; i < sizeof(vlasov_files) / sizeof(vlasov_files[0]); ++i) {
    remove(vlasov_files[i]);
  }
  TEST_CHECK(chdir(old) == 0);
  TEST_CHECK(rmdir(dir) == 0);
}

TEST_LIST = {{"scale_aware_comparison", test_scale_aware_comparison}, {NULL, NULL}};
