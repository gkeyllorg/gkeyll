#include <acutest.h>
#include <gkyl_tok_geo_wall_target_priv.h>

#include <math.h>
#include <stddef.h>

// Borrowed analytic coordinates only; no EFIT allocation/release is needed.
static struct gkyl_efit
wall(int n, double *r, double *z)
{
  return (struct gkyl_efit) {
    .limiter_status=1, .limiter_n=n, .limiter_R=r, .limiter_Z=z,
    .rdim=4.0, .zdim=4.0,
  };
}

static void
test_order_and_wrap(void)
{
  double r[]={1,3,3,1}, z[]={-1,-1,1,1};
  struct gkyl_efit e=wall(4,r,z);
  int forward[]={0,1}, reverse[]={1,0};
  int wrap[]={3,0}, reverse_wrap[]={0,3};
  double bottom[]={2,-1}, right[]={3,0}, left[]={1,0}, top[]={2,1};
  TEST_CHECK(tok_wall_target_valid(&e,2,forward));
  TEST_CHECK(tok_wall_target_valid(&e,2,reverse));
  TEST_CHECK(tok_wall_target_valid(&e,2,wrap));
  TEST_CHECK(tok_wall_target_valid(&e,2,reverse_wrap));
  TEST_CHECK(tok_wall_target_cap(&e,2,forward,bottom,right));
  TEST_CHECK(tok_wall_target_cap(&e,2,reverse,right,bottom));
  TEST_CHECK(tok_wall_target_cap(&e,2,wrap,left,bottom));
  TEST_CHECK(tok_wall_target_cap(&e,2,reverse_wrap,bottom,left));
  TEST_CHECK(!tok_wall_target_contains(&e,2,forward,left));
  TEST_CHECK(!tok_wall_target_contains(&e,2,wrap,right));
  TEST_CHECK(!tok_wall_target_contains(&e,2,wrap,top));

  // Reversing the underlying wall winding does not change membership when
  // the caller supplies the corresponding segment indices.
  double rr[]={1,3,3,1}, zz[]={1,1,-1,-1};
  struct gkyl_efit reversed=wall(4,rr,zz);
  int same_target[]={1,2};
  TEST_CHECK(tok_wall_target_cap(&reversed,2,same_target,bottom,right));
  TEST_CHECK(!tok_wall_target_contains(&reversed,2,same_target,top));
}

static void
test_ordered_arclength(void)
{
  double r[]={1,3,3,1}, z[]={-1,-1,1,1};
  struct gkyl_efit e=wall(4,r,z);
  int forward[]={0,1}, reverse[]={1,0}, wrap[]={3,0}, reverse_wrap[]={0,3};
  double bottom[]={2,-1}, right[]={3,0}, left[]={1,0}, join[]={3,-1};
  double s=-123.0;
  TEST_CHECK(tok_wall_target_coordinate(&e,2,forward,bottom,&s));
  TEST_CHECK(fabs(s-1.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,forward,right,&s));
  TEST_CHECK(fabs(s-3.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,reverse,bottom,&s));
  TEST_CHECK(fabs(s-3.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,reverse,right,&s));
  TEST_CHECK(fabs(s-1.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,wrap,left,&s));
  TEST_CHECK(fabs(s-1.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,wrap,bottom,&s));
  TEST_CHECK(fabs(s-3.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,reverse_wrap,left,&s));
  TEST_CHECK(fabs(s-3.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,reverse_wrap,bottom,&s));
  TEST_CHECK(fabs(s-1.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,forward,join,&s));
  TEST_CHECK(fabs(s-2.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,2,reverse,join,&s));
  TEST_CHECK(fabs(s-2.0)<1e-13);
  TEST_CHECK(tok_wall_target_coordinate(&e,1,forward,bottom,&s));
  TEST_CHECK(fabs(s-1.0)<1e-13);

  s=-123.0;
  TEST_CHECK(!tok_wall_target_coordinate(&e,2,forward,left,&s));
  TEST_CHECK(s==-123.0);
  TEST_CHECK(!tok_wall_target_coordinate(&e,0,forward,bottom,&s));
  TEST_CHECK(s==-123.0);
  TEST_CHECK(!tok_wall_target_coordinate(&e,2,forward,bottom,NULL));
}

static void
test_invalid_declarations(void)
{
  double r[]={1,3,3,1}, z[]={-1,-1,1,1};
  struct gkyl_efit e=wall(4,r,z);
  int good[]={0,1}, full[]={0,1,2,3}, duplicate[]={0,0};
  int disconnected[]={0,2}, negative[]={-1}, out_of_range[]={4};
  int repeated_turn[]={0,1,0};
  double p[]={2,-1};
  TEST_CHECK(!tok_wall_target_valid(&e,0,good));
  TEST_CHECK(!tok_wall_target_valid(&e,-1,good));
  TEST_CHECK(!tok_wall_target_valid(&e,4,full));
  TEST_CHECK(!tok_wall_target_valid(&e,5,full)); // Reject before reading list.
  TEST_CHECK(!tok_wall_target_valid(&e,1,NULL));
  TEST_CHECK(!tok_wall_target_valid(&e,2,duplicate));
  TEST_CHECK(!tok_wall_target_valid(&e,2,disconnected));
  TEST_CHECK(!tok_wall_target_valid(&e,1,negative));
  TEST_CHECK(!tok_wall_target_valid(&e,1,out_of_range));
  TEST_CHECK(!tok_wall_target_valid(&e,3,repeated_turn));
  TEST_CHECK(!tok_wall_target_contains(&e,2,disconnected,p));
  TEST_CHECK(!tok_wall_target_cap(&e,4,full,p,p));

  for (int k=0; k<6; ++k) {
    struct gkyl_efit invalid=e;
    switch (k) {
      case 0: invalid.limiter_status=0; break;
      case 1: invalid.limiter_status=-1; break;
      case 2: invalid.limiter_n=2; break;
      case 3: invalid.limiter_R=NULL; break;
      case 4: invalid.limiter_Z=NULL; break;
      case 5: invalid.rdim=-1.0; break;
    }
    TEST_CHECK(!tok_wall_target_valid(&invalid,2,good));
    TEST_CHECK(!tok_wall_target_cap(&invalid,2,good,p,p));
  }
  TEST_CHECK(!tok_wall_target_valid(NULL,2,good));
  TEST_CHECK(!tok_wall_target_contains(NULL,2,good,p));
}

static void
test_simple_arc_geometry(void)
{
  // Consecutive indices alone do not guarantee a geometrically simple arc.
  double crossed_r[]={1,3,1,3,4}, crossed_z[]={0,2,2,0,-1};
  struct gkyl_efit crossed=wall(5,crossed_r,crossed_z);
  int three[]={0,1,2};
  TEST_CHECK(!tok_wall_target_valid(&crossed,3,three));

  double back_r[]={1,3,2,2,1}, back_z[]={0,0,0,2,2};
  struct gkyl_efit back=wall(5,back_r,back_z);
  int two[]={0,1};
  TEST_CHECK(!tok_wall_target_valid(&back,2,two));
  int reverse[]={1,0};
  TEST_CHECK(!tok_wall_target_valid(&back,2,reverse));

  double straight_r[]={1,2,3,3,1}, straight_z[]={0,0,0,2,2};
  struct gkyl_efit straight=wall(5,straight_r,straight_z);
  TEST_CHECK(tok_wall_target_valid(&straight,2,two));
  TEST_CHECK(tok_wall_target_valid(&straight,2,reverse));

  double zero_r[]={1,1,3,3,1}, zero_z[]={0,0,0,2,2};
  struct gkyl_efit zero=wall(5,zero_r,zero_z);
  int one[]={0};
  TEST_CHECK(!tok_wall_target_valid(&zero,1,one));

  // An explicit repeated closing vertex must not let n<limiter_n disguise
  // the entire closed outline as a proper open target subarc.
  double closed_r[]={1,3,3,1,1}, closed_z[]={0,0,2,2,0};
  struct gkyl_efit closed=wall(5,closed_r,closed_z);
  int four[]={0,1,2,3};
  TEST_CHECK(!tok_wall_target_valid(&closed,4,four));
  TEST_CHECK(tok_wall_target_valid(&closed,3,three));
}

static void
test_material_cap_across_concavity(void)
{
  // A rectangular vessel with an exterior notch 2<R<2.2, 1<Z<=3.
  // The marked target follows the notch down, across its floor, then up.
  double r[]={1,4,4,2.2,2.2,2,2,1};
  double z[]={0,0,3,3,1,1,3,3};
  struct gkyl_efit e=wall(8,r,z);
  int target[]={3,4,5}, reverse[]={5,4,3}, only_right[]={3};
  double a[]={2.2,2}, b[]={2,2}, chord_midpoint[]={2.1,2};
  double unmarked_wall[]={3,0}, interior[]={2.1,0.5};
  TEST_CHECK(tok_wall_target_cap(&e,3,target,a,b));
  TEST_CHECK(tok_wall_target_cap(&e,3,reverse,b,a));
  TEST_CHECK(tok_wall_target_cap(&e,3,target,a,a));
  double sa=-1.0,sb=-1.0;
  TEST_CHECK(tok_wall_target_coordinate(&e,3,target,a,&sa));
  TEST_CHECK(tok_wall_target_coordinate(&e,3,target,b,&sb));
  TEST_CHECK(fabs(sa-1.0)<1e-13);
  TEST_CHECK(fabs(sb-3.2)<1e-13);
  TEST_CHECK(fabs((sb-sa)-2.2)<1e-13); // Wall arc, not the 0.2 m chord.
  // The straight chord's midpoint is in the exterior notch and not on the
  // target. A true cap result deliberately grants no chord/volume exemption.
  TEST_CHECK(!tok_wall_target_contains(&e,3,target,chord_midpoint));
  TEST_CHECK(!tok_wall_target_cap(&e,3,target,a,chord_midpoint));
  TEST_CHECK(!tok_wall_target_cap(&e,3,target,a,unmarked_wall));
  TEST_CHECK(!tok_wall_target_cap(&e,3,target,a,interior));
  TEST_CHECK(!tok_wall_target_cap(&e,1,only_right,a,b));
  double start[]={2.2,3}, end[]={2,3}, past_end[]={2,3.001};
  TEST_CHECK(tok_wall_target_cap(&e,3,target,start,end));
  TEST_CHECK(!tok_wall_target_contains(&e,3,target,past_end));
}

static void
test_membership_tolerance(void)
{
  double r[]={1,3,3,1}, z[]={-1,-1,1,1};
  struct gkyl_efit e=wall(4,r,z);
  int target[]={0};
  const double tol=4e-8; // 1e-8*max(1,rdim,zdim), in metres.
  double near[]={2,-1+0.5*tol}, far[]={2,-1+2*tol};
  double near_end[]={1-0.5*tol,-1}, past_end[]={1-2*tol,-1};
  TEST_CHECK(tok_wall_target_contains(&e,1,target,near));
  TEST_CHECK(!tok_wall_target_contains(&e,1,target,far));
  TEST_CHECK(tok_wall_target_contains(&e,1,target,near_end));
  TEST_CHECK(!tok_wall_target_contains(&e,1,target,past_end));
}

static void
test_nonfinite_inputs(void)
{
  double r[]={1,3,3,1}, z[]={-1,-1,1,1};
  struct gkyl_efit e=wall(4,r,z);
  int target[]={0};
  double good[]={2,-1};
  double bad[][2]={{NAN,-1},{2,NAN},{INFINITY,-1},{-INFINITY,-1},
    {2,INFINITY},{2,-INFINITY}};
  for (size_t i=0; i<sizeof(bad)/sizeof(bad[0]); ++i) {
    TEST_CHECK(!tok_wall_target_contains(&e,1,target,bad[i]));
    TEST_CHECK(!tok_wall_target_cap(&e,1,target,good,bad[i]));
    TEST_CHECK(!tok_wall_target_cap(&e,1,target,bad[i],good));
    double s=-123.0;
    TEST_CHECK(!tok_wall_target_coordinate(&e,1,target,bad[i],&s));
    TEST_CHECK(s==-123.0);
  }
  TEST_CHECK(!tok_wall_target_contains(&e,1,target,NULL));
  TEST_CHECK(!tok_wall_target_cap(&e,1,target,NULL,good));
  TEST_CHECK(!tok_wall_target_cap(&e,1,target,good,NULL));
  // Even an unmarked nonfinite wall vertex invalidates the supplied outline.
  r[3]=NAN;
  TEST_CHECK(!tok_wall_target_valid(&e,1,target));
  r[3]=1.0; z[2]=INFINITY;
  TEST_CHECK(!tok_wall_target_valid(&e,1,target));
  z[2]=1.0; e.rdim=NAN;
  TEST_CHECK(!tok_wall_target_valid(&e,1,target));
  e.rdim=4.0; e.zdim=INFINITY;
  TEST_CHECK(!tok_wall_target_valid(&e,1,target));
}

TEST_LIST = {
  { "target_order_and_wrap", test_order_and_wrap },
  { "target_ordered_arclength", test_ordered_arclength },
  { "target_invalid_declarations", test_invalid_declarations },
  { "target_simple_arc_geometry", test_simple_arc_geometry },
  { "target_material_cap_across_concavity", test_material_cap_across_concavity },
  { "target_membership_tolerance", test_membership_tolerance },
  { "target_nonfinite_inputs", test_nonfinite_inputs },
  { NULL, NULL },
};
