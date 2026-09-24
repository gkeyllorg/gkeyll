#pragma once

// Private header, not for direct use in user code

#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>

/**
 * Free wave geometry object.
 *
 * @param ref Reference counter for wave geometry
 */
void gkyl_wave_geom_free(const struct gkyl_ref_count *ref);

static void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int i=0; i<3; ++i) xp[i] = xc[i];
}

// Computes 1D geometry
static void
calc_geom_1d_from_nodes(const double *dx, const double *xlp, const double *xrp,
    evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  geo->kappa = fabs(xrp[0]-xlp[0])/dx[0];
  geo->lenr[0] = 1.0;

  geo->norm[0][0] = 1.0; geo->norm[0][1] = 0.0; geo->norm[0][2] = 0.0;
  geo->tau1[0][0] = 0.0; geo->tau1[0][1] = 1.0; geo->tau1[0][2] = 0.0;
  geo->tau2[0][0] = 0.0; geo->tau2[0][1] = 0.0; geo->tau2[0][2] = 1.0;
}

static void
calc_geom_1d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  double xlc[GKYL_MAX_CDIM], xrc[GKYL_MAX_CDIM];
  double xlp[GKYL_MAX_CDIM], xrp[GKYL_MAX_CDIM];

  xlc[0] = xc[0]-0.5*dx[0]; // left node
  xrc[0] = xc[0]+0.5*dx[0]; // right node

  // compute coordinates of left/right nodes
  mapc2p(0.0, xlc, xlp, ctx);
  mapc2p(0.0, xrc, xrp, ctx);

  calc_geom_1d_from_nodes(dx, xlp, xrp, mapc2p, ctx, geo);
}

// Computes 2D geometry
static void
calc_geom_2d_from_nodes(const double *dx, const struct gkyl_vec3 xll_p,
  const struct gkyl_vec3 xlr_p, const struct gkyl_vec3 xul_p,
  const struct gkyl_vec3 xur_p, evalf_t mapc2p, void *ctx,
  struct gkyl_wave_cell_geom *geo)
{
  // ll: lower-left; lr: lower-right
  // ul: upper-left; ur: upper-right

  // volume factor
  double area = 0.5*gkyl_vec3_len( gkyl_vec3_cross(gkyl_vec3_sub(xlr_p,xll_p), gkyl_vec3_sub(xul_p,xll_p)) )
    + 0.5*gkyl_vec3_len( gkyl_vec3_cross(gkyl_vec3_sub(xlr_p,xur_p), gkyl_vec3_sub(xul_p,xur_p)) );

  geo->kappa = area/(dx[0]*dx[1]);

  // face-area ratios for faces (a face is an edge in 2D)
  geo->lenr[0] = gkyl_vec3_len(gkyl_vec3_sub(xul_p, xll_p))/dx[1];
  geo->lenr[1] = gkyl_vec3_len(gkyl_vec3_sub(xlr_p, xll_p))/dx[0];

  // normal-tangent to left face
  struct gkyl_vec3 tau1_l = gkyl_vec3_norm(gkyl_vec3_sub(xul_p, xll_p));
  struct gkyl_vec3 tau2_l = gkyl_vec3_new(0.0, 0.0, 1.0); // ez
  struct gkyl_vec3 norm_l = gkyl_vec3_cross(tau1_l, tau2_l);
  
  for (int d=0; d<3; ++d) {
    geo->norm[0][d] = norm_l.x[d];
    geo->tau1[0][d] = tau1_l.x[d];
    geo->tau2[0][d] = tau2_l.x[d];
  }

  // normal-tangent to bottom face
  struct gkyl_vec3 tau1_b = gkyl_vec3_norm(gkyl_vec3_sub(xlr_p, xll_p));
  struct gkyl_vec3 tau2_b = gkyl_vec3_new(0.0, 0.0, -1.0); // -ez (ensures normal points into cell)
  struct gkyl_vec3 norm_b = gkyl_vec3_cross(tau1_b, tau2_b);

  for (int d=0; d<3; ++d) {
    geo->norm[1][d] = norm_b.x[d];
    geo->tau1[1][d] = tau1_b.x[d];
    geo->tau2[1][d] = tau2_b.x[d];
  }
}

static void
calc_geom_2d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  // ll: lower-left; lr: lower-right
  // ul: upper-left; ur: upper-right
  
  struct gkyl_vec3 xll_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xlr_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xul_p = gkyl_vec3_zeros();
  struct gkyl_vec3 xur_p = gkyl_vec3_zeros();

  struct gkyl_vec3 xll_c = gkyl_vec3_new(xc[0] - 0.5*dx[0],  xc[1] - 0.5*dx[1], 0.0);
  struct gkyl_vec3 xlr_c = gkyl_vec3_new(xc[0] + 0.5*dx[0],  xc[1] - 0.5*dx[1], 0.0);

  struct gkyl_vec3 xul_c = gkyl_vec3_new(xc[0] - 0.5*dx[0],  xc[1] + 0.5*dx[1], 0.0);
  struct gkyl_vec3 xur_c = gkyl_vec3_new(xc[0] + 0.5*dx[0],  xc[1] + 0.5*dx[1], 0.0);
  
  mapc2p(0.0, xll_c.x, xll_p.x, ctx);
  mapc2p(0.0, xlr_c.x, xlr_p.x, ctx);
  mapc2p(0.0, xul_c.x, xul_p.x, ctx);
  mapc2p(0.0, xur_c.x, xur_p.x, ctx);

  // need to set the final coordinate to 0.0
  xll_p.x[2] = xlr_p.x[2] = xul_p.x[2] = xur_p.x[2] = 0.0;

  calc_geom_2d_from_nodes(dx, xll_p, xlr_p, xul_p, xur_p, mapc2p, ctx, geo);
}

// Computes 3D geometry
static double
vol_tetra(const struct gkyl_vec3 p1,
          const struct gkyl_vec3 p2,
          const struct gkyl_vec3 p3,
          const struct gkyl_vec3 p4)
{
  struct gkyl_vec3 a = gkyl_vec3_sub(p1, p4);
  struct gkyl_vec3 b = gkyl_vec3_sub(p2, p4);
  struct gkyl_vec3 c = gkyl_vec3_sub(p3, p4);
  double vol = gkyl_vec3_dot(a, gkyl_vec3_cross(b, c)) / 6.;
  return vol > 0? vol : -vol;
}

/* The order of vertices is assumed to be (l: lower, u: upper):
[0], lll: (i,   j,   k  ) or (i-0.5, j-0.5, k-0.5)
[1], ull: (i+1, j,   k  ) or (i+0.5, j-0.5, k-0.5)
[2], uul: (i+1, j+1, k  ) or (i+0.5, j+0.5, k-0.5)
[3], lul: (i,   j+1, k  ) or (i-0.5, j+0.5, k-0.5)
[4], llu: (i,   j,   k+1) or (i-0.5, j-0.5, k+0.5)
[5], ulu: (i+1, j,   k+1) or (i+0.5, j-0.5, k+0.5)
[6], uuu: (i+1, j+1, k+1) or (i+0.5, j+0.5, k+0.5)
[7], luu: (i,   j+1, k+1) or (i-0.5, j+0.5, k+0.5)

       y
3----------2
|\     ^   |\
| \    |   | \
|  \   |   |  \
|   7------+---6
|   |  *-- |-- | -> x
0---+---\--1   |
 \  |    \  \  |
  \ |     \  \ |
   \|      z  \|
    4----------5
*/
static double
vol_hexa(const struct gkyl_vec3 *verts)
{
  // split the hexahedron into five tetrahedrons and add up their volumes
  // FIXME does this handle bad hexahedrons?
  double vol = 0;
  vol += vol_tetra(verts[1], verts[5], verts[4], verts[6]);
  vol += vol_tetra(verts[7], verts[4], verts[6], verts[3]);
  vol += vol_tetra(verts[4], verts[0], verts[1], verts[3]);
  vol += vol_tetra(verts[6], verts[1], verts[3], verts[4]);
  vol += vol_tetra(verts[1], verts[2], verts[6], verts[3]);

  return vol;
}

static inline double
triangle_area(const struct gkyl_vec3 p1,
              const struct gkyl_vec3 p2,
              const struct gkyl_vec3 p3)
{
  return 0.5*gkyl_vec3_len(
      gkyl_vec3_cross(gkyl_vec3_sub(p1, p2), gkyl_vec3_sub(p2, p3)));
}

// Points are in anti-clockwise order, i.e., p1-p3 and p2-p4 are diagonals.
static double
planar_quad_area_norm(const struct gkyl_vec3 p1,
                        const struct gkyl_vec3 p2,
                        const struct gkyl_vec3 p3,
                        const struct gkyl_vec3 p4,
                        struct gkyl_vec3 *norm)
{
  struct gkyl_vec3 v13 = gkyl_vec3_sub(p3, p1);
  struct gkyl_vec3 v24 = gkyl_vec3_sub(p4, p2);

  double area = 0.5 * gkyl_vec3_len(gkyl_vec3_cross(v13, v24));

  if (norm) {
    struct gkyl_vec3 norm_ = gkyl_vec3_norm(gkyl_vec3_cross(v13, v24));

    for (int d=0; d<3; ++d)
    {
      norm->x[d] = norm_.x[d];
    }
  }

  return area;
}

// ca * a + cb * b
static inline struct gkyl_vec3
gkyl_vec3_add_coeff(
  const double ca, struct gkyl_vec3 a, const double cb, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = {
    ca*a.x[0] + cb*b.x[0],
    ca*a.x[1] + cb*b.x[1],
    ca*a.x[2] + cb*b.x[2] } };
}

static double
quad_area_norm_tang(const struct gkyl_vec3 p1,
                    const struct gkyl_vec3 p2,
                    const struct gkyl_vec3 p3,
                    const struct gkyl_vec3 p4,
                    struct gkyl_vec3 *norm,
                    struct gkyl_vec3 *tau1,
                    struct gkyl_vec3 *tau2)
{
  // find Varignon parallelogram and use its normal as the quad's normal
  struct gkyl_vec3 pp1 = gkyl_vec3_add_coeff(0.5, p1, 0.5, p2);
  struct gkyl_vec3 pp2 = gkyl_vec3_add_coeff(0.5, p2, 0.5, p3);
  struct gkyl_vec3 pp3 = gkyl_vec3_add_coeff(0.5, p3, 0.5, p4);
  struct gkyl_vec3 pp4 = gkyl_vec3_add_coeff(0.5, p4, 0.5, p1);
  double area = planar_quad_area_norm(pp1, pp2, pp3, pp4, norm);

  // tau1 is v(p1->p2) projected on the normal plane, (tau1, tau2, norm)
  // completes a high-hand system
  struct gkyl_vec3 p12 = gkyl_vec3_sub(p2, p1);
  double d = gkyl_vec3_dot(p12, *norm);
  struct gkyl_vec3 t1 = gkyl_vec3_norm(gkyl_vec3_add_coeff(1, p12, -d, *norm));
  struct gkyl_vec3 t2 = gkyl_vec3_cross(*norm, t1);
  for (int d=0; d<3; ++d)
  {
    tau1->x[d] = t1.x[d];
    tau2->x[d] = t2.x[d];
  }

#if 1
  // project corner points onto the normal plane and use the projected planar
  // quad's area as the original quad's area
  d = gkyl_vec3_dot(*norm, p1);
  struct gkyl_vec3 proj1 = gkyl_vec3_add_coeff(1, p1, -d, *norm);
  d = gkyl_vec3_dot(*norm, p2);
  struct gkyl_vec3 proj2 = gkyl_vec3_add_coeff(1, p2, -d, *norm);
  d = gkyl_vec3_dot(*norm, p3);
  struct gkyl_vec3 proj3 = gkyl_vec3_add_coeff(1, p3, -d, *norm);
  d = gkyl_vec3_dot(*norm, p4);
  struct gkyl_vec3 proj4 = gkyl_vec3_add_coeff(1, p4, -d, *norm);
  
  area = planar_quad_area_norm(proj1, proj2, proj3, proj4, norm);
#else
  // accumulate triangle areas
  area += triangle_area(p1, pp1, pp4);
  area += triangle_area(p2, pp2, pp1);
  area += triangle_area(p3, pp3, pp2);
  area += triangle_area(p4, pp4, pp3);
#endif

  return area;
}

static void
calc_geom_3d_form_nodes(const double *dx, struct gkyl_vec3 verts[8],
    evalf_t mapc2p, void *ctx, struct gkyl_wave_cell_geom *geo)
{
  // compute cell volume and kappa
  double vol = vol_hexa(verts);
  double cell_vol_c = dx[0] * dx[1] * dx[2];
  geo->kappa = vol / cell_vol_c;

  // for each of the thee lower quad faces 'owned' by the present cell, compute
  // area, norm, tan1, tan2, and lenr
  int pt_idx[3][4] = { // indices of verices of each face in verts; see vol_hexa
    // lower-x face, lll, lul, luu, llu; v(p1->p2)=ey
    {0, 3, 7, 4},
    // lower-y face, lll, llu, ulu, ull; v(p1->p2)=ez
    {0, 4, 5, 1},
    // lower-z face, lll, ull, uul, lul; v(p1->p2)=ex
    {0, 1, 2, 3}
  };
  for (int face_idx = 0; face_idx < 3; ++face_idx)
  {
    int ip1 = pt_idx[face_idx][0];
    int ip2 = pt_idx[face_idx][1];
    int ip3 = pt_idx[face_idx][2];
    int ip4 = pt_idx[face_idx][3];

    struct gkyl_vec3 norm, tau1, tau2;
    double area= quad_area_norm_tang(
        verts[ip1], verts[ip2], verts[ip3], verts[ip4], &norm, &tau1, &tau2);
    double cell_area_c = cell_vol_c / dx[face_idx];
    geo->lenr[face_idx] = area / cell_area_c;
    for (int d=0; d<3; ++d)
    {
      geo->norm[face_idx][d] = norm.x[d];
      geo->tau1[face_idx][d] = tau1.x[d];
      geo->tau2[face_idx][d] = tau2.x[d];
    }
  }
}

static void
calc_geom_3d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx,
             struct gkyl_wave_cell_geom *geo)
{
  // get all vertices of the hexahedron; see vol_hexa for their order
  struct gkyl_vec3 verts_c[8] = {
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]-0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]-0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]+0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]+0.5*dx[1], xc[2]-0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]-0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]-0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]+0.5*dx[0], xc[1]+0.5*dx[1], xc[2]+0.5*dx[2]),
    gkyl_vec3_new(xc[0]-0.5*dx[0], xc[1]+0.5*dx[1], xc[2]+0.5*dx[2])
  };

  struct gkyl_vec3 verts[8]; // physical coordinate nodes
  for (int i=0; i<8; ++i)
    mapc2p(0.0, verts_c[i].x, verts[i].x, ctx);

  calc_geom_3d_form_nodes(dx, verts, mapc2p, ctx, geo);
}

