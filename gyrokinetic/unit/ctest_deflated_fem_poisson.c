#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_dg_bin_ops.h>


double calc_l2(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_range range_ext, struct gkyl_basis basis, struct gkyl_array* field1, struct gkyl_array* field2)
{
  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range_ext.volume);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *f1 = gkyl_array_cfetch(field1, lidx);
    const double *f2 = gkyl_array_cfetch(field2, lidx);
    double *diff_i = gkyl_array_fetch(diff, lidx);
    for(int i = 0; i <basis.num_basis; i++){
      diff_i[i] = f1[i] - f2[i];
    }
  }
  struct gkyl_array *l2_diff = gkyl_array_new(GKYL_DOUBLE, 1, range_ext.volume);
  gkyl_dg_calc_l2_range(&basis, 0, l2_diff, 0, diff, range);
  gkyl_array_scale_range(l2_diff, grid.cellVolume, &range);
  double l2[1];
  gkyl_array_reduce_range(l2, l2_diff, GKYL_SUM, &range);
  gkyl_array_release(diff);
  gkyl_array_release(l2_diff);
  return sqrt(l2[0]);
}


// functions for the charge density
void
rho_func_zdep_nd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 4*cos(z)*cos(2*x);
}

void
phi_func_zdep_nd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = cos(z)*cos(2*x);
}

void
rho_func_simplez_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 4*z*cos(2*x - M_PI/2);
}

void
phi_func_simplez_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*cos(2*x - M_PI/2);
}

void
rho_func_zind_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 4*cos(2*x - M_PI/2);
}

void
phi_func_zind_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = cos(2*x - M_PI/2);
}

void
rho_func_zind_dd_1x(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 4*cos(2*x - M_PI/2);
}

void
phi_func_zind_dd_1x(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  fout[0] = cos(2*x - M_PI/2);
}

void
rho_func_3x_dd_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = (4*cos(2*x)*sin(y) + 4*sin(2*x)*cos(y) + cos(2*x)*sin(y)) ;
  //fout[0] = 4*cos(2*y);
}

void
phi_func_3x_dd_dd(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = cos(2*x)*sin(y);
  //fout[0] = cos(2*y);
}



double
test_zdep_nd_nxnz(int nx, int ny){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, -M_PI}, upper[] = { 3*M_PI/4, M_PI };
  int cells[] = { nx, ny };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // project initial function on 2d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &rho_func_zdep_nd, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field_discont, "in_field.gkyl");

  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_dev, field);
  struct gkyl_array *field_discont_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_discont_dev, field_discont);
#else
  struct gkyl_array *field_dev = field;
  struct gkyl_array *field_discont_dev = field_discont;
#endif

  //smooth it
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &grid, &basis, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, 0, 0, 0, use_gpu);
  gkyl_fem_parproj_set_rhs(parproj, field_discont_dev, field_discont_dev);
  gkyl_fem_parproj_solve(parproj, field_dev);

  struct gkyl_poisson_bc poisson_bc = {0};
  poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, 2.0, 0); 
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *epsilon_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(epsilon_dev, epsilon);
  struct gkyl_array *phi_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *phi_dev = phi;
  struct gkyl_array *epsilon_dev = epsilon;
#endif

  struct gkyl_deflated_fem_poisson* deflated_fem_poisson = gkyl_deflated_fem_poisson_new(grid, basis_on_dev, basis, 
    local, local, epsilon_dev, 0, poisson_bc, NULL, use_gpu);
  gkyl_deflated_fem_poisson_advance(deflated_fem_poisson, field_dev, NULL, phi_dev);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(phi, phi_dev);
#endif
  //gkyl_grid_sub_array_write(&grid, &local, 0, phi, "out_field.gkyl");

  // project analytic solution
  struct gkyl_array *sol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_sol = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &phi_func_zdep_nd, 0);
  gkyl_proj_on_basis_advance(proj_sol, 0.0, &local, sol);
  gkyl_proj_on_basis_release(proj_sol);
  //gkyl_grid_sub_array_write(&grid, &local, 0, sol, "sol_field.gkyl");

  double l2 = calc_l2(grid, local,local_ext, basis, phi, sol);

  gkyl_deflated_fem_poisson_release(deflated_fem_poisson);
  gkyl_fem_parproj_release(parproj);
#ifdef GKYL_HAVE_CUDA 
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_discont_dev);
  gkyl_array_release(phi_dev);
  gkyl_array_release(epsilon_dev);
#endif
  gkyl_array_release(field_discont);
  gkyl_array_release(field);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(sol);
  return l2;
}

double
test_simplez_dd_nxnz(int nx, int ny){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, -M_PI }, upper[] = { M_PI, M_PI };
  int cells[] = { nx, ny };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // project initial function on 2d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &rho_func_simplez_dd, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field_discont, "in_field.gkyl");

  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_dev, field);
  struct gkyl_array *field_discont_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_discont_dev, field_discont);
#else
  struct gkyl_array *field_dev = field;
  struct gkyl_array *field_discont_dev = field_discont;
#endif

  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &grid, &basis, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, 0, 0, 0, use_gpu);
  gkyl_fem_parproj_set_rhs(parproj, field_discont_dev, field_discont_dev);
  gkyl_fem_parproj_solve(parproj, field_dev);


  struct gkyl_poisson_bc poisson_bc = {0};
  poisson_bc.lo_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, 2.0, 0); 
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *epsilon_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(epsilon_dev, epsilon);
  struct gkyl_array *phi_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *phi_dev = phi;
  struct gkyl_array *epsilon_dev = epsilon;
#endif
                                            
  struct gkyl_deflated_fem_poisson* deflated_fem_poisson = gkyl_deflated_fem_poisson_new(grid, basis_on_dev, basis, 
    local, local, epsilon_dev, 0, poisson_bc, NULL, use_gpu);
  gkyl_deflated_fem_poisson_advance(deflated_fem_poisson, field_dev, NULL, phi_dev);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(phi, phi_dev);
#endif
  //gkyl_grid_sub_array_write(&grid, &local, 0, phi, "out_field.gkyl");

  // project analytic solution
  struct gkyl_array *sol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_sol = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &phi_func_simplez_dd, 0);
  gkyl_proj_on_basis_advance(proj_sol, 0.0, &local, sol);
  gkyl_proj_on_basis_release(proj_sol);
  //gkyl_grid_sub_array_write(&grid, &local, 0, sol, "sol_field.gkyl");

  double l2 = calc_l2(grid, local,local_ext, basis, phi, sol);

  gkyl_deflated_fem_poisson_release(deflated_fem_poisson);
  gkyl_fem_parproj_release(parproj);
#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_discont_dev);
  gkyl_array_release(phi_dev);
  gkyl_array_release(epsilon_dev);
#endif
  gkyl_array_release(field_discont);
  gkyl_array_release(field);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(sol);
  return l2;
}

double
test_zind_dd_nxnz(int nx, int ny){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, -1 }, upper[] = { M_PI, 1 };
  int cells[] = { nx, ny };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // project initial function on 2d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &rho_func_zind_dd, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field_discont, "in_field.gkyl");

  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_dev, field);
  struct gkyl_array *field_discont_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_discont_dev, field_discont);
#else
  struct gkyl_array *field_dev = field;
  struct gkyl_array *field_discont_dev = field_discont;
#endif

  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &grid, &basis, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, 0, 0, 0, use_gpu);
  gkyl_fem_parproj_set_rhs(parproj, field_discont_dev, field_discont_dev);
  gkyl_fem_parproj_solve(parproj, field_dev);

  struct gkyl_poisson_bc poisson_bc = {0};
  poisson_bc.lo_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, 2.0, 0); 
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *epsilon_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(epsilon_dev, epsilon);
  struct gkyl_array *phi_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *phi_dev = phi;
  struct gkyl_array *epsilon_dev = epsilon;
#endif
                                            
  struct gkyl_deflated_fem_poisson* deflated_fem_poisson = gkyl_deflated_fem_poisson_new(grid, basis_on_dev, basis, 
    local, local, epsilon_dev, 0, poisson_bc, NULL, use_gpu);
  gkyl_deflated_fem_poisson_advance(deflated_fem_poisson, field_dev, NULL, phi_dev);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(phi, phi_dev);
#endif
  //gkyl_grid_sub_array_write(&grid, &local, 0, phi, "out_field.gkyl");

  // project analytic solution
  struct gkyl_array *sol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_sol = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &phi_func_zind_dd, 0);
  gkyl_proj_on_basis_advance(proj_sol, 0.0, &local, sol);
  gkyl_proj_on_basis_release(proj_sol);
  //gkyl_grid_sub_array_write(&grid, &local, 0, sol, "sol_field.gkyl");


  double l2 = calc_l2(grid, local,local_ext, basis, phi, sol);

  gkyl_deflated_fem_poisson_release(deflated_fem_poisson);
  gkyl_fem_parproj_release(parproj);
#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_discont_dev);
  gkyl_array_release(phi_dev);
  gkyl_array_release(epsilon_dev);
#endif
  gkyl_array_release(field_discont);
  gkyl_array_release(field);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(sol);
  return l2;
}


double
test_3x_dd_dd_nxnynz(int nx, int ny, int nz){
  // create the 2d field
  // create xz grid
  double lower[] = { -3*M_PI/4, -M_PI, -M_PI}, upper[] = { 3*M_PI/4, M_PI, M_PI };
  int cells[] = { nx, ny ,nz};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 ,1};
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 3, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // project initial function on 3d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &rho_func_3x_dd_dd, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  //gkyl_grid_sub_array_write(&grid, &local, 0, field_discont, "in_field.gkyl");

  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_dev, field);
  struct gkyl_array *field_discont_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_discont_dev, field_discont);
#else
  struct gkyl_array *field_dev = field;
  struct gkyl_array *field_discont_dev = field_discont;
#endif

  //smooth it
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &grid, &basis, GKYL_FEM_PARPROJ_DIRICHLET_SKIN, 0, 0, 0, use_gpu);
  gkyl_fem_parproj_set_rhs(parproj, field_discont_dev, field_discont_dev);
  gkyl_fem_parproj_solve(parproj, field_dev);

  //struct gkyl_poisson_bc poisson_bc;
  //poisson_bc.lo_type[0] = GKYL_POISSON_NEUMANN;
  //poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  //poisson_bc.lo_value[0].v[0] = 0.;
  //poisson_bc.up_value[0].v[0] = 0.;

  //poisson_bc.lo_type[1] = GKYL_POISSON_DIRICHLET;
  //poisson_bc.up_type[1] = GKYL_POISSON_DIRICHLET;
  //poisson_bc.lo_value[1].v[0] = 0.;
  //poisson_bc.up_value[1].v[0] = 0.;

  struct gkyl_poisson_bc poisson_bc = {0};
  poisson_bc.lo_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.up_type[0] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_type[1] = GKYL_POISSON_DIRICHLET;
  poisson_bc.up_type[1] = GKYL_POISSON_DIRICHLET;
  poisson_bc.lo_value[0].v[0] = 0.;
  poisson_bc.up_value[0].v[0] = 0.;
  poisson_bc.lo_value[1].v[0] = 0.;
  poisson_bc.up_value[1].v[0] = 0.;

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);
  gkyl_array_shiftc(epsilon, sqrt(pow(2,3)), 0); 
  gkyl_array_shiftc(epsilon, sqrt(pow(2,3)), basis.num_basis); 
  gkyl_array_shiftc(epsilon, sqrt(pow(2,3)), 2*basis.num_basis); 
  //gkyl_grid_sub_array_write(&grid, &local, 0, epsilon, "epsilon.gkyl");
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *epsilon_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);
  gkyl_array_copy(epsilon_dev, epsilon);
  struct gkyl_array *phi_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *phi_dev = phi;
  struct gkyl_array *epsilon_dev = epsilon;
#endif

  struct gkyl_deflated_fem_poisson* deflated_fem_poisson = gkyl_deflated_fem_poisson_new(grid, basis_on_dev, basis, 
    local, local, epsilon_dev, 0, poisson_bc, NULL, use_gpu);
  gkyl_deflated_fem_poisson_advance(deflated_fem_poisson, field_dev, NULL, phi_dev);
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(phi, phi_dev);
#endif
  //gkyl_grid_sub_array_write(&grid, &local, 0, phi, "out_field.gkyl");

  // project analytic solution
  struct gkyl_array *sol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj_sol = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &phi_func_3x_dd_dd, 0);
  gkyl_proj_on_basis_advance(proj_sol, 0.0, &local, sol);
  gkyl_proj_on_basis_release(proj_sol);
  //gkyl_grid_sub_array_write(&grid, &local, 0, sol, "sol_field.gkyl");

  double l2 = calc_l2(grid, local,local_ext, basis, phi, sol);

  gkyl_deflated_fem_poisson_release(deflated_fem_poisson);
  gkyl_fem_parproj_release(parproj);
#ifdef GKYL_HAVE_CUDA 
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_discont_dev);
  gkyl_array_release(phi_dev);
  gkyl_array_release(epsilon_dev);
#endif
  gkyl_array_release(field_discont);
  gkyl_array_release(field);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  gkyl_array_release(sol);
  return l2;
}



void test_zind_dd(){
  double l2s[6];
  int ny = 32;
  int i = 0;
  for(int nx = 4; nx < 129; nx*=2){
    l2s[i] = test_zind_dd_nxnz(nx,ny);
    if (i > 0)
      TEST_CHECK( l2s[i-1]>3.5*l2s[i] );
    i+=1;
  }
}

void test_simplez_dd(){
  double l2s[6];
  int ny = 32;
  int i = 0;
  for(int nx = 4; nx < 129; nx*=2){
    l2s[i] = test_simplez_dd_nxnz(nx,ny);
    if (i > 0)
      TEST_CHECK( l2s[i-1]>3.5*l2s[i] );
    i+=1;
  }
}

void test_zdep_nd(){
  // Expected results
  //double l2s[6] = { 1.4891748591339167, 0.4361776844752765, 0.1139546668200294, 0.0288104647768966, 0.0072320934639821, 0.0018431764053742};
  double l2s[6];
  int ny = 32;
  int i = 0;
  for(int nx = 4; nx < 129; nx*=2){
    l2s[i] = test_zdep_nd_nxnz(nx,ny);
    if (i > 0)
      TEST_CHECK( l2s[i-1]>3.0*l2s[i] );
    i+=1;
  }
}

void test_3x_dd_dd(){
  // Expected results
  //double l2s[6] = { 4.2333527815296019, 1.5169449008531153, 0.4618522366946782, 0.1324749882413162, 0.0438774243422054, 0.0212827374990034};
  double l2s[6];
  int ny = 32;
  int nz = 20;
  int i = 0;
  for(int nx = 4; nx < 129; nx*=2){
    l2s[i] = test_3x_dd_dd_nxnynz(nx,ny, nz);
    if (i > 0)
      TEST_CHECK( l2s[i-1]>2.0*l2s[i] );
    i+=1;
  }
}


TEST_LIST = {
  { "test_zind_dd", test_zind_dd},
  { "test_simplez_dd", test_simplez_dd},
  { "test_zdep_nd", test_zdep_nd},
  { "test_3x_dd_dd", test_3x_dd_dd},
  { NULL, NULL },
};
