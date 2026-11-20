#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#include <stdbool.h>

// Object type.
typedef struct gkyl_fem_gyroaverage gkyl_fem_gyroaverage;

// Object type of a batch of gyroaverage updaters on the GPU.
typedef struct gkyl_fem_gyroaverage_batch gkyl_fem_gyroaverage_batch;

// Inputs to create a new gyroaverage updater.
struct gkyl_fem_gyroaverage_inp {
  const struct gkyl_rect_grid *grid; // Grid of the plane perpendicular to B: (x) or (x,y).
  const struct gkyl_basis *basis; // Basis on that grid (piecewise linear).
  const struct gkyl_range *local; // Range the operator acts on (the whole plane).
  bool is_periodic[2]; // Whether each direction is periodic.
  // Extend the field beyond non-periodic boundaries by reflection, f(bnd - s) = f(bnd + s),
  // instead of linearly about its boundary value (see gkyl_fem_gyroaverage_inew).
  bool reflect_bnd;
  const struct gkyl_array *gxx; // Contravariant metric coefficient g^xx.
  const struct gkyl_array *gxy; // Contravariant metric coefficient g^xy (2D only).
  const struct gkyl_array *gyy; // Contravariant metric coefficient g^yy (2D only).
  const struct gkyl_array *bmag; // Magnetic field amplitude.
  double mass; // Species mass.
  double charge; // Species charge.
  int num_mu; // Number of magnetic moments to gyroaverage at.
  const double *mu; // Magnetic moments (m*vperp^2/(2*B)) to gyroaverage at.
  // Whether to run on the GPU. The geometry is given on the host either way (the operator is
  // assembled on the host), and the fields it is applied to are then on the device.
  bool use_gpu;
};

/**
 * Create a new updater that gyroaverages a field in the plane perpendicular to
 * the magnetic field, at a list of magnetic moments mu_n,
 *
 *   <f>_n(X) = 1/(2*pi) int_0^{2*pi} f(X + rho(X,mu_n,alpha)) dalpha.
 *
 * The gyroradius vector is mapped to computational coordinates with the
 * metric at the center of the cell containing X,
 *
 *   rho.grad(x) = rho*sqrt(g^xx)*cos(alpha),
 *   rho.grad(y) = rho*(g^xy*cos(alpha) + sqrt(g^xx*g^yy-(g^xy)^2)*sin(alpha))/sqrt(g^xx),
 *
 * with rho = sqrt(2*m*mu/B)/|q|. The operator is discretized in weak form on
 * the continuous piecewise linear (FEM) space: with phi_i the nodal basis,
 *
 *   B_ij = int phi_i(X) <phi_j>(X) dX,   M_ij = int phi_i phi_j dX,
 *
 * the gyroaverage of f is E*M^{-1}*B*M^{-1}*E^T*M_dg*f, where E interpolates
 * nodal values to the DG basis and M_dg is the DG mass matrix. The element
 * integrals in B are computed exactly, which keeps the operator accurate when
 * the gyroradius is smaller than a cell. In a non-periodic direction the
 * field is extended beyond the boundary linearly about its boundary value,
 * f(bnd - s) = 2*f(bnd) - f(bnd + s), so that the gyroaverage of a linear
 * function is that function up to the boundary; this is the odd reflection
 * of a field that vanishes at the boundary. The product of two such fields
 * is even about the boundary, so quantities quadratic in them are
 * gyroaveraged with the reflection f(bnd - s) = f(bnd + s) (reflect_bnd).
 *
 * @param inp Input parameters defined in gkyl_fem_gyroaverage_inp.
 * @return New updater pointer.
 */
struct gkyl_fem_gyroaverage *gkyl_fem_gyroaverage_inew(const struct gkyl_fem_gyroaverage_inp *inp);

/**
 * Gyroaverage a field at each magnetic moment.
 *
 * @param up Gyroaverage updater.
 * @param fin Input field (num_basis components).
 * @param fout Gyroaveraged field at each magnetic moment (num_mu*num_basis
 *             components, the coefficients for the n-th magnetic moment
 *             starting at n*num_basis).
 */
void gkyl_fem_gyroaverage_advance(
  struct gkyl_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Apply the adjoint of gkyl_fem_gyroaverage_advance, i.e. compute
 * sum_n G_n^dagger f_n, G_n being the gyroaverage at the n-th magnetic moment
 * and the adjoint taken in the L2 inner product, so that
 * sum_n int f_n G_n(g) dX = int (sum_n G_n^dagger f_n) g dX holds to machine
 * precision.
 *
 * @param up Gyroaverage updater.
 * @param fin Input field at each magnetic moment (num_mu*num_basis
 *            components, the coefficients for the n-th magnetic moment
 *            starting at n*num_basis).
 * @param fout Sum of the adjoint gyroaverages (num_basis components).
 */
void gkyl_fem_gyroaverage_advance_adjoint(
  struct gkyl_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Number of nonzero entries in the gyroaverage matrix B of a magnetic moment.
 *
 * @param up Gyroaverage updater.
 * @param mu_idx Index of the magnetic moment.
 * @return Number of nonzero entries.
 */
long gkyl_fem_gyroaverage_num_nonzeros(const struct gkyl_fem_gyroaverage *up, int mu_idx);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_gyroaverage_release(struct gkyl_fem_gyroaverage *up);

/**
 * Create a batch of gyroaverage updaters on the GPU. A batch applies an updater to each of
 * num_fld fields with a number of kernel launches that does not depend on num_fld. The
 * updaters must have the same grid, basis, range, periodicity and number of magnetic moments
 * (their geometry and magnetic moments may differ).
 *
 * @param ups Updater applied to each field (assembled on the host; the same updater may be
 *            given for several fields). They are not needed once the batch is created.
 * @param num_fld Number of fields.
 * @return New batch pointer.
 */
struct gkyl_fem_gyroaverage_batch *gkyl_fem_gyroaverage_batch_new(
  struct gkyl_fem_gyroaverage *const *ups, int num_fld
);

/**
 * Gyroaverage each field of a batch at each magnetic moment, on the GPU. The fields are
 * stored one after the other in the arrays: if size is the number of cells in the array of
 * one field, the s-th field is in the cells s*size to (s+1)*size-1.
 *
 * @param batch Batch of gyroaverage updaters.
 * @param fin Input fields (see gkyl_fem_gyroaverage_advance for each of them).
 * @param fout Gyroaveraged fields at each magnetic moment.
 */
void gkyl_fem_gyroaverage_batch_advance(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Apply the adjoint of the gyroaverage to each field of a batch, on the GPU. The fields are
 * stored as in gkyl_fem_gyroaverage_batch_advance.
 *
 * @param batch Batch of gyroaverage updaters.
 * @param fin Input fields at each magnetic moment (see
 *            gkyl_fem_gyroaverage_advance_adjoint for each of them).
 * @param fout Sums of the adjoint gyroaverages.
 */
void gkyl_fem_gyroaverage_batch_advance_adjoint(
  struct gkyl_fem_gyroaverage_batch *batch, const struct gkyl_array *fin, struct gkyl_array *fout
);

/**
 * Delete a batch.
 *
 * @param batch Batch to delete.
 */
void gkyl_fem_gyroaverage_batch_release(struct gkyl_fem_gyroaverage_batch *batch);
