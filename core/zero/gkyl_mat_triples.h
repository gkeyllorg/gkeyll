#pragma once

#include <stddef.h>
#include <stdbool.h>
#define COLMAJOR 0
#define ROWMAJOR 1

// A triple stores 'val' at (row, col)
struct gkyl_mtriple {
  size_t row, col;
  double val;  
};

/** Triples stores list of (i,j,val) for use in constructing sparse matrices */
typedef struct gkyl_mat_triples gkyl_mat_triples;
/** Iterator into triples */
typedef struct gkyl_mat_triples_iter gkyl_mat_triples_iter;

/**
 * Create new triples list in a matrix of shape (nr, nc).
 *
 * @param nr Number of rows
 * @param nc Number of cols
 * @return Pointer to new empty triples object.
 */
gkyl_mat_triples* gkyl_mat_triples_new(size_t nr, size_t nc);

/*
 * Set row-major/col-major ordering for triples.
 */
void gkyl_mat_triples_set_rowmaj_order(gkyl_mat_triples *tri);
void gkyl_mat_triples_set_colmaj_order(gkyl_mat_triples *tri);

bool gkyl_mat_triples_is_rowmaj(gkyl_mat_triples *tri);
bool gkyl_mat_triples_is_colmaj(gkyl_mat_triples *tri);

/**
 * Insert value 'val' in triples list at location (i,j)
 *
 * @return value inserted (val)
 */
GKYL_CU_DH double gkyl_mat_triples_insert(gkyl_mat_triples *tri, size_t i, size_t j, double val);

/**
 * Accumulate value 'val' in triples list at location (i,j). If an
 * element at this location exists it is incremented by 'val'. New
 * value at the location is returned.
 */
GKYL_CU_DH double gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val);

/**
 * Returns value in triples list at location (i,j)
 *
 * @return value at (i,j)
 */
double gkyl_mat_triples_get(const gkyl_mat_triples *tri, size_t i, size_t j);

/**
 * Return number of elements inserted.
 */
size_t gkyl_mat_triples_size(const gkyl_mat_triples *tri);

/**
 * Allocate and initialize a new iterator into triples. Release
 * calling gkyl_mat_triples_iter_release method.
 *
 * @param tri Triples to create a new iterator for.
 * @return New iterator
 */
gkyl_mat_triples_iter *gkyl_mat_triples_iter_new(const gkyl_mat_triples *tri);

/**
 * Initialize (set to the beginning) an existing iterator into
 * triples.
 *
 * @param iter Iterator to set back to the beginning.
 * @param tri Triples to create a new iterator for.
 */
void gkyl_mat_triples_iter_init(struct gkyl_mat_triples_iter *iter, const gkyl_mat_triples *tri);

/**
 * Move iterator to the next triple
 * 
 * @param iter Iterator to bump
 * @return True if more iterations are possible, false otherwise
 */
bool gkyl_mat_triples_iter_next(gkyl_mat_triples_iter *iter);

/**
 * Return the triple at the current location of the iterator.
 *
 * @param iter Iterator to get data from
 * @return triple at current iter location
 */
struct gkyl_mtriple gkyl_mat_triples_iter_at(const gkyl_mat_triples_iter *iter);

/**
 * Set the value of all triples to a given value.
 *
 * @param iter Iterator to release
 * @param val Value to set triples' value to.
 */
void gkyl_mat_triples_clear(struct gkyl_mat_triples *tri, double val);

/**
 * Release triples iterator
 *
 * @param iter Iterator to release
 */
void gkyl_mat_triples_iter_release(gkyl_mat_triples_iter *iter);

/**
 * Release triples
 *
 * @param tri Pointer to triples to release.
 */
void gkyl_mat_triples_release(gkyl_mat_triples *tri);
