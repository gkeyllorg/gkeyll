#ifndef __blas_lite__
#define __blas_lite__

#include <stdio.h>

#include "f2c.h"

// This file was generate automatically using cproto. It may hence look rather
// hairy
#ifdef __cplusplus
extern "C" {
#endif 
    int daxpy_(integer *n, doublereal *da, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
    doublereal dcabs1_(doublecomplex *z__);
    int dcopy_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
    doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
    int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);
    int dgemv_(char *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);
    int dger_(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda);
    doublereal dnrm2_(integer *n, doublereal *x, integer *incx);
    int drot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s);
    int dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx);
    int dswap_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy);
    int dsymv_(char *uplo, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy);
    int dsyr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, integer *incy, doublereal *a, integer *lda);
    int dsyr2k_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);
    int dtrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb);
    int dtrmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, integer *lda, doublereal *x, integer *incx);
    int dtrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb);
    doublereal dzasum_(integer *n, doublecomplex *zx, integer *incx);
    doublereal dznrm2_(integer *n, doublecomplex *x, integer *incx);
    integer idamax_(integer *n, doublereal *dx, integer *incx);
    integer izamax_(integer *n, doublecomplex *zx, integer *incx);
    int xerbla_(char *srname, integer *info);
    int zaxpy_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
    int zcopy_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
    void zdotc_(doublecomplex *ret_val, integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
    void zdotu_(doublecomplex *ret_val, integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
    int zdscal_(integer *n, doublereal *da, doublecomplex *zx, integer *incx);
    int zgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc);
    int zgemv_(char *trans, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);
    int zgerc_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);
    int zgeru_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);
    int zhemv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy);
    int zher2_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublecomplex *a, integer *lda);
    int zher2k_(char *uplo, char *trans, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *beta, doublecomplex *c__, integer *ldc);
    int zscal_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx);
    int zswap_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy);
    int ztrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);
    int ztrmv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx);
    int ztrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);
    int ztrsv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx);
    logical lsame_(char *ca, char *cb);
    int dsyrk_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, doublereal *c__, integer *ldc);
    int zherk_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublecomplex *a, integer *lda, doublereal *beta, doublecomplex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif

#endif // __blas_lite__
