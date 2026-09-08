#ifndef __dlapack_lite__
#define __dlapack_lite__

#include <stdio.h>

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif 
/* dlapack_lite.c */
    int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *iwork, integer *info);
    int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *ldc, doublereal *work, integer *info);
    int dgebak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *ldv, integer *info);
    int dgebal_(char *job, integer *n, doublereal *a, integer *lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);
    int dgebd2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *work, integer *info);
    int dgebrd_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *work, integer *lwork, integer *info);
    int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, integer *info);
    int dgehd2_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);
    int dgehrd_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dgelq2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);
    int dgelqf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dgelsd_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork, integer *iwork, integer *info);
    int dgeqr2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);
    int dgeqrf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *iwork, integer *info);
    int dgesv_(integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
    int dgetf2_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
    int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
    int dgetrs_(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
    int dhseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *info);
    int dlabad_(doublereal *small, doublereal *large);
    int dlabrd_(integer *m, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer *ldy);
    int dlacpy_(char *uplo, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb);
    int dladiv_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *p, doublereal *q);
    int dlae2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *rt1, doublereal *rt2);
    int dlaed0_(integer *icompq, integer *qsiz, integer *n, doublereal *d__, doublereal *e, doublereal *q, integer *ldq, doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, integer *info);
    int dlaed1_(integer *n, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, doublereal *work, integer *iwork, integer *info);
    int dlaed2_(integer *k, integer *n, integer *n1, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2, integer *indx, integer *indxc, integer *indxp, integer *coltyp, integer *info);
    int dlaed3_(integer *k, integer *n, integer *n1, doublereal *d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *q2, integer *indx, integer *ctot, doublereal *w, doublereal *s, integer *info);
    int dlaed4_(integer *n, integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam, integer *info);
    int dlaed5_(integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam);
    int dlaed6_(integer *kniter, logical *orgati, doublereal *rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *tau, integer *info);
    int dlaed7_(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *work, integer *iwork, integer *info);
    int dlaed8_(integer *icompq, integer *k, integer *n, integer *qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda, doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, integer *indxp, integer *indx, integer *info);
    int dlaed9_(integer *k, integer *kstart, integer *kstop, integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, integer *info);
    int dlaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, doublereal *z__, doublereal *ztemp, integer *info);
    int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);
    int dlahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, integer *ldz, integer *info);
    int dlahrd_(integer *n, integer *k, integer *nb, doublereal *a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, doublereal *y, integer *ldy);
    int dlaln2_(logical *ltrans, integer *na, integer *nw, doublereal *smin, doublereal *ca, doublereal *a, integer *lda, doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, doublereal *scale, doublereal *xnorm, integer *info);
    int dlals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, doublereal *work, integer *info);
    int dlalsa_(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *work, integer *iwork, integer *info);
    int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb, doublereal *rcond, integer *rank, doublereal *work, integer *iwork, integer *info);
    doublereal dlamch_(char *cmach);
    int dlamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
    int dlamc2_(integer *beta, integer *t, logical *rnd, doublereal *eps, integer *emin, doublereal *rmin, integer *emax, doublereal *rmax);
    doublereal dlamc3_(doublereal *a, doublereal *b);
    int dlamc4_(integer *emin, doublereal *start, integer *base);
    int dlamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, doublereal *rmax);
    int dlamrg_(integer *n1, integer *n2, doublereal *a, integer *dtrd1, integer *dtrd2, integer *index);
    doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer *lda, doublereal *work);
    doublereal dlanhs_(char *norm, integer *n, doublereal *a, integer *lda, doublereal *work);
    doublereal dlanst_(char *norm, integer *n, doublereal *d__, doublereal *e);
    doublereal dlansy_(char *norm, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *work);
    int dlanv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r, doublereal *rt2i, doublereal *cs, doublereal *sn);
    doublereal dlapy2_(doublereal *x, doublereal *y);
    doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z__);
    int dlarf_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);
    int dlarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k, doublereal *v, integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, doublereal *work, integer *ldwork);
    int dlarfg_(integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *tau);
    int dlarft_(char *direct, char *storev, integer *n, integer *k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, integer *ldt);
    int dlarfx_(char *side, integer *m, integer *n, doublereal *v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);
    int dlartg_(doublereal *f, doublereal *g, doublereal *cs, doublereal *sn, doublereal *r__);
    int dlas2_(doublereal *f, doublereal *g, doublereal *h__, doublereal *ssmin, doublereal *ssmax);
    int dlascl_(char *type__, integer *kl, integer *ku, doublereal *cfrom, doublereal *cto, integer *m, integer *n, doublereal *a, integer *lda, integer *info);
    int dlasd0_(integer *n, integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *info);
    int dlasd1_(integer *nl, integer *nr, integer *sqre, doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *iwork, doublereal *work, integer *info);
    int dlasd2_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info);
    int dlasd3_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma, doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, integer *idxc, integer *ctot, doublereal *z__, integer *info);
    int dlasd4_(integer *n, integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho, doublereal *sigma, doublereal *work, integer *info);
    int dlasd5_(integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *work);
    int dlasd6_(integer *icompq, integer *nl, integer *nr, integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl, doublereal *alpha, doublereal *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, doublereal *work, integer *iwork, integer *info);
    int dlasd7_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__, doublereal *z__, doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl, doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *c__, doublereal *s, integer *info);
    int dlasd8_(integer *icompq, integer *k, doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *work, integer *info);
    int dlasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *work, integer *iwork, integer *info);
    int dlasdq_(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *ldc, doublereal *work, integer *info);
    int dlasdt_(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml, integer *ndimr, integer *msub);
    int dlaset_(char *uplo, integer *m, integer *n, doublereal *alpha, doublereal *beta, doublereal *a, integer *lda);
    int dlasq1_(integer *n, doublereal *d__, doublereal *e, doublereal *work, integer *info);
    int dlasq2_(integer *n, doublereal *z__, integer *info);
    int dlasq3_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig, doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee);
    int dlasq4_(integer *i0, integer *n0, doublereal *z__, integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *tau, integer *ttype);
    int dlasq5_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2, logical *ieee);
    int dlasq6_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2);
    int dlasr_(char *side, char *pivot, char *direct, integer *m, integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *lda);
    int dlasrt_(char *id, integer *n, doublereal *d__, integer *info);
    int dlassq_(integer *n, doublereal *x, integer *incx, doublereal *scale, doublereal *sumsq);
    int dlasv2_(doublereal *f, doublereal *g, doublereal *h__, doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *csr, doublereal *snl, doublereal *csl);
    int dlaswp_(integer *n, doublereal *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx);
    int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, integer *ldw);
    int dorg2r_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);
    int dorgbr_(char *vect, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dorghr_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dorgl2_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *info);
    int dorglq_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dorgqr_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dorm2l_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *info);
    int dorm2r_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *info);
    int dormbr_(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
    int dorml2_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *info);
    int dormlq_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
    int dormql_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
    int dormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
    int dormtr_(char *side, char *uplo, char *trans, integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
    int dstedc_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
    int dsteqr_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz, doublereal *work, integer *info);
    int dsterf_(integer *n, doublereal *d__, doublereal *e, integer *info);
    int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
    int dsytd2_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);
    int dsytrd_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *work, integer *lwork, integer *info);
    int dtrevc_(char *side, char *howmny, logical *select, integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, doublereal *work, integer *info);
    integer ieeeck_(integer *ispec, real *zero, real *one);
    integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len);
    int dpotf2_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info);
    int dpotrf_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info);
    int dtrtri_(char *uplo, char *diag, integer *n, doublereal *a, integer *lda, integer *info);
    int dgbsv_(integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, integer *ldb, integer *info);
    int dgtsv_(integer *n, integer *nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer *info);

#ifdef __cplusplus
}
#endif

#endif // __dlapack_lite__
