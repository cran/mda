/* dcalcvar.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#define dbksl_ dtrsl_

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int calcvar_(nx, n, px, qr, qrank, qpivot, cov, tmpcov, work)
integer *nx, *n, *px;
doublereal *qr;
integer *qrank, *qpivot;
doublereal *cov, *tmpcov, *work;
{
    /* System generated locals */
    integer qr_dim1, qr_offset, cov_dim1, cov_offset, tmpcov_dim1, 
	    tmpcov_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal dsum;
    static integer i__, j, k;
    extern /* Subroutine */ int dbksl_();
    static integer km;

    /* Parameter adjustments */
    tmpcov_dim1 = *px;
    tmpcov_offset = tmpcov_dim1 + 1;
    tmpcov -= tmpcov_offset;
    cov_dim1 = *px;
    cov_offset = cov_dim1 + 1;
    cov -= cov_offset;
    --qpivot;
    qr_dim1 = *nx;
    qr_offset = qr_dim1 + 1;
    qr -= qr_offset;
    --work;

    /* Function Body */
    i__1 = *qrank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *qrank;
	for (j = 1; j <= i__2; ++j) {
	    tmpcov[i__ + j * tmpcov_dim1] = 0.;
	    cov[i__ + j * cov_dim1] = qr[i__ + j * qr_dim1];
/* L23002: */
	}
/* L23003: */
	tmpcov[i__ + i__ * tmpcov_dim1] = (float)1.;
/* L23000: */
    }
/* L23001: */
    info = 0;
    dbksl_(&cov[cov_offset], px, qrank, &tmpcov[tmpcov_offset], px, &info);
    i__1 = *qrank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *qrank;
	for (j = i__; j <= i__2; ++j) {
	    dsum = (float)0.;
	    km = max(i__,j);
	    k = km;
L23008:
	    if (! (k <= *qrank)) {
		goto L23010;
	    }
	    dsum += tmpcov[i__ + k * tmpcov_dim1] * tmpcov[j + k * 
		    tmpcov_dim1];
/* L23009: */
	    ++k;
	    goto L23008;
L23010:
	    tmpcov[i__ + j * tmpcov_dim1] = dsum;
	    tmpcov[j + i__ * tmpcov_dim1] = dsum;
/* L23006: */
	}
/* L23007: */
/* L23004: */
    }
/* L23005: */
    i__1 = *qrank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *qrank;
	for (j = 1; j <= i__2; ++j) {
	    cov[i__ + j * cov_dim1] = tmpcov[i__ + j * tmpcov_dim1];
/* L23013: */
	}
/* L23014: */
/* L23011: */
    }
/* L23012: */
    return 0;
} /* calcvar_ */

