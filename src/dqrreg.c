/* dqrreg.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int qrreg_(nx, n, px, p, nclass, x, xsc, in, y, qpivot, 
	qrank, beta, res, rss, cvar, var, varsc, scr1, work)
integer *nx, *n, *px, *p, *nclass;
doublereal *x, *xsc;
integer *in;
doublereal *y;
integer *qpivot, *qrank;
doublereal *beta, *res, *rss;
logical *cvar;
doublereal *var, *varsc, *scr1, *work;
{
    /* System generated locals */
    integer x_dim1, x_offset, xsc_dim1, xsc_offset, y_dim1, y_offset, 
	    res_dim1, res_offset, beta_dim1, beta_offset, var_dim1, 
	    var_offset, varsc_dim1, varsc_offset, i__1, i__2;

    /* Local variables */
    static integer ijob, info;
    static doublereal temp3;
    static integer i__, j, k;
    extern /* Subroutine */ int dqrsl_();
    static integer ii;
    extern /* Subroutine */ int dqrdca_();
    static integer nt;
    extern /* Subroutine */ int calcvar_();

    /* Parameter adjustments */
    --scr1;
    varsc_dim1 = *px;
    varsc_offset = varsc_dim1 + 1;
    varsc -= varsc_offset;
    var_dim1 = *px;
    var_offset = var_dim1 + 1;
    var -= var_offset;
    --qpivot;
    --in;
    xsc_dim1 = *n;
    xsc_offset = xsc_dim1 + 1;
    xsc -= xsc_offset;
    x_dim1 = *nx;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    res_dim1 = *nx;
    res_offset = res_dim1 + 1;
    res -= res_offset;
    beta_dim1 = *px;
    beta_offset = beta_dim1 + 1;
    beta -= beta_offset;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    --work;

    /* Function Body */
    ii = 0;
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	if (in[j] == 1) {
	    ++ii;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		xsc[i__ + ii * xsc_dim1] = x[i__ + j * x_dim1];
/* L23004: */
	    }
/* L23005: */
	}
/* L23000: */
    }
/* L23001: */
    nt = ii;
    ijob = 101;
    info = 1;
    temp3 = .01;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qpivot[i__] = i__;
/* L23006: */
    }
/* L23007: */
    dqrdca_(&xsc[xsc_offset], n, n, &nt, &scr1[1], &qpivot[1], &work[1], 
	    qrank, &temp3);
    *rss = (float)0.;
    i__1 = *nclass;
    for (k = 1; k <= i__1; ++k) {
	dqrsl_(&xsc[xsc_offset], n, n, qrank, &scr1[1], &y[k * y_dim1 + 1], &
		work[1], &work[1], &beta[k * beta_dim1 + 1], &work[1], &res[k 
		* res_dim1 + 1], &ijob, &info);
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    res[i__ + k * res_dim1] = y[i__ + k * y_dim1] - res[i__ + k * 
		    res_dim1];
	    *rss += res[i__ + k * res_dim1] * res[i__ + k * res_dim1];
/* L23010: */
	}
/* L23011: */
/* L23008: */
    }
/* L23009: */
    if (*cvar) {
	calcvar_(nx, n, px, &xsc[xsc_offset], qrank, &qpivot[1], &var[
		var_offset], &varsc[varsc_offset], &work[1]);
    }
    return 0;
} /* qrreg_ */

