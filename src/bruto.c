/* bruto.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static integer c__15 = 15;
static integer c__1 = 1;
static integer c__5 = 5;

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int bruto_(x, n, q, y, p, w, knot, nkmax, nk, wp, match, nef,
	 dfmax, cost, lambda, df, coef, type__, xrange, gcvsel, gcvbak, dfit, 
	maxit, nit, eta, resid, thresh, work, trace)
doublereal *x;
integer *n, *q;
doublereal *y;
integer *p;
doublereal *w, *knot;
integer *nkmax, *nk;
doublereal *wp;
integer *match, *nef;
doublereal *dfmax, *cost, *lambda, *df, *coef;
integer *type__;
doublereal *xrange, *gcvsel, *gcvbak, *dfit;
integer *maxit, *nit;
doublereal *eta, *resid, *thresh, *work;
logical *trace;
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, knot_dim1, knot_offset, 
	    coef_dim1, coef_offset, gcvsel_dim1, gcvsel_offset, gcvbak_dim1, 
	    gcvbak_offset, dfit_dim1, dfit_offset, eta_dim1, eta_offset, 
	    resid_dim1, resid_offset, match_dim1, match_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static logical select;
    extern /* Subroutine */ int bakssp_();

    /* Parameter adjustments */
    --w;
    xrange -= 3;
    --type__;
    --df;
    --lambda;
    --dfmax;
    --nef;
    match_dim1 = *n;
    match_offset = match_dim1 + 1;
    match -= match_offset;
    --nk;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    resid_dim1 = *n;
    resid_offset = resid_dim1 + 1;
    resid -= resid_offset;
    eta_dim1 = *n;
    eta_offset = eta_dim1 + 1;
    eta -= eta_offset;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    coef_dim1 = *nkmax * *p;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    knot_dim1 = *nkmax + 4;
    knot_offset = knot_dim1 + 1;
    knot -= knot_offset;
    --maxit;
    dfit_dim1 = *q;
    dfit_offset = dfit_dim1 + 1;
    dfit -= dfit_offset;
    gcvbak_dim1 = *q;
    gcvbak_offset = gcvbak_dim1 + 1;
    gcvbak -= gcvbak_offset;
    gcvsel_dim1 = *q;
    gcvsel_offset = gcvsel_dim1 + 1;
    gcvsel -= gcvsel_offset;
    --nit;
    --work;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    resid[i__ + j * resid_dim1] = y[i__ + j * y_dim1] - eta[i__ + j * 
		    eta_dim1];
/* L23002: */
	}
/* L23003: */
/* L23000: */
    }
/* L23001: */
    select = TRUE_;
    d__1 = *thresh * 10.;
    bakssp_(&select, &x[x_offset], n, q, &y[y_offset], p, &w[1], &knot[
	    knot_offset], nkmax, &nk[1], &wp[1], &match[match_offset], &nef[1]
	    , &dfmax[1], cost, &lambda[1], &df[1], &coef[coef_offset], &
	    type__[1], &xrange[3], &gcvsel[gcvsel_offset], &dfit[dfit_offset],
	     &maxit[1], &nit[1], &eta[eta_offset], &resid[resid_offset], &
	    d__1, &work[1], trace);
    select = FALSE_;
    bakssp_(&select, &x[x_offset], n, q, &y[y_offset], p, &w[1], &knot[
	    knot_offset], nkmax, &nk[1], &wp[1], &match[match_offset], &nef[1]
	    , &dfmax[1], cost, &lambda[1], &df[1], &coef[coef_offset], &
	    type__[1], &xrange[3], &gcvbak[gcvbak_offset], &dfit[dfit_offset],
	     &maxit[2], &nit[2], &eta[eta_offset], &resid[resid_offset], 
	    thresh, &work[1], trace);
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    eta[i__ + j * eta_dim1] = y[i__ + j * y_dim1] - resid[i__ + j * 
		    resid_dim1];
/* L23006: */
	}
/* L23007: */
/* L23004: */
    }
/* L23005: */
    return 0;
} /* bruto_ */

/* Subroutine */ int bakssp_(select, x, n, q, y, p, w, knot, nkmax, nk, wp, 
	match, nef, dfmax, cost, lambda, df, coef, type__, xrange, gcv, dfit, 
	maxit, nit, s, resid, thresh, work, trace)
logical *select;
doublereal *x;
integer *n, *q;
doublereal *y;
integer *p;
doublereal *w, *knot;
integer *nkmax, *nk;
doublereal *wp;
integer *match, *nef;
doublereal *dfmax, *cost, *lambda, *df, *coef;
integer *type__;
doublereal *xrange, *gcv, *dfit;
integer *maxit, *nit;
doublereal *s, *resid, *thresh, *work;
logical *trace;
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, knot_dim1, knot_offset, 
	    coef_dim1, coef_offset, gcv_dim1, gcv_offset, dfit_dim1, 
	    dfit_offset, s_dim1, s_offset, resid_dim1, resid_offset, 
	    match_dim1, match_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal sbar;
    extern /* Subroutine */ int sspl2_();
    static integer i__, j, k;
    static doublereal dfoff;
    extern doublereal wmean_();
    extern /* Subroutine */ int intpr_();
    static integer ntype;
    extern /* Subroutine */ int psspl2_();
    static doublereal ndfoff;
    extern /* Subroutine */ int dblepr_();
    static logical center;
    static doublereal gcvrat, ndf;
    static integer ier;
    static doublereal tol, rss, gcv0, gcv1;

    /* Parameter adjustments */
    --w;
    xrange -= 3;
    --type__;
    --df;
    --lambda;
    --dfmax;
    --nef;
    match_dim1 = *n;
    match_offset = match_dim1 + 1;
    match -= match_offset;
    --nk;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    resid_dim1 = *n;
    resid_offset = resid_dim1 + 1;
    resid -= resid_offset;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    coef_dim1 = *nkmax * *p;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    knot_dim1 = *nkmax + 4;
    knot_offset = knot_dim1 + 1;
    knot -= knot_offset;
    dfit_dim1 = *q;
    dfit_offset = dfit_dim1 + 1;
    dfit -= dfit_offset;
    gcv_dim1 = *q;
    gcv_offset = gcv_dim1 + 1;
    gcv -= gcv_offset;
    --work;

    /* Function Body */
    center = TRUE_;
    tol = .001;
    rss = 0.;
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	sbar = wmean_(n, &resid[j * resid_dim1 + 1], &w[1]);
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    resid[i__ + j * resid_dim1] -= sbar;
/* Computing 2nd power */
	    d__1 = resid[i__ + j * resid_dim1];
	    work[i__] = d__1 * d__1;
/* L23010: */
	}
/* L23011: */
	rss += wp[j] * wmean_(n, &work[1], &w[1]);
/* L23008: */
    }
/* L23009: */
    dfoff = 0.;
    i__1 = *q;
    for (k = 1; k <= i__1; ++k) {
	dfoff += df[k];
/* L23012: */
    }
/* L23013: */
/* Computing 2nd power */
    d__1 = 1 - (dfoff * *cost + 1) / *n;
    gcv1 = rss / (d__1 * d__1);
    gcvrat = 1.;
    *nit = 0;
L23014:
    if (*nit < *maxit && gcvrat > *thresh) {
	gcv0 = gcv1;
	++(*nit);
	i__1 = *q;
	for (k = 1; k <= i__1; ++k) {
	    gcv[k + *nit * gcv_dim1] = gcv1;
	    if (! (*select) && type__[k] == 1) {
		goto L23016;
	    }
	    if (type__[k] > 1) {
		psspl2_(&x[k * x_dim1 + 1], n, p, &knot[k * knot_dim1 + 1], &
			nk[k], &xrange[(k << 1) + 1], &coef[k * coef_dim1 + 1]
			, &coef[k * coef_dim1 + 1], &s[s_offset], &c__0, &
			type__[k]);
		i__2 = *p;
		for (j = 1; j <= i__2; ++j) {
		    sbar = wmean_(n, &s[j * s_dim1 + 1], &w[1]);
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			resid[i__ + j * resid_dim1] = resid[i__ + j * 
				resid_dim1] + s[i__ + j * s_dim1] - sbar;
/* L23024: */
		    }
/* L23025: */
/* L23022: */
		}
/* L23023: */
	    }
	    ndfoff = dfoff - df[k];
	    if (*select) {
		ntype = 0;
	    } else {
		ntype = type__[k];
	    }
	    sspl2_(&x[k * x_dim1 + 1], &resid[resid_offset], &w[1], n, p, &
		    knot[k * knot_dim1 + 1], &nk[k], &wp[1], &match[k * 
		    match_dim1 + 1], &nef[k], &ndfoff, &dfmax[k], cost, &
		    lambda[k], &ndf, &gcv1, &coef[k * coef_dim1 + 1], &s[
		    s_offset], &ntype, &center, &xrange[(k << 1) + 1], &work[
		    1], &tol, &ier);
	    gcv[k + *nit * gcv_dim1] = gcv1;
	    if (*select) {
		dfit[k + *nit * dfit_dim1] = ndf;
		df[k] = ndf;
		dfoff = ndfoff + ndf;
		type__[k] = ntype;
	    }
	    if (type__[k] > 1) {
		i__2 = *p;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			resid[i__ + j * resid_dim1] -= s[i__ + j * s_dim1];
/* L23034: */
		    }
/* L23035: */
/* L23032: */
		}
/* L23033: */
	    }
L23016:
	    ;
	}
/* L23017: */
	gcvrat = (d__1 = gcv1 - gcv0, abs(d__1)) / gcv0;
	if (*trace) {
	    intpr_("outer iteration", &c__15, nit, &c__1, 15L);
	    dblepr_("gcv  ", &c__5, &gcv1, &c__1, 5L);
	    dblepr_("ratio", &c__5, &gcvrat, &c__1, 5L);
	}
	goto L23014;
    }
/* L23015: */
    return 0;
} /* bakssp_ */

/* Subroutine */ int pbruto_(x, n, q, ybar, p, knot, nkmax, nk, coef, type__, 
	xrange, eta, work)
doublereal *x;
integer *n, *q;
doublereal *ybar;
integer *p;
doublereal *knot;
integer *nkmax, *nk;
doublereal *coef;
integer *type__;
doublereal *xrange, *eta, *work;
{
    /* System generated locals */
    integer x_dim1, x_offset, knot_dim1, knot_offset, coef_dim1, coef_offset, 
	    eta_dim1, eta_offset, work_dim1, work_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int psspl2_();

    /* Parameter adjustments */
    xrange -= 3;
    --type__;
    --nk;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    work_dim1 = *n;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    eta_dim1 = *n;
    eta_offset = eta_dim1 + 1;
    eta -= eta_offset;
    --ybar;
    coef_dim1 = *nkmax * *p;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    knot_dim1 = *nkmax + 4;
    knot_offset = knot_dim1 + 1;
    knot -= knot_offset;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    eta[i__ + j * eta_dim1] = ybar[j];
/* L23040: */
	}
/* L23041: */
/* L23038: */
    }
/* L23039: */
    i__1 = *q;
    for (k = 1; k <= i__1; ++k) {
	if (type__[k] == 1) {
	    goto L23042;
	}
	psspl2_(&x[k * x_dim1 + 1], n, p, &knot[k * knot_dim1 + 1], &nk[k], &
		xrange[(k << 1) + 1], &coef[k * coef_dim1 + 1], &coef[k * 
		coef_dim1 + 1], &work[work_offset], &c__0, &type__[k]);
	i__2 = *p;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		eta[i__ + j * eta_dim1] += work[i__ + j * work_dim1];
/* L23048: */
	    }
/* L23049: */
/* L23046: */
	}
/* L23047: */
L23042:
	;
    }
/* L23043: */
    return 0;
} /* pbruto_ */

