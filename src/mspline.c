/* mspline.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b9 = 16.;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__0 = 0;
static doublereal c_b167 = 0.;
static doublereal c_b170 = 1.;

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int sspl_(x, y, w, n, ldy, p, knot, nk, method, tol, wp, ssy,
	 dfoff, dfmax, cost, lambda, df, cv, gcv, coef, s, lev, xwy, hs, sg, 
	abd, p1ip, ier)
doublereal *x, *y, *w;
integer *n, *ldy, *p;
doublereal *knot;
integer *nk, *method;
doublereal *tol, *wp, *ssy, *dfoff, *dfmax, *cost, *lambda, *df, *cv, *gcv, *
	coef, *s, *lev, *xwy, *hs, *sg, *abd, *p1ip;
integer *ier;
{
    /* System generated locals */
    integer y_dim1, y_offset, coef_dim1, coef_offset, s_dim1, s_offset, 
	    xwy_dim1, xwy_offset, hs_dim1, hs_offset, sg_dim1, sg_offset;

    /* Local variables */
    extern /* Subroutine */ int sgram_(), sslvr2_(), stxwx2_(), fmm_();

    /* Parameter adjustments */
    --lev;
    --w;
    --x;
    s_dim1 = *ldy;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --ssy;
    --wp;
    y_dim1 = *ldy;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    p1ip -= 5;
    abd -= 5;
    sg_dim1 = *nk;
    sg_offset = sg_dim1 + 1;
    sg -= sg_offset;
    hs_dim1 = *nk;
    hs_offset = hs_dim1 + 1;
    hs -= hs_offset;
    xwy_dim1 = *nk;
    xwy_offset = xwy_dim1 + 1;
    xwy -= xwy_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;

    /* Function Body */
    sgram_(&sg[sg_dim1 + 1], &sg[(sg_dim1 << 1) + 1], &sg[sg_dim1 * 3 + 1], &
	    sg[(sg_dim1 << 2) + 1], &knot[1], nk);
    stxwx2_(&x[1], &y[y_offset], &w[1], n, ldy, p, &knot[1], nk, &xwy[
	    xwy_offset], &hs[hs_dim1 + 1], &hs[(hs_dim1 << 1) + 1], &hs[
	    hs_dim1 * 3 + 1], &hs[(hs_dim1 << 2) + 1]);
    if (*method == 1) {
	sslvr2_(&x[1], &y[y_offset], &w[1], n, ldy, p, &knot[1], nk, method, 
		tol, &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, &coef[
		coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset], &hs[
		hs_dim1 + 1], &hs[(hs_dim1 << 1) + 1], &hs[hs_dim1 * 3 + 1], &
		hs[(hs_dim1 << 2) + 1], &sg[sg_dim1 + 1], &sg[(sg_dim1 << 1) 
		+ 1], &sg[sg_dim1 * 3 + 1], &sg[(sg_dim1 << 2) + 1], &abd[5], 
		&p1ip[5], ier);
    } else {
	fmm_(&x[1], &y[y_offset], &w[1], n, ldy, p, &knot[1], nk, method, tol,
		 &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, &coef[
		coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset], &hs[
		hs_offset], &sg[sg_offset], &abd[5], &p1ip[5], ier);
	if (*method > 2 && *df > *dfmax) {
	    *df = *dfmax;
	    fmm_(&x[1], &y[y_offset], &w[1], n, ldy, p, &knot[1], nk, &c__2, 
		    tol, &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, &
		    coef[coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset]
		    , &hs[hs_offset], &sg[sg_offset], &abd[5], &p1ip[5], ier);
	}
    }
    return 0;
} /* sspl_ */

/* Subroutine */ int fmm_(xs, ys, ws, n, ldy, nvar, knot, nk, method, tol, wp,
	 ssy, dfoff, cost, lambda, df, cv, gcv, coef, s, lev, xwy, hs, sg, 
	abd, p1ip, ier)
doublereal *xs, *ys, *ws;
integer *n, *ldy, *nvar;
doublereal *knot;
integer *nk, *method;
doublereal *tol, *wp, *ssy, *dfoff, *cost, *lambda, *df, *cv, *gcv, *coef, *s,
	 *lev, *xwy, *hs, *sg, *abd, *p1ip;
integer *ier;
{
    /* System generated locals */
    integer ys_dim1, ys_offset, coef_dim1, coef_offset, s_dim1, s_offset, 
	    xwy_dim1, xwy_offset, hs_dim1, hs_offset, sg_dim1, sg_offset, 
	    i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(), pow_dd(), d_sign();

    /* Local variables */
    static doublereal a, b, c__, d__, e;
    static integer i__;
    static doublereal p, q, r__, u, v, w, x, ratio, t1, t2, fu, fv, fw, fx, 
	    ax, bx;
    extern /* Subroutine */ int sslvr2_();
    static doublereal xm, targdf, eps;
    static integer i23010, i23039;
    static doublereal tol1, tol2;

    /* Parameter adjustments */
    --lev;
    --ws;
    --xs;
    s_dim1 = *ldy;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --ssy;
    --wp;
    ys_dim1 = *ldy;
    ys_offset = ys_dim1 + 1;
    ys -= ys_offset;
    p1ip -= 5;
    abd -= 5;
    sg_dim1 = *nk;
    sg_offset = sg_dim1 + 1;
    sg -= sg_offset;
    hs_dim1 = *nk;
    hs_offset = hs_dim1 + 1;
    hs -= hs_offset;
    xwy_dim1 = *nk;
    xwy_offset = xwy_dim1 + 1;
    xwy -= xwy_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;

    /* Function Body */
    ax = 1e-10;
    bx = (float)1.5;
    t1 = (float)0.;
    t2 = (float)0.;
    targdf = *df;
    i__1 = *nk - 3;
    for (i__ = 3; i__ <= i__1; ++i__) {
	t1 += hs[i__ + hs_dim1];
/* L23004: */
    }
/* L23005: */
    i__1 = *nk - 3;
    for (i__ = 3; i__ <= i__1; ++i__) {
	t2 += sg[i__ + sg_dim1];
/* L23006: */
    }
/* L23007: */
    ratio = t1 / t2;
    c__ = ((float)3. - sqrt(5.)) * (float).5;
    eps = 1.;
L10:
    eps /= 2.;
    tol1 = eps + 1.;
    if (tol1 > 1.) {
	goto L10;
    }
    eps = sqrt(eps);
    a = ax;
    b = bx;
    v = a + c__ * (b - a);
    w = v;
    x = v;
    e = (float)0.;
    d__1 = x * (float)6. - (float)2.;
    *lambda = ratio * pow_dd(&c_b9, &d__1);
    sslvr2_(&xs[1], &ys[ys_offset], &ws[1], n, ldy, nvar, &knot[1], nk, 
	    method, tol, &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, &
	    coef[coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset], &hs[
	    hs_dim1 + 1], &hs[(hs_dim1 << 1) + 1], &hs[hs_dim1 * 3 + 1], &hs[(
	    hs_dim1 << 2) + 1], &sg[sg_dim1 + 1], &sg[(sg_dim1 << 1) + 1], &
	    sg[sg_dim1 * 3 + 1], &sg[(sg_dim1 << 2) + 1], &abd[5], &p1ip[5], 
	    ier);
    i23010 = *method;
    goto L23010;
L23012:
/* Computing 2nd power */
    d__1 = targdf - *df;
    fx = d__1 * d__1 + 3.;
    goto L23011;
L23013:
    fx = *gcv;
    goto L23011;
L23014:
    fx = *cv;
    goto L23011;
L23010:
    if (i23010 == 2) {
	goto L23012;
    }
    if (i23010 == 3) {
	goto L23013;
    }
    if (i23010 == 4) {
	goto L23014;
    }
L23011:
    fv = fx;
    fw = fx;
L20:
    xm = (a + b) * (float).5;
    tol1 = eps * abs(x) + *tol / 3.;
    tol2 = tol1 * 2.;
    if ((d__1 = x - xm, abs(d__1)) <= tol2 - (b - a) * (float).5) {
	goto L90;
    }
    if (abs(e) <= tol1) {
	goto L40;
    }
    r__ = (x - w) * (fx - fv);
    q = (x - v) * (fx - fw);
    p = (x - v) * q - (x - w) * r__;
    q = (q - r__) * (float)2.;
    if (q > (float)0.) {
	p = -p;
    }
    q = abs(q);
    r__ = e;
    e = d__;
/* L30: */
    if (abs(p) >= (d__1 = q * (float).5 * r__, abs(d__1))) {
	goto L40;
    }
    if (p <= q * (a - x)) {
	goto L40;
    }
    if (p >= q * (b - x)) {
	goto L40;
    }
    d__ = p / q;
    u = x + d__;
    if (u - a < tol2) {
	d__1 = xm - x;
	d__ = d_sign(&tol1, &d__1);
    }
    if (b - u < tol2) {
	d__1 = xm - x;
	d__ = d_sign(&tol1, &d__1);
    }
    goto L50;
L40:
    if (x >= xm) {
	e = a - x;
    }
    if (x < xm) {
	e = b - x;
    }
    d__ = c__ * e;
L50:
    if (abs(d__) >= tol1) {
	u = x + d__;
    }
    if (abs(d__) < tol1) {
	u = x + d_sign(&tol1, &d__);
    }
    d__1 = u * (float)6. - (float)2.;
    *lambda = ratio * pow_dd(&c_b9, &d__1);
    sslvr2_(&xs[1], &ys[ys_offset], &ws[1], n, ldy, nvar, &knot[1], nk, 
	    method, tol, &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, &
	    coef[coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset], &hs[
	    hs_dim1 + 1], &hs[(hs_dim1 << 1) + 1], &hs[hs_dim1 * 3 + 1], &hs[(
	    hs_dim1 << 2) + 1], &sg[sg_dim1 + 1], &sg[(sg_dim1 << 1) + 1], &
	    sg[sg_dim1 * 3 + 1], &sg[(sg_dim1 << 2) + 1], &abd[5], &p1ip[5], 
	    ier);
    i23039 = *method;
    goto L23039;
L23041:
/* Computing 2nd power */
    d__1 = targdf - *df;
    fu = d__1 * d__1 + 3.;
    goto L23040;
L23042:
    fu = *gcv;
    goto L23040;
L23043:
    fu = *cv;
    goto L23040;
L23039:
    if (i23039 == 2) {
	goto L23041;
    }
    if (i23039 == 3) {
	goto L23042;
    }
    if (i23039 == 4) {
	goto L23043;
    }
L23040:
    if (fu > fx) {
	goto L60;
    }
    if (u >= x) {
	a = x;
    }
    if (u < x) {
	b = x;
    }
    v = w;
    fv = fw;
    w = x;
    fw = fx;
    x = u;
    fx = fu;
    goto L20;
L60:
    if (u < x) {
	a = u;
    }
    if (u >= x) {
	b = u;
    }
    if (fu <= fw) {
	goto L70;
    }
    if (w == x) {
	goto L70;
    }
    if (fu <= fv) {
	goto L80;
    }
    if (v == x) {
	goto L80;
    }
    if (v == w) {
	goto L80;
    }
    goto L20;
L70:
    v = w;
    fv = fw;
    w = u;
    fw = fu;
    goto L20;
L80:
    v = u;
    fv = fu;
    goto L20;
L90:
    if (*method == 2) {
	sslvr2_(&xs[1], &ys[ys_offset], &ws[1], n, ldy, nvar, &knot[1], nk, &
		c__1, tol, &wp[1], &ssy[1], dfoff, cost, lambda, df, cv, gcv, 
		&coef[coef_offset], &s[s_offset], &lev[1], &xwy[xwy_offset], &
		hs[hs_dim1 + 1], &hs[(hs_dim1 << 1) + 1], &hs[hs_dim1 * 3 + 1]
		, &hs[(hs_dim1 << 2) + 1], &sg[sg_dim1 + 1], &sg[(sg_dim1 << 
		1) + 1], &sg[sg_dim1 * 3 + 1], &sg[(sg_dim1 << 2) + 1], &abd[
		5], &p1ip[5], ier);
    }
    return 0;
} /* fmm_ */

/* Subroutine */ int stxwx2_(x, z__, w, k, ldy, pz, xknot, n, y, hs0, hs1, 
	hs2, hs3)
doublereal *x, *z__, *w;
integer *k, *ldy, *pz;
doublereal *xknot;
integer *n;
doublereal *y, *hs0, *hs1, *hs2, *hs3;
{
    /* System generated locals */
    integer z_dim1, z_offset, y_dim1, y_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal work[16];
    static integer i__, j, mflag, ileft;
    static doublereal vnikx[4]	/* was [4][1] */;
    static integer pp;
    extern /* Subroutine */ int bsplvd_(), interv_();
    static doublereal eps;

    /* Parameter adjustments */
    --w;
    --x;
    z_dim1 = *ldy;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --hs3;
    --hs2;
    --hs1;
    --hs0;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    --xknot;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hs0[i__] = 0.;
	hs1[i__] = 0.;
	hs2[i__] = 0.;
	hs3[i__] = 0.;
	i__2 = *pz;
	for (j = 1; j <= i__2; ++j) {
	    y[i__ + j * y_dim1] = 0.;
/* L23068: */
	}
/* L23069: */
/* L23066: */
    }
/* L23067: */
    eps = 1e-10;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n + 1;
	interv_(&xknot[1], &i__2, &x[i__], &ileft, &mflag);
	if (mflag == 1) {
	    if (x[i__] <= xknot[ileft] + eps) {
		--ileft;
	    } else {
		return 0;
	    }
	}
	bsplvd_(&xknot[1], &c__4, &x[i__], &ileft, work, vnikx, &c__1);
	j = ileft - 3;
	i__2 = *pz;
	for (pp = 1; pp <= i__2; ++pp) {
	    y[j + pp * y_dim1] += w[i__] * z__[i__ + pp * z_dim1] * vnikx[0];
/* L23076: */
	}
/* L23077: */
/* Computing 2nd power */
	d__1 = vnikx[0];
	hs0[j] += w[i__] * (d__1 * d__1);
	hs1[j] += w[i__] * vnikx[0] * vnikx[1];
	hs2[j] += w[i__] * vnikx[0] * vnikx[2];
	hs3[j] += w[i__] * vnikx[0] * vnikx[3];
	j = ileft - 2;
	i__2 = *pz;
	for (pp = 1; pp <= i__2; ++pp) {
	    y[j + pp * y_dim1] += w[i__] * z__[i__ + pp * z_dim1] * vnikx[1];
/* L23078: */
	}
/* L23079: */
/* Computing 2nd power */
	d__1 = vnikx[1];
	hs0[j] += w[i__] * (d__1 * d__1);
	hs1[j] += w[i__] * vnikx[1] * vnikx[2];
	hs2[j] += w[i__] * vnikx[1] * vnikx[3];
	j = ileft - 1;
	i__2 = *pz;
	for (pp = 1; pp <= i__2; ++pp) {
	    y[j + pp * y_dim1] += w[i__] * z__[i__ + pp * z_dim1] * vnikx[2];
/* L23080: */
	}
/* L23081: */
/* Computing 2nd power */
	d__1 = vnikx[2];
	hs0[j] += w[i__] * (d__1 * d__1);
	hs1[j] += w[i__] * vnikx[2] * vnikx[3];
	j = ileft;
	i__2 = *pz;
	for (pp = 1; pp <= i__2; ++pp) {
	    y[j + pp * y_dim1] += w[i__] * z__[i__ + pp * z_dim1] * vnikx[3];
/* L23082: */
	}
/* L23083: */
/* Computing 2nd power */
	d__1 = vnikx[3];
	hs0[j] += w[i__] * (d__1 * d__1);
/* L23070: */
    }
/* L23071: */
    return 0;
} /* stxwx2_ */

/* Subroutine */ int sslvr2_(x, y, w, n, ldy, p, knot, nk, method, tol, wp, 
	ssy, dfoff, cost, lambda, df, cv, gcv, coef, sz, lev, xwy, hs0, hs1, 
	hs2, hs3, sg0, sg1, sg2, sg3, abd, p1ip, info)
doublereal *x, *y, *w;
integer *n, *ldy, *p;
doublereal *knot;
integer *nk, *method;
doublereal *tol, *wp, *ssy, *dfoff, *cost, *lambda, *df, *cv, *gcv, *coef, *
	sz, *lev, *xwy, *hs0, *hs1, *hs2, *hs3, *sg0, *sg1, *sg2, *sg3, *abd, 
	*p1ip;
integer *info;
{
    /* System generated locals */
    integer y_dim1, y_offset, coef_dim1, coef_offset, sz_dim1, sz_offset, 
	    xwy_dim1, xwy_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal work[16], sumw, tssy;
    static integer i__, j;
    extern /* Subroutine */ int dpbfa_();
    static integer icoef, mflag, ileft;
    extern /* Subroutine */ int dpbsl_();
    static doublereal b0, b1, b2, b3, vnikx[4]	/* was [4][1] */;
    extern /* Subroutine */ int sinrp2_();
    static doublereal xv;
    extern doublereal bvalue_();
    extern /* Subroutine */ int bsplvd_();
    static integer ld4;
    static logical fittoo;
    extern /* Subroutine */ int interv_();
    static integer ilo;
    static doublereal eps, rss;

    /* Parameter adjustments */
    --lev;
    --w;
    --x;
    sz_dim1 = *ldy;
    sz_offset = sz_dim1 + 1;
    sz -= sz_offset;
    --ssy;
    --wp;
    y_dim1 = *ldy;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    p1ip -= 5;
    abd -= 5;
    --sg3;
    --sg2;
    --sg1;
    --sg0;
    --hs3;
    --hs2;
    --hs1;
    --hs0;
    xwy_dim1 = *nk;
    xwy_offset = xwy_dim1 + 1;
    xwy -= xwy_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;

    /* Function Body */
    fittoo = *method != 2;
    ilo = 1;
    eps = 1e-11;
    ld4 = 4;
    if (fittoo) {
	i__1 = *nk;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *p;
	    for (j = 1; j <= i__2; ++j) {
		coef[i__ + j * coef_dim1] = xwy[i__ + j * xwy_dim1];
/* L23088: */
	    }
/* L23089: */
/* L23086: */
	}
/* L23087: */
    }
    i__1 = *nk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ << 2) + 4] = hs0[i__] + *lambda * sg0[i__];
/* L23090: */
    }
/* L23091: */
    i__1 = *nk - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 1 << 2) + 3] = hs1[i__] + *lambda * sg1[i__];
/* L23092: */
    }
/* L23093: */
    i__1 = *nk - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 2 << 2) + 2] = hs2[i__] + *lambda * sg2[i__];
/* L23094: */
    }
/* L23095: */
    i__1 = *nk - 3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 3 << 2) + 1] = hs3[i__] + *lambda * sg3[i__];
/* L23096: */
    }
/* L23097: */
    dpbfa_(&abd[5], &ld4, nk, &c__3, info);
    if (*info != 0) {
	return 0;
    }
    if (fittoo) {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    dpbsl_(&abd[5], &ld4, nk, &c__3, &coef[j * coef_dim1 + 1]);
/* L23102: */
	}
/* L23103: */
	icoef = 1;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xv = x[i__];
	    i__2 = *p;
	    for (j = 1; j <= i__2; ++j) {
		sz[i__ + j * sz_dim1] = bvalue_(&knot[1], &coef[j * coef_dim1 
			+ 1], nk, &c__4, &xv, &c__0);
/* L23106: */
	    }
/* L23107: */
/* L23104: */
	}
/* L23105: */
    }
    sinrp2_(&abd[5], &ld4, nk, &p1ip[5]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xv = x[i__];
	i__2 = *nk + 1;
	interv_(&knot[1], &i__2, &xv, &ileft, &mflag);
	if (mflag == -1) {
	    ileft = 4;
	    xv = knot[4] + eps;
	}
	if (mflag == 1) {
	    ileft = *nk;
	    xv = knot[*nk + 1] - eps;
	}
	j = ileft - 3;
	bsplvd_(&knot[1], &c__4, &xv, &ileft, work, vnikx, &c__1);
	b0 = vnikx[0];
	b1 = vnikx[1];
	b2 = vnikx[2];
	b3 = vnikx[3];
/* Computing 2nd power */
	d__1 = b0;
/* Computing 2nd power */
	d__2 = b1;
/* Computing 2nd power */
	d__3 = b2;
/* Computing 2nd power */
	d__4 = b3;
	lev[i__] = (p1ip[(j << 2) + 4] * (d__1 * d__1) + p1ip[(j << 2) + 3] * 
		(float)2. * b0 * b1 + p1ip[(j << 2) + 2] * (float)2. * b0 * 
		b2 + p1ip[(j << 2) + 1] * (float)2. * b0 * b3 + p1ip[(j + 1 <<
		 2) + 4] * (d__2 * d__2) + p1ip[(j + 1 << 2) + 3] * (float)2. 
		* b1 * b2 + p1ip[(j + 1 << 2) + 2] * (float)2. * b1 * b3 + 
		p1ip[(j + 2 << 2) + 4] * (d__3 * d__3) + p1ip[(j + 2 << 2) + 
		3] * (float)2. * b2 * b3 + p1ip[(j + 3 << 2) + 4] * (d__4 * 
		d__4)) * w[i__];
/* L23108: */
    }
/* L23109: */
    rss = 0.;
    *df = 0.;
    sumw = 0.;
    *gcv = 0.;
    *cv = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*df += lev[i__];
/* L23114: */
    }
/* L23115: */
    if (fittoo) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sumw += w[i__];
	    i__2 = *p;
	    for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
		d__1 = y[i__ + j * y_dim1] - sz[i__ + j * sz_dim1];
		rss += w[i__] * wp[j] * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = (y[i__ + j * y_dim1] - sz[i__ + j * sz_dim1]) / (1 - 
			lev[i__]);
		*cv += w[i__] * wp[j] * (d__1 * d__1);
/* L23120: */
	    }
/* L23121: */
/* L23118: */
	}
/* L23119: */
	tssy = 0.;
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    tssy += wp[j] * ssy[j];
/* L23122: */
	}
/* L23123: */
/* Computing 2nd power */
	d__1 = 1. - ((*dfoff + *df - 1) * *cost + 1) / sumw;
	*gcv = (rss + tssy) / sumw / (d__1 * d__1);
	*cv = (*cv + tssy) / sumw;
    }
    return 0;
} /* sslvr2_ */

/* Subroutine */ int sinrp2_(abd, ld4, nk, p1ip)
doublereal *abd;
integer *ld4, *nk;
doublereal *p1ip;
{
    /* System generated locals */
    integer abd_dim1, abd_offset, p1ip_dim1, p1ip_offset, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j;
    static doublereal c0, c1, c2, c3, wjm1[1], wjm2[2], wjm3[3];

    /* Parameter adjustments */
    p1ip_dim1 = *ld4;
    p1ip_offset = p1ip_dim1 + 1;
    p1ip -= p1ip_offset;
    abd_dim1 = *ld4;
    abd_offset = abd_dim1 + 1;
    abd -= abd_offset;

    /* Function Body */
    wjm3[0] = 0.;
    wjm3[1] = 0.;
    wjm3[0] = 0.;
    wjm2[0] = 0.;
    wjm2[1] = 0.;
    wjm1[0] = 0.;
    i__1 = *nk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *nk - i__ + 1;
	c0 = 1. / abd[j * abd_dim1 + 4];
	if (j <= *nk - 3) {
	    c1 = abd[(j + 3) * abd_dim1 + 1] * c0;
	    c2 = abd[(j + 2) * abd_dim1 + 2] * c0;
	    c3 = abd[(j + 1) * abd_dim1 + 3] * c0;
	} else {
	    if (j == *nk - 2) {
		c1 = 0.;
		c2 = abd[(j + 2) * abd_dim1 + 2] * c0;
		c3 = abd[(j + 1) * abd_dim1 + 3] * c0;
	    } else {
		if (j == *nk - 1) {
		    c1 = 0.;
		    c2 = 0.;
		    c3 = abd[(j + 1) * abd_dim1 + 3] * c0;
		} else {
		    if (j == *nk) {
			c1 = 0.;
			c2 = 0.;
			c3 = 0.;
		    }
		}
	    }
	}
	p1ip[j * p1ip_dim1 + 1] = 0. - (c1 * wjm3[0] + c2 * wjm3[1] + c3 * 
		wjm3[2]);
	p1ip[j * p1ip_dim1 + 2] = 0. - (c1 * wjm3[1] + c2 * wjm2[0] + c3 * 
		wjm2[1]);
	p1ip[j * p1ip_dim1 + 3] = 0. - (c1 * wjm3[2] + c2 * wjm2[1] + c3 * 
		wjm1[0]);
/* Computing 2nd power */
	d__1 = c0;
/* Computing 2nd power */
	d__2 = c1;
/* Computing 2nd power */
	d__3 = c2;
/* Computing 2nd power */
	d__4 = c3;
	p1ip[j * p1ip_dim1 + 4] = d__1 * d__1 + d__2 * d__2 * wjm3[0] + c1 * (
		float)2. * c2 * wjm3[1] + c1 * (float)2. * c3 * wjm3[2] + 
		d__3 * d__3 * wjm2[0] + c2 * (float)2. * c3 * wjm2[1] + d__4 *
		 d__4 * wjm1[0];
	wjm3[0] = wjm2[0];
	wjm3[1] = wjm2[1];
	wjm3[2] = p1ip[j * p1ip_dim1 + 2];
	wjm2[0] = wjm1[0];
	wjm2[1] = p1ip[j * p1ip_dim1 + 3];
	wjm1[0] = p1ip[j * p1ip_dim1 + 4];
/* L23124: */
    }
/* L23125: */
    return 0;
} /* sinrp2_ */

/* Subroutine */ int suff2_(n, p, ny, match, y, w, ybar, wbar, ssy, work)
integer *n, *p, *ny, *match;
doublereal *y, *w, *ybar, *wbar, *ssy, *work;
{
    /* System generated locals */
    integer y_dim1, y_offset, ybar_dim1, ybar_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int pack_();
    static doublereal tsum;
    static integer i__, j;
    extern /* Subroutine */ int unpack_();

    /* Parameter adjustments */
    --work;
    --w;
    --match;
    --wbar;
    --ssy;
    ybar_dim1 = *p + 1;
    ybar_offset = ybar_dim1 + 1;
    ybar -= ybar_offset;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;

    /* Function Body */
    pack_(n, p, &match[1], &w[1], &wbar[1]);
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__] = y[i__ + j * y_dim1] * w[i__];
/* L23136: */
	}
/* L23137: */
	pack_(n, p, &match[1], &work[1], &ybar[j * ybar_dim1 + 1]);
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (wbar[i__] > 0.) {
		ybar[i__ + j * ybar_dim1] /= wbar[i__];
	    } else {
		ybar[i__ + j * ybar_dim1] = 0.;
	    }
/* L23138: */
	}
/* L23139: */
	unpack_(n, p, &match[1], &ybar[j * ybar_dim1 + 1], &work[1]);
	tsum = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = y[i__ + j * y_dim1] - work[i__];
	    tsum += w[i__] * (d__1 * d__1);
/* L23142: */
	}
/* L23143: */
	ssy[j] = tsum;
/* L23134: */
    }
/* L23135: */
    return 0;
} /* suff2_ */

/* Subroutine */ int sspl0_(x, y, w, n, p, knot, nk, method, tol, wp, match, 
	nef, center, dfoff, dfmax, cost, lambda, df, cv, gcv, coef, s, lev, 
	xrange, work, ier)
doublereal *x, *y, *w;
integer *n, *p;
doublereal *knot;
integer *nk, *method;
doublereal *tol, *wp;
integer *match, *nef;
logical *center;
doublereal *dfoff, *dfmax, *cost, *lambda, *df, *cv, *gcv, *coef, *s, *lev, *
	xrange, *work;
integer *ier;
{
    /* System generated locals */
    integer y_dim1, y_offset, s_dim1, s_offset, i__1;

    /* Local variables */
    static doublereal temp;
    static integer nefp1;
    extern /* Subroutine */ int sspl1_();
    static integer i__;
    static doublereal xdiff;
    extern /* Subroutine */ int namat_();
    static integer n2;
    static doublereal xmiss, sigtol;
    extern /* Subroutine */ int sknotl_();

    /* Parameter adjustments */
    --w;
    --x;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    --knot;
    --match;
    --lev;
    --coef;
    --xrange;
    --work;

    /* Function Body */
    if (*nef == 0) {
	xmiss = 1e20;
	sigtol = 1e-5;
	namat_(&x[1], &match[1], n, nef, &work[1], &work[*n + 1], &xmiss, &
		sigtol);
	xrange[1] = work[1];
	xrange[2] = work[*nef];
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[match[i__]] = x[i__];
/* L23146: */
	}
/* L23147: */
    }
    xdiff = xrange[2] - xrange[1];
    i__1 = *nef;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__] = (work[i__] - xrange[1]) / xdiff;
/* L23148: */
    }
/* L23149: */
    if (*nk == 0) {
	sknotl_(&work[1], nef, &knot[1], nk);
	*nk += -4;
    }
    if (*dfmax > (doublereal) (*nk)) {
	*dfmax = (doublereal) (*nk);
    }
    if (*cost > 0.) {
	temp = (*n - (doublereal) (*center)) / *cost - *dfoff;
	if (*dfmax > temp) {
	    *dfmax = temp;
	}
    }
    nefp1 = *nef + 1;
    n2 = nefp1 * ((*p << 1) + 2) + 1;
    sspl1_(&x[1], &y[y_offset], &w[1], n, p, &knot[1], nk, method, tol, &wp[1]
	    , &match[1], nef, &nefp1, center, dfoff, dfmax, cost, lambda, df, 
	    cv, gcv, &coef[1], &s[s_offset], &lev[1], &work[1], &work[nefp1 + 
	    1], &work[nefp1 * (*p + 1) + 1], &work[nefp1 * (*p + 2) + 1], &
	    work[n2], &work[n2 + *p * *nk], &work[n2 + (*p + 4) * *nk], &work[
	    n2 + (*p + 8) * *nk], &work[n2 + (*p + 12) * *nk], &work[n2 + (*p 
	    + 16) * *nk], &work[n2 + (*p + 16) * *nk + *p], ier);
    return 0;
} /* sspl0_ */

/* Subroutine */ int sspl1_(x, y, w, n, p, knot, nk, method, tol, wp, match, 
	nef, nefp1, center, dfoff, dfmax, cost, lambda, df, cv, gcv, coef, s, 
	lev, xin, yin, win, sout, xwy, hs, sg, abd, p1ip, ssy, work, ier)
doublereal *x, *y, *w;
integer *n, *p;
doublereal *knot;
integer *nk, *method;
doublereal *tol, *wp;
integer *match, *nef, *nefp1;
logical *center;
doublereal *dfoff, *dfmax, *cost, *lambda, *df, *cv, *gcv, *coef, *s, *lev, *
	xin, *yin, *win, *sout, *xwy, *hs, *sg, *abd, *p1ip, *ssy, *work;
integer *ier;
{
    /* System generated locals */
    integer y_dim1, y_offset, coef_dim1, coef_offset, s_dim1, s_offset, 
	    yin_dim1, yin_offset, sout_dim1, sout_offset, xwy_dim1, 
	    xwy_offset, hs_dim1, hs_offset, sg_dim1, sg_offset, i__1, i__2;

    /* Local variables */
    static doublereal sbar;
    extern /* Subroutine */ int sspl_(), suff2_();
    static integer i__, j;
    extern doublereal wmean_();
    static doublereal tdfoff;
    extern /* Subroutine */ int unpack_();

    /* Parameter adjustments */
    --work;
    --match;
    --w;
    --x;
    --ssy;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    p1ip -= 5;
    abd -= 5;
    sg_dim1 = *nk;
    sg_offset = sg_dim1 + 1;
    sg -= sg_offset;
    hs_dim1 = *nk;
    hs_offset = hs_dim1 + 1;
    hs -= hs_offset;
    xwy_dim1 = *nk;
    xwy_offset = xwy_dim1 + 1;
    xwy -= xwy_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;
    --lev;
    sout_dim1 = *nefp1;
    sout_offset = sout_dim1 + 1;
    sout -= sout_offset;
    --win;
    yin_dim1 = *nefp1;
    yin_offset = yin_dim1 + 1;
    yin -= yin_offset;
    --xin;

    /* Function Body */
    suff2_(n, nef, p, &match[1], &y[y_offset], &w[1], &yin[yin_offset], &win[
	    1], &ssy[1], &work[1]);
    if (*center) {
	if (*cost > 0.) {
	    tdfoff = *dfoff - 1 / *cost;
	}
    }
    sspl_(&xin[1], &yin[yin_offset], &win[1], nef, nefp1, p, &knot[1], nk, 
	    method, tol, &wp[1], &ssy[1], &tdfoff, dfmax, cost, lambda, df, 
	    cv, gcv, &coef[coef_offset], &sout[sout_offset], &lev[1], &xwy[
	    xwy_offset], &hs[hs_offset], &sg[sg_offset], &abd[5], &p1ip[5], 
	    ier);
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	unpack_(n, nef, &match[1], &sout[j * sout_dim1 + 1], &s[j * s_dim1 + 
		1]);
	if (*center) {
	    sbar = wmean_(nef, &sout[j * sout_dim1 + 1], &win[1]);
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s[i__ + j * s_dim1] -= sbar;
/* L23166: */
	    }
/* L23167: */
	}
/* L23162: */
    }
/* L23163: */
    if (*center) {
	*df += -1;
    }
    return 0;
} /* sspl1_ */

/* Subroutine */ int namat_(x, match, n, nef, work, iwork, xmiss, tol)
doublereal *x;
integer *match, *n, *nef;
doublereal *work;
integer *iwork;
doublereal *xmiss, *tol;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal xend;
    static integer i__, index;
    extern /* Subroutine */ int sortdi_();
    static doublereal xstart;

    /* Parameter adjustments */
    --iwork;
    --work;
    --match;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__] = x[i__];
	iwork[i__] = i__;
/* L23170: */
    }
/* L23171: */
    sortdi_(&work[1], n, &iwork[1], &c__1, n);
    xstart = x[iwork[1]];
    index = *n;
    xend = x[iwork[*n]];
L23172:
    if (xend >= *xmiss && index > 1) {
	--index;
	xend = x[iwork[index]];
	goto L23172;
    }
/* L23173: */
    *tol *= xend - xstart;
    index = 1;
    work[1] = xstart;
    i__ = 1;
L23174:
    if (! (i__ <= *n)) {
	goto L23176;
    }
L23177:
    if (x[iwork[i__]] - xstart < *tol) {
	match[iwork[i__]] = index;
	++i__;
	if (i__ > *n) {
	    goto L10;
	}
	goto L23177;
    }
/* L23178: */
    xstart = x[iwork[i__]];
    ++index;
    match[iwork[i__]] = index;
    work[index] = xstart;
/* L23175: */
    ++i__;
    goto L23174;
L23176:
L10:
    if (xstart >= *xmiss) {
	*nef = index - 1;
    } else {
	*nef = index;
    }
    return 0;
} /* namat_ */

/* Subroutine */ int simfit_(x, y, w, n, p, dfoff, cost, wp, gcv, coef, s, 
	type__, center, work)
doublereal *x, *y, *w;
integer *n, *p;
doublereal *dfoff, *cost, *wp, *gcv, *coef, *s;
integer *type__;
logical *center;
doublereal *work;
{
    /* System generated locals */
    integer y_dim1, y_offset, s_dim1, s_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal gcvc, gcvl, xbar, sumw;
    static integer i__, j;
    static doublereal dcent, sx, sy, sxx, sxy, syy;

    /* Parameter adjustments */
    --w;
    --x;
    --work;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    coef -= 3;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;

    /* Function Body */
    dcent = 1 - (doublereal) (*center);
    sumw = 0.;
    gcvc = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sumw += w[i__];
/* L23183: */
    }
/* L23184: */
    if (*type__ != 1) {
	sx = (float)0.;
	sxx = (float)0.;
	gcvl = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sx += w[i__] * x[i__];
/* L23187: */
	}
/* L23188: */
	xbar = sx / sumw;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sxx += w[i__] * (x[i__] - xbar) * x[i__];
/* L23189: */
	}
/* L23190: */
    }
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	sy = 0.;
	syy = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sy += w[i__] * y[i__ + j * y_dim1];
/* L23193: */
	}
/* L23194: */
	work[j] = sy / sumw;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    syy += w[i__] * (y[i__ + j * y_dim1] - work[j]) * y[i__ + j * 
		    y_dim1];
/* L23195: */
	}
/* L23196: */
	gcvc += wp[j] * syy;
	if (*type__ != 1) {
	    sxy = (float)0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sxy += w[i__] * (x[i__] - xbar) * y[i__ + j * y_dim1];
/* L23199: */
	    }
/* L23200: */
	    coef[(j << 1) + 2] = sxy / sxx;
	    gcvl += wp[j] * (syy - sxy * coef[(j << 1) + 2]);
	}
/* L23191: */
    }
/* L23192: */
    if (*type__ == 0) {
/* Computing 2nd power */
	d__1 = 1 - (*dfoff * *cost + dcent) / sumw;
	gcvc /= sumw * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = 1 - (dcent + (*dfoff + 1) * *cost) / sumw;
	gcvl /= sumw * (d__1 * d__1);
	if (gcvc <= gcvl) {
	    *type__ = 1;
	    *gcv = gcvc;
	} else {
	    *type__ = 2;
	    *gcv = gcvl;
	}
    } else {
	if (*type__ == 1) {
/* Computing 2nd power */
	    d__1 = 1 - (*dfoff * *cost + dcent) / sumw;
	    *gcv = gcvc / (sumw * (d__1 * d__1));
	} else {
/* Computing 2nd power */
	    d__1 = 1 - (dcent + (*dfoff + 1) * *cost) / sumw;
	    *gcv = gcvl / (sumw * (d__1 * d__1));
	}
    }
    if (*type__ == 1) {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    coef[(j << 1) + 1] = work[j] * dcent;
	    coef[(j << 1) + 2] = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s[i__ + j * s_dim1] = coef[(j << 1) + 1];
/* L23211: */
	    }
/* L23212: */
/* L23209: */
	}
/* L23210: */
    } else {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    coef[(j << 1) + 1] = work[j] * dcent - xbar * coef[(j << 1) + 2];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s[i__ + j * s_dim1] = coef[(j << 1) + 1] + coef[(j << 1) + 2] 
			* x[i__];
/* L23215: */
	    }
/* L23216: */
/* L23213: */
	}
/* L23214: */
    }
    return 0;
} /* simfit_ */

/* Subroutine */ int sspl2_(x, y, w, n, p, knot, nk, wp, match, nef, dfoff, 
	dfmax, cost, lambda, df, gcv, coef, s, type__, center, xrange, work, 
	tol, ier)
doublereal *x, *y, *w;
integer *n, *p;
doublereal *knot;
integer *nk;
doublereal *wp;
integer *match, *nef;
doublereal *dfoff, *dfmax, *cost, *lambda, *df, *gcv, *coef, *s;
integer *type__;
logical *center;
doublereal *xrange, *work, *tol;
integer *ier;
{
    /* System generated locals */
    integer y_dim1, y_offset, s_dim1, s_offset, i__1, i__2;

    /* Local variables */
    static doublereal coef1, coef2;
    extern /* Subroutine */ int sspl0_();
    static integer i__, j;
    static doublereal cv;
    static integer method;
    extern /* Subroutine */ int simfit_();
    static doublereal gcv1;

    /* Parameter adjustments */
    --match;
    --w;
    --x;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    --wp;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    --knot;
    --coef;
    --xrange;
    --work;

    /* Function Body */
    if (*type__ == 3) {
	method = 1;
	sspl0_(&x[1], &y[y_offset], &w[1], n, p, &knot[1], nk, &method, tol, &
		wp[1], &match[1], nef, center, dfoff, dfmax, cost, lambda, df,
		 &cv, gcv, &coef[1], &s[s_offset], &work[1], &xrange[1], &
		work[*n + 1], ier);
	return 0;
    }
    if (*type__ > 0) {
	simfit_(&x[1], &y[y_offset], &w[1], n, p, dfoff, cost, &wp[1], gcv, &
		coef[1], &s[s_offset], type__, center, &work[1]);
	*df = (doublereal) (*type__) - (doublereal) (*center);
	return 0;
    }
    method = 3;
    sspl0_(&x[1], &y[y_offset], &w[1], n, p, &knot[1], nk, &method, tol, &wp[
	    1], &match[1], nef, center, dfoff, dfmax, cost, lambda, df, &cv, 
	    gcv, &coef[1], &s[s_offset], &work[1], &xrange[1], &work[*n + 1], 
	    ier);
    gcv1 = *gcv;
    simfit_(&x[1], &y[y_offset], &w[1], n, p, dfoff, cost, &wp[1], gcv, &work[
	    1], &work[(*p << 1) + 1], type__, center, &work[(*n + 2) * *p + 1]
	    );
    if (*gcv <= gcv1) {
	*df = (doublereal) (*type__) - (doublereal) (*center);
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    coef1 = work[(j - 1 << 1) + 1];
	    coef2 = work[(j - 1 << 1) + 2];
	    if (*type__ == 1) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[i__ + j * s_dim1] = coef1;
/* L23227: */
		}
/* L23228: */
	    } else {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[i__ + j * s_dim1] = coef1 + coef2 * x[i__];
/* L23229: */
		}
/* L23230: */
	    }
	    coef[(j - 1 << 1) + 1] = coef1;
	    coef[(j - 1 << 1) + 2] = coef2;
/* L23223: */
	}
/* L23224: */
    } else {
	*type__ = 3;
	*gcv = gcv1;
    }
    return 0;
} /* sspl2_ */

/* Subroutine */ int psspl2_(x, n, p, knot, nk, xrange, coef, coefl, s, order,
	 type__)
doublereal *x;
integer *n, *p;
doublereal *knot;
integer *nk;
doublereal *xrange, *coef, *coefl, *s;
integer *order, *type__;
{
    /* System generated locals */
    integer coef_dim1, coef_offset, s_dim1, s_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal ytemp;
    extern /* Subroutine */ int psspl_();
    static integer i23231;

    /* Parameter adjustments */
    --x;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;
    --xrange;
    coefl -= 3;

    /* Function Body */
    i23231 = *type__;
    goto L23231;
L23233:
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	if (*order >= 1) {
	    ytemp = 0.;
	} else {
	    ytemp = coefl[(j << 1) + 1];
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s[i__ + j * s_dim1] = ytemp;
/* L23238: */
	}
/* L23239: */
/* L23234: */
    }
/* L23235: */
    goto L23232;
L23240:
    if (*order >= 1) {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    if (*order == 1) {
		ytemp = coefl[(j << 1) + 2];
	    } else {
		ytemp = 0.;
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s[i__ + j * s_dim1] = ytemp;
/* L23247: */
	    }
/* L23248: */
/* L23243: */
	}
/* L23244: */
    } else {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s[i__ + j * s_dim1] = coefl[(j << 1) + 1] + coefl[(j << 1) + 
			2] * x[i__];
/* L23251: */
	    }
/* L23252: */
/* L23249: */
	}
/* L23250: */
    }
    goto L23232;
L23253:
    psspl_(&x[1], n, p, &knot[1], nk, &xrange[1], &coef[coef_offset], &s[
	    s_offset], order);
    goto L23232;
L23231:
    if (i23231 == 1) {
	goto L23233;
    }
    if (i23231 == 2) {
	goto L23240;
    }
    if (i23231 == 3) {
	goto L23253;
    }
L23232:
    return 0;
} /* psspl2_ */

/* Subroutine */ int psspl_(x, n, p, knot, nk, xrange, coef, s, order)
doublereal *x;
integer *n, *p;
doublereal *knot;
integer *nk;
doublereal *xrange, *coef, *s;
integer *order;
{
    /* System generated locals */
    integer coef_dim1, coef_offset, s_dim1, s_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd();

    /* Local variables */
    static doublereal ends[2], xdif, endv[2], xmin;
    static integer i__, j, where;
    static doublereal xends[2], stemp;
    extern doublereal bvalue_();
    static doublereal xcs;
    static integer i23270;

    /* Parameter adjustments */
    --x;
    s_dim1 = *n;
    s_offset = s_dim1 + 1;
    s -= s_offset;
    coef_dim1 = *nk;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --knot;
    --xrange;

    /* Function Body */
    if (*order > 2 || *order < 0) {
	return 0;
    }
    xdif = xrange[2] - xrange[1];
    xmin = xrange[1];
    xends[0] = 0.;
    xends[1] = 1.;
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	if (*order == 0) {
	    endv[0] = bvalue_(&knot[1], &coef[j * coef_dim1 + 1], nk, &c__4, &
		    c_b167, &c__0);
	    endv[1] = bvalue_(&knot[1], &coef[j * coef_dim1 + 1], nk, &c__4, &
		    c_b170, &c__0);
	}
	if (*order <= 1) {
	    ends[0] = bvalue_(&knot[1], &coef[j * coef_dim1 + 1], nk, &c__4, &
		    c_b167, &c__1);
	    ends[1] = bvalue_(&knot[1], &coef[j * coef_dim1 + 1], nk, &c__4, &
		    c_b170, &c__1);
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xcs = (x[i__] - xmin) / xdif;
	    where = 0;
	    if (xcs < 0.) {
		where = 1;
	    }
	    if (xcs > 1.) {
		where = 2;
	    }
	    if (where > 0) {
		i23270 = *order;
		goto L23270;
L23272:
		stemp = endv[where - 1] + (xcs - xends[where - 1]) * ends[
			where - 1];
		goto L23271;
L23273:
		stemp = ends[where - 1];
		goto L23271;
L23274:
		stemp = 0.;
		goto L23271;
L23270:
		if (i23270 == 0) {
		    goto L23272;
		}
		if (i23270 == 1) {
		    goto L23273;
		}
		if (i23270 == 2) {
		    goto L23274;
		}
L23271:
		;
	    } else {
		stemp = bvalue_(&knot[1], &coef[j * coef_dim1 + 1], nk, &c__4,
			 &xcs, order);
	    }
	    if (*order > 0) {
		d__1 = (doublereal) (*order);
		s[i__ + j * s_dim1] = stemp / pow_dd(&xdif, &d__1);
	    } else {
		s[i__ + j * s_dim1] = stemp;
	    }
/* L23262: */
	}
/* L23263: */
/* L23256: */
    }
/* L23257: */
    return 0;
} /* psspl_ */

