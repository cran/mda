/* dmarss.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

union {
    struct {
	logical tracec;
    } _1;
    struct {
	logical trace;
    } _2;
} _BLNK__;

#define _BLNK__1 (_BLNK__._1)
#define _BLNK__2 (_BLNK__._2)

/* Table of constant values */

static integer c__11 = 11;
static integer c__1 = 1;
static integer c__12 = 12;
static integer c__21 = 21;
static integer c__3 = 3;
static integer c__19 = 19;
static integer c__4 = 4;
static integer c__43 = 43;
static integer c__10 = 10;
static integer c__7 = 7;
static integer c__13 = 13;
static integer c__25 = 25;
static integer c__5 = 5;

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int marss_(nx, n, p, nclass, y, x, w, tagx, maxorder, mmax, 
	penalty, thresh, forwstep, interms, prune, bx, fullin, lenb, bestgcv, 
	bestin, flag__, cut, dir, res, alpha, beta, scrat, iscrat, trace)
integer *nx, *n, *p, *nclass;
doublereal *y, *x, *w;
integer *tagx, *maxorder, *mmax;
doublereal *penalty, *thresh;
logical *forwstep;
integer *interms;
logical *prune;
doublereal *bx;
integer *fullin, *lenb;
doublereal *bestgcv;
integer *bestin, *flag__;
doublereal *cut, *dir, *res, *alpha, *beta, *scrat;
integer *iscrat;
logical *trace;
{
    /* System generated locals */
    integer tagx_dim1, tagx_offset, flag_dim1, flag_offset, y_dim1, y_offset, 
	    x_dim1, x_offset, bx_dim1, bx_offset, cut_dim1, cut_offset, 
	    dir_dim1, dir_offset, res_dim1, res_offset, beta_dim1, 
	    beta_offset;

    /* Local variables */
    extern /* Subroutine */ int marsnew1_();
    static integer len10, len11, len12, len13, len14, n1, n2, n3, n4, n5, n6, 
	    n7, n8, n9, n10, n11, n12, n13, n14, n15, len1, len2, len3, len4, 
	    len5, len6, len7, len8, len9;

    /* Parameter adjustments */
    --w;
    tagx_dim1 = *nx;
    tagx_offset = tagx_dim1 + 1;
    tagx -= tagx_offset;
    x_dim1 = *nx;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --alpha;
    res_dim1 = *nx;
    res_offset = res_dim1 + 1;
    res -= res_offset;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    beta_dim1 = *mmax;
    beta_offset = beta_dim1 + 1;
    beta -= beta_offset;
    dir_dim1 = *mmax;
    dir_offset = dir_dim1 + 1;
    dir -= dir_offset;
    cut_dim1 = *mmax;
    cut_offset = cut_dim1 + 1;
    cut -= cut_offset;
    flag_dim1 = *mmax;
    flag_offset = flag_dim1 + 1;
    flag__ -= flag_offset;
    --bestin;
    --fullin;
    bx_dim1 = *nx;
    bx_offset = bx_dim1 + 1;
    bx -= bx_offset;
    --scrat;
    --iscrat;

    /* Function Body */
    _BLNK__1.tracec = *trace;
    len1 = *n * *mmax;
    len2 = *mmax;
    len3 = *mmax * *mmax;
    len4 = *mmax * *nclass;
    len5 = *nclass;
    len6 = *mmax;
    len7 = *mmax;
    len8 = *nclass;
    len9 = *n;
    len10 = *n * *mmax;
    len11 = *mmax * *mmax;
    len12 = *mmax * *nclass;
    len13 = *mmax * *mmax;
    len14 = *mmax * *mmax;
    n1 = 1;
    n2 = n1 + len1;
    n3 = n2 + len2;
    n4 = n3 + len3;
    n5 = n4 + len4;
    n6 = n5 + len5;
    n7 = n6 + len6;
    n8 = n7 + len7;
    n9 = n8 + len8;
    n10 = n9 + len9;
    n11 = n10 + len10;
    n12 = n11 + len11;
    n13 = n12 + len12;
    n14 = n13 + len13;
    n15 = n14 + len14;
    marsnew1_(nx, n, p, nclass, &y[y_offset], &x[x_offset], &w[1], &tagx[
	    tagx_offset], maxorder, mmax, &bx[bx_offset], bestgcv, &bestin[1],
	     &fullin[1], lenb, &flag__[flag_offset], &cut[cut_offset], &dir[
	    dir_offset], &res[res_offset], &alpha[1], &beta[beta_offset], 
	    penalty, thresh, forwstep, interms, prune, &scrat[1], &scrat[n2], 
	    &scrat[n3], &scrat[n4], &scrat[n5], &scrat[n6], &scrat[n7], &
	    scrat[n8], &scrat[n9], &scrat[n10], &scrat[n11], &scrat[n12], &
	    scrat[n13], &scrat[n14], &scrat[n15], &iscrat[1], &iscrat[*mmax + 
	    1], &iscrat[(*mmax << 1) + 1], &iscrat[*mmax * 3 + 1]);
    return 0;
} /* marss_ */

/* Subroutine */ int marsnew1_(nx, n, p, nclass, y, x, w, tagx, maxorder, 
	mmax, bx, bestgcv, bestin, fullin, lenb, flag__, cut, dir, res, alpha,
	 beta, penalty, thresh, forwstep, interms, prune, bxorth, bxorthm, 
	cov, covsy, ybar, scr1, scr5, scr6, temp, bxsc, r__, betasc, varsc, 
	var, work, termlen, in, tempin, qpivot)
integer *nx, *n, *p, *nclass;
doublereal *y, *x, *w;
integer *tagx, *maxorder, *mmax;
doublereal *bx, *bestgcv;
integer *bestin, *fullin, *lenb, *flag__;
doublereal *cut, *dir, *res, *alpha, *beta, *penalty, *thresh;
logical *forwstep;
integer *interms;
logical *prune;
doublereal *bxorth, *bxorthm, *cov, *covsy, *ybar, *scr1, *scr5, *scr6, *temp,
	 *bxsc, *r__, *betasc, *varsc, *var, *work;
integer *termlen, *in, *tempin, *qpivot;
{
    /* System generated locals */
    integer flag_dim1, flag_offset, tagx_dim1, tagx_offset, cov_dim1, 
	    cov_offset, covsy_dim1, covsy_offset, x_dim1, x_offset, bx_dim1, 
	    bx_offset, bxorth_dim1, bxorth_offset, y_dim1, y_offset, cut_dim1,
	     cut_offset, dir_dim1, dir_offset, beta_dim1, beta_offset, 
	    bxsc_dim1, bxsc_offset, r_dim1, r_offset, res_dim1, res_offset, 
	    betasc_dim1, betasc_offset, varsc_dim1, varsc_offset, var_dim1, 
	    var_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static logical cvar;
    static integer imax, jmax, kmax;
    static doublereal prevcrit, temp1, temp2, temp7;
    static integer i__, j, k;
    static doublereal dofit;
    static integer qrank, itemp[4];
    extern /* Subroutine */ int qrreg_();
    static doublereal rtemp[4], tolbx;
    static integer i1, i2;
    extern /* Subroutine */ int intpr_();
    static integer kc, ii, ik;
    static logical go;
    static integer kk, jo, nt;
    extern /* Subroutine */ int dblepr_(), addtrm_();
    static integer nterms, nterms2;
    static doublereal cmm, gcv, rss, doftemp, stopfac, critmax, gcvnull;
    extern /* Subroutine */ int orthreg_();
    static logical newform;
    static integer minterm;
    static doublereal rssfull, temprss, rsstemp, rssnull;

    /* Parameter adjustments */
    --temp;
    --w;
    tagx_dim1 = *nx;
    tagx_offset = tagx_dim1 + 1;
    tagx -= tagx_offset;
    x_dim1 = *nx;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --scr6;
    --ybar;
    --alpha;
    res_dim1 = *nx;
    res_offset = res_dim1 + 1;
    res -= res_offset;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    --qpivot;
    --tempin;
    --termlen;
    var_dim1 = *mmax;
    var_offset = var_dim1 + 1;
    var -= var_offset;
    varsc_dim1 = *mmax;
    varsc_offset = varsc_dim1 + 1;
    varsc -= varsc_offset;
    betasc_dim1 = *mmax;
    betasc_offset = betasc_dim1 + 1;
    betasc -= betasc_offset;
    r_dim1 = *mmax;
    r_offset = r_dim1 + 1;
    r__ -= r_offset;
    bxsc_dim1 = *n;
    bxsc_offset = bxsc_dim1 + 1;
    bxsc -= bxsc_offset;
    --scr5;
    --scr1;
    covsy_dim1 = *mmax;
    covsy_offset = covsy_dim1 + 1;
    covsy -= covsy_offset;
    cov_dim1 = *mmax;
    cov_offset = cov_dim1 + 1;
    cov -= cov_offset;
    --bxorthm;
    bxorth_dim1 = *n;
    bxorth_offset = bxorth_dim1 + 1;
    bxorth -= bxorth_offset;
    beta_dim1 = *mmax;
    beta_offset = beta_dim1 + 1;
    beta -= beta_offset;
    dir_dim1 = *mmax;
    dir_offset = dir_dim1 + 1;
    dir -= dir_offset;
    cut_dim1 = *mmax;
    cut_offset = cut_dim1 + 1;
    cut -= cut_offset;
    flag_dim1 = *mmax;
    flag_offset = flag_dim1 + 1;
    flag__ -= flag_offset;
    --fullin;
    --bestin;
    bx_dim1 = *nx;
    bx_offset = bx_dim1 + 1;
    bx -= bx_offset;
    --work;

    /* Function Body */
    tolbx = (float).01;
    stopfac = (float)10.;
    prevcrit = (float)1e10;
    if (*interms == 1) {
	dofit = 0.;
    } else {
	dofit = 0.;
	i__1 = *lenb;
	for (j = 2; j <= i__1; ++j) {
	    dofit += fullin[j];
/* L23002: */
	}
/* L23003: */
	nterms = *interms;
    }
    if (*forwstep) {
	fullin[1] = 1;
	i__1 = *mmax;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    fullin[i__] = 0;
/* L23006: */
	}
/* L23007: */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    w[i__] = 1.;
/* L23008: */
	}
/* L23009: */
	i__1 = *mmax;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    termlen[i__] = 0;
	    i__2 = *p;
	    for (j = 1; j <= i__2; ++j) {
		flag__[i__ + j * flag_dim1] = 0;
		cut[i__ + j * cut_dim1] = 0.;
/* L23012: */
	    }
/* L23013: */
/* L23010: */
	}
/* L23011: */
	nterms = 1;
	nterms2 = 2;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    bx[i__ + bx_dim1] = 1.;
	    bxorth[i__ + bxorth_dim1] = (float)1. / sqrt((doublereal) (*n));
/* L23014: */
	}
/* L23015: */
	bxorthm[1] = 1 / sqrt((doublereal) (*n));
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *mmax;
	    for (j = 1; j <= i__2; ++j) {
		bx[i__ + j * bx_dim1] = (float)0.;
/* L23018: */
	    }
/* L23019: */
/* L23016: */
	}
/* L23017: */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    bx[i__ + bx_dim1] = 1.;
/* L23020: */
	}
/* L23021: */
	i__1 = *nclass;
	for (k = 1; k <= i__1; ++k) {
	    ybar[k] = (float)0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ybar[k] += y[i__ + k * y_dim1] / *n;
/* L23024: */
	    }
/* L23025: */
/* L23022: */
	}
/* L23023: */
	if (*interms == 1) {
	    rssnull = (float)0.;
	    i__1 = *nclass;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
		    d__1 = y[i__ + k * y_dim1] - ybar[k];
		    rssnull += d__1 * d__1;
/* L23030: */
		}
/* L23031: */
/* L23028: */
	    }
/* L23029: */
	} else {
	    rssnull = (float)0.;
	    i__1 = *nclass;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
		    d__1 = res[i__ + k * res_dim1];
		    rssnull += d__1 * d__1;
/* L23034: */
		}
/* L23035: */
/* L23032: */
	    }
/* L23033: */
	}
	rss = rssnull;
	cmm = dofit + 1 + *penalty * (dofit * (float).5);
/* Computing 2nd power */
	d__1 = (float)1. - cmm / *n;
	gcvnull = rssnull / *n / (d__1 * d__1);
	if (_BLNK__2.trace) {
	    dblepr_("initial rss=", &c__11, &rssnull, &c__1, 12L);
	}
	if (_BLNK__2.trace) {
	    dblepr_("initial gcv=", &c__11, &gcvnull, &c__1, 12L);
	}
	*lenb = 1;
	ii = *interms - 1;
	go = TRUE_;
L23040:
	if (ii < *mmax - 1 && rss / rssnull > *thresh && go) {
	    ii += 2;
	    i__1 = nterms;
	    for (i1 = 1; i1 <= i__1; ++i1) {
		i__2 = nterms;
		for (i2 = 1; i2 <= i__2; ++i2) {
		    cov[i1 + i2 * cov_dim1] = 0.;
/* L23044: */
		}
/* L23045: */
/* L23042: */
	    }
/* L23043: */
	    i__1 = nterms;
	    for (j = 1; j <= i__1; ++j) {
		cov[j + j * cov_dim1] = (float)0.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    cov[j + j * cov_dim1] += (bxorth[i__ + j * bxorth_dim1] - 
			    bxorthm[j]) * (bxorth[i__ + j * bxorth_dim1] - 
			    bxorthm[j]);
/* L23048: */
		}
/* L23049: */
/* L23046: */
	    }
/* L23047: */
	    i__1 = *nclass;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = nterms;
		for (j = 1; j <= i__2; ++j) {
		    covsy[j + k * covsy_dim1] = (float)0.;
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			covsy[j + k * covsy_dim1] += (y[i__ + k * y_dim1] - 
				ybar[k]) * bxorth[i__ + j * bxorth_dim1];
/* L23054: */
		    }
/* L23055: */
/* L23052: */
		}
/* L23053: */
/* L23050: */
	    }
/* L23051: */
	    i__1 = *mmax;
	    for (ik = 1; ik <= i__1; ++ik) {
		tempin[ik] = fullin[ik];
/* L23056: */
	    }
/* L23057: */
	    addtrm_(nx, &bx[bx_offset], &tempin[1], &bxorth[bxorth_offset], &
		    bxorthm[1], p, n, nclass, &rss, &prevcrit, &cov[
		    cov_offset], &covsy[covsy_offset], &y[y_offset], &ybar[1],
		     &x[x_offset], &tagx[tagx_offset], &w[1], &termlen[1], 
		    mmax, &tolbx, &nterms, &flag__[flag_offset], maxorder, &
		    scr1[1], &scr5[1], &scr6[1], &imax, &jmax, &kmax, &
		    critmax, &newform, &bxsc[bxsc_offset], &r__[r_offset], &
		    betasc[betasc_offset], &temp[1]);
	    doftemp = dofit;
	    doftemp += 1;
	    if (imax > 1 && newform) {
		doftemp += 1;
	    }
	    temprss = rss - critmax;
	    cmm = doftemp + 1 + *penalty * (doftemp * (float).5);
/* Computing 2nd power */
	    d__1 = (float)1. - cmm / *n;
	    gcv = temprss / *n / (d__1 * d__1);
	    go = FALSE_;
	    if (critmax / rss > *thresh && gcv / gcvnull < stopfac) {
		go = TRUE_;
		dofit = doftemp;
		rss -= critmax;
		kk = tagx[imax + jmax * tagx_dim1];
/* L256: */
		itemp[0] = jmax;
		itemp[1] = imax;
		itemp[2] = kmax;
		rtemp[0] = critmax;
		rtemp[1] = x[kk + jmax * x_dim1];
		rtemp[2] = rss;
		rtemp[3] = gcv;
		if (_BLNK__2.trace) {
		    intpr_("adding term ", &c__12, &ii, &c__1, 12L);
		}
		if (_BLNK__2.trace) {
		    intpr_("var, sp index, parent", &c__21, itemp, &c__3, 21L)
			    ;
		}
		if (_BLNK__2.trace) {
		    dblepr_("critmax cut rss gcv", &c__19, rtemp, &c__4, 19L);
		}
		prevcrit = critmax;
		i__1 = *p;
		for (j = 1; j <= i__1; ++j) {
		    flag__[ii + j * flag_dim1] = flag__[kmax + j * flag_dim1];
		    flag__[ii + 1 + j * flag_dim1] = flag__[kmax + j * 
			    flag_dim1];
		    cut[ii + j * cut_dim1] = cut[kmax + j * cut_dim1];
		    cut[ii + 1 + j * cut_dim1] = cut[kmax + j * cut_dim1];
		    dir[ii + j * dir_dim1] = dir[kmax + j * dir_dim1];
		    dir[ii + 1 + j * dir_dim1] = dir[kmax + j * dir_dim1];
/* L23068: */
		}
/* L23069: */
		termlen[ii] = termlen[kmax] + 1;
		termlen[ii + 1] = termlen[kmax] + 1;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    temp[i__] = x[tagx[i__ + jmax * tagx_dim1] + jmax * 
			    x_dim1];
/* L23070: */
		}
/* L23071: */
		temp1 = temp[imax];
		fullin[ii] = 1;
		if (imax > 1 && newform) {
		    fullin[ii + 1] = 1;
		}
		flag__[ii + jmax * flag_dim1] = 1;
		flag__[ii + 1 + jmax * flag_dim1] = 1;
		cut[ii + jmax * cut_dim1] = temp1;
		cut[ii + 1 + jmax * cut_dim1] = temp1;
		dir[ii + jmax * dir_dim1] = 1.;
		dir[ii + 1 + jmax * dir_dim1] = -1.;
		if (fullin[ii + 1] == 0) {
		    termlen[ii + 1] = *maxorder + 1;
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (x[i__ + jmax * x_dim1] - temp1 > 0.) {
			bx[i__ + ii * bx_dim1] = bx[i__ + kmax * bx_dim1] * (
				x[i__ + jmax * x_dim1] - temp1);
		    }
		    if (temp1 - x[i__ + jmax * x_dim1] >= 0.) {
			bx[i__ + (ii + 1) * bx_dim1] = bx[i__ + kmax * 
				bx_dim1] * (temp1 - x[i__ + jmax * x_dim1]);
		    }
/* L23076: */
		}
/* L23077: */
		if (nterms == 1) {
		    temp1 = (float)0.;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			temp1 += bx[i__ + (bx_dim1 << 1)] / *n;
/* L23084: */
		    }
/* L23085: */
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			bxorth[i__ + (bxorth_dim1 << 1)] = bx[i__ + (bx_dim1 
				<< 1)] - temp1;
/* L23086: */
		    }
/* L23087: */
		} else {
		    orthreg_(n, n, &nterms, &bxorth[bxorth_offset], &fullin[1]
			    , &bx[ii * bx_dim1 + 1], &bxorth[nterms2 * 
			    bxorth_dim1 + 1]);
		}
		if (fullin[ii + 1] == 1) {
		    i__1 = nterms + 1;
		    orthreg_(n, n, &i__1, &bxorth[bxorth_offset], &fullin[1], 
			    &bx[(ii + 1) * bx_dim1 + 1], &bxorth[(nterms2 + 1)
			     * bxorth_dim1 + 1]);
		} else {
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			bxorth[i__ + (nterms2 + 1) * bxorth_dim1] = 0.;
/* L23090: */
		    }
/* L23091: */
		}
		bxorthm[nterms2] = (float)0.;
		bxorthm[nterms2 + 1] = (float)0.;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    bxorthm[nterms2] += bxorth[i__ + nterms2 * bxorth_dim1] / 
			    *n;
		    bxorthm[nterms2 + 1] += bxorth[i__ + (nterms2 + 1) * 
			    bxorth_dim1] / *n;
/* L23092: */
		}
/* L23093: */
		temp1 = (float)0.;
		temp2 = (float)0.;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		    d__1 = bxorth[i__ + nterms2 * bxorth_dim1];
		    temp1 += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = bxorth[i__ + (nterms2 + 1) * bxorth_dim1];
		    temp2 += d__1 * d__1;
/* L23094: */
		}
/* L23095: */
		if (temp1 > (float)0.) {
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			bxorth[i__ + nterms2 * bxorth_dim1] /= sqrt(temp1);
/* L23098: */
		    }
/* L23099: */
		}
		if (temp2 > (float)0.) {
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			bxorth[i__ + (nterms2 + 1) * bxorth_dim1] /= sqrt(
				temp2);
/* L23102: */
		    }
/* L23103: */
		}
		*lenb += 2;
		nterms += 2;
		nterms2 += 2;
	    }
	    goto L23040;
	}
/* L23041: */
	rtemp[0] = rss / rssnull;
	rtemp[1] = critmax / rss;
	rtemp[2] = gcv / gcvnull;
	if (_BLNK__2.trace) {
	    dblepr_("stopping forw step; rss crit and gcv ratios", &c__43, 
		    rtemp, &c__3, 43L);
	}
	if (_BLNK__2.trace) {
	    if (rss / rssnull <= *thresh) {
		d__1 = rss / rssnull;
		dblepr_("rss ratio=", &c__10, &d__1, &c__1, 10L);
	    }
	    if (critmax / rss <= *thresh) {
		d__1 = critmax / rss;
		dblepr_("crit ratio=", &c__11, &d__1, &c__1, 11L);
	    }
	    dblepr_("critmax", &c__7, &critmax, &c__1, 7L);
	    dblepr_("rss", &c__3, &rss, &c__1, 3L);
	    if (gcv / gcvnull > stopfac) {
		d__1 = gcv / gcvnull;
		dblepr_("gcv ratio=", &c__10, &d__1, &c__1, 10L);
	    }
	}
    }
    dofit = -1.;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bestin[i__] = fullin[i__];
	dofit += fullin[i__];
/* L23114: */
    }
/* L23115: */
    if (_BLNK__2.trace) {
	intpr_("aft forw step", &c__13, &nterms, &c__1, 13L);
    }
    qrreg_(nx, n, mmax, lenb, nclass, &bx[bx_offset], &bxsc[bxsc_offset], &
	    bestin[1], &y[y_offset], &qpivot[1], &qrank, &beta[beta_offset], &
	    res[res_offset], &rss, &cvar, &var[var_offset], &varsc[
	    varsc_offset], &scr1[1], &work[1]);
    nt = (integer) (dofit + 1);
    if (qrank < nt) {
	i__1 = nt;
	for (i__ = qrank + 1; i__ <= i__1; ++i__) {
	    bestin[qpivot[i__]] = 0;
	    fullin[qpivot[i__]] = 0;
	    dofit += -1;
/* L23120: */
	}
/* L23121: */
    }
    cvar = TRUE_;
    rssfull = rss;
    cmm = dofit + 1 + *penalty * (dofit * (float).5);
/* Computing 2nd power */
    d__1 = (float)1. - cmm / *n;
    *bestgcv = rss / *n / (d__1 * d__1);
    rtemp[0] = *bestgcv;
    rtemp[1] = rssfull;
    rtemp[2] = dofit;
    if (_BLNK__2.trace) {
	dblepr_("full model: gcv rss dofit", &c__25, rtemp, &c__3, 25L);
    }
    if (_BLNK__2.trace) {
	intpr_("terms", &c__5, &fullin[1], lenb, 5L);
    }
    if (*prune) {
	i__1 = *mmax;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tempin[i__] = bestin[i__];
/* L23128: */
	}
/* L23129: */
L23130:
	if (dofit > 0.) {
	    jo = 1;
	    rsstemp = (float)1e100;
	    minterm = 0;
	    i__1 = *lenb;
	    for (ii = 2; ii <= i__1; ++ii) {
		if (tempin[ii] == 1) {
		    ++jo;
		    temp7 = (float)0.;
		    i__2 = *nclass;
		    for (kc = 1; kc <= i__2; ++kc) {
/* Computing 2nd power */
			d__1 = beta[jo + kc * beta_dim1];
			temp7 += d__1 * d__1 / var[jo + jo * var_dim1];
/* L23136: */
		    }
/* L23137: */
		    if (temp7 < rsstemp) {
			minterm = ii;
			rsstemp = temp7;
		    }
		}
/* L23132: */
	    }
/* L23133: */
	    rss += rsstemp;
	    dofit += -1;
	    cmm = dofit + (float)1. + *penalty * (dofit * (float).5);
/* Computing 2nd power */
	    d__1 = (float)1. - cmm / *n;
	    gcv = rss / *n / (d__1 * d__1);
	    tempin[minterm] = 0;
/* L100: */
	    if (gcv < *bestgcv) {
		*bestgcv = gcv;
		i__1 = *mmax;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    bestin[i__] = tempin[i__];
/* L23142: */
		}
/* L23143: */
	    }
	    if (dofit > 0.) {
		cvar = TRUE_;
		qrreg_(nx, n, mmax, lenb, nclass, &bx[bx_offset], &bxsc[
			bxsc_offset], &tempin[1], &y[y_offset], &qpivot[1], &
			qrank, &beta[beta_offset], &res[res_offset], &rss, &
			cvar, &var[var_offset], &varsc[varsc_offset], &scr1[1]
			, &work[1]);
	    }
	    goto L23130;
	}
/* L23131: */
	qrreg_(nx, n, mmax, lenb, nclass, &bx[bx_offset], &bxsc[bxsc_offset], 
		&bestin[1], &y[y_offset], &qpivot[1], &qrank, &beta[
		beta_offset], &res[res_offset], &rss, &cvar, &var[var_offset],
		 &varsc[varsc_offset], &scr1[1], &work[1]);
/* L101: */
	if (_BLNK__2.trace) {
	    intpr_("best model", &c__10, &bestin[1], lenb, 10L);
	}
	if (_BLNK__2.trace) {
	    dblepr_(" gcv=", &c__4, bestgcv, &c__1, 5L);
	}
    }
    return 0;
} /* marsnew1_ */

/* Subroutine */ int addtrm_(nx, bx, tempin, bxorth, bxorthm, p, n, nclass, 
	rss, prevcrit, cov, covsy, y, ybar, x, tagx, w, termlen, mmax, tolbx, 
	nterms, flag__, maxorder, scr1, scr5, scr6, imax, jmax, kmax, critmax,
	 newform, bxsc, r__, betasc, scrat)
integer *nx;
doublereal *bx;
integer *tempin;
doublereal *bxorth, *bxorthm;
integer *p, *n, *nclass;
doublereal *rss, *prevcrit, *cov, *covsy, *y, *ybar, *x;
integer *tagx;
doublereal *w;
integer *termlen, *mmax;
doublereal *tolbx;
integer *nterms, *flag__, *maxorder;
doublereal *scr1, *scr5, *scr6;
integer *imax, *jmax, *kmax;
doublereal *critmax;
logical *newform;
doublereal *bxsc, *r__, *betasc, *scrat;
{
    /* System generated locals */
    integer flag_dim1, flag_offset, tagx_dim1, tagx_offset, cov_dim1, 
	    cov_offset, covsy_dim1, covsy_offset, x_dim1, x_offset, bx_dim1, 
	    bx_offset, bxorth_dim1, bxorth_offset, y_dim1, y_offset, 
	    bxsc_dim1, bxsc_offset, r_dim1, r_offset, betasc_dim1, 
	    betasc_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double log(), sqrt();

    /* Local variables */
    static doublereal sumb;
    static integer iendspan;
    static doublereal crittemp;
    static logical tnewform;
    static doublereal temp1, temp2;
    static integer i__, j, k, m, v, i1, k0;
    static doublereal sumbx;
    static integer kc;
    static doublereal sumbx2;
    static integer ii, jk, kk, nm, mm;
    static doublereal st, su;
    static integer kk1, nterms2, jjj;
    static doublereal tem, critadd;
    static integer minspan;
    extern /* Subroutine */ int orthreg_();
    static doublereal critnew;
    static integer nterms21;
    static doublereal scr2;

    /* Parameter adjustments */
    tagx_dim1 = *nx;
    tagx_offset = tagx_dim1 + 1;
    tagx -= tagx_offset;
    x_dim1 = *nx;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --scrat;
    --w;
    --scr6;
    --ybar;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    betasc_dim1 = *mmax;
    betasc_offset = betasc_dim1 + 1;
    betasc -= betasc_offset;
    r_dim1 = *mmax;
    r_offset = r_dim1 + 1;
    r__ -= r_offset;
    bxsc_dim1 = *n;
    bxsc_offset = bxsc_dim1 + 1;
    bxsc -= bxsc_offset;
    --scr5;
    --scr1;
    flag_dim1 = *mmax;
    flag_offset = flag_dim1 + 1;
    flag__ -= flag_offset;
    --termlen;
    covsy_dim1 = *mmax;
    covsy_offset = covsy_dim1 + 1;
    covsy -= covsy_offset;
    cov_dim1 = *mmax;
    cov_offset = cov_dim1 + 1;
    cov -= cov_offset;
    --bxorthm;
    bxorth_dim1 = *n;
    bxorth_offset = bxorth_dim1 + 1;
    bxorth -= bxorth_offset;
    --tempin;
    bx_dim1 = *nx;
    bx_offset = bx_dim1 + 1;
    bx -= bx_offset;

    /* Function Body */
    *critmax = 0.;
    *jmax = 0;
    *imax = 0;
    *kmax = 0;
    i__1 = *nterms;
    for (m = 1; m <= i__1; ++m) {
	nm = 0;
	i__2 = *n;
	for (jjj = 1; jjj <= i__2; ++jjj) {
	    if (bx[jjj + m * bx_dim1] > 0.) {
		++nm;
	    }
/* L23152: */
	}
/* L23153: */
	tem = -(1. / (*n * nm)) * log(.94999999999999996);
	minspan = (integer) (log(tem) / log(2.) * -1. / (float)2.5);
	tem = .05 / *n;
	iendspan = (integer) (3. - log(tem) / log(2.));
	if (termlen[m] < *maxorder) {
	    i__2 = *p;
	    for (v = 1; v <= i__2; ++v) {
		if (flag__[m + v * flag_dim1] == 0) {
		    tnewform = TRUE_;
		    mm = 1;
L23162:
		    if (mm <= *nterms && tnewform) {
			++mm;
			if (tempin[mm] == 1) {
			    tnewform = FALSE_;
			    if (flag__[mm + v * flag_dim1] != 1) {
				tnewform = TRUE_;
				goto L9911;
			    }
			    i__3 = *p;
			    for (j = 1; j <= i__3; ++j) {
				if (j != v) {
				    if (flag__[mm + j * flag_dim1] != flag__[
					    m + j * flag_dim1]) {
					tnewform = TRUE_;
					goto L9911;
				    }
				}
/* L23168: */
			    }
/* L23169: */
			}
L9911:
			goto L23162;
		    }
/* L23163: */
		    if (tnewform) {
			nterms2 = *nterms + 1;
			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    scrat[i__] = x[i__ + v * x_dim1] * bx[i__ + m * 
				    bx_dim1];
/* L23176: */
			}
/* L23177: */
			if (*nterms > 1) {
			    orthreg_(n, n, nterms, &bxorth[bxorth_offset], &
				    tempin[1], &scrat[1], &bxorth[nterms2 * 
				    bxorth_dim1 + 1]);
			} else {
			    tem = 0.;
			    i__3 = *n;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				tem += scrat[i__] / *n;
/* L23180: */
			    }
/* L23181: */
			    i__3 = *n;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				bxorth[i__ + (bxorth_dim1 << 1)] = scrat[i__] 
					- tem;
/* L23182: */
			    }
/* L23183: */
			}
			bxorthm[nterms2] = (float)0.;
			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    bxorthm[nterms2] += bxorth[i__ + nterms2 * 
				    bxorth_dim1] / *n;
/* L23184: */
			}
/* L23185: */
			temp1 = (float)0.;
			i__3 = *n;
			for (i__ = 1; i__ <= i__3; ++i__) {
/* Computing 2nd power */
			    d__1 = bxorth[i__ + nterms2 * bxorth_dim1];
			    temp1 += d__1 * d__1;
/* L23186: */
			}
/* L23187: */
			if (temp1 > *tolbx) {
			    i__3 = *n;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				bxorth[i__ + nterms2 * bxorth_dim1] /= sqrt(
					temp1);
/* L23190: */
			    }
/* L23191: */
			} else {
			    i__3 = *n;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				bxorth[i__ + nterms2 * bxorth_dim1] = 0.;
/* L23192: */
			    }
/* L23193: */
			    tnewform = FALSE_;
			}
			i__3 = nterms2;
			for (i1 = 1; i1 <= i__3; ++i1) {
			    cov[i1 + nterms2 * cov_dim1] = (float)0.;
			    cov[nterms2 + i1 * cov_dim1] = (float)0.;
/* L23194: */
			}
/* L23195: */
			cov[nterms2 + nterms2 * cov_dim1] = 1.;
			i__3 = *nclass;
			for (kc = 1; kc <= i__3; ++kc) {
			    covsy[nterms2 + kc * covsy_dim1] = (float)0.;
			    i__4 = *n;
			    for (i__ = 1; i__ <= i__4; ++i__) {
				covsy[nterms2 + kc * covsy_dim1] += (y[i__ + 
					kc * y_dim1] - ybar[kc]) * bxorth[i__ 
					+ nterms2 * bxorth_dim1];
/* L23198: */
			    }
/* L23199: */
/* L23196: */
			}
/* L23197: */
			critnew = (float)0.;
			i__3 = *nclass;
			for (kc = 1; kc <= i__3; ++kc) {
			    temp1 = 0.;
			    i__4 = *n;
			    for (i__ = 1; i__ <= i__4; ++i__) {
				temp1 += y[i__ + kc * y_dim1] * bxorth[i__ + 
					nterms2 * bxorth_dim1];
/* L23202: */
			    }
/* L23203: */
/* Computing 2nd power */
			    d__1 = temp1;
			    critnew += d__1 * d__1;
/* L23200: */
			}
/* L23201: */
			if (critnew > *critmax) {
			    *jmax = v;
			    *critmax = critnew;
			    *imax = 1;
			    *kmax = m;
			}
		    }
		    if (tnewform) {
			nterms2 = *nterms + 1;
			nterms21 = *nterms + 2;
		    } else {
			nterms2 = *nterms;
			nterms21 = *nterms + 1;
			critnew = (float)0.;
		    }
		    i__3 = *nclass;
		    for (kc = 1; kc <= i__3; ++kc) {
			covsy[nterms21 + kc * covsy_dim1] = 0.;
/* L23208: */
		    }
/* L23209: */
		    i__3 = nterms21;
		    for (ii = 1; ii <= i__3; ++ii) {
			cov[ii + nterms21 * cov_dim1] = 0.;
			cov[nterms21 + ii * cov_dim1] = 0.;
/* L23210: */
		    }
/* L23211: */
		    i__3 = *nclass;
		    for (kc = 1; kc <= i__3; ++kc) {
			scr6[kc] = 0.;
/* L23212: */
		    }
/* L23213: */
		    i__3 = nterms21;
		    for (ii = 1; ii <= i__3; ++ii) {
			scr1[ii] = 0.;
/* L23214: */
		    }
/* L23215: */
		    scr2 = 0.;
		    su = 0.;
		    st = 0.;
		    sumbx2 = 0.;
		    sumb = (float)0.;
		    sumbx = (float)0.;
		    k = *n - 1;
L23216:
		    if (! (k > 0)) {
			goto L23218;
		    }
		    i__3 = nterms2;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			kk = tagx[k + v * tagx_dim1];
			kk1 = tagx[k + 1 + v * tagx_dim1];
			scr1[i__] += (bxorth[kk1 + i__ * bxorth_dim1] - 
				bxorthm[i__]) * bx[kk1 + m * bx_dim1];
			cov[i__ + nterms21 * cov_dim1] += (x[kk1 + v * x_dim1]
				 - x[kk + v * x_dim1]) * scr1[i__];
			cov[nterms21 + i__ * cov_dim1] = cov[i__ + nterms21 * 
				cov_dim1];
/* L23219: */
		    }
/* L23220: */
/* Computing 2nd power */
		    d__1 = bx[kk1 + m * bx_dim1];
		    scr2 += d__1 * d__1 * x[kk1 + v * x_dim1];
/* Computing 2nd power */
		    d__1 = bx[kk1 + m * bx_dim1];
		    sumbx2 += d__1 * d__1;
		    sumb += bx[kk1 + m * bx_dim1];
		    sumbx += bx[kk1 + m * bx_dim1] * x[kk1 + v * x_dim1];
		    su = st;
		    st = sumbx - sumb * x[kk + v * x_dim1];
		    cov[nterms21 + nterms21 * cov_dim1] = cov[nterms21 + 
			    nterms21 * cov_dim1] + (x[kk1 + v * x_dim1] - x[
			    kk + v * x_dim1]) * (scr2 * 2 - sumbx2 * (x[kk + 
			    v * x_dim1] + x[kk1 + v * x_dim1])) + (su * su - 
			    st * st) / *n;
		    crittemp = critnew;
		    i__3 = *nclass;
		    for (kc = 1; kc <= i__3; ++kc) {
			scr6[kc] += (y[kk1 + kc * y_dim1] - ybar[kc]) * bx[
				kk1 + m * bx_dim1];
			covsy[nterms21 + kc * covsy_dim1] += (x[kk1 + v * 
				x_dim1] - x[kk + v * x_dim1]) * scr6[kc];
			temp1 = covsy[nterms21 + kc * covsy_dim1];
			temp2 = cov[nterms21 + nterms21 * cov_dim1];
			i__4 = nterms2;
			for (jk = 1; jk <= i__4; ++jk) {
			    temp1 -= covsy[jk + kc * covsy_dim1] * cov[jk + 
				    nterms21 * cov_dim1];
			    temp2 -= cov[jk + nterms21 * cov_dim1] * cov[jk + 
				    nterms21 * cov_dim1];
/* L23223: */
			}
/* L23224: */
			if (cov[nterms21 + nterms21 * cov_dim1] > 0.) {
			    if (temp2 / cov[nterms21 + nterms21 * cov_dim1] > 
				    *tolbx) {
				critadd = temp1 * temp1 / temp2;
			    } else {
				critadd = (float)0.;
			    }
			} else {
			    critadd = 0.;
			}
			crittemp += critadd;
			if (crittemp > *rss * (float)1.01) {
			    crittemp = (float)0.;
			}
			if (crittemp > *prevcrit * 2) {
			    crittemp = (float)0.;
			}
/* L23221: */
		    }
/* L23222: */
		    if (k > 1) {
			k0 = tagx[k - 1 + v * tagx_dim1];
		    }
		    if (crittemp > *critmax && k % minspan == 0 && k >= 
			    iendspan && k <= *n - iendspan && bx[kk1 + m * 
			    bx_dim1] > 0. && ! (k > 1 && x[kk + v * x_dim1] ==
			     x[k0 + v * x_dim1])) {
			*jmax = v;
			*critmax = crittemp;
			*imax = k;
			*kmax = m;
			*newform = tnewform;
		    }
/* L23217: */
		    --k;
		    goto L23216;
L23218:
		    ;
		}
/* L9999: */
/* L23158: */
	    }
/* L23159: */
	}
/* L23150: */
    }
/* L23151: */
    return 0;
} /* addtrm_ */

