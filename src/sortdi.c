/* sortdi.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int sortdi_(v, lenv, a, ii, jj)
doublereal *v;
integer *lenv, *a, *ii, *jj;
{
    static integer i__, j, k, l, m, t, ij, il[20], iu[20], tt;
    static doublereal vt, vtt;


/*     puts into a the permutation vector which sorts v into */
/*     increasing order.  only elements from ii to jj are considered. */
/*     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements */
/*     v is returned sorted */
/*     this is a modification of cacm algorithm #347 by r. c. singleton, 
*/
/*     which is a modified hoare quicksort. */

    /* Parameter adjustments */
    --v;
    --a;

    /* Function Body */
    m = 1;
    i__ = *ii;
    j = *jj;
L10:
    if (i__ >= j) {
	goto L80;
    }
L20:
    k = i__;
    ij = (j + i__) / 2;
    t = a[ij];
    vt = v[ij];
    if (v[i__] <= vt) {
	goto L30;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    v[ij] = v[i__];
    v[i__] = vt;
    vt = v[ij];
L30:
    l = j;
    if (v[j] >= vt) {
	goto L50;
    }
    a[ij] = a[j];
    a[j] = t;
    t = a[ij];
    v[ij] = v[j];
    v[j] = vt;
    vt = v[ij];
    if (v[i__] <= vt) {
	goto L50;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    v[ij] = v[i__];
    v[i__] = vt;
    vt = v[ij];
    goto L50;
L40:
    a[l] = a[k];
    a[k] = tt;
    v[l] = v[k];
    v[k] = vtt;
L50:
    --l;
    if (v[l] > vt) {
	goto L50;
    }
    tt = a[l];
    vtt = v[l];
L60:
    ++k;
    if (v[k] < vt) {
	goto L60;
    }
    if (k <= l) {
	goto L40;
    }
    if (l - i__ <= j - k) {
	goto L70;
    }
    il[m - 1] = i__;
    iu[m - 1] = l;
    i__ = k;
    ++m;
    goto L90;
L70:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L90;
L80:
    --m;
    if (m == 0) {
	return 0;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L90:
    if (j - i__ > 10) {
	goto L20;
    }
    if (i__ == *ii) {
	goto L10;
    }
    --i__;
L100:
    ++i__;
    if (i__ == j) {
	goto L80;
    }
    t = a[i__ + 1];
    vt = v[i__ + 1];
    if (v[i__] <= vt) {
	goto L100;
    }
    k = i__;
L110:
    a[k + 1] = a[k];
    v[k + 1] = v[k];
    --k;
    if (vt < v[k]) {
	goto L110;
    }
    a[k + 1] = t;
    v[k + 1] = vt;
    goto L100;
} /* sortdi_ */

/* Subroutine */ int interv_(xt, lxt, x, left, mflag)
doublereal *xt;
integer *lxt;
doublereal *x;
integer *left, *mflag;
{
    /* Initialized data */

    static integer ilo = 1;

    static integer istep, middle, ihi;

/* omputes  left = max( i ; 1 .le. i .le. lxt  .and.  xt(i) .le. x )  . */

/* ******  i n p u t  ****** */
/* xt.....a double precision sequence, of length  lxt , assumed to be nond
ecreasing*/
/*  lxt.....number of terms in the sequence  xt . */
/*  x.....the point whose location with respect to the sequence  xt  is */
/*        to be determined. */

/* ******  o u t p u t  ****** */
/*  left, mflag.....both integers, whose value is */

/*   1     -1      if               x .lt.  xt(1) */
/*   i      0      if   xt(i)  .le. x .lt. xt(i+1) */
/*  lxt     1      if  xt(lxt) .le. x */

/*        in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0 */
/*        indicates that  x  lies outside the halfopen interval */
/*        xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the */
/*        interval is due to the decision to make all pp functions cont- 
*/
/*        inuous from the right. */

/* ******  m e t h o d  ****** */
/*  the program is designed to be efficient in the common situation that 
*/
/*  it is called repeatedly, with  x  taken from an increasing or decrea- 
*/
/*  sing sequence. this will happen, e.g., when a pp function is to be */
/*  graphed. the first guess for  left  is therefore taken to be the val- 
*/
/*  ue returned at the previous call and stored in the  l o c a l  varia- 
*/
/*  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec- 
*/
/*  essary since the present call may have nothing to do with the previ- 
*/
/*  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left = */
/*  ilo  and are done after just three comparisons. */
/*     otherwise, we repeatedly double the difference  istep = ihi - ilo 
*/
/*  while also moving  ilo  and  ihi  in the direction of  x , until */
/*                      xt(ilo) .le. x .lt. xt(ihi) , */
/*  after which we use bisection to get, in addition, ilo+1 = ihi . */
/*  left = ilo  is then returned. */

    /* Parameter adjustments */
    --xt;

    /* Function Body */
/*     save ilo  (a valid fortran statement in the new 1977 standard) */
    ihi = ilo + 1;
    if (ihi < *lxt) {
	goto L20;
    }
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    if (*lxt <= 1) {
	goto L90;
    }
    ilo = *lxt - 1;
    ihi = *lxt;

L20:
    if (*x >= xt[ihi]) {
	goto L40;
    }
    if (*x >= xt[ilo]) {
	goto L100;
    }

/*              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x . 
*/
/* L30: */
    istep = 1;
L31:
    ihi = ilo;
    ilo = ihi - istep;
    if (ilo <= 1) {
	goto L35;
    }
    if (*x >= xt[ilo]) {
	goto L50;
    }
    istep <<= 1;
    goto L31;
L35:
    ilo = 1;
    if (*x < xt[1]) {
	goto L90;
    }
    goto L50;
/*              **** now x .ge. xt(ihi) . increase  ihi  to capture  x . 
*/
L40:
    istep = 1;
L41:
    ilo = ihi;
    ihi = ilo + istep;
    if (ihi >= *lxt) {
	goto L45;
    }
    if (*x < xt[ihi]) {
	goto L50;
    }
    istep <<= 1;
    goto L41;
L45:
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    ihi = *lxt;

/*           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval. 
*/
L50:
    middle = (ilo + ihi) / 2;
    if (middle == ilo) {
	goto L100;
    }
/*     note. it is assumed that middle = ilo in case ihi = ilo+1 . */
    if (*x < xt[middle]) {
	goto L53;
    }
    ilo = middle;
    goto L50;
L53:
    ihi = middle;
    goto L50;
/* **** set output and return. */
L90:
    *mflag = -1;
    *left = 1;
    return 0;
L100:
    *mflag = 0;
    *left = ilo;
    return 0;
L110:
    *mflag = 1;
    *left = *lxt;
    return 0;
} /* interv_ */

/* Subroutine */ int bsplvd_(t, k, x, left, a, dbiatx, nderiv)
doublereal *t;
integer *k;
doublereal *x;
integer *left;
doublereal *a, *dbiatx;
integer *nderiv;
{
    /* System generated locals */
    integer a_dim1, a_offset, dbiatx_dim1, dbiatx_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer jlow, kp1mm, i__, j, m, mhigh, jp1mid;
    static doublereal fkp1mm;
    static integer il;
    static doublereal factor;
    static integer ideriv;
    extern /* Subroutine */ int bsplvb_();
    static integer ldummy, kp1;
    static doublereal sum;

/* alls bsplvb */
/* alculates value and deriv.s of all b-splines which do not vanish at x 
*/

/* ******  i n p u t  ****** */
/*  t     the knot array, of length left+k (at least) */
/*  k     the order of the b-splines to be evaluated */
/*  x     the point at which these values are sought */
/*  left  an integer indicating the left endpoint of the interval of */
/*        interest. the  k  b-splines whose support contains the interval 
*/
/*               (t(left), t(left+1)) */
/*        are to be considered. */
/*  a s s u m p t i o n  - - -  it is assumed that */
/*               t(left) .lt. t(left+1) */
/*        division by zero will result otherwise (in  b s p l v b ). */
/*        also, the output is as advertised only if */
/*               t(left) .le. x .le. t(left+1) . */
/*  nderiv   an integer indicating that values of b-splines and their */
/*        derivatives up to but not including the  nderiv-th  are asked */
/*        for. ( nderiv  is replaced internally by the integer in (1,k) */
/*        closest to it.) */

/* ******  w o r k   a r e a  ****** */
/*  a     an array of order (k,k), to contain b-coeff.s of the derivat- */
/*        ives of a certain order of the  k  b-splines of interest. */

/* ******  o u t p u t  ****** */
/*  dbiatx   an array of order (k,nderiv). its entry  (i,m)  contains */
/*        value of  (m-1)st  derivative of  (left-k+i)-th  b-spline of */
/*        order  k  for knot sequence  t , i=m,...,k; m=1,...,nderiv. */

/* ******  m e t h o d  ****** */
/*  values at  x  of all the relevant b-splines of order k,k-1,..., */
/*  k+1-nderiv  are generated via  bsplvb  and stored temporarily */
/*  in  dbiatx .  then, the b-coeffs of the required derivatives of the */
/*  b-splines of interest are generated by differencing, each from the */
/*  preceding one of lower order, and combined with the values of b- */
/*  splines of corresponding order in  dbiatx  to produce the desired */
/*  values. */

    /* Parameter adjustments */
    --t;
    a_dim1 = *k;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    dbiatx_dim1 = *k;
    dbiatx_offset = dbiatx_dim1 + 1;
    dbiatx -= dbiatx_offset;

    /* Function Body */
/* Computing MAX */
    i__1 = min(*nderiv,*k);
    mhigh = max(i__1,1);
/*     mhigh is usually equal to nderiv. */
    kp1 = *k + 1;
    i__1 = kp1 - mhigh;
    bsplvb_(&t[1], &i__1, &c__1, x, left, &dbiatx[dbiatx_offset]);
    if (mhigh == 1) {
	goto L99;
    }
/*     the first column of  dbiatx  always contains the b-spline values */
/*     for the current order. these are stored in column k+1-current */
/*     order  before  bsplvb  is called to put values for the next */
/*     higher order on top of it. */
    ideriv = mhigh;
    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	jp1mid = 1;
	i__2 = *k;
	for (j = ideriv; j <= i__2; ++j) {
	    dbiatx[j + ideriv * dbiatx_dim1] = dbiatx[jp1mid + dbiatx_dim1];
/* L11: */
	    ++jp1mid;
	}
	--ideriv;
	i__2 = kp1 - ideriv;
	bsplvb_(&t[1], &i__2, &c__2, x, left, &dbiatx[dbiatx_offset]);
/* L15: */
    }

/*     at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j) for */
/*     i=j,...,k and j=1,...,mhigh ('=' nderiv). in particular, the */
/*     first column of  dbiatx  is already in final form. to obtain cor- 
*/
/*     responding derivatives of b-splines in subsequent columns, gene- */
/*     rate their b-repr. by differencing, then evaluate at  x. */

    jlow = 1;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *k;
	for (j = jlow; j <= i__2; ++j) {
/* L19: */
	    a[j + i__ * a_dim1] = 0.;
	}
	jlow = i__;
/* L20: */
	a[i__ + i__ * a_dim1] = 1.;
    }
/*     at this point, a(.,j) contains the b-coeffs for the j-th of the */
/*     k  b-splines of interest here. */

    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	kp1mm = kp1 - m;
	fkp1mm = (doublereal) kp1mm;
	il = *left;
	i__ = *k;

/*        for j=1,...,k, construct b-coeffs of  (m-1)st  derivative of
 */
/*        b-splines from those for preceding derivative by differencin
g */
/*        and store again in  a(.,j) . the fact that  a(i,j) = 0  for 
*/
/*        i .lt. j  is used.sed. */
	i__2 = kp1mm;
	for (ldummy = 1; ldummy <= i__2; ++ldummy) {
	    factor = fkp1mm / (t[il + kp1mm] - t[il]);
/*           the assumption that t(left).lt.t(left+1) makes denomi
nator */
/*           in  factor  nonzero. */
	    i__3 = i__;
	    for (j = 1; j <= i__3; ++j) {
/* L24: */
		a[i__ + j * a_dim1] = (a[i__ + j * a_dim1] - a[i__ - 1 + j * 
			a_dim1]) * factor;
	    }
	    --il;
/* L25: */
	    --i__;
	}

/*        for i=1,...,k, combine b-coeffs a(.,i) with b-spline values 
*/
/*        stored in dbiatx(.,m) to get value of  (m-1)st  derivative o
f */
/*        i-th b-spline (of interest here) at  x , and store in */
/*        dbiatx(i,m). storage of this value over the value of a b-spl
ine */
/*        of order m there is safe since the remaining b-spline deriva
t- */
/*        ive of the same order do not use this value due to the fact 
*/
/*        that  a(j,i) = 0  for j .lt. i . */
/* L30: */
	i__2 = *k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = (float)0.;
	    jlow = max(i__,m);
	    i__3 = *k;
	    for (j = jlow; j <= i__3; ++j) {
/* L35: */
		sum = a[j + i__ * a_dim1] * dbiatx[j + m * dbiatx_dim1] + sum;
	    }
/* L40: */
	    dbiatx[i__ + m * dbiatx_dim1] = sum;
	}
    }
L99:
    return 0;
} /* bsplvd_ */

doublereal bvalue_(t, bcoef, n, k, x, jderiv)
doublereal *t, *bcoef;
integer *n, *k;
doublereal *x;
integer *jderiv;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static doublereal fkmj;
    static integer i__, j, mflag, jcmin, jcmax, jdrvp1;
    static doublereal aj[20];
    static integer jc;
    static doublereal dm[20];
    static integer jj;
    static doublereal dp[20];
    extern /* Subroutine */ int interv_();
    static integer km1, imk, kmj, ilo, nmi;

/* alls  interv */

/* alculates value at  x  of  jderiv-th derivative of spline from b-repr. 
*/
/*  the spline is taken to be continuous from the right. */

/* ******  i n p u t ****** */
/*  t, bcoef, n, k......forms the b-representation of the spline  f  to */
/*        be evaluated. specifically, */
/*  t.....knot sequence, of length  n+k, assumed nondecreasing. */
/*  bcoef.....b-coefficient sequence, of length  n . */
/*  n.....length of  bcoef  and dimension of s(k,t), */
/*        a s s u m e d  positive . */
/*  k.....order of the spline . */

/*  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed 
*/
/*        arbitrarily by the dimension statement for  aj, dm, dm  below, 
*/
/*        but is  n o w h e r e  c h e c k e d  for. */

/*  x.....the point at which to evaluate . */
/*  jderiv.....integer giving the order of the derivative to be evaluated 
*/
/*        a s s u m e d  to be zero or positive. */

/* ******  o u t p u t  ****** */
/*  bvalue.....the value of the (jderiv)-th derivative of  f  at  x . */

/* ******  m e t h o d  ****** */
/*     the nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo- 
*/
/*  cated with the aid of  interv . the  k  b-coeffs of  f  relevant for 
*/
/*  this interval are then obtained from  bcoef (or taken to be zero if */
/*  not explicitly available) and are then differenced  jderiv  times to 
*/
/*  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval. */
/*  precisely, with  j = jderiv, we have from x.(12) of the text that */

/*     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) ) */

/*  where */
/*                   / bcoef(.),                     ,  j .eq. 0 */
/*                   / */
/*    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1) */
/*                   / ----------------------------- ,  j .gt. 0 */
/*                   /    (t(.+k-j) - t(.))/(k-j) */

/*     then, we use repeatedly the fact that */

/*    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) ) */
/*  with */
/*                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1) */
/*    a(.,x)  =    --------------------------------------- */
/*                 (x - t(.))      + (t(.+m-1) - x) */

/*  to write  (d**j)f(x)  eventually as a linear combination of b-splines 
*/
/*  of order  1 , and the coefficient for  b(i,1,t)(x)  must then */
/*  be the desired number  (d**j)f(x). (see x.(17)-(19) of text). */

/*     dimension t(n+k) */
/* urrent fortran standard makes it impossible to specify the length of */
/*  t  precisely without the introduction of otherwise superfluous */
/*  additional arguments. */
    /* Parameter adjustments */
    --t;
    --bcoef;

    /* Function Body */
    ret_val = (float)0.;
    if (*jderiv >= *k) {
	goto L99;
    }

/*  *** find  i  s.t.  1 .le. i .lt. n+k  and  t(i) .lt. t(i+1) and */
/*      t(i) .le. x .lt. t(i+1) . if no such i can be found,  x  lies */
/*      outside the support of  the spline  f  and  bvalue = 0. */
/*      (the asymmetry in this choice of  i  makes  f  rightcontinuous) */
    if (*x != t[*n + 1] || t[*n + 1] != t[*n + *k]) {
	goto L700;
    }
    i__ = *n;
    goto L701;
L700:
    i__1 = *n + *k;
    interv_(&t[1], &i__1, x, &i__, &mflag);
    if (mflag != 0) {
	goto L99;
    }
L701:
/*  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i). */
    km1 = *k - 1;
    if (km1 > 0) {
	goto L1;
    }
    ret_val = bcoef[i__];
    goto L99;

/*  *** store the k b-spline coefficients relevant for the knot interval 
*/
/*     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dm(j) = x - t(i+1-j), 
*/
/*     dp(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable 
*/
/*     from input to zero. set any t.s not obtainable equal to t(1) or */
/*     to t(n+k) appropriately. */
L1:
    jcmin = 1;
    imk = i__ - *k;
    if (imk >= 0) {
	goto L8;
    }
    jcmin = 1 - imk;
    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
/* L5: */
	dm[j - 1] = *x - t[i__ + 1 - j];
    }
    i__1 = km1;
    for (j = i__; j <= i__1; ++j) {
	aj[*k - j - 1] = (float)0.;
/* L6: */
	dm[j - 1] = dm[i__ - 1];
    }
    goto L10;
L8:
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
/* L9: */
	dm[j - 1] = *x - t[i__ + 1 - j];
    }

L10:
    jcmax = *k;
    nmi = *n - i__;
    if (nmi >= 0) {
	goto L18;
    }
    jcmax = *k + nmi;
    i__1 = jcmax;
    for (j = 1; j <= i__1; ++j) {
/* L15: */
	dp[j - 1] = t[i__ + j] - *x;
    }
    i__1 = km1;
    for (j = jcmax; j <= i__1; ++j) {
	aj[j] = (float)0.;
/* L16: */
	dp[j - 1] = dp[jcmax - 1];
    }
    goto L20;
L18:
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
/* L19: */
	dp[j - 1] = t[i__ + j] - *x;
    }

L20:
    i__1 = jcmax;
    for (jc = jcmin; jc <= i__1; ++jc) {
/* L21: */
	aj[jc - 1] = bcoef[imk + jc];
    }

/*               *** difference the coefficients  jderiv  times. */
    if (*jderiv == 0) {
	goto L30;
    }
    i__1 = *jderiv;
    for (j = 1; j <= i__1; ++j) {
	kmj = *k - j;
	fkmj = (doublereal) kmj;
	ilo = kmj;
	i__2 = kmj;
	for (jj = 1; jj <= i__2; ++jj) {
	    aj[jj - 1] = (aj[jj] - aj[jj - 1]) / (dm[ilo - 1] + dp[jj - 1]) * 
		    fkmj;
/* L23: */
	    --ilo;
	}
    }

/*  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative, */
/*     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv). */
L30:
    if (*jderiv == km1) {
	goto L39;
    }
    jdrvp1 = *jderiv + 1;
    i__2 = km1;
    for (j = jdrvp1; j <= i__2; ++j) {
	kmj = *k - j;
	ilo = kmj;
	i__1 = kmj;
	for (jj = 1; jj <= i__1; ++jj) {
	    aj[jj - 1] = (aj[jj] * dm[ilo - 1] + aj[jj - 1] * dp[jj - 1]) / (
		    dm[ilo - 1] + dp[jj - 1]);
/* L33: */
	    --ilo;
	}
    }
L39:
    ret_val = aj[0];

L99:
    return ret_val;
} /* bvalue_ */

/* Subroutine */ int bsplvb_(t, jhigh, index, x, left, biatx)
doublereal *t;
integer *jhigh, *index;
doublereal *x;
integer *left;
doublereal *biatx;
{
    /* Initialized data */

    static integer j = 1;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal term;
    static integer i__;
    static doublereal saved, deltal[20], deltar[20];
    static integer jp1;

/* alculates the value of all possibly nonzero b-splines at  x  of order 
*/

/*               jout  =  dmax( jhigh , (j+1)*(index-1) ) */

/*  with knot sequence  t . */

/* ******  i n p u t  ****** */
/*  t.....knot sequence, of length  left + jout  , assumed to be nonde- */
/*        creasing.  a s s u m p t i o n . . . . */
/*                       t(left)  .lt.  t(left + 1)   . */
/*   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1) */
/*  jhigh, */
/*  index.....integers which determine the order  jout = max(jhigh, */
/*        (j+1)*(index-1))  of the b-splines whose values at  x  are to */
/*        be returned.  index  is used to avoid recalculations when seve- 
*/
/*        ral columns of the triangular array of b-spline values are nee- 
*/
/*        ded (e.g., in  bvalue  or in  bsplvd ). precisely, */
/*                     if  index = 1 , */
/*        the calculation starts from scratch and the entire triangular */
/*        array of b-spline values of orders 1,2,...,jhigh  is generated 
*/
/*        order by order , i.e., column by column . */
/*                     if  index = 2 , */
/*        only the b-spline values of order  j+1, j+2, ..., jout  are ge- 
*/
/*        nerated, the assumption being that  biatx , j , deltal , deltar 
*/
/*        are, on entry, as they were on exit at the previous call. */
/*           in particular, if  jhigh = 0, then  jout = j+1, i.e., just */
/*        the next column of b-spline values is generated. */

/*  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im- 
*/
/*        posed arbitrarily by the dimension statement for  deltal  and */
/*        deltar  below, but is  n o w h e r e  c h e c k e d  for . */

/*  x.....the point at which the b-splines are to be evaluated. */
/*  left.....an integer chosen (usually) so that */
/*                  t(left) .le. x .le. t(left+1)  . */

/* ******  o u t p u t  ****** */
/*  biatx.....array of length  jout , with  biatx(i)  containing the val- 
*/
/*        ue at  x  of the polynomial of order  jout  which agrees with */
/*        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left), 
*/
/*        t(left+1)) . */

/* ******  m e t h o d  ****** */
/*  the recurrence relation */

/*                       x - t(i)              t(i+j+1) - x */
/*     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x) 
*/
/*                     t(i+j)-t(i)            t(i+j+1)-t(i+1) */

/*  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x), 
*/
/*  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),..., */
/*  b(left,j)(x), storing the new values in  biatx  over the old. the */
/*  facts that */
/*            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1) */
/*  and that */
/*            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j) */
/*  are used. the particular organization of the calculations follows al- 
*/
/*  gorithm  (8)  in chapter x of the text. */

/*     dimension biatx(jout), t(left+jout) */
/* urrent fortran standard makes it impossible to specify the length of */
/*  t  and of  biatx  precisely without the introduction of otherwise */
/*  superfluous additional arguments. */
    /* Parameter adjustments */
    --t;
    --biatx;

    /* Function Body */
/*     save j,deltal,deltar (valid in fortran 77) */

    switch ((int)*index) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    j = 1;
    biatx[1] = 1.;
    if (j >= *jhigh) {
	goto L99;
    }

L20:
    jp1 = j + 1;
    deltar[j - 1] = t[*left + j] - *x;
    deltal[j - 1] = *x - t[*left + 1 - j];
    saved = 0.;
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
	term = biatx[i__] / (deltar[i__ - 1] + deltal[jp1 - i__ - 1]);
	biatx[i__] = saved + deltar[i__ - 1] * term;
/* L26: */
	saved = deltal[jp1 - i__ - 1] * term;
    }
    biatx[jp1] = saved;
    j = jp1;
    if (j < *jhigh) {
	goto L20;
    }

L99:
    return 0;
} /* bsplvb_ */

/* Subroutine */ int dpbfa_(abd, lda, n, m, info)
doublereal *abd;
integer *lda, *n, *m, *info;
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal ddot_();
    static integer j, k;
    static doublereal s, t;
    static integer ik, jk, mu;


/*     dpbfa factors a double precision symmetric positive definite */
/*     matrix stored in band form. */

/*     dpbfa is usually called by dpbco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                the matrix to be factored.  the columns of the upper */
/*                triangle are stored in the columns of abd and the */
/*                diagonals of the upper triangle are stored in the */
/*                rows of abd .  see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. m + 1 . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */
/*                0 .le. m .lt. n . */

/*     on return */

/*        abd     an upper triangular matrix  r , stored in band */
/*                form, so that  a = trans(r)*r . */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  if the leading minor of order  k  is not */
/*                     positive definite. */

/*     band storage */

/*           if  a  is a symmetric positive definite band matrix, */
/*           the following program segment will set up the input. */

/*                   m = (band width above diagonal) */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-m) */
/*                      do 10 i = i1, j */
/*                         k = i-j+m+1 */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas ddot */
/*     fortran max0,dsqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = abd_dim1 + 1;
    abd -= abd_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	ik = *m + 1;
/* Computing MAX */
	i__2 = j - *m;
	jk = max(i__2,1);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (k = mu; k <= i__2; ++k) {
	    i__3 = k - mu;
	    t = abd[k + j * abd_dim1] - ddot_(&i__3, &abd[ik + jk * abd_dim1],
		     &c__1, &abd[mu + j * abd_dim1], &c__1);
	    t /= abd[*m + 1 + jk * abd_dim1];
	    abd[k + j * abd_dim1] = t;
	    s += t * t;
	    --ik;
	    ++jk;
/* L10: */
	}
L20:
	s = abd[*m + 1 + j * abd_dim1] - s;
/*     ......exit */
	if (s <= 0.) {
	    goto L40;
	}
	abd[*m + 1 + j * abd_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* dpbfa_ */

/* Subroutine */ int dpbsl_(abd, lda, n, m, b)
doublereal *abd;
integer *lda, *n, *m;
doublereal *b;
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;

    /* Local variables */
    extern doublereal ddot_();
    static integer k;
    static doublereal t;
    extern /* Subroutine */ int daxpy_();
    static integer kb, la, lb, lm;


/*     dpbsl solves the double precision symmetric positive definite */
/*     band system  a*x = b */
/*     using the factors computed by dpbco or dpbfa. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                the output from dpbco or dpbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  technically this indicates */
/*        singularity but it is usually caused by improper subroutine */
/*        arguments.  it will not occur if the subroutines are called */
/*        correctly and  info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dpbco(abd,lda,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call dpbsl(abd,lda,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran min0 */

/*     internal variables */


/*     solve trans(r)*y = b */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = abd_dim1 + 1;
    abd -= abd_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = ddot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	b[k] = (b[k] - t) / abd[*m + 1 + k * abd_dim1];
/* L10: */
    }

/*     solve r*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	b[k] /= abd[*m + 1 + k * abd_dim1];
	t = -b[k];
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* dpbsl_ */

