/* dorthreg.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int orthreg_(nx, n, p, x, in, y, res)
integer *nx, *n, *p;
doublereal *x;
integer *in;
doublereal *y, *res;
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static doublereal beta, temp1, temp2;
    static integer i__, j;

    /* Parameter adjustments */
    --res;
    --y;
    --in;
    x_dim1 = *nx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	res[i__] = y[i__];
/* L23000: */
    }
/* L23001: */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	if (in[j] == 1) {
	    temp1 = 0.;
	    temp2 = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp1 += res[i__] * x[i__ + j * x_dim1];
		temp2 += x[i__ + j * x_dim1] * x[i__ + j * x_dim1];
/* L23006: */
	    }
/* L23007: */
	    beta = temp1 / temp2;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		res[i__] -= beta * x[i__ + j * x_dim1];
/* L23008: */
	    }
/* L23009: */
	}
/* L23002: */
    }
/* L23003: */
    return 0;
} /* orthreg_ */

