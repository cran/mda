/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* S compatibility library - mapping some internal functions in S to R
 *
 * $Id: S_compat.c,v 1.1 1998/01/07 23:40:13 bates Exp $
 */

#include "S.h"


/* Not needed if we could include "Linpack.h" but that doesn't work
 * right now.  Linpack.h needs Fortran.h needs Platform.h which is not 
 * available.  Platform.h should be copied to $RHOME/include
 */

extern longint
F77_CALL(dqrdc2)(double*, longint*, longint*, longint*, double*,
		 longint*, double*, longint*, double*);

void
F77_NAME(dqrdca) (double *x, longint *ldx, longint *n, longint *p,
		  double *qraux, longint *jpvt, double *work,
		  longint *rank, double *tol)
{
  F77_CALL(dqrdc2) (x, ldx, n, p, tol, rank, qraux, jpvt, work);
}

/* #include<Linpack.h> */
/* #include<S_compat.h>*/

void dblepr_(char *label, int nchar, double *data, int ndata);
void intpr_(char *label, int nchar, int *data, int ndata);

void dblepr_(char *label, int nchar, double *data, int ndata)
{
  int k;
  for(k=0;k<nchar; k++)
    putchar(label[k]);
  putchar('\n');

  for(k=0;k<ndata; k++)
    printf("[%d] %lf\n", k, data[k]);
}
  
void intpr_(char *label, int nchar, int *data, int ndata)
{
  int k;
  for(k=0;k<nchar; k++)
    putchar(label[k]);
  putchar('\n');

  for(k=0;k<ndata; k++)
    printf("[%d] %d\n", k, data[k]);
}
  


