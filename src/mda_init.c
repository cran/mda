// Automatically generated, editing not advised.
#ifndef R_MDA_H
#define R_MDA_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("mda", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(bruto)(
double *x,
int *n,
int *q,
double *y,
int *p,
double *w,
double *knot,
int *nkmax,
int *nk,
double *wp,
int *match,
int *nef,
double *dfmax,
double *cost,
double *lambda,
double *df,
double *coef,
int *type,
double *xrange,
double *gcvsel,
double *gcvbak,
double *dfit,
int *maxit,
int *nit,
double *eta,
double *resid,
double *thresh,
double *work,
int *iwork,
int *trace
);
 
static R_NativePrimitiveArgType bruto_t[] = {
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP
};
void F77_SUB(marss)(
int *nx,
int *n,
int *p,
int *nclass,
double *y,
double *x,
double *w,
int *tagx,
int *maxorder,
int *mmax,
double *penalty,
double *thresh,
int *forwstep,
int *interms,
int *prune,
double *bx,
int *fullin,
int *lenb,
double *bestgcv,
int *bestin,
int *flag,
double *cut,
double *dir,
double *res,
double *beta,
double *scrat,
int *iscrat,
int *trace
);
 
static R_NativePrimitiveArgType marss_t[] = {
INTSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP
};
void F77_SUB(sspl0)(
double *x,
double *y,
double *w,
int *n,
int *p,
double *knot,
int *nk,
int *method,
double *tol,
double *wp,
int *match,
int *nef,
int *icen,
double *dfoff,
double *dfmax,
double *cost,
double *lambda,
double *df,
double *cv,
double *gcv,
double *coef,
double *s,
double *lev,
double *xrange,
double *work,
int *iwork,
int *ier
);
 
static R_NativePrimitiveArgType sspl0_t[] = {
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP
};
void F77_SUB(pbruto)(
double *x,
int *n,
int *q,
double *ybar,
int *p,
double *knot,
int *nkmax,
int *nk,
double *coef,
int *type,
double *xrange,
double *eta,
double *work
);
 
static R_NativePrimitiveArgType pbruto_t[] = {
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP
};
void F77_SUB(psspl2)(
double *x,
int *n,
int *p,
double *knot,
int *nk,
double *xrange,
double *coef,
double *coefl,
double *s,
int *order,
int *type
);
 
static R_NativePrimitiveArgType psspl2_t[] = {
REALSXP,
INTSXP,
INTSXP,
REALSXP,
INTSXP,
REALSXP,
REALSXP,
REALSXP,
REALSXP,
INTSXP,
INTSXP
};

static R_FortranMethodDef fMethods[] = {
FDEF(bruto) ,
FDEF(marss) ,
FDEF(sspl0) ,
FDEF(pbruto) ,
FDEF(psspl2) ,
{NULL, NULL, 0}
};

void R_init_mda(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
