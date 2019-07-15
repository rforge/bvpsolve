#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#include "bvpSolve.h"

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */

/*
   ToDo:
   - consider replacing SEXP with REALSXP, INTSXP, STRSXP (character), VEXSXP (lists) etc.
   - unlock
*/

/* .C calls */
/* extern void unlock_solver();*/

extern void fsub(int *n, double *x, double *z, double *f,
                 double * RPAR, int * IPAR);
                 
extern void dfsub(int * n, double *x, double *z, double * df,
      double *RPAR, int *IPAR);  
                 
extern void gsub(int *i, int *n, double *z, double *g,
      double *RPAR, int *IPAR);  
      
extern void dgsub(int *i, int *n, double *z, double *dg,
      double *RPAR, int *IPAR);  

static const R_CMethodDef CEntries[] = {
/*    {"unlock_solver", (DL_FUNC) &unlock_solver, 0},*/
    {"fsub",          (DL_FUNC) &fsub,          6},
    {"dfsub",         (DL_FUNC) &dfsub,         6},
    {"gsub",          (DL_FUNC) &gsub,          6},
    {"dgsub",         (DL_FUNC) &dgsub,         6},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP call_acdc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_colnew(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_bvptwp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_colmodsys(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"call_acdc",       (DL_FUNC) &call_acdc,      31},
    {"call_colnew",     (DL_FUNC) &call_colnew,    26},
    {"call_bvptwp",     (DL_FUNC) &call_bvptwp,    30},
    {"call_colmodsys",  (DL_FUNC) &call_colmodsys, 28},
    {NULL, NULL, 0}
};

/* .Fortran calls */
void F77_NAME(appsln)(double *, double *, double *, int *);
void F77_NAME(sysappsln)(double *, double *, double *, int *);    
void F77_NAME(mappsln)(double *, double *, double *, int *);

void F77_NAME(initbnd)(void *);
void F77_NAME(funbnd)(int *, double *, double *, double *,
                 double * , int * );

void F77_NAME(dfbnd)(int * , double *, double *, double *,
                  double *, int *);  

void F77_NAME(gbnd)(int *, int *, double *, double *,
                 double *, int *);  

void F77_NAME(dgbnd)(int *, int *, double *, double *,
                  double *, int *);  
 
R_FortranMethodDef FEntries[] = {
    {"mappsln",       (DL_FUNC) &F77_SUB(mappsln),      4},
    {"appsln",        (DL_FUNC) &F77_SUB(appsln),       4},
    {"sysappsln",     (DL_FUNC) &F77_SUB(sysappsln),    4},
    {"initbnd",       (DL_FUNC) &F77_SUB(initbnd),      1},
    {"funbnd",        (DL_FUNC) &F77_SUB(funbnd),       6},
    {"dfbnd",         (DL_FUNC) &F77_SUB(dfbnd),        6},
    {"gbnd",          (DL_FUNC) &F77_SUB(gbnd),         6},
    {"dgbnd",         (DL_FUNC) &F77_SUB(dgbnd),        6},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_bvpSolve(DllInfo *dll) {

  R_registerRoutines(dll, CEntries, CallEntries, FEntries, NULL);

  // the following two lines protect against accidentially finding entry points

  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}
