#include<R_ext/Rdynload.h>
#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#include "de.h"

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/*
    .C	R_CMethodDef
    .Call	R_CallMethodDef
    .Fortran	R_FortranMethodDef
    .External	R_ExternalMethodDef
*/

/* .Call calls */
extern SEXP call_mebdfi(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_gambim(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_dop(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[]  = {
  {"call_mebdfi", (DL_FUNC) &call_mebdfi, 27},
  {"call_gambim", (DL_FUNC) &call_gambim, 25},
  {"call_dop",    (DL_FUNC) &call_dop   , 18},
  {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_deTestSet(DllInfo *dll) {

  // register entry points
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);

  R_useDynamicSymbols(dll, TRUE); // enable dynamic searching

}
