#include <R.h>
#include <Rdefines.h>

/* give name to data types */
typedef void deriv_func    (int *, double *,double *, double *, double *, int *);
typedef void bound_func    (int *, int *, double *,double *, double *, int *);
typedef void jac_func      (int *, double *, double *,double *, double *, int *);
typedef void jacbound_func (int *, int *, double *, double *, double *, int *);

deriv_func * derfun;

/* global variables */
SEXP X, Y, J;
extern SEXP bvp_gparms;

/* bvp globals */
extern SEXP twp_deriv_func;
extern SEXP twp_jac_func;
extern SEXP twp_bound_func;
extern SEXP twp_jacbound_func;
extern SEXP bvp_envir;

extern SEXP bvp_envir;
/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initibvpparms;*/
void Initbvpparms(int *, double *);
void initParms(SEXP Initfunc, SEXP Parms);
typedef void init_func (void (*)(int *, double *));

/* forcings */
void updatedeforc(double *);
void Initdeforc(int *, double *);

int initForcings(SEXP list);

SEXP getListElement(SEXP list, const char *str);

long int nforc;  /* the number of forcings */

/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
double * tvec;
double * fvec;
int    * ivec;
int    fmethod;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
int    * findex;
double * intpol;
int    * maxindex;

double * forcings;
