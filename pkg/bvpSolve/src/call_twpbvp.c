#include <time.h>
#include <string.h>
#include "twpbvp.h"

/* definition of the calls to the fortran functions -
      subroutine twpbvp(ncomp, nlbc, aleft, aright,
     *       nfxpnt, fixpnt, ntol, ltol, tol,
     *       linear, givmsh, giveu, nmsh,
     *       xx, nudim, u, nmax,
     *       lwrkfl, wrk, lwrkin, iwrk,
     *       fsub, dfsub, gsub, dgsub, rpar, ipar, iflbvp)
*/

void F77_NAME(twpbvp)(int*, int*, double *, double *,
         int *, double *, int *, int *, double *,
         int *, int *, int *, int *,
         double *, int *, double *, int *,
         int *, double *, int*, int*,
         void (*)(int *, double *, double *, double *, double *, int *), /* fsub(n,x,u,f,rp,ip)   */
		     void (*)(int *, double *, double *, double *, double *, int *), /* dfsub(n,x,u,df,rp,ip) */
			   void (*)(int *, int *, double *, double *, double *, int *),    /* gsub(i,n,u,g,rp,ip)   */
		     void (*)(int *, int *, double *, double *, double *, int *),    /* dgsub(i,n,u,dg,rp,ip) */
         double *, int *, int *);

/* interface between fortran function calls and R functions
   Fortran code calls twp_derivs(m, x, y, ydot)
   R code called as twp_deriv_func(x, y) and returns ydot
   Note: passing of parameter values and "..." is done in R-function bvpcol*/

static void twp_derivs (int *n,  double *x, double *y, double *ydot,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                             REAL(X)[0]   = *x;
  for (i = 0; i < *n ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(twp_deriv_func,X,Y)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvp_envir));      incr_N_Protect();

  for (i = 0; i < *n ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
                                                my_unprotect(2);
}

static void forc_twp (int *neq, double *x, double *y,
                         double *ydot, double *rpar, int *ipar)
{
  updatedeforc(x);
  derfun(neq, x, y, ydot, rpar, ipar);
}

/* interface between fortran call to jacobian and R function */
static void twp_jac (int *n,  double *x, double *y, double *pd,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                           REAL(X)[0]   = *x;
  for (i = 0; i < *n; i++) REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(twp_jac_func,X,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvp_envir));      incr_N_Protect();

  for (i = 0; i < *n * *n; i++)  pd[i] = REAL(ans)[i];
                                                my_unprotect(2);
}

/* interface between fortran call to boundary condition and corresponding R function */

static void twp_bound (int *ii, int *n, double *y, double *gout,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < *n ; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(twp_bound_func,J,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvp_envir));        incr_N_Protect();
  /* only one element returned... */
  gout[0] = REAL(ans)[0];
                                                  my_unprotect(2);
}
/*interface between fortran call to jacobian of boundary and corresponding R function */

static void twp_jacbound (int *ii, int *n, double *y, double *dg,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                           INTEGER(J)[0] = *ii;
  for (i = 0; i < *n; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(twp_jacbound_func,J,Y));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvp_envir));          incr_N_Protect();

  for (i = 0; i < *n ; i++)  dg[i] = REAL(ans)[i];
                                                    my_unprotect(2);
}

/* number of eqs, order of eqs, summed order of eqns, from, to,
boundary points, settings, number of tolerances, tolerances,
mesh points, initial value of continuation parameter */

SEXP call_bvptwp(SEXP Ncomp, SEXP Nlbc, SEXP Fixpnt, SEXP Aleft, SEXP Aright,
		SEXP Tol, SEXP Linear, SEXP Givmesh, SEXP Givu, SEXP Nmesh,
		SEXP Nmax, SEXP Lwrkfl, SEXP Lwrkin, SEXP Xguess, SEXP Yguess,
    SEXP Rpar, SEXP Ipar, SEXP func, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP Initfunc, SEXP Parms, SEXP flist, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE;

  int  j, ii, ncomp, nlbc, nmax, lwrkfl, lwrkin, nx, *ipar, isForcing;
  double aleft, aright, *wrk, *tol, *fixpnt, *u, *xx, *rpar;
  int *ltol, *iwrk, ntol, iflag, nfixpnt, linear, givmesh, givu, nmesh, isDll;

  deriv_func    *derivs;
  jac_func      *jac;
  jacbound_func *jacbound;
  bound_func    *bound;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();

  ncomp  = INTEGER(Ncomp)[0];    /* number of equations */
  nlbc   = INTEGER(Nlbc)[0];     /* number of left boundary conditions */
  nmax   = INTEGER(Nmax)[0];     /* max number of mesh points */
  lwrkfl = INTEGER(Lwrkfl)[0];   /* length of double workspace */
  lwrkin = INTEGER(Lwrkin)[0];   /* length of integer workspace */
  linear = INTEGER(Linear)[0];   /* true if linear problem */
  givu   = INTEGER(Givu)[0];     /* true if initial trial solution given */
  givmesh = INTEGER(Givmesh)[0]; /* true if initial mesh given */
  nmesh  = INTEGER(Nmesh)[0];    /* size of mesh */

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

/* is function a dll ?*/
  if (inherits(func, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  ntol = LENGTH(Tol);
  tol   =(double *) R_alloc(ntol, sizeof(double));
    for (j = 0; j < ntol;j++) tol[j] = REAL(Tol)[j];

  ltol   =(int *) R_alloc(ntol, sizeof(int));
    for (j = 0; j < ntol;j++) ltol[j] = j+1;

  nfixpnt =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(nfixpnt, sizeof(double));
   for (j = 0; j < nfixpnt;j++) fixpnt[j] = REAL(Fixpnt)[j];

  xx   =(double *) R_alloc(nmax, sizeof(double));
   for (j = 0; j < nmesh; j++) xx[j] = REAL(Xguess)[j];
   for (j = nmesh; j < nmax; j++) xx[j] = 0;

  ii = nmax*ncomp;
  u   =(double *) R_alloc(ii, sizeof(double));
   for (j = 0; j < nmesh*ncomp; j++) u[j] = REAL(Yguess)[j];
   for (j = nmesh*ncomp; j < nmax*ncomp; j++) u[j] = 0;

  wrk = (double *) R_alloc(lwrkfl, sizeof(double));
  iwrk= (int *)    R_alloc(lwrkin, sizeof(int));

  ii = LENGTH(Ipar);
  ipar = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) ipar[j] = INTEGER(Ipar)[j];

  ii = LENGTH(Rpar);
  rpar = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rpar[j] = REAL(Rpar)[j];

/* initialise global R-variables... */
  if (isDll == 0) {
    PROTECT(X  = NEW_NUMERIC(1));               incr_N_Protect();
    PROTECT(J = NEW_INTEGER(1));                incr_N_Protect();
    PROTECT(Y = allocVector(REALSXP,ncomp));    incr_N_Protect();
  }

  
  isForcing = initForcings(flist);
  initParms(Initfunc, Parms);

   bvp_envir = rho;

  if (isDll) {   /* DLL addresses passed to fortran */
      derivs = (deriv_func *) R_ExternalPtrAddr(func);
      jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
      bound = (bound_func *) R_ExternalPtrAddr(boundfunc);
      jacbound = (jacbound_func *) R_ExternalPtrAddr(jacboundfunc);

	  /* here overruling derivs if forcing */
      if (isForcing) {
        derfun = (deriv_func *) R_ExternalPtrAddr(func);
        derivs = (deriv_func *) forc_twp;
      }

  } else {      /* interface functions between fortran and R */
      derivs = (deriv_func *) twp_derivs;
      twp_deriv_func = func;

      jac = twp_jac;
      twp_jac_func = jacfunc;

      bound = twp_bound;
      twp_bound_func = boundfunc;

      jacbound = twp_jacbound;
      twp_jacbound_func = jacboundfunc;
    }

/* Call the fortran function -
      subroutine twpbvp(ncomp, nlbc, aleft, aright,
     *       nfxpnt, fixpnt, ntol, ltol, tol,
     *       linear, givmsh, giveu, nmsh,
     *       xx, nudim, u, nmax,
     *       lwrkfl, wrk, lwrkin, iwrk,
     *       fsub, dfsub, gsub, dgsub, iflbvp)
  error("ncomp, nlbc, nfixpnt, ntol, atol, %f, %i, %i, %i, %i",tol[0],ncomp,nlbc,nfixpnt,ntol);
  error("till here %i, %i, %i, %i",ncomp,nlbc,nfixpnt,ntol);

*/
	  F77_CALL(twpbvp) (&ncomp, &nlbc, &aleft, &aright, &nfixpnt, fixpnt,
        &ntol, ltol, tol, &linear, &givmesh, &givu, &nmesh, xx, &ncomp,
        u, &nmax, &lwrkfl, wrk, &lwrkin, iwrk,
        derivs,jac,bound,jacbound, rpar, ipar,&iflag);
/*
C....   iflag - The Mode Of Return From twpbvp
C....         =  0  For Normal Return
C....         =  1  If The Expected No. Of Subintervals Exceeds Storage
C....               Specifications.
C....         = -1  If There Is An Input Data Error.
*/
	  if (iflag == -1)
     {
	   unprotect_all();
     error("One of the input parameters is invalid.\n");
     }
    else if (iflag == 1)
	{
	  unprotect_all();
	  error("The Expected No. Of Subintervals Exceeds Storage Specifications.\n");
	}

  else
	{
/*  warning("nmesh, nmax, ncomp, iflag %i, %i, %i, %i, %i, %i",nmesh,nmax,ncomp,iflag);
   error("TILL HERE");
*/
   nx = nmesh;

    PROTECT(yout = allocVector(REALSXP,(ncomp+1)*(nx)));incr_N_Protect();
	  for (j = 0; j < nx; j++)       REAL(yout)[j]    = xx[j];
    for (j = 0; j < ncomp*nx; j++) REAL(yout)[nx+j] =  u[j];
   }
  PROTECT(ISTATE = allocVector(INTSXP, 3));incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  INTEGER(ISTATE)[1] = nmax;
  INTEGER(ISTATE)[2] = nmesh;
  setAttrib(yout, install("istate"), ISTATE);

/*     ii = ncomp+7;

  ii = ispace[6];
  PROTECT(RWORK = allocVector(REALSXP, ii));incr_N_Protect();
  for (k = 0;k<ii;k++) REAL(RWORK)[k] = fspace[k];
  setAttrib(yout, install("rstate"), RWORK);
  }
                    ####   termination   ####                            */
  unprotect_all();
  return(yout);
}
