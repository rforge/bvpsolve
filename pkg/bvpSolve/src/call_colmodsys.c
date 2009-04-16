#include <time.h>
#include <string.h>
#include "colmod.h"

/* definition of the calls to the fortran functions -
      Subroutine Colmod(Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)
     
     
          Call Appsln(Xx,Z,Fspace,Ispace)               */
void F77_NAME(colmod)(int*, int*, double *, double *, double *, int *, int *,
         double *, double *, int *, double *, int *, double *, double *,
         void (*)(double *, double *, double *, double *),        /* fsub  */
		     void (*)(double *, double *, double *, int *, double *), /* dfsub */
			   void (*)(int *, double *, double *, double *),           /* gsub  */
		     void (*)(int *, double *, double *, double *),           /* dgsub */
         void (*)(double *, double *, double *, double *));       /* guess */

void F77_NAME(appsln)(double *, double *, double *, int *);

/* interface between fortran function calls and R functions
   Fortran code calls colsys_derivs(x, y, ydot, eps)
   R code called as colsys_deriv_func(time, y) and returns ydot
   Note: passing of parameter values and "..." is done in R-function bvpcol*/
static void colsys_guess (double *x, double *y,
                        double *ydot, double *eps)
{
}

static void colsys_derivs (double *x, double *y,
                        double *ydot, double *eps)
{
  int i;
  SEXP R_fcall, ans;
                                REAL(EPS)[0] = *eps;
                                REAL(X)[0]   = *x;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang4(colsys_deriv_func,X,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvpcolmod_envir));      incr_N_Protect();

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

/* interface between fortran call to jacobian and R function */
static void colsys_jac (double *x, double *y, double *pd, int *neq, double *eps)
{
  int i;
  SEXP R_fcall, ans;
                              REAL(EPS)[0] = *eps;
                              REAL(X)[0]   = *x;
  for (i = 0; i < mstar; i++) REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang4(colsys_jac_func,X,Y,EPS));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvpcolmod_envir));      incr_N_Protect();

  for (i = 0; i < *neq * mstar; i++)  pd[i] = REAL(ans)[i];
                                                my_unprotect(2);
}

/* interface between fortran call to boundary condition and corresponding R function */

static void colsys_bound (int *ii, double *y, double *gout, double *eps)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(EPS)[0]  = *eps;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(colsys_bound_func,J,Y,EPS));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvpcolmod_envir));        incr_N_Protect();
  /* only on e element returned... */
   gout[0] = REAL(ans)[0];
                                                  my_unprotect(2);
}
/*interface between fortran call to jacobian of boundary and corresponding R function */

static void colsys_jacbound (int *ii, double *y, double *dg, double *eps)
{
  int i;
  SEXP R_fcall, ans;
                             REAL(EPS)[0]  = *eps;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < mstar; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(colsys_jacbound_func,J,Y,EPS));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, bvpcolmod_envir));          incr_N_Protect();

  for (i = 0; i < mstar ; i++)  dg[i] = REAL(ans)[i];
                                                    my_unprotect(2);
}

/* give name to data types */
typedef void deriv_func2(double *, double *,double *, double *);
typedef void bound_func2 (int *, double *, double *,double *);
typedef void jac_func2  ( double *, double *, double *, int *, double *);
typedef void jacbound_func2   (int *, double *, double *, double *);
/* typedef void init_func (void (*)(int *, double *));*/



/* MAIN C-FUNCTION, CALLED FROM R-code

      Subroutine Colmod(Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)             */

/* number of eqs, order of eqs, summed order of eqns, from, to,
boundary points, settings, number of tolerances, tolerances,
mesh points, initial value of continuation parameter */

SEXP call_colmodsys(SEXP Ncomp, SEXP Mstar, SEXP M, SEXP Xout, SEXP Aleft, SEXP Aright,
		SEXP Zeta, SEXP Ipar, SEXP Ltol, SEXP Tol, SEXP Fixpnt, SEXP Eps,
		SEXP Epsmin, SEXP func, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP Guess, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, RWORK;

  int  j, k, ii, nx, ncomp;
  double *aleft, *aright, *zeta, *guess, *fspace, *tol, *fixpnt, *z;
  double eps, epsmin, xout;
  int *m, *ispace, *ipar, ltol, iflag;

  deriv_func2    *derivs;
  jac_func2      *jac;
  jacbound_func2 *jacbound;
  bound_func2    *bound;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();

  ncomp = INTEGER(Ncomp)[0];   /* number of equations -global variable */
  n_eq  = INTEGER(Ncomp)[0];   /* number of equations -global variable */
  mstar = INTEGER(Mstar)[0];   /* number of variables */

  m  = (int *) R_alloc(n_eq, sizeof(int));  /* order of diff eqns */
  for (j = 0; j < n_eq; j++) m[j] = INTEGER(M)[j];

  ii = LENGTH(Aleft);
  aleft  =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) aleft[j] = REAL(Aleft)[j];

  ii = LENGTH(Aright);
  aright =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) aright[j] = REAL(Aright)[j];

  ii = LENGTH(Zeta);
  zeta   =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) zeta[j] = REAL(Zeta)[j];

  ii = LENGTH(Ipar);    /* length of ipar */
  ipar  = (int *)    R_alloc(ii, sizeof(int));
  for (j = 0; j < ii;j++) ipar[j] = INTEGER(Ipar)[j];

  ltol = INTEGER(Ltol)[0];
  tol   =(double *) R_alloc(ltol, sizeof(double));
    for (j = 0; j < ltol;j++) tol[j] = REAL(Tol)[j];

  ii =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) fixpnt[j] = REAL(Fixpnt)[j];

  ii = LENGTH(Guess);
  guess = (double*) R_alloc(ii, sizeof(double));
    for (j = 0; j<ii;j++) guess[j] = REAL(Guess)[j];

  ii = ipar[5];
  ispace = (int *) R_alloc(ii, sizeof(int));

  ii = ipar[4];
  fspace = (double *) R_alloc(ii, sizeof(double));

  eps    = REAL(Eps)[0];
  epsmin = REAL(Epsmin)[0];


/* initialise global R-variables... */

  PROTECT(X  = NEW_NUMERIC(1));               incr_N_Protect();
  PROTECT(EPS = NEW_NUMERIC(1));              incr_N_Protect();
  PROTECT(J = NEW_INTEGER(1));                incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,mstar));    incr_N_Protect();
/*  PROTECT(bvp_gparms = parms);                incr_N_Protect();*/

      bvpcolmod_envir = rho;

      /* interface function between fortran and R passed to Fortran*/
      derivs = (deriv_func2 *) colsys_derivs;
      colsys_deriv_func = func;

      jac = colsys_jac;
      colsys_jac_func = jacfunc;

      bound = colsys_bound;
      colsys_bound_func = boundfunc;

      jacbound = colsys_jacbound;
      colsys_jacbound_func = jacboundfunc;

/* Call the fortran function -
      Subroutine Colmod(Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)             */

	  F77_CALL(colmod) (&ncomp, m, aleft, aright, zeta, ipar, &ltol,
        tol, fixpnt, ispace, fspace, &iflag, &eps, &epsmin,
        derivs,jac,bound,jacbound,colsys_guess);

/*             Call Appsln(Xx,Z,Fspace,Ispace)
C....   Iflag - The Mode Of Return From Colmod.
C....         =  1  For Normal Return
C....         =  0  If The Collocation Matrix Is Singular For The Final
C....               Continuation Problem.
C....         = -1  If The Expected No. Of Subintervals Exceeds Storage
C....               Specifications.
C....         = -2  If The Nonlinear Iteration Has Not Converged For The
C....               Final Continuation Problem.
C....         = -3  If There Is An Input Data Error.

	  if (iflag == -1)
     {
      warning("an excessive amount of work (> maxsteps ) was done, but integration was successful - increase maxsteps");
     }  else             */

  if (iflag == 0)
	{
	  unprotect_all();
	  error("The collocation matrix is singular for the final continuation problem\n");
	}
  else if (iflag == -1)
	{
	  unprotect_all();
	  error("The Expected No. Of Subintervals Exceeds Storage Specifications.\n");
	}
  else if (iflag == -2)
	{
	  unprotect_all();
	  error("The Nonlinear Iteration Has Not Converged For The Final Continuation Problem.\n");
	}
  else  if (iflag == -3)
	{
	  unprotect_all();
	  error("Illegal input to colmod\n");
	}
  else
	{
    nx = LENGTH(Xout);
    z  =(double *) R_alloc(mstar, sizeof(double));

    PROTECT(yout = allocMatrix(REALSXP,mstar+1,nx));incr_N_Protect();
	  for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(appsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */
  ii = ncomp+7;
  PROTECT(ISTATE = allocVector(INTSXP, ii+1));incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  for (k = 0;k<ii;k++) INTEGER(ISTATE)[k+1] = ispace[k];
  ii = ispace[6];
  PROTECT(RWORK = allocVector(REALSXP, ii));incr_N_Protect();
  for (k = 0;k<ii;k++) REAL(RWORK)[k] = fspace[k];
  setAttrib(yout, install("istate"), ISTATE);
  setAttrib(yout, install("rstate"), RWORK);
  }
/*                       ####   termination   ####                            */    
  unprotect_all();
  return(yout);
}
