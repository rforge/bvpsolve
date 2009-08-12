\name{mebdfi}
\alias{mebdfi}
\title{Solver for Differential Algebraic Equations (DAE) up to index 3}
\description{
  Solves either:
  \itemize{
    \item a system of ordinary differential equations (ODE) of the form
      \deqn{y' = f(t,y,...)} or
    \item a system of differential algebraic equations (DAE) of the form
      \deqn{F(t,y,y') = 0}
  }
  using a combination of backward differentiation formula (BDF) and a
  direct linear system solution method (dense or banded).
  
  The \R function \code{mebdfi} provides an interface to the Fortran DAE
  solver of the same name, written by T.J. Abdulla and J.R. Cash.
    
  The system of DE's is written as an \R function (which may, of course,
  use \code{\link{.C}}, \code{.Fortran}, \code{\link{.Call}}, etc., to
  call foreign code) or be defined in compiled code that has been
  dynamically loaded.
}
\usage{
mebdfi(y, times, func = NULL, parms, dy = NULL, res = NULL,
  nind=c(length(y),0,0), rtol = 1e-6, atol = 1e-8, jacfunc = NULL,
  jacres = NULL, jactype = "fullint", verbose = FALSE,
  tcrit = NULL,  hini = 0, ynames = TRUE, maxord = 7, bandup = NULL,
  banddown = NULL, maxsteps = 5000, dllname = NULL,
  initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL,
  forcings=NULL, initforc = NULL, fcontrol=NULL,  ...)
}

\arguments{
  \item{y }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times} = \code{NULL}.
  }
  \item{func }{cannot be used if the model is a DAE system. If an ODE
    system, \code{func} should be an \R-function that computes the
    values of the derivatives in the ODE system (the \emph{model
    definition}) at time t.
  
    \code{func} must be defined as: \code{yprime = func(t, y,
    parms,...)}.  \cr \code{t} is the current time point in the
    integration, \code{y} is the current estimate of the variables in
    the ODE system.  If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}, unless
    \code{ynames} is FALSE.  \code{parms} is a vector or list of
    parameters.  \code{...} (optional) are any other arguments passed to
    the function.
    
    The return value of \code{func} should be a list,
    whose first element is a vector containing the derivatives of
    \code{y} with respect to \code{time}, and whose next elements are
    global values that are required at each point in \code{times}.
    The derivatives should be specified in the same order as the specification
    of the state variables \code{y}.

    Note that it is not possible to define \code{func} as a compiled
    function in a dynamically loaded shared library. Use \code{res}
    instead.
  }
  \item{parms }{vector or list of parameters used in \code{func},
    \code{jacfunc}, or \code{res}
  }
  \item{dy }{the initial derivatives of the state variables of the DE
    system.  Ignored if an ODE.
  }
  \item{res }{if a DAE system: either an \R-function that computes the
    residual function F(t,y,y') of the DAE system (the model
    defininition) at time \code{t}, or a character string giving the
    name of a compiled function in a dynamically loaded shared library.
      
    If \code{res} is a user-supplied R-function, it must be called as:
    \code{F = func(t, y, dy, parms, ...)}.
    
    Here \code{t} is the current time point in the integration, \code{y}
    is the current estimate of the variables in the ODE system,
    \code{dy} are the corresponding rates of change.  If the initial
    \code{y} or \code{dy} have a names attribute, the names will be
    available inside \code{func}, unless \code{ynames} is \code{FALSE}.
    \code{parms} is a vector of parameters.
    
    The return value of \code{res} should be a list, whose first element
    is a vector containing the residuals of the DAE system,
    i.e. delta = F(t,y,y'), and whose next elements contain output
    variables that are required at each point in \code{times}.
  
    If \code{res} is a string, then \code{dllname} must give the name of
    the shared library (without extension) which must be loaded before
    \code{mebdfi()} is called (see package vignette \code{"compiledCode"}
    for more information).
  }
  \item{nind }{if a DAE system: a three-valued vector with the number of
    variables of index 1, 2, 3 respectively.
    The equations must be defined such that the index 1 variables precede
    the index 2 variables which in turn precede the index 3 variables.
    The sum of the variables of different index should equal N,
    the total number of variables.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations. Only used in case
    the system is an ODE (y' = f(t,y)), specified by \code{func}. The \R
    calling sequence for \code{jacfunc} is identical to that of
    \code{func}.
    
    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix dydot/dy, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and Fortran store matrices).
    
    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of lsode.
  }
  \item{jacres }{ \code{jacres} and not \code{jacfunc} should be used if
    the system is specified by the residual function F(t,y,y'),
    i.e. \code{jacres} is used in conjunction with \code{res}.
    
    If \code{jacres} is an \R-function, the calling sequence for
    \code{jacres} is identical to that of \code{res}, but with extra
    parameter \code{cj}.  Thus it should be called as: \code{jacres = 
    func(t, y, dy, parms, cj, ...)}.  Here \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system, \eqn{y'} are the corresponding rates of
    change and \code{cj} is a scalar, which is normally proportional to
    the inverse of the stepsize. If the initial y or dy have a names
    attribute, the names will be available inside jacres, unless
    \code{ynames} is \code{FALSE}. \code{parms} is a vector of
    parameters (which may have a names attribute).
    
    If the Jacobian is a full matrix, \code{jacres} should return the
    matrix dG/dy + cj*dG/dyprime, where the ith row is the sum of the
    derivatives of \eqn{G_i} with respect to \eqn{y_j} and the scaled
    derivatives of \eqn{G_i} with respect to \eqn{dy_j}.
    
    Ff the Jacobian is banded, \code{jacres} should return only the
    nonzero bands of the Jacobian, rotated rowwise. See details for the
    calling sequence when \code{jacres} is a string.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by the user.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will output
    the settings of vector \bold{istate} and \bold{rstate} - see details.
  }
  \item{tcrit }{the Fortran routine \code{mebdfi} overshoots its targets
    (times points in the vector \code{times}), and interpolates values
    for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined by the solver
  }
  \item{ynames }{if FALSE: names of state variables are not passed to
    function \code{func} ; this may speed up the simulation especially
    for large models.
  }
  \item{maxord }{the maximum order to be allowed. Reduce \code{maxord}
    to save storage space ( <= 5)
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver; will be recalculated to be at least 500 and a multiple of
    500; the solver will give a warning if more than 500 steps are
    taken, but it will continue till \code{maxsteps} steps.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{func}/or \code{res} and
    \code{jacfunc}.  See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (fortran) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \file{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    See "Details" of \link{vode} or vignette \code{compiledCode}.
  }
  \item{... }{additional arguments passed to \code{func},
    \code{jacfunc}, \code{res} and \code{jacres}, allowing this to be a
    generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the next elements of the return from \code{func} or
  \code{res}, plus an additional column (the first) for the time value.
  There will be one row for each element in \code{times} unless the
  Fortran routine `mebdfi' returns with an unrecoverable error.  If
  \code{y} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## Coupled chemical reactions including an equilibrium
## modeled as (1) an ODE and (2) as a DAE
##
## The model describes three chemical species A,B,D:
## subjected to equilibrium reaction D <- > A + B
## D is produced at a constant rate, prod
## B is consumed at 1s-t order rate, r
## Chemical problem formulation 1: ODE
## =======================================================================

## Dissociation constant
K <- 1 

## parameters
pars <- c(
        ka   = 1e6,     # forward rate
        r    = 1,
        prod = 0.1)


Fun_ODE <- function (t, y, pars)
{
  with (as.list(c(y, pars)), {
    ra  <- ka*D        # forward rate
    rb  <- ka/K *A*B   # backward rate

    ## rates of changes
    dD  <- -ra + rb + prod
    dA  <-  ra - rb
    dB  <-  ra - rb - r*B
    return(list(dy = c(dA, dB, dD),
                CONC = A+B+D))
  })
}

## =======================================================================
## Chemical problem formulation 2: DAE
## 1. get rid of the fast reactions ra and rb by taking
## linear combinations   : dD+dA = prod (res1) and
##                         dB-dA = -r*B (res2)
## 2. In addition, the equilibrium condition (eq) reads:
## as ra = rb : ka*D = ka/K*A*B = >      K*D = A*B
## =======================================================================

Res_DAE <- function (t, y, yprime, pars)
{
  with (as.list(c(y, yprime, pars)), {
    ## residuals of lumped rates of changes
    res1 <- -dD - dA + prod
    res2 <- -dB + dA - r*B
    
    ## and the equilibrium equation
    eq   <- K*D - A*B

    return(list(c(res1, res2, eq),
                CONC = A+B+D))
  })
}

times <- seq(0, 100, by = 1)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = 2*3/K)
require(mebdfi)

## ODE model solved with mebdfi
ODE <- as.data.frame(mebdfi(y = y, times = times, func = Fun_ODE,
                     parms = pars, atol = 1e-10, rtol = 1e-10))

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0) 
## DAE model solved with mebdfi
DAE <- as.data.frame(mebdfi(y = y, dy = dy, times = times, hini=1e-10,
         res = Res_DAE, parms = pars, atol = 1e-10, rtol = 1e-10))

## ================
## plotting output
## ================
opa <- par(mfrow = c(2,2))

for (i in 2:5)
{
plot(ODE$time,ODE[,i],xlab = "time",
     ylab = "conc",main = names(ODE)[i],type = "l")
points(DAE$time,DAE[,i],col = "red")
}
legend("bottomright",lty = c(1,NA),pch = c(NA,1),
       col = c("black","red"),legend = c("ODE","DAE"))      
       
# difference between both implementations:
max(abs(ODE-DAE))

par(mfrow = opa)

\dontrun{
## =======================================================================
## same DAE model, now with the Jacobian
## =======================================================================
jacres_DAE <- function (t, y, yprime, pars, cj)
{
    with (as.list(c(y, yprime, pars)), {
##    res1 = -dD - dA + prod
      PD[1,1] <- -1*cj      # d(res1)/d(A)-cj*d(res1)/d(dA)
      PD[1,2] <- 0          # d(res1)/d(B)-cj*d(res1)/d(dB)
      PD[1,3] <- -1*cj      # d(res1)/d(D)-cj*d(res1)/d(dD)
##     res2 = -dB + dA - r*B
      PD[2,1] <- 1*cj
      PD[2,2] <- -r -1*cj
      PD[2,3] <- 0
##    eq = K*D - A*B
      PD[3,1] <- -B
      PD[3,2] <- -A
      PD[3,3] <- K
      return(PD)
   })
}

PD <- matrix(nc = 3, nr = 3, 0)

DAE2 <- as.data.frame(mebdfi(y = y, dy = dy, times = times,
          res = Res_DAE, jacres = jacres_DAE, jactype = "fullusr",
          parms = pars, atol = 1e-10, rtol = 1e-10))
         
max(abs(DAE-DAE2))

## See \dynload subdirectory for a FORTRAN implementation of this model

## =======================================================================
## The chemical model as a DLL, with production a forcing function
## =======================================================================
times <- seq(0, 100, by = 2)

pars <- c(K = 1, ka   = 1e6, r    = 1)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = 2*3/pars["K"])

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0)

# production increases with time
prod <- matrix(nc=2,data=c(seq(0,100,by=10),0.1*(1+runif(11)*1)))

ODE_dll <- as.data.frame(mebdfi(y=y,dy=dy,times=times,res="chemres",
          dllname="deSolve", initfunc="initparms",
          initforc="initforcs", parms=pars, forcings=prod,
          atol=1e-10,rtol=1e-10,nout=2, outnames=c("CONC","Prod")))

opa <- par(mfrow = c(1,2))
plot(ODE_dll$time,ODE_dll$Prod,xlab = "time",
     ylab = "/day",main = "production rate",type = "l")
plot(ODE_dll$time,ODE_dll$D,xlab = "time",
     ylab = "conc",main = "D",type = "l")
par(mfrow = opa)
}
}
\references{
J. R. Cash, The integration of stiff initial value problems
         in O.D.E.S using modified extended backward differentiation
         formulae, Comp. and Maths. with applics., 9, 645-657, (1983).
J.R. Cash and S. Considine, an MEBDF code for stiff
         initial value problmes, ACM Trans Math Software, 142-158,
         (1992).
J.R. Cash, Stable recursions with applications to the
         numerical solution of stiff systems, Academic Press,(1979).
}
\details{
  The mebdfi solver uses modified extended backward differentiation
  formulas of orders one through eight (specified with \code{maxord})
  to solve either:
  \itemize{
    \item an ODE system of the form \deqn{y' = f(t,y,...)}  for y = Y,
    or
    \item a DAE system of the form \deqn{F(t,y,y') = 0} for y = Y and y'
    = YPRIME.
  }
  

  The recommended value of \code{maxord} is eight, unless it is believed
  that there are severe stability problems in which case \code{maxord} =
  4 or 5 should be tried instead.

  ODEs are specified in \code{func}, DAEs are specified in \code{res}.
  
  If a DAE system, Values for Y \emph{and} YPRIME at the initial time
  must be given as input. Ideally,these values should be consistent,
  that is, if T, Y, YPRIME are the given initial values, they should
  satisfy F(T,Y,YPRIME) = 0.  \cr However, if consistent values are not
  known, in many cases mebdfi can solve for them: when \code{estini} = 1,
  y' and algebraic variables (their number specified with \code{nalg})
  will be estimated, when \code{estini} = 2, y will be estimated.
  
  The form of the \bold{Jacobian} can be specified by
  \code{jactype}. This is one of:
  
  \describe{
    \item{jactype = "fullint":}{a full Jacobian, calculated internally
      by \code{mebdfi}, the default,
    }
    \item{jactype = "fullusr":}{a full Jacobian, specified by user
      function \code{jacfunc} or \code{jacres},
    }
    \item{jactype = "bandusr":}{a banded Jacobian, specified by user
      function \code{jacfunc} or \code{jacres}; the size of the bands
      specified by \code{bandup} and \code{banddown},
    }
    \item{jactype = "bandint":}{a banded Jacobian, calculated by
      \code{mebdfi}; the size of the bands specified by \code{bandup} and
      \code{banddown}.
    }
  }
  
  If \code{jactype} = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc}.

  If jactype = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc} or \code{jacres}.
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver.  If the request for
  precision exceeds the capabilities of the machine, mebdfi will return
  an error code. See \code{\link{lsoda}} for details.
  
  \bold{res and jacres} may be defined in compiled C or Fortran code, as
  well as in an R-function. See package vignette \code{"compiledCode"}
  for details.  Examples
  in Fortran are in the \file{dynload} subdirectory of the
  \code{deSolve} package directory.

  The output will have the attributes \bold{istate} and \bold{rstate},
  two vectors with several
  useful elements. Their values can be displayed using R-function
  \code{attributes}, or by calling \code{\link{diagnostics}}

  If \code{verbose} = TRUE, the settings of istate and rstate will be
  written to the screen.
  See vignette("deSolve") for an explanation of each element in these vectors.

}
\seealso{
  \itemize{
    \item \code{\link{daspk}} another DAE solver,
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\note{

}
\keyword{math}

