\name{bvpcolmod}
\alias{bvpcolmod}
\title{
  Solver boundary value problems of ordinary differential equations -
  Automatic Continuation With Collocation method
}
\description{
  Solves `Stiff' Boundary Value Problems For Ordinary
  Differential Equations By Using Continuation.
  This is an implementation of the fortran code colmod written by
  J.R. Cash and M.H. Wright.
}

\note{
  colmod is a revised version of the package colnew by
  Bader and Ascher (1987), which in turn is a modification of the package
  colsys by Ascher, Christiansen and Russell (1981). Colmod has been
  adapted to allow an automatic continuation strategy to
  be Used (Cash et al., 1995).
  The mesh selection algorithm used in colmod differs from that used in colnew
}

\usage{
bvpcolmod(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
  ynames = NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
  posbound =NULL, leftbc=NULL, islin=FALSE, colp= NULL,
  epsini=NULL, eps=epsini, atol=1e-8, nmax=1000, verbose=FALSE, ...)
}
\arguments{
  \item{yini }{either a vector with the initial (state) values for the ODE
    system, or \code{NULL}. If \code{yini} is a vector, use \code{NA} for an
    initial value which is not available.
  }
  \item{x }{ sequence of the independent variable for which output is wanted;
    the first value of \code{x} must be the initial value (at which
    \code{yini} is defined), the final value the end condition (at which
    \code{yend} is defined).
  }
  \item{func }{an \R-function that computes the values of the derivatives in
    the ODE system (the model definition) at x.

    If \code{eps} is NULL then \code{func} must be defined as:
    \code{yprime = func(x, y, parms,...)}.  \code{x} is the current point
    of the independent variable in the integration, \code{y} is the current
    estimate of the (state) variables in the ODE system.  If the initial
    values vector \code{yini} has a names attribute, the names will be
    available inside \code{func}.  \code{parms} is a vector or list of
    parameters; ... (optional) are any other arguments passed to the function.

    If \code{eps} is given a value then \code{func} must be defined as:
    \code{yprime = func(x, y, parms, eps,...)}.
      
    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{x}, and whose next elements are global values that are required at
    each point in \code{x}.
  }
  \item{yend }{either a vector with the final (state) values for the ODE
    system, or \code{NULL}; if \code{yend} is a vector use \code{NA}
    for an end value which is not available.
  }
  \item{parms }{parameters passed to \code{func}.
  }
  \item{guess }{guess for the value(s) of the unknown initial conditions;
    i.e. one value for each \code{NA} in \code{yini}.
    The length of \code{guess} should thus equal the number of NAs
    in \code{yini}. If not provided, a value = 0 is assumed for each
    \code{NA} and a warning printed.
  }
  \item{ynames }{only if \code{yini} and \code{yend} not provided or their
    elements not given names: the names of the y's.
  }
  \item{jacfunc }{jacobian (optional) - an \R-function that evaluates the
    jacobian of \code{func} at point \code{x},
    \code{jacfunc} must be defined as: \code{jac = func(x, y, parms, eps,...)}
    (where \code{eps} should not be present, if it is not defined as
    an input argument to \code{bvpcolmod}. It should return the partial
    derivatives of func with respect to y,
    i.e. df(i,j) = dfi/dyj. See last example.

    If \code{jacfunc} is \code{NULL}, then a numerical approximation using
    differences is used.
  }
  \item{bound }{boundary function (optional) - only if \code{yini} and
    \code{yend} are not available. An \R function that evaluates
    the i-th boundary element at point \code{x}. It should be defined as:
    \code{bound = func(i, y, parms, eps,...)} (where \code{eps} should
    not be present, if it is not defined as an input argument to
    \code{bvpcolmod}. It should return the i-th boundary condition.
    See last example.
   }
  \item{jacbound }{jacobian of the boundary function (optional) - only if
    \code{bound} is defined. An \R function that evaluates the gradient of
    the i-th boundary element with respect to the state variables, at
    point \code{x}.
    It should be defined as: \code{jacbound = func(i, x, y, parms, eps,...)}
    (where \code{eps} should not be present, if it is not defined as
    an input argument to \code{bvpcolmod}. It should return the gradient
    of the i-th boundary condition. See last example.

    If \code{jacbound} is \code{NULL}, then a numerical approximation using
    differences is used.
   }
  \item{posbound }{only if \code{bound} is given: the position (in the
    mesh) of the boundary conditions -corresponds to fortran input
    "Zeta" in colmod FORTRAN code.
  }
  \item{leftbc}{only if \code{yini} and \code{yend} are not available and
    \code{posbound} is not specified: the number of left boundary conditions.
  }
  \item{islin }{set to \code{TRUE} if the problem is linear - this will
    speed up the simulation.
  }
  \item{colp }{number of collocation points per subinterval.
  }
  \item{epsini }{the initial value of the continuation parameter. If
    \code{NULL} and \code{eps} is given a value, then \code{epsini} takes
    the default starting value of 0.5. For many singular perturbation type
    problems, the choice of 0.1 < \code{eps} < 1 represents a (fairly) easy
    problem. The user should attempt to specify an initial
    problem that is not `too' challenging.
    \code{epsini} must be initialised strictly less than 1 and greater than 0.
  }
  \item{eps }{the desired value of precision for which the user would like
    to solve the problem. \code{eps} must be less than or equal to
    \code{epsini}.
  }
  \item{nmax }{maximal number of subintervals.
  }
  \item{atol }{ error tolerance, a scalar.
  }
  \item{verbose }{when \code{TRUE}, writes messages (warnings) to the screen.
  }
  \item{... }{additional arguments passed to the model functions.
  }
}
\value{
  A matrix with up to as many rows as elements in times and as many columns
  as elements in \code{yini} plus the number of "global" values returned
  in the second element of the return from \code{func}, plus an additional
  column (the first) for the x-value.

  There will be one row for each element in \code{x} unless the solver returns
  with an unrecoverable error.

  If \code{yini} has a names attribute, it will be used to label the
  columns of the output value.

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\details{
  This is an implementation of the Automatic Continuation With Collocation
  method called "colmod" to solve boundary value problems of ordinary
  differential equations.
  
  The following description is adapted from the documentation of the colmod
  source code.

  A boundary value problem does not have all initial values of
  the state variable specified. Rather some conditions are specified at
  the end of the integration interval.

  The ODEs and boundary conditions are made available through the
  user-provided routines, \code{func} and vectors \code{yini} and \code{yend}
  or (optionally) \code{bound}.

  The corresponding partial derivatives are optionally available through the
  user-provided routines, \code{jacfunc} and \code{jacbound}. Default is that
  they are automatically generated by \R, using numerical differences.

  The user-requested tolerance is provided through \code{atol}.

  If the function terminates because the maximum
  number of subintervals was exceeded, then it is recommended that
  'the program be run again with a larger value for this maximum.'

  The type of problems which \code{bvpcolmod} is designed to solve typically
  involve a small positive parameter 0 < eps << 1.
  As eps becomes progressively smaller, the problem normally becomes
  increasingly difficult to approximate numerically (for example, due
  to the appearance of narrow transition layers in the profile of the
  analytic solution).

  The idea of continuation is to solve a chain of problems in which the
  parameter eps decreases monotonically towards some desired value.
  That is, a sequence of problems is attempted to be solved:

  epsini  > eps1  > eps2  > eps3  >  .....  > eps  > 0

  where \code{epsini} is a user provided starting value and \code{eps} is a
  user desired final value for the parameter.
}

\seealso{
%%  \code{\link{bvpshoot}} for the shooting method

%%  \code{bvptwp} for a MIRK formula
}

\references{
  U. Ascher, J. Christiansen and R. D. Russell,  (1981)
  collocation software for boundary-value odes,
  acm trans. math software 7, 209-222.

  G. Bader and U. Ascher, (1987)
  a new basis implementation for a mixed order
  boundary value ode solver,
  siam j. scient. stat. comput. 8, 487-483.

  J. R. Cash, G. Moore and R. W. Wright, (1995)
  an automatic continuation strategy for the solution of
  singularly perturbed linear two-point boundary value problems,
  j. comp. phys. 122, 266-279.
}

\keyword{math}

\examples{
#################################################################
# Example 1: simple standard problem
# solve the BVP ODE:
# d2y/dt^2=-3py/(p+t^2)^2
# y(t= -0.1)=-0.1/sqrt(p+0.01)
# y(t=  0.1)= 0.1/sqrt(p+0.01)
# where p = 1e-5
#
# analytical solution y(t) = t/sqrt(p + t^2).
#
# The problem is rewritten as a system of 2 ODEs:
# dy=y2
# dy2=-3p*y/(p+t^2)^2
################################################################################

#--------------------------------
# Derivative function
#--------------------------------
fun <- function(t,y,pars)
{ dy1 <-y[2]
  dy2 <- - 3*p*y[1]/(p+t*t)^2
  return(list(c(dy1,
             dy2))) }


# parameter value
p    <-1e-5

# initial and final condition; second conditions unknown
init <- c(-0.1/sqrt(p+0.01), NA)
end  <- c(0.1/sqrt(p+0.01), NA)

# Solve bvp
sol  <- as.data.frame(bvpcolmod(yini=init,x=seq(-0.1,0.1,by=0.001),
        func=fun, yend=end, guess=1))
plot(sol$time,sol[,2],type="l")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="p")

#################################################################
# Example 1b: simple
# solve d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y = sqrt(x)*cos(x),
# on the interval [1,6] and with boundary conditions:
# y(1)=1, y(6)=-0.5
#
# Write as set of 2 odes
# dy/dx = y2
# dy2/dx  = - 1/x*dy/dx - (1-1/(4x^2)y + sqrt(x)*cos(x)
#################################################################

f2 <- function(x,y,parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2]-(1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy,dy2))
}

x    <- seq(1,6,0.1)
sol  <- bvpcolmod(yini=c(1,NA),yend=c(-0.5,NA),x=x,func=f2,guess=1)
plot(sol)

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x),add=TRUE,type="l")

################################################################################
# Example 3. Uses continuation
# Test problem 24
################################################################################

Prob24<- function(t,y,pars,ks)      #eps is called ks here
{
A <- 1+t*t
AA <- 2*t
ga <- 1.4
 list(c(y[2],(((1+ga)/2 -ks*AA)*y[1]*y[2]-y[2]/y[1]-
              (AA/A)*(1-(ga-1)*y[1]^2/2))/(ks*A*y[1])))

}

ks  <- 0.001
ini <- c(0.9129,NA)
end <- c(0.375,NA)

# bvpcolmod works with eps NOT too small...

solcol  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess=c(0),eps=ks,verbose=FALSE)
plot(solcol[,1],solcol[,2],type="l",main="test problem 24",
     lwd=2,col="blue")

solcol2  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess=c(0),eps=0.01)
lines(solcol2,lty=2,lwd=2)

solcol3  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess=c(0),eps=0.1)
lines(solcol3,lty=3,col="green",lwd=2)

legend("topright",col=c("blue","black","green2"),lty=1:3,
       lwd=2,legend=c("0.001","0.01","0.1"),title="eps")

#################################################################
# Example 4  - solved with specification of boundary, and jacobians
# d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
# y(0)=y'(0)=0
# y(1)=1, y'(1)=0
#
# dy/dx  = y2
# dy2/dx = y3    (=d2y/dx2)
# dy3/dx = y4    (=d3y/dx3)
# dy4/dx = R*(y2*y3 -y*y4)
#################################################################

f2<- function(x,y,parms,R)
{
  list(c(y[2],y[3],y[4],R*(y[2]*y[3]-y[1]*y[4]) ))
}

df2 <- function(x,y,parms,R)
{
 df <- matrix(nr=4,nc=4,byrow=TRUE,data=c(
             0,1,0,0,
             0,0,1,0,
             0,0,0,1,
             -1*R*y[4],R*y[3],R*y[2],-R*y[1]))
}

g2 <- function(i,y,parms,R)
{
  if ( i ==1) return(y[1])
  if ( i ==2) return(y[2])
  if ( i ==3) return(y[1]-1)
  if ( i ==4) return(y[2])
}

dg2 <- function(i,y,parms,R)
{
  if ( i ==1) return(c(1,0,0,0))
  if ( i ==2) return(c(0,1,0,0))
  if ( i ==3) return(c(1,0,0,0))
  if ( i ==4) return(c(0,1,0,0))
}

init <- c(1,NA)
R    <- 0.01
sol  <- as.data.frame(bvpcolmod(x=seq(0,1,by=0.01),
        posbound=c(0,0,1,1),func=f2, guess=1,R=R,
        bound=g2, jacfunc=df2, jacbound=dg2))
plot(sol[,1:2])

}

