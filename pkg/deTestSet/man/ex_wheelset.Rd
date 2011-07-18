\name{wheelset}
\alias{wheelset}
\title{Wheel Set problem, mechanics, Index 2 IDE}
\description{
Describes the motion of a simple wheelset on a rail track. 

It is a differential algebraic equation of index 2, 17 equations.

}
\usage{
wheelset (times = seq(0, 10, by = 0.001), yini = NULL, dyini = NULL, 
  parms = list(), method = "mebdfi", maxsteps = 1e5, ...) 
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{dyini }{the initial derivatives of the state variables of the DE
    system.  
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{parms }{list of parameters that overrule the default parameter values
  }
  \item{method }{the solver to use; only \code{mebdfi} available for now
  }
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
  }
  \item{... }{additional arguments passed to the solver .
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{yini}, plus an additional column (the first) 
  for the time value.

  There will be one row for each element in \code{times} unless the
  solver returns with an unrecoverable error.  If
  \code{yini} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>

Francesca Mazzia
}
\examples{
out <- wheelset()
plot(out, lwd=2)
max(abs(out[nrow(out), -1] - reference("wheelset")))
}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\note{
This model is implemented in R.
}

\details{

  The default parameters are:
       MR = 16.08, G = 9.81, V = 30., RN0 = 0.1, LI1 = 0.0605,
       LI2 = 0.366, MA = 0.0, HA = 0.2, MU = 0.12 , XL = 0.19, 
       CX = 6400., CZ = 6400. ,
       E = 1.3537956, GG = 0.7115218, SIGMA = 0.28, GM = 7.92e10,
       C11 = 4.72772197, C22 = 4.27526987, C23 = 1.97203505,
       DELTA0 = 0.0262, AR = 0.1506, RS = 0.06, EPS = 0.00001, 
       B1 = 0.0, B2 = 4.0  
}
\keyword{utilities}
