\name{diagnose}
\alias{diagnose}
\title{Print Diagnostic Characteristics of ODE and DAE Solvers}
\description{
  Prints several diagnostics of the simulation to the screen, e.g.
  number of steps taken, the last step size,...
}
\usage{
diagnose(obj)
}
\arguments{
  \item{obj }{is the output matrix as produced by twpbvp.
  }
}
\value{
  The integer and real vector with diagnostic values; for function \code{lsodar}
  also the root information.
  
  See tables 2 and 3 in vignette("deSolve") for what these vectors contain.
}
\details{
  When the integration output is saved as a data.frame, then the required
  attributes are lost and method \code{diagnostics} will not work anymore.
}
\examples{
f2 <- function(x,y,parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2]-(1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy,dy2))
}

x    <- seq(1,6,0.1)
sol  <- bvptwp(yini=c(y=1,dy=NA),yend=c(-0.5,NA),x=x,func=f2,guess=1)
plot(sol,which="y")
diagnose(sol)
}
\keyword{ utilities }