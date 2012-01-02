\name{rkMethodPlot}
\alias{rkMethodPlot}
\title{Plots the steps in runge-kutta methods}

\description{
...
}
\usage{
rkMethodPlot (rk, ...) 
}

\arguments{
  \item{rk }{A list containing the runge-kutta coefficients,
    implicit or explicit, e.g. matrix \code{A}, vectors \code{b1}, \code{b2},
   code{c}. The list can be of type \code{rkMethod}, as defined in package
   \code{deSolve}.
  }
  \item{... }{arguments passed to the plotting function.
  }
}
\value{
  Returns nothing
}
\author{Karline Soetaert 
}
\examples{
# This to plot all runge kutta methods
#RKS <- rkMethod()
#for (i in 4:21)   rkMethodPlot( rkMethod(RKS[i]))

## -----------------------------------------------------------------------------
## Figures A and B: Cash-Karp and Radau 5 steps
## -----------------------------------------------------------------------------

par(mfrow=c(2,2))

rkMethodPlot( rkMethod("rk45ck"), main="Cash-Karp")
writelabel("A")

rkMethodPlot( rkMethod("irk5"), main="Radau5")
writelabel("B")

rkMethodPlot( rkMethod("rk45dp6"), main="Dopri")
writelabel("C")
                                                                                   
rkMethodPlot( rkMethod("irk6l"), main="Lobatto")
writelabel("D")

legend("bottomright", pch = c(16, 16, 1, NA), pt.cex = c(1.5, 1.5, 1), 
       legend = c(expression(y[0]), expression(y[1]), "intermediary", "k"),
       col = c("grey", "black", "black", "black"), lty = c(NA, NA, NA, 1), 
       lwd = c(1, 1, 1, 2))

}
\seealso{
\link{stability.bruteforce} for plotting stability regions of Runge-Kutta methods.
}
\keyword{plot}

