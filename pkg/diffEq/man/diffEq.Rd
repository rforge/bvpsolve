\name{diffEq-package}
\alias{diffEq-package}
\alias{diffEq}
\docType{package}
\title{
  Functions to create the figures of the book 
  "solving differential equations in R" by Karline Soetaert, Jeff R. Cash and 
  Francesca Mazzia.
  
  Published by Springer
}
\description{
  R package diffEq contains the functions for generating figures relating to 
  differential equations
}
\details{
  \tabular{ll}{
    Package: \tab diffEq\cr
    Type:    \tab Package\cr
    Version: \tab 1.0\cr
    License: \tab GNU Public License 2 or above\cr
  }

 
}
\author{
  Karline Soetaert (Maintainer),
  
  Jeff Cash,
  
  Francesca Mazzia
}
\seealso{
  \code{\link{rkMethodPlot}} for plotting the steps of runge-kutta methods
  
  \link{stability.multistep}, \link{stability.bruteforce} for plotting stability regions

  \link{Coefficients} for the BDF, AdamsMoulton, AdamsBashford coefficients.
  
  \code{rkMethod} from package \code{deSolve} for the Runge-Kutta coefficients
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(rkMethodPlot)
example(stability.multistep)
example(BDF)

## open the directory with R sourcecode examples
browseURL(paste(system.file(package = "diffEq"), "/doc/examples", sep = ""))

## show package vignette with how to use the package
## + source code of the vignette
vignette("diffEq")
edit(vignette("diffEq"))

}
}

\keyword{ package }