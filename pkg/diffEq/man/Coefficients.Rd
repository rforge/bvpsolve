\name{Coefficients}
\alias{Coefficients}
\alias{BDF}
\alias{AdamsMoulton}
\alias{AdamsBashforth}
\docType{data}

\title{The coefficients of multistep methods}

\description{
Coefficients \code{alpha} and \code{beta} of the Adams-Bashforth, 
Adams-Moulton and Backward differentiation formulae.

sum_{j=0}^k alpha_j y_{n-j} = h sum_(j=0)^k beta_j f(x_{n-j},y_{n-j})

For the BDF methods, the \code{angle} of the stability-region (the \code{alpha}
of the A(alpha) stability, in radians is also given.
}
\usage{
BDF
AdamsMoulton
AdamsBashforth
}

\author{Karline Soetaert 
}
\examples{
## =============================================================================
## Stability properties
## =============================================================================

BDF

stability.multistep(alpha = BDF$alpha[3,], beta = BDF$beta[3,], 
  xlim = c(-7,7), ylim = c(-7,7))

stability.multistep(alpha = AdamsMoulton$alpha[3,], beta = AdamsMoulton$beta[3,], 
  xlim = c(-7,7), ylim = c(-7,7) )

stability.multistep(alpha = AdamsBashforth$alpha[3,], beta = AdamsBashforth$beta[3,] )

## =============================================================================
## Running a BDF
## =============================================================================
# test model
ode1  <- function (t, y)  return(cos(t)*y )

h     <- 0.01
times <- seq(from = 0, to = 20, by = h)
yout  <- vector (length = length(times))
yout[1] <- 1

# 3rd order BDF
Alpha <- BDF$alpha [3,2:4]
Beta  <- BDF$beta[3,]

bdf <- function(y, t, h, f, ys) {   

  rootfun <- function(ynext) 
    - ynext - sum(Alpha * ys) + Beta * h * f(t + h, ynext)

  y <- multiroot(f=rootfun, start=y)$root

  ys[2:3] <- ys[1:2]
  ys[1]   <- y 

  list (y = y, ys=ys)                                               
}

# Spinup uses Euler...
Euler <- function(y, t, h, f) {
  fn    <- f(t, y) 
  ynext <- y + h * fn

  list (y = ynext, fn = fn)
}

for (i in 2:3)
  yout[i] <- Euler(yout[i-1], times[i-1], h, ode1)$y

ys <- rev(yout[1:3])

# BDF steps
for (i in 4:length(times)){
  step    <- bdf (y=yout[i-1], t=times[i-1], h=h, f=ode1, ys=ys) 
  yout[i] <- step$y
  ys      <- step$ys
} 

}
\seealso{
\link{stability.multistep} for plotting stability regions
}
\keyword{datasets}

