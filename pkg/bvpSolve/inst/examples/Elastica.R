## =============================================================================
## Elastica equation
## as from Jeff Cash's webpage 
## http://www.ma.ic.ac.uk/~jcash/BVP_software
## =============================================================================

Elastica <- function (t, y, pars) {

  list( c(cos(y[3]),
          sin(y[3]),
          y[4],
          y[5]*cos(y[3]),
          0))
}

require(bvpSolve)

Sol <- bvptwp(func=Elastica,
              yini = c(x=0, y=0, p=NA,   k=0, F=NA),
              yend = c(x=NA,y=0, p=-pi/2,k=NA,F=NA),
              x = seq(0,0.5,len=16),
              guess=c(0,0) )

Sol
plot(Sol)
