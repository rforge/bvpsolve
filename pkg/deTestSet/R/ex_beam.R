## =============================================================================
##
## Beam problem 
##
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 80
##
## =============================================================================


beam <- function(times = seq(0,5,by=0.01), yini = NULL, 
    ...) {

### check input 

    # there are no parameters...

    if (is.null(yini)) yini <- rep(0,80)

    checkini(80,yini)

### solve 
    out <- ode(func="beamf", parms=NULL, dllname="deTestSet", y = yini, 
           times=times, initfunc=NULL,  ...)

  return(out)
}

