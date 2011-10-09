## =============================================================================
##
## The ring modulator  - 
##
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 15
##
## =============================================================================

ring <- function(times = seq(0, 1e-3, by = 1e-6),
                 yini = NULL, dyini = NULL,
                 parms = list(), method = "mebdfi", maxsteps = 1e6, ...) {

# initial conditions of state variables

    # parameter values
    parameter <- c(c=1.6e-8 , cs=2e-12 , cp=1e-8  , r=25e3 , rp=50,
          lh=4.45   , ls1=2e-3 , ls2=5e-4 , ls3=5e-4,
          rg1=36.3  , rg2=17.3 , rg3=17.3 , ri=50 , rc=600,
          gamma=40.67286402e-9 , delta=17.7493332)

    parameter <- overrulepar(parameter, parms, 16)

### check input 
    if (is.null(yini)) 
      yini <- rep(0,15)
    if (is.null(dyini)) 
      dyini <- rep(0,15)

    checkini(15, yini, dyini)
    
### solve
   if (! is.function(method))
     if (method %in% c("mebdfi", "daspk"))
     return( mebdfi(y = yini, dy = dyini, times = times, res = "ringres",
          dllname = "deTestSet", initfunc = "ringpar",
          parms = parameter,  maxsteps = maxsteps, ...))

   out <- dae(y = yini, times = times,
          func = "ringfunc", 
          dllname = "deTestSet", initfunc = "ringpar", parms = parameter,
          method = method,  maxsteps = maxsteps, ...)

  return(out)
}

# -------------------------------------------------------
