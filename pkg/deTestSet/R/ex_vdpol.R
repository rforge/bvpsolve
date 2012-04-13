## =============================================================================
##
## van der Pol problem
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 2
##
## =============================================================================

vdpol <- function(times = 0:2000, yini = NULL, parms = list(), 
		  atol = 1e-6, rtol = 1e-6, printmescd = TRUE, ...) {

### derivative function
  Vdpol <- function(t,y,mu) {
    list(c(
    y[2],
    mu * (1 - y[1]^2) * y[2] - y[1]
  ))
}

### check input 
    parameter <- c(mu = 1000)

    parameter <- overrulepar(parameter, parms, 1)

    if (is.null(yini))  yini <- c(y1 = 2, y2 = 0) 
    checkini(2, yini)
	prob <- vdpolprob()

### solve
    out <- ode(func = Vdpol, parms = parameter, y = yini, times = times,
			atol=atol,rtol=rtol, ...)
	
	if (nrow(out) > 0) 
	if (printmescd & ( out[nrow(out),1] == prob$t[2] )) { 
		ref = reference("vdpol")
		mescd = min(-log10(abs(out[nrow(out),-1] - ref)/(atol/rtol+abs(ref))))
		printM(prob$fullnm)
		cat('Solved with ')
		printM(attributes(out)$type)
		cat('Using rtol = ')
		cat(rtol)
		cat(', atol=')
		printM(atol)
		printM("Mixed error significant digits:")
		printM(mescd)}
	

    return(out)
}




vdpolprob <- function(){ 
	fullnm <- 'Problem VANDERPOL'
	problm <- 'vdpolm'
	type   <- 'ODE'
	neqn   <- 2
	t <- matrix(1,2)
	t[1]   <- 0
	t[2]   <- 2000
	numjac <- TRUE
	mljac  <- neqn
	mujac  <- neqn	
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}
