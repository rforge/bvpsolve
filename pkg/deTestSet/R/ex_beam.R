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


beam <- function(times = seq(0, 5, by = 0.05), yini = NULL, method = gamd, 
		printmescd = TRUE,atol=1e-6,rtol=1e-6, ...) {

### check input 

    # there are no parameters...
	prob <- beamprob()
	
    if (is.null(yini)) yini <- rep(0, 80)

    checkini(80, yini)

### solve 
    out <- ode(func = "beamfunc", parms = NULL, dllname = "deTestSet", y = yini,
           times = times, initfunc = NULL, method=method, atol=atol,rtol=rtol,...)
   
   if (printmescd & ( out[nrow(out),1] == prob$t[2] )) { 
	   ref = reference("beam")
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

beamprob <- function(){
	fullnm <- "Beam"
	problm <- 'beam'
	type <- 'ODE'
	neqn <- 80
	ndisc <- 0
	t <- matrix(1,2)
	t[1] <- 0
	t[2] <- 5
	numjac <- TRUE
	mljac <- neqn
	mujac <- neqn
	return(list(fullnm=fullnm, problm=problm,type=type,neqn=neqn,ndisc=ndisc,
					t=t,numjac=numjac,mljac=mljac,mujac=mujac))
}

