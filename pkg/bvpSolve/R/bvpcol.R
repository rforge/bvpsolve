
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using collocation method "colnew"
##==============================================================================

bvpcol<- function(yini=NULL, x, func, yend=NULL, parms=NULL, ynames = NULL,
    xguess=NULL, yguess=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
    leftbc = NULL, islin=FALSE, nmax=1000, atol=1e-8, ncomp=NULL,
    colp=NULL, posbound=NULL, order = NULL, fullOut = TRUE, 
    dllname=NULL, initfunc=dllname, rpar = NULL, ipar = NULL, nout = 0,
     forcings=NULL, initforc = NULL, fcontrol=NULL, verbose = FALSE,...)   {

  rho <- environment(func)

  aleft  <- x[1]
  aright <- x[length(x)]
  fixpt  <- x[-c(1,length(x))]
  
##---------------------
## error checking...
##---------------------
  if (is.null(yini)   && is.null(bound))
    stop("either 'yini' and 'yend' or 'bound' should be inputted")
  if (!is.null(bound) && is.null(posbound)&& is.null(leftbc))
    stop("if 'bound' is given, 'posbound' or 'leftbc' should also be given")
  if (!is.null(yini)  && is.null(yend))
    stop("if 'yini' is given, 'yend' should also be given")
  if (!is.null(yini)  && !is.null(bound))
    stop("either 'yini' or bound should be given, not both")

##---------------------
## Boundary conditions (specified by yini and yend)
##---------------------
  mstar <- NULL
  y     <- NULL
  guess <- NULL
  Restart <- "bvpSolve" %in% class(yguess)    

  if (! is.null(yini)) {  # yini is either a vector - or a function
    y      <- yini
    mstar  <- length(y)       # summed order of the differential equations
    if (! is.null(order)) {
      ncomp <- length (order)
      if (ncomp > mstar)
        stop ("length of 'order' can not exceed the length of 'yini'")
      if (ncomp <= 0)
        stop ("length of 'order' should be > 0")
      if (sum(order) != mstar)
        stop ("summed values of 'order' should be = length of 'yini' = number of variables")  
    } else {
      ncomp  <- length(yini) 
      order <- rep(1,ncomp)
    }

    Y      <- y
    inix    <- which (is.na(y))  # NA when initial condition not known
    nas     <- length(inix)
    leftbc  <- length(which (!is.na(y)))   # NA when initial condition not known

    if ( ! is.null(yguess) & ! Restart)
      guess <- yguess[inix,1]

    if (nas > 0 )  {
      guess <- rep(0,nas)
    }

    if (nas > 0)
      y[inix] <- guess

    inix   <- which (!is.na(yini))
    finalx <- which (!is.na(yend))
    ll <- length(finalx) + length(inix)
    if (ll != mstar)
      stop (paste("number of boundary conditions wrong: should be ",mstar," but is", ll))
    zeta   <- c(rep(aleft,length(inix)),rep(aright,length(finalx)))

  } else   {      # the conditions specified by boundary function 'bound'
    if (is.null(leftbc) & is.null(posbound))
       stop("'leftbc' or 'posbound' should be inputted if 'bound' is given")
    if (! is.null(order)) {
      ncomp <- length (order)
      mstar <- sum(order)
      if (ncomp > mstar)
        stop ("length of 'order' can not exceed the length of 'yini'")
      if (ncomp <= 0)
        stop ("length of 'order' should be > 0")
    }  
    if (! is.null(posbound))
       mstar   <- length (posbound)    # summed order of the differential equations
    if (! is.null(yguess) & ! Restart) {
      guess <- yguess[,1]   
      if (is.null(mstar))
       mstar <- nrow(yguess)
      else if (mstar != nrow(yguess))
       stop ("number of rows in 'yguess' should = number of variables")  
    }  
   if (! is.null(ynames) & is.null(ncomp)) ncomp <- length(ynames) 
    if (is.null(mstar))
       mstar <- ncomp
    if (is.null(mstar))
      stop ("don't know the number of equations - provide 'ncomp' ")
       
    nas    <- mstar
    if (is.null(guess)){
      if (verbose) warning("estimates for unknown initial conditions not given ('xguess','yguess'); assuming 0's")
      guess <- rep(0,mstar)
      }
   Y  <- y <- guess
         
   if (is.null (posbound) && ! is.null(leftbc)) {
     posbound <- NULL
     if (leftbc >0)          posbound <- c(posbound,rep(aleft,leftbc))
     if (mstar - leftbc >0)  posbound <- c(posbound,rep(aright,mstar-leftbc))
   } 
   zeta   <- posbound
  }
  Ynames <- ynames
  if (is.null(Ynames)) Ynames <- attr(y,"names")
  if (is.null(Ynames) & ! is.null(yini))   Ynames <- names(yini)
  if (is.null(Ynames) & ! is.null(yend))   Ynames <- names(yend)
  if (is.null(Ynames) & is.matrix(yguess)) Ynames <- rownames(yguess)
  if (is.null(Ynames) & is.vector(yguess)) Ynames <- names(yguess)

  if (is.null(ncomp)) ncomp <- mstar 
  if (is.null(order)) order <- rep(1,ncomp)

  # Functions
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  # The functions are in a DLL
  if (is.character(func)) {
    if (sum(duplicated (c(func,initfunc,jacfunc))) >0)
      stop("func, initfunc, or jacfunc cannot be the same")

    if (! is.null(initfunc))  # KS: ADDED THAT to allow absence of initfunc
      if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
        stop(paste("'initfunc' not loaded ",initfunc))

    # Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA

    funcname <- func
    ## get the pointer and put it in func
    if (is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("bvp function not loaded",funcname))

    ## the Jacobian
    if (is.loaded(jacfunc, PACKAGE = dllname))  {
      JacFunc <- getNativeSymbolInfo(jacfunc, PACKAGE = dllname)$address
    } else
      stop(paste("jacobian function not loaded ",jacfunc))

    ## the boundary
    if (is.loaded(bound, PACKAGE = dllname))  {
      Bound <- getNativeSymbolInfo(bound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary function not loaded ",bound))

    ## the boundary Jacobian
    if (is.loaded(jacbound, PACKAGE = dllname))  {
      JacBound <- getNativeSymbolInfo(jacbound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary jac function not loaded ",jacbound))
    Nglobal <-  0
    Nmtot <- NULL

    if (! is.null(forcings))
      flist <- checkforcings(forcings,x,dllname,initforc,FALSE,fcontrol)

  } else {      # The functions are R-code
    
    Func    <- function(x, state)  {
        attr(state,"names") <- Ynames
        func   (x, state, parms, ...)[1]
      }
    Func2   <- function(x, state)  {
        attr(state,"names") <- Ynames
        func   (x, state, parms, ...)
      }
    if (! is.null(jacfunc))
      JacFunc <- function(x, state)  {
          attr(state,"names") <- Ynames
          jacfunc(x, state, parms, ...)
        }
    if (! is.null(bound))
        Bound  <- function(i, state)  {
          attr(state,"names") <- Ynames
          bound   (i, state, parms, ...)
        }
    if (! is.null(jacbound))
        JacBound   <- function(ii, state)  {
          attr(state,"names") <- Ynames
          jacbound(ii, state, parms, ...)
        }

## function evaluation            
  if (is.null(y) & ! is.null(ncomp)) y<-runif(ncomp) 
  tmp <- eval(Func2(x[1], y), rho)
  if (!is.list(tmp))
     stop("Model function must return a list\n")
  NN  <- length(tmp[[1]])    # number of differential equations
  if (NN > 20)
    stop ("number of equations must be <= 20 for bvpcol")
  if (NN != ncomp)
    stop(paste("The number of function evaluations returned by func() (",
        NN, ") must equal than the number of variables (",ncomp, ")", sep = ""))

## in case jacobian function is not defined...
  if ( is.null(jacfunc)) {
    JAC      <- matrix(nr=ncomp,nc=mstar)
    perturbfac  <- 1e-8

    JacFunc <- function (x, state)  {
      state2 <- state
      tmp2   <- unlist(eval(Func(x, state), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- unlist(eval(Func(x, state), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }
  }

## in case boundary function is not defined...
  if ( is.null(bound)) {  
    iini <- sum(!is.na(Y))  
    iend <- sum(!is.na(yend))
    iibb <- c(which( !is.na(Y)), which(!is.na(yend)))
    bb   <- c(Y[!is.na(Y)], yend[!is.na(yend)])

    Bound  <- function(i, state) {
    if (is.function(yini)) {
      y   <-yini(state, parms, ...)
      bb  <- c(y[!is.na(y)], yend[!is.na(yend)])
    }
    return(state[iibb[i]]-bb[i]) }        # too simple for now..
  }

  if ( is.null(jacbound)) {
    JacBound <- function (ii, state)  {
    BJAC        <- numeric(mstar)
    perturbfac  <- 1e-8
      state2 <- state
      tmp2     <- Bound(ii, state)
      for (i in 1:mstar)  {
        dstate   <-max(perturbfac,state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- Bound(ii, state)
        BJAC[i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(BJAC)
    }
  }

  Nglobal <- if (length(tmp) > 1)
        length(unlist(tmp[-1]))
  else 0
  Nmtot <- attr(unlist(tmp[-1]), "names")
  }
  
## legal values
  if (is.null(colp))
    colp=max(max(order) + 1, 5 - max(order))
  if (colp > 7)
    stop ("colp must be smaller than 8")

  iset   <- rep(0,14)
  tol    <- rep(atol,ncomp)
  iset[1] <- 1-islin    # linear/nonlinear problem
  iset[2] <- colp       # no of colocation points per subinterval
  iset[3] <- 0          # no of subintervals in original mesh
  iset[4] <- ncomp      # no of tolerances
  iset[14] <- fullOut
  
  if (verbose)
    iset[7] <- 0 else iset[7]<-1

  GuessFunc <- function(x) y     # dummy function
  if (Restart) {   # previous solution used - same mesh 
    ATT <- attributes(yguess)
    if (ATT$name != "bvpcol")
       stop ("can only continuate solution if previous solution in 'yguess' was obtained with 'bvpcol', not with ", ATT$name)
    rwork <- ATT$rstate
    if ( length (rwork) == 0)
       stop ("attributes(yguess)$rstate should be present, if continuation is requested")
    iwork <- ATT$istate[-1]   
    if ( length (iwork) == 0)
       stop ("attributes(yguess)$istate should be present, if continuation is requested")
             
    iset[9]<-2
    iset[3]<- iwork[1]
  } else {
  rwork <- 0.
  iwork <- 0
  
  Xguess <- xguess
  Yguess <- yguess
# check dimensions
  if (! is.null(yguess) && ! is.null(yguess))
     if (length(xguess) != ncol(yguess))
       stop("yguess should have as many columns as length of xguess")
       
  if (! is.null (xguess))  {
    givmesh <- TRUE
    nmesh   <- length(xguess)
    
      rr <- range(Xguess) - range(x)
      if (rr[1] <0)
        stop ("minimum of 'xguess' (",Xguess[1], "), should be <= minimum of 'x' (", x[1], ")")
      if (rr[2] <0)
        stop ("maximum of 'xguess' (",Xguess[nmesh], "), should be >= maximum of 'x' (", x[nmesh], ")")
 
  }
  if (! is.null(yguess))  {
    if (is.null(xguess))
      stop ("xguess must be provided if yguess is")
    givu <- TRUE
    Yguess <- yguess # ncomp,nmesh
    if (length(Yguess) != nmesh*mstar) stop ("xguess and yguess not compatible")
    GuessFunc <- function(x) {
      zz <- NULL
      for (i in 1:ncomp) 
        zz <- rbind(zz,approx(Xguess,Yguess[i,],x, rule=2)$y)
      zz  
    }
  }
   
  if (! is.null(Xguess))   {
   iset[9]<-1
   iset[3]<-length(xguess)
  }
  }
  
# work spaces
  kd      <- ncomp*colp
  kdm     <- kd + mstar
  iset[6] <- nmax*(kdm + 3)  # length of ispace

  nrec    <- 0 # should be number of right hand side boundary conditions...
  nsizef <- 4 + 3*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar
  iset[5] <- nmax*nsizef     # length of fspace
  iset[5] = max(iset[5], length(rwork))     # length of fspace

  iset[10] <- 0
  if (is.null(ipar)) ipar <- 1
  if (is.null(rpar)) rpar <- 1
  
  ii <- which (zeta %in% c(aleft, aright))
  fixpnt   <- zeta[-ii]  # points that have to be excluded
  iset[11] <- length(fixpnt)
  if(is.null(initfunc))
    initpar <- NULL # parameter initialisation not needed if function is not a DLL
  else
    initpar <- as.double(parms)

  storage.mode(y) <- storage.mode(x) <- "double"

  out <- .Call("call_colnew",as.integer(ncomp),as.double(x),
            as.double(aleft),as.double(aright),as.double(zeta),
            as.integer(mstar),as.integer(order),
            as.integer(iset),as.double(rwork), as.integer(iwork),
            as.double(tol),
            as.double(fixpnt), as.double (rpar), as.integer (ipar), 
            Func, JacFunc, Bound, JacBound, GuessFunc, ModelInit, initpar,
            flist, type = as.integer(2), rho,  PACKAGE="bvpSolve")

#  attr(out,"istate") <- NULL
#  attr(out,"rstate") <- NULL
     nm <- c("x",
          if (!is.null(Ynames)) Ynames else as.character(1:mstar))
# if there are other variables...
    if (Nglobal > 0) {
       if (!is.character(func)) {                  # if a DLL: already done...    
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2)) {
          y <- out[-1,i]
          names(y) <- nm[-1]
          out2[, i] <- unlist(Func2(out[1, i], y)[-1])  # KS: Func2 rather than func
          }
        out <- rbind(out,out2)
        }  # end !is.character func
        nm <- c(nm,
                if (!is.null(Nmtot)) Nmtot else
                                     as.character((ncomp+1) : (ncomp + Nglobal)))
  } 
  class(out) <- c("bvpSolve","matrix")  # a boundary value problem
  attr(out,"name") <- "bvpcol"
  dimnames(out) <- list(nm,NULL)
  t(out)
}
