
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using collocation method "colmod"
##==============================================================================

bvpcolmod<- function(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
    ynames = NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
    posbound=NULL, leftbc = NULL, islin=FALSE, colp=NULL, epsini=NULL,
    eps =epsini, atol=1e-8, nmax=1000, verbose=FALSE, ...)   {

  rho <- environment(func)

  aleft  <- x[1]
  aright <- x[length(x)]
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

  if (! is.null(yini)) {  # yini is either a vector - or a function
    y <- yini
    mstar   <- length(y)       # summed order of the differential equations
    Y       <- y
    inix    <- which (is.na(y))  # NA when initial condition not known
    lini    <- length(inix)
    if (lini > 0 & is.null(guess))  {
      warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
      guess <- rep(0,lini)
    }

    if (lini != length(guess))
      stop("length of 'guess' should be equal to number of NAs in y")
    if (lini > 0)
      y[inix] <- guess
    inix   <- which (!is.na(yini))
    finalx <- which (!is.na(yend))
    zeta   <- c(rep(aleft,length(inix)),rep(aright,length(finalx)))
  } else   {      # the conditions specified by boundary function 'bound'
    mstar   <- length (posbound)    # summed order of the differential equations
    lini    <- mstar
    if (is.null(guess)){
      warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
      guess <- rep(0,mstar)
      }
   Y  <- y <- guess
         
    if (is.null (posbound) && ! is.null(leftbc))
      posbound <- c(rep(aleft,leftbc),rep(aright,mstar-leftbc))
    zeta   <- posbound
  }
  Ynames <- attr(y,"names")
  if (is.null(Ynames) ) Ynames <- ynames

  if (!is.null(eps))  {   #parms AND eps
    Func    <- function(x,state,eps)  {
      attr(state,"names") <- Ynames
      func  (x, state, parms, eps, ...)[1]
    }
      
    Func2   <- function(x, state, eps)  {
      attr(state,"names") <- Ynames
      func   (x, state, parms, eps, ...)
    }
    if (! is.null(jacfunc))
      JacFunc <- function(x, state, eps) {
        attr(state,"names") <- Ynames
        jacfunc(x, state, parms, eps, ...)
      }
    if (! is.null(bound))
      Bound  <- function(i, state, eps) {
        attr(state,"names") <- Ynames
        bound  (i, state, parms, eps, ...)
      }
    if (! is.null(jacbound))
      JacBound   <- function(ii, state, eps) {
        attr(state,"names") <- Ynames
        jacbound(ii, state, parms, eps, ...)
      }
  } else {        #  NO eps
    Func    <- function(x, state, eps)  {
      attr(state,"names") <- Ynames
      func   (x, state, parms, ...)[1]
    }
    Func2   <- function(x, state, eps)  {
      attr(state,"names") <- Ynames
      func   (x, state, parms, ...)
    }
    if (! is.null(jacfunc))
      JacFunc <- function(x, state, eps)  {
        attr(state,"names") <- Ynames
        jacfunc(x, state, parms, ...)
      }
    if (! is.null(bound))
      Bound  <- function(i, state, eps)  {
        attr(state,"names") <- Ynames
        bound   (i, state, parms, ...)
      }
    if (! is.null(jacbound))
      JacBound   <- function(ii, state, eps)  {
        attr(state,"names") <- Ynames
        jacbound(ii, state, parms, ...)
      }
  }

## function evaluation

  tmp <- eval(Func(x[1], y, eps), rho)
  if (!is.list(tmp))
     stop("Model function must return a list\n")
  ncomp  <- length(tmp[[1]])    # number of differential equations
  if (ncomp != mstar)
    stop(paste("The number of function evaluations returned by func() (",
        ncomp, "must equal than the length of the initial conditions vector (",
       length(y), ")", sep = ""))

## in case jacobian function is not defined...

  if ( is.null(jacfunc)) {
    JAC      <- matrix(nr=ncomp,nc=mstar)
    perturbfac  <- 1e-8

    JacFunc <- function (x, state, eps)  {
      state2 <- state
      tmp2   <- unlist(eval(Func(x, state, eps), rho))
      for (i in 1:mstar) {
        dstate   <- max (perturbfac, state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- unlist(eval(Func(x, state, eps), rho))
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

    Bound  <- function(i, state, eps) {
    if (is.function(yini)) {
      y   <-yini(state, parms, ...)
      bb  <- c(y[!is.na(y)], yend[!is.na(yend)])
    }
    return(state[iibb[i]]-bb[i]) }        # too simple for now..
  }

  if ( is.null(jacbound)) {
    JacBound <- function (ii, state, eps)  {
    BJAC        <- numeric(mstar)
    perturbfac  <- 1e-8
      state2 <- state
      tmp2     <- Bound(ii, state, eps)
      for (i in 1:mstar)  {
        dstate   <-max(perturbfac,state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- Bound(ii, state, eps)
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

## the order of the differential equations
## for now i cannot be inputted: order 1 assumed
  order <- rep(1,ncomp)

## legal values
  if (is.null(colp))
    colp=max(max(order) + 1, 5 - max(order))
  if (colp > 7)
    stop ("colp must be smaller than 8")

  ipar   <- rep(0,13)
  itol   <- ncomp
  tol    <- rep(atol,ncomp)
  ipar[1] <- 1          # nonlinear problem
  ipar[2] <- colp       # no of colocation points per subinterval
  ipar[3] <- 0          # no of subintervals in original mesh
  ipar[4] <- 1          # no of tolerances
  if (verbose)
    ipar[7] <- 0 else ipar[7]<-1
  if (! is.null(epsini))
    ipar[13]<-1

#  if (! is.null(xguess))
#   {
#   ipar[9]<-2
#   ipar[3]<-length(xguess)
#   }

# work spaces
  kd      <- ncomp*colp
  kdm     <- kd + mstar
  ipar[6] <- nmax*(kdm + 3)  # length of ispace
  nrec    <- 0 # should be number of right hand side boundary conditions...
  if (islin) {
    nsizef <- 7+3*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar
  } else {
    nsizef <- 8+4*mstar+ (5+kd)*kdm+(2*mstar-nrec)*2*mstar +kd
  }
  ipar[5] <- nmax*nsizef     # length of fspace

  if (is.null(eps))
    eps    <- 1e-10
  if (is.null(epsini))
    epsini <- eps
  if (epsini<eps)
    stop (" eps must be smaller of equal to epsini")

  fixpnt  <- 0
  storage.mode(y) <- storage.mode(x) <- "double"

  out <- .Call("call_colmodsys",as.integer(ncomp),as.integer(mstar),
            as.integer(order),as.double(x),as.double(aleft),as.double(aright),
            as.double(zeta),as.integer(ipar),as.integer(itol),as.double(tol),
            as.double(fixpnt),as.double(epsini),as.double(eps),
            Func, JacFunc, Bound, JacBound, as.double(guess),
            rho, PACKAGE="bvpSolve")

  attr(out,"istate") <- NULL
  attr(out,"rstate") <- NULL
  nm <- c("time",
          if (!is.null(attr(y,"names"))) names(y) else as.character(1:ncomp))
  dimnames(out) <- list(nm,NULL)
  t(out)
}
