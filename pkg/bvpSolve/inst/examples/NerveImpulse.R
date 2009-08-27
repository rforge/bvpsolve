## =============================================================================
##   This is a nerve impulse model considered in Example 7.1 of
##   R. Seydel, From Equilibrium to Chaos, Elsevier, 1988.  The
##   differential equations
##
##      y1' = 3*t*(y1 + y2 - 1/3*y1^3 + lambda)
##      y2' = -t*(y1 - 0.7 + 0.8*y2)/3
##
##   are to be solved subject to periodic boundary conditions
##
##      y1(0) = y1(1)
##      y2(0) = y2(1)
##
##
##   The parameter lambda has the value -1.3.  The
##   period T is unknown, i.e the length of the interval is not known.
##
##   An extra equation to estimate T is: -T(y1(0)-0.7+0.8*y2(0))/3=1
##
## =============================================================================


require(bvpSolve)

nerve <- function (t,y,T)
  list(c( 3*T*(y[1] + y[2] - 1/3*(y[1]^3) - 1.3),
        (-1/3)*T*(y[1] - 0.7 + 0.8*y[2]) ))

res<- function (Y,yini,T)
  c(Y[1]-yini[1],
    Y[2]-yini[2],
    T*(-1/3)*(yini[1] - 0.7 + 0.8*yini[2]) - 1)
  
init <- c(y1=NA,y2=NA)
sol  <- bvpshoot(yini=init,x=seq(0,1,by=0.01),
        func=nerve, guess=c(0.5,0.5), yend=res, extra=2*pi)
attributes(sol)$root
plot(sol,type="l",lwd=2)

## =============================================================================
## Second solution: augment the number of state variables...
## =============================================================================

nerve2 <- function (t,y,p)
  list(c( 3*y[3]*(y[1] + y[2] - 1/3*(y[1]^3) - 1.3),
        (-1/3)*y[3]*(y[1] - 0.7 + 0.8*y[2]) ,
         0))

res2<- function (Y,yini,T)
  c(Y[1]-yini[1],
    Y[2]-yini[2],
    Y[3]*(-1/3)*(yini[1] - 0.7 + 0.8*yini[2]) - 1)

init <- c(y1=NA,y2=NA,T=NA)
sol2 <- bvpshoot(yini=init,x=seq(0,1,by=0.01),
        func=nerve2, guess=c(0.5,0.5,2*pi), yend=res2)
        
attributes(sol)$root
plot(sol,type="l",lwd=2)

yini <- NULL
## =============================================================================
## Third solution: augment + bvptwp - does not work...
## =============================================================================

bound <- function(i,y,p) {
 if (i==1) {
  yini <<- y
  return(y[3]*(-1/3)*(y[1] - 0.7 + 0.8*y[2]) - 1)
 }
 if (i ==2) return(y[1]-yini[1])
 if (i ==3) return(y[2]-yini[2])
 
}

Sol  <- bvptwp(bound=bound, x=seq(0,1,by=0.01), guess=c(y=0.5,dy=0.5,T=2*pi),
        func=nerve2,leftbc=1)
plot(sol,type="l",lwd=2)
Sol  <- bvptwp(bound=bound, x=seq(0,1,by=0.01),
        func=nerve2,leftbc=1,xguess=xguess,yguess=yguess)
