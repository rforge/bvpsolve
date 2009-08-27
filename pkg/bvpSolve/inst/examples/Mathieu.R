## =============================================================================
## Find the 4th eigenvalue of Mathieu's equation:
## y''+(lam-10cos2t)y=0   on the interval [0,pi]
## y(0)=1, y'(0)=0  and y'(pi)=0
##
## 2nd order problem is rewritten as:
## dy=y2
## dy2= -(lam-10cos(2t))*y
## =============================================================================

require(bvpSolve)

mathieu<- function(t,y,lambda=15)
{
 list(c(y[2],-(lambda-10*cos(2*t))*y[1]))
}

#----------------------
# Solution method 1
#  **  shooting  **
#----------------------
x = seq(0,pi,by=0.01)

init <- c(1,0)
sol  <- bvpshoot(yini=init,yend=c(NA,0),x=x,
        func=mathieu, guess=NULL,  extra=15)
plot(sol[,1:2])

#----------------------
# Solution method 2
# multiroot + bvptwp
#----------------------

cost <- function(X)
{  sol<- bvptwp(yini=c(1,NA), yend=c(NA,0),x=c(0,pi),parms=X,
        func=mathieu,guess=0)
  return(sol[2,3])  # y2[0]=0
}

# find the root
lam <- multiroot(f=cost,start=15)

# solve the mode with this root...
Sol<- bvptwp(yini=c(1,NA), yend=c(NA,0),x=x,parms=lam$root,
        func=mathieu,atol=1e-10,guess=1)
lines(Sol,col="red")

#----------------------
# Solution method 3
# augmented equations...
#----------------------


mathieu2<- function(t,y,parms)
{
 list(c(y[2],
        -(y[3]-10*cos(2*t))*y[1],
        0 ))
}

# Solvable with bvpshoot
init <- c(y=1,dy=0,lambda=NA)
sol  <- bvpshoot(yini=init,yend=c(NA,0, NA),x=x,
        func=mathieu2, guess=15)
plot(sol)

# Not solvable with bvptwp

jac <- function(x,y,p) {
  df <- matrix(nr=3,nc=3,0)
  df[1,2] <- 1
  df[2,1] <- -(y[3]-10*cos(2*x))
  df[2,3] <- -y[1]
  df
}

# only works for bvptwp if yguess is not 0!...
sol  <- bvptwp(yini=init,yend=c(NA,0, NA),x=x,
        func=mathieu2, jacfunc =jac, xguess = c(0,1,2*pi),
        yguess = matrix(nr=3,rep(1,9)),guess=1)
plot(sol, type="l",lwd=2)

  
