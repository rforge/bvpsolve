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

jac <- function(x,y,p) {
  df <- matrix(nr=3,nc=3,0)
  df[1,2] <- 1
  df[2,1] <- -(y[3]-10*cos(2*x))
  df[2,3] <- -y[1]
  df
}
xguess <-  c(0,1,2*pi)
yguess <- matrix(nr=3,rep(1,9))
rownames(yguess) <- c("y", "dy", "lambda")

# only works for bvptwp if yguess is not 0!...
print(system.time(
sol  <- bvptwp(yini=init,yend=c(NA,0, NA),x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

plot(sol, type="l",lwd=2)

print(system.time(
sol2 <- bvpshoot(yini=init,yend=c(NA,0, NA),x=x,
        func=mathieu2, jacfunc =jac, guess=c(15))
))

plot(sol2, type="l",lwd=2)
  
# including bound
bound <- function(i,y,parms){
  if (i ==1) return(y[1]-1)
  if (i ==2) return(y[2])
  if (i ==3) return(y[2])
}

print(system.time(
sol3  <- bvptwp(bound = bound, leftbc = 2,x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

print(system.time(
sol4  <- bvpshoot(bound = bound, leftbc = 2,x=x,
        func=mathieu2, jacfunc =jac, guess=c(y=1,dy=0,lambda=15))
))
        
## and jacbound
jacbound <- function(i,y,parms){
  if (i ==1) return(c(1,0,0))
  else return(c(0,1,0))
}

print(system.time(
sol5  <- bvptwp(bound = bound, jacbound = jacbound, leftbc = 2, x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

print(system.time(
sol6  <- bvpshoot(bound = bound, jacbound = jacbound, leftbc = 2,x=x,
        func=mathieu2, jacfunc =jac, guess=c(y=1,dy=1,lam=15))
))

par(mfrow=c(2,3))
plot(sol,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol2,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol3,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol4,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol5,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol6,which="y", mfrow=NULL,type="l",lwd=2)
par(mfrow=c(1,1))
