################################################################################
#   This is the example for MUSN in U. Ascher, R. Mattheij, and R. Russell,
#   Numerical Solution of Boundary Value Problems for Ordinary Differential
#   Equations, SIAM, Philadelphia, PA, 1995.  MUSN is a multiple shooting
#   code for nonlinear BVPs.  The problem is
#
#      u' =  0.5*u*(w - u)/v
#      v' = -0.5*(w - u)
#      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
#      z' =  0.5*(w - u)
#      y' = -100*(y - w)
#
#   The interval is [0 1] and the boundary conditions are
#
#      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
#
# note: there are two solutions...
################################################################################
require(bvpSolve)

#--------------------------------
# Derivatives
#--------------------------------
musn <- function(t,Y,pars)
{
  with (as.list(Y),
  {
   du=0.5*u*(w-u)/v
   dv=-0.5*(w-u)
   dw=(0.9-1000*(w-y)-0.5*w*(w-u))/z
   dz=0.5*(w-u)
   dy=-100*(y-w)
   return(list(c(du,dv,dw,dz,dy)))
  })
}

x<- seq(0,1,by=0.05)
#----------------------
# Solution method 1
#  ** shooting **
#----------------------

# Residual function for yend...
res  <- function (Y,yini,pars)  with (as.list(Y), w-y)

# Initial values; NA= not available
init <- c(u=1,v=1,w=1,z=-10,y=NA)

print(system.time(
sol   <-bvpshoot(yini= init, x=x,func=musn,
           yend=res,guess=1,atol=1e-10,rtol=0)
))

# second solution...
sol2  <-bvpshoot(yini= init, x=x,func=musn,
           yend=res,guess=0.9,atol=1e-10,rtol=0)

pairs(sol)

# check the solution by simple integration...
yini <- sol[1,-1]
out <- as.data.frame(ode(fun=musn,y=yini,parms=0,times=x,atol=1e-10,rtol=0))
out$w[nrow(out)]-out$y[nrow(out)]

#----------------------
# Solution method 2
#  ** collocation **
#----------------------
# Does not work unless these initial conditions are used
Cost <- function(x)
{
  Sol <- bvptwp(yini= init, x=c(0,1),func=musn, xguess=seq(0,1,len=5),
              yguess=matrix(nr=5,(rep(c(1,1,1,-10,0.91),5)),byrow=TRUE),
           yend=c(NA,NA,NA,NA,x),guess=1,atol=1e-10)
  Sol[nrow(Sol),4]-x
}

print(system.time(
ye  <- multiroot(f=Cost, start=1)
))


# this does not work...
bound <- function(i,y,parms)
{
 if (i == 1) y[1]-1
 else if (i == 2) y[2]-1
 else if (i == 3) y[3]-1
 else if (i == 4) y[4]+10
 else if (i == 5) y[5]-y[3]
}

 yguess=matrix(nr=5,(rep(c(1,1,1,-10,0.91),5)),byrow=TRUE)
 colnames(yguess) <- c("u","v","w","z","y")
# does not work
Sol <- bvptwp(bound=bound,x=x,func=musn, xguess=seq(0,1,len=5),
              leftbc = 4,yguess=yguess)

# but this does...
print(system.time(
Sol <- bvptwp(yini= init, x=c(0,1),func=musn, xguess=seq(0,1,len=5),
              yguess=matrix(nr=5,(rep(c(1,1,1,-10,0.91),5)),byrow=TRUE),
              yend=c(NA,NA,NA,NA,ye$root),guess=1,atol=1e-10)
))
pairs(Sol)

# does not work...
print(system.time(
Sol2 <- bvpcolmod(yini= init, x=c(0,1),func=musn,
                   yend=c(NA,NA,NA,NA,ye$root),guess=1,atol=1e-10)
))