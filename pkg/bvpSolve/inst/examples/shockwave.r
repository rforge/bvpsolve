#==============================================================================
# "Shock Wave", a test problem in Ascher, Mattheij, and
# Russell, Numerical Solution of Boundary Value Problems
# for Ordinary Differential Equations", Classics in Applied
# Mathematics Series, SIAM, Philadelphia, 1995.
#==============================================================================
require(bvpSolve)

# Declare global problem dependent parameters.
eps     <- 0.01
gamma   <- 1.4

neqns   <- 2   #number of differential equations
leftbc  <- 1   #number of boundary conditions at left end of the interval

fsub <- function (t,y,pars)
{ term1 <- 1/(eps*(1.+t*t))
  term2 <- 0.5 +0.5*gamma - 2.0*eps*t
  f1    <- y[2]
  f2    <- term1/y[1]*(term2*y[1]*y[2] - y[2]/y[1] -
           (2.0*t/(1.0+t^2))*(1.0-0.5*(gamma-1.0)*y[1]^2) )
  return(list(c(f1,f2)))
}

A      <- 0.
B      <- 1.
x     <- seq(A,B,len=100)

# bvpshoot does not work for eps=0.01
#print(system.time(Sol <- bvpshoot (yini=c(0.9129,NA), yend=c(0.375,NA),
#       x=x,func=fsub,guess=0,atol=1e-10,rtol=1e-10)))
#plot(Sol)

# bvptwp does not work either...
#print(system.time(Sol <- bvptwp (yini=c(0.9129,NA), yend=c(0.375,NA),
#       x=x,func=fsub,guess=0,atol=1e-10)))
#plot(Sol)


#bvpcolmod works but is more complex.
fsub2<- function (t,y,pars,eps)
{ term1 <- 1/(eps*(1.+t*t))
  term2 <- 0.5 +0.5*gamma - 2.0*eps*t
  f1    <- y[2]
  f2    <- term1/y[1]*(term2*y[1]*y[2] - y[2]/y[1] -
           (2.0*t/(1.0+t^2))*(1.0-0.5*(gamma-1.0)*y[1]^2) )
  return(list(c(f1,f2)))
}
# Solve with simple interface
print(system.time(Sol <- bvpcolmod (yini=c(0.9129,NA), yend=c(0.375,NA),
        x=x,func=fsub2,guess=0.5,eps= 0.001)))
plot(Sol, type="l", lwd=2)

# More complex but faster solution methods with analytical jacobian
df <- function(t,y,pars,eps)       # jacobian
{
  pd <-matrix(nr=neqns,nc=neqns,data=0.0)
  term1 <- 1/(eps*(1.+t*t))
  term2 <- 0.5 +0.5*gamma - 2.0*eps*t

  pd[1,1] = 0.0
  pd[1,2] = 1.0

  pd[2,1] = term1*( 2.0*y[2]/(y[1]^3) + 2.0*t/
	        ((1.0+t^2)*y[1]^2) + (t/(1.0+t^2))*(gamma-1.0) )

  pd[2,2] = (term1/y[1])*( term2*y[1] - 1.0/y[1] )
  return(pd)
}

print(system.time(Sol <- bvpcolmod (yini=c(0.9129,NA), yend=c(0.375,NA),
        x=x,func=fsub2,guess=0.5,eps= 0.001,jacfunc=df)))
        
lines(Sol,col="blue",lty=2)

