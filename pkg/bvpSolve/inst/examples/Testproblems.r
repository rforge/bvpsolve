#===============================================================================
# The test problems from
# http://www.ma.ic.ac.uk/~jcash/BVP_software
#
# See also the \inst\Test directory
#===============================================================================

###################################################
### chunk number 1: preliminaries
###################################################
library("bvpSolve")


#===============================================================================
### Test problem 1
#===============================================================================

Prob1 <- function(t, y, pars) {
   list(c( y[2] , y[1]/xi ))
}

xi <-0.1
print(system.time(
  shoot  <- bvpshoot(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
            func=Prob1, guess=0)))
print(system.time(
  twp <- bvptwp(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
            func=Prob1, guess=0)))
print(system.time(
  mod <- bvpcolmod(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
            func=Prob1, guess=0)))

xi <-0.01
shoot2  <- bvpshoot(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
            func=Prob1, guess=0)

xi <-0.001
shoot3  <- bvpshoot(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob1, guess=0)
col2  <- bvpcolmod(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob1, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 1",xlab="x",ylab="y")
lines(shoot2[,1],shoot2[,2],lty=2)
lines(shoot3[,1],shoot3[,2],lty=3)
# exact solution
curve(exp(-x/sqrt(xi))-exp((x-2)/sqrt(xi))/(1-exp(-2/sqrt(xi))),
      0,1,add=TRUE,type="p")


#===============================================================================
### Test problem 2
#===============================================================================

Prob2 <- function(t, y, pars) {
  list(c( y[2], y[2]/xi ))
}

xi <-0.2
shoot <- bvpshoot(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob2, guess=0)

xi <-0.1
twp   <- bvptwp(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob2, guess=0, atol=1e-10)
xi <- 0.01
mod   <- bvpcolmod(yini=c(1,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob2, guess=0, atol=1e-10)


plot(shoot[,1],shoot[,2],type="l",main="test problem 2",xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)
curve((1-exp((x-1)/xi))/(1-exp(-1/xi)),0,1,type="p",add=TRUE)


#===============================================================================
### Test problem 3
#===============================================================================

Prob3 <- function(x, y, pars) {
  list(c( y[2],
         1/xi*(-(2+cos(pi*x))*y[2]+y[1]-
          (1+xi*pi*pi)*cos(pi*x)-(2+cos(pi*x))*pi*sin(pi*x))
      ))
}

xi <-0.1
shoot  <- bvpshoot(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
           func=Prob3, guess=0)

xi <-0.01
twp    <- bvptwp(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
           func=Prob3, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 3",xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
curve(cos(pi*x),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 4
#===============================================================================

Prob4 <- function(t, y, pars) {
  list(c( y[2], (-y[2]+(1+xi)*y[1])/xi ))
}

yini <- c(1+exp(-2),NA)

xi   <-0.5
yend <- c(1+exp(-2*(1+xi)/xi),NA)
shoot  <- bvpshoot(yini=yini,yend=yend,x=seq(-1,1,by=0.01),
          func=Prob4, guess=0)

xi <-0.1
yend <- c(1+exp(-2*(1+xi)/xi),NA)
twp    <- bvptwp(yini=yini,yend=yend,x=seq(-1,1,by=0.01),
          func=Prob4, guess=0)

xi <-0.01
yend <- c(1+exp(-2*(1+xi)/xi),NA)
mod    <- bvpcolmod(yini=yini,yend=yend,x=seq(-1,1,by=0.01),
          func=Prob4, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 4",ylim=c(0,1.2),
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)
curve(exp(x-1)+exp(-(1+xi)*(1+x)/xi),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 5
#===============================================================================

Prob5 <- function(x, y, pars) {
  list(c( y[2],
          x*y[2]+y[1]-(1+pi*pi)*cos(pi*x)+pi*x*sin(pi*x) ))
}

xi <-0.1
shoot  <- bvpshoot(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
           func=Prob5, guess=0)

twp    <- bvptwp(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
           func=Prob5, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 5",
  xlab="x",ylab="y")
curve(cos(pi*x),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 6
#===============================================================================

Prob6 <- function(t, y, pars) {
  list(c( y[2],
          1/xi*(-t*y[2]-xi*pi*pi*cos(pi*t)-pi*t*sin(pi*t)) ))
}

xi <-0.1
shoot  <- bvpshoot(yini=c(-2,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
           func=Prob6, guess=0)

xi <-0.01
twp  <- bvptwp(yini=c(-2,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
           func=Prob6, guess=0)

xi <-0.001
mod  <- bvpcolmod(yini=c(-2,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
           func=Prob6, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 6",
     ylim =c(-2,2),xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x)+erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 7
#===============================================================================

prob7 <- function(x, y, pars) {
  list(c(y[2],
         1/xi*(-x*y[2]+y[1]-(1+xi*pi*pi)*cos(pi*x)-pi*x*sin(pi*x)))
     )
}

x  <- seq(-1,1,by=0.01)

xi <- 0.01
twp  <- bvptwp(yini=c(-1,NA),yend=c(1,NA),x=x,func=prob7, guess=0)

xi <- 0.001
mod  <- bvpcolmod(yini=c(-1,NA),yend=c(1,NA),x=x,func=prob7, guess=0)

xi <- 0.0001
twp2  <- bvptwp(yini=c(-1,NA),yend=c(1,NA),x=x,func=prob7, guess=0)


plot(twp[,1],twp[,2],type="l",main="test problem 7",
  xlab="x",ylab="y")
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x)+x+(x*erf(x/sqrt(2*xi))+sqrt(2*xi/pi)*exp(-x^2/2/xi))/
         (erf(1/(2*xi))+sqrt(2*xi/pi)*exp(-1/2/xi)),-1,1,type="p",add=TRUE)

plot(twp[,1],twp[,3],type="l",main="test problem 7",
  xlab="x",ylab="y'")
lines(mod[,1],mod[,3],lty=2)
lines(twp2[,1],twp2[,3],lty=3)


#===============================================================================
### Test problem 8
#===============================================================================

prob8 <- function(x, y, pars) {
  list(c( y[2], -1/xi*y[2]))
}

x  <- seq(0,1,by=0.01)
xi <- 0.2
shoot  <- bvpshoot(yini=c(1,NA),yend=c(2,NA),x=x,func=prob8,guess=0)

xi <- 0.1
twp  <- bvptwp(yini=c(1,NA),yend=c(2,NA),x=x,func=prob8, guess=0)

xi <- 0.01
mod  <- bvpcolmod(yini=c(1,NA),yend=c(2,NA),x=x,func=prob8, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 8",
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)

curve(2-exp(-1/xi)-exp(-x/xi)/(1-exp(-1/xi)),0,1,add=TRUE,type="p")


#===============================================================================
### Test problem 9
#===============================================================================

Prob9 <- function(x, y, pars) {
  list(c( y[2], -1/(xi+x^2)*(4*x*y[2]+2*y[1]) ))
}

xi <-0.05
twp  <- bvptwp(yini=c(1/(1+xi),NA),yend=c(1/(1+xi),NA),x=seq(-1,1,by=0.01),
           func=Prob9, guess=0)

xi <-0.02
twp2  <- bvptwp(yini=c(1/(1+xi),NA),yend=c(1/(1+xi),NA),x=seq(-1,1,by=0.01),
           func=Prob9, guess=0)

xi <-0.01
mod  <- bvpcolmod(yini=c(1/(1+xi),NA),yend=c(1/(1+xi),NA),x=seq(-1,1,by=0.01),
           func=Prob9, guess=0)


plot(twp[,1],twp[,2],type="l",main="test problem 9", ylim=c(0,100),
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)

curve(1/(xi+x^2),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 10
#===============================================================================

Prob10 <- function(x, y, pars) {
  list(c( y[2], -1/xi*x*y[2] ))
}

xi <-0.1
shoot  <- bvpshoot(yini=c(0,NA),yend=c(2,NA),x=seq(-1,1,by=0.01),
        func=Prob10, guess=0,atol=1e-10)

xi <- 0.05
twp  <- bvptwp(yini=c(0,NA),yend=c(2,NA),x=seq(-1,1,by=0.01),
           func=Prob10, guess=0)

xi <- 0.01
mod  <- bvpcolmod(yini=c(0,NA),yend=c(2,NA),x=seq(-1,1,by=0.01),
           func=Prob10, guess=0)


plot(shoot[,1],shoot[,2],type="l",main="test problem 10",
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(1+erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 11
#===============================================================================

Prob11 <- function(x, y, pars) {
  list(c(y[2], 1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x)) ))
}

xi <-0.1
# Shooting
print(system.time(
shoot  <- bvpshoot(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob11, guess=0,atol=1e-10)
))

print(system.time(
twp  <- bvptwp(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob11, guess=0,atol=1e-10)
))
print(system.time(
mod  <- bvpcolmod(yini=c(-1,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob11, guess=0,atol=1e-10)
))


plot(shoot[,1],shoot[,2],type="l",main="test problem 11",
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)

curve(cos(pi*x),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 12
#===============================================================================

Prob12 <- function(x, y, pars) {
  list(c(y[2],1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x))))
}

xi <-0.01
shoot  <- bvpshoot(yini=c(-1,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob12, guess=0,atol=1e-10)

xi <-0.0025
twp  <- bvptwp(yini=c(-1,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob12, guess=0,atol=1e-10)

xi <-0.0001
mod  <- bvpcolmod(yini=c(-1,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob12, guess=0,atol=1e-10)


plot(shoot[,1],shoot[,2],type="l",main="test problem 12",
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)
curve(cos(pi*x)+exp((x-1)/sqrt(xi)),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 13
#===============================================================================

Prob13 <- function(x, y, pars)  {
  list(c( y[2], 1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x)) ))
}

xi <-0.01
shoot  <- bvpshoot(yini=c(0,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob13, guess=0,atol=1e-10)

xi <-0.0025
twp  <- bvptwp(yini=c(0,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob13, guess=0,atol=1e-10)

xi <-0.0001
mod  <- bvpcolmod(yini=c(0,NA),yend=c(-1,NA),x=seq(-1,1,by=0.01),
        func=Prob13, guess=0,atol=1e-10)


plot(shoot[,1],shoot[,2],type="l",main="test problem 13",
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)
curve(cos(pi*x)+exp(-(x+1)/sqrt(xi)),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 14
#===============================================================================

Prob14 <- function(x, y, pars)  {
  list(c( y[2], 1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x))))
}

xi <-0.01
shoot  <- bvpshoot(yini=c(0,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob14, guess=0,atol=1e-10)

xi <-0.0025
twp  <- bvptwp(yini=c(0,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob14, guess=0,atol=1e-10)

xi <-0.0001
mod  <- bvpcolmod(yini=c(0,NA),yend=c(0,NA),x=seq(-1,1,by=0.01),
        func=Prob14, guess=0,atol=1e-10)


plot(shoot[,1],shoot[,2],type="l",main="test problem 14", ylim=c(-1,1),
  xlab="x",ylab="y")
lines(mod[,1],mod[,2],lty=2)
lines(twp[,1],twp[,2],lty=3)
curve(cos(pi*x)+exp((x-1)/sqrt(xi))+exp(-(x+1)/sqrt(xi)),-1,1,type="p",add=TRUE)


#===============================================================================
### Test problem 15
#===============================================================================

Prob15 <- function(x, y, pars)  {
  list(c( y[2], 1/xi*x*y[1] ))
}

xi <-0.003
# Shooting
print(system.time(
shoot  <- bvpshoot(yini=c(1,NA),yend=c(1,NA),x=seq(-1,1,by=0.01),
           func=Prob15, guess=0,atol=1e-10)
))

xi <- 0.005

print(system.time(
twp  <- bvptwp(yini=c(1,NA),yend=c(1,NA),x=seq(-1,1,by=0.01),
           func=Prob15, guess=0)
))

xi <- 0.01
print(system.time(
mod  <- bvpcolmod(yini=c(1,NA),yend=c(1,NA),x=seq(-1,1,by=0.01),
           func=Prob15, guess=0)
))


plot(shoot[,1],shoot[,2],type="l",main="test problem 15",
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)


#===============================================================================
### Test problem 16
#===============================================================================

Prob16 <- function(x, y, pars) {
  list(c( y[2], -1/xi^2*pi^2*y[1]/4 ))
}

xi <-0.11
# Shooting
print(system.time(
shoot  <- bvpshoot(yini=c(0,NA),yend=c(sin(pi/2/xi),NA),x=seq(0,1,by=0.01),
           func=Prob16, guess=0,atol=1e-10)
))


plot(shoot[,1],shoot[,2],type="l",main="test problem 16",
  xlab="x",ylab="y")
curve(sin(pi*x/2/xi),0,1,type="p",add=TRUE)


#===============================================================================
### Test problem 17
#===============================================================================

Prob17 <- function(x, y, pars)  {
  list(c( y[2], -3*xi*y[1]/(xi+x^2)^2 ))
}

xseq<-seq(-0.1,0.1,by=0.001)

xi <-0.01
twp1  <- bvptwp(yini=c(-0.1/sqrt(xi+0.01),NA),
                   yend=c(0.1/sqrt(xi+0.01),NA),x=xseq,
                   func=Prob17, guess=0,atol=1e-10)

xi <- 0.001
twp2  <- bvptwp(yini=c(-0.1/sqrt(xi+0.01),NA),
                   yend=c(0.1/sqrt(xi+0.01),NA),x=xseq,
                   func=Prob17, guess=0,atol=1e-8)

xi <- 0.0001
twp3  <- bvptwp(yini=c(-0.1/sqrt(xi+0.01),NA),
                   yend=c(0.1/sqrt(xi+0.01),NA),x=xseq,
                   func=Prob17, guess=0,atol=1e-8)


plot(twp1[,1],twp1[,2],type="l",main="test problem 17",ylim=c(-1,1),
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=2)
curve(x/sqrt(xi+x^2),-0.1,0.1,type="p",add=TRUE)


#===============================================================================
### Test problem 18
#===============================================================================

Prob18 <- function(x, y, pars) {
  list(c( y[2], -1/xi*y[2]) )
}

xseq<-seq(0,1,by=0.01)
xi <-0.2
shoot  <- bvpshoot(yini=c(1,NA),yend=c(exp(-1/xi),NA),x=xseq,
           func=Prob18, guess=0,atol=1e-10)

xi <- 0.1
twp  <- bvptwp(yini=c(1,NA),yend=c(exp(-1/xi),NA),x=xseq,
           func=Prob18, guess=0,atol=1e-10)

xi <- 0.01
mod  <- bvpcolmod(yini=c(1,NA),yend=c(exp(-1/xi),NA),x=xseq,
           func=Prob18, guess=0,atol=1e-10)


plot(shoot[,1],shoot[,2],type="l",main="test problem 18",ylim=c(0,1),
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)
curve(exp(-x/xi),0,1,type="p",add=TRUE)


#===============================================================================
### Test problem 19
#===============================================================================

Prob19 <- function(t, y, pars, ksi) {
  pit = pi*t
  list(c(y[2],(pi/2*sin(pit/2)*exp(2*y[1])-exp(y[1])*y[2])/ksi))
}

xi <-0.05
shoot  <- bvpshoot(yini=c(0,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
        func=Prob19, guess=0,ksi=xi)

xi <- 0.03
twp  <- bvptwp(yini=c(0,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob19, guess=0,ksi=xi,atol=1e-15)

xi <- 0.005
print(system.time(

mod  <- bvpcolmod(yini=c(0,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob19, guess=0,ksi=xi, atol=1e-10)
))

# same, using continuation
print(system.time(
solcol3  <- bvpcolmod(yini=c(0,NA),yend=c(0,NA),x=seq(0,1,by=0.01),
           func=Prob19, guess=0,eps=xi)
))


plot(shoot[,1],shoot[,2],type="l",main="test problem 19", ylim=c(-0.7,0),
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)


#===============================================================================
### Test problem 20
#===============================================================================

Prob20 <- function(x, y, pars, xi) {
  list(c( y[2] , 1/xi *(1-y[2]^2) ))
}

xi <-0.5
ini <- c(1+xi * log(cosh(0.745/xi)),NA)
end <- c(1+xi * log(cosh(0.255/xi)),NA)
mod1  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob20, guess=0, eps=xi)

xi <-0.3
ini <- c(1+xi * log(cosh(0.745/xi)),NA)
end <- c(1+xi * log(cosh(0.255/xi)),NA)
mod2  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob20, guess=0, eps=xi)

xi <-0.01
ini <- c(1+xi * log(cosh(0.745/xi)),NA)
end <- c(1+xi * log(cosh(0.255/xi)),NA)
mod3  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob20, guess=0, eps=xi)


plot(mod1[,1],mod1[,2],type="l",main="test problem 20", ylim=c(1,1.8),
  xlab="x",ylab="y")
lines(mod2[,1],mod2[,2],lty=2)
lines(mod3[,1],mod3[,2],lty=3)
curve(1+xi * log(cosh((x-0.745)/xi)),0,1,add=TRUE,type="p")


#===============================================================================
### Test problem 21
#===============================================================================

Prob21 <- function(x, y, pars, xi) {
  list(c( y[2], 1/xi *(y[1]+y[1]^2-exp(-2*x/sqrt(xi))) ))
}

ini <- c(1,NA)

xi <-0.2
end <- c(exp(-1/sqrt(xi)),NA)
shoot  <- bvpshoot(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob21, guess=0, xi=xi)

xi <-0.1
end <- c(exp(-1/sqrt(xi)),NA)
twp  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob21, guess=0, xi=xi)

xi <-0.01
end <- c(exp(-1/sqrt(xi)),NA)
mod  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob21, guess=0, eps=xi)


plot(shoot[,1],shoot[,2],type="l",main="test problem 21", ylim=c(0,1),
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)
curve(exp(-x/sqrt(xi)),0,1,add=TRUE,type="p")


#===============================================================================
### Test problem 22
#===============================================================================

Prob22 <- function(t, y, pars, xi) {
  list(c( y[2], -1/xi *(y[2]+y[1]^2) ))
}

ini <- c(0,NA)
end <- c(1/2,NA)

xi <-0.1
shoot  <- bvpshoot(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob22, guess=0, xi=xi)

xi <-0.05
twp  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob22, guess=0, xi=xi)

xi <- 0.01
mod  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob22, guess=0,xi=xi)


plot(shoot[,1],shoot[,2],type="l",main="test problem 22", ylim=c(0,1),
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)


#===============================================================================
### Test problem 23
#===============================================================================

Prob23 <- function(t, y, pars, xi) {
  list(c( y[2], sinh(y[1]/xi)/xi) )
}

ini <- c(0,NA)
end <- c(1,NA)

xi <- 1/5
twp1  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob23, guess=c(0),xi=xi)

xi <- 1/7
twp  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob23, guess=0,xi=xi)

xi <- 1/9
mod  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob23, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 23",
  xlab="x",ylab="y")
lines(twp[,1],twp[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)


#===============================================================================
### Test problem 24
#===============================================================================

Prob24 <- function(t, y, pars, xi) {
  A <- 1+t*t
  AA <- 2*t
  ga <- 1.4
  list(c(y[2],
       (((1+ga)/2 -xi*AA)*y[1]*y[2]-y[2]/y[1]-
       (AA/A)*(1-(ga-1)*y[1]^2/2))/(xi*A*y[1])  ))
}

ini <- c(0.9129,NA)
end <- c(0.375,NA)

xi <-0.03
mod1  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess=0 ,eps=xi)

xi <-0.01
mod2  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess=0 ,eps=xi)

xi <-0.0015
mod3  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24, guess= c(0) ,eps=xi, verbose=TRUE)


plot(mod1[,1],mod1[,2],type="l",main="test problem 24", ylim=c(0.2,1.4),
  xlab="x",ylab="y")
lines(mod2[,1],mod2[,2],lty=2)
lines(mod3[,1],mod3[,2],lty=3)


#===============================================================================
### Test problem 25
#===============================================================================

Prob25 <- function(t, y, pars, xi) {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(-1/3,NA)
end <- c(1/3,NA)
xi  <-0.1
twp1  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob25, guess=0, xi=xi)

xi <-0.01
twp2  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob25, guess=0,xi=xi)

xi <- 0.001
mod  <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob25, guess=0,eps=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 25",
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(mod[,1],mod[,2],lty=3)


#===============================================================================
### Test problem 26
#===============================================================================

Prob26 <- function(t, y, pars, xi)  {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(1,NA)
end <- c(-1/3,NA)

xi <-0.1
twp1  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob26, guess=0,xi=xi)

xi <-0.02
twp2  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob26, guess=0,xi=xi)

xi <-0.005
twp3  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob26, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 26",
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=3)


#===============================================================================
### Test problem 27
#===============================================================================

Prob27 <- function(t, y, pars, xi) {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(1,NA)
end <- c(1/3,NA)

xi   <-0.1
twp1 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob27, guess=0,xi=xi)

xi   <-0.02
twp2 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob27, guess=0,xi=xi)

xi <-0.005
twp3 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob27, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 27", ylim=c(0,1),
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=3)


#===============================================================================
### Test problem 28
#===============================================================================

Prob28 <- function(t, y, pars, xi)  {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(1,NA)
end <- c(3/2,NA)

xi   <-0.1
twp1 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob28, guess=0,xi=xi)

xi   <-0.02
twp2 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob28, guess=0,xi=xi)

xi <-0.005
twp3 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob28, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 28", ylim=c(0.4,1.5),
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=3)


#===============================================================================
### Test problem 29
#===============================================================================

Prob29 <- function(t, y, pars, xi)  {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(0,NA)
end <- c(3/2,NA)

xi   <-0.1
twp1 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob29, guess=0,xi=xi)

xi   <-0.02
twp2 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob29, guess=0,xi=xi)

xi <-0.005
twp3 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob29, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 29",
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=3)


#===============================================================================
### Test problem 30
#===============================================================================

Prob30 <- function(t, y, pars, xi)  {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(-7/6,NA)
end <- c(3/2,NA)

xi   <-0.1
twp1 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob30, guess=0,xi=xi)

xi   <-0.02
twp2 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob30, guess=0,xi=xi)

# xi <-0.005     # this small value may give problems...
xi <-0.01
twp3 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob30, guess=0,xi=xi)


plot(twp1[,1],twp1[,2],type="l",main="test problem 30",
  xlab="x",ylab="y")
lines(twp2[,1],twp2[,2],lty=2)
lines(twp3[,1],twp3[,2],lty=3)


#===============================================================================
### Test problem 31
#===============================================================================

Prob31 <- function(t, Y, pars)  {
  with (as.list(Y), {
    dy    <- sin(Tet)
    dTet  <- M
    dM    <- -Q/xi
    T <- 1/cos (Tet) +xi*Q*tan(Tet)
    dQ    <- 1/xi*((y-1)*cos(Tet)-M*T)
    list(c( dy, dTet, dM, dQ))
  })
}

ini <- c(y=0,Tet=NA,M=0,Q=NA)
end <- c(y=0,Tet=NA,M=0,Q=NA)

xi <-0.1
twp  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob31, guess=c(0,0),atol=1e-10)

xi <- 0.05
mod <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob31, guess=c(0,0),atol=1e-10)

xi <- 0.01
mod2 <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob31, guess=c(0,0),atol=1e-10)


plot(twp[,1],twp[,4],type="l",main="test problem 31",
  xlab="x",ylab="y3")
lines(mod[,1],mod[,4],lty=2)
lines(mod2[,1],mod2[,4],lty=3)


#===============================================================================
### Test problem 32
#===============================================================================

Prob32 <- function(t, y, pars, xi) {
  list(c( y[2], y[3], y[4], 1/xi*(y[2]*y[3]-y[1]*y[4])))

}
ini <- c(0,0,NA,NA)
end <- c(1,0,NA,NA)

xi <-0.01
twp <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob32, guess=c(0,0), xi=xi)

xi <-0.002
twp2 <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob32, guess=c(0,0), xi=xi)

xi <-0.0001
mod <- bvpcolmod(yini=ini,yend=end,x=seq(0,1,by=0.01),
           func=Prob32, guess=c(0,0), xi=xi)


plot(twp[,1],twp[,3],type="l",main="test problem 32",
  xlab="x",ylab="y'")
lines(twp2[,1],twp2[,3],lty=2)
lines(mod[,1],mod[,3],lty=3)

