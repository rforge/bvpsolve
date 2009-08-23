## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## Can be solved with bvpshoot
## =============================================================================

Reactor<-function(x,y,parms)
{
  list(c(y[2],Pe*(y[2]+R*(y[1]^n))))
}

Pe <- 1
R  <- 2
n  <- 2

yini <- function (x,parms) c(x,Pe*(x-1))

sol<-bvpshoot(func=Reactor, yend=c(y=NA,dy=0), yini=yini,
              x=seq(0,1,by=0.01),extra=1)
plot(sol)

