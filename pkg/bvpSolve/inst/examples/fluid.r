## =============================================================================
## fluid injection problem
## f'''- R[(f')^2 -f*f''] + RA = 0
## h'' + R*f*h' + 1 = 0 
## O'' + P*f*O' = 0
## A is unknown
##
## rewritten as:
## df=f1                   #f
## df1=f2                  #f1=f'
## df2=R(f1^2-f*f2)-A      #f2=f''
## dh=h1
## dh1= -Rfh1-1
## dO=O1
## dO1 = O2
## dO2 = -P*f*O1
## dA = 0
## =============================================================================
require(bvpSolve)

fluid<-function(t,y,pars,R)
{
 P    <- 0.7*R
 with(as.list(y),
 {
  df = f1                   #f'
  df1= f2                   #f''
  df2= R*(f1^2-f*f2)-A      #f'''
  dh = h1
  dh1= -R*f*h1-1
  dO = O1
  dO1= -P*f*O1              # the constant to be estimated
  dA = 0
  return(list(c(df,df1,df2,dh,dh1,dO,dO1,dA)))
 })
}


R      <- 100
P      <- 0.7*R
times  <- seq(0,1,by=0.01)
yini   <- c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA)
yend   <- c(1,0,NA,0,NA,1,NA,NA)

Soltwp1 <- bvptwp(func=fluid, x=times, parms=NULL, R=R,
                 yini=c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA),
                    yend=c(1,0,NA,0,NA,1,NA,NA))
plot(Soltwp1, which="f1", type="l", lwd=2)

R      <- 1000
P      <- 0.7*R
times  <- seq(0,1,by=0.01)
yini   <- c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA)
yend   <- c(1,0,NA,0,NA,1,NA,NA)

Soltwp2 <- bvptwp(func=fluid, x=times, parms=NULL, R=R,
                 yini=c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA),
                    yend=c(1,0,NA,0,NA,1,NA,NA))
lines(Soltwp2[,"x"], Soltwp2[,"f1"], col="red", lwd =2)

R      <- 10000
P      <- 0.7*R
times  <- seq(0,1,by=0.01)
yini   <- c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA)
yend   <- c(1,0,NA,0,NA,1,NA,NA)

print(system.time(
Soltwp3 <- bvptwp(func=fluid, x=times, parms=NULL, R=R,
                 yini=c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA),
                    yend=c(1,0,NA,0,NA,1,NA,NA))
))
                    
lines(Soltwp3[,"x"], Soltwp3[,"f1"], col="darkblue", lwd =2)

# For more extreme value we use the previous solution as initial guess..
R      <- 20000
P      <- 0.7*R
times  <- seq(0,1,by=0.01)
yini   <- c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA)
yend   <- c(1,0,NA,0,NA,1,NA,NA)

print(system.time(
Soltwp4 <- bvptwp(func=fluid, x=times, parms=NULL, R=R,
                 yini=c(f=0,f1=0,f2=NA,h=0,h1=NA,O=0,O1=NA,A=NA),
                 xguess = Soltwp3[,1], yguess=t(Soltwp3[,-1]),
                    yend=c(1,0,NA,0,NA,1,NA,NA))
))
                    
lines(Soltwp4[,"x"], Soltwp4[,"f1"], col="darkgreen", lwd =2)
