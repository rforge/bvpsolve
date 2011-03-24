## =============================================================================
##
## The medical AKZO problem 
##                     
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 400
##
## =============================================================================
#
# use ReacTran to solve this..
require(ReacTran)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# The grid properties
Len <- 10 
dx  <- 0.05
x   <- seq (from = dx/2, to = Len, by = dx)
N   <- length(x)

# initial condition of state variables
yini <- c(rep(0,N),rep(1,N))

# model parameters
k   <- 100      # reaction rate
D   <- 1        # diffusion coefficient

# derivative function
medAKZO <- function(t,y,parms) {

    # two state variable vectors
    u <- y[1:N]
    v <- y[-(1:N)]
    
    # boundary condition for u  
    phi <- if (t <= 5) 2 else 0

    # rate of change: only u is transported
    du <- tran.1D(C=u, D=D, C.up=phi, dx=dx)$dC - k*u*v
    dv <- - k*u*v

    # return the rate of changes, concatenated, as a list
    list(c(du,dv))
}

# -------------------------------------------------------
# solve the model
# -------------------------------------------------------

# time sequence
times <- seq(0,20, by = 0.01)

print(system.time(
  out  <- ode.1D(func=medAKZO, times=times, parms=NULL, y=yini, nspec=2)
))

# -------------------------------------------------------
# plot
# -------------------------------------------------------

# remove time column
Out <- out[,-1]

par (mfrow=c(2,2))  # 2 rows, 2 columns
image(Out[,1:N],x=times, y = x, col=femmecol(), main="u", zlim=c(0,2))
image(Out[,-(1:N)],x=times, y = x,col=femmecol(), main="v", zlim=c(0,2))

emptyplot()
colorlegend(zlim=c(0,2), digit=2, posx=c(0.5,0.53))

plot(out,which=1,type="l",lwd=2, mfrow=NULL, col="darkblue", main="u")
for ( i in 1:10) lines(times,out[,10*i],col="lightblue")
mtext(outer=TRUE,side=3,"Chemical AKZO",line=-1.5, cex=1.5)
