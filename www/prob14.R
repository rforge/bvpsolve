library(bvpSolve)
library(deSolve)
Prob14 <- function(x, y, pars)  {
  list(1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x)))
}


xi    <- 0.01
mod1 <- bvpshoot(yini = c(0, NA), yend = c(0, NA), 
                 order = 2, x = seq(-1, 1, by = 0.01), func = Prob14, guess = 0)
xi   <- 0.0025
mod2 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
               order = 2, x = seq(-1, 1, by = 0.01), func = Prob14)
xi    <- 0.0001
mod3 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
               order = 2, x = seq(-1, 1, by = 0.01), func = Prob14)

png("prob14.png", width=600, height=200)
par(lwd=2)
par(mar=c(3,3,1.5,0)+.1)
plot(mod1, mod2, mod3, which = 1, lty = c(1, 2,3 ), main = "test problem 14")
mtext("x", side=1, line=2)
mtext("y", side=2, line=2)
#curve(cos(pi*x)+exp((x-1)/sqrt(xi))+exp(-(x+1)/sqrt(xi)),
#      -1, 1, type = "p", add = TRUE)
dev.off()
