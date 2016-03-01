library(fBasics)


Vasicek = function(r_0, a, b, sigma, dt) { r_0*exp(-a*dt) + b*(1-exp(-a*dt)) + sigma*(sqrt( (1-exp(-2*a*dt))/(2*a) ))*rnorm(1, 0, 1)}

####################### Vasicek Parameters #################################
## r_0 = initial interest rate
## a = speed of reversion
## b = long term mean
## sigma = volatility
## T = time to maturity
## N = number of timesteps in path

r_0 = .07
a = .5
b = .045
sigma = .01
T = 1
N = 100


####################### Monte Carlo Parameters ###########################
## MC_N = number of paths

MC_N = 500





dt = T/N
MC = rep(0, MC_N)
MC[1] = r_0
Paths = matrix(0, MC_N, N)
minVal = r_0
maxVal = r_0

## Alter screen layout for graphs
nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), respect=TRUE)
#layout.show(nf)


par(mar=c(3,3,2,2)) 
for ( j in 1:MC_N ) {
  
  Paths[j,1] = r_0

  for ( i in 2:N )
  {
    Paths[j,i] = Vasicek(Paths[j,i-1], a, b, sigma, dt)
    if ( Paths[j,i] < minVal ) { minVal = Paths[j,i] }
    if ( Paths[j,i] > maxVal ) { maxVal = Paths[j,i] }
  }
  MC[j] = Paths[j,i]
}


for ( j in 1:MC_N )
{
  if ( j == 1 ) { plot(1:N, Paths[j,],type='l', ylim=c(minVal, maxVal)) }
  else { 
  	if ( j %% 2 == 0 ) { lines(1:N, Paths[j,], type='l', col=1)  }
  	else { lines(1:N, Paths[j,], type='l', col=2)  }
  }
}
title("Simulated Interest Rate Paths")


yhist = hist(MC, plot=FALSE,)
par(mar=c(3,0,2,2)) 
barplot(yhist$counts, axes=FALSE, xlim=c(0, max(yhist$counts)), space=0, horiz=TRUE) 
hist(MC, freq=FALSE, col="yellow", main="Histogram and Density of Monte Carlo Simulation")
lines(density(MC), lwd=3, lty=3)
abline(v=mean(MC), lwd=8)
legend("topright", "Density Function", lwd=4, lty=3)

mean(MC)
