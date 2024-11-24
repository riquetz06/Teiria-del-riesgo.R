#TALLER DE R CON APLICACIONES ACTUARIALES
#ACT. Y M. EN F. RIQUET ZEQUEIRA FERNÁNDEZ

IMPORTAR BASES DE DATOS DE EXCEL      

install.packages("RODBC")
library(RODBC)
data= odbcConnectExcel(file.choose())
sqlTables(data)
mydat=sqlFetch(data, "nombre hoja") 
odbcClose(data)
mydat

COMANDOS BÁSICOS PÁRA ANÁLISIS ESTADÍSTICO

ejemplo<-scan()
max(ejemplo)
pie(ejemplo)
hist(ejemplo)
normal<-rnorm(250)
hist(normal)
hist(normal, breaks=50, freq=F)
> hist(normal, breaks=50, freq=F, main= "HISTOGRAMA DIST NORMAL",
+ xlab="numeros del eje x", ylab="números del eje y",
+ xlim=c(-3,3), ylim=c(0,0.6)
+ , col="51")
curve(dnorm, add=T)


APLICACIONES A LA TEORÍA DEL RIESGO (UTILIZANDO EL PAQUETE ACTUAR EN R)

####################################
#  F\'ormula de De Pril en R v1.1  #
####################################
I <- 5 # Montos de reclamaciones
J <- 3 # \'Indice m\'aximo para tasas de muerte
R <- 20 # Valor m\'aximo para r en g(r)
#
n <- array(1:15, dim=c(5,3))
n[1,1]<-1
n[2,1]<-3
n[3,1]<-5
n[4,1]<-2
n[5,1]<-2
n[1,2]<-3
n[2,2]<-5
n[3,2]<-3
n[4,2]<-2
n[5,2]<-3
n[1,3]<-1
n[2,3]<-4
n[3,3]<-4
n[4,3]<-6
n[5,3]<-4
#
q <- array(1:3, dim=c(3))
q[1]<-0.02
q[2]<-0.025
q[3]<-0.03
#...............................
# Funci\'on h(i,k)
#...............................
h <- function(i,k) {
aux <- 0
for (j in 1:J) {
aux <- aux+n[i,j]*(q[j]/(1-q[j]))^k
}
aux <- i*((-1)^(k-1))*aux
return(aux)
}
#...............................
# C\'alculo de la densidad de S
#...............................
gc <- array(1:R, dim=c(R))
gc0 <- g(0)
#
g <- function(r) {
if (r==0) {
aux <- 1
for (i in 1:I) {
    for (j in 1:J) {
    aux <- aux*((1-q[j])^n[i,j])
    }
}
return(aux)
}
else
{
aux <- 0
for (i in 1:min(r,I)) {
    for (k in 1:floor(r/i)) {
        if (r-i*k==0) { aux <- aux + gc0*h(i,k) }
        else {aux <- aux + gc[r-i*k]*h(i,k)}
        }
}
aux <- aux/r
gc[r] <- aux
return(aux)
}
}
#...............................
# Asignaci\'on en el arreglo "gc" y graficaci\'on.
#...............................
for (i in 1:R) {
gc[i] <- g(i)
}
# Nota: Se omite en la gr\'afica el valor de la densidad en cero "gc0".
barplot(gc,main="Funci?n de densidad de S",xlab="r", ylab="g(r)")
#################################################################
# Fin de c\'odigo
#################################################################
####################################
#  F\'ormula de Panjer en R v1.0   #
#  [Caso Poisson]                  #
####################################
#
R <- 30 # Valor m\'aximo para r en g(r)
#
#...............................
# c\'alculo de p_k=P(N=k) (Caso Poisson)
#...............................
a <- 0
b <- 3  #lambda
p0 <- 2.7172^{-b}
p <- array(1:R, dim=c(R))
p[1] <- (a+b)*p0
for (k in 2:R) {
	p[k] <- (a+b/k)*p[k-1]
		}
#...............................
# c\'alculo de f_r=P(Y=r), r>=1
#...............................
#
f <- array(1:R, dim=c(R))
f[1] <- 0.1
f[2] <- 0.2
f[3] <- 0.3
f[4] <- 0.4
for (i in 5:R) { f[i] <- 0 }
#................................
# C\'alculo de la densidad de S
#................................
g0 <- p0
g <- array(1:R, dim=c(R))
g[1] <- (a+b)*f[1]*g0
for (r in 2: R) {
	aux <- 0
	for (i in 1:{r-1}) {
	    aux <- aux + (a+b*i/r)*f[i]*g[r-i]
	               }
	aux <- aux + (a+b)*f[r]*g0
	g[r] <- aux
		}
#...............................
# Graficaci\'on
#...............................
# Nota: Se omite en la gr\'afica el valor de la densidad en cero "g0".
barplot(g,main="Funci?n de densidad de S",xlab="r", ylab="g(r)")
#
#################################################################
# Fin de c\'odigo
#################################################################


#####################################
# Distribución del monto de reclamo #
#####################################

## Basic example: no reinsurance, exponential claim severity and wait
## times, premium rate computed with expected value principle and
## safety loading of 20%.
adjCoef(mgfexp, premium = 1.2, upper = 1)
## Same thing, giving function h.
h <- function(x) 1/((1 - x) * (1 + 1.2 * x))
adjCoef(h = h, upper = 1)
## Example 11.4 of Klugman et al. (2008)
mgfx <- function(x) 0.6 * exp(x) + 0.4 * exp(2 * x)
adjCoef(mgfx(x), mgfexp(x, 4), prem = 7, upper = 0.3182)
## Proportional reinsurance, same assumptions as above, reinsurer’s
## safety loading of 30%.
mgfx <- function(x, y) mgfexp(x * y)
p <- function(x) 1.3 * x - 0.1
h <- function(x, a) 1/((1 - a * x) * (1 + x * p(a)))
R1 <- adjCoef(mgfx, premium = p, upper = 1, reins = "proportional",
from = 0, to = 1, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "p",
from = 0, to = 1, n = 101)
R1(seq(0, 1, length = 10)) # evaluation for various retention rates
R2(seq(0, 1, length = 10)) # same
plot(R1) # graphical representation
plot(R2, col = "green", add = TRUE) # smoother function
## Excess-of-loss reinsurance
p <- function(x) 1.3 * levgamma(x, 2, 2) - 0.1
mgfx <- function(x, l)
mgfgamma(x, 2, 2) * pgamma(l, 2, 2 - x) +
exp(x * l) * pgamma(l, 2, 2, lower = FALSE)
h <- function(x, l) mgfx(x, l) * mgfexp(-x * p(l))
R1 <- adjCoef(mgfx, upper = 1, premium = p, reins = "excess-of-loss",
from = 0, to = 10, n = 11)
R2 <- adjCoef(h = h, upper = 1, reins = "e",
from = 0, to = 10, n = 101)
plot(R1)
plot(R2, col = "green", add = TRUE)

















#################
# Convoluciones #
#################

## Convolution method (example 9.5 of Klugman et al. (2008))
fx <- c(0, 0.15, 0.2, 0.25, 0.125, 0.075,
0.05, 0.05, 0.05, 0.025, 0.025)
pn <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.15, 0.06, 0.03, 0.01)
Fs <- aggregateDist("convolution", model.freq = pn,
model.sev = fx, x.scale = 25)
summary(Fs)
c(Fs(0), diff(Fs(25 * 0:21))) # probability mass function
plot(Fs)
## Recursive method
Fs <- aggregateDist("recursive", model.freq = "poisson",
model.sev = fx, lambda = 3, x.scale = 25)
plot(Fs)
Fs(knots(Fs)) # cdf evaluated at its knots
diff(Fs) # probability mass function
## Recursive method (high frequency)
## Not run: Fs <- aggregateDist("recursive", model.freq = "poisson",
model.sev = fx, lambda = 1000)
## End(Not run)
Fs <- aggregateDist("recursive", model.freq = "poisson",
model.sev = fx, lambda = 250, convolve = 2, maxit = 1500)
plot(Fs)
## Normal Power approximation
Fs <- aggregateDist("npower", moments = c(200, 200, 0.5))
Fs(210)
## Simulation method
model.freq <- expression(data = rpois(3))
model.sev <- expression(data = rgamma(100, 2))
Fs <- aggregateDist("simulation", nb.simul = 1000,
model.freq, model.sev)
mean(Fs)
plot(Fs)
## Evaluation of ruin probabilities using Beekman’s formula with
## Exponential(1) claim severity, Poisson(1) frequency and premium rate
## c = 1.2.
fx <- discretize(pexp(x, 1), from = 0, to = 100, method = "lower")
phi0 <- 0.2/1.2
Fs <- aggregateDist(method = "recursive", model.freq = "geometric",
model.sev = fx, prob = phi0)
1 - Fs(400) # approximate ruin probability
u <- 0:100
plot(u, 1 - Fs(u), type = "l", main = "Ruin probability")




####################################
#  Distribución de la pérdida      #
#  [Caso Gamma]                    #
####################################
> n <- 100000
> sample <- rgamma(n, rate=0.5, shape=5)
> library(MASS)
> truehist(sample, main="Histogram of sample")
> truehist(log(sample), main="Histogram of log(sample)")
> fitdistr(sample, "Weibull")
> # Fit distributions
> fit.w <- fitdistr(sample, "Weibull")
> fit.g <- fitdistr(sample, "gamma")
> # Caclulate points for density curves
> x <- seq(from=0, to=40, by=1)
> y.w <- dweibull(x, shape=fit.w$estimate[1],
+ scale=fit.w$estimate[2])
> y.g <- dgamma(x, shape=fit.g$estimate[1],
+ rate=fit.g$estimate[2])
> # Draw histogram and density curves
> truehist(sample, main="Histogram of sample")
> lines(x, y.w)
> lines(x, y.g, col="red")
> legend(30, 0.1, legend=c("Gamma distn", "Weibull distn"),
+ lwd=1, col=c("red", "black"))
> # Draw cumulative density functions
> plot(ecdf(sample), cex=0)
> curve(pweibull(x, shape=fit.w$estimate[1],
+ scale=fit.w$estimate[2]), add=T, col=2)
> curve(pgamma(x, shape=fit.g$estimate[1],
+ rate=fit.g$estimate[2]), add=T, col=3)
> legend(30, 0.6, lwd=1, col=1:3,
+ legend=c("Empirical CDF", "Weibull CDF", "Gamma CDF"))
RESERVAS PARA SEGUROS GENERALES

####################################
#            Chain Ladder          #
####################################

Cargar paquete Chain Ladder:
> library(ChainLadder)
> require(ChainLadder)
> data(package="ChainLadder")
> ## Sample triangle

Ejemplo 1:

> RAA
> plot(RAA)
> plot(RAA, lattice=TRUE)
> ?plot.triangle
> raa.inc <- cum2incr(RAA)
## Show first origin period and its incremental development
> raa.inc[1,]
> raa.cum <- incr2cum(raa.inc)
> ## Show first origin period and its cumulative development
> raa.cum[1,]


Ejemplo 2: (Importar triángulo de siniestros ocurridos desde Excel).
> file.choose() #indica la ruta de nuestro archivo .csv#
> myCSVfile <- "C:\\Users\\Riquet\\Downloads\\Libro1.csv" 
> tri <- read.csv(file=myCSVfile, header = FALSE)
> library(ChainLadder)
> tri <- as.triangle(as.matrix(tri))

Ejemplo 3: (DEMOS).
> install.packages("ROBDC")
> library(RODBC)
> demo(DatabaseExamples)




