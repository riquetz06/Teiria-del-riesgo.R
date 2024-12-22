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

COMANDOS BÁSICOS PÁRA ANALISIS ESTADÍSTICO

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

#Loss Distributions for Actuaries Course
lnClaims < - rlnorm(10000, meanlog = 2, sdlog = 0.5)
expClaims < - rexp(10000, rate = 1)
hist(lnClaims, main = "Histogtram of 10,000 Simulated Claims Values from  the Log Normal Distribution")
hist(expClaims, main = "Histogtram of 10,000 Simulated Claims Values from  the Log Exponential Distribution")
install.packages("MASS")
library(MASS)
fitdistr(lnClaims,"lognormal")
fitdistr(lnClaims,"exponential")

#Data Visualization
set.seed(123)
lnClaims<-rlnorm(10000, meanlog = 2, sdlog = 0.5)
lnfit<-fitdistr(lnClaims, "lognormal")
lnfitmean<-lnfit$estimate["meanlog"]
lnfitsd<-lnfit$estimate["sdlog"]
hist(lnClaims, main = "Lognormal Loss Distribution mean of 2 and sd 0.5", freq = FALSE)
curve(dlnorm(x,lnfitmean,lnfitsd), add = TRUE, col = "red")

lnClaims2<-rlnorm(10000, meanlog = 2, sdlog = 0.5)
lnfit2<-fitdistr(lnClaims2, "lognormal")
lnfitmean2<-lnfit2$estimate["meanlog"]
lnfitsd2<-lnfit2$estimate["sdlog"]
hist(lnClaims2, main = "Lognormal Loss Distribution mean of 2 and sd 0.5", freq = FALSE)
curve(dlnorm(x,lnfitmean2,lnfitsd2), add = TRUE, col = "red")  
lines(density(lnClaims2), col="blue" )
n=length(lnClaims2)
qqplot(qlnorm(ppoints(n), lnfitmean2,lnfitsd2), lnClaims2, xlab="Theorical quantiles", ylab="Sample quantiles",
       main ="Comparision of actual vs Fitted")
abline(0,1)

###############
# Reinsurance #
###############
ClaimsX<-rlnorm(10000, meanlog = 2, sdlog = 0.5) #Synthetic data
RetentionM = 15 # Retention Limit M = 15
ClaimsY <- pmin(ClaimsX,RetentionM)
library(ggplot2)
datY <- data.frame(loss = ClaimsY)
head(datY)
ggplot(datY, aes(x=loss))+geom_histogram(aes(y=..density..),binswidth=0.5, bins=50, colour ="black", fill="white")+
  geom_density(alpha=0.2,fill="#00FFFF")+geom_vline(aes(xintercept=mean(loss)),colour="red", linetype="dashed",size=1)+
  xlim(0,20)+ylim(0, 0.3)

