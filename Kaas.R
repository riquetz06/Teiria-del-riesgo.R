###
### Chapter 1 
###

set.seed(2525); n <- 1e6; X <- rnorm(n)
rbind(c(mean(X[X>1]), mean(X*(X>1)), mean(pmax(X-1,0))),
      c(dnorm(1)/pnorm(-1), dnorm(1), dnorm(1)-pnorm(-1)))
## Result:
##          [,1]      [,2]       [,3]
## [1,] 1.524587 0.2416943 0.08316332
## [2,] 1.525135 0.2419707 0.08331547



###
### Chapter 2
###

x <- 3.5; mu <- 1; sig <- 1; gam <- 1; z <- (x-mu)/sig
1-pbinom(x, 1000, 0.001)                      ##  0.01892683
1-ppois(x,1)                                  ##  0.01898816
1-pnorm(z)                                    ##  0.00620967
1-pnorm(sqrt(9/gam^2 + 6*z/gam + 1) - 3/gam)  ##  0.02275013
1-pgamma(x-(mu-2*sig/gam), 4/gam^2, 2/gam/sig)##  0.02122649



###
### Chapter 3
###

SparseVec <- function (freq)
{if (any(freq<0)) stop("negative frequency")
 M <- length(freq)
 mu <- sum((1:M)*freq); sigma2 <- sum((1:M)^2*freq)
 ##mean and variance of the compound r.v.; see (3.4)
 MM <- ceiling(mu + 10 * sqrt(sigma2)) + 6
 fs <- dpois(0:(MM-1), freq[1])  ##density of S_1 = 1*N_1
 for (j in 2:M)
 {MMM <- trunc((MM-1)/j)
  fj <- rep(0, MM)  ##construct the density of j*N_j
  fj[(0:MMM)*j+1] <- dpois(0:MMM, freq[j])
  fs <- convolve(fs, rev(fj), type="o") }
 ##fs is the density of S_j = 1*N_1 + ... + j*N_j, j=2..M
 return(fs)  }
f <- SparseVec(c(1,2,1)); f[1:7] * exp(4)



Panjer.Poisson <- function (p, lambda) 
{ if (sum(p)>1||any(p<0)) stop("p parameter not a density")
  if (lambda * sum(p) > 727) stop("Underflow")
  cumul <- f <- exp(-lambda * sum(p))
  r <- length(p)
  s <- 0
  repeat
  { s <- s+1
    m <- min(s, r)
    last <- lambda / s * sum(1:m * head(p,m) * rev(tail(f,m)))
    f <- c(f,last)
    cumul <- cumul + last
    if (cumul > 0.99999999) break  }
    return(f)  }
Panjer.Poisson(c(0.25,0.5,0.25), 4) * exp(4)



n <- 64; p <- rep(0, n); p[2:3] <- 0.5; lab <- 1
f <- Re(fft(exp(lab*(fft(p)-1)), inverse=TRUE))/n



set.seed(1); n <- 2000; r <- 2; p <- 0.5
hh <- rpois(n, rgamma(n,r,p/(1-p)))
n.j <- tabulate(1+hh); j <- 0:max(hh)
rbind(n.j, round(dnbinom(j,r,p)*n))



y.bar <- sum(j*n.j/n); y2.bar <- sum(j^2*n.j/n)
p0 <- y.bar/(y2.bar-y.bar^2); r0 <- p0 * y.bar/(1-p0)



g <- function (r) {sum(dnbinom(j,r,r/(r+y.bar),log=T)*n.j)}
r <- optimize(g, c(r0/2, 2*r0), maximum=T, tol=1e-12)$maximum
p <- r/(r+y.bar)



h <- function (x) {-sum(dnbinom(j,x[1],x[2],log=T)*n.j)}
optim(c(r0,p0), h, control=list(reltol=1e-14))



set.seed(2525); y <- rgamma(2000, shape=5, rate=1)
aux <- log(mean(y)) - mean(log(y))
f <- function(x) log(x) - digamma(x) - aux
alpha <- uniroot(f, c(1e-8,1e8))$root  ##  5.049 
beta <- alpha/mean(y)                  ##  1.024


require(statmod)
library(tweedie); set.seed(2525)
y <- rinvgauss(2000, mean=5, shape=3)
alpha <- 1/(mean(y)*mean(1/y)-1); beta <- alpha/mean(y)
alpha/beta; alpha^2/beta  ## 4.9448 3.0446


set.seed(1); n <- 2000; q <- 1.5; alpha <- 1; beta <- 2
gam <- (beta-alpha)/beta * q
y <- rbinom(n,1,gam) * rexp(n)/alpha + rexp(n)/beta



f <- function (y, q, alpha, beta){
  q * alpha * exp(-alpha*y) + (1-q) * beta * exp(-beta*y)}
h <- function (x) {-sum(log(f(y, x[1], x[2], x[3])))}
optim(c(0.8, 0.9, 1.8), h)



set.seed(2525); x0 <- 100; alpha <- 2; n <- 2000
y <- x0*exp(rexp(n)/alpha)
x0.hat <- min(y); alpha.hat <- 1/mean(log(y/x0.hat))



y <- rep(0,64); y[2:7] <- 1/6; Re(fft(fft(y)^10,TRUE))/64



###
### Chapter 4 
###

lab <- 1; EX <- 1; theta <- 0.3; u <- 7.5; alpha <- 2
n <- 400; nSim <- 10000; set.seed(2)
c. <- (1+theta)*lab*EX
N <- rep(Inf, nSim)
for (k in 1:nSim){
  Wi <- rexp(n)/lab; Ti <- cumsum(Wi)
  Xi <- rgamma(n, shape=alpha)/alpha ## severity has mean EX=1
  Si <- cumsum(Xi); Ui <- u + Ti*c. - Si
  ruin <- !all(Ui>=0)
  if (ruin) N[k] <- min(which(Ui<0))}
N <- N[N<Inf]; length(N); mean(N); sd(N); max(N)
##  745  30.78792  24.71228  255



f <- function (r) exp(r)/2 + exp(2*r)/2 - 1 - 1.8*r
R <- uniroot(f, lower=0.00001, upper=1)$root ## 0.2105433



psi <- function (u, theta, x, y, px)
{ if (px<0||px>1||theta<=0||x<=0||y<=0||u<0) stop("bad params")
  mu <- x*px + y*(1-px)
  ss <- 0
  for (k in 0:(u/x))
  { n <- 0:((u-k*x)/y)
    h <- 1/mu/(1+theta) * (u - k*x - n*y)
    tt <- sum((-h)^(k+n)*exp(h)*(1-px)^n/factorial(n))
    ss <- ss + px^k / factorial(k) * tt }
  return(1 - theta/(1+theta)*ss) }
psi(2.5,0.5,1,2,0.5) ## 0.2475216



mu <- 1.5; sig2 <- 0.25; gam <- 0; 
mu3 <- gam * sig2^1.5 + 3 * mu * sig2 + mu^3
mu3.diff <- function(x)
{  y <- mu + sig2/(mu-x); px <- sig2 / (sig2+(x-mu)^2)
   px*x^3 + (1-px)*y^3 - mu3   }
x <- uniroot(mu3.diff, lower=0, upper=mu*0.9999999)$root
psi(2.5, 0.5, x, mu + sig2/(mu-x), sig2 / (sig2+(x-mu)^2))



set.seed(1); S <- exp(rnorm(1000,0,1))
levels <- c(0.5,0.95,0.975,0.99,0.995)
quantile(S, levels, type=7)
##      50%        95%      97.5%        99%      99.5% 
##0.9652926  5.7200919  7.4343584 10.0554305 11.5512498 



###
### Chapter 6 
###

BM.frac <- c(1.2,1,.9,.8,.7,.6,.55,.5,.45,.4,.375,.35,.325,.3)
Next <- rbind(  ##  see Table 6.1
   c( 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,14),
   c( 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 7, 8, 8, 9),
   c( 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 4, 4, 5),
   c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
FillP <- function (p)
{ PP <- matrix(0,nrow=14,ncol=14)
  for (k1 in 1:4) for (b in 1:14) 
  {j <- Next[k1,b]; PP[b,j] <- PP[b,j] + p[k1]}
  return(PP)}



b <- c(0,0); lbs <- c(0.05*(1-0.0001), 0.05*(1+0.0001)) 
for (i in 1:2)
{ pp <- dpois(0:2, lbs[i])
  P <- FillP(c(pp[1],pp[2],pp[3],1-sum(pp)))
  for (k in 1:10) P <- P %*% P ## i.e., P <- P^(2^10)
  stst <- P[14,] ## bottom row is near steady state db
  b[i] <- sum(stst*BM.frac)} ## b(lambda)
(log(b[2])-log(b[1])) / (log(lbs[2])-log(lbs[1]))
## = 0.1030403315



pp <- dpois(0:2, 0.05)
P <- FillP(c(pp[1],pp[2],pp[3],1-sum(pp)))
stst <- eigen(t(P))$vectors[,1]
stst <- stst/sum(stst); stst <- Re(stst)



P.prime <- FillP(c(-pp[1],pp[1]-pp[2],pp[2]-pp[3],pp[3]))
IP <- diag(14)-P
lP <- stst %*% P.prime         ## system (6.21) is g%IP = lP
IP[,14] <- 1; lP[14] <- 0      ## replace last eqn by sum(g)=0
g <- lP %*% solve(IP)          ## g = lP (IP)^{-1}
0.05*sum(g*BM.frac)/sum(stst*BM.frac)     ## = 0.1030403317



TMax <- 50; NSim <- 10000; FinalBM <- numeric(NSim)
lbs <- c(0.05*(1-0.1), 0.05*(1+0.1)); b <- c(0,0); 
for (ii in 1:2) ## just as in Method 1
{ for (n in 1:NSim)
  { cn1 <- rpois(TMax,lbs[ii]); cn1 <- pmin(cn1, 3) + 1
    BM <- 14; for (i in 1:TMax) BM <- Next[cn1[i],BM]
    FinalBM[n] <- BM  
  }
  print(table(FinalBM)/NSim); b[ii] <- mean(BM.frac[FinalBM]) 
}
(log(b[2])-log(b[1])) / (log(lbs[2])-log(lbs[1]))



###
### Chapter 7
###

nsim <- 4000; n <- 20; alpha <- rep(1, n)
mu.Y <- rep(20/100,n); Sigma.Y <- diag((13.5/100)^2, n)
L <- lower.tri(Sigma.Y, diag=TRUE)
   ## TRUE (=1) below, and on, diagonal
mu.Z <- as.vector(L %*% mu.Y)     ## means of (Z1,..Zn)
Sigma.Z <- L %*% Sigma.Y %*% t(L) ## and variance matrix
library(MASS) 
S <- exp(-mvrnorm(nsim, mu.Z, Sigma.Z)) %*% alpha



w <- rnorm(nsim) %o% sqrt(diag(Sigma.Z)) + rep(1,nsim) %o% mu.Z
S.upp <- exp(-w) %*% alpha



b <- rep(1,n) ## for the choice V = b'Y = Y1+...+Yn
C <- rbind(L,b)
mu.Zk <- (C %*% mu.Y)[1:n]; mu.V <- (C %*% mu.Y)[n+1]
Omega <- C %*% Sigma.Y %*% t(C)
sd.Zk <- sqrt(diag(Omega))[1:n]; sd.V <- sqrt(diag(Omega)[n+1])
rho.ZkV <- pmin(cov2cor(Omega)[n+1,1:n],1)



V <- rnorm(nsim, mu.V, sd.V)
EZkV <- rep(1,nsim) %o% mu.Zk + ## nsim x n matrix
           (V-mu.V) %o% (rho.ZkV * sd.Zk / sd.V)
VarZkV <- sd.Zk^2*(1-rho.ZkV^2) ## n vector
S.low <- exp(-EZkV + 0.5 * rep(1,nsim) %o% VarZkV) %*% alpha



ZkV <- EZkV + rnorm(nsim) %o% sqrt(VarZkV)
S.imp <- exp(-ZkV) %*% alpha



y <- (-0.5+1:nsim)/nsim
plot(sort(S), y, type="l", yaxp=c(0,1,1), xaxp=c(2,10,4),
     ylab="", xlab="", lwd=1.5)
lines(sort(S.upp), y, lty="solid")
lines(sort(S.imp), y, lty="dashed")
lines(sort(S.low), y, lty="dotted")
mu.S <- sum(alpha*exp(-(mu.Z-diag(Sigma.Z)/2)))
lines(c(min(S.upp),mu.S,mu.S,max(S.upp)), c(0,0,1,1))


mu.S <- sum(alpha*exp(-(mu.Z-diag(Sigma.Z)/2)))
s <- 0 ## see (2.11) in Vanduffel et al. (2008)
for (i in 1:n) for (j in 1:n)
{s <- s + alpha[i] * alpha[j] * 
          exp(- mu.Z[i] - mu.Z[j] + Sigma.Z[i,i]/2
              + Sigma.Z[j,j]/2 + Sigma.Z[i,j])}
sd.S <- sqrt(s-mu.S^2)



lines(qnorm(y, mu.S, sd.S), y, lty="longdash")



###
### Chapter 8
###

J <- 3; K <- 5; X <- scan(n=J*K)
 99.3  93.7 103.9  92.5 110.6
112.5 108.3 118.0  99.4 111.8
129.2 140.9 108.3 105.0 116.6
j <- rep(1:J, each=K); j <- as.factor(j)
X.bar <- mean(X); Xj.bar <- tapply(X, j, mean)
MSB <- sum((Xj.bar-X.bar)^2) * K / (J-1)
MSW <- sum((X-rep(Xj.bar,each=K))^2)/J/(K-1)
MSB/MSW; qf(0.95, J-1, J*(K-1))  ## 4.588 and 3.885



anova(lm(X~j))



J <- 10; K <- 5; j <- rep(1:J, each=K); j <- as.factor(j)
m <- 100; a <- 100; s2 <- 64;
set.seed(6345789) 
w <- 0.50 + runif(J*K)
X <- m + rep(rnorm(J, 0, sqrt(a)), each=K) +
         rnorm(J*K, 0, sqrt(s2/w))



anova(lm(X~j,weight=w))



w.js <- tapply(w, j, sum); w.ss <- sum(w.js)
z.j <- 1 / (1 + s2/(a*w.js)); z.s <- sum(z.j)
X.jw <- tapply(X*w, j, sum)/w.js
X.ww <- sum(X.jw * w.js) / w.ss
X.zw <- sum(X.jw * z.j) / z.s
pr.j <- z.j * X.jw + (1-z.j) * X.zw #(8.34)



m.tilde <- X.ww
SSW <- sum(w*(X-X.jw[j])^2)
s2.tilde <- SSW/J/(K-1)
SSB <- sum(w.js*(X.jw-X.ww)^2)
a.tilde <- (SSB - (J-1)*s2.tilde) / (w.ss - sum(w.js^2)/w.ss)



z.j.tilde <- 1 / (1 + s2.tilde / (a.tilde * w.js))
z.s.tilde <- sum(z.j.tilde)
X.zw.tilde <- sum(X.jw * z.j.tilde)/ z.s.tilde
pr.j.tilde <- z.j.tilde * X.jw + (1-z.j.tilde) * X.zw.tilde



a.hat <- a.tilde
repeat {
  a.hat.old <- a.hat
  z.j.hat <- 1/(1+s2.tilde/(w.js*a.hat))
  X.zw.hat <- sum(z.j.hat * X.jw) / sum(z.j.hat)
  a.hat <- sum(z.j.hat*(X.jw-X.zw.hat)^2)/(J-1)
  if (abs((a.hat-a.hat.old)/a.hat.old) < 1e-6) break}



###
### Chapter 9
###

y <- c(10,15,20,35); w <- c(300,500,700,100)
i <- c(1,1,2,2); j <- c(1,2,1,2); beta <- c(1,1)
for (iter in 1:20){
  alpha <- sqrt(tapply(w*y^2/beta[j],i,sum)/
                tapply(w*beta[j],i,sum))
  beta <- sqrt(tapply(w*y^2/alpha[i],j,sum)/
               tapply(w*alpha[i],j,sum))}



n <- scan(n=54) ## read 54 numbers into vector n
 1  8 10  8  5 11 14 12 11 10  5 12 13 12 15 13 12 24
12 11  6  8 16 19 28 11 14  4 12  8 18  3 17  6 11 18
12  3 10 18 10 13 12 31 16 16 13 14  8 19 20  9 23 27
expo <- scan(n=54) * 7 ## number of policies times 7 
10 22 30 11 15 20 25 25 23 28 19 22 19 21 19 16 18 29
25 18 20 13 26 21 27 14 16 11 23 26 29 13 26 13 17 27
20 18 20 29 27 24 23 26 18 25 17 29 11 24 16 11 22 29



sex <- as.factor(rep(1:2, each=27, len=54))
region <- as.factor(rep(1:3, each=9, len=54))
type <- as.factor(rep(1:3, each=3, len=54))
job <- as.factor(rep(1:3, each=1, len=54))
AnnClFr <- round(1000 * n/expo)
data.frame(expo, n, sex, region, type, job, AnnClFr)[1:10,]



xt <- xtabs(AnnClFr ~ sex+region+type+job)
ftable(xt, row.vars=1:2, col.vars=3:4)



glm(n/expo ~ sex+region+type+job,
    fam=poisson(link=log), wei=expo)


glm(n ~ sex+region+type+job+offset(log(expo)), 
    fam=poisson(link=log))

g <- glm(n/expo ~ 1+region+type+region:type, poisson, wei=expo)
anova(g, test="Chisq")


rm(list = ls())
fn <- "http://www1.fee.uva.nl/ke/act/people/kaas/Cars.txt"

Cars <- read.table(fn, header=TRUE)
Cars$A <- as.factor(Cars$A); Cars$R <- as.factor(Cars$R)
Cars$M <- as.factor(Cars$M); Cars$U <- as.factor(Cars$U)
Bminus1 <- Cars$B - 1; Bis14 <- as.numeric(Cars$B==14)
Cars$B <- as.factor(Cars$B); Cars$WW <- as.factor(Cars$WW)
ActualWt <- c(650,750,825,875,925,975,1025,1075,1175,1375,1600)
W <- log(ActualWt/650)[Cars$WW]
str(Cars)
attach(Cars)
ftable(xtabs(cbind(Expo, nCl, TotCl, TotPrem) ~ R+A+M+U))

100 * tapply(nCl, R, sum) / tapply(Expo, R, sum) 
## 7.649282  9.581532 12.680255
100 * tapply(nCl,list(R,A),sum) / tapply(Expo,list(R,A),sum)

g1 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14, quasipoisson, wei=Expo)
g2 <- glm(TotCl/Expo~R+A+U+W+Bminus1+Bis14+M, quasipoisson, wei=Expo)
g3 <- glm(TotCl/Expo~R+A+U+W+B, quasipoisson, wei=Expo)
anova(g1,g2)
anova(g1,g3)

glm(TotCl/Expo ~ (R+A+U+W+Bminus1+Bis14)^2, quasipoisson,
    wei=Expo)

set.seed(1); y <- rpois(10, 7+2*(1:10))
g <- glm(y~I(1:10), poisson(link=identity))
2-AIC(g)/2 == logLik(g)

l <- list(Use=U,Age=A,Area=R,Mile=M)
round(ftable(100*tapply(TotCl,l,sum)/tapply(TotPrem,l,sum)))

detach(Cars)

###
### Chapter 10 
###

Xij <- c(232,106,35,16,2, 258,115,56,27, 221,82,4, 359,71, 349)
  i <- c(  1,  1, 1, 1,1,   2,  2, 2, 2,   3, 3,3,   4, 4,   5)
  j <- c(  1,  2, 3, 4,5,   1,  2, 3, 4,   1, 2,3,   1, 2,   1)



TT <- trunc(sqrt(2*length(Xij)))
i <- rep(1:TT,TT:1); j <- sequence(TT:1)



CL <- glm(Xij~as.factor(i)+as.factor(j), family=poisson)
coefs <- exp(coef(CL)) ##exponents of parameter estimates
alpha.glm <- coefs[1] * c(1, coefs[2:TT])
beta.glm <- c(1, coefs[(TT+1):(2*TT-1)])



coefs



Ri <- tapply(Xij, i, sum); Cj <- tapply(Xij, j, sum)
alpha <- beta <- numeric(TT)
aa <- alpha[1] <- Ri[1]
bb <- beta[TT] <-  Cj[TT] / Ri[1]
for (n in 2:TT) {
   aa <- aa + (alpha[n] <- Ri[n]/(1-bb))
   bb <- bb + (beta[TT-n+1] <- Cj[TT-n+1] / aa)}
pred <- alpha %*% t(beta)



Xij.mat.cum <- Xij.mat <- matrix(0, nrow=TT, ncol=TT)
for (k in 1:length(Xij)) Xij.mat[i[k],j[k]] <- Xij[k]
for (k in 1:TT) Xij.mat.cum[k,] <- cumsum(Xij.mat[k,])



i.mat <- row(Xij.mat); j.mat <- col(Xij.mat);
future <- i.mat + j.mat - 1 > TT  
t(Xij.mat)[!t(future)] ## equals the vector Xij



Xij <- scan(n=55)
357848  766940  610542  482940 527326 574398 146342 139950 227229 67948
352118  884021  933894 1183289 445745 320996 527804 266172 425046
290507 1001799  926219 1016654 750816 146923 495992 280405
310608 1108250  776189 1562400 272482 352053 206286
443160  693190  991983  769488 504851 470639
396132  937085  847498  805037 705960
440832  847631 1131398 1063269
359480 1061648 1443370
376686  986608
344014



n <- length(Xij); TT <- trunc(sqrt(2*n))
i <- rep(1:TT, TT:1); i <- as.factor(i) ## row nrs
j <- sequence(TT:1); j <- as.factor(j)  ## col nrs
Orig.CL <- glm(Xij~i+j, quasipoisson)
coefs <- exp(as.numeric(coef(Orig.CL)))
alpha <- c(1, coefs[2:TT]) * coefs[1]
beta <- c(1, coefs[(TT+1):(2*TT-1)])
Orig.fits <- alpha %*% t(beta)
future <- row(Orig.fits) + col(Orig.fits) - 1 > TT
Orig.reserve <- sum(Orig.fits[future]) ## 18680856



Prs.resid <- (Xij-fitted(Orig.CL))/sqrt(fitted(Orig.CL))



p <- 2*TT-1; phi.P <- sum(Prs.resid^2)/(n-p)



Adj.Prs.resid <- Prs.resid * sqrt(n/(n-p))



set.seed(6345789)



nBoot <- 1000; payments <- reserves <- numeric(nBoot)
for (boots in 1:nBoot){ ## Start of bootstrap-loop



Ps.Xij <- sample(Adj.Prs.resid, n, replace=TRUE) 



Ps.Xij <- Ps.Xij * sqrt(fitted(Orig.CL)) + fitted(Orig.CL)
Ps.Xij <- pmax(Ps.Xij, 0) ## Set `observations' < 0 to 0



Ps.CL <- glm(Ps.Xij~i+j, quasipoisson)
coefs <- exp(as.numeric(coef(Ps.CL)))
Ps.alpha <- c(1, coefs[2:TT]) * coefs[1]
Ps.beta <- c(1, coefs[(TT+1):(2*TT-1)])



Ps.fits <- Ps.alpha %*% t(Ps.beta)
Ps.reserve <- sum(Ps.fits[future])



Ps.totpayments <- phi.P * rpois(1, Ps.reserve/phi.P)



reserves[boots] <- Ps.reserve
payments[boots] <- Ps.totpayments
}  ## Curly bracket indicates end of bootstrap-loop



PEbs <- sqrt(phi.P*Orig.reserve + sd(reserves)^2) ## 2882413
sd(reserves)^2 / (phi.P * Orig.reserve)           ## 7.455098



payments <- payments/1e6 ## expressed in millions
quantile(payments, c(0.5,0.75,0.9,0.95,0.99))
##       50%      75%      90%      95%      99% 
##  18.56828 20.67234 22.35558 23.61801 26.19600 
mean(payments)           ## 18.75786
sd(payments)             ## 2.873488
100 * sd(payments) / mean(payments)  ## 15.31885 = c.v. in %
pp <- (payments-mean(payments))/sd(payments)
sum(pp^3)/(nBoot-1)      ## 0.2468513 estimates the skewness
sum(pp^4)/(nBoot-1) - 3  ## 0.2701999 estimates the kurtosis



hist(payments,breaks=21,prob=TRUE)
lines(density(payments), lty="dashed")
curve(dnorm(x, mean = mean(payments), sd = sd(payments)),
   lty="dotted", add=TRUE)



Xij.1 <- xtabs(Xij~i+j) ## full square matrix
ii <- row(Xij.1); jj <- col(Xij.1); Xij.1 <- as.vector(Xij.1)
future <- as.numeric(ii+jj-1 > TT)
ii <- as.factor(ii); jj <- as.factor(jj)  ## are now vectors
Full.CL <- glm(Xij.1~ii+jj, fam=quasipoisson, wei=1-future)
Sig <- vcov(Full.CL); X <- model.matrix(Full.CL)
Cov.eta <- X%*%Sig%*%t(X)
mu.hat <- fitted(Full.CL)*future
pe2 <- phi.P * sum(mu.hat) + t(mu.hat) %*% Cov.eta %*% mu.hat
cat("Total reserve =", sum(mu.hat), "p.e. =", sqrt(pe2), "\n")
## Total reserve = 18680856 p.e. = 2945659



Xij <- scan(n=36)
156  37   6   5   3   2   1   0
154  42   8   5   6   3   0
178  63  14   5   3   1
198  56  13  11   2
206  49   9   5
250  85  28
252  44
221
TT <- trunc(sqrt(2*length(Xij)))
i <- rep(1:TT, TT:1); j <- sequence(TT:1)
ni <- c(28950,29754,31141,32443,34700,36268,37032,36637)



Ri <- tapply(Xij, i, sum); Cj <- tapply(Xij, j, sum)
alpha <- beta <- numeric(TT)
aa <- alpha[1] <- Ri[1]
bb <- beta[TT] <- Cj[TT] / Ri[1]
for (n in 2:TT) {
aa <- aa + (alpha[n] <- Ri[n]/(1-bb))
bb <- bb + (beta[TT-n+1] <- Cj[TT-n+1] / aa)}
Mi <- ni * alpha[1] / ni[1]
BF <- Mi %*% t(beta); CL <- alpha %*% t(beta)
future <- row(BF) + col(BF) - 1 > TT
rowSums(BF * future) ## 0.0  0.0  0.5  2.6  6.4 13.2 26.3 76.9
rowSums(CL * future) ## 0.0  0.0  0.6  3.1  7.0 19.2 32.1 90.0



Expo <- ni[i] ## Exposures with each element of Xij
CL <- glm(Xij~as.factor(i)+as.factor(j), poisson)
CLoff <- glm(Xij~offset(log(Expo))+as.factor(j), poisson)



fitted(CL)
alpha[i]*beta[j]
alpha*beta
alpha%o%beta
alpha%*%beta
outer(alpha,beta)
alpha%*%t(beta)



###
### Chapter 11
###

w <- c(1,2,2,1,1,1,1,4,2,4,2,3,2,1,1,2)
Y <- c(0,1,0,8,0,0,0,30,0,1,1,38,0,0,0,26) / w
gender <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
income <- c(1,2,5,20,1,2,5,20,1,2,5,20,1,2,5,20)



lm(Y ~ gender+income, weights=w)
glm(Y ~ gender+income, weights=w, family=quasipoisson)



X <- cbind(rep(1,10),1:10)
y <- c(14,0,8,8,16,16,32,18,28,22)
z <- log(y+(y==0))
W <- diag(10)
b <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z 
cat("Start:", b[1], b[2], "\n")



for (it in 1:5){
 eta <- as.vector(X %*% b)
 mu <- exp(eta)                  ## eta = g(mu) = log(mu)
 W <- diag(mu)           ## (g'(mu))^(-2)/V(mu) = mu^2/mu
 z <- X %*% b + (y-mu)/mu   ## d eta/d mu = g'(mu) = 1/mu
 S <- solve(t(X) %*% W %*% X)
 b <- S %*% t(X) %*% W %*% z
 cat("it =", it, b[1], b[2], "\n")}



coef(glm(y~I(1:10),quasipoisson))



dTweedie <- function (y, power, mu, phi)
{ if (power==2) s <- dgamma(y, 1/phi, 1/(phi*mu)) else
  if (power==1) s <- dpois(y/phi, mu/phi) else
  { lambda <- mu^(2-power)/phi/(2-power)
    if (y==0) s <- exp(-lambda) else 
    { alpha <- (2-power)/(power-1)
      beta <- 1 / (phi * (power-1) * mu^(power-1))
      k <- max(10, ceiling(lambda + 7*sqrt(lambda)))
      s <- sum(dpois(1:k,lambda) * dgamma(y,alpha*(1:k),beta))
  } }
  return(s) }


require(statmod)###If "FALSE" results, download it from CRAN first
TT <- 10; i <- rep(1:TT, each=TT); j <- rep(1:TT, TT)
past <- i + j - 1 <= TT; n <- sum(past)
Expo <- c(100, 110, 115, 120, 130, 135, 130, 140, 130, 120)
Runoff <- c(30, 30, 10, 20, 5, 3, 1, 0.5, 0.3, 0.2)
Off <- rep(Expo, each=TT) * rep(Runoff, TT); lOff <- log(Off)
##note that future values are input as 0.01; they get weight 0 anyway
Xij <- scan(n=100)
4289.93 3093.71 1145.72 1387.58 293.92 189.17  42.36 11.41 4.31 12.39
3053.09 2788.81  682.44 1475.69 253.31 100.58  79.35 15.48 8.06  0.01
4388.93 2708.67  688.42 2049.57 353.20 266.43 109.42 47.90 0.01  0.01
4144.15 2045.63 1642.27 1310.97 548.97 159.87  69.86  0.01 0.01  0.01
2912.73 4078.56 1652.28 2500.94 394.99 220.89   0.01  0.01 0.01  0.01
5757.18 5200.83 1177.65 2486.30 580.29   0.01   0.01  0.01 0.01  0.01
4594.18 3928.15 1236.01 2729.68   0.01   0.01   0.01  0.01 0.01  0.01
3695.03 3688.23 1300.97    0.01   0.01   0.01   0.01  0.01 0.01  0.01
3967.13 4240.97    0.01    0.01   0.01   0.01   0.01  0.01 0.01  0.01
4933.06    0.01    0.01    0.01   0.01   0.01   0.01  0.01 0.01  0.01
round(xtabs(Xij~i+j)) ## produces a table of the input values
y <- Xij[past]
Tweedie.logL <- function(pow)
{ gg <- glm(Xij~i+j+offset(lOff), tweedie(pow,0), wei=as.numeric(past))
  reserve <- sum(fitted.values(gg)[!past])
  dev <- deviance(gg); phi.hat <- dev/n
  mu <- fitted.values(gg)[past]; hat.logL <- 0
  for (ii in 1:length(y))
  { hat.logL <- hat.logL + log(dTweedie(y[ii], pow, mu[ii], phi.hat)) }
  cat("Power =", round(pow,3), "\tphi =", round(phi.hat,2), 
      "\tRes. =", round(reserve), "\tlogL =", round(hat.logL,3), "\n")
  hat.logL   }
for (pow in c(1,1.25,1.5,1.75,2)) Tweedie.logL(pow)
oo <- optimize(Tweedie.logL, c(1.01,1.99), maximum=T, tol=1e-4)


## add the following 2 lines:
lm. <- lm(Y ~ gender+income, weights=w)
model.matrix(lm.)

mean(lm.$residuals^2*w) *
 solve(t(model.matrix(lm.))%*%diag(w)%*%model.matrix(lm.)) /
    vcov(lm.)



Size  <- c(59,60,62,56,63,59,62,60)
Score <- c(69,72,76,78,81,84,86,88)
Def   <- c( 6,13,18,28,52,53,61,60)



hist(rtweedie(1000, power=1.001, mu=1, phi=1), breaks=41)
sum(dtweedie((1:1999-.5)/100, 1.5, 1, 1)) / 100 +
  dtweedie(0, 1.5, 1, 1)
