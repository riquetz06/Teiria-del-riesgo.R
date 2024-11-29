#Convoluciones por método recursivo
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
