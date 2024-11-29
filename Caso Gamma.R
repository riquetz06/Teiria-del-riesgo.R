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
> # Calculate points for density curves
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
