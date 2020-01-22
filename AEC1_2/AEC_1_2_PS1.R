# Plotting the Likelihood and Loglikelihood of a Bernoulli distribution

data = c(0, 1, 1, 1, 0, 1, 1, 0, 1, 1)

log.likelihood <- function(data, theta){
  sum(dbinom(x = data, size = 1, prob = theta, log = T))
}

theta = seq(0, 1, 0.01)
lls <- vector(mode = "numeric", length = length(theta))
for(i in 1:length(theta)){
  lls[i] <- log.likelihood(data, theta[i])
  }

plot(theta, lls, type = "l")

plot(theta, exp(lls), type = "l")

lls_finite <- lls[2:100]
var(exp(lls_finite))
