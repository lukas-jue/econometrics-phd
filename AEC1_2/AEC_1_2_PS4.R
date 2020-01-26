# Advanced Econometrics I - Part 2
# Problem Set 4 - Due 30th January 2020

library(MASS)
library(bbmle)
library(margins) # replicates STATA's margins command
library(dplyr)
library(stargazer)

set.seed(1)

# set parameters and mean, variance
b0 <- 0.5
b1 <- 3
b2 <- 2
means <- c(0, 0)
cov_mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
n <- 10000

# generate data
multnorm <- mvrnorm(n, mu = means, Sigma = cov_mat)
x1 <- multnorm[,1]
x2 <- multnorm[,2]
e <- rnorm(n, mean = 0, sd = 1)

# compute y_star
y_star <- b0 + b1*x1 + b2*x2 + e

# translate y_star into y
y <- vector(length = n)
for (i in 1:n) {
  y[i] <- ifelse(y_star[i] > 0, y[i] <- 1, y[i] <- 0)
}

#########################################################
# Exercise 1: Estimate beta using LPM
#########################################################

model_lpm <- lm(y ~ x1 + x2)
summary(model_lpm)
# use heteroskedasticity-robust standard errors, because LPM errors are heteroskedastic by design
lmtest::coeftest(model_lpm, vcov. = sandwich::vcovHC, type = "HC0") # HC0: White Estimator

plot(x1, y)
  abline(model_lpm$coefficients[c(1, 2)])
plot(x2, y)
  abline(model_lpm$coefficients[c(1, 3)])
  
e_hat <- y - model_lpm$coefficients[1] - model_lpm$coefficients[2] * x1 - model_lpm$coefficients[3] * x2
y_hat_lpm <- model_lpm$coefficients[1] + model_lpm$coefficients[2] * x1 + model_lpm$coefficients[3] * x2

plot(y_hat, e_hat) # residuals vs. fitted
plot(x1, e_hat)
plot(x2, e_hat)
plot(y_hat_lpm, y)

var_y <- y_hat * (1 - y_hat)
par(mfrow = c(2, 2))
plot(x1, var_y)
plot(x2, var_y)
hist(var_y)

# Graphical test for heteroskedasticity
par(mfrow = c(2, 2))
plot(model_lpm)
# top left plot: Residuals clearly depend on fitted values
# => heteroskedasticity is present

# verify by Breusch-Pagan Test
lmtest::bptest(model_lpm)

sq_error <- e_hat^2
model_bp <- lm(sq_error ~ x1 + x2)
summary(model_bp)


#########################################################
# Exercise 2: Estimate beta using Probit
#########################################################

model_probit <- glm(y ~ x1 + x2, family = binomial(link = "probit"))
summary(model_probit)

# slide 130, 132
xb_probit <- model_probit$coefficients[1] + model_probit$coefficients[2] * x1 + model_probit$coefficients[3] * x2

y_hat_probit <- pnorm(xb_probit) # phi
# y_hat_probit <- (1/sqrt(2*pi)) * exp(-(xb_probit)^2 / 2) # should be the same to pnorm, but isn't
z_sq <- xb_probit^2

plot(x1, phi)
plot(x2, phi)

y_hat_logit <- exp(xb_probit) / (1 - exp(xb_probit))^2
plot(x1, y_hat_probit)


#########################################################
# ML estimation of beta (LPM)

LL <- function(beta0, beta1, beta2, mu, sigma) {
  e = y - x1 * beta1 - x2 * beta2 - beta0
  #
  e = suppressWarnings(dnorm(e, mu, sigma, log = TRUE))
  #
  -sum(e)
}

model_mle_lpm <- mle2(LL, start = list(beta0 = 1, beta1 = 1, beta2 = 1, mu = 0, sigma = 1))
summary(model_mle_lpm)

#########################################################
# ML estimation of beta (Probit) 
# Results don't seem too convincing, likely contains errors

LL_probit <- function(beta0, beta1, beta2, mu, sigma) {
  z = x1 * beta1 + x2 * beta2 + beta0
  e = y - (1/sqrt(2*pi)) * exp(-(z)^2 / 2)
  #
  e = suppressWarnings(dnorm(e, mu, sigma, log = TRUE))
  #
  -sum(e)
}

model_mle_probit <- mle2(LL_probit, start = list(beta0 = 1, beta1 = 1, beta2 = 1, mu = 0, sigma = 1))
summary(model_mle_probit)


#########################################################
# Exercise 3
#########################################################

# a) Compare P(y = 1 | x) btw. LPM and Probit
plot(y_hat_lpm, y_hat_probit)

# b) Obtain the Average Partial Effects (APEs) of x1 and x2
model_lpm %>% 
  margins() %>% 
  summary()

model_probit %>% 
  margins() %>% 
  summary()


par(mfrow = c(1, 2))
model_lpm %>% 
  margins() %>% 
  plot(main = "LPM")

model_probit %>% 
  margins() %>% 
  plot(main = "Probit")

# c) Obtain the APE of x1 by hand
xb_probit <- model_probit$coefficients[1] + model_probit$coefficients[2] * x1 + model_probit$coefficients[3] * x2
APE_b1 <- model_probit$coefficients[2] * 1/length(x1) * sum(dnorm(xb_probit))

# d) Obtain the APE of x1, evaluated at x2 = 0; 1; 2; 4 using the margins command
model_probit %>% 
  margins(at = list(x2 = c(0, 1, 2, 4))) %>% 
  summary()

# e) Calculate - by hand - the APE of x1 evaluated at x2 = 0; 1; 2; 4.
# NEED TO FIX!
xb <- model_probit$coefficients[1] + model_probit$coefficients[2] * x1 + model_probit$coefficients[3] * x2
model_probit$coefficients[2] * model_probit$coefficients[3] * 2



#########################################################
# Exercise 4: Estimate Logit and compare to Probit
#########################################################

model_logit <- glm(y ~ x1 + x2, family = binomial(link = "logit"))
summary(model_logit)

model_logit %>% 
  margins() %>% 
  summary()
# result: almost identical AME

par(mfrow = c(1, 2))
model_logit %>% 
  margins() %>% 
  plot(main = "Logit")

model_probit %>% 
  margins() %>% 
  plot(main = "Probit")

#########################################################
# Exercise 5: Run DGP again with e ~ N(0,2) and compare to benchmark probit
#########################################################

e2 <- rnorm(n, mean = 0, sd = 2)

# compute y_star
y_star2 <- b0 + b1*x1 + b2*x2 + e2

# translate y_star into y
y2 <- vector(length = n)
for (i in 1:n) {
  y2[i] <- ifelse(y_star2[i] > 0, y2[i] <- 1, y2[i] <- 0)
}


model_lpm2 <- lm(y2 ~ x1 + x2)
model_probit2 <- glm(y2 ~ x1 + x2, family = binomial(link = "probit"))
model_logit2 <- glm(y2 ~ x1 + x2, family = binomial(link = "logit"))


summary(model_probit2)
stargazer(model_lpm, model_lpm2, model_probit, model_probit2, model_logit, model_logit2, type = "text")

# note: OLS not very sensitive to increase in error variance
# However, Probit and Logit models decrease their coefficients by ~50%
# Is there a direct connection to Var(e)?


#########################################################
# Exercise 6: Run DGP again with e ~ N(0,2) and compare to benchmark probit
#########################################################