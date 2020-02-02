# Advanced Econometrics I - Part 2
# Problem Set 5 - Due 6th February 2020

library(MASS)
library(bbmle)
library(margins) # replicates STATA's margins command
library(dplyr)
library(stargazer)
library(ggplot2)
library(latex2exp)

#########################################################
# Exercise 2: DGP
#########################################################

set.seed(1)

# set parameters
b0 <- 5
b1 <- 3
b2 <- 2
g0 <- 1
g1 <- 2
g2 <- 1.5

# specify joint distribution for X
means_x <- c(0, 0)
cov_mat_x <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)

# specify joint distribution for u1, v2
means_e <- c(0, 0)
rho <- 0
cov_mat_e <- matrix(c(1, 2*rho, 2*rho, 4), nrow = 2, ncol = 2)

n <- 10000

# generate data
multnorm_x <- mvrnorm(n, mu = means_x, Sigma = cov_mat_x)
x1 <- multnorm_x[,1]
x2 <- multnorm_x[,2]

multnorm_e <- mvrnorm(n, mu = means_e, Sigma = cov_mat_e)
u1 <- multnorm_e[,1]
v2 <- multnorm_e[,2]

# compute y1_star and y2
y2 <- g0 + g1*x1 + g2*x2 + v2
y1_star <- b0 + b1*x1 + b2*y2 + u1
y1 <- vector(length = n)
for (i in 1:n) {
  y1[i] <- ifelse(y1_star[i] > 0, y1[i] <- 1, y1[i] <- 0)
}


#########################################################
# Exercise 3: rho = 0, no covariance btw. x1, x2 (exogeneity)
#########################################################

# (a) Perform the two-step approach described on slide 187.

# 1. OLS of y2 on x1 and x2 to obtain residuals v2_hat
model_ols <- lm(y2 ~ x1 + x2)
v2_hat <- model_ols$residuals

# 2. Probit of y1 on x1, x2, v2_hat
model_probit <- glm(y1 ~ x1 + y2 + v2_hat, family = binomial(link = "probit"))
model_lpm <- lm(y1 ~ x1 + x2 + v2_hat)


# (b) Test the H0 that y2 is exogenous.

theta1_hat <- cov(v2, u1) / var(v2)
theta1 <- 0 # H0
t <- (theta1_hat - theta1) / (sd(v2) / sqrt(n))
t.test(v2, mu = 0) # can't reject, so v2 is likely exogenous

# (c) Are the coeffcients for y2 and x1 in the probit model as expected?

summary(model_probit)
rbind(coefficients(model_probit), c(b0, b1, b2, rho))
# b0, b1, b2 almost identical to true parameters
# what can we expect for v2_hat?

# (d) Does the procedure work without an identifying variable x2. Why?

model_ols2 <- lm(y2 ~ x1)
v2_hat2 <- model_ols2$residuals
model_probit2 <- glm(y1 ~ x1 + y2 + v2_hat2, family = binomial(link = "probit"))
summary(model_probit2)
# it doesn't work, because of singularities
# all other determinants of y2 would be captured by the constant if x2 weren't included


#########################################################
# Exercise 4: rho = 0.7, endogeneity
#########################################################

set.seed(1)

rho2 <- 0.7
cov_mat_e2 <- matrix(c(1, 2*rho2, 2*rho2, 4), nrow = 2, ncol = 2)
multnorm_e2 <- mvrnorm(n, mu = means_e, Sigma = cov_mat_e2)
u12 <- multnorm_e2[,1]
v22 <- multnorm_e2[,2]

# compute y1_star2 and y22
y22 <- g0 + g1*x1 + g2*x2 + v22
y1_star2 <- b0 + b1*x1 + b2*y22 + u12
y12 <- vector(length = n)
for (i in 1:n) {
  y12[i] <- ifelse(y1_star2[i] > 0, y12[i] <- 1, y12[i] <- 0)
}


# (a) Regress u1 on v2 and calculate the residuals. Is the variance of
# these residuals as suggested on page 185?
model_ols3 <- lm(u1 ~ v2 - 1)
e1 <- residuals(model_ols3)
mean(e1) # close to 0 as expected
# var(e1) = var(u1) - theta1^2 * var(v2) = 1 - rho2^2
var(u1) - coefficients(model_ols3)[1]^2 * var(v2) # holds (= var(e1))
# last part of the equation doesn't hold; why?
1 - rho2^2

# (b) Perform the two-step approach again.

# 1. OLS of y22 on x1 and x2 to obtain residuals v2_hat22
model_ols4 <- lm(y22 ~ x1 + x2)
v2_hat22 <- model_ols4$residuals

# 2. Probit of y1 on x1, x2, v2_hat22
model_probit22 <- glm(y12 ~ x1 + y22 + v2_hat22, family = binomial(link = "probit"), maxit = 1000)
model_lpm22 <- lm(y12 ~ x1 + x2 + v2_hat22)


# (c) Test the H0 that y2 is exogenous. Is the coeffcient for v2 as
# expected?



# (d) Compare the coeffcient estimates for y2 and x1 with equation (1).
# Can you explain the difference?



# (e) What is the APE of y2 on y1?



# (f) Run an LPM model using 2SLS and compare the marginal effect
# of y2 to the APE above?