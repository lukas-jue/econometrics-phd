################################################################
### Advanced Econometrics 1 – PhD Course
### Winter Term 2019
### Problem Set 2
################################################################

################################################################
# Question 1 – Homogeneous IV Setup
################################################################
# see slides chapter 3 for notation

# 0. Data Generating Process
set.seed(123)
v <- rnorm(200)
e <- rnorm(200)
z <- rnorm(200)

d <- 0.25 * z + v
u <- 0.5 * v + e   # create error term which is highly correlated to regressor (model endogeneity problem)
y <- 2 * d + u

# 1. Derive OLS estimator from regression of y on d w/o constant
ols_y_d <- lm(y ~ 0 + d)
summary(ols_y_d)
beta_hat <- coefficients(ols_y_d)["d"] # biased estimator for influence of d on y (d correlated to error term u)

# 2. Derive the IV-estimator where d is instrumented by z
cor(y, z) # check relevance condition
cor(d, z)
cor(d, u) # check if error term is correlated to regressor (if yes, we have an endogeneity problem)
cor(z, u) # check if instrument z is exogenous

first_iv <- lm(d ~ z)
d_hat <- fitted.values(first_iv)
pi_hat <- as.numeric(coefficients(first_iv)["z"])

reduced_form <- lm(y ~ d_hat)
phi_hat_1 <- as.numeric(coefficients(reduced_form)["d_hat"]) # unbiased estimator for influence of d on y through usage of instrument z
print("true beta = 2.000")
paste("biased beta hat = ", beta_hat)
paste("unbiased beta_hat_iv = : " ,phi_hat_1)

# 3. Perform a Monte Carlo study to to study the behavior of beta_OLS and beta_IV for a sample size N = 200
library(MonteCarlo)

# first, write code above as a function
simulate_iv <- function(n){
  v <- rnorm(n)
  e <- rnorm(n)
  z <- rnorm(n)
  d <- 0.25 * z + v
  u <- 0.5 * v + e
  y <- 2 * d + u
  ols_y_d <- lm(y ~ 0 + d)
  beta_hat <- coefficients(ols_y_d)["d"]
  first_iv <- lm(d ~ z)
  d_hat <- fitted.values(first_iv)
  pi_hat <- as.numeric(coefficients(first_iv)["z"])
  reduced_form <- lm(y ~ d_hat)
  phi_hat_1 <- as.numeric(coefficients(reduced_form)["d_hat"])
  
  return(phi_hat_1)
}

# single Monte Carlo simulation
n <- 1000 #MC runs
mc_beta_iv <- vector(length = n)
for (i in 1:n) {
  mc_beta_iv[i] <- simulate_iv(200)
}
hist(mc_beta_iv, xlim = c(-2, 4), breaks = 100)
plot(density(mc_beta_iv), xlim = c(0, 4))

# multiple MC simulations to compare distributions of beta_iv's
n_grid <- c(2, 5, 10, 50, 100, 500, 1000)

mc_beta_iv_list <- list()
for (i in 1:length(n_grid)) {
  vector <- vector()
  for (k in 1:n_grid[i]) {
    vector[[k]] <- simulate_iv(n_grid[i])
  }
  mc_beta_iv_list[[i]] <- vector
}

dim_graphs <- 3  #round(sqrt(length(mc_beta_iv_list)))
par(mfcol = c(dim_graphs, dim_graphs))
for (i in 1:length(mc_beta_iv_list)) {
  hist(mc_beta_iv_list[[i]], breaks = 50, xlim = c(0, 4), main = paste("n = " , length(mc_beta_iv_list[[i]])))
  abline(v = 2, col = "red")
}


# 4. Repeat 3. but change the data generating process (dgp) of u: What does this dgp imply? u = 0.5v + 0.5z + e
# Answer: Estimates for beta_iv are upward biased. No valid instrument anymore, since error u is now correlated with instrument z
simulate_iv_4 <- function(n){
  v <- rnorm(n)
  e <- rnorm(n)
  z <- rnorm(n)
  d <- 0.25 * z + v
  u <- 0.5 * v + 0.5 * z + e # only change is here. Error is now correlated with instrument z. Before u <- 0.5 * v + e
  y <- 2 * d + u
  ols_y_d <- lm(y ~ 0 + d)
  beta_hat <- coefficients(ols_y_d)["d"]
  first_iv <- lm(d ~ z)
  d_hat <- fitted.values(first_iv)
  pi_hat <- as.numeric(coefficients(first_iv)["z"])
  reduced_form <- lm(y ~ d_hat)
  phi_hat_1 <- as.numeric(coefficients(reduced_form)["d_hat"])
  
  return(phi_hat_1)
}

mc_beta_iv_list <- list()
for (i in 1:length(n_grid)) {
  vector <- vector()
  for (k in 1:n_grid[i]) {
    vector[[k]] <- simulate_iv_4(n_grid[i])
  }
  mc_beta_iv_list[[i]] <- vector
}

dim_graphs <- 3
par(mfcol = c(dim_graphs, dim_graphs))
for (i in 1:length(mc_beta_iv_list)) {
  hist(mc_beta_iv_list[[i]], breaks = 50, xlim = c(0, 4), main = paste("n = " , length(mc_beta_iv_list[[i]])))
  abline(v = 2, col = "red")
}

# 5. Repeat 3 but change the dgp of d: What does this dgp imply? d = a/N * z + v
# What is a? Assumed constant. If a/n is large, then less variance for smaller n MC simulations of beta_iv

simulate_iv_5 <- function(n){
  v <- rnorm(n)
  e <- rnorm(n)
  z <- rnorm(n)
  a <- 100
  d <- a/n * z + v # only change here. Before: d <- 0.25 * z + e
  u <- 0.5 * v + e
  y <- 2 * d + u
  ols_y_d <- lm(y ~ 0 + d)
  beta_hat <- coefficients(ols_y_d)["d"]
  first_iv <- lm(d ~ z)
  d_hat <- fitted.values(first_iv)
  pi_hat <- as.numeric(coefficients(first_iv)["z"])
  reduced_form <- lm(y ~ d_hat)
  phi_hat_1 <- as.numeric(coefficients(reduced_form)["d_hat"])
  
  return(phi_hat_1)
}

mc_beta_iv_list <- list()
for (i in 1:length(n_grid)) {
  vector <- vector()
  for (k in 1:n_grid[i]) {
    vector[[k]] <- simulate_iv_5(n_grid[i])
  }
  mc_beta_iv_list[[i]] <- vector
}

dim_graphs <- 3
par(mfcol = c(dim_graphs, dim_graphs))
for (i in 1:length(mc_beta_iv_list)) {
  hist(mc_beta_iv_list[[i]], breaks = 50, xlim = c(0, 4), main = paste("n = " , length(mc_beta_iv_list[[i]])))
  abline(v = 2, col = "red")
}

