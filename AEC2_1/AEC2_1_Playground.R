
# AEC 2.1: Univariate Time Series Econometrics ----------------------------

# Playground


# Wiener Process ----------------------------------------------------------

library(e1071)

wp <- rwiener(end = 10)
plot(wp)
var(wp)
mean(wp)


# White Noise -------------------------------------------------------------

wn <- rnorm(1000)
plot(wn)
var(wn)

# example of stationarity without ergodicity
# note (characterization of ergodicity):
# - as t->inf, mean(x_t) -> (a.s.) mu
# - all measurable, strictly stationary, and ergodic mappings f_t of a process x_t also converges a.s. to E(f(x_t))

# define x_t = \epsilon_t + delta, where delta ~ N(0,1)
# delta shifts the white noise process eps_t by delta, which is normally distributed, but time-independent

eps_t <- rnorm(1000)
delta <- rnorm(1)

x_t <- eps_t + delta
plot(x_t)

# this process is stationary, but not ergodic:
# E(x_t) -> (a.s.) mu + delta != mu


# Conditional Heteroskedasticity ------------------------------------------

# Consider x_t = sqrt(h_t) + e_t, where e_t ~ iid(0,1)
# and h_t = f(x_{t-1}, x_{t-2}, ...)
n <- 1000

library(MASS)
df <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
e_t <- df[,1]
h_t <- abs(df[,2])

e_t <- abs(rnorm(n))
h_t <- lag(e_t, 1) + lag(e_t, 2) + lag(e_t, 3)
x_t <- sqrt(h_t) + e_t
plot(x_t, h_t) # see heteroskedasticity

# (G)ARCH model accounting for conditional heteroskedasticity
x_t_sq <- x_t^2
arch_model <- lm(h_t ~ x_t_sq + lag(x_t_sq, 1) + lag(x_t_sq, 2) + lag(x_t_sq, 3))
summary(arch_model)

# variance gamma(0)
gamma0 <- arch_model$coefficients[1] / (1 - arch_model$coefficients[2] - arch_model$coefficients[3] - arch_model$coefficients[4])
