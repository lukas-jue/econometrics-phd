# Advanced Econometrics I - Part 2
# Problem Set 6 - Due 13th February 2020
 
# TO DO
  # (d) Correlated random effects

# Use the data-set PSID80 92:dta. We are interested in the binary de-
# pendent variable inlf, which indicates whether a woman was in the
# labor force. The x vector of explanatory variables consist of marr;
# nwif einc; ch0 2; ch3 5; ch6 17; and educ (see the variable labels for
# the meaning) along with a full set of year dummies.

library(dplyr)
library(haven)
library(plm)
library(pglm)
library(stargazer)
library(PanelCount)

df <- read_dta(file = "AEC1_2/PSID80_92.dta")

# 1. Estimate the following static models and compare the partial effects of
# the coeffcient on the marriage dummy:
# (a) LPM with fixed effects

lpm_fe <- plm(inlf ~ marr + nwifeincb + ch0_2 + ch3_5 + ch6_17 + educ, 
            data = df,
            index = c("id", "year"), 
            model = "within")

# (b) Pooled probit

probit_pooled <- glm(inlf ~ marr + nwifeincb + ch0_2 + ch3_5 + ch6_17 + educ, 
                     data = df,
                     family = binomial("probit"))

# (c) Random effects probit

probit_re <- pglm(inlf ~ marr + nwifeincb + ch0_2 + ch3_5 + ch6_17 + educ,
                  data = df,
                  model = ("random"),
                  effect = "individual",
                  index = c("id", "year"),
                  family = binomial("probit"))

summary(probit_re)

# (d) Correlated random effects

# lpm_cre <- CRE

# (e) Fixed effects logit

logit_fe <- glm(inlf ~ marr + nwifeincb + ch0_2 + ch3_5 + ch6_17 + educ, 
                     data = df,
                     family = binomial("logit"))

stargazer(lpm_fe, probit_pooled, logit_fe, type = "text")

# 2. Now define a lag of inlf, say inlf1 and estimate the following dynamic
# models.

df <- df %>%
  group_by(id) %>%
  mutate(inlf_lag = dplyr::lag(inlf, n = 1, default = NA))

# (a) Estimate a dynamic linear probability model How strong is the
# estimated persistence of being in the labor force?

# (b) Now estimate the model with fixed effects.

# (c) Next, estimate the linear model using the Arellano-Bond approach.
# Is the relationship of the coeffcient estimates for 2a, 2b and 2c as
# expected?

# (d) Estimate a dynamic probit model using the same explanatory vari-
#   ables as in part 2a, i.e. unobserved heterogeneity is not taken
# account of. Compute the average partial effect of inlf1.

# (e) Using the Wooldrige-approach, estimate a dynamic probit model
# that allows for unobserved heterogeneity. How does the magnitude
# of the coeffcient on inlf1 compare with that in 2d?