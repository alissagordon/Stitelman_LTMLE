library(tidyverse)
library(SuperLearner)

## DGP 1
n <- 500

W <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, plogis(1.5*W - 1))

# Time 1 covariates
L1 <- rbinom(n, 1, plogis(1.5*W + 3*A - W*A - 2))
C1 <- rbinom(n, 1, 0.1)

# Time 2 covariates (only for those not censored at time 1)
L2 <- rbinom(n, 1, plogis(1.5*W + 2*A - 2*W*A + L1 - 1))
C2 <- rbinom(n, 1, 0.1)
L2_obs <- ifelse(C1==1, NA, L2)

# Outcome Y only observed if C1 and C2 are both 0
Y <- rbinom(n, 1, plogis(
  1.5*W + 2*A - 1.5*W*A*L1 + 3*L1 - 1.5*L1*A - 2 - L2*W + 2*L2
))
Y_obs <- ifelse(C1 == 1 | C2 == 1, NA, Y)

# Final data frame with correct temporal ordering
dat <- data.frame(W, A, L1, C1, L2_obs, C2, Y_obs)


# Implementation of Stitelman LTMLE 

# Targeting Y
# Set intervention level
abar <- 1

# Define SL library
SL.library <- c("SL.glm.interaction", "SL.mean")

# Fit g_A: P(A = 1 | W)
gA_fit <- SuperLearner(Y = dat$A, X = data.frame(W = dat$W),
                       family = binomial(), SL.library = SL.library)
gA_pred <- predict(gA_fit, newdata = data.frame(W = dat$W), type = "response")$pred

# Fit g_C1: P(C1 = 1 | A, W, L1)
gC1_fit <- SuperLearner(Y = dat$C1, X = dat[, c("A", "W", "L1")],
                        family = binomial(), SL.library = SL.library)
gC1_pred <- predict(gC1_fit, newdata = dat[, c("A", "W", "L1")], type = "response")$pred

# Fit g_C2: P(C2 = 1 | A, W, L1, L2)
gC2_fit <- SuperLearner(Y = dat$C2[!is.na(dat$L2_obs)], X = dat[!is.na(dat$L2_obs), c("A", "W", "L1", "L2_obs")],
                        family = binomial(), SL.library = SL.library)
gC2_pred <- predict(gC2_fit, newdata = dat[, c("A", "W", "L1", "L2_obs")], type = "response")$pred

# Now construct the clever covariate:
# H_Y = I(A==1 & C1==0 & C2==0) / [gA * (1 - gC1) * (1 - gC2)]

clever_covariate_y_g <- (dat$A == 1 & dat$C1 == 0 & dat$C2 == 0) /
                    (gA_pred * (1 - gC1_pred) * (1 - gC2_pred))

# Fit initial Q_Y
qY_fit <- SuperLearner(Y = dat$Y_obs[!is.na(dat$Y_obs)],
                       X = dat[!is.na(dat$Y_obs), c("L2_obs", "L1", "A", "W")],
                       family = binomial(), SL.library = SL.library)

# Predictions for different histories of L(t) and updates
L11 = dat
L11$C1 = 0
L11$C2 = 0
L11$L1 = 1
L11$L2_obs = 1
L11$A = 1
qy.fit.L11 = predict(qY_fit,type = "link",newdata = L11)$pred
qy.fit.L11.update = glm(Y ~ offset(qy.fit.L11)+1, weights = clever_covariate_y_g, family = "quasibinomial", data = L11)
qy.star.L11 = predict(qy.fit.L11.update, type = "response", newdata=dat)

# Update P(Y = 1, C1 = 0, C2 = 0, A = 1 | L2 = 1, L1 = 0, W)
L01 = dat
L01$C1 = 0
L01$C2 = 0
L01$L1 = 1
L01$L2_obs = 0
L01$A = 1
qy.fit.L01 = predict(qY_fit,type = "link",newdata = L01)$pred
qy.fit.L01.update = glm(Y ~ offset(qy.fit.L01)+1, weights = clever_covariate_y_g, family = "quasibinomial", data = L01)
qy.star.L01 = predict(qy.fit.L01.update, type = "response", newdata=dat)

# Update P(Y = 1, C1 = 0, C2 = 0, A = 1 | L2 = 0, L1 = 1, W)
L10 = dat
L10$C1 = 0
L10$C2 = 0
L10$L1 = 0
L10$L2_obs = 1
L10$A = 1
qy.fit.L10 = predict(qY_fit,type = "link",newdata = L10)$pred
qy.fit.L10.update = glm(Y ~ offset(qy.fit.L10)+1, weights = clever_covariate_y_g, family = "quasibinomial", data = L10)
qy.star.L10 = predict(qy.fit.L10.update, type = "response", newdata=dat)

# Update P(Y = 1, C1 = 0, C2 = 0, A = 1 | L2 = 0, L1 = 0, W)
L00 = dat
L00$C1 = 0
L00$C2 = 0
L00$L1 = 0
L00$L2_obs = 0
L00$A = 1
qy.fit.L00 = predict(qY_fit,type = "link",newdata = L00)$pred
qy.fit.L00.update = glm(Y ~ offset(qy.fit.L00)+1, weights = clever_covariate_y_g, family = "quasibinomial", data = L00)
qy.star.L00 = predict(qy.fit.L00.update, type = "response", newdata=dat)

# Targeting L2

# Update P(L2 = 1| C1 = 0, L1 = 1, A1 = 1, W)
Cl2q1 = qy.star.L11 - qy.star.L01
Cl2q0 = qy.star.L10 - qy.star.L00
Cl2g = as.numeric(dat$A == abar & dat$C1 == 0)/(gA_pred*(1-gC1_pred))


qL2_fit <- SuperLearner(Y = dat$L2_obs[!is.na(dat$L2_obs)], X = dat[!is.na(dat$L2_obs), c("A", "W", "L1")],
                        family = binomial(), SL.library = SL.library)
#Q_L2_init <- predict(qL2_fit, newdata = dat[, c("A", "W", "L1")], type = "response")$pred
#logit_Q_L2 <- qlogis(Q_L2_init)

#Make predictions for L1=1 and L0=0
L1=dat
L1$A=abar
L1$C1=0
L1$L1=1
ql2.fit.L1 = predict(qL2_fit,type = "link",newdata = L1)$pred
ql2.fit.L1.update = glm(L2_obs ~ offset(ql2.fit.L1)+Cl2q1, weights = Cl2g, family = "quasibinomial", data = L1)
ql2.star.L1 = predict(ql2.fit.L1.update, type = "response", newdata=dat)


L0=dat
L0$A=abar
L0$C1=0
L0$L1=0
ql2.fit.L0 = predict(qL2_fit,type = "link",newdata = L0)$pred
ql2.fit.L0.update = glm(L2_obs ~ offset(ql2.fit.L0)+Cl2q0, weights = Cl2g, family = "quasibinomial", data = L0)
ql2.star.L0 = predict(ql2.fit.L0.update, type = "response", newdata=dat)

C_Q <- qy.star.L01*(1-ql2.star.L1)+qy.star.L11*ql2.star.L1- (qy.star.L00*(1-ql2.star.L0)+qy.star.L10*ql2.star.L0)
Cl2g = as.numeric(dat$A == abar)/(gA_pred)

qL1_fit <- SuperLearner(Y = dat$L1, X = dat[, c("A", "W")],
                        family = binomial(), SL.library = SL.library)
L <- dat
L$A <- abar

ql1.fit = predict(qL1_fit, type="link", newdata=L)$pred
ql1.fit.update = glm(L1~offset(ql1.fit)+C_Q, weights = Cl2g, family= "quasibinomial", data=dat)
ql1.star = predict(ql1.fit.update, type="response")

ori_estimate=mean(qy.star.L11*ql2.star.L1*ql1.star+qy.star.L10*ql2.star.L0*(1-ql1.star)+qy.star.L01*(1-ql2.star.L1)*ql1.star+ qy.star.L00*(1-ql2.star.L0)*(1-ql1.star))
