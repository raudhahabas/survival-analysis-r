# =============================================================================
# Mortality Modelling in Insurance: survival analysis in R
#
# Sections:
#   1. Log-rank test comparing two groups           (Question 2c)
#   2. Maximum likelihood fit of a parametric hazard (Question 3b)
#   3. Graphical check of the parametric model       (Question 3c)
#   4. Cox proportional hazards model                (Question 4b-4e)
#
# Data: built-in datasets from the KMsurv and survival packages.
# =============================================================================

library(survival)
library(KMsurv)
library(ggplot2)


# -----------------------------------------------------------------------------
# 1. Log-rank test: do Group a and Group b have the same hazard?
# -----------------------------------------------------------------------------
times = c(10, 15, 20, 20, 25, 30, 18, 22, 25, 35, 40, 45)
status = c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1)  # 1 = death, 0 = censored
group = factor(c(rep("a", 6), rep("b", 6)))

survobj = Surv(times, status)
survdiff(survobj ~ group)
survdiff(survobj ~ group)$chisq  # 1.99 < 3.84, so no evidence of a difference


# -----------------------------------------------------------------------------
# 2. MLE for the hazard mu(t) = 2*theta*exp(theta*t) / (1 + exp(theta*t))
#    Channing House males, left-truncated at entry age, right-censored.
# -----------------------------------------------------------------------------
data("channing")
mcdata = subset(channing, gender == 1 & age > ageentry)

loglk = function(theta){
  ti = mcdata$age # observed age (exit age)
  ai = mcdata$ageentry # entry age (truncation)
  delta = mcdata$death # death indicator

  log_mu = log(2*theta) + theta*ti - log(1+exp(theta*ti))

  # A(ti) - A(ai): cumulative hazard between entry and exit
  cum_haz_diff = 2*log((1+exp(theta*ti))/(1+exp(theta*ai)))

  loglkvalue = sum(delta*log_mu - cum_haz_diff)
  return(loglkvalue)
}

optimize(loglk, interval = c(1e-8, 0.1), maximum = TRUE)  # theta_hat = 0.00334


# -----------------------------------------------------------------------------
# 3. Graphical check: if the model is right, log(2*exp(A(t)/2) - 1) = theta*t,
#    so the transformed Nelson-Aalen estimate should be linear in t.
# -----------------------------------------------------------------------------
# fit NA estimate with left-truncation
survobj = Surv(mcdata$ageentry, mcdata$age, mcdata$death)
fitobj = survfit(survobj ~ 1)

sumobj = summary(fitobj)
tau = sumobj$time
na_est = sumobj$cumhaz

# compute f(A) = log(...)
y_vals = log(2*exp(na_est/2) - 1)

# plot transformed NA estimate against tau
par(mar = c(5, 5.5, 4, 2))
plot(tau, y_vals,
     xlab = "t (age in months)",
     ylab = expression(log(2*e^{hat(A)(t)/2} - 1)),
     pch = 16)
abline(lm(y_vals ~ tau), col = "blue")


# -----------------------------------------------------------------------------
# 4. Cox proportional hazards model on the mgus2 dataset
# -----------------------------------------------------------------------------

# 4b. Fit the model and compare two patients
coxfit = coxph(Surv(futime, death) ~ age + sex + hgb + pstat, data = mgus2)
summary(coxfit)

# covariates
xA = c(age = 64, sexM = 1, hgb = 13.5, pstat = 1)
xB = c(age = 71, sexM = 0, hgb = 7.8, pstat = 1)

# estimate coefficients
beta_hat = coef(coxfit)

# hazard ratio of patient A relative to patient B
hr = exp(sum(beta_hat*(xA - xB)))
hr  # 0.479: patient A's risk of death is about 52% lower

# 4c. Breslow estimate: P(patient A dies in year 3 | survived 2 years)
patA = data.frame(age = 64, sex = "M", hgb = 13.5, pstat = 1)
survA = survfit(coxfit, newdata = patA, censor = FALSE)

# survival estimates at t = 24 and t = 36 months
S_24 = summary(survA, times = 24)$surv
S_36 = summary(survA, times = 36)$surv

prob = 1 - S_36/S_24
prob  # about 6.3%

# 4d. Likelihood ratio test: is pstat a significant covariate?
coxfit_full = coxph(Surv(futime, death) ~ age + sex + hgb + pstat, data = mgus2)
coxfit_red = coxph(Surv(futime, death) ~ age + sex + hgb, data = mgus2)

test_stat = 2*(coxfit_full$loglik[2] - coxfit_red$loglik[2])
test_stat  # 5.22 > 3.84, so pstat is significant at the 5% level

# 4e. Check the proportional hazards assumption with stratified models.
#     PH holds if the log cumulative hazard curves stay a constant distance apart.
cox_sex = coxph(Surv(futime, death) ~ age + hgb + pstat + strata(sex), data = mgus2)
cox_pstat = coxph(Surv(futime, death) ~ age + sex + hgb + strata(pstat), data = mgus2)

surv_sex = survfit(cox_sex)
surv_pstat = survfit(cox_pstat)

# stratified by sex
n1 = surv_sex$strata[1]
n2 = surv_sex$strata[2]

P_sex = ggplot(NULL) +
  geom_step(data = data.frame(time = surv_sex$time[1:n1],
                              lch = log(surv_sex$cumhaz[1:n1])),
            aes(x = time, y = lch), col = "red") +
  geom_step(data = data.frame(time = surv_sex$time[(n1+1):(n1+n2)],
                              lch = log(surv_sex$cumhaz[(n1+1):(n1+n2)])),
            aes(x = time, y = lch), col = "blue") +
  labs(x = "Time (months)", y = "log cumulative hazard",
       title = "Stratified by sex") +
  annotate("text", x = 380, y = c(1.0, 1.8),
           label = c("Female", "Male"), col = c("red", "blue"))

print(P_sex)

# stratified by pstat
m1 = surv_pstat$strata[1]
m2 = surv_pstat$strata[2]

P_pstat = ggplot(NULL) +
  geom_step(data = data.frame(time = surv_pstat$time[1:m1],
                              lch = log(surv_pstat$cumhaz[1:m1])),
            aes(x = time, y = lch), col = "black") +
  geom_step(data = data.frame(time = surv_pstat$time[(m1+1):(m1+m2)],
                              lch = log(surv_pstat$cumhaz[(m1+1):(m1+m2)])),
            aes(x = time, y = lch), col = "blue") +
  labs(x = "Time (months)", y = "log cumulative hazard",
       title = "Stratified by pstat") +
  annotate("text", x = 380, y = c(0.8, 1.6),
           label = c("pstat = 0", "pstat = 1"), col = c("black", "blue"))

print(P_pstat)
