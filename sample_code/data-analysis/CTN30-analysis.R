library(tidyverse)
library(survival)
library(survminer)
# import data and covariate imputation
setwd("~/Dropbox/research/Clinical-trial/covariate-adaptive/Covariate-adaptive/")
source("R/ICAD-AIPW.R")
source("R/ICAD.R")
source("R/ICAD-km.R")
CTN30 <- readRDS("data/CTN30.rds")
CTN30$age[which(is.na(CTN30$age))] <- median(CTN30$age, na.rm = T)
CTN30$BL[which(is.na(CTN30$BL))] <- FALSE
CTN30 <- mutate(CTN30, strata = interaction(CTN30$heroin_history, CTN30$chronic_pain))
pi <- 0.5

# average lab results -------------------
# handling missing data: missing if one missed more than 2 weeks
Y <- apply(CTN30[,8:12], 1, function(x){mean(x, na.rm = T)})
missing <- (is.na(CTN30[,8]) & is.na(CTN30[,9])) + is.na(CTN30[,10]) + is.na(CTN30[,11]) + is.na(CTN30[,12])
Y <- ifelse(missing > 2, NA, Y)
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(3)

# retention outcome -------------------------------
Y <- CTN30$complete1
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))
ICAD(Y,A, Strata, W, pi = pi, family = "binomial") %>% round(2)

# survival analysis for restricted mean survival time (RMST) using the AIPW esitmator-------
# Y: time to first two consecutive negative lab result
# M: time to first missing visit
event_table <- matrix(NA, nrow = nrow(CTN30), ncol = 4)
for(j in 1:4){
  temp_event <- CTN30[,7+j] + CTN30[8+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN30[,8:12], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 5)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 5)

# survival analysis: Kaplan-Meier estimator --------
# Y: time to first two consecutive negative lab result
# M: time to first missing visit
event_table <- matrix(NA, nrow = nrow(CTN30), ncol = 4)
for(j in 1:4){
  temp_event <- CTN30[,7+j] + CTN30[8+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN30[,8:12], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)

# KM estimator of the treatment arm
fit_1 <- ICAD_km(E[A==1], C[A==1], Strata[A==1], pi, 5)
CTN30_1 <- data.frame(E = E[A==1], C = C[A==1])
fit_CTN30_1 <- survfit(Surv(E, C) ~ 1, data = CTN30_1)
fit_CTN30_1s <- fit_CTN30_1 
fit_CTN30_1$lower <- c(qnorm(0.025, fit_1$S, fit_1$std_iid), NA)
fit_CTN30_1$upper <- c(qnorm(0.975, fit_1$S, fit_1$std_iid), NA)
fit_CTN30_1s$lower <- c(qnorm(0.025, fit_1$S, fit_1$std), NA)
fit_CTN30_1s$upper <- c(qnorm(0.975, fit_1$S, fit_1$std), NA)

variance_reduction_1 <- round(1-(fit_1$std/fit_1$std_iid)^2,2)
p <- ggsurvplot(list(fit1 = fit_CTN30_1, fit2 =fit_CTN30_1s), data = CTN30_1, 
           combine = T, conf.int = T, xlim = c(0,5), conf.int.style = "step",
           legend.title = "", legend.labs = c("Stratification Ignored", "Stratification not Ignored"),
           break.x.by = 1, font.legend = 14, palette = c("blue","black"),
           title = "NIDA-CTN-0030 Treatment group")
ggsave("data-analysis/CTN30-treatment.png", print(p), width = 10, height = 7)

# KM estimator of the control arm
fit_0 <- ICAD_km(E[A==0], C[A==0], Strata[A==0], pi, 5)
CTN30_0 <- data.frame(E = E[A==0], C = C[A==0])
fit_CTN30_0 <- survfit(Surv(E, C) ~ 1, data = CTN30_0)
fit_CTN30_0s <- fit_CTN30_0
fit_CTN30_0$lower <- c(qnorm(0.025, fit_0$S, fit_0$std_iid), NA)
fit_CTN30_0$upper <- c(qnorm(0.975, fit_0$S, fit_0$std_iid), NA)
fit_CTN30_0s$lower <- c(qnorm(0.025, fit_0$S, fit_0$std), NA)
fit_CTN30_0s$upper <- c(qnorm(0.975, fit_0$S, fit_0$std), NA)

variance_reduction_0 <- round(1-(fit_0$std/fit_0$std_iid)^2,2)
p <- ggsurvplot(list(fit1 = fit_CTN30_0, fit2 =fit_CTN30_0s), data = CTN30_0, 
           combine = T, conf.int = T, xlim = c(0,5), conf.int.style = "step",
           legend.title = "", legend.labs = c("Stratification Ignored", "Stratification not Ignored"),
           break.x.by = 1, font.legend = 14, palette = c("blue","black"),
           title = "NIDA-CTN-0030 Control group")
ggsave("CTN30-control.png", print(p), width = 10, height = 7)
