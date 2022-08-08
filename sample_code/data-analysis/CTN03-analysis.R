library(tidyverse)
# import data and imputation for covariates
setwd("~/Dropbox/research/Clinical-trial/covariate-adaptive/Covariate-adaptive/")
source("R/ICAD.R")
CTN03 <- readRDS("data/CTN03.rds")
CTN03$baseline.OPI[which(is.na(CTN03$baseline.OPI))] <- median(CTN03$baseline.OPI, na.rm = T)
CTN03$COWS[which(is.na(CTN03$COWS))] <- median(CTN03$COWS, na.rm = T)
CTN03$ARSW[which(is.na(CTN03$ARSW))] <- median(CTN03$ARSW, na.rm = T)
CTN03$VAS[which(is.na(CTN03$VAS))] <- median(CTN03$VAS, na.rm = T)
CTN03 <- CTN03[!is.na(CTN03$strata),]
CTN03$strata <- as.factor(CTN03$strata)
pi <- 0.5

# lab outcome -------------------------------
Y <- CTN03$outcome.OPI
A <- as.numeric(CTN03$arm == "7-day taper")
Strata <- CTN03$strata
W <- select(CTN03, sex, baseline.OPI, COWS, ARSW, VAS)
n <- nrow(CTN03)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

# inference
ICAD(Y,A, Strata, W, pi = pi, family = "binomial") %>% round(3)


# retention outcome -------------------------------
Y <- CTN03$complete
A <- as.numeric(CTN03$arm == "28-day taper")
Strata <- CTN03$strata
W <- select(CTN03, sex, baseline.OPI, COWS, ARSW, VAS)
n <- nrow(CTN03)
missing_proportion <- 1 - sum(is.na(Y))/length(Y)
n_strata <- length(unique(CTN03$strata))

# inference
ICAD(Y,A, Strata, W, pi = pi, family = "binomial") %>% round(3)
