
library(tidyverse)
library(survival)
library(survminer)

source("../sample_code/R/ICAD-AIPW.R")
source("../sample_code/R/ICAD.R")
source("../sample_code/R/ICAD-km.R")

setwd("./CTN-0030-analysis")

CTN30 <- read_csv("CTN30.csv")

CTN30$age[which(is.na(CTN30$age))] <- median(CTN30$age, na.rm = T)
CTN30$BL[which(is.na(CTN30$BL))] <- FALSE
CTN30 <- mutate(CTN30, strata = interaction(CTN30$heroin_history, CTN30$chronic_pain))
pi <- 0.5

# average lab results -------------------
# handling missing data: missing if one missed more than 2 weeks
Y <- apply(CTN30[,8:12], 1, function(x){mean(x, na.rm = T)})
missing <- (is.na(CTN30[,8]) & is.na(CTN30[,9])) + is.na(CTN30[,10]) + is.na(CTN30[,11]) + is.na(CTN30[,12])
missing <- missing %>% as.vector()
Y <- ifelse(missing > 2, NA, Y)
# 2 is Enhanced
A <- as.numeric(CTN30$arm1 == 2)
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

I <- ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(3)

setwd("..")
