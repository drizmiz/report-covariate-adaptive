
library(tidyverse)
library(survival)
library(survminer)
library(cowplot)

source("../sample_code/R/ICAD-AIPW.R")
source("../sample_code/R/ICAD.R")
source("../sample_code/R/ICAD-km.R")

CTN44 <- read_csv("./CTN-0044-analysis/CTN44.csv")
CTN44$`0`[which(is.na(CTN44$`0`))] <- TRUE
CTN44$gender[which(is.na(CTN44$gender))] <- 1
pi <- 0.5

# average lab results -------------------
# handling missing data: missing if one missed more than 6 weeks

Y <- apply(CTN44[,7:30], 1, function(x){mean(x, na.rm = T)})
missing <- (is.na(CTN44[,seq(7, 30, by = 2)]) * is.na(CTN44[,seq(8, 30, by = 2)])) %>% apply(1,sum)
Y <- ifelse(missing > 6, NA, Y)
A <- as.numeric(CTN44$arm == "TES")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, gender, `0`)
I <- ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(3)

# survival analysis-------
# Y: time to first two consecutive negative lab result
# M: time to first missing visits
event_table <- matrix(NA, nrow = nrow(CTN44), ncol = 23)
for(j in 1:23){
  temp_event <- CTN44[,6+j] + CTN44[7+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN44[,7:30], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN44$arm == "TES")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, `0`)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 7)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 7)

# survival analysis: Kaplan-Meier estimator --------
# Do the following first: run trace("ggsurvplot_df", edit = T) and change
# all "lineytpe = "dashed"" (lines 108 and 110) to 
# "linetype = ifelse(df$strata == unique(df$strata)[1], "dashed", "dotted")"
event_table <- matrix(NA, nrow = nrow(CTN44), ncol = 23)
for(j in 1:23){
  temp_event <- CTN44[,6+j] + CTN44[7+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN44[,7:30], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN44$arm == "TES")
Strata <- as.factor(CTN44$strata)

# KM estimator of the treatment arm
fit_1 <- ICAD_km(E[A==1], C[A==1], Strata[A==1], pi, 14)
CTN44_1 <- data.frame(E = E[A==1], C = C[A==1])
fit_CTN44_1 <- survfit(Surv(E, C) ~ 1, data = CTN44_1)
fit_CTN44_1s <- fit_CTN44_1 
fit_CTN44_1$lower <- c(qnorm(0.025, fit_1$S, fit_1$std_iid), NA)
fit_CTN44_1$upper <- c(qnorm(0.975, fit_1$S, fit_1$std_iid), NA)
fit_CTN44_1s$lower <- c(qnorm(0.025, fit_1$S, fit_1$std), NA)
fit_CTN44_1s$upper <- c(qnorm(0.975, fit_1$S, fit_1$std), NA)

variance_reduction_1 <- round(1-(fit_1$std/fit_1$std_iid)^2,2)
p <- ggsurvplot(list(fit1 = fit_CTN44_1, fit2 =fit_CTN44_1s), data = CTN44_1, legend = "none",
                combine = T, conf.int = T, xlim = c(0,12), conf.int.style = "step",#  legend = c(0.8,0.9),
                #legend.title = "",  legend.labs = c("Stratification ignored", "Stratification not ignored"),
                break.x.by = 1, font.legend = 10, palette = c("black","black"),  conf.int.fill = 1,
                title = "NIDA-CTN-0044 Treatment group", 
                font.title = 8, font.xlab = 7, font.ylab = 7, font.tickslab = 6)
d_text <- data.frame(x = (1:12) + 0.2, y = 0, text = paste0(round(variance_reduction_1, 2) * 100, "%")[1:12])
p1 <- p$plot + geom_text(aes(x = x, y = y, label = text), size = 2.5, data = d_text) +
  annotate(geom="text", x=2.2, y=0.1, label="Variance Reduction", size = 3) + xlab("Visit")
ggsave("CTN44-treatment.tex", p1, width = 4.4, height = 2.2)

# KM estimator of the control arm
fit_0 <- ICAD_km(E[A==0], C[A==0], Strata[A==0], pi, 12)
CTN44_0 <- data.frame(E = E[A==0], C = C[A==0])
fit_CTN44_0 <- survfit(Surv(E, C) ~ 1, data = CTN44_0)
fit_CTN44_0s <- fit_CTN44_0
fit_CTN44_0$lower <- c(qnorm(0.025, fit_0$S, fit_0$std_iid), NA)
fit_CTN44_0$upper <- c(qnorm(0.975, fit_0$S, fit_0$std_iid), NA)
fit_CTN44_0s$lower <- c(qnorm(0.025, fit_0$S, fit_0$std), NA)
fit_CTN44_0s$upper <- c(qnorm(0.975, fit_0$S, fit_0$std), NA)

variance_reduction_0 <- round(1-(fit_0$std/fit_0$std_iid)^2,2)
p <- ggsurvplot(list(fit1 = fit_CTN44_0, fit2 =fit_CTN44_0s), data = CTN44_0, legend = "none",
                combine = T, conf.int = T, xlim = c(0,12), conf.int.style = "step", #legend = c(0.8,0.9),
                legend.title = "", legend.labs = c("Stratification ignored", "Stratification not ignored"),
                break.x.by = 1, font.legend = 10, palette = c("black","red"),
                title = "NIDA-CTN-0044 Control group", 
                font.title = 12, font.xlab = 10, font.ylab = 10)
d_text <- data.frame(x = (1:12) + 0.2, y = 0, text = paste0(round(variance_reduction_0, 2) * 100, "%")[1:12])
p0 <- p$plot + geom_text(aes(x = x, y = y, label = text), size = 3.5, data = d_text) +
  annotate(geom="text", x=2.2, y=0.1, label="Variance Reduction", size = 4) + xlab("Visit")
ggsave("CTN44-control.tex", p0, width = 6.6, height = 3.3)

pp <- plot_grid(p1, p0, ncol = 1, labels = c("A", "B"))
save_plot(filename = "CTN44-combined.tex",plot = pp, base_asp = 1)

# p1 <- ggdraw() + draw_image("CTN44-treatment.pdf")
# p0 <- ggdraw() + draw_image("CTN44-control.pdf")
