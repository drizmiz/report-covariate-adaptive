
library(boot)
library(tidyverse)

source("../sample_code/R/ICAD-AIPW.R")
source("../sample_code/R/ICAD.R")
source("../sample_code/R/ICAD-km.R")

simple_rand <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian") {
  d <- data.frame(Y, A, Strata, W)
  d1 <- data.frame(Y, A = 1, Strata, W)
  d0 <- data.frame(Y, A = 0, Strata, W)
  glm.fit <- glm(Y ~ ., data = d, family = family)
  p1 <- predict(glm.fit, d1, type = "response")
  p0 <- predict(glm.fit, d0, type = "response")
  est <- mean(p1) - mean(p0)
  IF <- (A-pi)/pi/(1-pi) * (Y - A * p1  - (1-A) * p0) + p1 - p0 - est
  return(c(est, var(IF)/N))
}

N = 300
cov_N = 2

for(xxx in 1:10) {

S = rbinom(N, 1, 0.5)

W = matrix(nrow = cov_N, ncol = N)

for (j in 1:N) {
  for (i in 1:cov_N) {
    W[i, j] = rnorm(1)
  }
}

pi = 0.2

bsize = 10
block = rep(0, bsize)
block[sample(1:bsize, bsize * pi)] = 1

beta_0 = 1
beta_W = c(-0.2, 0.8, -0.5, -1.0)
beta_A = 2.5

W_cross = W[1,] * S

W = W %>% t()

prob_0 = inv.logit(beta_0 + beta_W %*% rbind(W %>% t(), W_cross, S) + beta_A * 0)
prob_1 = inv.logit(beta_0 + beta_W %*% rbind(W %>% t(), W_cross, S) + beta_A * 1)

est0 = c()
est1 = c()
est_ua = c()
var0 = c()
var1 = c()
var_ua = c()
rl = c()

B = 500
for (a in 1:B) {
  Y_0 = c()
  Y_1 = c()
  
  for (i in 1:N) {
    Y_0[i] = rbinom(1, 1, prob_0[i])
    Y_1[i] = rbinom(1, 1, prob_1[i])
  }
  rl[a] = Y_1 %>% mean() - Y_0 %>% mean()
  
  A = c()
  strata_0 = sum(S == 0)
  strata_1 = sum(S == 1)
  A[which(S == 0)] = rep_len(block, strata_0)
  A[which(S == 1)] = rep_len(block, strata_1)
  
  Y = c()
  Y[which(A == 0)] = Y_0[which(A == 0)]
  Y[which(A == 1)] = Y_1[which(A == 1)]
  
  ic = ICAD(Y, A, S, W, pi = pi, family = "binomial")
  est0[a] = ic[3, 1]
  var0[a] = ic[3, 2]
  est_ua[a] = ic[1, 1]
  var_ua[a] = ic[1, 2]
  
  A = rbinom(N, 1, pi)
  
  Y[which(A == 0)] = Y_0[which(A == 0)]
  Y[which(A == 1)] = Y_1[which(A == 1)]
  
  sr = simple_rand(Y, A, S, W, pi = pi, family = "binomial")
  est1[a] = sr[1]
  var1[a] = sr[2]
}

mean(var_ua)
mean(var0)
mean(var1)
x1 = x1 %>% append(var(est_ua))
x2 = x2 %>% append(var(est0))
x3 = x3 %>% append(var(est1))
}
