
library(boot)
library(tidyverse)

source("../sample_code/R/ICAD-AIPW.R")
source("../sample_code/R/ICAD.R")
source("../sample_code/R/ICAD-km.R")

simple_rand <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian") {
  d <- data.frame(Y, A, Strata, W)
  d1 <- data.frame(Y, A = 1, Strata, W)
  d0 <- data.frame(Y, A = 0, Strata, W)
  M <- !is.na(Y)
  propensity.fit <- glm(M ~ ., data = data.frame(M, A, Strata, W), family = "binomial")
  propensity_score <- predict(propensity.fit, type = "response")
  glm.fit <- glm(Y~ ., data = d, family = family, weights = 1/propensity_score)
  X <- model.matrix(propensity.fit)
  p <- ncol(X)
  n <- nrow(d)
  pA <- predict(glm.fit, d, type = "response")
  p1 <- predict(glm.fit, d1, type = "response")
  p0 <- predict(glm.fit, d0, type = "response")
  est <- mean(p1) - mean(p0)
  if(family == "gaussian"){
    h_beta_A <- h_beta_0 <- h_beta_1 <- X
    h_beta_1[,"A"] <- 1
    h_beta_0[,"A"] <- 0
  }else if(family == "binomial"){
    h_beta_A <- h_beta_0 <- h_beta_1 <- (pA - pA^2) %*% t(rep(1, p)) * X
    h_beta_1[,"A"] <- pA - pA^2
    h_beta_0[,"A"] <- 0
  }
  e_beta_A <- (propensity_score - propensity_score^2) %*% t(rep(1, p)) * X
  MY <- ifelse(is.na(Y), 0, Y)
  c_1 <- colMeans((A-pi)/pi/(1-pi) * (M/propensity_score - 1) * h_beta_A) %*%
    solve(t(h_beta_A) %*% ((M/propensity_score) %*% t(rep(1, p)) * X)/n)
  c_2 <- colMeans(h_beta_1 - h_beta_0) %*% 
    solve(t(h_beta_A) %*% ((M/propensity_score) %*% t(rep(1, p)) * X)/n) %*%
    (t(e_beta_A) %*% (((M * (MY - pA)/propensity_score^2) %*% t(rep(1, p))) * X)/n) %*%
    solve(t(e_beta_A) %*% X/n)
  IF <- ((A-pi)/pi/(1-pi) - X %*% t(c_1)) * (M/propensity_score * (MY - pA)) - X %*% t(c_2) * (M - propensity_score) + p1 - p0 - est
  return(c(est, var(IF)/N))
}

N = 300
cov_N = 2

x1 = c()
x2 = c()
x3 = c()

for(xxx in 1:8) {
  
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
  
  alpha_0 = -3
  alpha_W = c(-0.7, -0.7, 0, -1.5)
  alpha_A = 0.5
  
  beta_0 = 1
  beta_W = c(-0.2, 0.8, 0, -1.0)
  beta_A = 2.5
  
  W_cross = W[1,] * S
  
  W = W %>% t()
  
  Mprob_0 = inv.logit(alpha_0 + alpha_W %*% rbind(W %>% t(), W_cross, S) + alpha_A * 0)
  Mprob_1 = inv.logit(alpha_0 + alpha_W %*% rbind(W %>% t(), W_cross, S) + alpha_A * 1)
  
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
    M_0 = c()
    M_1 = c()
    
    for (i in 1:N) {
      M_0[i] = rbinom(1, 1, Mprob_0[i])
      M_1[i] = rbinom(1, 1, Mprob_1[i])
      Y_0[i] = rbinom(1, 1, prob_0[i])
      Y_1[i] = rbinom(1, 1, prob_1[i])
    }
    rl[a] = Y_1 %>% mean() - Y_0 %>% mean()
    
    A = c()
    strata_0 = sum(S == 0)
    strata_1 = sum(S == 1)
    A[which(S == 0)] = rep_len(block, strata_0)
    A[which(S == 1)] = rep_len(block, strata_1)
    
    M = c()
    M[which(A == 0)] = M_0[which(A == 0)]
    M[which(A == 1)] = M_1[which(A == 1)]
    
    Y = c()
    Y[which(A == 0)] = Y_0[which(A == 0)]
    Y[which(A == 1)] = Y_1[which(A == 1)]
    Y[which(M == 1)] = NA
    
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

vr = (x3 - x2)/x3
