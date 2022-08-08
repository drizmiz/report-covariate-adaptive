#' Inference under Covariate-Adaptive Design (ICAD) with time-to-event outcomes using the AIPW estimator
#'
#' @param E time of event, that is, min(Y, M)
#' @param C censoring indicator, that is I(Y <= M)
#' @param A treatment indicator
#' @param Strata stratification variable (factor)
#' @param W baseline variables
#' @param pi target proportion of people getting treatment
#' @param tau restriction time
#' 
#' @import tidyverse
#'
#' @return point estimate and its estimated variance
#' @export
#'
#' @examples
#' expit <- function(x) {1/(1+exp(-x))}
#' n <- 20000
#' tau <- 6
#' A <- rbinom(n, 1, 0.5)
#' W <- rbinom(n, 1, 0.5)
#' S_t <- matrix(NA, nrow = n, ncol = tau)
#' for(i in 1:ncol(S_t)){
#'    S_t[,i] <- (1 - expit(- 1 - A - W))^i
#' }
#' S_t <- cbind(S_t, rep(0, n))
#' Y <- apply(S_t, 1, function(s){which(s < runif(1))[1]})
#' M_t <- matrix(NA, nrow = n, ncol = tau)
#' for(i in 1:ncol(M_t)){
#'    M_t[,i] <- (1 - expit(- 1 - A - W))^i
#' }
#' M_t <- cbind(M_t, rep(0, n))
#' M <- apply(M_t, 1, function(s){which(s < runif(1))[1]})
#' E <- pmin(Y, M)
#' C <- Y <= M
#' ICAD_tte(E,C,A,W, tau = 6)
#' ## the true RMST for tau = 6 is 
#' sum((1-expit(-3))^(1:5) - (1-expit(-1))^(1:5))/2
#' 
ICAD_tte <- function(E, C, A, Strata, W = NULL, pi = 0.5, tau = NULL){
  n <- length(E)
  if(is.null(W)){
    covariates <- Strata
  }else{
    covariates <- cbind(Strata, W)
  }
  d <- data.frame(E, C, A, covariates)
  R <- map(1:(tau-1), ~(E >= .)) %>% unlist %>% matrix(nrow = n)
  X <- model.matrix(lm(E ~ . -C, d))
  d <- data.frame(E, C, X[,-1])
  hazard <- matrix(NA, nrow = n, ncol = tau - 1)
  hazard_1 <- matrix(NA, nrow = n, ncol = tau - 1)
  hazard_0 <- matrix(NA, nrow = n, ncol = tau - 1)
  censor <- matrix(NA, nrow = n, ncol = tau - 1)
  for(t in 1:(tau-1)){
    d_t <- d[E >= t, ]
    Y_t <- (d_t$E == t) & (d_t$C == 1)
    M_t <- (d_t$E == t) & (d_t$C == 0)
    fit_event_t <- glm(Y_t ~ . - E - C, data = cbind(d_t, Y_t), family = "binomial")
    if(length(fit_event_t$coefficients) > fit_event_t$rank){
      fit_event_t <- glm(Y_t ~ 0 + . - E - C, data = cbind(d_t, Y_t), family = "binomial")
    }
    hazard[,t] <- predict(fit_event_t, data.frame(A, X[,-(1:2)]), type = "response")
    hazard_1[,t] <- predict(fit_event_t, data.frame(A = 1, X[,-(1:2)]), type = "response")
    hazard_0[,t] <- predict(fit_event_t, data.frame(A = 0, X[,-(1:2)]), type = "response")
    fit_censor_t <- glm(M_t ~ . - E - C, data = cbind(d_t, M_t), family = "binomial")
    if(length(fit_censor_t$coefficients) > fit_censor_t$rank){
      fit_censor_t <- glm(M_t ~ 0 + . - E - C, data = cbind(d_t, M_t), family = "binomial")
    }
    censor[,t] <- predict(fit_censor_t, data.frame(A, X[,-(1:2)]), type = "response")
  }
  hazard <- bounded(hazard)
  hazard_1 <- bounded(hazard_1)
  hazard_0 <- bounded(hazard_0)
  censor <- bounded(censor)
  survival_t <- apply(1 - hazard, 1, cumprod) %>% t
  survival_t1 <- apply(1 - hazard_1, 1, cumprod) %>% t
  survival_t0 <- apply(1 - hazard_0, 1, cumprod) %>% t
  cumcensor_t <- apply(1 - censor, 1, cumprod) %>% t
  cumcensor_t <- cbind(1 - 1e-5, cumcensor_t)
  phi <- 0
  for(t in 1:(tau-1)){
    for(m in 1:t){
      phi <- phi - (A-pi)/(pi*(1-pi)) * R[,m] / cumcensor_t[,m] * survival_t[,t] / survival_t[,m] * (((E == m) & (C == 1)) - hazard[,m])
    }
    phi <- phi + survival_t1[,t] - survival_t0[,t]
  }
  
  # calculate c_1(t) and c_2(t)
  dev_phi_beta <- map(1:(tau-1), function(j){
    dev_phi_beta_t <- 0
    for(t in j:(tau-1)){
      if(j > 1){
        for(m in 1:(j-1)){
          dev_phi_beta_t <- dev_phi_beta_t + (A-pi)/(pi*(1-pi)) * R[,m] / cumcensor_t[,m] * hazard[,j] *
            survival_t[,t] / survival_t[,m] * (((E == m) & (C == 1)) - hazard[,m])
        }
      }
      dev_phi_beta_t <- dev_phi_beta_t + 
        (A-pi)/(pi*(1-pi)) * R[,j] / cumcensor_t[,j] * hazard[,j] * survival_t[,t] / survival_t[,j] +
        survival_t1[,t] * hazard_1[,j] - survival_t0[,t] * hazard_0[,j]
    }
    colMeans((dev_phi_beta_t %*% t(rep(1,ncol(X)))) * X)
  })
  dev_phi_alpha <- map(1:(tau-1), function(j){
    dev_phi_alpha_t <- 0
    for(t in j:(tau-1)){
      for(m in j:t){
        dev_phi_alpha_t <- dev_phi_alpha_t - (A-pi)/(pi*(1-pi)) * R[,m] / cumcensor_t[,m] * censor[,j] *
          survival_t[,t] / survival_t[,m] * (((E == m) & (C == 1)) - hazard[,m])
      }
    }
    colMeans((dev_phi_alpha_t %*% t(rep(1,ncol(X)))) * X)
  })
  c_1 <- map(1:(tau-1), ~dev_phi_beta[[.]] %*% solve(t((R[,.] * hazard[,.] * (1-hazard[,.])) %*% t(rep(1, ncol(X))) * X) %*% X))
  c_2 <- map(1:(tau-1), ~dev_phi_alpha[[.]] %*% solve(t((R[,.] * censor[,.] * (1-censor[,.])) %*% t(rep(1, ncol(X))) * X) %*% X))
  
  # construct the AIPW estimator
  est <- phi
  for(t in 1:(tau-1)){
    est <- est - (X %*% t(c_1[[t]])) * R[,t] * (((E == t) & (C == 1)) - hazard[,t]) - 
      (X %*% t(c_2[[t]])) * R[,t] * (((E == t) & (C == 0)) - censor[,t])
  }
  
  # calculate variance
  tilde_V <- var(est)/n
  temp <- 0
  for(t in 1:(tau-1)){
    for(m in 1:t){
      temp <- temp + (A - pi)/pi/(1-pi) * R[,m] / cumcensor_t[,m] * survival_t[,t] / survival_t[,m] * (((E == m) & (C == 1)) - hazard[,m])
    }
  }
  E.S <- map_dbl(unique(Strata), ~mean(temp[Strata == .]))
  p.S <- map_dbl(unique(Strata), ~mean(Strata == .))
  V <- tilde_V - (1-2*pi)^2/pi/(1-pi) * sum(p.S * E.S^2)
  
  return(c(est = mean(est), var = V))
}


bounded <- function(x, r = 1e-2){
  xx <- x
  xx[x < r] <- r
  xx[x > 1-r] <- 1-r
  return(xx)
}
