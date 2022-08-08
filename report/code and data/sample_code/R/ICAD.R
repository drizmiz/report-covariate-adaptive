# Inference under Covariate-Adaptive Design (ICAD) with continuous or binary outcomes
#'
#' @param Y outcome
#' @param A binary treatment indicator
#' @param Strata stratification variable (factor)
#' @param W baseline variables
#' @param pi target proportion of people getting treatment
#' @param family "gaussian" for continuous outcomes and "binomial" for binary outcomes
#'
#' @return the unadjusted, adjusted and DR-WLS estimators with variance and 95% confidence intervals.
#' @export
#'
#' @examples
ICAD <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian"){
  if(any(is.na(A))){
    stop("No missing treatment is allowed.")
  }
  if(any(is.na(Strata))){
    stop("Strata must be provided for asymptotic efficiency.")
  }
  if(!is.null(W) & any(is.na(W))){
    stop("Please impute missing covariates first.") 
  }
  if(!any(is.na(Y))){
    message("No missing outcome presents. The unadjusted, adjusted-for-strata and adjusted-for-all estimator are cacluated.")
    return(rbind(
      `unadjusted` = eCAD(Y, A, Strata, W, pi, family, method = "unadjusted"),
      `adjusted-for-strata` = eCAD(Y, A, Strata, W, pi, family, method = "adjusted-for-strata"),
      `adjusted-for-all` = eCAD(Y, A, Strata, W, pi, family, method = "adjusted-for-all")
    ))
  }else{
    message("There are missing outcomes. The unadjusted, adjusted-for-strata, adjusted-for-all and DR-WLS estimator are calculated.")
    return(rbind(
      unadjsed = eCAD(Y[!is.na(Y)], A[!is.na(Y)], Strata[!is.na(Y)], W[!is.na(Y),], pi, family, method = "unadjusted"),
      `adjusted-for-strata` = eCAD(Y[!is.na(Y)], A[!is.na(Y)], Strata[!is.na(Y)], W[!is.na(Y),], pi, family, method = "adjusted-for-strata"),
      `adjusted-for-all` = eCAD(Y[!is.na(Y)], A[!is.na(Y)], Strata[!is.na(Y)], W[!is.na(Y),], pi, family, method = "adjusted-for-all"),
      drwls = eCAD(Y, A, Strata, W, pi, family, method = "DR-WLS")
    )) 
  }
}

eCAD <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian", method){
  n <- length(Y)
  if(method == "unadjusted"){
    est <- mean(Y[A==1]) - mean(Y[A==0])
    IF <- (A-pi)/pi/(1-pi) * (Y - A * mean(Y[A==1]) - (1-A) * mean(Y[A==0]))
  }
  if(method == "adjusted-for-strata"){
    d <- data.frame(Y, A, Strata)
    d1 <- data.frame(Y, A = 1, Strata)
    d0 <- data.frame(Y, A = 0, Strata)
    glm.fit <- glm(Y~ ., data = d, family = family)
    p1 <- predict(glm.fit, d1, type = "response")
    p0 <- predict(glm.fit, d0, type = "response")
    est <- mean(p1) - mean(p0)
    IF <- (A-pi)/pi/(1-pi) * (Y - A * p1  - (1-A) * p0) + p1 - p0 - est
  }
  if(method == "adjusted-for-all"){
    d <- data.frame(Y, A, Strata, W)
    d1 <- data.frame(Y, A = 1, Strata, W)
    d0 <- data.frame(Y, A = 0, Strata, W)
    glm.fit <- glm(Y~ ., data = d, family = family)
    p1 <- predict(glm.fit, d1, type = "response")
    p0 <- predict(glm.fit, d0, type = "response")
    est <- mean(p1) - mean(p0)
    IF <- (A-pi)/pi/(1-pi) * (Y - A * p1  - (1-A) * p0) + p1 - p0 - est
  }
  if(method == "DR-WLS"){
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
  }
  E.S <- map_dbl(unique(Strata), ~mean((A[Strata == .]-pi) * IF[Strata == .]))
  p.S <- map_dbl(unique(Strata), ~mean(Strata == .))
  var <- var(IF) - sum(p.S * E.S^2)/pi/(1-pi)
  return(c(est = est, var = var/n, 
           CI.lower = qnorm(0.025, mean = est, sd =  sqrt(var/n)),
           CI.upper = qnorm(0.975, mean = est, sd =  sqrt(var/n))))
}
