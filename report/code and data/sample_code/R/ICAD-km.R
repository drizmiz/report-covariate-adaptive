#' Inference under Covariate-Adaptive Design (ICAD) with time-to-event outcomes using the Kaplan-Meier estimator
#'
#' @param U time of event, that is, min(Y, M)
#' @param delta censoring indicator, that is I(Y <= M)
#' @param strata stratification variable (factor)
#' @param pi target proportion of people getting treatment
#' @param tau restriction time
#'
#' @return a list with survival functions (S), standard deviation using correct variance formula (std)
#' and standard deviation ingorning stratification (std_iid)
#' @export
#'
#' @examples
ICAD_km <- function(U, delta, strata, pi, tau){
  d <- data.frame(U, delta)
  dlambda <- map_dbl(1:tau, function(t){
    sum((d$U == t) & (d$delta == 1))/sum(d$U >= t)
  })
  S <- cumprod(1-dlambda)
  dR <- map(1:tau, function(t){
    (((d$U == t) & (d$delta == 1)) - (d$U >= t) * dlambda[t])/mean(d$U >= t)/(1- dlambda[t])
  }) %>% unlist %>% matrix(ncol = tau)
  R <- t(apply(dR, 1, cumsum))
  
  ER_S <- matrix(NA, nrow = length(unique(strata)), ncol = tau)
  for(t in 1:tau){
    for(s in 1:length(unique(strata))){
      ER_S[s,t] <- mean(R[strata == unique(strata)[s],t])
    }
  }
  p.S <- map_dbl(unique(strata), ~mean(strata == .))
  std <- sqrt(apply(R, 2, var)/pi - (1-pi)/pi * map_dbl(1:ncol(ER_S), ~sum(p.S * ER_S[,.]^2))) * S / sqrt(nrow(d))
  std_iid <- sqrt(apply(R, 2, var)/pi) * S / sqrt(nrow(d))
  # stdS <- sqrt(cumsum(dlambda/(1-dlambda)/map_dbl(1:tau, ~sum(d$U >= .))) * S^2 / pi)
  return(list(S = S, std = std, std_iid = std_iid))
}
