#' Log likelihood of the Conway-Maxwell-Poisson Distribution
#'
#' @description This function calculates the log likelihood of the distribution as described by Benson and Friel (2021)
#'
#' @param y Count value
#' @param mu Location parameter
#' @param nu Shape parameter
#' @param r Number of acceptances
#'
#' @return Estimation of the log likelihood of the distribution
#' @export
#'
#' @examples llk_cmp(10, 5, 2)
llk_cmp <- function(y, mu, nu, r = 1000){
  #Estimation of the log-likelihood based on Benson & Friel
  cmp_sam <- mapply(com_sampler, mu, nu, n = r, ndraws = TRUE)
  n_r <- unlist(cmp_sam[2,])
  log_Bf <- unlist(cmp_sam[3,])
  log_Z <- rep(0, length(nu))
  log_Z[nu >= 1] <- mu[nu >= 1]
  log_qf <- log_cmp(y, mu, nu)
  log_M <- log(n_r/r)

  loglik <- unname(c(log_qf - log_Z + log_M - log_Bf))

  return(loglik)
}
