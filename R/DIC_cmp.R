#' DIC of the regression model
#'
#' @description This function used an approach similar to the presented by Benson & Friel (2021) to calculate the BIC. We select S a sample size of the posterior samples to speed up computation
#'
#' @param fit An object from the mcmc_cmp_mh
#' @param S Number of iterations used to calculate the DIC
#'
#' @return Vector of approximated DIC
#' @export
#'
DIC_cmp <- function(fit, S = 100){
  X <- fit$X
  y <- fit$y
  l <- ncol(fit$posterior_b[[1]])
  n <- nrow(y)
  J <- ncol(y)

  selected_S <- sample(1:l, size = S)
  ll1 <- matrix(nrow = n, ncol = S)
  ll2 <- matrix(nrow = n, ncol = S)

  llk <- list()

  for(j in 1:J){
    llk[[j]] <- matrix(nrow = n, ncol = S)
    for(s in 1:S){
      mu <- fit$fitted_mu[[j]][,s]
      nu <- fit$fitted_nu[[j]][,s]

      llk[[j]][,s] <- llk_cmp(y[,j], mu, nu)
    }
  }

  post_b <- fit$posterior_b

  b_est <- purrr::map(post_b, rowMeans)
  beta_est <- fit$estimation_beta
  gamma_est <- fit$estimation_gamma

  Xtbeta <- purrr::map2(beta_est, X, function(a,b) c(b%*%a))
  Xtgamma <- purrr::map2(gamma_est, X, function(a,b) c(b%*%a))

  mu_est <- lapply(purrr::map2(Xtbeta, b_est, `+`), exp)
  nu_est <- lapply(Xtgamma, exp)

  llk_esp <- list()

  for(j in 1:J){
    mu_esp <- mu_est[[j]]
    nu_esp <- nu_est[[j]]
    llk_esp[[j]] <- llk_cmp(y[,j], mu_esp, nu_esp)
  }

  d_bar <- do.call(cbind, llk_esp)

  d_bar_post <- lapply(llk, rowMeans)

  p_dic <- 2*(d_bar - do.call(cbind, d_bar_post))

  DIC_cmp <- colSums(-2*d_bar + p_dic)

  return(DIC_cmp)
}
