#' Rootograms plots - Multivariate CMP
#'
#' @param fit An element from `mcmc_cmp`
#' @param type Wheter to do a bar plot or a rootogram
#' @param S Optional. Indicates the number of posterior samples used (Default 100)
#' @export
#' @import bayesplot
#' @import ggplot2
#' @import cowplot
#'
#' @return No return value, called for plotting only
#'
#' @examples
#' \donttest{
#'   n = 50; J = 2
#'   X = list(matrix(rnorm(3*n), ncol = 3), matrix(rnorm(3*n), ncol = 3))
#'   beta <- list(c(1,0.1, 1), c(0, 0.5, -0.5))
#'   mu <- exp(prod_list(X, beta))
#'   y = matrix(rpois(n = length(mu), lambda = mu), nrow = n)
#'   fit <- mcmc_cmp(y, X, S = 1000, nburn = 1000, scale_cov_b = 0.8,
#'   scale_cov_beta = 0.04, scale_cov_gamma = 0.06)
#'   fitting_plots(fit)
#' }
fitting_plots <- function(fit, type = "rootogram", S = 100){
  n <- nrow(fit$posterior_b[[1]])
  plots <- list()
  posterior_S <- sample(1:ncol(fit$posterior_b[[1]]), S)
  J <- length(fit$posterior_beta)
  for(j in 1:J){
    mu_est <- fit$fitted_mu[[j]][,posterior_S]
    nu_est <- fit$fitted_nu[[j]][,posterior_S]

    y_est <- mapply(com_sampler, mu_est, nu_est)
    y_est <- matrix(y_est, nrow = n)

    if(type == "rootogram") plot(bayesplot::ppc_rootogram(fit$y[,j], t(y_est)) +
                                   labs(title = paste0("Response Variable #: ", j)))
    if(type == "bar") plot(bayesplot::ppc_bars(fit$y[,j], t(y_est)) +
                             labs(title = paste0("Response Variable #: ", j)))
  }
}


