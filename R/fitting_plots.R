#' Rootograms plots - Multivariate CMP
#'
#' @param fit An element from `mcmc_cmp`
#' @param S Optional. Indicates the number of posterior samples used (Default 1000)
#'
#' @return A rootogram plot for each response
#' @export
#' @import bayesplot
#'
#' @examples
#' \dontrun{
#'   n = 50; J = 2
#'   X = list(matrix(rnorm(3*n), ncol = 3), matrix(rnorm(3*n), ncol = 3))
#'   beta <- list(c(1,0.1, 1), c(0, 0.5, -0.5))
#'   mu <- exp(prod_list(X, beta))
#'   y = matrix(rpois(n = length(mu), lambda = mu), nrow = n)
#'   fit <- mcmc_cmp(y, X, S = 10000, nburn = 1000, scale_cov_b = 0.8,
#'   scale_cov_beta = 0.04, scale_cov_gamma = 0.06)
#'
#'   fitting_plot(fit)
#' }
fitting_plots <- function(fit, S = 1000){
  n <- nrow(fit$posterior_b[[1]])
  plots <- list()
  posterior_S <- sample(1:ncol(fit$posterior_b[[1]]), S)
  for(j in 1:length(fit$posterior_beta)){
    mu_est <- fit$fitted_mu[[j]][,posterior_S]
    nu_est <- fit$fitted_nu[[j]][,posterior_S]

    y_est <- mapply(com_sampler, mu_est, nu_est)
    y_est <- matrix(y_est, nrow = n)

    bayesplot::ppc_rootogram(fit$y[,j], t(y_est))
  }
}
