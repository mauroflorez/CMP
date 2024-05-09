#' MCMC Algorithm for Conway-Maxwell-Poisson Regression Model for Multivariate Correlated Count Data
#'
#' @description MCMC Algorithm to estimate the parameters in the regression model for multivariate correlated count data
#'
#' @param y Matrix of observations
#' @param X Covariates list, each element is the design matrix for each column of y
#' @param S Number of MCMC samples to be drawn
#' @param nburn Number of MCMC samples to burn-in
#' @param initial_beta List with initial value of \eqn{beta} for each response
#' @param initial_gamma List with initial value of \eqn{gamma} for each response
#' @param initial_b Initital value of \eqn{b}.
#' @param prior_mean_beta Prior mean for \eqn{beta}. (Default zero vector)
#' @param prior_var_beta Prior covariance matrix for \eqn{beta} (Default \eqn{I})
#' @param prior_mean_gamma Prior mean for \eqn{beta}. (Default zero vector)
#' @param prior_var_gamma Prior covariance matrix for \eqn{gamma} (Default \eqn{I})
#' @param v_0 Prior degrees of freedom of random effects
#' @param R_0 Prior covariance matrix of random effects
#' @param intercept Logical value indicating whether include the intercept
#' @param scale_b Covariance matrix for RW proposals of the random effects (Default \eqn{I})
#' @param scale_beta List with initial values for the scale matrices of \eqn{beta} (Default \eqn{I})
#' @param scale_gamma List with initial values for the scale matrices of \eqn{gamma} (Default \eqn{I})
#' @param scale_cov_b Scale parameter for the RW of random effects. (Default \eqn{2.4/sqrt(2)})
#' @param scale_cov_beta Scale parameter for the covariance of the proposals.
#' @param scale_cov_gamma Scale parameter for the covariance of the proposals.
#' @param inc_burn logical: include burned samples in the return
#' @param re_chain logical: If the posterior samples for the r.e are include. False return just the mean
#' @param progress output of algorithm: "acc_rates" update each 10 interactions the acceptance ratios, "bar" show a bar of progress
#' @param way How to calculate the MCMC updates, based on Chib (2001)
#' @param random_seed Random seed
#' @param ... Additional parameters of the MCMC algorithm
#'
#' @return A list:
#' \item{posterior_b}{List with posterior values of the random effects}
#' \item{estimation_beta}{Estimation of beta parameters}
#' \item{posterior_beta}{List with posterior values of beta}
#' \item{estimation_gamma}{Estimation of gamma parameters}
#' \item{posterior_gamma}{List with posterior values of gamma}
#' \item{posterior_D}{Values of covariance matrix D}
#' \item{fitted_mu}{Posterior of location parameters for each response}
#' \item{fitted_nu}{Posterior of shape parameters for ecah response}
#' \item{accept_rate_b}{Acceptance rate of Random Effects}
#' \item{accept_rate_beta}{Acceptance rate of beta}
#' \item{accept_rate_gamma}{Acceptance rate of gamma}
#' \item{scale_beta}{Estimated Scale matrix for beta parameters}
#' \item{scale_gamma}{Estimated Scale matrix for gamma parameters}
#' \item{X}{List of covariates used}
#' \item{y}{Matrix of observed counts}
#' @export
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
#' }
mcmc_cmp <- function(y, X, S = 10000, nburn = 5000, initial_beta, initial_gamma, initial_b,
                     prior_mean_beta, prior_var_beta, prior_mean_gamma, prior_var_gamma,
                     v_0, R_0, intercept = FALSE, scale_b, scale_beta, scale_gamma,
                     scale_cov_b, scale_cov_beta, scale_cov_gamma,
                     inc_burn = FALSE, re_chain = TRUE, progress = "acc_rates", way = 2, random_seed,...){
  ### -------------------------- Check Values ---------------------------- ###

  #Set seed
  if(!missing(random_seed)) set.seed(random_seed)

  #Include intercept to each design_matrix if necessary
  design_matrix <- X
  if(intercept) design_matrix <- purrr::map(X, function(x) cbind(1, x))

  #Number of variables
  J <- ncol(y)
  #Number of observations
  n <- nrow(y)
  #Dimension of each design matrix
  k_j <- purrr::map_dbl(design_matrix, base::ncol)

  ###--------------------------- Default Values ---------------------------###

  if(missing(prior_mean_beta))
    prior_mean_beta <- purrr::map(seq_len(ncol(y)), function(j) rep(0, ncol(X[[j]])))

  if(missing(prior_var_beta))
    prior_var_beta <- purrr::map(seq_len(ncol(y)), function(j) 0.01*diag(ncol(X[[j]])))

  if(missing(prior_mean_gamma))
    prior_mean_gamma <- purrr::map(seq_len(ncol(y)), function(j) rep(0, ncol(X[[j]])))

  if(missing(prior_var_gamma))
    prior_var_gamma <- purrr::map(seq_len(ncol(y)), function(j) 0.01*diag(ncol(X[[j]])))

  if(missing(v_0))
    v_0 <- 2*J

  if(missing(R_0))
    R_0 <-  diag(J)

  #Random Effects
  if(missing(initial_b))
    initial_b <- matrix(c(rep(0, n*J)), nrow = n, ncol = J)

  # Extract initial values for beta and gamma
  if(missing(initial_beta) || missing(initial_gamma)){
    initial_beta <- if(missing(initial_beta)) purrr::map(k_j, function(x) rep(0, x)) else initial_beta
    initial_gamma <- if(missing(initial_gamma)) purrr::map(k_j, function(x) rep(0, x)) else initial_gamma
  }
  # Default scale parameter
  if(missing(scale_cov_b))
    scale_cov_b = 2.4/sqrt(J)
  # Jumps for RW
  if(missing(scale_b))
    scale_b <- scale_cov_b^2 * diag(J)

  ### -------------------- Robust Adaptive Metropolis (RAM) ----------------###

  #Default Values

  if(!missing(scale_beta)){
    if(!is.list(scale_beta)) {
      stop("scale_beta should be a list!")
    } else if(!all(purrr::map_lgl(scale_beta, is.matrix))){
      #If elements are vector or a scalar then pass to a diagonal matrix
      scale_beta = purrr::map2(scale_beta, k_j, function(x,y) diag(x, y))
      S_beta = purrr::map2(scale_beta, k_j, function(x,y) diag(x, y))
    } else{
      S_beta = scale_beta
    }
  }
  else {
    if(missing(scale_cov_beta)) scale_cov_beta <- purrr::map(k_j, function(x) 2.4/sqrt(x))
    S_beta = purrr::map2(scale_cov_beta, k_j, function(x,y) x^2*diag(y))
  }

  if(missing(scale_cov_gamma))
    scale_cov_beta <- purrr::map(k_j, function(x) 2.4/sqrt(x))

  if(!missing(scale_gamma)){
    if(!is.list(scale_gamma)) {
      stop("scale_gamma should be a list!")
    } else if(!all(purrr::map_lgl(scale_gamma, is.matrix))){
      #If elements are vector or a scalar then pass to a diagonal matrix
      scale_gamma = purrr::map2(scale_gamma, k_j, function(x,y) diag(x, y))
      S_gamma = purrr::map(scale_gamma, function(x) t(chol(x)))
    } else{
      S_gamma = scale_gamma
    }
  } else {
    if(missing(scale_cov_gamma)) scale_cov_gamma <- purrr::map(k_j, function(x) 2.4/sqrt(x))
    S_gamma = purrr::map2(scale_cov_gamma, k_j, function(x,y) x^2*diag(y))
  }

  ###------------------------- Auxiliary Functions -------------------------###

  #Sample MVN
  sample_mvn <- function(mean, cov) {
    drop(mvnfast::rmvn(1, mu = mean, sigma = cov))
  }

  #Evaluate MVN
  evaluate_mvn <- function(x, mean, cov){
    mvnfast::dmvn(X = x, mu = mean, sigma = cov)
  }
  ## Prior densities to evaluate in acceptance ratios
  prior_b <- function(x, D) mvnfast::dmvn(x, mu = rep(0,J), sigma = D)

  prior_beta <- function(beta){
    purrr::pmap_dbl(list(beta, prior_mean_beta, prior_var_beta), evaluate_mvn)
  }

  prior_gamma <- function(gamma){
    purrr::pmap_dbl(list(gamma, prior_mean_gamma, prior_var_gamma), evaluate_mvn)
  }

  #Unnormalized density COM-Poisson
  log_dcmp <- function(y, mu, nu) nu*(y*log(mu) - lfactorial(y))

  #Matrix to list
  matrix2list <- function(x){
    lapply(seq_len(ncol(x)), function(i) x[, i])
  }

  #Acceptance ratios
  accept_rule = function(y, mu_new, nu_new, mu, nu, log_prior, param){
    log_qf <- log_dcmp(y, mu, nu)  # Compute the log-density of the target distribution
    log_qf_star <- log_dcmp(y, mu_new, nu_new)  # Compute the log-density of the proposal distribution
    # Generate an auxiliary data set
    y_aux <- matrix(mapply(com_sampler, c(mu_new), c(nu_new)), nrow = n, byrow = F)
    # Compute the log-density of the target and proposal distributions for the auxiliary data set
    log_qf_aux <- log_dcmp(y_aux, mu, nu)
    log_qf_aux_star <- log_dcmp(y_aux, mu_new, nu_new)

    # Evaluate the specified version of the algorithm using a switch statement
    switch(param,
           r_eff = {
             # Use rowSums() to sum over the rows of the log-density matrices
             loglike <- log_prior + rowSums(log_qf_star - log_qf) + rowSums(log_qf_aux - log_qf_aux_star)
             accept <- pmin(1, exp(loglike))
           },
           way_1 = {
             # Use sum() to sum over all the elements of the log-density matrices
             loglike <- log_prior + sum(log_qf_star - log_qf) + sum(log_qf_aux - log_qf_aux_star)
             accept <-  min(1, exp(loglike))
           },
           way_2 = {
             # Use colSums() to sum over the columns of the log-density matrices
             loglike <- log_prior + colSums(log_qf_star - log_qf) + colSums(log_qf_aux - log_qf_aux_star)
             accept <- pmin(1, exp(loglike))
           },
    )
    return(accept)  # Return the acceptance probability of the proposal distribution
  }

  ### -------------------------- Initial Values --------------------------###

  D_current <- solve(drop(stats::rWishart(n = 1, df = v_0, Sigma = R_0)))
  b_current <- initial_b
  beta_current <- initial_beta
  gamma_current <- initial_gamma

  mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
  nu_current <- exp(prod_list(design_matrix, gamma_current))

  #Log - likelihood

  loglik_b <- log(apply(b_current, 1, prior_b, D = D_current))
  loglik_beta <- log(prior_beta(beta_current))
  loglik_gamma <- log(prior_gamma(gamma_current))

  ### ---------------------------- Storage ------------------------------###

  #Set values for the number of acceptances
  accept_beta <- rep(0, J)
  accept_gamma <- rep(0, J)
  accept_b <- rep(0, n)
  #Define matrices to save the posterior draws into

  post_b <- list()
  post_beta <- list()
  post_gamma <- list()
  post_D <- matrix(ncol = nburn+S, nrow = length(c(D_current)))

  for(j in 1:J){
    post_b[[j]] <- matrix(nrow = n, ncol = nburn+S)
    post_beta[[j]] <- matrix(nrow = nburn+S, ncol = k_j[j])
    post_gamma[[j]] <- matrix(nrow = nburn+S, ncol = k_j[j])

    post_b[[j]][,1] <- b_current[,j]
    post_beta[[j]][1,] <- beta_current[[j]]
    post_gamma[[j]][1,] <- gamma_current[[j]]
  }

  post_D[,1] = c(D_current)


  ##############################################################################
  ############################### MCMC - Model #################################

  for(s in 2:(S+nburn)){

    ############################################################################
    #                            Block 1: Sample b
    ############################################################################

    #Proposal
    b_star <- t(apply(b_current, 1, function(x)
      drop(mvnfast::rmvn(1, mu = x, sigma = scale_b))))
    #Proposal \mu and \nu
    mu_star <- exp(prod_list(design_matrix, beta_current) + b_star)
    nu_star <- nu_current
    #Evaluate the priors

    loglik_b_star <- apply(b_star, 1, function(x) mvnfast::dmvn(x, mu = rep(0,J),
                                                                sigma = D_current,
                                                                log = T))
    log_prior <- loglik_b_star - loglik_b
    #Which i<n accepted?
    accept_prob_b <- accept_rule(y = y,
                                 mu_new = mu_star, nu_new = nu_star,
                                 mu = mu_current, nu = nu_current,
                                 log_prior = log_prior, param = "r_eff")
    ind_accept <- which(stats::runif(nrow(y)) < accept_prob_b)

    #Accept new values
    b_current[ind_accept,] <- b_star[ind_accept,]
    #Update Loglik
    loglik_b[ind_accept] <- loglik_b_star[ind_accept]
    #Update values
    mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
    #Save Acceptance Ratios
    accept_b[ind_accept] <- accept_b[ind_accept] + 1

    ############################################################################
    #                            Block 3: Sample Beta
    ############################################################################

    #Proposals
    beta_star <- purrr::map2(beta_current, S_beta, sample_mvn)

    mu_star <- exp(prod_list(design_matrix, beta_star) + b_current)
    nu_star <- nu_current
    #Log Priors
    loglik_beta_star = log(prior_beta(beta_star))
    log_prior <- loglik_beta_star - loglik_beta
    #Which variable j<J accepted?
    accept_prob_beta <- accept_rule(y = y,
                                    mu_new = mu_star, nu_new = nu_star,
                                    mu = mu_current, nu = nu_current,
                                    log_prior = log_prior, param = way)

    var_accept <- which(stats::runif(J) < accept_prob_beta)
    #Accept
    beta_current[var_accept] <- beta_star[var_accept]
    loglik_beta[var_accept] <- loglik_beta_star[var_accept]
    #Update
    mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
    accept_beta[var_accept] <- accept_beta[var_accept] + 1

    ############################################################################
    #                         Block 4: Sample Gamma
    ############################################################################

    #Proposals
    gamma_star <- purrr::map2(gamma_current, S_gamma, sample_mvn)
    mu_star <- mu_current
    nu_star <- exp(prod_list(design_matrix, gamma_star))
    #Log Priors
    loglik_gamma_star <- log(prior_gamma(gamma_star))
    log_prior <- loglik_gamma_star - loglik_gamma
    #Which variable j<J accepted?
    accept_prob_gamma <- accept_rule(y = y,
                                     mu_new = mu_star, nu_new = nu_star,
                                     mu = mu_current, nu = nu_current,
                                     log_prior = log_prior, param = way)

    var_accept <- which(stats::runif(J) < accept_prob_gamma)
    #Accept
    gamma_current[var_accept] <- gamma_star[var_accept]
    loglik_gamma[var_accept] <- loglik_gamma_star[var_accept]
    #Update
    nu_current <- exp(prod_list(design_matrix, gamma_current))
    accept_gamma[var_accept] <- accept_gamma[var_accept] + 1

    ############################################################################
    #                         Block 5: Sample D
    ############################################################################

    scale <- solve(solve(R_0) + crossprod(b_current))
    #Sampling
    D_current <- solve(drop(rWishart(n = 1, df = n + v_0, Sigma = scale)))

    #############################  Save Simulations ############################

    for(j in 1:J){
      post_b[[j]][,s] <- b_current[,j]
      post_beta[[j]][s,] <- beta_current[[j]]
      post_gamma[[j]][s,] <- gamma_current[[j]]
    }
    post_D[,s] = c(D_current)

    ### ------------------------- Output ---------------------------- ###

    if(progress == "bar"){
      if(s == 2){
        pb <- progress::progress_bar$new(format = "Completing [:bar] :percent || ETA: :eta]",
                               total = S+nburn,    # Current bar character
                               clear = FALSE,    # If TRUE, clears the bar when finish
                               width = 80)
      } else {
        pb$tick()
        Sys.sleep(1 / 100)
      }
    } else if(progress == "acc_rates" && s%%100 == 0){
      cat("Progress: ", round(s/(nburn+S)*100,2), "% ---- Iteration: ", s,"/", nburn+S, "\n")
      cat("Beta Acceptance Ratio: ", accept_beta/s, "\n")
      cat("Gamma Acceptance Ratio: ", accept_gamma/s, "\n")
      cat("R.E Acceptance Ratio: ", mean(accept_b)/s, "\n")
      cat("#----------------------------------------------# \n")
    }

  }

  #Return
  if(inc_burn == FALSE){
    post_b <- purrr::map(post_b, ~ .x[, -seq_len(nburn)])
    post_beta <- purrr::map(post_beta, ~ .x[-seq_len(nburn), ])
    post_gamma <- purrr::map(post_gamma, ~ .x[-seq_len(nburn), ])
  }

  #Estimations
  est_b <- purrr::map(post_b, rowMeans)
  est_beta <- purrr::map(post_beta, colMeans)
  est_gamma <- purrr::map(post_gamma, colMeans)

  #Posterior fitted values

  Xtbeta <- purrr::map2(post_beta, X, function(a,b) b%*%t(a))
  Xtgamma <- purrr::map2(post_gamma, X, function(a,b) b%*%t(a))

  fitted_mu <- lapply(purrr::map2(Xtbeta, post_b, `+`), exp)
  fitted_nu <- lapply(Xtgamma, exp)

  if(re_chain == FALSE){
    post_b <- est_b #Just the mean of the R.E
  }

  return(list(posterior_b = post_b, estimation_beta = est_beta, posterior_beta = post_beta,
              estimation_gamma = est_gamma, posterior_gamma = post_gamma,
              posterior_D = post_D, fitted_mu = fitted_mu, fitted_nu = fitted_nu,
              accept_rate_b = accept_b/(nburn+S), accept_rate_beta = accept_beta/(nburn + S),
              accept_rate_gamma = accept_gamma/(nburn + S), Scale_beta = S_beta, Scale_gamma = S_gamma,
              X = design_matrix, y = y))
}
