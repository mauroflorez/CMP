% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_cmp.R
\name{mcmc_cmp}
\alias{mcmc_cmp}
\title{MCMC Algorithm for Conway-Maxwell-Poisson Regression Model for Multivariate Correlated Count Data}
\usage{
mcmc_cmp(
  y,
  X,
  S = 10000,
  nburn = 5000,
  initial_beta,
  initial_gamma,
  initial_b,
  prior_mean_beta,
  prior_var_beta,
  prior_mean_gamma,
  prior_var_gamma,
  v_0,
  R_0,
  intercept = FALSE,
  scale_b,
  scale_beta,
  scale_gamma,
  scale_cov_b,
  scale_cov_beta,
  scale_cov_gamma,
  inc_burn = FALSE,
  re_chain = TRUE,
  way = 2,
  random_seed,
  ...
)
}
\arguments{
\item{y}{Matrix of observations}

\item{X}{Covariates list, each element is the design matrix for each column of y}

\item{S}{Number of MCMC samples to be drawn}

\item{nburn}{Number of MCMC samples to burn-in}

\item{initial_beta}{List with initial value of \eqn{beta} for each response}

\item{initial_gamma}{List with initial value of \eqn{gamma} for each response}

\item{initial_b}{Initital value of \eqn{b}.}

\item{prior_mean_beta}{Prior mean for \eqn{beta}. (Default zero vector)}

\item{prior_var_beta}{Prior covariance matrix for \eqn{beta} (Default \eqn{I})}

\item{prior_mean_gamma}{Prior mean for \eqn{beta}. (Default zero vector)}

\item{prior_var_gamma}{Prior covariance matrix for \eqn{gamma} (Default \eqn{I})}

\item{v_0}{Prior degrees of freedom of random effects}

\item{R_0}{Prior covariance matrix of random effects}

\item{intercept}{Logical value indicating whether include the intercept}

\item{scale_b}{Covariance matrix for RW proposals of the random effects (Default \eqn{I})}

\item{scale_beta}{List with initial values for the scale matrices of \eqn{beta} (Default \eqn{I})}

\item{scale_gamma}{List with initial values for the scale matrices of \eqn{gamma} (Default \eqn{I})}

\item{scale_cov_b}{Scale parameter for the RW of random effects. (Default \eqn{2.4/sqrt(2)})}

\item{scale_cov_beta}{Scale parameter for the covariance of the proposals.}

\item{scale_cov_gamma}{Scale parameter for the covariance of the proposals.}

\item{inc_burn}{logical: include burned samples in the return}

\item{re_chain}{logical: If the posterior samples for the r.e are include. False return just the mean}

\item{way}{How to calculate the MCMC updates, based on Chib (2001)}

\item{random_seed}{Random seed}

\item{...}{Additional parameters of the MCMC algorithm}
}
\value{
A list:
\item{posterior_b}{List with posterior values of the random effects}
\item{estimation_beta}{Estimation of beta parameters}
\item{posterior_beta}{List with posterior values of beta}
\item{estimation_gamma}{Estimation of gamma parameters}
\item{posterior_gamma}{List with posterior values of gamma}
\item{posterior_D}{Values of covariance matrix D}
\item{fitted_mu}{Posterior of location parameters for each response}
\item{fitted_nu}{Posterior of shape parameters for ecah response}
\item{accept_rate_b}{Acceptance rate of Random Effects}
\item{accept_rate_beta}{Acceptance rate of beta}
\item{accept_rate_gamma}{Acceptance rate of gamma}
\item{scale_beta}{Estimated Scale matrix for beta parameters}
\item{scale_gamma}{Estimated Scale matrix for gamma parameters}
\item{X}{List of covariates used}
\item{y}{Matrix of observed counts}
}
\description{
MCMC Algorithm to estimate the parameters in the regression model for multivariate correlated count data
}
\examples{
\donttest{
  n = 50; J = 2
  X = list(matrix(rnorm(3*n), ncol = 3), matrix(rnorm(3*n), ncol = 3))
  beta <- list(c(1,0.1, 1), c(0, 0.5, -0.5))
  mu <- exp(prod_list(X, beta))
  y = matrix(rpois(n = length(mu), lambda = mu), nrow = n)
  fit <- mcmc_cmp(y, X, S = 10000, nburn = 1000, scale_cov_b = 0.8,
  scale_cov_beta = 0.04, scale_cov_gamma = 0.06)
}
}
