#' Log density of the normalized component of the Conway-Maxwell-Poisson
#'
#' @param y Value
#' @param mu Location parameter
#' @param nu Shape parameter
#'
#' @return the log of the normalized component
log_cmp <- function(y, mu, nu){nu*(y*log(mu) - lfactorial(y))}
