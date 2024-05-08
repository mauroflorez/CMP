#' Rejection Sampler - COM-Poisson
#'
#' @description Sampler for the Conway-Maxwell-Poisson as described in Algorithm 2 - Benson & Friel (2021)
#' @param mu Location parameter
#' @param nu Shape parameter
#' @param n Number of draws (default = 1)
#' @param ndraws Optional: Return the number of draws required to generate the n samples.
#'
#' @return A list or numeric if ndraws = FALSE:
#' \item{sample}{Values sampled from the distribution}
#' \item{drawsa}{Number of draws required in the rejection sampler}
#' \item{log_Bf}{Log of the boundary of the rejection sampler}
#' @export
#'
#' @examples com_sampler(2, 0.2, n = 10, ndraws = TRUE)
#' com_sampler(1, 2)
com_sampler = function(mu, nu, n = 1, ndraws = FALSE){
  sample <- c()
  draws <- 0
  if(nu >= 1){
    #Enveloping Bound
    #B_fg = (mu^floor(mu)/factorial(floor(mu)))^(nu-1)
    log_B_fg = (nu-1)*(floor(mu)*log(mu)-lfactorial(floor(mu)))
    while(length(sample) < n){
      y_sam = stats::rpois(1,mu) #Proposal
      #alpha = (mu^y_sam/factorial(y_sam))^nu/(B_fg*(mu^(y_sam)/factorial(y_sam)))
      log_alpha = nu*(y_sam*log(mu)-lfactorial(y_sam))-(log_B_fg + y_sam*log(mu) - lfactorial(y_sam))
      if(runif(1) <= exp(log_alpha)){
        sample <- c(sample, y_sam)
      }
      draws <- draws + 1
    }
  } else{
    #Enveloping Bound
    p = (2*nu)/(2*nu*mu + 1 + nu)
    ratio = floor(mu/((1-p)^(1/nu)))
    #B_fg = (1/p)*(mu^(nu*ratio))/((1-p)^ratio*(factorial(ratio))^nu)
    log_B_fg = -log(p) + nu*ratio*log(mu) - (ratio*log(1-p) + nu*lfactorial(ratio))
    while(length(sample) < n){
      u0 = stats::runif(1)
      y_sam = floor(log(u0)/log(1-p)) #Proposal
      #alpha = ((mu^y_sam/factorial(y_sam))^nu)/(B_fg*(1-p)^(y_sam)*p)
      log_alpha = nu*(y_sam*log(mu) - lfactorial(y_sam)) -
        (log_B_fg + y_sam*log(1-p) + log(p))
      if(runif(1) <= exp(log_alpha)){
        sample <- c(sample, y_sam)
      }
      draws <- draws + 1
    }
  }
  if(ndraws){
    return(list("samples" = sample, "draws" = draws, "log_Bf" = log_B_fg))
  }
  else{
    return(sample)
  }
}
