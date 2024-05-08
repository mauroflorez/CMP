#' DIC of the regression model
#'
#' @param fit An object from the mcmc_cmp_mh
#' @param X Covariates list, each element is the design matrix for each column of y
#' @param y Matrix of observations
#' @param S Number of iterations used to calculate the DIC
#'
#' @return Approximated DIC
#' @export
#'
DIC_cmp <- function(fit, X, y, S = 100){
  l <- length(fit$posterior_b)
  n <- nrow(y)
  selected_S <- sample(1:l, size = S)
  ll1 <- matrix(nrow = n, ncol = S)
  ll2 <- matrix(nrow = n, ncol = S)

  X1 <- X[[1]]
  X2 <- X[[2]]

  y1 <- y[,1]
  y2 <- y[,2]

  for(s in 1:S){
    b_est <- fit$posterior_b[[s]]

    beta1_est <- fit$posterior_beta[[1]][selected_S[s],]
    beta2_est <- fit$posterior_beta[[2]][selected_S[s],]

    gamma1_est <- fit$posterior_gamma[[1]][selected_S[s],]
    gamma2_est <- fit$posterior_gamma[[2]][selected_S[s],]

    mu1 <- exp(X1%*%beta1_est + b_est[,1])
    mu2 <- exp(X2%*%beta2_est + b_est[,2])

    nu1 <- exp(X1%*%gamma1_est)
    nu2 <- exp(X2%*%gamma2_est)

    ll1[,s] <- llk_cmp(y1, mu1, nu1)
    ll2[,s] <- llk_cmp(y2, mu2, nu2)
  }

  beta1_esp <- fit$estimation_beta[[1]]
  beta2_esp <- fit$estimation_beta[[2]]

  gamma1_esp <- fit$estimation_gamma[[1]]
  gamma2_esp <- fit$estimation_gamma[[2]]

  #Random Effects
  post_b <- fit$posterior_b
  post_b1 = matrix(nrow = n, ncol = length(post_b))
  post_b2 = matrix(nrow = n, ncol = length(post_b))

  for(i in 1:length(post_b)){
    post_b1[,i] = post_b[[i]][,1]
    post_b2[,i] = post_b[[i]][,2]
  }

  b1_esp = apply(post_b1,1,mean)
  b2_esp = apply(post_b2,1,mean)

  mu_esp <- exp(cbind(X1%*%beta1_esp + b1_esp, X2%*%beta2_esp + b2_esp))
  nu_esp <- exp(cbind(X1%*%gamma1_esp, X2%*%gamma2_esp))

  d_bar <- cbind(llk_cmp(y1, mu_esp[,1], nu_esp[,1]),
                 llk_cmp(y2, mu_esp[,2], nu_esp[,2]))

  p_dic <- 2*(d_bar - cbind(rowMeans(ll1), rowMeans(ll2)))

  DIC_cmp <- colSums(-2*d_bar + p_dic)

  return(DIC_cmp)
}
