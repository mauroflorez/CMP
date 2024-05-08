#' Product of lists between matrices
#'
#' @param X Data
#' @param beta Parameters
#'
#' @return A list with the products element-wise
#' @export
prod_list = function(X, beta){
  matrix(unlist(purrr::map2(X, beta, `%*%`)), ncol = length(beta))
}
