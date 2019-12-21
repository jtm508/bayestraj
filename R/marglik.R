#' marg_lik
#'
#' Calculate log marginal likelihood under unit-information Zellner's g-prior
#'
#' @param y: Vector, outcomes
#' @param X: Matrix, design matrix
#'
#' @export

marg_lik = function(y,X) {
  n = dim(X)[1]
  p = dim(X)[2] #-1 #subtract 1 for intercept
  g = n
  b.ols = as.vector(solve(crossprod(X),crossprod(X,y)))
  R2 = 1 - crossprod(y - X %*% b.ols) / crossprod(y - mean(y))
  return(lgamma((n - 1) / 2) - 
           (n - 1) * log(sqrt(pi)) - 
           log(sqrt(n)) - 
           (n - 1) * log(sqrt(sum((y - mean(y))^2))) +
           (n - 1 - p) / 2 * log(1 + g) - 
           (n - 1) / 2 * log(1 + g * (1 - R2))
  )
}