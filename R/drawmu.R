#' drawmu
#'
#' Draw mu from posterior distribution
#'
#' @param beta: beta: Matrix, coefficients for outcome in each row
#' @param Sigma: Matrix, prior covariance of beta
#' @param Lambda: Matrix, prior covariance of mu
#' @param K: Number of groups
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @export

drawmu = function(beta, Sigma, Lambda, K) {
  SigInv = solve(Lambda)
  A = solve(K * Sigma + SigInv)
  B = K * solve(Sigma) %*% colMeans(beta)
  mu = as.vector(rmvnorm(1, A %*% B, A))
  return(mu)
}