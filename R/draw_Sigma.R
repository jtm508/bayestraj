#' drawbeta
#'
#' Draw beta from posterior distribution
#'
#' @param X: Matrix, design matrix
#' @param Y: Vector, outcomes
#' @param sigma: Float, variance of Y
#' @param mu: Vector, prior mean of beta
#' @param Sigma: Matrix, prior covariance of beta
#'
#' @importFrom MCMCpack riwish
#'
#' @export

drawSigma = function(beta, mu, nu, V, K) {
  diff = sweep(beta, 2, mu)
  S = crossprod(diff)
  Sigma = riwish(nu + K, V + S)
  return(Sigma)
}
  