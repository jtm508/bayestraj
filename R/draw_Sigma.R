#' drawSigma
#'
#' Draw Sigma (covariance) from posterior distribution
#'
#' @param X: Matrix, design matrix
#' @param Y: Vector, outcomes
#' @param mu: Vector, prior mean of beta
#' @param nu: Inverse Wishart spread hypeparameter
#' @param V: Inverse Wishart location hyperparameter
#' @param K: Number of groups
#'
#' @importFrom MCMCpack riwish
#'
#' @export

drawSigma = function(beta, mu, nu, V, K) {
  diff = t(t(beta) - mu)
  S = crossprod(diff)
  Sigma = riwish(nu + K, V + S)
  return(Sigma)
}
  