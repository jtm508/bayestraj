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
#' @importFrom mvtnorm rmvnorm
#'
#' @export

drawbeta = function(X,Y,sigma,mu,Sigma) {
  SigInv = solve(Sigma)
  A = solve(crossprod(X) / sigma + SigInv)
  B = crossprod(X,Y) / sigma + SigInv %*% mu
  beta = rmvnorm(1,A %*% B,A)
  return(beta)
}