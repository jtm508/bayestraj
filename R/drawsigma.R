#' drawsigma
#'
#' Draw sigma from posterior distribution
#'
#' @param X: Matrix, design matrix
#' @param Y: Vector, outcomes
#' @param beta: Matrix, coefficients for Y in each row
#' @param index: Vector, group membership for each row of X and Y
#' @param nu: Float, prior parameter
#' @param sigma0: Float, prior parameter
#' @param N: Integer, sample size
#'
#' @importFrom MCMCpack rinvgamma
#'
#' @export

drawsigma = function(X, Y, beta, index, nu, sigma0, N) {
  SSR = sum((Y - rowSums(X * beta[index,]))^2)
  sigma = rinvgamma(1, 0.5 * (nu + N), 0.5 * (nu * sigma0 + SSR))
  return(sigma)
}