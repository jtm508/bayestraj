#' drawintercept
#'
#' Estimate the trajectory model with model averaging
#'
#' @param X: Matrix, design matrix. Does not include intercept column
#' @param y: Vector, outcomes
#' @param beta: Vector, regression coefficients
#' @param sigmasq: Float, regression variance
#' @param n: Integer, sample size
#'
#' @export

drawintercept = function(X,y,beta,sigmasq,n) {
  mu = mean(y-X%*%beta)
  sd = sqrt(sigmasq/n)
  intercept = rnorm(1,mu,sd)
  return(intercept)
}