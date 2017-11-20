#' drawbeta
#'
#' Draw beta from posterior distribution
#'
#' @param X: Matrix, design matrix
#' @param Y: Vector, outcomes
#' @param N: Integer, number of tandems
#' @param id: Vector, tandem for each row of X and Y
#' @param g2: Vector, group membership for series 2
#' @param pi1: Vector, group probabilities for series 1
#' @param pi2: Vector, transition probabilities for series 2
#' @param beta: Matrix, likelihood coefficients for series 1
#' @param sigma: Float, likelihood variance
#' @param K: Number of groups in series 1
#'
#' @export
#' 

drawgroup2 = function(X,Y,N,id,g1,pi1,pi2,beta,sigma,K) {
  logpi1 = log(pi1)
  logpi2 = log(pi2)
  
  #log-likelihood of each observation under each group
  ll = dnorm(Y, mean = X %*% t(beta), sd = sqrt(sigma), log = TRUE)
  
  for (i in 1:N) {
    denom = colSums(ll[id==i,]) + logpi2[g1[i],] + logpi1[g1[i]]
    numer = logsumexp(denom)
    prob = exp(denom - numer)
    g[i] = sample(c(1:K), 1, prob=prob)
  }
  return(g)
}