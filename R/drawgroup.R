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

#drawgroup_old = function(X,Y,id,g2,pi1,pi2,beta,sigma,K) {
#  
#  #likelihood of each observation under each group
#  likelihood = dnorm(Y,mean=X %*% t(beta),sd=sqrt(sigma))
#  for (i in 1:N) {
#    denom = apply(likelihood[id1==i,],2,prod) * pi2[,g2[i]] * pi1
#    prob = denom/sum(denom)
#    g[i] = sample(c(1:K),1,prob=prob)
#  }
#  return(g)
#}

drawgroup = function(X,Y,N,id,g2,pi1,pi2,beta,sigma,K) {
  logpi1 = log(pi1)
  logpi2 = log(pi2)
  
  #log-likelihood of each observation under each group
  ll = dnorm(Y, mean = X %*% t(beta), sd = sqrt(sigma), log = TRUE)
  
  g = rep(NA,N)
  for (i in 1:N) {
    denom = colSums(ll[id==i,]) + logpi2[,g2[i]] + logpi1
    numer = logsumexp(denom)
    prob = exp(denom - numer)
    g[i] = sample(c(1:K), 1, prob=prob)
  }
  return(g)
}