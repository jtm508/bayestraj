#' drawgroup2
#'
#' Draw second series' groups from posterior distribution
#'
#' @param X: Matrix, design matrix
#' @param y: Vector, outcomes
#' @param N: Integer, number of tandems
#' @param id: Vector, tandem for each row of X and Y
#' @param g1: Vector, group membership for series 1
#' @param pi1: Vector, group probabilities for series 1
#' @param pi1_2: Vector, transition probabilities for series 2
#' @param beta: Matrix, likelihood coefficients for series 1
#' @param sigma: Float or vector, likelihood variance
#' @param K: Number of groups in series 1
#'
#' @export

drawgroup2 = function(X,y,N,id,g1,pi1,pi1_2,beta,sigma,K) {
  logpi1 = log(pi1)
  logpi1_2 = log(pi1_2)
  
  #log-likelihood of each observation under each group
  ll = log_lik_obs(X,y,beta,sigma)
  llSum = rowsum(ll,group=id)
  
  #log-likelihood of each observation under each group
  numer = llSum + logpi1_2[g1,] + logpi1[g1]
  maxNumer = apply(numer,1,max)
  denom = maxNumer + log(rowSums(exp(numer-maxNumer))) #logsumexp trick
  prob = exp(numer - denom)
  
  #get probabilities 
  u = runif(N)
  cumProb = prob %*% upper.tri(diag(ncol(prob)), diag = TRUE)
  g = rowSums(u > cumProb) + 1
  
  return(g)
}