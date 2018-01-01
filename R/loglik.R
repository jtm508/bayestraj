#' log_lik_obs
#'
#' Compute log likelihood of each observation
#'
#' @param X: Matrix, design matrix
#' @param y: Vector, outcomes
#' @param beta: Matrix, likelihood coefficients
#' @param sigma: Float or Vector, likelihood variance
#'
#' @export
log_lik_obs = function(X,y,beta,sigma) {
  if (length(sigma) == 1) {
    ll = dnorm(y, mean = X %*% t(beta), sd = sqrt(sigma), log = TRUE)
  }
  else {
    ll = t(dnorm(t(matrix(y,nrow=length(y),ncol=length(sigma))),
                 mean = t(X %*% t(beta)),
                 sd = sqrt(sigma),
                 log = TRUE))
  }
  return(ll)
}

#' log_lik
#'
#' Compute log likelihood of data for latent trajectory model
#'
#' @param X: Matrix, design matrix
#' @param y: Vector, outcomes
#' @param pi: Vector, group probabilities
#' @param beta: Matrix, likelihood coefficients
#' @param sigma: Float, likelihood variance
#' @param id: Vector, id for each observation
#'
#' @export

log_lik = function(X,y,pi,beta,sigma,id) {
  ll.obs = log_lik_obs(X,y,beta,sigma) #each observation under each class
  ll.id = rowsum(ll.obs,group=id) #sum over time periods
  ll.id = exp(t(t(ll.id) + log(pi))) #add logpi
  ll.sum = log(rowSums(ll.id)) #sum over classes
  return(sum(ll.sum)) #sum over individuals
}



#' log_lik_dual
#'
#' Compute log likelihood of data for dual latent trajectory model
#'
#' @param X1: Matrix, design matrix group 1
#' @param X2: Matrix, design matrix group 2
#' @param y1: Vector, outcomes group 1
#' @param y2: Vector, outcomes group 2
#' @param pi1: Vector, group probabilities group 1
#' @param pi1_2: Vector, transition probabilities from group 1 to group 2
#' @param beta1: Matrix, likelihood coefficients group 1
#' @param beta2: Matrix, likelihood coefficients group 2
#' @param sigma1: Float or vector, likelihood variance group 1
#' @param sigma2: Float or vector, likelihood variance group 2
#' @param id1: Vector, id for each observation group 1
#' @param id2: Vector, id for each observation group 2
#'
#' @export

log_lik_dual = function(X1,X2,y1,y2,pi1,pi1_2,beta1,beta2,sigma1,sigma2,id1,id2) {
  ll.obs.1 = log_lik_obs(X1,y1,beta1,sigma1) #each observation under each class
  ll.id.1 = rowsum(ll.obs.1,group=id1) #sum over time periods
  ll.id.1 = exp(t(t(ll.id.1) + log(pi1))) #add logpi
  
  ll.obs.2 = log_lik_obs(X2,y2,beta2,sigma2) #each observation under each class 
  ll.id.2 = rowsum(ll.obs.2,group=id2) #sum over time periods
  
  #is there a way to vectorize this?
  K1 = dim(beta1)[1]
  for (k in 1:K1) {
    trans = rowSums(exp(t(t(ll.id.2) + log(pi1_2[k,])))) #transition likelihood summed over transitions from class k
    ll.id.1[,k] = ll.id.1[,k] * trans #multiply by transition likelihood
  }
  
  ll.sum = log(rowSums(ll.id.1)) #sum over classes in group 1
  return(sum(ll.sum)) #sum over individuals
}