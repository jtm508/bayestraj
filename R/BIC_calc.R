#' BIC_calc
#'
#' Computes BIC for latent trajectory model:
#' BIC = log_likelihood - 0.5 x num_parameters x log(num_observations)
#'
#' @param X: Matrix, design matrix
#' @param y: Vector, outcomes
#' @param pi: Vector, group probabilities
#' @param beta: Matrix, likelihood coefficients
#' @param sigma: Float, likelihood variance
#' @param id: Vector, id for each observation
#' @param z: Matrix, K x dim(X)[2] matrix, indicating which variables to include in each group.
#'
#' @export

BIC_calc = function(X,y,pi,beta,sigma,id,z) {
  ll = log_lik(X,y,pi,beta,sigma,id)
  n = length(y)
  #if pi is K dimensional, it only has k-1 parameters because they must sum to 1
  k = sum(z) + length(pi)-1 + length(sigma)
  BIC = ll - 0.5*k*log(n)
  return(BIC)
}

#' BIC_calc_dual
#'
#' Compute BIC for dual latent trajectory model
#' BIC = log_likelihood - 0.5 x num_parameters x log(num_observations)
#'
#' @param X1: Matrix, design matrix group 1
#' @param X2: Matrix, design matrix group 2
#' @param y1: Vector, outcomes group 1
#' @param y2: Vector, outcomes group 2
#' @param pi1: Vector, group probabilities group 1
#' @param pi1_2: Vector, transition probabilities from group 1 to group 2
#' @param beta1: Matrix, likelihood coefficients group 1
#' @param beta2: Matrix, likelihood coefficients group 2
#' @param sigma1: Float, likelihood variance group 1
#' @param sigma2: Float, likelihood variance group 2
#' @param id1: Vector, id for each observation group 1
#' @param id2: Vector, id for each observation group 2
#' @param z1: Matrix, K1 x dim(X1)[2] matrix, indicating which variables to include in each group.
#' @param z2: Matrix, K2 x dim(X2)[2] matrix, indicating which variables to include in each group.
#'
#' @export
#' 
BIC_calc_dual = function(X1,X2,y1,y2,pi1,pi1_2,beta1,beta2,sigma1,sigma2,id1,id2,z1,z2,constrain=FALSE) {
  ll = log_lik_dual(X1,X2,y1,y2,pi1,pi1_2,beta1,beta2,sigma1,sigma2,id1,id2)
  n = length(y1) + length(y2)
  #For group memberships, use the number of joint probabilities. Don't need marginal, transitions
  #etc. because these are all determined from the joint distribution. Again subtract 1 because they 
  #must sumto 1
  k = sum(z1) + sum(z2)+ length(pi1_2)-1 + length(sigma1) + length(sigma2)
  if (constrain == TRUE)
    #these are not sepearate parameters in the constrained model
    k = k - sum(z2) - length(sigma2)
  BIC = ll - 0.5*k*log(n)
  return(BIC)
}