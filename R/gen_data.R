#' gen_data
#'
#' Generate data for a dual trajectory model
#' @param N: Integer, number of tandems
#' @param T1: Integer, number of time periods in series 1
#' @param T2: Integer, number of time periods in series 2
#' @param pi1: Vector, probability of being assigned to each group in series 1
#' @param pi2: Matrix, transition matrix. Row j is the transition probabilities for when the series 1 group is j
#' @param beta1: List, coefficients for each group in series 1
#' @param beta2: List, coefficients for each group in series 2
#' @param sigma1: Float, variance for outcomes in series 1
#' @param sigma2: Float, variance for outcomes in series 2
#'
#'
#' @return List: \cr
#'    X1: Matrix, dataset for series 1 \cr
#'    X2: Matrix, dataset for series 2 \cr
#'    Y1: Vector, outcomes for series 1 \cr
#'    Y2: Vector, outcomes for series 2 \cr
#' @export

gen_data = function(N,T1,T2,pi1,pi2,beta1,beta2,sigma1,sigma2) {
  if(length(pi1)!=dim(pi2)[1])
    stop("Incompatible dimensions for pi1 and pi2")

  if(length(pi1)!=length(beta1))
    stop("Incompatible dimensions for pi1 and beta1")

  if(dim(pi2)[2]!=length(beta2))
    stop("Incompatible dimensions for pi2 and beta2")

  #generate IDs
  id1 = rep(1:N,each=T1)
  id2 = rep(1:N,each=T2)

  #dimension of design matrix
  d1 = length(beta1[[1]])
  d2 = length(beta2[[1]])

  #number of random uniform covariates.
  ncov1 = d1 - 3
  ncov2 = d2 - 3

  #time dimension
  time1 = rep(1:T1,N)
  time2 = rep(1:T2,N)

  #generate design matrices
  X1 = cbind(rep(1,N*T1),
             matrix(runif(N*T1*ncov1),N*T1,ncov1),
             time1,
             time1^2)

  X2 = cbind(rep(1,N*T2),
             matrix(runif(N*T2*ncov2),N*T2,ncov2),
             time2,
             time2^2)

  #generate group memberships
  K1 = length(beta1)
  K2 = length(beta2)
  g1 = sample(c(1:K1),N,replace=TRUE,prob=pi1)
  g2 = numeric(N)
  for (i in 1:N) {
    g2[i] = sample(c(1:K2),1,replace=TRUE,prob=pi2[g1[i],])
  }

  #genereate outcomes
  Y1 = numeric(N*T1)
  for (i in 1:K1){
    index = id1 %in% which(g1==i)
    Y1[index] = X1[index,]%*%beta1[[i]] + rnorm(n=length(which(index)),sd=sqrt(sigma1))
  }

  Y2 = numeric(N*T2)
  for (i in 1:K2){
    index = id2 %in% which(g2==i)
    Y2[index] = X2[index,]%*%beta2[[i]] + rnorm(n=length(which(index)),sd=sqrt(sigma2))
  }

  return(list(X1=X1,X2=X2,Y1=Y1,Y2=Y2))
}
