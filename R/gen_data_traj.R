#' gen_data_traj
#'
#' Generate data for a single trajectory model
#' @param N: Integer, number of tandems
#' @param T1: Integer, number of time periods
#' @param pi1: Vector, probability of being assigned to each group
#' @param beta: Matrix, coefficients for each group
#' @param sigma: Float, variance for outcomes in series 1
#' @param scale: Boolean, TRUE to scale design matrix
#'
#'
#' @return List: \cr
#'    X1: Matrix, dataset \cr
#'    Y1: Vector, outcomes \cr
#' @export

gen_data_traj = function(N,T1,pi,beta,sigma,zeta=matrix(),scale=FALSE) {
  #generate IDs
  id = rep(1:N,each=T1)
  
  #dimension of design matrix
  d = dim(beta)[2]
  
  #number of random uniform covariates.
  ncov = d - 3
  
  #time dimension
  time = rep(1:T1,N)
  
  #generate design matrices
  X = cbind(rep(1,N*T1),
             matrix(runif(N*T1*ncov)*10,N*T1,ncov),
             time,
             time^2,deparse.level=FALSE)
  
  
  if (scale == TRUE) {
    X[,-1] = scale(X[,-1])
  }
  
  #generate group memberships
  K = dim(beta)[1]
  if (dim(zeta)[2] == 1) {
    Z = matrix()
    g = sample(c(1:K),N,replace=TRUE,prob=pi)
  }
  else {
    Z = cbind(1,rnorm(N),rnorm(N))
    u = Z %*% zeta
    pi = exp(u)/rowSums(exp(u))
    g = rep(0,N)
    for (i in 1:N){
      g[i] = sample(c(1:K),1,prob=as.vector(pi[i,]))
    }
  }
  
  #genereate outcomes
  index = g[id]
  
  Y = rowSums(X*beta[index,]) + rnorm(N*T1,sd=sqrt(sigma))
  
  #replace 1st column with id
  X[,1]=id
  
  return(list(X=X,Y=Y,g=g,Z=Z))
}