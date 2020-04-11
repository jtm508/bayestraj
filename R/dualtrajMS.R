#' dualtraj
#'
#' Estimate the dual trajectory model using a Gibbs Sampler
#'
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param y1: Vector, outcomes for series 1
#' @param y2: Vector, outcomes for series 2
#' @param K1: Integer, number of latent classes in series 1
#' @param K2: Integer, number of latent classes in series 2
#' @param time_index: Integer, column of X corresponding to time
#' @param iterations: Integer, number of MCMC iterations
#' @param thin: Integer, store every 'thin' iteration
#' @param dispIter: Integer, frequency of printing the iteration number
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#' @importFrom mvtnorm rmvnorm
#'
#' @export

dualtrajMS = function(X1,X2,y1,y2,K1,K2,time_index,iterations,thin=1,dispIter=10) {
  #extract ids from design matrices
  id1 = X1[,1]
  id2 = X2[,1]
  
  #relabel id's to begin at 1
  uniqueIDs = sort(unique(c(id1,id2)))
  id1 = match(id1,uniqueIDs)
  id2 = match(id2,uniqueIDs)
  
  #number of pairs
  N = length(unique(id1))
  
  #length of data
  N1 = length(id1)
  N2 = length(id2)
  
  #number of non-time covariates
  ncov = time_index - 2
  
  #replace id with intercept column
  X1[,1]=1
  X2[,1]=1
  
  #dimensions of theta
  d1 = dim(X1)[2]
  d2 = dim(X2)[2]
  
  #hyperparameters
  alpha1 = rep(1,K1)
  alpha2 = rep(1,K2)
  
  #initialize parameters
  c1 = sample(c(1:K1),N,replace=TRUE)
  c2 = sample(c(1:K2),N,replace=TRUE)
  pi1 = as.vector(rdirichlet(1,alpha1))
  pi1_2 = matrix(nrow=K1,ncol=K2)
  for (j in 1:K1)
    pi1_2[j,] = rdirichlet(1,alpha2)
  sigma1 = rep(1,K1)
  sigma2 = rep(1,K2)
  beta1=matrix(0,nrow=K1,ncol=d1,byrow=TRUE)
  beta2=matrix(0,nrow=K2,ncol=d2,byrow=TRUE)
  z1 = matrix(1,nrow=K1,ncol=d1)
  z2 = matrix(1,nrow=K2,ncol=d2)
  marg.lik1.c = rep(0,K1)
  marg.lik2.c = rep(0,K2)
  
  #initialize storage
  c1Store = matrix(nrow=iterations/thin, ncol=N)
  c2Store = matrix(nrow=iterations/thin, ncol=N)
  pi1Store = matrix(nrow=iterations/thin, ncol=K1)
  pi1_2Store = array(NA,dim=c(iterations/thin,K1,K2))
  beta1Store = list()
  z1Store = list()
  for (j in 1:K1) {
    beta1Store[[j]] = matrix(nrow=iterations/thin, ncol=d1)
    z1Store[[j]] = matrix(nrow=iterations/thin, ncol=d1)
  }
  beta2Store = list()
  z2Store = list()
  for (j in 1:K2) {
    beta2Store[[j]] = matrix(nrow=iterations/thin, ncol=d2)
    z2Store[[j]] = matrix(nrow=iterations/thin, ncol=d1)
  }
  sigma1Store = matrix(nrow=iterations/thin, ncol=K1)
  sigma2Store = matrix(nrow=iterations/thin, ncol=K2)
  
  for (q in 1:iterations) {
    if (q %% dispIter == 0) {
      print(q)
    }
    
    #draw groups
    c1 = drawgroup_dual(X1,y1,N,id1,c2,pi1,pi1_2,beta1,sigma1,K1)
    c2 = drawgroup2(X2,y2,N,id2,c1,pi1,pi1_2,beta2,sigma2,K2)
    
    #reindex according to new groups
    index1 = c1[id1]
    index2 = c2[id2]
    
    #draw group/transition probabilities
    pi1 = drawpi(c1, alpha1, K1)
    for (j in 1:K1)
      pi1_2[j,] = drawpi(c2[c1==j], alpha2, K2)
    
    #draw group parameters
    for (j in 1:K1) {
      #recalculate marginal likelihood for new group memberships
      marg.lik1.c[j] = marg_lik(y1[index1==j], X1[index1==j,z1[j,]==1,drop=FALSE])
      #sample z
      if (ncov > 0) {
        for (b in resamp(2:(1+ncov))) {
          zlist = drawz(z1,j,b,y1,X1,index1,marg.lik1.c[j])
          z1[j,b] = zlist[[1]]
          marg.lik1.c[j] = zlist[[2]]
        }
      }
      #draw polynomial degree
      z1[j,time_index:d1] = drawpoly(y1[index1==j], X1[index1==j,,drop=FALSE],time_index,z1[j,])
      
      #some calculations
      X.temp = X1[index1==j,z1[j,]==1,drop=FALSE]
      y.temp = y1[index1==j]
      n = length(y.temp)
      g = n
      A = crossprod(X.temp)
      beta.ols = as.vector(solve(A,crossprod(X.temp,y.temp)))
      SSR = sum((y.temp - X.temp %*% beta.ols)^2)
      
      #sample sigma^2
      sigma1[j] = rinvgamma(1,n/2,SSR/2)
      
      #sample beta
      beta1[j,] = rep(0,d1)
      beta1[j,z1[j,]==1] = as.vector(rmvnorm(1, beta.ols, g / (g + 1)  * sigma1[j] * solve(A)))
    }
    
    for (j in 1:K2) {
      #recalculate marginal likelihood for new group memberships
      marg.lik2.c[j] = marg_lik(y2[index2==j], X2[index2==j,z2[j,]==1,drop=FALSE])
      if (ncov > 0) {
        for (b in resamp(2:(1+ncov))) {
          zlist = drawz(z2,j,b,y2,X2,index2,marg.lik2.c[j])
          z2[j,b] = zlist[[1]]
          marg.lik2.c[j] = zlist[[2]]
        }
      }
      #draw polynomial degree
      z2[j,time_index:d2] = drawpoly(y2[index2==j], X2[index2==j,,drop=FALSE],time_index,z2[j,])
      
      #some calculations
      X.temp = X2[index2==j,z2[j,]==1,drop=FALSE]
      y.temp = y2[index2==j]
      n = length(y.temp)
      g = n
      A = crossprod(X.temp)
      beta.ols = as.vector(solve(A, crossprod(X.temp, y.temp)))
      SSR = sum((y.temp - X.temp %*% beta.ols)^2)
      
      #sample sigma^2
      sigma2[j] = rinvgamma(1,n/2,SSR/2)
      
      #sample beta
      beta2[j,] = rep(0,d2)
      beta2[j,z2[j,]==1] = as.vector(rmvnorm(1, beta.ols, g / (g + 1)  * sigma2[j] * solve(A)))
    }
    
    if (q %% thin == 0) {
      store = q/thin
      pi1Store[store,] = pi1
      pi1_2Store[store,,] = pi1_2
      for (j in 1:K1) {
        beta1Store[[j]][store,] = beta1[j,]
        z1Store[[j]][store,] = z1[j,]
      }
      for (j in 1:K2) {
        beta2Store[[j]][store,] = beta2[j,]
        z2Store[[j]][store,] = z2[j,]
      }
      c1Store[store,] = c1
      c2Store[store,] = c2
      sigma1Store[store,] = sigma1
      sigma2Store[store,] = sigma2
    }
  }
  
  #calculate other class probabilities
  
  #marginal probability of group2
  pi2 = matrix(nrow=iterations/thin,ncol=K2)
  for (i in 1:K2) {
    pi2[,i] = rowSums(pi1Store * pi1_2Store[,,i])
  }
  
  #joint probability
  pi12 = array(NA,dim=c(iterations/thin,K1,K2))
  for (i in 1:K1) {
    for (j in 1:K2) {
      pi12[,i,j] = pi1Store[,i] * pi1_2Store[,i,j]
    }
  }
  
  #group 1 probability conditional on group 2
  pi2_1 = array(NA,dim=c(iterations/thin,K1,K2))
  for (i in 1:K1) {
    for (j in 1:K2) {
      pi2_1[,i,j] = pi12[,i,j] / pi2[,j]
    }
  }
  
  return(list(beta1  = beta1Store,
              beta2  = beta2Store,
              c1     = c1Store,
              c2     = c2Store,
              pi1    = pi1Store,
              pi2    = pi2,
              pi12   = pi12,
              pi1_2  = pi1_2Store,
              pi2_1  = pi2_1,
              sigma1 = sigma1Store,
              sigma2 = sigma2Store,
              z1     = z1Store,
              z2     = z2Store))
  
}

#fuction to ensure proper sampling of z
resamp = function(x) {
  if(length(x)==1) 
    x 
  else 
    sample(x)
}