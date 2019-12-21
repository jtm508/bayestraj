#' trajMS
#'
#' Estimate the trajectory model with model averaging
#'
#' @param X: Matrix, design matrix for series 1. 1st column should be the id. Last columns of X correspond to time and polynomials of time.
#' @param y: Vector, outcomes for series 1
#' @param K: Integer, number of latent classes in series 1
#' @param time_index: Integer, column of X corresponding to time
#' @param iterations: Integer, number of MCMC iterations
#' @param thin: Integer, store every 'thin' iteration
#' @param dispIter: Integer, frequency of printing the iteration number
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#' @importFrom mvtnorm rmvnorm
#'
#' @export

trajMS = function(X,y,K,time_index,iterations,thin=1,dispIter=10) {
  #extract ids from design matrices
  id = X[,1]
  
  #relabel id's to begin at 1
  uniqueIDs = sort(unique(c(id)))
  id = match(id,uniqueIDs)
  
  #number of pairs
  N = length(unique(id))
  
  #length of data
  N_ = length(id)
  
  #number of non-time covariates
  ncov = time_index - 2
  
  #replace id with intercept column
  X[,1]=1
  
  #dimensions of theta
  d = dim(X)[2]
  
  #hyperparameters
  alpha = rep(1,K)
  
  #initialize parameters
  c = sample(c(1:K),N,replace=TRUE)
  pi = as.vector(rdirichlet(1,alpha))
  sigma = rep(1,K)
  beta=matrix(0,nrow=K,ncol=d,byrow=TRUE)
  z = matrix(1,nrow=K,ncol=d)
  z[,1] = -1 #intercept column. Code as -1 so its not double counted
  marg.lik.c = rep(0,K)
  
  #initialize storage
  cStore = matrix(nrow=iterations/thin, ncol=N)
  piStore = matrix(nrow=iterations/thin, ncol=K)
  betaStore = list()
  zStore = list()
  for (j in 1:K) {
    betaStore[[j]] = matrix(nrow=iterations/thin, ncol=d)
    zStore[[j]] = matrix(nrow=iterations/thin, ncol=d)
  }
  sigmaStore = matrix(nrow=iterations/thin, ncol=K)
  
  for (q in 1:iterations) {
    if (q %% dispIter == 0) {
      print(q)
    }
    
    #draw groups
    c = drawgroup(X,y,N,id,pi,beta,sigma,K)
    #print(table(c))
    
    #reindex according to new groups
    index = c[id]
    
    #draw group/transition probabilities
    pi = drawpi(c, alpha, K)
    
    #draw group parameters
    for (j in 1:K) {
      #recalculate marginal likelihood for new group memberships
      marg.lik.c[j] = marg_lik(y[index==j], X[index==j,z[j,] %in% c(-1,1),drop=FALSE])
      #sample z
      if (ncov > 0) {
        for (b in resamp(2:(1+ncov))) {
          zlist = drawz(z,j,b,y,X,index,marg.lik.c[j])
          z[j,b] = zlist[[1]]
          marg.lik.c[j] = zlist[[2]]
        }
      }
      #draw polynomial degree
      z[j,time_index:d] = drawpoly(y[index==j], X[index==j,,drop=FALSE],time_index,z[j,])
      
      #some calculations
      X.temp = X[index==j,z[j,] %in% c(-1,1),drop=FALSE]
      y.temp = y[index==j]
      n = length(y.temp)
      g = n
      A = crossprod(X.temp)
      beta.ols = as.vector(solve(A,crossprod(X.temp,y.temp)))
      SSR = sum((y.temp - X.temp %*% beta.ols)^2)
      
      #sample sigma^2
      X.tilde = X[index==j,z[j,]==1,drop=FALSE] #does not include intercept
      y.tilde = y[index==j] - beta[j,1] #subtract intercept
      A.tilde = crossprod(X.tilde)
      beta.ols = as.vector(solve(A.tilde,crossprod(X.tilde,y.tilde)))
      sigma[j] = rinvgamma(1,n/2,SSR/2 + (crossprod(beta.ols,A.tilde) %*% beta.ols) / (2 * (g + 1)))
      
      #sample beta
      beta[j,-1] = 0 #do not overwrite intercept
      beta[j,z[j,]==1] = as.vector(rmvnorm(1, g / (g + 1) * beta.ols, g / (g + 1)  * sigma[j] * solve(A.tilde)))
      beta[j,1] = drawintercept(X.tilde,y.temp,beta[j,z[j,]==1],sigma[j],n)
    }
    
    
    if (q %% thin == 0) {
      store = q/thin
      piStore[store,] = pi
      for (j in 1:K) {
        betaStore[[j]][store,] = beta[j,]
        zStore[[j]][store,] = z[j,]
      }
      cStore[store,] = c
      sigmaStore[store,] = sigma
    }
  }
  
  return(list(beta  = betaStore,
              c     = cStore,
              pi    = piStore,
              sigma = sigmaStore,
              z     = zStore))
  
}

#fuction to ensure proper sampling of z
resamp = function(x) {
  if(length(x)==1) 
    x 
  else 
    sample(x)
}