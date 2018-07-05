#' dualtraj_consrained
#'
#' Estimate the dual trajectory model using a Gibbs Sampler
#'
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param y1: Vector, outcomes for series 1
#' @param y2: Vector, outcomes for series 2
#' @param K: Integer, number of latent classes
#' @param z: Matrix, K x dim(X)[2] indicator matrix indicating which variables to inlcude in each group.
#' @param iterations: Integer, number of MCMC iterations
#' @param thin: Integer, store every 'thin' iteration
#' @param dispIter: Integer, frequency of printing the iteration number
#' @param ll: Boolean, Set to TRUE to display the maximum log-likelihood over all the draws.
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#' @importFrom mvtnorm rmvnorm
#'
#' @export

dualtraj_constrained = function(X1,X2,y1,y2,K,z,iterations,thin=1,dispIter=10,ll=FALSE) {

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

#replace id with intercept column
X1[,1]=1
X2[,1]=1

#get complete datasets:
X_ = rbind(X1,X2)
y_ = c(y1,y2)
N_ = N1+N2

#dimensions of theta
#assert that dim(X1)[2]=dim(X2)[2]
d = dim(X1)[2]

#hyperparameters
alpha1 = rep(1,K)
alpha2 = rep(1,K)
nu0 = 0.001
sigma0 = 1

#initialize parameters
c1 = sample(c(1:K),N,replace=TRUE)
c2 = sample(c(1:K),N,replace=TRUE)
pi1 = as.vector(rdirichlet(1,alpha1))
pi1_2 = matrix(nrow=K,ncol=K)
for (j in 1:K)
  pi1_2[j,] = rdirichlet(1,alpha2)
sigma = 1
Sigma = list()
for (j in 1:K)
  Sigma[[j]] = 100*diag(sum(z[j,]))
beta=matrix(0,nrow=K,ncol=d,byrow=TRUE)
mu = list()
for (j in 1:K)
  mu[[j]] = rep(0,sum(z[j,]))

#initialize storage
c1Store = matrix(nrow=iterations/thin, ncol=N)
c2Store = matrix(nrow=iterations/thin, ncol=N)
pi1Store = matrix(nrow=iterations/thin, ncol=K)
pi1_2Store = array(NA,dim=c(iterations/thin,K,K))
betaStore = list()
for (j in 1:K)
  betaStore[[j]] = matrix(nrow=iterations/thin, ncol=sum(z[j,]))
sigmaStore = rep(NA,iterations/thin)

maxll = -Inf

for (q in 1:iterations) {
  if (q %% dispIter == 0) {
    print(q)
  }
  
  #draw groups
  c1 = drawgroup(X1,y1,N,id1,c2,pi1,pi1_2,beta,sigma,K)
  c2 = drawgroup2(X2,y2,N,id2,c1,pi1,pi1_2,beta,sigma,K)
  
  #reindex according to new groups
  index1 = c1[id1]
  index2 = c2[id2]
  index = c(index1,index2)
  
  #draw group/transition probabilities
  pi1 = drawpi(c1, alpha1, K)
  for (j in 1:K)
    pi1_2[j,] = drawpi(c2[c1==j], alpha2, K)
  
  #draw beta
  for (j in 1:K) {
    beta[j,] = rep(0,d)
    beta[j,z[j,]==1] = drawbeta(X_[index==j,z[j,]==1], y_[index==j], sigma, mu[[j]], Sigma[[j]])
  }
  
  #draw variances
  sigma = drawsigma(X_, y_, beta, index, nu0, sigma0, N_)
  
  if (ll == TRUE) {
    ll.c = log_lik_dual(X1,X2,y1,y2,pi1,pi1_2,beta,beta,sigma,sigma,id1,id2)
    maxll = max(maxll,ll.c)
  }
  
  
  if (q %% thin == 0) {
    store = q/thin
    pi1Store[store,] = pi1
    pi1_2Store[store,,] = pi1_2
    for (j in 1:K)
      betaStore[[j]][store,] = beta[j,z[j,]==1]
    c1Store[store,] = c1
    c2Store[store,] = c2
    sigmaStore[store] = sigma
  }
}

#calculate other class probabilities

#marginal probability of group2
pi2 = matrix(nrow=iterations/thin,ncol=K)
for (i in 1:K) {
  pi2[,i] = rowSums(pi1Store * pi1_2Store[,,i])
}

#joint probability
pi12 = array(NA,dim=c(iterations/thin,K,K))
for (i in 1:K) {
  for (j in 1:K) {
    pi12[,i,j] = pi1Store[,i] * pi1_2Store[,i,j]
  }
}

#group 1 probability conditional on group 2
pi2_1 = array(NA,dim=c(iterations/thin,K,K))
for (i in 1:K) {
  for (j in 1:K) {
    pi2_1[,i,j] = pi12[,i,j] / pi2[,j]
  }
}

if (ll == TRUE) {
  print(maxll)
}

#return draws
return(list(beta   = betaStore,
            c1     = c1Store,
            c2     = c2Store,
            pi1    = pi1Store,
            pi2    = pi2,
            pi12   = pi12,
            pi1_2  = pi1_2Store,
            pi2_1  = pi2_1,
            sigma = sigmaStore))
}
