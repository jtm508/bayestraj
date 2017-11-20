#' dualtraj
#'
#' Estimate the dual trajectory model using a Gibbs Sampler
#'
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param Y1: Vector, outcomes for series 1
#' @param Y2: Vector, outcomes for series 2
#' @param K1: Integer, number of latent classes in series 1
#' @param K2: Integer, number of latent classes in series 2
#' @param iterations: Integer, number of MCMC iterations
#' @param thin: Integer, store every 'thin' iteration
#' @param dispIter: Integer, frequency of printing the iteration number
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#' @importFrom mvtnorm rmvnorm
#'
#' @export

dualtraj = function(X1,X2,Y1,Y2,g1,g2,K1,K2,iterations,thin=1,dispIter=10) {

#extract ids from design matrices
id1 = X1[,1]
id2 = X2[,1]

#number of tandems
N = length(unique(id1))

#length of data
N1 = length(id1)
N2 = length(id2)

#replace id with intercept column
X1[,1]=1
X2[,1]=1

#dimensions of theta
d1 = dim(X1)[2]
d2 = dim(X2)[2]

#hyperparameters
alpha1 = rep(1,K1)
alpha2 = rep(1,K2)
nu0 = 0.001
sigma0 = 1
nu1 = d1+3
nu2 = d2+3
V1 = nu1*diag(d1)
V2 = nu2*diag(d2)
Lambda1 = 10000*diag(d1)
Lambda2 = 10000*diag(d2)

#initialize parameters
pi1 = rdirichlet(1,alpha1)
pi2 = matrix(nrow=K1,ncol=K2)
for (j in 1:K1)
  pi2[j,] = rdirichlet(1,alpha2)
sigma1 = 1
sigma2 = 1
Sigma1 = diag(d1)
Sigma2 = diag(d2)
beta1=matrix(0,nrow=K1,ncol=d1,byrow=TRUE)
beta2=matrix(0,nrow=K2,ncol=d2,byrow=TRUE)
mu1 = rep(0,d1)
mu2 = rep(0,d2)

#initialize storage
pi1Store = matrix(nrow=iterations/thin, ncol=K1)
pi2Store = matrix(nrow=iterations/thin, ncol=K1*K2)
beta1Store = list()
for (j in 1:K1)
  beta1Store[[j]] = matrix(nrow=iterations/thin, ncol=d1)
beta2Store = list()
for (j in 1:K2)
  beta2Store[[j]] = matrix(nrow=iterations/thin, ncol=d2)
sigma1Store = rep(NA,iterations/thin)
sigma2Store = rep(NA,iterations/thin)

for (q in 1:iterations) {
  if (q %% dispIter == 0)
    print(q)
  
  #draw groups
  g1 = drawgroup(X1,Y1,N,id1,g2,pi1,pi2,beta1,sigma1)
  g2 = drawgroup2(X1,Y1,N,id2,g1,pi1,pi2,beta1,sigma1)
  
  #reindex according to new groups
  index1 = g1[id1]
  index2 = g2[id2]
  
  #draw group/transition probabilities
  pi1 = drawpi(g1, alpha1, K1)
  for (j in 1:K1)
    pi2[j,] = drawpi(g2[g1==j], alpha2, K2)
  
  #draw beta
  for (j in 1:K1)
    beta1[j,] = drawbeta(X1[index1==j,], Y1[index1==j], sigma1, mu1, Sigma1)
  
  for (j in 1:K2) 
    beta2[j,] = drawbeta(X2[index2==j,], Y2[index2==j], sigma2, mu2, Sigma2)
  
  #draw variances
  sigma1 = drawsigma(X1, Y1, beta1, index1, nu0, sigma0, N1)
  sigma2 = drawsigma(X2, Y2, beta2, index2, nu0, sigma0, N2)
  
  
  #draw prior means
  mu1 = drawmu(beta1, Sigma1, Lambda1, K1)
  mu2 = drawmu(beta2, Sigma2, Lambda2, K2)
  
  #draw prior covariances)
  Sigma1 = drawSigma(beta1, mu1, nu1, V1, K1)
  Sigma2 = drawSigma(beta2, mu2, nu2, V2, K2)
  
  if (q %% thin ==0) {
    store = q/thin
    pi1Store[store,] = pi1
    pi2Store[store,] = as.vector(t(pi2))
    for (j in 1:K1)
      beta1Store[[j]][store,] = beta1[j,]
    for (j in 1:K2)
      beta2Store[[j]][store,] = beta2[j,]
    sigma1Store[store] = sigma1
    sigma2Store[store] = sigma2
  }
}
return(list(beta1  = beta1Store,
            beta2  = beta2Store,
            pi1    = pi1Store,
            pi2    = pi2Store,
            sigma1 = sigma1Store,
            sigma2 = sigma2Store))
}
