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
#' @param z1: Matrix, K1 x dim(X1)[2] indicator matrix indicating which variables to include in each group.
#' @param z2: Matrix, K2 x dim(X2)[2] indicator matrix indicating which variables to include in each group.
#' @param iterations: Integer, number of MCMC iterations
#' @param thin: Integer, store every 'thin' iteration
#' @param dispIter: Integer, frequency of printing the iteration number
#' @param ll: Boolean, Set to TRUE to display the maximum log-likelihood over all the draws.
#' @param lambda: Numeric, prior for beta coefficients are N(0,lambda*I) where I is identity matrix
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#' @importFrom mvtnorm rmvnorm
#'
#' @export

dualtraj = function(X1,X2,y1,y2,K1,K2,z1,z2,iterations,thin=1,dispIter=10,ll=FALSE,lambda=100) {

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

#dimensions of theta
d1 = dim(X1)[2]
d2 = dim(X2)[2]

#hyperparameters
alpha1 = rep(1,K1)
alpha2 = rep(1,K2)
nu0 = 0.001
sigma0 = 1
Sigma1 = list()
for (j in 1:K1)
  Sigma1[[j]] = lambda*diag(sum(z1[j,])) #Beta1 prior. Scale identity matrix by lambda coefficient to make non-informative
Sigma2 = list()
for (j in 1:K2)
  Sigma2[[j]] = lambda*diag(sum(z2[j,]))

#initialize parameters
c1 = sample(c(1:K1),N,replace=TRUE) #randomly assign groups
c2 = sample(c(1:K2),N,replace=TRUE)
pi1 = as.vector(rdirichlet(1,alpha1)) #randomly assign group membership probabilities
pi1_2 = matrix(nrow=K1,ncol=K2)
for (j in 1:K1)
  pi1_2[j,] = rdirichlet(1,alpha2)
sigma1 = 1 #arbitrarily set to 1
sigma2 = 1
beta1=matrix(0,nrow=K1,ncol=d1,byrow=TRUE) #set to 0
beta2=matrix(0,nrow=K2,ncol=d2,byrow=TRUE)
mu1 = list() #set to 0
for (j in 1:K1)
  mu1[[j]] = rep(0,sum(z1[j,]))
mu2 = list()
for (j in 1:K2)
  mu2[[j]] = rep(0,sum(z2[j,]))

#initialize storage
c1Store = matrix(nrow=iterations/thin, ncol=N)
c2Store = matrix(nrow=iterations/thin, ncol=N)
pi1Store = matrix(nrow=iterations/thin, ncol=K1)
pi1_2Store = array(NA,dim=c(iterations/thin,K1,K2))
beta1Store = list()
for (j in 1:K1)
  beta1Store[[j]] = matrix(nrow=iterations/thin, ncol=sum(z1[j,]))
beta2Store = list()
for (j in 1:K2)
  beta2Store[[j]] = matrix(nrow=iterations/thin, ncol=sum(z2[j,]))
sigma1Store = rep(NA,iterations/thin)
sigma2Store = rep(NA,iterations/thin)

maxll = -Inf

#MCMC
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
  
  #draw beta
  for (j in 1:K1) {
    beta1[j,] = rep(0,d1)
    beta1[j,z1[j,]==1] = drawbeta(X1[index1==j,z1[j,]==1], y1[index1==j], sigma1, mu1[[j]], Sigma1[[j]])
  }
  
  for (j in 1:K2) {
    beta2[j,] = rep(0,d2)
    beta2[j,z2[j,]==1] = drawbeta(X2[index2==j,z2[j,]==1], y2[index2==j], sigma2, mu2[[j]], Sigma2[[j]])
  }
  
  #draw variances
  sigma1 = drawsigma(X1, y1, beta1, index1, nu0, sigma0, N1)
  sigma2 = drawsigma(X2, y2, beta2, index2, nu0, sigma0, N2)
  
  if (ll == TRUE) {
    ll.c = log_lik_dual(X1,X2,y1,y2,pi1,pi1_2,beta1,beta2,sigma1,sigma2,id1,id2)
    maxll = max(maxll,ll.c)
  }
  
  #store results
  if (q %% thin == 0) {
    store = q/thin
    pi1Store[store,] = pi1
    pi1_2Store[store,,] = pi1_2
    for (j in 1:K1)
      beta1Store[[j]][store,] = beta1[j,z1[j,]==1]
    for (j in 1:K2)
      beta2Store[[j]][store,] = beta2[j,z2[j,]==1]
    c1Store[store,] = c1
    c2Store[store,] = c2
    sigma1Store[store] = sigma1
    sigma2Store[store] = sigma2
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

if (ll == TRUE) {
  print(maxll)
}

#return draws
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
            sigma2 = sigma2Store))
}
