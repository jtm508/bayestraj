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
#'
#' @importFrom MCMCpack rdirichlet riwish rinvgamma
#'
#' @export

dualtraj = function(X1,X2,Y1,Y2,K1,K2,iterations) {

#extract ids from design matrices
id1 = X1[,1]
id2 = X2[,2]

#number of tandems
N = length(unique(id1))

#replace id with intercept column
X1[,1]=1
X2[,1]=1

for (q in 1:iterations) {
  if (q %% 10 == 0)
    print(q)
}
}
