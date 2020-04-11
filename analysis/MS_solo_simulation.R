library(BayesTraj)
set.seed(3)

N=1000 #number of units
T=9 #time periods
pi=c(0.5,0.2,0.3) #group membership probabilities
#coefficients
beta=matrix(c(110,5,-0.5,
              111,-2,0,
              118,0,0),nrow=3,ncol=3,byrow=TRUE)
sigma=4 #variance of outcomes

set.seed(1)
data = gen_data(N=N,
                T=T,
                pi=pi,
                beta=beta,
                sigma=sigma,
                poly = 2 #degree of polynomial
)
#unpack data
X=data$X
y=data$Y
g=data$g
K = length(pi)

#Estimate model with model selection
iter = 5000
thin = 1
z = matrix(1,nrow=K,ncol=dim(X)[2])
model = trajMS(X=X, #data matrix
               y=y, #outcomes
               K=K, #number of groups
               time_index=2, #column of X corresponding to time
               iterations=iter, #number of iterations
               thin=thin, #thinning
               dispIter=1000) #Print a message every 1000 iterations

#print results
burn = 0.9
summary = summary_single_MS(model,X,y,burn)
print(summary$estimates)