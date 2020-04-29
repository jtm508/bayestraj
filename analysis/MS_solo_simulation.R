library(BayesTraj)
set.seed(3)

sim_number = 2

N=1000 #number of units
T=15 #time periods
if (sim_number==1) {
  pi=c(0.5,0.2,0.3) #group membership probabilities
  #coefficients
  beta=matrix(c(110,5,-0.5,
                111,-2,0,
                118,0,0),nrow=3,ncol=3,byrow=TRUE)
  sigma=2^2 #variance of outcomes
  poly = 2 #degree of polynomial
} else if (sim_number==2) {
  pi=c(0.4,0.2,0.3,0.1) #group membership probabilities
  #coefficeints
  beta=matrix(c(5,1,-0.1,0.01,
                4,0.3,0.1,0,
                7,-0.3,0,0,
                6,0,0,0),nrow=4,ncol=4,byrow=TRUE)
  #variance
  sigma=1^2
  poly = 3 #degree of polynomial
} 

data = gen_data(N=N,
                T=T,
                pi=pi,
                beta=beta,
                sigma=sigma,
                poly = poly 
)
#unpack data
X=data$X
y=data$Y
g=data$g
K = length(pi)

#Estimate model with model selection
iter = 25000
thin = 1
model = trajMS(X=X, #data matrix
               y=y, #outcomes
               K=K, #number of groups
               time_index=2, #column of X corresponding to time
               iterations=iter, #number of iterations
               thin=thin, #thinning
               dispIter=100) #Print a message every 1000 iterations

#print results
burn = 0.8
summary = summary_single_MS(model,X,y,burn)
print(summary$estimates)