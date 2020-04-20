library(BayesTraj)

sim_number = 1

N=1000 #number of units
T=9 #time periods
if (sim_number==1) {
  pi=c(0.25,0.25,0.5) #group membership probabilities
  #coefficeints
  beta=matrix(c(20,1,-0.3,
                25,-1,0.2,
                30,2,0.1),nrow=3,ncol=3,byrow=TRUE)
  #variance
  sigma=1.5^2
} else if (sim_number==2) {
  pi=c(0.3,0.4,0.3) #group membership probabilities
  #coefficeints
  beta=matrix(c(5,1,-0.1,
                4,0.3,0.1,
                5,-1,0.1),nrow=3,ncol=3,byrow=TRUE)
  #variance
  sigma=1^2
} 

#generate data
set.seed(1)
data = gen_data(N=N,
                T=T,
                pi=pi,
                beta=beta,
                sigma=sigma,
                poly = 2)
#unpack data
X=data$X
y=data$Y
g=data$g
K = length(pi)

#estimate model
iter = 10000
thin = 1
z = matrix(1,nrow=K,ncol=dim(X)[2])
model = traj(X=X,
             y=y,
             K=K,
             z=z,
             iterations=iter,
             thin=thin,
             dispIter=1000,
             ll=TRUE)
burn = 0.9
n1 = dim(model$c1)[1]

#print output
summary = summary_single(model,X,y,z,burn)
print(summary$estimates)
print(summary$log.likelihood)
print(summary$BIC)