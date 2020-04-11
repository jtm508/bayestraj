library(BayesTraj)

N=1000 #number of units
T=9 #time periods
pi=c(0.5,0.2,0.3) #group membership probabilities
#coefficeints
beta=matrix(c(110,5,-0.5,
              111,-2,0.1,
              118,3,0.1),nrow=3,ncol=3,byrow=TRUE)
#standard deviation
sigma=2

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
iter = 5000
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