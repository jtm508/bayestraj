library(BayesTraj)
set.seed(2)

N=1000 #number of units
T1=8 #series 1 time periods
T2=9 #series 2 time periods
pi1=c(0.5,0.2,0.3) #group membership probabilities
#transition probabilities
pi1_2=matrix(c(0.3,0.3,0.4,
               0.2,0.5,0.3,
               0.7,0.2,0.1),
             nrow=3,ncol=3,byrow=TRUE)
#series 1 coefficients
beta1=matrix(c(100,4,10,-0.5,0.1,
               80,8,-1,0.5,0,
               120,-6,20,-2,0),nrow=3,ncol=5,byrow=TRUE)
#series 2 coefficients
beta2=matrix(c(90,6,0.4,0,0,
               50,0,1,-0.5,0,
               100,-3,-30,3,-0.3),nrow=3,ncol=5,byrow=TRUE)
sigma1=16 #series 1 variance
sigma2=36 #series 2 variance

#generate data
data = gen_data_dual(N=N,
                 T1=T1,
                 T2=T2,
                 pi1=pi1,
                 pi2=pi1_2,
                 beta1=beta1,
                 beta2=beta2,
                 sigma1=sigma1,
                 sigma2=sigma2,
                 poly=3)
#unpack data
X1=data$X1
X2=data$X2
y1=data$Y1
y2=data$Y2
g1=data$g1
g2=data$g2
K1 = length(pi1)
K2 = dim(pi1_2)[2]

#estimate model with model selection#
model = dualtrajMS(X1=X1,
                   X2=X2,
                   y1=y1,
                   y2=y2,
                   K1=3,
                   K2=3,
                   time_index=3,
                   iterations=5000,
                   thin=1,
                   dispIter=10)

#print results
burn = 0.9
summary = summary_dual_MS(model,X1,X2,y1,y2,burn)
print(summary$estimates)