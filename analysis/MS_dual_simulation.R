library(BayesTraj)
set.seed(2)

sim_number = 1

N=1000 #number of units
T1=8 #series 1 time periods
T2=9 #series 2 time periods

if (sim_number==1) {
  pi1=c(0.5,0.2,0.3) #group membership probabilities
  #transition probabilities
  pi1_2=matrix(c(0.3,0.3,0.4,
                 0.2,0.5,0.3,
                 0.7,0.2,0.1),
               nrow=3,ncol=3,byrow=TRUE)
  #joint probability calculations
  pi12 = matrix(nrow=dim(pi1_2)[1],ncol=dim(pi1_2)[2])
  for (i in 1:dim(pi1_2)[1]) {
    for (j in 1:dim(pi1_2)[2]) {
      pi12[i,j] = pi1[i] * pi1_2[i,j]
    }
  }
  pi2 = colSums(pi12)
  pi2_1 = matrix(nrow=length(pi2),ncol=length(pi1))
  for (i in 1:length(pi2)) {
    for (j in 1:length(pi1)) {
      pi2_1[i,j] = pi12[j,i] / pi2[i]
    }
  }
  #series 1 coefficients
  beta1=matrix(c(100,10,-0.5,0.1,
                 80,-1,0.5,0,
                 120,20,-2,0),nrow=3,ncol=4,byrow=TRUE)
  #series 2 coefficients
  beta2=matrix(c(90,0.4,0,0,
                 50,1,-0.5,0,
                 100,-30,3,-0.3),nrow=3,ncol=4,byrow=TRUE)
  sigma1=4^2 #series 1 variance
  sigma2=6^2 #series 2 variance
  poly = 3 #degree of polynomial
} else if (sim_number==2) {
  pi1=c(0.2,0.15,0.65) #group membership probabilities
  #transition probabilities
  pi1_2=matrix(c(0.4,0.3,0.4,
                 0.5,0.1,0.2,
                 0.3,0.5,0.2),
               nrow=3,ncol=3,byrow=TRUE)
  #series 1 coefficients
  beta1=matrix(c(25,10,0,
                 10,-1,0.5,
                 30,20,0),nrow=3,ncol=3,byrow=TRUE)
  #series 2 coefficients
  beta2=matrix(c(25,1,0,
                 25,0.5,-0.2,
                 20,0,0),nrow=3,ncol=3,byrow=TRUE)
  sigma1=2^2 #series 1 variance
  sigma2=3^2 #series 2 variance
  poly = 2 #degree of polynomial
}

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
                 poly=poly)
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
                   time_index=2,
                   iterations=10000,
                   thin=1,
                   dispIter=100)

#print results
burn = 0.9
summary = summary_dual_MS(model,X1,X2,y1,y2,burn)
print(summary$estimates)