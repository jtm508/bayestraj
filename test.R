N=1000
T1=5
T2=6
pi1=c(0.5,0.2,0.3)
pi2=matrix(c(0.3,0.3,0.4,
             0.2,0.5,0.3,
             0.7,0.2,0.1),
           nrow=3,ncol=3,byrow=TRUE)
beta1=matrix(c(1,1,6,-3,
               3,2,2,-1,
               2,1,9,-3),nrow=3,ncol=4,byrow=TRUE)
beta2=matrix(c(4,1,5,-2,
               6,2,5,-4,
               5,1,15,-2),nrow=3,ncol=4,byrow=TRUE)
sigma1=2
sigma2=4

data = gen_data2(N=N,
                T1=T1,
                T2=T2,
                pi1=pi1,
                pi2=pi2,
                beta1=beta1,
                beta2=beta2,
                sigma1=sigma1,
                sigma2=sigma2)
X1=data$X1
X2=data$X2
Y1=data$Y1
Y2=data$Y2
g1=data$g1
g2=data$g2
K1 = length(pi1)
K2 = dim(pi2)[2]

model = dualtraj(X1=X1,
                 X2=X2,
                 Y1=Y1,
                 Y2=Y2,
                 g1=g1,
                 g2=g2,
                 K1=3,
                 K2=3,
                 iterations=10000,
                 thin=10,
                 dispIter=100)