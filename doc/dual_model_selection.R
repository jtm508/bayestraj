## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(BayesTraj)

## ------------------------------------------------------------------------
N=1000 #number of units
T1=9 #Time periods for Group 1
T2=9 #Time periods for Group 2
pi1=c(0.5,0.2,0.3) #Group 1 membership probabilities
#Transition Matrix
pi1_2=matrix(c(0.3,0.3,0.4,
               0.2,0.5,0.3,
               0.7,0.2,0.1),
             nrow=3,ncol=3,byrow=TRUE)
K1 = length(pi1) #Number of groups in series 1
K2 = dim(pi1_2)[2] #Number of groups in series 2
#Coefficients for Series 1
beta1=matrix(c(100,10,-0.5,0.1,
               80,-1,0.5,0,
               120,20,-2,0),nrow=3,ncol=4,byrow=TRUE)
#Coefficients for Series 2
beta2=matrix(c(90,0.4,0,0,
               50,1,-0.5,0,
               100,-30,3,-0.3),nrow=3,ncol=4,byrow=TRUE)
sigma1=16 #standard deviation of Series 1 outcomes
sigma2=36 #standard deviation of Series 2 outcomes

set.seed(1)
data = gen_data_dual(N=N,
                T1=T1,
                T2=T2,
                pi1=pi1,
                pi2=pi1_2,
                beta1=beta1,
                beta2=beta2,
                sigma1=sigma1,
                sigma2=sigma2,
                poly = 3) #degree of polynomial

## ------------------------------------------------------------------------
X1=data$X1
X2=data$X2
y1=data$Y1
y2=data$Y2

## ------------------------------------------------------------------------
print(head(X1,18))

## ------------------------------------------------------------------------
print(head(y1,18))

## ------------------------------------------------------------------------
iter = 5000
thin = 1
model = dualtrajMS(X1=X1, #data matrix Series 1
                 X2=X2, #data matrix Series 2
                 y1=y1, #outcomes Series 1
                 y2=y2, #outcomes Series 2
                 K1=K1, #number of groups Series 1
                 K2=K2, #number of groups Series 2
                 time_index=2, #column of X corresponding to time
                 iterations=iter, #number of iterations
                 thin=thin, #thinning
                 dispIter=1000) #Print a message every 1000 iterations

## ------------------------------------------------------------------------
head(model$beta1[[1]]) #Series 1 group 1's coefficients
head(model$beta1[[2]]) #Series 1 group 2's coefficients
head(model$beta1[[3]]) #Series 1 group 3's coefficients
head(model$sigma1) #Series 1 variance - NOT THE STANDARD DEVIATION
model$c1[1:6,1:10] #Series 1 unit-level group memberships
head(model$pi1) #Series 1  group-membership probabilities
head(model$pi2) #Series 2  group-membership probabilities
model$pi1_2[1,,] #Transition probabilities from Series 1 Group 1.
model$pi12[1,,] #Joint probability of both Series group memberships

## ------------------------------------------------------------------------
burn = 0.9
summary = summary_dual_MS(model,X1,X2,y1,y2,burn)

## ------------------------------------------------------------------------
print(summary$estimates)

## ------------------------------------------------------------------------
plot(model$beta1[[3]][1000:5000,4],type='l')

## ------------------------------------------------------------------------
print(summary$log.likelihood)

## ------------------------------------------------------------------------
plot(model$beta1[[1]][1000:5000,1],type='l')

