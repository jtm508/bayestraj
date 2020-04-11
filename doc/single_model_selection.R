## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(BayesTraj)

## ------------------------------------------------------------------------
N=1000 #number of units
T=9 #time periods
pi=c(0.5,0.2,0.3) #group membership probabilities
K = length(pi) #number of groups
#coefficients
beta=matrix(c(110,5,-0.5,
              111,-2,0,
              118,0,0),nrow=3,ncol=3,byrow=TRUE)
sigma=2 #standard deviation of outcomes

set.seed(1)
data = gen_data(N=N,
                T=T,
                pi=pi,
                beta=beta,
                sigma=sigma,
                poly = 2 #degree of polynomial
                )

## ------------------------------------------------------------------------
X=data$X
y=data$Y

## ------------------------------------------------------------------------
print(head(X,18))

## ------------------------------------------------------------------------
print(head(y,18))

## ------------------------------------------------------------------------
iter = 5000
thin = 1
model = trajMS(X=X, #data matrix
               y=y, #outcomes
               K=K, #number of groups
               time_index=2, #column of X corresponding to time
               iterations=iter, #number of iterations
               thin=thin, #thinning
               dispIter=1000) #Print a message every 1000 iterations

## ------------------------------------------------------------------------
head(model$beta[[1]]) #group 1's coefficients
head(model$beta[[2]]) #group 2's coefficients
head(model$beta[[3]]) #group 3's coefficients
head(model$sigma) #variance - NOT THE STANDARD DEVIATION
model$c[1:6,1:10] #unit-level group memberships
head(model$pi) #group-membership probabilities

## ------------------------------------------------------------------------
burn = 0.9
summary = summary_single_MS(model,X,y,burn)

## ------------------------------------------------------------------------
print(summary$estimates)

## ------------------------------------------------------------------------
plot(model$beta[[2]][1000:5000,2],type='l')

## ------------------------------------------------------------------------
print(summary$log.likelihood)

## ------------------------------------------------------------------------
plot(model$beta[[1]][1000:5000,1],type='l')

