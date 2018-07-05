#' summary_dual
#'
#' Summarize output of dual trajectory model
#'
#' @param model: List, output from dualtraj or dualtrajMS 
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param y1: Vector, outcomes for series 1
#' @param y2: Vector, outcomes for series 2
#' @param z1: Matrix, K1 x dim(X1)[2] indicator matrix indicating which variables to inlcude in each group.
#' @param z2: Matrix, K2 x dim(X2)[2] indicator matrix indicating which variables to inlcude in each group.
#' @param burn: float, fraction of draws to keep for post burn-in period.
#'
#' @export
 
summary_dual = function(model,X1,X2,y1,y2,z1,z2,burn) {
  n1 = dim(model$c1)[1]
  K1 = length(model$beta1)
  K2 = length(model$beta2)
  
  df = data.frame(A= numeric(0), B= numeric(0), C= numeric(0), D= numeric(0),E= numeric(0))
  #group 1
  beta1.e = matrix(0,nrow=K1,ncol=max(do.call(rbind,lapply(model$beta1,dim))[,2]))
  for (k in 1:K1) {
    A = cbind(colMeans(tail(model$beta1[[k]],n1*burn)),
              apply(tail(model$beta1[[k]],n1*burn), 2, sd),
              t(apply(tail(model$beta1[[k]],n1*burn),2,quantile,probs=c(0.025,0.5,0.975))))
    beta1.e[k,z1[k,]==1] = A[,1]
    rownames(A) = sprintf(paste("beta1_",k,"[%d]",sep=''),seq(1:dim(model$beta1[[k]])[2]))
    df = rbind(df,A)
  }
  #apply(tail(model$beta1[[1]],n1*burn),2,quantile,probs=c(0.025,0.5,0.975))
  
  A = as.data.frame(t(c(mean(tail(sqrt(model$sigma1),n1*burn)),
                        sd(tail(sqrt(model$sigma1),n1*burn)),
                        quantile(tail(sqrt(model$sigma1),n1*burn),probs=c(0.025,0.5,0.975)))))
  sigma1.e = as.numeric(A[1])^2
  rownames(A) = "sigma1"
  df = rbind(df,A)
  
  A = cbind(colMeans(tail(100*model$pi1,n1*burn)),
            apply(tail(100*model$pi1,n1*burn),2,sd),
            t(apply(tail(100*model$pi1,n1*burn),2,quantile,probs=c(0.025,0.5,0.975))))
  pi1.e = A[,1]/100
  rownames(A) = sprintf("pi1[%d]",seq(1:K1))
  df = rbind(df,A)
  
  #group 2
  beta2.e = matrix(0,nrow=K2,ncol=max(do.call(rbind,lapply(model$beta2,dim))[,2]))
  for (k in 1:K2) {
    A = cbind(colMeans(tail(model$beta2[[k]],n1*burn)),
              apply(tail(model$beta2[[k]],n1*burn), 2, sd),
              t(apply(tail(model$beta2[[k]],n1*burn),2,quantile,probs=c(0.025,0.5,0.975))))
    beta2.e[k,z2[k,]==1] = A[,1]
    rownames(A) = sprintf(paste("beta2_",k,"[%d]",sep=''),seq(1:dim(model$beta2[[k]])[2]))
    df = rbind(df,A)
  }
  
  A = as.data.frame(t(c(mean(tail(sqrt(model$sigma2),n1*burn)),
                        sd(tail(sqrt(model$sigma2),n1*burn)),
                        quantile(tail(sqrt(model$sigma2),n1*burn),probs=c(0.025,0.5,0.975)))))
  sigma2.e = as.numeric(A[1])^2
  rownames(A) = "sigma2"
  df = rbind(df,A)
  
  A = cbind(colMeans(tail(100*model$pi2,n1*burn)),
            apply(tail(100*model$pi2,n1*burn),2,sd),
            t(apply(tail(100*model$pi2,n1*burn),2,quantile,probs=c(0.025,0.5,0.975))))
  rownames(A) = sprintf("pi2[%d]",seq(1:K2))
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975))))
  pi1_2.e = matrix(A[,1],nrow=K1,ncol=K2,byrow=TRUE)/100
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K1) {
    for (k2 in 1:K2) {
      rn = c(rn,paste("1->2,",k1,"->",k2,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), mean)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), sd)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975)))
  colnames(A) = colnames(df)
  rn = c()
  for (k2 in 1:K2) {
    for (k1 in 1:K1) {
      rn = c(rn,paste("2->1,",k2,"->",k1,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #joint
  A = cbind(as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975))))
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K1) {
    for (k2 in 1:K2) {
      rn = c(rn,paste(k1,k2,sep=','))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  colnames(df)[1:2] = c("Estimate","Standard Deviation")
  
  id1 = X1[,1]
  id2 = X2[,1]
  X1[,1]=1
  X2[,1]=1
  ll = log_lik_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta1.e,beta2.e,sigma1.e,sigma2.e,id1,id2)
  BIC = BIC_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta1.e,beta2.e,sigma1.e,sigma2.e,id1,id2,z1,z2)
  
  return(list(estimates = df,
              log.likelihood = ll,
              BIC = BIC))
}

#' summary_dual_constrained
#'
#' Summarize output of dual trajectory model
#'
#' @param model: List, output from dualtraj or dualtrajMS 
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param y1: Vector, outcomes for series 1
#' @param y2: Vector, outcomes for series 2
#' @param z: Matrix, K x dim(X1)[2] indicator matrix indicating which variables to inlcude
#' @param burn: float, fraction of draws to keep for post burn-in period.
#'
#' @export

summary_dual_constrained = function(model,X1,X2,y1,y2,z,burn) {
  n = dim(model$c1)[1]
  K = length(model$beta)
  
  df = data.frame(A= numeric(0), B= numeric(0), C= numeric(0), D= numeric(0),E= numeric(0))
  #group 1
  beta.e = matrix(0,nrow=K,ncol=max(do.call(rbind,lapply(model$beta,dim))[,2]))
  for (k in 1:K) {
    A = cbind(colMeans(tail(model$beta[[k]],n*burn)),
              apply(tail(model$beta[[k]],n*burn), 2, sd),
              t(apply(tail(model$beta[[k]],n*burn),2,quantile,probs=c(0.025,0.5,0.975))))
    beta.e[k,z[k,]==1] = A[,1]
    rownames(A) = sprintf(paste("beta_",k,"[%d]",sep=''),seq(1:dim(model$beta[[k]])[2]))
    df = rbind(df,A)
  }
  
  A = as.data.frame(t(c(mean(tail(sqrt(model$sigma),n*burn)),
                        sd(tail(sqrt(model$sigma),n*burn)),
                        quantile(tail(sqrt(model$sigma),n*burn),probs=c(0.025,0.5,0.975)))))
  sigma.e = as.numeric(A[1])^2
  rownames(A) = "sigma"
  df = rbind(df,A)
  
  A = cbind(colMeans(tail(100*model$pi1,n*burn)),
            apply(tail(100*model$pi1,n*burn),2,sd),
            t(apply(tail(100*model$pi1,n*burn),2,quantile,probs=c(0.025,0.5,0.975))))
  pi1.e = A[,1]/100
  rownames(A) = sprintf("pi1[%d]",seq(1:K))
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(t(apply(100*model$pi1_2[(n*(1-burn)+1):n,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi1_2[(n*(1-burn)+1):n,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi1_2[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi1_2[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi1_2[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.975))))
  pi1_2.e = matrix(A[,1],nrow=K,ncol=K,byrow=TRUE)/100
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      rn = c(rn,paste("1->2,",k1,"->",k2,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(apply(100*model$pi2_1[(n*(1-burn)+1):n,,], c(2,3), mean)),
            as.vector(apply(100*model$pi2_1[(n*(1-burn)+1):n,,], c(2,3), sd)),
            as.vector(apply(100*model$pi2_1[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.025)),
            as.vector(apply(100*model$pi2_1[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.5)),
            as.vector(apply(100*model$pi2_1[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.975)))
  colnames(A) = colnames(df)
  rn = c()
  for (k2 in 1:K) {
    for (k1 in 1:K) {
      rn = c(rn,paste("2->1,",k2,"->",k1,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #joint
  A = cbind(as.vector(t(apply(100*model$pi12[(n*(1-burn)+1):n,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi12[(n*(1-burn)+1):n,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi12[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi12[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi12[(n*(1-burn)+1):n,,], c(2,3), quantile,probs=0.975))))
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      rn = c(rn,paste(k1,k2,sep=','))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  colnames(df)[1:2] = c("Estimate","Standard Deviation")
  
  id1 = X1[,1]
  id2 = X2[,1]
  X1[,1]=1
  X2[,1]=1
  ll = log_lik_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta.e,beta.e,sigma.e,sigma.e,id1,id2)
  BIC = BIC_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta.e,beta.e,sigma.e,sigma.e,id1,id2,z,z,constrain=TRUE)
  
  return(list(estimates = df,
              log.likelihood = ll,
              BIC = BIC))
}


#' summary_dual_MS
#'
#' Summarize output of dual trajectory model with Bayesian model averaging.
#'
#' @param model: List, output from dualtraj or dualtrajMS 
#' @param X1: Matrix, design matrix for series 1. 1st column should be the id.
#' @param X2: Matrix, design matrix for series 2. 1st column should be the id.
#' @param y1: Vector, outcomes for series 1
#' @param y2: Vector, outcomes for series 2
#' @param burn: float, fraction of draws to keep for post burn-in period.
#'
#' @export

summary_dual_MS = function(model,X1,X2,y1,y2,burn) {
  n1 = dim(model$c1)[1]
  K1 = length(model$beta1)
  K2 = length(model$beta2)
  
  df = data.frame(A= numeric(0), B= numeric(0), C= numeric(0), D= numeric(0),E= numeric(0),F= numeric(0))
  #group 1
  beta1.e = matrix(0,nrow=K1,ncol=dim(X1)[2])
  sigma1.e = rep(NA,K1)
  for (k in 1:K1) {
    A = cbind(colMeans(tail(model$beta1[[k]],n1*burn)),
              apply(tail(model$beta1[[k]],n1*burn), 2, sd),
              t(apply(tail(model$beta1[[k]],n1*burn),2,quantile,probs=c(0.025,0.5,0.975))),
              colMeans(model$z1[[k]]))
    beta1.e[k,] = A[,4]
    rownames(A) = sprintf(paste("beta1_",k,"[%d]",sep=''),seq(1:dim(model$beta1[[k]])[2]))
    df = rbind(df,A)
    
    A = as.data.frame(t(c(mean(tail(sqrt(model$sigma1[,k]),n1*burn)),
                          sd(tail(sqrt(model$sigma1[,k]),n1*burn)),
                          quantile(tail(sqrt(model$sigma1[,k]),n1*burn),probs=c(0.025,0.5,0.975)),
                          NA)))
    sigma1.e[k] = as.numeric(A[1])^2
    rownames(A) = paste("sigma1_",k,sep='')
    df = rbind(df,A)
  }
  
  A = cbind(colMeans(tail(100*model$pi1,n1*burn)),
            apply(tail(100*model$pi1,n1*burn),2,sd),
            t(apply(tail(100*model$pi1,n1*burn),2,quantile,probs=c(0.025,0.5,0.975))),
            NA)
  pi1.e = A[,1]/100
  rownames(A) = sprintf("pi1[%d]",seq(1:K1))
  df = rbind(df,A)
  
  #group 2
  beta2.e = matrix(0,nrow=K2,ncol=dim(X2)[2])
  sigma2.e = rep(NA,K2)
  for (k in 1:K2) {
    A = cbind(colMeans(tail(model$beta2[[k]],n1*burn)),
              apply(tail(model$beta2[[k]],n1*burn), 2, sd),
              t(apply(tail(model$beta2[[k]],n1*burn),2,quantile,probs=c(0.025,0.5,0.975))),
              colMeans(model$z2[[k]]))
    beta2.e[k,] = A[,4]
    rownames(A) = sprintf(paste("beta2_",k,"[%d]",sep=''),seq(1:dim(model$beta2[[k]])[2]))
    df = rbind(df,A)
    
    A = as.data.frame(t(c(mean(tail(sqrt(model$sigma2[,k]),n1*burn)),
                          sd(tail(sqrt(model$sigma2[,k]),n1*burn)),
                          quantile(tail(sqrt(model$sigma2[,k]),n1*burn),probs=c(0.025,0.5,0.975)),
                          NA)))
    sigma2.e[k] = as.numeric(A[1])^2
    rownames(A) = paste("sigma2_",k,sep='')
    df = rbind(df,A)
  }
  
  A = cbind(colMeans(tail(100*model$pi2,n1*burn)),
            apply(tail(100*model$pi2,n1*burn),2,sd),
            t(apply(tail(100*model$pi2,n1*burn),2,quantile,probs=c(0.025,0.5,0.975))),
            NA)
  rownames(A) = sprintf("pi2[%d]",seq(1:K2))
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi1_2[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975))),
            NA)
  pi1_2.e = matrix(A[,1],nrow=K1,ncol=K2,byrow=TRUE)/100
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K1) {
    for (k2 in 1:K2) {
      rn = c(rn,paste("1->2,",k1,"->",k2,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #transitions
  A = cbind(as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), mean)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), sd)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5)),
            as.vector(apply(100*model$pi2_1[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975)),
            NA)
  colnames(A) = colnames(df)
  rn = c()
  for (k2 in 1:K2) {
    for (k1 in 1:K1) {
      rn = c(rn,paste("2->1,",k2,"->",k1,sep=''))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  #joint
  A = cbind(as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), mean))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), sd))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.025))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.5))),
            as.vector(t(apply(100*model$pi12[(n1*(1-burn)+1):n1,,], c(2,3), quantile,probs=0.975))),
            NA)
  colnames(A) = colnames(df)
  rn = c()
  for (k1 in 1:K1) {
    for (k2 in 1:K2) {
      rn = c(rn,paste(k1,k2,sep=','))
    }
  }
  rownames(A) = rn
  df = rbind(df,A)
  
  colnames(df)[1:2] = c("Estimate","Standard Deviation")
  colnames(df)[6] = "Inclusion Prob."
  
  id1 = X1[,1]
  id2 = X2[,1]
  X1[,1]=1
  X2[,1]=1
  z1 = (beta1.e != 0)*1
  z2 = (beta2.e != 0)*1
  
  ll = log_lik_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta1.e,beta2.e,sigma1.e,sigma2.e,id1,id2)
  BIC = BIC_dual(X1,X2,y1,y2,pi1.e,pi1_2.e,beta1.e,beta2.e,sigma1.e,sigma2.e,id1,id2,z1,z2)
  
  return(list(estimates = df,
              log.likelihood = ll,
              BIC = BIC))
}