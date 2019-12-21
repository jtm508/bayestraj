library(BayesTraj)

N=1000
T1=9
T2=9
pi1=c(0.5,0.2,0.3)
pi1_2=matrix(c(0.3,0.3,0.4,
             0.49,0.50,0.01,
             0.7,0.2,0.1),
           nrow=3,ncol=3,byrow=TRUE)
pi12 = matrix(nrow=dim(pi1_2)[1],ncol=dim(pi1_2)[2])
for (i in 1:dim(pi1_2)[1]) {
  for (j in 1:dim(pi1_2)[2]) {
    pi12[i,j] = pi1[i] * pi1_2[i,j]
  }
}
beta1=matrix(c(110,5,-0.5,
               111,-2,0.1,
               118,3,0.1),nrow=3,ncol=3,byrow=TRUE)
beta2=matrix(c(110,6,-0.6,
               111,-3,0.1,
               112,2,0.7),nrow=3,ncol=3,byrow=TRUE)
sigma1=2
sigma2=4

set.seed(3)
data = gen_data2(N=N,
                T1=T1,
                T2=T2,
                pi1=pi1,
                pi2=pi1_2,
                beta1=beta1,
                beta2=beta2,
                sigma1=sigma1,
                sigma2=sigma2,
                poly = 2)
X1=data$X1
X2=data$X2
y1=data$Y1
y2=data$Y2
g1=data$g1
g2=data$g2
K1 = length(pi1)
K2 = dim(pi1_2)[2]
print(table(g1,g2))
print(min(y1))
print(min(y2))
print(max(y1))
print(max(y2))

#create dataset for stata
stata1 = data.frame(X1,y1,'gen1')
colnames(stata1) = c('id','t','t2','y','gen')
stata2 = data.frame(X2,y2,'gen2')
colnames(stata2) = c('id','t','t2','y','gen')
stata = rbind(stata1,stata2)
write.csv(stata,"C:/Users/jtm50/Dropbox/Emma and Justin's Secret Folder/DualTrajectory/Data/stata_simulation.csv",
          row.names=FALSE)

iter = 10000
thin = 1
z1 = matrix(1,nrow=K1,ncol=dim(X1)[2])
z2 = matrix(1,nrow=K2,ncol=dim(X2)[2])
model = dualtraj(X1=X1,
                 X2=X2,
                 y1=y1,
                 y2=y2,
                 K1=K1,
                 K2=K2,
                 z1=z1,
                 z2=z2,
                 iterations=iter,
                 thin=thin,
                 dispIter=100)
burn = 0.9
n1 = dim(model$c1)[1]

summary = summary_dual(model,X1,X2,y1,y2,z1,z2,burn)
print(summary$estimates)
print(summary$log.likelihood)
print(summary$BIC)

library(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

dev.off()

#1: 3, 1, 2
#2: 2, 3, 1

df = model$pi1_2[(n1*(1-burn)+1):n1,,]
dim(df) = c(n1*burn,K1*K2)
df = as.data.frame(df)
names(df) = c('pi11', 'pi21', 'pi31',
              'pi12', 'pi22', 'pi32',
              'pi13', 'pi23', 'pi33')

p11 = ggplot(df, aes(x=pi11)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[3,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p12 = ggplot(df, aes(x=pi12)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[3,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p13 = ggplot(df, aes(x=pi13)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[3,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p21 = ggplot(df, aes(x=pi21)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[1,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p22 = ggplot(df, aes(x=pi22)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[1,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p23 = ggplot(df, aes(x=pi23)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[1,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p31 = ggplot(df, aes(x=pi31)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[2,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p32 = ggplot(df, aes(x=pi32)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[2,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p33 = ggplot(df, aes(x=pi33)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi1_2[2,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()

multiplot(p11, p21, p31,
          p12, p22, p32,
          p13, p23, p33,
          cols=3)

df = model$pi12[(n1*(1-burn)+1):n1,,]
dim(df) = c(n1*burn,K1*K2)
df = as.data.frame(df)
names(df) = c('pi11', 'pi21', 'pi31',
              'pi12', 'pi22', 'pi32',
              'pi13', 'pi23', 'pi33')

p11 = ggplot(df, aes(x=pi11)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[3,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p12 = ggplot(df, aes(x=pi12)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[3,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p13 = ggplot(df, aes(x=pi13)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[3,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p21 = ggplot(df, aes(x=pi21)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[1,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p22 = ggplot(df, aes(x=pi22)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[1,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p23 = ggplot(df, aes(x=pi23)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[1,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p31 = ggplot(df, aes(x=pi31)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[2,3], linetype = "dashed") + 
  labs(y = "Density") + theme_bw()
p32 = ggplot(df, aes(x=pi32)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[2,2], linetype = "dashed") + 
  labs(y = "") + theme_bw()
p33 = ggplot(df, aes(x=pi33)) + geom_density(color="darkblue", fill="lightblue") + geom_vline(xintercept = pi12[2,1], linetype = "dashed") + 
  labs(y = "") + theme_bw()

multiplot(p11, p21, p31,
          p12, p22, p32,
          p13, p23, p33,
          cols=3)


#trace plots to show convergence
par(mfrow = c(2,2))
matplot(c(1:50), model$beta1[[1]][1:50,1], type='l', xlab='Iteration', ylab='Intecept', col=1:5)
matplot(c(1:50), model$beta1[[1]][1:50,2], type='l', xlab='Iteration', ylab='Age', col=1:5)
matplot(c(1:50), model$beta1[[1]][1:50,3], type='l', xlab='Iteration', ylab='Age^2', col=1:5)
matplot(c(1:50), sqrt(model$sigma1[1:50]), type='l', xlab='Iteration', ylab='sigma', col=1:5)

#leave these out, just mention they are similar
par(mfrow = c(3,3))
matplot(c(1:50), model$pi12[1:50,2,1], type='l', xlab='Iteration', ylab='pi_11')
matplot(c(1:50), model$pi12[1:50,2,2], type='l', xlab='Iteration', ylab='pi_12')
matplot(c(1:50), model$pi12[1:50,2,3], type='l', xlab='Iteration', ylab='pi_13')
matplot(c(1:50), model$pi12[1:50,1,1], type='l', xlab='Iteration', ylab='pi_21')
matplot(c(1:50), model$pi12[1:50,1,2], type='l', xlab='Iteration', ylab='pi_22')
matplot(c(1:50), model$pi12[1:50,1,3], type='l', xlab='Iteration', ylab='pi_23')
matplot(c(1:50), model$pi12[1:50,3,1], type='l', xlab='Iteration', ylab='pi_31')
matplot(c(1:50), model$pi12[1:50,3,2], type='l', xlab='Iteration', ylab='pi_32')
matplot(c(1:50), model$pi12[1:50,3,3], type='l', xlab='Iteration', ylab='pi_33')
