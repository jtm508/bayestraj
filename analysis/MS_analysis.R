library(ggplot2)  
library(reshape2)
library(BayesTraj)

#The loop is for seraching for a seed which finds a global mode. Ignore unless redoing this search
#for (i in 1:50){
#print(i)
#set.seed(i)
set.seed(1)

#Son data
son_jtm = read.csv("C:/Users/jtm50/Dropbox/Emma and Justin's Secret Folder/DualTrajectory/Data/offspring_f.csv")
son_jtm = son_jtm[son_jtm['male']==1,]
son_jtm = son_jtm[c('pid','age','loginc')]
son_jtm['age'] = son_jtm['age'] - 25
son_jtm = son_jtm[complete.cases(son_jtm), ] # drop missing observations

#Father data
father_jtm = read.csv("C:/Users/jtm50/Dropbox/Emma and Justin's Secret Folder/DualTrajectory/Data/father.csv")
father_jtm = father_jtm[c('pid','agef','logincf')]
father_jtm['agef'] = father_jtm['agef'] - 25
father_jtm = father_jtm[complete.cases(father_jtm), ] # drop missing observations

#only keep common pid's
son_pid = as.numeric(unlist(unique(son_jtm['pid'])))
father_pid = as.numeric(unlist(unique(father_jtm['pid'])))
common_id = intersect(son_pid,father_pid)
son_jtm = son_jtm[as.numeric(unlist(son_jtm['pid'])) %in% common_id,]
father_jtm = father_jtm[as.numeric(unlist(father_jtm['pid'])) %in% common_id,]

#Format Data
Y_son = as.vector(son_jtm$loginc)
X_son = as.matrix(son_jtm[c('pid','age')])
X_son = cbind(X_son,X_son[,2]^2,X_son[,2]^3)
rm(son_jtm)

Y_father = as.vector(father_jtm$logincf)
X_father = as.matrix(father_jtm[c('pid','agef')])
X_father = cbind(X_father,X_father[,2]^2,X_father[,2]^3)
rm(father_jtm)

iter=25000
K = 5
K1 = K
K2 = K
tryCatch({
model = dualtrajMS(X1=X_father,
                 X2=X_son,
                 y1=Y_father,
                 y2=Y_son,
                 K1=K1,
                 K2=K2,
                 time_index=2,
                 iterations=iter,
                 thin=1,
                 dispIter=100)

burn = 0.8
summary = summary_dual_MS(model,X_father,X_son,Y_father,Y_son,burn)
print(summary$estimates)
print(summary$log.likelihood)
print(summary$BIC)
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#}

n1 = dim(model$c1)[1]

#function for putting multiple plots in one image
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

#posterior predictive distributions for plotting trajectories
sampleTraj = function(beta,bounds,poly,log=FALSE,div=FALSE) {
  #generate covariates for each time period
  cov = rep(1,(bounds[2]-bounds[1]+1))
  for (i in 1:poly) {
    cov = cbind(cov,(0:30)^i)
  }
  if (log == TRUE) {
    cov = cbind(cov,(log(cov[,2]+1)))
  }
  if (div == TRUE) {
    cov = cbind(cov,(1/(cov[,2]+1)))
  }
  
  #simulate y using cov*beta + n(0,sqrt(sigma))
  yTraj= list()
  for (i in 1:length(beta)) {
    yTraj[[i]] = t(cov %*% t(tail(beta[[i]],n1*burn)))
  }
  return(yTraj)
}

sampleTraj = function(beta,bounds,poly,log=FALSE,div=FALSE) {
  #generate covariates for each time period
  cov = rep(1,(bounds[2]-bounds[1]+1))
  for (i in 1:poly) {
    cov = cbind(cov,(0:30)^i)
  }
  if (log == TRUE) {
    cov = cbind(cov,(log(cov[,2]+1)))
  }
  if (div == TRUE) {
    cov = cbind(cov,(1/(cov[,2]+1)))
  }
  
  #simulate y using cov*beta + n(0,sqrt(sigma))
  yTraj= list()
  for (i in 1:length(beta)) {
    ppd = t(cov %*% t(tail(beta[[i]],n1*burn)))
    df = exp(as.data.frame(t(apply(ppd,2,quantile,probs=c(0.025,0.5,0.975)))))
    colnames(df) = c('lower','median','upper')
    df$time = seq.int(nrow(df)) + 24
    df$group = as.factor(i)
    yTraj[[i]] = df
  }
  #append datasets
  df_all = do.call(rbind,yTraj)
  
  #reorder groups by mean value
  df_group = aggregate(df_all$median,list(df_all$group),mean)
  df_group = df_group[order(df_group$x,decreasing=TRUE),]
  df_group$order = seq.int(nrow(df_group))
  df_group = df_group[order(df_group$Group.1),]
  group_order = df_group$order
  df_all$group = as.factor(group_order[df_all$group])
  df_all = df_all[order(df_all$group),]
  return(df_all)
}

#get trajectory samples
yTraj1 = sampleTraj(model$beta1,bounds=c(0,30),poly=3,log=FALSE,div=FALSE)
yTraj2 = sampleTraj(model$beta2,bounds=c(0,30),poly=3,log=FALSE,div=FALSE)

#function to plot the trajectories
trajPlot = function(yTraj,bounds,title) {
  order = order(unlist(lapply(yTraj,mean)))
  print(length(order))
  obj.list = list()
  cols = c()
  for (j in 1:length(order)){
    obj.list[[j]] = as.data.frame(t(apply(yTraj[[order[j]]],2,quantile,probs=c(0.025,0.5,0.975))))
    cols = c(cols,c(paste("A",j,sep=""),paste("B",j,sep=""),paste("C",j,sep="")))
  }
  predictiveIntervals = exp(do.call(cbind, obj.list))  / 1000
  colnames(predictiveIntervals) = cols
  predictiveIntervals$index = seq(bounds[1],bounds[2]) + 25
  
  p = ggplot(data = predictiveIntervals)
  p = p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank())
  p = p + scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  p = p + expand_limits(y=0)
  col_list = c('dodgerblue4','darkgreen','darkorange3','goldenrod4','darkred','turquoise2','palegreen1')
  pch_list = c(15,16,17,18,8,4,3)
  for (j in 1:length(order)){
    p = p +
      geom_ribbon(aes_string(x = 'index', ymin=paste("A",j,sep=""), ymax=paste("C",j,sep="")), 
                  fill='gray',alpha = 0.5) +
      geom_line(aes_string(x = 'index', y = paste("B",j,sep="")), color = col_list[j]) + 
      geom_point(aes_string(x = 'index', y = paste("B",j,sep="")), 
                 pch = pch_list[j],size=2,fill='blue',color=col_list[j])
  }
  p = p +
    labs(title = title) + 
    xlab("Age") + 
    ylab("Income (Thousands)")
  p
}

trajPlot = function(yTraj,title,names=NULL) {
  col_list = c('dodgerblue4','darkgreen','darkorange3','goldenrod4','darkred','turquoise2','palegreen1')
  pch_list = c(15,16,17,18,8,4,3)
  p = ggplot(yTraj, aes(x=time, y=median, group=group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), fill='gray',alpha = 0.5) +
    geom_line(color=col_list[yTraj$group]) +
    geom_point(pch=pch_list[yTraj$group],color=col_list[yTraj$group],size=2) +
    theme(legend.title = element_blank(),legend.position="bottom",panel.background = element_rect(fill = 'white', colour = 'white')) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    expand_limits(y=0) +
    labs(title = title) + 
    xlab("Age") + 
    ylab("Income (Thousands)")
  p
}
trajPlot(yTraj1,'Father Trajectories')

#function to plot the trajectories
trajPlot = function(yTraj,title,names=NULL) {
  group_num = yTraj$group
  if (!is.null(names)) {
    yTraj$group = names[yTraj$group]
  }
  col_list = c('dodgerblue4','darkgreen','darkorange3','goldenrod4','darkred','turquoise2','palegreen1')
  pch_list = c(15,16,17,18,8,4,3)
  p = ggplot(yTraj, aes(x=time, y=median, group=group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), fill='gray',alpha = 0.5) +
    geom_line(aes(color=group)) +
    #geom_point(aes(linetype = group),pch=pch_list[yTraj$group],color=col_list[yTraj$group],size=2) +
    geom_point(aes(shape=group,color=group)) +
    theme(legend.title = element_blank(),legend.position="bottom",panel.background = element_rect(fill = 'white', colour = 'white')) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    expand_limits(y=0) +
    labs(title = title) + 
    xlab("Age") + 
    ylab("Income (Thousands)")
  p
}
trajPlot(yTraj1,'Father Trajectories',c('Test','sdfsd','sdfs','tada','sdf'))

#combine plots
multiplot(trajPlot(yTraj1,'Father Trajectories',c('Test','sdfsd','sdfs','tada','sdf')),
          trajPlot(yTraj2,'Son Trajectories'),
          cols=2)

#group membership plots
multiDensity = function(x) {
  data =  melt(as.data.frame(x))
  ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
}
multiDensity(tail(model$pi1,n1*burn))
multiDensity(tail(model$pi2,n1*burn))

#joint distribution plots
joint_dist_plot = function(pi1,pi2,pi12) {
  df1 = as.data.frame(tail(pi1,n1*burn))
  names(df1) = c('pi1_1','pi1_2','pi1_3','pi1_4','pi1_5')
  
  p1_1 = ggplot(df1, aes(x=pi1_1)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p1_2 = ggplot(df1, aes(x=pi1_2)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p1_3 = ggplot(df1, aes(x=pi1_3)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p1_4 = ggplot(df1, aes(x=pi1_4)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p1_5 = ggplot(df1, aes(x=pi1_5)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  
  df2 = as.data.frame(tail(pi2,n1*burn))
  names(df2) = c('pi2_1','pi2_2','pi2_3','pi2_4','pi2_5')
  
  p2_1 = ggplot(df2, aes(x=pi2_1)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p2_2 = ggplot(df2, aes(x=pi2_2)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p2_3 = ggplot(df2, aes(x=pi2_3)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p2_4 = ggplot(df2, aes(x=pi2_4)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  p2_5 = ggplot(df2, aes(x=pi2_5)) + geom_density(color="darkblue", fill="lightblue") + xlim(0,0.35)
  
  df12 = aperm(pi12[(n1*(1-burn)+1):n1,,],c(1,3,2))
  dim(df12) = c(n1*burn,K1*K2)
  df12 = as.data.frame(tail(df12,n1*burn))
  names(df12) = c('pi11','pi12','pi13','pi14','pi15',
                  'pi21','pi22','pi23','pi24','pi25',
                  'pi31','pi32','pi33','pi34','pi35',
                  'pi41','pi42','pi43','pi44','pi45',
                  'pi51','pi52','pi53','pi54','pi55')
  p11 = ggplot(df12, aes(x=pi11)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p12 = ggplot(df12, aes(x=pi12)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p13 = ggplot(df12, aes(x=pi13)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p14 = ggplot(df12, aes(x=pi14)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p15 = ggplot(df12, aes(x=pi15)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p21 = ggplot(df12, aes(x=pi21)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p22 = ggplot(df12, aes(x=pi22)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p23 = ggplot(df12, aes(x=pi23)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p24 = ggplot(df12, aes(x=pi24)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p25 = ggplot(df12, aes(x=pi25)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p31 = ggplot(df12, aes(x=pi31)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p32 = ggplot(df12, aes(x=pi32)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p33 = ggplot(df12, aes(x=pi33)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p34 = ggplot(df12, aes(x=pi34)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p35 = ggplot(df12, aes(x=pi35)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p41 = ggplot(df12, aes(x=pi41)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p42 = ggplot(df12, aes(x=pi42)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p43 = ggplot(df12, aes(x=pi43)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p44 = ggplot(df12, aes(x=pi44)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p45 = ggplot(df12, aes(x=pi45)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p51 = ggplot(df12, aes(x=pi51)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p52 = ggplot(df12, aes(x=pi52)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p53 = ggplot(df12, aes(x=pi53)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p54 = ggplot(df12, aes(x=pi54)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  p55 = ggplot(df12, aes(x=pi55)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,0.12)
  
  blank = plot.new()
  multiplot(p1_1, p1_2, p1_3, p1_4, p1_5, blank,
            p11, p21, p31, p41, p51, p2_1,
            p12, p22, p32, p42, p52, p2_2,
            p13, p23, p33, p43, p53, p2_3,
            p14, p24, p34, p44, p54, p2_4,
            p15, p25, p35, p45, p55, p2_4,
            cols=6)
}

joint_dist_plot(model$pi1,model$pi2,model$pi12)

#transition probabilit plots
trans_plot = function(pi) {
  df1_2 = aperm(pi[(n1*(1-burn)+1):n1,,],c(1,3,2))
  dim(df1_2) = c(n1*burn,K1*K2)
  df1_2 = as.data.frame(tail(df1_2,n1*burn))
  names(df1_2) = c('pi1_1','pi1_2','pi1_3','pi1_4','pi1_5',
                   'pi2_1','pi2_2','pi2_3','pi2_4','pi2_5',
                   'pi3_1','pi3_2','pi3_3','pi3_4','pi3_5',
                   'pi4_1','pi4_2','pi4_3','pi4_4','pi4_5',
                   'pi5_1','pi5_2','pi5_3','pi5_4','pi5_5')
  xl = 0.5
  p11 = ggplot(df1_2, aes(x=pi1_1)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p12 = ggplot(df1_2, aes(x=pi1_2)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p13 = ggplot(df1_2, aes(x=pi1_3)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p14 = ggplot(df1_2, aes(x=pi1_4)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p15 = ggplot(df1_2, aes(x=pi1_5)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p21 = ggplot(df1_2, aes(x=pi2_1)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p22 = ggplot(df1_2, aes(x=pi2_2)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p23 = ggplot(df1_2, aes(x=pi2_3)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p24 = ggplot(df1_2, aes(x=pi2_4)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p25 = ggplot(df1_2, aes(x=pi2_5)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p31 = ggplot(df1_2, aes(x=pi3_1)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p32 = ggplot(df1_2, aes(x=pi3_2)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p33 = ggplot(df1_2, aes(x=pi3_3)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p34 = ggplot(df1_2, aes(x=pi3_4)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p35 = ggplot(df1_2, aes(x=pi3_5)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p41 = ggplot(df1_2, aes(x=pi4_1)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p42 = ggplot(df1_2, aes(x=pi4_2)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p43 = ggplot(df1_2, aes(x=pi4_3)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p44 = ggplot(df1_2, aes(x=pi4_4)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p45 = ggplot(df1_2, aes(x=pi4_5)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p51 = ggplot(df1_2, aes(x=pi5_1)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p52 = ggplot(df1_2, aes(x=pi5_2)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p53 = ggplot(df1_2, aes(x=pi5_3)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p54 = ggplot(df1_2, aes(x=pi5_4)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  p55 = ggplot(df1_2, aes(x=pi5_5)) + geom_density(color="darkblue", fill="lightpink") + xlim(0,xl)
  
  multiplot(p11, p21, p31, p41, p51,
            p12, p22, p32, p42, p52,
            p13, p23, p33, p43, p53,
            p14, p24, p34, p44, p54,
            p15, p25, p35, p45, p55,
            cols=5)
}
trans_plot(model$pi1_2)
trans_plot(model$pi2_1)

options(scipen=999)
plot(model$beta2[[2]][0:1000,4],ylab='Coefficient',xlab='MCMC Iteration',main='Cubic Coefficients MCMC Plot')
plot(model$beta2[[2]][0:1000,4],type='l',ylab='Coefficient',xlab='MCMC Iteration',main='Cubic Coefficients Trace Plot')

par(mfrow = c(3,1))

#trace plots to show convergence
age = cbind(model$beta2[[1]][,2],model$beta1[[2]][,2],model$beta1[[3]][,2])

matplot(c(1:iter), age, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)

agesq = cbind(model$beta1[[1]][,3],model$beta1[[2]][,3],model$beta1[[3]][,3])

matplot(c(1:iter), agesq, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)

age3 = cbind(model$beta1[[1]][,4],model$beta1[[2]][,4],model$beta1[[3]][,4])

matplot(c(1:iter), age3, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)

par(mfrow = c(3,1))

age = cbind(model$beta2[[1]][,2],model$beta1[[2]][,2],model$beta1[[3]][,2])

matplot(c(1:iter), age, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)

agesq = cbind(model$beta2[[1]][,3],model$beta1[[2]][,3],model$beta1[[3]][,3])

matplot(c(1:iter), agesq, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)

age3 = cbind(model$beta2[[1]][,4],model$beta1[[2]][,4],model$beta1[[3]][,4])

matplot(c(1:iter), age3, type='l', xlab='Iteration', ylab='Age Coefficient', col=1:3)
legend('bottomright', inset=.05, legend=c(1,2,3), 
       pch=1, horiz=TRUE, col=1:3)