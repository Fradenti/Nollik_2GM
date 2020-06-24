# processing the simulations
library(reshape2)
library(knitr)
library(tidyverse)
library(mltools)
library(pROC)
library(MLmetrics)
library(modEvA)
library(latex2exp)

##################################################### BAYESIAN FDR -- > pick the threshold
vannucci <- function(k,p,alpha){
  fdr <- sum((1-p)*(p>k))/sum(p>k) 
  fdr  - alpha
}
choose_threshold <- function(prob,AL){
  uniroot(function(x) vannucci(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
}
# beware: my function picks as threshold the last NON-valid obs
# Michele's, Simon's ecc's one picks the first valid one

acc <- function(preds,actual) mean(preds==actual)
##################################################### All the performance indexes from BNPT testing
compute_stat2 <- function(A,P,mppi){ #actual, #pred
  PRECISION   = Precision(y_true = A,y_pred = P,positive = 1)
  RECALL      = Recall(y_true = A,y_pred = P,positive = 1)
  SPECIFICITY = Specificity(y_true = A, y_pred = P,positive = 1)
  s=c(MCC = mcc(preds = P,actuals = A), 
      F1  = 2/(1/RECALL+1/PRECISION),
      SPECIFICITY=SPECIFICITY,
      ACC=acc(P,A),
      PRE = PRECISION,
      AUC = modEvA::AUC(obs = A, pred =mppi, simplif = T)
  )
  return(s)
}
#####################################################
Extract_summaries_MCMC <- function(a){
  aa <- matrix(NA,length(a$Kcon),9)
  aa[,1] <- a$RHOcon
  aa[,2] <- a$M1con
  aa[,3] <- a$S1con
  aa[,4] <- a$M2con
  aa[,5] <- a$S2con
  aa[,6] <- a$THcon
  aa[,7] <- a$S0
  aa[,8] <- a$M0
  aa[,9] <- a$a
  return(as.mcmc(aa))  
}
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
f_post_i <- function(x,res1, nsim, i){
  (1-(res1$RHOcon[nsim-i]))*dnorm(x,(res1$M0[nsim-i]),sqrt(res1$S0[nsim-i])) +  
    (res1$RHOcon[nsim-i]) * wd_mixmom_theta(z = x,theta = (res1$THcon[nsim-i]),
                                            m1 =    (res1$M1con[nsim-i]),
                                            m2 =    (res1$M2con[nsim-i]),
                                            s1 =    (res1$S1con[nsim-i]),
                                            s2 =    (res1$S2con[nsim-i]),
                                            k  =    (res1$Kcon[nsim-i]),
                                            a  =    (res1$a[nsim-i]))
}
pif0_post_i <- function(x,res1, nsim,i){
  (1-(res1$RHOcon[nsim-i]))*dnorm(x,res1$M0[nsim-i],sqrt((res1$S0[nsim-i])))
}
#########################################################################################
locdfr_post_i   <- function(x,res1, nsim,i) pif0_post_i(x,res1, nsim,i)/f_post_i(x,res1, nsim, i)
locdfr_post_k_i <- function(x,res1,kkk, nsim,i) pif0_post_i(x,res1, nsim,i)/f_post_i(x, res1, nsim,i)-kkk
#####################################################################################################

##################################################### Function for Nice Summary -- 1 single scenario, 1 dataset
SUMMA2 <- function(MA,TZ,TH){
  KKK <- matrix(NA,6,ncol(MA))
  ROC <- array(NA,c(1001,2,ncol(MA)))
  for(i in 1:ncol(MA)){
    #mppi      <- colMeans(t(MA[,i]$Lcon))
    mppi      <- 1-locdfr_post(MA[,i]$data,MA[,i]) 
    thr       <- choose_threshold(mppi,TH)
    cat(thr)
    P         <- ifelse(mppi>thr,1,0)
    KKK[,i]   <- compute_stat2(A = TZ,P = P,mppi = mppi)
    ss  <- modEvA::AUC(obs = TZ, pred =mppi,interval = .001,plot = F)$th
    ROC[,,i] <- cbind(ss$sensitivity,1-ss$specificity) ## Errore here
  }
  return(list(KKK,ROC)) 
}

############################################################################################
# ERRATO
# post_process <- function(res1,ZZ,thre, truez=NULL,tit="",mppi){
#   th <- cbind(uniroot(function(x) locdfr_post_k(x,res1,thre), c(-10,0)  )$root,
#               uniroot(function(x) locdfr_post_k(x,res1,thre), c( 0,10)  )$root)
#   curve(locdfr_post(x,res1),-10,10,main=tit)
#   points( 1-mppi ~ ZZ,cex=.5,pch="+")
#   if(is.null(truez)){
#     truez <- c(rep(0,length(ZZ)*.9),rep(1,length(ZZ)*.1))}
#   postz <- ifelse( ZZ<th[1] | ZZ>th[2], 1, 0 )
#   plot(ZZ,col=postz+1,pch=truez+1)
#   abline(h=th)
#   return(list(th,PRED=postz,ACTUAL=truez,table(truez,postz),mppi=apply(res1$Lcon,1,mean)))
# }
















RESULTS2 <- function(OUT,data=NULL, tz,BFDR){
  L <- dim(OUT)[2]
  LIST <- list()
  thresholds <- matrix(NA,L,2)
  results    <- matrix(NA,L,6)
  fdr        <- array(NA,dim=c(101,2,L))
  tz         <- c(rep(0,900),rep(1,100))
  for(i in 1:dim(OUT)[2]){
    data           <- OUT[,i]$data
    #mppi          <- colMeans(t(OUT[,i]$Lcon))
    mppi           <- 1 - locdfr_post(OUT[,i]$data,OUT[,i]) 
    thr            <- choose_threshold(mppi,BFDR)
    #A              <- post_process(OUT[,i],data,thr,tit = i,mppi = mppi) # threshold trovata con BFDR
    #LIST[[i]]      <- A
    #thresholds[i,] <- A[[1]]
    predicted      <- as.numeric(mppi>thr)
    LIST[[i]]      <- table(predicted,tz)
    
    results[i,]    <- as.numeric(compute_stat2(A = tz, P = predicted, mppi = mppi))
    ff             <- curve(locdfr_post(x,OUT[,i]),-10,10,main=i)
#    abline(v=A[[1]])
    fdr[,,i]       <- cbind(ff$x,ff$y) 
    cat(i)
  }
  
  ### ggplot
  plotfdr       <- matrix(NA,101,L+1)
  plotfdr[,1:2] <- fdr[,,1]
  for(i in 1:L){
    plotfdr[,i+1] <-  1- fdr[,2,i]
  }
  ggplotfdr <- melt(as.tibble(plotfdr),id=1)
  p1        <- ggplot(data=ggplotfdr)+
    geom_line(aes(x=V1,y=value,group=variable),alpha=.2,col="blue")+theme_bw()+xlab(TeX("$z$-score"))+ylab(TeX("$\\hat{P}_{1}(z)$"))
  
  #####################################################
  LIST[[L+1]] <- p1
  LIST[[L+2]] <- results
  #####################################################
  return(LIST)
}












post_processS3S4 <- function(res1,ZZ,thre, TZ=NULL,mppi){
  th <- cbind(uniroot(function(x) locdfr_post_k(x,res1,thre), c(-10,0)  )$root,
              uniroot(function(x) locdfr_post_k(x,res1,thre), c( 0,10)  )$root)
  
  curve(locdfr_post(x,res1),-10,10)
  points( 1-mppi ~ ZZ,cex=.5,pch="+")
  
  truez <- TZ
  postz <- ifelse( ZZ<th[1] | ZZ>th[2], 1, 0 )
  plot(ZZ,col=postz+1,pch=truez+1)
  abline(h=th)
  return(list(th,PRED=postz,ACTUAL=truez,table(truez,postz),mppi=apply(res1$Lcon,1,mean)))
}


RESULTS2_S3S4 <- function(OUTS3,data=NULL, tz=NULL,BFDR){
  L <- dim(OUTS3)[2]
  LIST <- list()
  thresholds <- matrix(NA,L,2)
  results    <- matrix(NA,L,6)
  fdr        <- array(NA,dim=c(101,2,L))
  OUT <- list()
  for(k in 1:L){
    OUT[[k]] <- OUTS3[,k][[1]]
  }
  
  for(i in 1:length(OUT)){
    data           <- OUT[[i]]$data
    TZ             <- OUTS3[,i][[2]]
    mppi2          <- colMeans(t(OUT[[i]]$Lcon))
    mppi           <- 1-locdfr_post(OUT[[i]]$data,OUT[[i]]) 
    thr            <- choose_threshold(mppi,BFDR)
    #A              <- post_processS3S4(res1 = OUT[[i]],ZZ = data,thre = thr,TZ = TZ,mppi = mppi)
    LIST[[i]]      <- thr
    predicted      <- as.numeric(mppi>thr)
    results[i,]    <- as.numeric(compute_stat2(A =TZ,P =predicted,mppi = mppi))
    ff             <- curve(locdfr_post(x,OUT[[i]]), -10, 10,main=i)
    #abline(v=A[[1]])
    fdr[,,i]       <- cbind(ff$x,ff$y) 
    cat(i)
  }
  
  plotfdr <- matrix(NA,101,L+1)
  plotfdr[,1:2] <- fdr[,,1]
  for(i in 1:L){
    plotfdr[,i+1] <-   1-fdr[,2,i]
  }
  ggplotfdr <- melt(as.tibble(plotfdr),id=1)
  p1 <- ggplot(data=ggplotfdr)+
    geom_line(aes(x=V1,y=value,group=variable),alpha=.2,col="blue")+theme_bw()+xlab(TeX("$z$-score"))+ylab(TeX("$\\hat{P}_{1}(z)$"))
  ##################################################################  
  LIST[[L+1]] <- p1
  LIST[[L+2]] <- results
  ##################################################################  
  return(LIST)
}





# 
# GIMME_MEAN_of_MEANS <- function(a){
# mat <- matrix(NA,50,9)
#     for( i in 1:50){
#     mat[i,] <- colMeans(Extract_summaries_MCMC(a[,i]))
#     }
# return(mat)
# }
# 
# GIMME_MEAN_of_MEANS <- function(a,S3S4=F){
#   mat <- matrix(NA,50,9)
#   if(S3S4==T){
# 
#   for( i in 1:50){
#     dd <- a[,i][[1]]
#     mat[i,] <- colMeans(Extract_summaries_MCMC(dd))  
#   }
#   }else{
#   for( i in 1:50){
#     mat[i,] <- colMeans(Extract_summaries_MCMC(a[,i]))
#   }}
#   return(mat)
# }
# 
# 
# apply(matMean,2,function(x) c(mean(x),sd(x)))
# 
# 
# 
#   # mat1 <- apply(RS1[[52]],2,function(x) cbind(mean(x),sd(x)))
#   # vec1 <- numeric(12)
#   # vec1[c(1,3,5,7,9,11)] <- mat1[1,]
#   # vec1[c(1,3,5,7,9,11)+1] <- mat1[2,]
#   mat1 <- apply(RS4[[52]],2,function(x) cbind(mean(x),sd(x)))
#   tmat1 <- round(as.data.frame(t(mat1)),4)
#   tmat1[,2] <- paste0("(",tmat1[,2],")")
# tmat1
# 
# 
# matMean <- GIMME_MEAN_of_MEANS(S4,S3S4 = T)
# mm <- apply(matMean,2,function(x) c(mean(x),sd(x)))
# tmat1 <- round(as.data.frame(t(mm[,c(1,6,9)])),4)
# tmat1[,2] <- paste0("(",tmat1[,2],")")
# tmat1