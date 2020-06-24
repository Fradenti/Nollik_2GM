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
    mppi      <- MA[,i]$Lcon
    #mppi      <- 1-locdfr_post(MA[,i]$data,MA[,i]) 
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
















RESULTS2bnp <- function(OUT,data=NULL, tz,BFDR){ # uso mppi invece the 1-locfdr per label swithcing
  L <- dim(OUT)[2]
  LIST <- list()
  thresholds <- matrix(NA,L,2)
  results    <- matrix(NA,L,6)
  fdr        <- array(NA,dim=c(101,2,L))
  tz         <- c(rep(0,900),rep(1,100))
  for(i in 1:dim(OUT)[2]){
    data            <- OUT[,i]$data
    mppi            <- OUT[,i]$Lcon
  #  mppi           <- 1 - locdfr_post(OUT[,i]$data,OUT[,i]) 
    thr             <- choose_threshold(mppi,BFDR)
    #A              <- post_process(OUT[,i],data,thr,tit = i,mppi = mppi) # threshold trovata con BFDR
    #LIST[[i]]      <- A
    #thresholds[i,] <- A[[1]]
    predicted      <- as.numeric(mppi>thr)
    LIST[[i]]      <- table(predicted,tz)
    
    results[i,]    <- as.numeric(compute_stat2(A = tz, P = predicted, mppi = mppi))
    cat(i)
  }
  
  #####################################################
  return(results)
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


RESULTS2_S3S4bnp <- function(OUTS3,data=NULL, tz=NULL,BFDR){
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
    mppi           <- (OUT[[i]]$Lcon)
    #mppi          <- 1-locdfr_post(OUT[[i]]$data,OUT[[i]]) 
    thr            <- choose_threshold(mppi,BFDR)
    #A             <- post_processS3S4(res1 = OUT[[i]],ZZ = data,thre = thr,TZ = TZ,mppi = mppi)
    LIST[[i]]      <- thr
    predicted      <- as.numeric(mppi>thr)
    results[i,]    <- as.numeric(compute_stat2(A =TZ,P =predicted,mppi = mppi))
  }
  
  ##################################################################  
  return(results)
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

Extract_summaries_MCMC_bnp <- function(x){
  aa <- matrix(NA,10000,2)
  aa[,1] <- x$RHOcon
  aa[,2] <- x$a
  return(as.mcmc(aa))  
}

GIMME_MEAN_of_MEANS_bnp <- function(a,S3S4=F){
  mat <- matrix(NA,50,2)
  if(S3S4==T){

  for( i in 1:50){
    dd <- a[,i][[1]]
    mat[i,] <- colMeans(cbind(dd$a,dd$RHOcon))
  }
  }else{
  for( i in 1:50){
    mat[i,] <- colMeans(Extract_summaries_MCMC_bnp(a[,i]))
  }}
  return(mat)
}
# 
# 
  apply(matMean,2,function(x) c(mean(x),sd(x)))
# 
# 
# 
#   # mat1 <- apply(RS1[[52]],2,function(x) cbind(mean(x),sd(x)))
#   # vec1 <- numeric(12)
#   # vec1[c(1,3,5,7,9,11)] <- mat1[1,]
#   # vec1[c(1,3,5,7,9,11)+1] <- mat1[2,]
x=a
  #a     <- RESULTS2bnp(S2,BFDR = .05)
  a     <- RESULTS2bnp(S2,BFDR = .05)
  mat1  <- apply(a,2,function(x) cbind(mean(x),sd(x)))
  tmat1 <- round(as.data.frame(t(mat1)),4)
  tmat1[,2] <- paste0("(",tmat1[,2],")")
  tmat1
# 
# 
# matMean <- GIMME_MEAN_of_MEANS(S4,S3S4 = T)
# mm <- apply(matMean,2,function(x) c(mean(x),sd(x)))
# tmat1 <- round(as.data.frame(t(mm[,c(1,6,9)])),4)
# tmat1[,2] <- paste0("(",tmat1[,2],")")
# tmat1
  
  
  interp <- function(x,S){
    grid   <- S$DENS[,1]
    DENS   <- S$DENS
    minimo <- which.max(x<grid)
    mean(c(DENS[minimo,5],DENS[minimo-1,5]))
  }
  # d <- sapply(a$data,function(x) (interp(x,S)))
  # plot(d~a$Lcon)
  
  RESULTS2_S3S4bnp <- function(OUTS3,data=NULL, tz=NULL,BFDR){
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
      mppi           <- sapply(data,function(x) (interp(x,OUT[[i]])))
      #mppi          <- 1-locdfr_post(OUT[[i]]$data,OUT[[i]]) 
      thr            <- choose_threshold(mppi,BFDR)
      #A             <- post_processS3S4(res1 = OUT[[i]],ZZ = data,thre = thr,TZ = TZ,mppi = mppi)
      LIST[[i]]      <- thr
      predicted      <- as.numeric(mppi>thr)
      results[i,]    <- as.numeric(compute_stat2(A =TZ,P =predicted,mppi = mppi))
    cat(i)
      }
    
    ##################################################################  
    return(results)
  }
  
  x=a
  #a     <- RESULTS2bnp(S2,BFDR = .05)
  a     <- RESULTS2_S3S4bnp(S3,BFDR = .05)
  mat1  <- apply(a,2,function(x) cbind(mean(x),sd(x)))
  tmat1 <- round(as.data.frame(t(mat1)),4)
  tmat1[,2] <- paste0("(",tmat1[,2],")")
  tmat1
  