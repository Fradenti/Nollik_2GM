
library(reshape2)
library(knitr)
library(tidyverse)
library(mltools)
library(pROC)
library(MLmetrics)
library(modEvA)
library(latex2exp)
#load("C://Users/frali/Desktop/Application to golub data/golub.RData")
vannucci <- function(k,p,alpha){
  fdr <- sum((1-p)*(p>k))/sum(p>k) 
  fdr  - alpha
}
choose_threshold <- function(prob,AL){
  uniroot(function(x) vannucci(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
}
# beware: my function picks as threshold the last NON-valid obs
# Michele's, Simon's ecc's one picks the first valid one
##############################################################################################################
wd_mixmom_theta <- function(z,theta,m1,m2,k,s1=1,s2=1,a=1, logscale=F){
  g <- function(z)   {
    (1-exp(-(z/a)^(2*k)))*
      ( (1-theta) * dnorm(z,m1,sqrt(s1)) + theta * dnorm(z,m2,sqrt(s2)) )
  }
  Kost <- integrate(g,-40,40)$v
  RR <-   log(g(z)) - log(Kost)
  if(logscale){
    return(RR)
  }else{
    return(exp(RR))
  }
}
################################################################################################################
f_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon)) * dnorm(x,mean(res1$M0[,j]),sqrt(mean(res1$S0[,j]))) +  
    mean(res1$RHOcon)  * wd_mixmom_theta(z = x,theta = mean(res1$THcon[,j]),
                                         m1 =    mean(res1$M1con[,j]),
                                         m2 =    mean(res1$M2con[,j]),
                                         s1 =    mean(res1$S1con[,j]),
                                         s2 =    mean(res1$S2con[,j]),
                                         k  =    round(mean(res1$Kcon)),
                                         a  =    mean(res1$a[,j]))
}
################################################################################################################
f1_postj <- function(x,res1,j){
  mean(res1$RHOcon)  * wd_mixmom_theta(z = x,theta = mean(res1$THcon[,j]),
                                       m1 =    mean(res1$M1con[,j]),
                                       m2 =    mean(res1$M2con[,j]),
                                       s1 =    mean(res1$S1con[,j]),
                                       s2 =    mean(res1$S2con[,j]),
                                       k  =    round(mean(res1$Kcon)),
                                       a  =    mean(res1$a[,j]))
}
f0_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon)) * dnorm(x,mean(res1$M0[,j]),sqrt(mean(res1$S0[,j]))) 
}
locfdr_j <- function(x,res1,j){
  1-f1_postj(x,res1,j)/(f0_postj(x,res1,j)+f1_postj(x,res1,j))
}

res1 <- readRDS("~/Git_repo_Ahead/Nollik2GM/Grouped_model/New_limma_pm3_th50_BI100k.RDS")
res1$HP

plot(ts(res1$a))

plot(res1$RHOcon)
colMeans(res1$THcon)
apply(res1$THcon,2,sd)
colMeans(res1$a)
apply(res1$a,2,sd)

mean(res1$RHOcon)
sd(res1$RHOcon)

###############################################################################################

PLOT.MAKER.HIER <- function(z,zg,res1,j){
  
  efron <- locfdr::locfdr(z[zg==j],plot = F)
  
  mppi2 <- 1 - locfdr_j(res1$data[zg==j],res1,j)  # RAO BLACKWELLIZED VERSION (DO), 1-stima fdr
  thr2  <- choose_threshold(prob = mppi2,.05)
  cat("soglia:", thr2, "\n Numero relevant:", sum(mppi2>thr2),
      "\n Relevant Efron:", sum(efron$fdr<.2))
  D <- cbind(id=1:length(res1$data[zg==j]),
             y=res1$data[zg==j],
             mppi=mppi2,
             thr=thr2,
             mmpiefron = 1-efron$fdr)
  D <- as_tibble(D)
  D <- D %>% mutate(norm=dnorm(y))
  D1 <- D %>% mutate(THR=thr2)
  # New plot
  P1 <- ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
    geom_line(data=D,aes(x=y,y=mppi/2),alpha=(.8)) +
    geom_line(data=D,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
    #geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
    scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
    geom_hline(yintercept=thr2/2,alpha=1,linetype=3)+
  #  geom_line(data=D,aes(x=y,y=norm),linetype=5)+
    ylab("Density")+xlab("Observed z-score") 
  
  a   <- res1$Lcon
  xxx <- seq(min(res1$data[zg==j])-2,max(res1$data[zg==j])+2,length.out = length(res1$data[zg==j]))
  f_j  <- f_postj(x = xxx,res1 = res1,j =  j)
  f1_j <- f1_postj(x = xxx,res1 = res1,j = j)
  f0_j <- f0_postj(x = xxx,res1 = res1,j = j)
  
  D <- as_tibble(cbind(xxx,f_j,f1_j,f0_j,
                       y=res1$data[zg==j], norm=dnorm(xxx)))
  D2 <- D
  P2 <- ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1,bins = 30)+theme_bw()+
    geom_line(data=D,aes(x=xxx,y=f0_j),col="blue") +
    geom_line(data=D,aes(x=xxx,y=f1_j),col="red") +
    geom_line(data=D,aes(x=xxx,y=f_j) ,col="black",lty=2)+
    geom_line(data=D,aes(x=xxx,y=norm) ,col="black",lty=3)+
    ylab("Density")+xlab("Observed z-score") 
  
  return(list(D1,P1,D2,P2))
}

z <- res1$data
zg <- rep(1:3,rep(1899,3))
out1 <- PLOT.MAKER.HIER(z,zg,res1 = res1,j = 1)
out2 <- PLOT.MAKER.HIER(z,zg,res1,2)
out3 <- PLOT.MAKER.HIER(z,zg,res1,3)

out1[[2]]
out2[[2]]
out3[[2]]



DD1 <- rbind(out1[[1]],out2[[1]],out3[[1]])
DD1 <- DD1 %>% mutate(group=rep(c("Ubi1","Ubi4","Ubi6"),rep(1899,3)))

KAPPA <- 3

PP1 <- 
  ggplot(data=DD1)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=DD1,aes(x=y,y=mppi/KAPPA),alpha=(.8)) +
  geom_line(data=DD1,aes(x=y,y=mmpiefron/KAPPA),alpha=(.8),linetype=4) +
  scale_y_continuous(sec.axis = sec_axis(~.*KAPPA, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(data=DD1,aes(yintercept=THR/KAPPA),alpha=1,linetype=3)+
#  geom_line(data=DD1,aes(x=y,y=norm),linetype=5)+
  facet_wrap(~group)+
  ylab("Density")+xlab("Observed z-score") 
PP1
ggsave("~/Git_repo_Ahead/Nollik2GM/Grouped_model/Hier_Proteomics_1_th50_bi100k.png",width = 15,height = 5)

  



  DD2 <- rbind(out1[[3]],out2[[3]],out3[[3]])
DD2 <- DD2 %>% mutate(group=rep(c("Ubi1","Ubi4","Ubi6"),rep(1899,3)))
PP2 <- 
  ggplot(data=DD2)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1,bins = 40)+theme_bw()+
  geom_line(data=DD2,aes(x=xxx,y=norm) ,col="black",lty=3)+
  geom_line(data=DD2,aes(x=xxx,y=f_j),col="blue") +
  geom_line(data=DD2,aes(x=xxx,y=f1_j),col="red") +
  geom_line(data=DD2,aes(x=xxx,y=f0_j) ,col="black")+
  ylab("Density")+xlab("Observed z-score")+
  facet_wrap(~group)
ggsave("Grouped_model/Hier_Proteomics_2_1_th50_bi100.png",width = 15,height = 5)

PP2




rawp1 <- (2*pnorm(-abs(z[zg==1])))
rawp4 <- (2*pnorm(-abs(z[zg==2])))
rawp6 <- (2*pnorm(-abs(z[zg==3])))

library(sgof)
BH(rawp1)
BH(rawp4)
BH(rawp6)
