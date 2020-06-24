Rcpp::sourceCpp('~/Git_repo_Ahead/Nollik2GM/Grouped_model/Hier_W3_v3.cpp')    
source('~/Git_repo_Ahead/Nollik2GM/Grouped_model/Hier_W3_v3.R')    

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



################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

burn=1000
######################################################
#ESEMPIO SCEMO Hierarchical
###################################################
HP <- list(
  KAPPA = 2,
  aa    = 20,    bb   = 57,
  a_rho = 1,    b_rho = 9,
  ax    = 1,     bx   = 1, 
  m01   = -3,    m02  = 3, 
  V01   = 1,     V02  = 1, 
  a01   = 2,     b01  = 5,
  a02   = 2,     b02  = 5,
  a_phi0 = 10, b_phi0 = 10, 
  m_phi0 = 0,  V_phi0 = 1/100
)  


set.seed(1234)
z1 <- c(rnorm(900,0,.5),rnorm(50,5),rnorm(50,-5))
z2 <- c(rnorm(900,0,.5),rnorm(100,-3))
z3 <- c(rnorm(900,0,.5),rnorm(50,5),rnorm(50,-5))
z4 <- c(rnorm(900,0,0.5),rnorm(50,3),rnorm(50,-3))
z5 <- c(rnorm(900,0,0.5),rnorm(50,2),rnorm(50,-4))
z  <- c(z1,z2,z3,z4,z5)
zg <- rep(1:5,rep(1000,5))
plot(z,col=zg)

G <- length(unique(zg))
res1 <- Nollik_w4_hier_v3(z = z, 
                          zg = zg,  
                          hyperpar = HP, 
                          NSIM = 5000,
                          BURN = 1,
                          THIN = 1,
                          fixeda = 0,
                          verbose.step = 10,
                          SIG1 = array(1*diag(2),c(2,2,G)), 
                          SIG2 = array(1*diag(2),c(2,2,G)),  
                          sigA = rep(.5,G),
                          batch = 50)
# problema con a
###############################################################################################
res1 <- Nollik_w4(z = z1,
                          hyperpar = HP,
                          NSIM = 10000,
                          BURN = 10000,
                          THIN = 1,cores = 1,
                          verbose.step = 500,
                          SIG1 = .5*diag(2),
                          SIG2 = .5*diag(2),
                          sigA = rep(.5,1))

###############################################################################################
plot(res1$RHOcon)
plot(res1$a)
plot(res1$RHOcon)
mean(res1$THcon)
###############################################################################################


pairs(cbind(res1$Lcon[,5000],res1$data),pch="+")
plot(res1$M1con[,5])
plot(res1$a[,1])



plot.prob <- function(g){
hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
ef  <- locfdr::locfdr(z[zg==g],plot = F)
pef <- 1-ef$fdr
p1  <- apply(res1$Lcon[zg==g,],1,mean)
points(p1~z[zg==g],col=2)  
points(pef~z[zg==g],col=4)
}

plot.prob(1)
plot.prob(2)
plot.prob(3)
plot.prob(4)
plot.prob(5)



plot(res1$M1con)
plot(res1$M2con)
plot(res1$M0)
plot(res1$S0)
plot(res1$a)


x <- seq(-7,7,length.out = 1000)
f1 <- f_postj(x,res1,1)
f2 <- f_postj(x,res1,2)
f3 <- f_postj(x,res1,3)
f4 <- f_postj(x,res1,4)
f5 <- f_postj(x,res1,5)
EFF <- list(f1,f2,f3,f4,f5)

plot.dens <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFF[[g]]~x,col=2,pch=".")  
}
plot.dens(1)
plot.dens(2)
plot.dens(3)
plot.dens(4)
plot.dens(5)



EFFE1 <- list()
EFFE0 <- list()

for(g in 1:G){
  
  EFFE0[[g]] <- f0_postj(x,res1,g)
  EFFE1[[g]] <- f1_postj(x,res1,g)
}



plot.dens2 <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFFE0[[g]]~x,col=2,pch=".")  
  points(EFFE1[[g]]~x,col=4,pch=".")  
}
plot.dens2(1)
plot.dens2(2)
plot.dens2(3)
plot.dens2(4)
plot.dens2(5)





# ggsave("/home/fra/SIMU_STEF/alon_hist_zstat_corretto_thin50.png",height = 7,width = 7)
# 

PLOT.MAKER.HIER <- function(z,zg,res1,j){

efron <- locfdr::locfdr(z[zg==j],plot = F)
mppi2 <- 1 - locfdr_j(res1$data[zg==j],res1,1)  # RAO BLACKWELLIZED VERSION (DO), 1-stima fdr
thr2  <- choose_threshold(prob = mppi2,.05)

f0 <- res1$data[mppi2<= thr2]
f1 <- res1$data[mppi2> thr2]
postz <- mppi2> thr2
A <- postz
D <- cbind(id=1:length(res1$data[zg==j]),
           y=res1$data[zg==j],
           postz=A,
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
  geom_line(data=D,aes(x=y,y=norm),linetype=5)+
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


out1 <- PLOT.MAKER.HIER(z,zg,res1,1)
out2 <- PLOT.MAKER.HIER(z,zg,res1,2)
out3 <- PLOT.MAKER.HIER(z,zg,res1,3)
out4 <- PLOT.MAKER.HIER(z,zg,res1,4)
out5 <- PLOT.MAKER.HIER(z,zg,res1,5)


DD1 <- rbind(out1[[1]],out2[[1]],out3[[1]],out4[[1]],out5[[1]])
DD1 <- DD1 %>% mutate(group=zg)
PP1 <- 
  ggplot(data=DD1)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=DD1,aes(x=y,y=mppi/2),alpha=(.8)) +
  geom_line(data=DD1,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
  scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(data=DD1,aes(yintercept=THR/2),alpha=1,linetype=3)+
  geom_line(data=DD1,aes(x=y,y=norm),linetype=5)+
  facet_wrap(~group)+
  ylab("Density")+xlab("Observed z-score") 
PP1
  

out2[[4]]




DD2 <- rbind(out1[[3]],out2[[3]],out3[[3]],out4[[3]],out5[[3]])
DD2 <- DD2 %>% mutate(group=zg)
PP2 <- 
  ggplot(data=DD2)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1,bins = 50)+theme_bw()+
  geom_line(data=DD2,aes(x=xxx,y=f_j),col="blue") +
  geom_line(data=DD2,aes(x=xxx,y=f1_j),col="red") +
  geom_line(data=DD2,aes(x=xxx,y=f_j) ,col="black",lty=2)+
  geom_line(data=DD2,aes(x=xxx,y=norm) ,col="black",lty=3)+
  ylab("Density")+xlab("Observed z-score")+
  facet_wrap(~group)

PP2
