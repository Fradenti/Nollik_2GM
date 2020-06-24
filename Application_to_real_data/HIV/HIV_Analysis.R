      library(dglars)
      data(alon)
      library(readr)
      library(tidyverse)
      library(reshape2)
      library(viridis)
      library(knitr)
      library(coda)
      require(multtest)
      library(sgof)
      library(golubEsets)
      
      Rcpp::sourceCpp('Param_version/weight_W3/W3.cpp')
      source('Param_version/weight_W3/W3.R')
      
      library(locfdr)
      data(hivdata)
      hist(hivdata,freq = F)
      curve(dnorm,add=T)
      locfdr(hivdata)
      teststat <- hivdata
      qplot(teststat,col=I(4),fill=I("gray"))+theme_bw()
      zstat <- teststat
      qplot(zstat,col=I(4),fill=I("gray"))+theme_bw()
      ##############################################################################
      plt = ggplot(data.frame(zstat), aes(sample = zstat)) + stat_qq() + theme_bw()
      rawp = sort(2 * (1 - pnorm(abs(zstat))))
      plot(rawp,2*pnorm(sort(-abs(zstat))))
      re <- BH(rawp,.05)
      sum(re$Adjusted.pvalues<.1)
      sum(re$Adjusted.pvalues<.15)
      ### FIRST TRIAL
      plot(sort(rawp),pch="O",xlim=c(1,30),ylim=c(0,.001))
      n <- length(rawp)
      lines((1:n)/n*.1)
      lines((1:n)/n*.05)
      
      plot((rawp*n)/(1:n),pch=".")
      points(re$Adjusted.pvalues,pch=".")
      points(p.adjust(rawp,"BH"),pch=".",col=2)
      
      plot((sort(rawp))*(n/(1:(n))),pch="o",xlim=c(0,500))
      points(re$Adjusted.pvalues,pch="o",col=2)
    points(p.adjust(rawp,"BH"),pch="o",col=3)
    plot(re$Adjusted.pvalues,pch="o",col=2,xlim=c(0,100))
    points((sort(rawp))*(n/(1:(n))),pch="o",xlim=c(0,500))
    
    #####################################
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
    
    
    plot.PT <- function(x,TH=.2){
      require("tidyverse")
      mppi   <- (x$Lcon)
      thr    <- choose_threshold(mppi,TH)
      postz  <- ifelse(mppi>thr,1,0)
      D      <- cbind(id=1:length(x$data),y=x$data,mppi=mppi,postz=postz)
      D      <- as.tibble(D)
      
      p1    <- ggplot()+geom_line(data=D,aes(x=y, y=mppi),col="grey")+
        geom_point(data=D,aes(x=y, y=mppi,col=as.factor(postz+1)),alpha=.5,pch="+")+
        theme_bw()+
        geom_hline(yintercept=thr,alpha=1,col=4,linetype=2)+ylab(unname(TeX("$P\\left(\\lambda_i=1|data\\right)$")))+xlab("z")+
        scale_color_manual(guide=F,values=c(2,3))
      p1
    }
    traces <- function(x){
      as.mcmc(cbind(RHO = (x$RHOcon), 
                    TH  = (x$THcon),
                    M1  = (x$M1), 
                    S1  = (x$S1), 
                    M2  = (x$M2),
                    S2  = (x$S2), 
                    M0  = (x$M0),
                    S0  = (x$S0),
                    a   = (x$a)))
    }
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
      m_phi0 = 0,  V_phi0 = 1/100)
    
      set.seed(1234567891)
      
      res1 <- Nollik_w3(z = zstat, 
                        hyperpar = HP, 
                        NSIM = 10000,
                        BURN = 100000,
                        THIN = 50,
                        verbose.step = 1000,
                        SIG1 = diag(c(.5,.5)), 
                        SIG2 = diag(c(.5,.5)),  
                        sigA = .5,
                        optThresh = .4,
                        cores = 3)

#save(res1,file="/home/fra/SIMU_STEF/application_dati_seri/golub_alon/environment_ALON_21AGO_zstat_thinning50.RData")
#load("/home/fra/SIMULAZIONI per articoli/SIMU_STEF/application_dati_seri/golub_alon/W3Alon/environment_ALON_21AGO_zstat_thinning50.RData")
#### Postprocessing Results 1
#saveRDS(object = res1,file = "~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv_th50_Bi100k.RDS")
    
#Proportions of observations from the alternative process / negative component
res1 <- readRDS("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv_th50_Bi100k.RDS")
res1 <- readRDS("Application_to_real_data/HIV/Hiv_th50_Bi100k.RDS")

mean(res1$RHOcon); sd(res1$RHOcon)
#[1] 0.07933
#[1] 0.0111
mean(res1$THcon); sd(res1$THcon)
#[1] 0.1221
#[1] 0.05118
mean(res1$a); sd(res1$a)
#[1] 2.0579
#[1] 0.3097
mean(res1$M0); sd(res1$M0)


plot(traces(x = res1)[,3])
plot(traces(x = res1)[,7:9])


acf(res1$M1con)
acf(res1$M2con)
acf(res1$S1con)
acf(res1$S2con)
####################################################################
pp=posterior_probability_of_alterative=res1$Lcon
# pp24=posterior_probability_of_alterative=apply(res1$Lcon,1,function(x) quantile(x,c(.05,.95)))


z=zstat
d=as_tibble(cbind(ind=1:length(zstat),z=zstat,pp=posterior_probability_of_alterative))
ggplot(data=d)+geom_point(aes(x=ind,y=z,fill=pp),size=2,pch=21)+theme_bw()+scale_fill_gradient(low="red",high = "blue")


# mppi   = pp
# plot(mppi~zstat)
# efron  = locfdr::locfdr(teststat)
# thr    = choose_threshold(prob = mppi,.15)
# curve(locdfr_post(x,res1),min(res1$data),max(res1$data))
# 
# plot(mppi~res1$data)
# curve(1-locdfr_post(x,res1),min(res1$data),max(res1$data),add=T,col=2)
# 
# 
# 
# f0 <- res1$data[mppi<= thr]
# f1 <- res1$data[mppi> thr]
# postz <- mppi> thr
# y0 <- density(f0,bw = 1,n = 600)
# y1 <- density(f1,bw = 1,n = 600)
# A <- postz
# D <- cbind(id=1:length(res1$data),y=res1$data,postz=A,mppi=mppi,thr=thr,mmpiefron = 1-efron$fdr)
# D <- as_tibble(D)
# f <- cbind(y1=y1$y,y0=y0$y,x1=y1$x,x0=y0$x)
# f <- as_tibble(f)
# 
# 
# ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1)+theme_bw()+ 
#   geom_line(data=D,aes(x=y,y=mppi/4),alpha=(.5),col="red") +
#   geom_line(data=D,aes(x=y,y=mmpiefron/4),alpha=(.5),col="blue") +
#   geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
#   scale_y_continuous(sec.axis = sec_axis(~.*4, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
#   geom_hline(yintercept=thr/4,alpha=1,col=2,linetype=2)+
#   ylab("Density")+xlab("Observed z-score") 
# ggsave("/home/fra/SIMU_STEF/application_dati_seri/golub_alon/Alon_output_zstat_correttothin50_soglia15.png",height = 7,width = 7)

efron <- locfdr::locfdr(zstat)
# sum(zstat< -2.471104 | zstat >  2.227792)
# sum(efron$fdr<(1-.84))
#efron <- locfdr::locfdr(teststat)
plot(efron$fdr~res1$data)
mppi2 <- 1-locdfr_post(res1$data,res1)  # RAO BLACKWELLIZED VERSION (DO), 1-stima fdr

plot(mppi2~zstat)
#points(res1$Lcon~zstat,col=2)
thr2    = choose_threshold(prob = mppi2,.05)
#[1] 0.8422263
signif = sum(mppi2 > thr2)

sum(efron$fdr<.2)
#signif = sum(mppi2 > .8)
sum(signif)# 144

f0 <- res1$data[mppi2<= thr2]
f1 <- res1$data[mppi2> thr2]
postz <- mppi2> thr2
y0 <- density(f0,bw = 1,n = 600)
y1 <- density(f1,bw = 1,n = 600)
A <- postz
D <- cbind(id=1:length(res1$data),y=res1$data,postz=A,mppi=mppi2,thr=thr2,mmpiefron = 1-efron$fdr)
D <- as_tibble(D)
f <- cbind(y1=y1$y,y0=y0$y,x1=y1$x,x0=y0$x)
f <- as_tibble(f)

# 
# ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
#   geom_line(data=D,aes(x=y,y=mppi2/4),alpha=(.5),col="red") +
#   geom_line(data=D,aes(x=y,y=mmpiefron/4),alpha=(.5),col="blue") +
#   geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
#   scale_y_continuous(sec.axis = sec_axis(~.*4, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
#   geom_hline(yintercept=thr2/4,alpha=1,col=2,linetype=2)+
#   ylab("Density")+xlab("Observed z-score") 
# ggsave("/home/fra/SIMU_STEF/application_dati_seri/golub_alon/Alon_output_zstat_corretto_RAOBLACKWELLIZED_thin50_soglia10_efron_nohigherDF.png",height = 7,width = 7)


# a <- res1$Lcon
# xxx <- seq(min(res1$data),max(res1$data),length.out = length(res1$data))
# fx  <- f_post(xxx,res1)
# f1x <- pif1_post(xxx,res1)
# f0x <- pif0_post(xxx,res1)
# 
# D <- as_tibble(cbind(xxx,fx,f1x,f0x,y=res1$data, norm=dnorm(xxx)))
# 
# ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1,bins = 50)+theme_bw()+ xlim(-10,10)+
#   geom_line(data=D,aes(x=xxx,y=f0x),col="blue") +
#   geom_line(data=D,aes(x=xxx,y=f1x),col="red") +
#   geom_line(data=D,aes(x=xxx,y=fx) ,col="black")+
#   geom_line(data=D,aes(x=xxx,y=norm) ,col="black",lty=3)+
#   ylab("Density")+xlab("Observed z-score") 
# ggsave("/home/fra/SIMU_STEF/alon_hist_zstat_corretto_thin50.png",height = 7,width = 7)
# 


plot(efron$fdr~res1$data)
mppi2 <- 1-locdfr_post(res1$data,res1)  # RAO BLACKWELLIZED VERSION (DO), 1-stima fdr
thr2    = choose_threshold(prob = mppi2,.05)
signif = sum(mppi2 > thr2)

thr2_1    = choose_threshold(prob = res1$Lcon,.05)
signif2   = sum(res1$Lcon > thr2_1)

sum(efron$fdr<.2)
sum(signif)

f0 <- res1$data[mppi2<= thr2]
f1 <- res1$data[mppi2> thr2]
postz <- mppi2> thr2
y0 <- density(f0,bw = 1,n = 600)
y1 <- density(f1,bw = 1,n = 600)
A <- postz
D <- cbind(id=1:length(res1$data),y=res1$data,postz=A,mppi=mppi2,thr=thr2,mmpiefron = 1-efron$fdr)
D <- as_tibble(D)
f <- cbind(y1=y1$y,y0=y0$y,x1=y1$x,x0=y0$x,nn=dnorm(y0$x))
f <- as_tibble(f)

D <- D %>% mutate(norm=dnorm(y))

# New plot
ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=D,aes(x=y,y=mppi2/2),alpha=(.8)) +
  geom_line(data=D,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
  #geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
  scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(yintercept=thr2/2,alpha=1,linetype=3)+
  geom_line(data=D,aes(x=y,y=norm),linetype=5)+
  ylab("Density")+xlab("Observed z-score") 
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv1_th50_bi100k.png",height = 4.98, width = 7.8)

# 7.8 x 4.98 
# New plot
ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=D,aes(x=y,y=mppi2/2),alpha=(.8)) +
  geom_line(data=D,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
  #geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
  scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(yintercept=thr2/2,alpha=1,linetype=3)+
  ylab("Density")+xlab("Observed z-score") 
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv2_th50_bi100k.png",height = 4.98, width = 7.8)


# New plot
ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=D,aes(x=y,y=mppi2/2),alpha=(.8)) +
#  geom_line(data=D,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
  #geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
  scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(yintercept=thr2/2,alpha=1,linetype=3)+
  geom_line(data=D,aes(x=y,y=norm),linetype=5)+
  ylab("Density")+xlab("Observed z-score") 
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv3.png",height = 4.98, width = 7.8)



ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),col=1, fill="lightgray")+theme_bw()+ 
  geom_line(data=D,aes(x=y,y=mppi2/2),alpha=(.8)) +
  #  geom_line(data=D,aes(x=y,y=mmpiefron/2),alpha=(.8),linetype=4) +
  #geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
  scale_y_continuous(sec.axis = sec_axis(~.*2, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(yintercept=thr2/2,alpha=1,linetype=3)+
  #geom_line(data=D,aes(x=y,y=norm),linetype=5)+
  ylab("Density")+xlab("Observed z-score") 
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv4.png",height = 4.98, width = 7.8)


a <- res1$Lcon
xxx <- seq(min(res1$data),max(res1$data),length.out = length(res1$data))
fx  <- f_post(xxx,res1)
f1x <- pif1_post(xxx,res1)
f0x <- pif0_post(xxx,res1)

D <- as_tibble(cbind(xxx,fx,f1x,f0x,y=res1$data, norm=dnorm(xxx)))

ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),alpha=.1,col=1,bins = 50)+theme_bw()+
  geom_line(data=D,aes(x=xxx,y=norm) ,col="black",lty=3)+
  geom_line(data=D,aes(x=xxx,y=fx),col="blue") +
  geom_line(data=D,aes(x=xxx,y=f1x),col="red") +
  geom_line(data=D,aes(x=xxx,y=f0x) ,col="black")+
  ylab("Density")+xlab("Observed z-score")
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/HIV/Hiv_2_th50_bi100k.png",height = 4.98,width = 7.8)
