    library(reshape2)
    library(knitr)
    library(tidyverse)
    library(mltools)
    library(sgof)
    library(coda)
    library(pROC)
    library(MLmetrics)
    library(modEvA)
    library(latex2exp)
    library("phyloseq"); packageVersion("phyloseq")
    library(edgeR)
    library(DESeq2)
    library(microbiomeSeq)
    library(edgeR)
    library(latex2exp)
    library(DESeq2)
    Rcpp::sourceCpp("BNP_version/BNP_W4/W4_BNP.cpp")
    source("BNP_version/BNP_W4/W4_BNP.R")
    
    #############################################################################
    vannucci <- function(k,p,alpha){
      fdr <- sum((1-p)*(p>k))/sum(p>k) 
      fdr  - alpha
    }
    choose_threshold <- function(prob,AL){
      uniroot(function(x) vannucci(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
    }
    #############################################################################
    plot.PT <- function(x,TH=.2){
      require("tidyverse")
      mppi   <- rowMeans(x$Lcon)
      thr    <- choose_threshold(mppi,TH)
      postz  <- ifelse(mppi>thr,1,0)
      D      <- cbind(id=1:length(x$data),y=x$data,mppi=mppi,postz=postz)
      D      <- as.tibble(D)
      p1     <- ggplot()+geom_line(data=D,aes(x=y, y=mppi),col="grey")+
        geom_point(data=D,aes(x=y, y=mppi,col=as.factor(postz+1)),alpha=.5)+
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
    #############################################################################
    
    
    
    
    #############################################################################
    data(pitlatrine)
    physeq <- pitlatrine
    table(physeq@sam_data$Country)
    physeq
    library(phyloseq)
    #############################################################################
    pscp = transform_sample_counts(physeq, function(x){x/sum(x)})
    hist(log10(apply(otu_table(pscp), 2, var)),
         xlab="log10(variance)", breaks=50,
         main="A large fraction of OTUs have very low variance")
    varianceThreshold = 1e-7
    keepOTUs = names(which(apply(otu_table(pscp), 2, var) > varianceThreshold))
    physeqB  = prune_taxa(keepOTUs, physeq)
    physeqB
    
    diagdds = phyloseq_to_deseq2(physeqB, ~ Country)
    # calculate geometric means prior to estimate size factors
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
    
    res <- results(diagdds)
    hist(res$stat,freq=F,breaks = 20)
    curve(dnorm,add=T)
    hist(res$log2FoldChange,freq=F)
    curve(dnorm,add=T)
    #######################################################################################
    #######################################################################################
    
    mean(res$stat)
    z1    <- scale(res$stat)
    zstat <- scale(res$stat,center = F)
    
    hist(z1,breaks = 100)
    hist(zstat,breaks = 100,add=T)
    
    hist(zstat,breaks = 100,freq=F)
    curve(dnorm,add=T)
    
    aa <- efron <- locfdr::locfdr(zstat,df = 18)
    plot(aa$fdr~zstat)
    
    sum((aa$fdr)<.2)
    sum((aa$fdr)<.2)/length(aa$fdr)
    
    rawp = sort(2 * (1 - pnorm(abs(zstat))))
    re <- BH(rawp,.05)
    sum(re$Adjusted.pvalues<.05)
    sum(re$Adjusted.pvalues<.15)
    BH(rawp,.05)
    
    
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
    
    set.seed(19922598)
    hist(zstat)
    
    zstat <- readRDS("BNP_version/Pitlatrine_data_scaled_noncentered.RDS")
    HP <- list(
      KAPPA = 2,
      aa    = 20,    bb   = 57,
      a_rho = 1,    b_rho = 9,
      ax    = 1,     bx   = 1, 
      m01   = 0, V01   = 100, 
      a01   = 3,     b01  = 1,
      a_phi0 = 10, b_phi0 = 10, 
      m_phi0 = 0,  V_phi0 = 1/100,
      alpha_a=1,alpha_b=1
    )
    
    res1 <- Nollik_w4_BNP(z = zstat,   
                          hyperpar = HP, 
                          NSIM = 10000,
                          BURN = 100000,
                          J = 30,
                          THIN = 50,
                          verbose.step = 50,
                          sigA = .25, 
                          fixedA = T,
                          setseed = 123456,
                          optThresh = .4 )
saveRDS(object = res1,file = "~/Documents/GitHub/Nollik2GM/Application_to_real_data/Pitlatrine/Pitlatrine_analysis_notcentered_BI100k_thin50_W4.RDS")
# ###############################################################################################0
res1 <- readRDS(file = "~/Documents/GitHub/Nollik2GM/Application_to_real_data/Pitlatrine/Pitlatrine_analysis_notcentered_BI100k_thin50_W4.RDS")

zstat <- z <- res1$data
length(res1$data)
mean(res1$RHOcon);sd(res1$RHOcon)
# [1] 0.09860348
# [1] 0.01698395
mean(res1$a);sd(res1$a)
# [1] 1.998894
# [1] 0.2787727
####################################################################################
###############################################################################################
posterior_probability_of_alterative=res1$Lcon
plot(posterior_probability_of_alterative~zstat)
# d=as.tibble(cbind(ind=1:length(zstat),z=zstat,pp=posterior_probability_of_alterative))
# ggplot(data=d)+geom_point(aes(x=ind,y=z,fill=pp),size=2,pch=21,alpha=.5)+theme_bw()+scale_fill_gradient(low="red",high = "blue")

zstat <- res1$data
plot(zstat,1-posterior_probability_of_alterative,pch="+")
points(zstat,aa$fdr,col=2,pch="x")
###############################################################################################################
###############################################################################################################
mppi <- posterior_probability_of_alterative
plot(mppi~zstat)
thr = choose_threshold(prob = mppi,.05)
thr  #0.8130104
signif = sum(mppi > thr)
#plot.PT(x=res1,.05)+geom_point(aes(x=zstat,y=1-efron$fdr),alpha=.2,col=4)
sum(aa$fdr<.2) # 101
sum(signif)    # 66

plot(res1$a)
plot(res1$RHOcon)
hist(zstat)
plot(res1$M0)
plot(res1$S0)
plot(res1$RHOcon)
plot(t(res1$MuSigArray[,,1]))
plot(t(res1$MuSigArray[,,10000]))
sum(res1$Picon[[1]][10000,])

#View(cbind(res1$RHOcon,res1$S0,res1$M0,res1$a)[7100:7300,])



y <- res1$data
f0 <- res1$data[mppi<= thr]
f1 <- res1$data[mppi> thr]
postz <- mppi> thr
y0 <- density(f0,bw = 1,n = 600)
y1 <- density(f1,bw = 1,n = 600)
A <- postz
D <- cbind(id=1:length(res1$data),y=res1$data,postz=A,mppi=mppi,thr=thr,mmpiefron = 1-efron$fdr)
D <- as_tibble(D)
f <- cbind(y1=y1$y,y0=y0$y,x1=y1$x,x0=y0$x)
f <- as_tibble(f)

KAPPA=1/.7
ggplot(data=D)+geom_histogram(aes(x=y,y=..density..),fill="lightgrey",col=1,bins = 30)+theme_bw()+ 
  geom_line(data=D,aes(x=y,y=mppi/KAPPA),alpha=.8) +
  geom_line(data=D,aes(x=y,y=mmpiefron/KAPPA),alpha=(.8),linetype=4) +
  #  geom_point(data=D,aes(x=y,y=rep(-.02,length(res1$data))),col=as.factor(postz+1),alpha=(.8),pch="|") +
  scale_y_continuous(sec.axis = sec_axis(~.*KAPPA, unname(TeX("$P\\left(z_i=1|data\\right)$"))))  +
  geom_hline(yintercept=thr/KAPPA,alpha=1,linetype=3)+
  ylab("Density")+xlab("Observed z-score") 
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/Pitlatrine/PitLatrine_1_noncentered_BI100_th50_BNPW4.png",width = 7.8,height =  4.98 )

D <- as_tibble(cbind(xxx=res1$DENS[,1],fx=res1$DENS[,4],
                     f1x=res1$DENS[,3],
                     f0x=res1$DENS[,2],
                     norm=dnorm(res1$DENS[,1])))
ggplot()+
  geom_histogram(aes(x=res1$data,y=..density..),alpha=.1,col=1,bins = 15)+theme_bw()+
  geom_line(data=D,aes(x=xxx,y=norm) ,col="black",lty=3)+
  geom_line(data=D,aes(x=xxx,y=fx),col="blue") +
  geom_line(data=D,aes(x=xxx,y=f1x),col="red") +
  geom_line(data=D,aes(x=xxx,y=f0x) ,col="black")+
  ylab("Density")+xlab("Observed z-score")+xlim(-5,5)
ggsave("~/Git_repo_Ahead/Nollik2GM/Application_to_real_data/Pitlatrine/PitLatrine_2_noncentered_BI100_th50_W4.png",width = 7.8,height =  4.98 )



plot(t(res1$MuSigArray[,,10000]))
