zz <- res1$data
zg <- rep(1:3,rep(1899,3))
G <- 3
res2 <- Nollik_w3_hier_v3(z = zz, 
                          zg = zg,  
                          hyperpar = HP, 
                          NSIM = 10000,
                          BURN = 100000,
                          THIN = 30,
                          verbose.step = 50,
                          SIG1 = array(.5*diag(2),c(2,2,G)), 
                          SIG2 = array(.5*diag(2),c(2,2,G)),  
                          sigA = rep(.5,G))
#######################################################################
plot(res2$RHOcon)
colMeans(res1$THcon)
###############################################################################################

z <- zz
res1 <- res2
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


plot(res1$M1con)
plot(res1$M2con)
plot(res1$M0)
plot(res1$a)

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

x <- seq(-7,7,length.out = 1000)
f1 <- f_postj(x,res1,1)
f2 <- f_postj(x,res1,2)
f3 <- f_postj(x,res1,3)
EFF <- list(f1,f2,f3)

plot.dens <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFF[[g]]~x,col=1,pch=".")  
}
plot.dens(1)
plot.dens(2)
plot.dens(3)


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



plot.dens2 <- function(g){
  hist(z[zg==g],breaks = 40,freq=F)
  points(EFFE0[[g]]~x,col=2,pch=".")  
  points(EFFE1[[g]]~x,col=4,pch=".")  
}
plot.dens2(1)
plot.dens2(2)
plot.dens2(3)
