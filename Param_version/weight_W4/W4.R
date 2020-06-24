library(mombf)
library(fOptions)
#library(CharFun)
library(mvtnorm)
library(coda)
library(MCMCpack)
library(tidyverse)
#############################GIBBS FUNCTION#################################################
Nollik_w4                 =function(z,hyperpar,
                                      NSIM, BURN, THIN, 
                                      verbose.step=1000,
                                      SIG1=diag(2),SIG2=diag(2),sigA=1,
                                      batch=50,
                                      optThresh=.44,
                                      cores,setseed=NA){
  
  n = length(z)
  ########### Hyperparameters extraction
  
  # min_kN <- hyperpar$min_kN
  # max_kN <- hyperpar$max_kN
  KAPPA  <- hyperpar$KAPPA
  a_rho  <- hyperpar$a_rho
  b_rho  <- hyperpar$b_rho 
  ax     <- hyperpar$ax
  bx     <- hyperpar$bx
  m01    <- hyperpar$m01 # media a priori hyperprior su m1
  m02    <- hyperpar$m02 # media a priori hyperprior su m2
  V01    <- hyperpar$V01  # varianza a priori hyperprior su m1
  V02    <- hyperpar$V01 # varianza a priori hyperprior su m2
  a01    <- hyperpar$a01 # shape IG hyperprior su s1
  b01    <- hyperpar$b01 # scale IG hyperprior su s1
  a02    <- hyperpar$a02 # shape IG hyperprior su s2
  b02    <- hyperpar$b02 # scale IG hyperprior su s2
  a_phi0 <- hyperpar$a_phi0 
  b_phi0 <- hyperpar$b_phi0 
  m_phi0 <- hyperpar$m_phi0 
  V_phi0 <- hyperpar$V_phi0 # location and shape M0S0
  aa     <- hyperpar$aa
  bb     <- hyperpar$bb #parameters of IG on a

  ###### container
  Tcon = Kcon = RHOcon = numeric(NSIM) 
  ACCa = ACC1 = ACC2 =  numeric(batch)
  Lcon = matrix(NA,n,NSIM)
  THETAcon = numeric(NSIM)
  Xcon = matrix(NA,n,NSIM)
  A = S0 = S1 = S2 = M0 = M1 = M2 = numeric(NSIM)
  
  ################ INITIALIZATION
  if(is.na(setseed)){setseed <- sample(1:1000,1)}
  set.seed(setseed)
  ind    = 0
  k      = KAPPA #sample(min_kN:max_kN,1)
  a      = rinvgamma(1,aa,bb)
  rho    = rbeta(1,a_rho,b_rho)
  
  lambda = sample(0:1,n,T)
  x      = sapply(lambda,function(s) ifelse(s==0,-1,sample(0:1,1,T)))
  
  theta  = rbeta(1,ax,bx)
  s1     = 1/rgamma(1,shape = a01, scale = 1/b01) 
  s2     = 1/rgamma(1,shape = a02, scale = 1/b02)  
  m1     = -abs(rnorm(1, m01, sqrt(V01*s1)))
  m2     =  abs(rnorm(1, m02, sqrt(V02*s2)))
  s0     = 1/rgamma(1,shape = a_phi0, scale = 1/b_phi0)
  m0     = rnorm(1, m_phi0, sqrt(V_phi0*s0))
  
  
  poss_lab_01 = c(0,1)
  # poss_lab_k  = c(min_kN:max_kN)
  t1 <- Sys.time()
  for(i in 1:(NSIM*THIN + BURN)){
    
    #################################################################################################################################à
    LX     = update_lambda_x_cpp_parallel(z = z,rho = rho,theta = theta,k = k,n = n,
                                          m1 = m1,m2 = m2,m0 = m0,
                                          s1 = s1,s2 = s2,s0 = s0,
                                          possible_label123 = 1:3,a = a,
                                          cores=cores  )
    lambda = LX[,1]
    x      = LX[,2]
    #################################################################################################################################à
    rho    = update_rho_cpp(lambda = lambda,
                            a_rho = a_rho,
                            b_rho = b_rho,
                            n = n)
    theta  = update_theta_cpp(x =x, ax=ax, bx=bx)
    #################################################################################################################################à
    k      = KAPPA #update_k_cpp(z = z,min_kN=min_kN,max_kN=max_kN, lambda = lambda,x = x,m1 = m1,m2 = m2,s1 = s1,s2 = s2,a = a,possible_label_k = poss_lab_k)
    #################################################################################################################################à
    m1s1   = update_m1s1_cpp(z = z,
                             x = x, 
                             lambda = lambda,
                             m0j = m01, V0j = V01,
                             a0j = a01, b0j = b01,
                             prevMj = m1, prevSj = s1,
                             SIG = SIG1, k = k ,a = a)
    m1 = m1s1[1]; s1 = m1s1[2]; 
    ACC1[i-ind*batch] = m1s1[3]
    ###################################################################################################################################
    m2s2   = update_m2s2_cpp(z = z,
                             x = x, lambda = lambda,
                             m0j = m02, V0j = V02,
                             a0j = a02, b0j = b02,
                             prevMj = m2, prevSj = s2,
                             SIG = SIG2, k = k,a = a)
    m2 = m2s2[1];  s2 = m2s2[2]; 
    ACC2[i-ind*batch] = m2s2[3]
    #################################################################################################################################à
    z0     = z[lambda==0]
    m0s0   = update_m0s0_cpp(z0,m_phi0,V_phi0,a_phi0,b_phi0)
    m0     = m0s0[1]; s0 = m0s0[2];
    #################################################################################################################################à
    z1     = z[lambda==1];
    n11    = sum(x==0); n12 = sum(x==1);
    Avec   = MH_step_a_cpp2(previous_a = a,m1 = m1,s1 = s1,m2 = m2,s2 = s2,aa = aa,bb = bb,sigma_a = sigA,sub_z = z1,k = k,n11 = n11, n12 = n12)
    a      = Avec[1]; 
    ACCa[i-ind*batch] =Avec[2];
    
    if( ((i) %% batch) == 0){
      diag(SIG1) = exp(ifelse( mean(ACC1) < optThresh, (log(diag(SIG1)) - min(0.01,1/sqrt(i)) ) , (log(diag(SIG1)) + min(0.01,1/sqrt(i)) ) ))
      diag(SIG2) = exp(ifelse( mean(ACC2) < optThresh, (log(diag(SIG2)) - min(0.01,1/sqrt(i) )) , (log(diag(SIG2)) + min(0.01,1/sqrt(i)) ) ))
      sigA       = exp(ifelse( mean(ACCa) < .44,#optThresh, 
                               (log(sigA)       - min(0.01,1/sqrt(i) )) , (log(sigA)       + min(0.01,1/sqrt(i)) ) ))
      ind <- ind + 1
    }
    
    if (i > BURN & ((i - BURN) %% THIN == 0)) {
      rr                      <- floor((i - BURN)/THIN)
      RHOcon[rr] = rho
      Lcon[,rr]  = lambda
      Kcon[rr]   = k
      THETAcon[rr] =theta
      M1[rr]     = m1
      M2[rr]     = m2
      S1[rr]     = s1
      S2[rr]     = s2
      S0[rr]     = s0
      M0[rr]     = m0
      Xcon[,rr]  = x
      A[rr] =a  
    }
    
    if (i%%(verbose.step*THIN) == 0) {
      cat(paste("Sampling iteration: ", i, " out of ",NSIM*THIN + BURN,"\n",
                sep = "" ))}
  }
  t2 <- Sys.time()
  return(list(Kcon   = as.mcmc(Kcon),
              RHOcon = as.mcmc(RHOcon),
              Lcon   = apply(Lcon,1,mean), 
              M1con  = as.mcmc(M1), 
              M2con  = as.mcmc(M2),
              S1con  = as.mcmc(S1), 
              S2con  = as.mcmc(S2), 
              Xcon   = apply(Xcon,1,mean), 
              THcon  = as.mcmc(THETAcon), 
              S0     = as.mcmc(S0), 
              M0     = as.mcmc(M0), 
              a      = as.mcmc(A),
              data   = z, 
              HP     = hyperpar,
              time   = difftime(t2, t1, units='mins')))
}




###############################################################################################################
wd_mixmom_theta <- function(z,theta,m1,m2,k,s1=1,s2=1,a=1, logscale=F){
  g <- function(z)   {
    (exp(-(a/z)^(2*k)))*
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

SIMU2 <- function(z,nsim=10000,burn=10000,thin=1,COR=3,HP,SS=1234){
  res1 <- Nollik_w4(z = z,   hyperpar = HP, NSIM = nsim,
                                    BURN = burn,
                                    THIN = thin,
                                    verbose.step = 5,
                                    SIG1 = diag(c(.5,.5)), 
                                    SIG2 = diag(c(.5,.5)),  
                                    sigA = .5,setseed = SS,
                                    optThresh = .4,
                                    cores = COR )
  return(res1)
}

f_post <- function(x,res1){
  (1-mean(res1$RHOcon))*dnorm(x,mean(res1$M0),sqrt(mean(res1$S0))) +  
    mean(res1$RHOcon) * wd_mixmom_theta(z = x,theta = mean(res1$THcon),
                                        m1 =    mean(res1$M1con),
                                        m2 =    mean(res1$M2con),
                                        s1 =    mean(res1$S1con),
                                        s2 =    mean(res1$S2con),
                                        k  =    mean(res1$Kcon),
                                        a  =    mean(res1$a))
}

pif0_post <- function(x,res1){
  (1-mean(res1$RHOcon))*dnorm(x,mean(res1$M0),sqrt(mean(res1$S0)))
}

pif1_post <- function(x,res1){
  mean(res1$RHOcon) * wd_mixmom_theta(z = x,theta = mean(res1$THcon),
                                      m1 =    mean(res1$M1con),
                                      m2 =    mean(res1$M2con),
                                      s1 =    mean(res1$S1con),
                                      s2 =    mean(res1$S2con),
                                      k  =    mean(res1$Kcon),
                                      a  =    mean(res1$a))
}



locdfr_post <- function(x,res1) pif0_post(x,res1)/f_post(x,res1)
locdfr_post_k <- function(x,res1,kkk) pif0_post(x,res1)/f_post(x,res1)-kkk





######## media di funzioni e non funzione di media
f_post_mcmc <- function(x,res1,NSIM){
  NN <- matrix(NA,length(x),NSIM)
  for(bb in 1:NSIM){
    NN[,bb] <- (1-res1$RHOcon[bb])*dnorm(x,(res1$M0[bb]),sqrt((res1$S0[bb]))) +  
      (res1$RHOcon[bb]) * wd_mixmom_theta(z = x,theta = (res1$THcon[bb]),
                                          m2 =    (res1$M2con[bb]),
                                          m1 =    (res1$M1con[bb]),
                                          s1 =    (res1$S1con[bb]),
                                          s2 =    (res1$S2con[bb]),
                                          k  =    (res1$Kcon[bb]),
                                          a  =    (res1$a[bb]))
  }
  LL <- list(mat=NN, mean=rowMeans(NN))
  return(LL)
}

pif0_post_mcmc <- function(x,res1, NSIM){
  NN <- matrix(NA,length(x),NSIM)
  for(bb in 1:NSIM){
    NN[,bb] <- (1-(res1$RHOcon[bb]))*dnorm(x,(res1$M0[bb]),sqrt((res1$S0[bb])))
  }
  LL <- list(mat=NN, mean=rowMeans(NN))
  return(LL)
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
