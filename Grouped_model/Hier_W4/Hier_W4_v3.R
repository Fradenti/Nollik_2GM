# v3.0: x variabile
library(mombf)
library(fOptions)
library(mvtnorm)
library(coda)
library(MCMCpack)
library(tidyverse)
#############################GIBBS FUNCTION#################################################
Nollik_w4_hier_v3              =function(z,zg,hyperpar,
                                         NSIM, BURN, THIN, 
                                         verbose.step=1000,
                                         SIG1 = array(diag(2),c(2,2,G)),
                                         SIG2 = array(diag(2),c(2,2,G)),
                                         sigA = rep(1,G),
                                         batch     = 50,
                                         fixeda    = 0,
                                         optThresh = .44,
                                         cores){
  
  n       = length(z)
  zg_freq = table(zg)
  G       = max(unique(zg))
  ########### Hyperparameters extraction
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
  Kcon     = numeric(NSIM) 
  ACCa     = ACC1     = ACC2 =  matrix(NA,(batch),G)
  Lcon     = matrix(NA,n,NSIM)
  THcon    = RHOcon = numeric(NSIM)
  Xcon     = matrix(NA,n,NSIM)
  A        = THcon = S0 = S1 = S2 = M0 = M1 = M2 = matrix(NA,NSIM,G)
  
  ind    = 0
  k      = KAPPA #sample(min_kN:max_kN,1)
  a      = rinvgamma(G,aa,bb) #rep(3,G)
  rho    = rbeta(1,a_rho,b_rho)
  
  lambda = sample(0:1,n,T,prob = c(1-rho,rho))
  x      = sapply(lambda,function(s) ifelse(s==0,-1,sample(0:1,1,T)))
  
  theta  = rbeta(G,ax,bx)
  s1     = 1/rgamma(G,shape = a01, scale = 1/b01) 
  s2     = 1/rgamma(G,shape = a02, scale = 1/b02)  
  s0     = 1/rgamma(G,shape = a_phi0, scale = 1/b_phi0)
  m1     = rnorm(G, m01, sqrt(V01*s1))
  m2     = rnorm(G, m02, sqrt(V02*s2))
  m0     = rnorm(G, m_phi0, sqrt(V_phi0*s0))
  
  poss_lab_01 = c(0,1)
  
  for(i in 1:(NSIM*THIN + BURN)){
    
    #################################################################################################################################à
    LX     = update_lambda_x_HIER4(z = z, 
                                   zg = zg,
                                   G = G,
                                   rho = rho,
                                   theta = theta,
                                   k = k,n = n,
                                   m1 = m1,m2 = m2,m0 = m0,
                                   s1 = s1,s2 = s2,s0 = s0,
                                   possible_label123 = 1:3,
                                   a = a  )
    lambda = LX[,1]
    x      = LX[,2]
    #################################################################################################################################à
    rho    = update_rho_cpp(lambda = lambda,
                            a_rho = a_rho,
                            b_rho = b_rho,
                            n = n)
    theta  = update_theta_cpp_nest(x =x,z_g = zg,nG = G, ax=ax, bx=bx)
    #################################################################################################################################à
    k      = KAPPA #update_k_cpp(z = z,min_kN=min_kN,max_kN=max_kN, lambda = lambda,x = x,m1 = m1,m2 = m2,s1 = s1,s2 = s2,a = a,possible_label_k = poss_lab_k)
    #################################################################################################################################à
    for(g in 1:G){
      m1s1   = update_m1s1_cpp(z = z[zg==g],
                               x = x[zg==g], 
                               lambda = lambda[zg==g],
                               m0j = m01, V0j = V01,
                               a0j = a01, b0j = b01,
                               prevMj = m1[g], 
                               prevSj = s1[g],
                               SIG = SIG1[,,g], 
                               k = k ,a = a[g])
      m1[g] = m1s1[1]; 
      s1[g] = m1s1[2]; 
      ACC1[i-ind*batch,g] = m1s1[3]
      ###################################################################################################################################
      m2s2   = update_m2s2_cpp(z = z[zg==g],
                               x = x[zg==g], 
                               lambda = lambda[zg==g],
                               m0j = m02, V0j = V02,
                               a0j = a02, b0j = b02,
                               prevMj = m2[g], prevSj = s2[g],
                               SIG = SIG2[,,g], 
                               k = k,a = a[g])
      m2[g] = m2s2[1];  
      s2[g] = m2s2[2]; 
      ACC2[i-ind*batch,g] = m2s2[3]
      #################################################################################################################################à
      z0     = z[lambda==0 & zg==g]
      m0s0   = update_m0s0_cpp(sub_z0 = z0,
                               m_phi0 = m_phi0,
                               V_phi0,
                               a_phi0 = a_phi0,b_phi0 = b_phi0)
      m0[g]  = m0s0[1]; 
      s0[g]  = m0s0[2];
      #################################################################################################################################à
      
      if(fixeda>0){
        a[g] = fixeda;
      }else{
        z1  = z[lambda==1 & zg==g];
        n11 = sum(x==0 & zg==g); 
        n12 = sum(x==1 & zg==g);
        Avec   = MH_step_a_cpp2(previous_a = a[g],
                                m1 = m1[g],
                                s1 = s1[g],
                                m2 = m2[g],
                                s2 = s2[g],
                                aa = aa,
                                bb = bb,
                                sigma_a = sigA[g],
                                sub_z = z1,
                                k = k,
                                n11 = n11, n12 = n12)
       a[g] = Avec[1]; 
       ACCa[i-ind*batch,g] =Avec[2];
      }
      
      if( ((i) %% batch) == 0){
        diag(SIG1[,,g]) = exp(ifelse( mean(ACC1[,g]) < optThresh, (log(diag(SIG1[,,g])) - min(0.01,1/sqrt(i)) ) , (log(diag(SIG1[,,g])) + min(0.01,1/sqrt(i)) ) ))
        diag(SIG2[,,g]) = exp(ifelse( mean(ACC2[,g]) < optThresh, (log(diag(SIG2[,,g])) - min(0.01,1/sqrt(i) )) , (log(diag(SIG2[,,g])) + min(0.01,1/sqrt(i)) ) ))
        sigA[g]         = exp(ifelse( mean(ACCa[,g]) < .44,#optThresh, 
                                    (log(sigA[g])       - min(0.01,1/sqrt(i) )) , (log(sigA[g])       + min(0.01,1/sqrt(i)) ) ))
        if(g==G) {ind <- ind + 1}
      }
      # 
    }
    
    
    if (i > BURN & ((i - BURN) %% THIN == 0)) {
      rr                      <- floor((i - BURN)/THIN)
      RHOcon[rr] = rho
      Lcon[,rr]  = lambda
      Kcon[rr]   = k
      THcon[rr,]  =theta
      M1[rr,]     = m1
      M2[rr,]     = m2
      S1[rr,]     = s1
      S2[rr,]     = s2
      S0[rr,]     = s0
      M0[rr,]     = m0
      Xcon[,rr]   = x
      A[rr,]      = a  
    }
    
    if (i%%(verbose.step*THIN) == 0) {
      cat(paste("Sampling iteration: ", i, " out of ",NSIM*THIN + BURN,"\n",
                sep = "" ))}
  }
  
  return(list(Kcon   = as.mcmc(Kcon),
              RHOcon = as.mcmc(RHOcon),
              Lcon   = as.mcmc(Lcon), 
              M1con  = as.mcmc(M1), 
              M2con  = as.mcmc(M2),
              S1con  = as.mcmc(S1), 
              S2con  = as.mcmc(S2), 
              Xcon   = as.mcmc(Xcon), 
              THcon  = as.mcmc(THcon), 
              S0     = as.mcmc(S0), 
              M0     = as.mcmc(M0), 
              a      = as.mcmc(A),
              data   = z, 
              HP     = hyperpar))
}


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
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

plot.dens <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFF[[g]]~x,col=2,pch=".")  
}


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




plot.dens2 <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFFE0[[g]]~x,col=2,pch=".")  
  points(EFFE1[[g]]~x,col=4,pch=".")  
}
