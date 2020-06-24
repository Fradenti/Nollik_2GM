library(mombf) 
library(fOptions)
#library(CharFun)
library(mvtnorm)
library(coda)
library(MCMCpack)
library(tidyverse)
#############################GIBBS FUNCTION#################################################
Nollik_w4_BNP                      = function(z,hyperpar,J,
                                              NSIM, BURN, THIN, 
                                              verbose.step=1000,
                                              SIG1=array(diag(c(.5,.5)),dim = c(2,2,J)),
                                              sigA=.5,
                                              batch=50,
                                              optThresh=.44,
                                              fixedA=T,
                                              setseed=NA){
  
  n           = length(z)
  poss_lab_01 = c(0,1)
  poss_lab_0J = 0:J
  ########### Hyperparameters extraction
  ygrid  <- seq(min(z)-5,max(z)+5,length.out = 1000)
  DENS   <- matrix(0,1000,3); DENS[,1] <- ygrid
  ##################################################################
  KAPPA  <- hyperpar$KAPPA
  ##################################################################
  a_rho  <- hyperpar$a_rho
  b_rho  <- hyperpar$b_rho 
  ##################################################################
  m01    <- hyperpar$m01 # media a priori hyperprior su m1
  V01    <- hyperpar$V01
  a01    <- hyperpar$a01 # shape IG hyperprior su s1
  b01    <- hyperpar$b01 # scale IG hyperprior su s1
  ##################################################################
  a_phi0 <- hyperpar$a_phi0 
  b_phi0 <- hyperpar$b_phi0 
  m_phi0 <- hyperpar$m_phi0 
  V_phi0 <- hyperpar$V_phi0 # location and shape M0S0
  ##################################################################
  aa     <- hyperpar$aa
  bb     <- hyperpar$bb #parameters of IG on a
  ##################################################################
  a_alpha<- hyperpar$alpha_a 
  b_alpha<- hyperpar$alpha_b 
  ##################################################################
  ##################################################################
  
  ###### containers
  ACCa       = numeric(batch)
  ACC1       = matrix(NA,J,(batch))
  Xcon       = Lcon     = matrix(NA,n,NSIM)
  Picon      = matrix(NA,J,NSIM)
  MuSig      = matrix(NA,2,J)
  Alphacon   = Kcon = RHOcon =  A = S0 = M0 = numeric(NSIM)
  MuSigArray = array(NA,dim = c(2,J,NSIM))
  ##################################################################
  ##################################################################
  
  # Initialization
  if(is.na(setseed)){ setseed= sample(1:1000,1)}
  set.seed(setseed)
  ind    = 0
  k      = KAPPA #sample(min_kN:max_kN,1)
  ##################################################################
  a      = rinvgamma(1,aa,bb)
  ##################################################################
  rho    = rbeta(1,a_rho,b_rho)
  ##################################################################
  if(fixedA) alpha=1 else{  alpha  = rgamma(1,a_alpha,b_alpha) }
  
  uis    = rbeta(J,1,alpha); 
  uis[J] <- 1
  pis    = SB_given_u2(uis)  

  lambda = sample(0:1,n,T,prob = c(1-rho,rho))
  x      = sapply(lambda,function(s) ifelse(s==0,0,sample(1:J,1,T,prob = pis)))
  ##################################################################
  MuSig[2,] <- rinvgamma(J,shape = a01, scale = b01) 
  MuSig[1,] <- rnorm(J, m01, sqrt(V01*MuSig[2,]))
  s0        <- rinvgamma(1,shape = a_phi0, scale = b_phi0)
  m0        <- rnorm(1, m_phi0, sqrt(V_phi0*s0))
  ###################################################################################################################################
  ###################################################################################################################################
  t1 <- Sys.time()
  for(i in 1:(NSIM*THIN + BURN)){
    
    #################################################################################################################################à
    LX     = update_lambda_gamma_cpp(z = z,
                                          rho = rho,
                                          pis=pis,k = k,n = n,
                                          MuSig2 = MuSig,
                                          m0 = m0,s0 = s0,
                                          J = J,
                                          possible_label0J = poss_lab_0J, 
                                          a = a)
    lambda = LX[,1]
    x      = LX[,2]
    #################################################################################################################################à
    rho    = update_rho_cpp(lambda = lambda,
                            a_rho = a_rho,
                            b_rho = b_rho,
                            n = n)
    ##################################################################################################################################
    u_pis = UPD_Pi_v_z2(zj = x, J = J, alpha = alpha) - 1e-8
    
    # if(  !is.finite(sum(log(1-u_pis[-J])) ))
    # {cat(paste("Not enough groups: fix alpha"))
    #   break;}
    pis   = SB_given_u2(u_pis)
    ###################################################################################################################################
    if(fixedA & i == 1){
      alpha <- 1
    }else if (fixedA==F){
      alpha  <- rgamma(1,  a_alpha + (J-1),  b_alpha -  sum(log(1-u_pis[-J])) )+1e-5
    }
    #################################################################################################################################à
    #k      = KAPPA #update_k_cpp(z = z,min_kN=min_kN,max_kN=max_kN, lambda = lambda,x = x,m1 = m1,m2 = m2,s1 = s1,s2 = s2,a = a,possible_label_k = poss_lab_k)
    #################################################################################################################################à
    m1s1   = update_msBNP_cpp(z = z,
                              x = x, 
                              lambda = lambda,
                              m0j = m01, V0j = V01,
                              a0j = a01, b0j = b01,
                              prevMuSig2 = MuSig, J = J,  
                              SIG = SIG1, k = k ,a = a)
    
    MuSig = m1s1[1:2,] 
    ACC1[,i-ind*batch] = m1s1[3,]
    ###################################################################################################################################
    z0     = z[lambda==0]
    m0s0   = update_m0s0_cpp(z0,m_phi0,V_phi0,a_phi0,b_phi0)
    m0     = m0s0[1]; s0 = m0s0[2];
    #################################################################################################################################à
    z1     = z[lambda==1];
    nj     = table(factor(x,levels=1:J))
    Avec   = MH_step_a_cpp2(previous_a = a, MuSig2 = MuSig,
                            aa = aa,bb = bb,
                            sigma_a = sigA, sub_z = z1,k = k,J = J,nj = nj)
    Avec
    a                 = Avec[1]; 
    ACCa[i-ind*batch] = Avec[2];
    ##########################################################################################
    if( ((i) %% batch) == 0){
      for(jj in 1:J){
        diag(SIG1[,,jj]) = exp(ifelse( mean(ACC1[jj,]) < optThresh, 
                                       (log(diag(SIG1[,,jj])) - min(0.01,1/sqrt(i)) ) , 
                                       (log(diag(SIG1[,,jj])) + min(0.01,1/sqrt(i)) ) ))
      }
      
      sigA       = exp(ifelse( mean(ACCa) < .44,#optThresh, 
                               (log(sigA) - min(0.01,1/sqrt(i) )) , 
                               (log(sigA) + min(0.01,1/sqrt(i)) ) ))
      ind <- ind + 1
    }
    ##########################################################################################
    if (i > BURN & ((i - BURN) %% THIN == 0)) {
      rr                      <- floor((i - BURN)/THIN)
      RHOcon[rr]   = rho
      Lcon[,rr]    = lambda
      #Kcon[rr]     = k
      Picon[,rr]    = pis
      MuSigArray[,,rr] = MuSig
      S0[rr]       = s0
      M0[rr]       = m0
      Xcon[,rr]    = x
      A[rr]        = a
      Alphacon[rr] = alpha
      DENS[,2]     = DENS[,2] +  (1-rho) * dnorm(ygrid,m0,sqrt(s0))/NSIM
      DENS[,3]     = DENS[,3] +     rho  * wd_mixmom_cppBNP_vec(z = ygrid,MuSig2 = MuSig,J = J,pis = pis,k = k,a = a)/NSIM
    }
    ##########################################################################################    
    if (i%%(verbose.step*THIN) == 0) {
      cat(paste("\n Sampling iteration: ", i, " out of ",NSIM*THIN + BURN,"---- alpha=",round(alpha,3),"\n",
                "---- M0=",round(m0,3),"---- S0=",round(s0,3),"---- a=",round(a,3),"\n"," ---- SigA=", round(sigA,5), "------\n",
                sep = "" ))
      }
    ##########################################################################################
  }
  t2 <- Sys.time()
  
  ##########################################################################################
  DENS <- cbind(DENS,DENS[,2]+DENS[,3])
  DENS <- cbind(DENS,1-DENS[,2]/DENS[,4])
  ##########################################################################################
  
  return(list(Kcon   = as.mcmc(Kcon),
              RHOcon = as.mcmc(RHOcon),
              Lcon   = apply(Lcon,1,mean), 
              MuSigArray=MuSigArray,
              Xcon   = Xcon, 
              Picon  = mcmc.list(as.mcmc(t(Picon))), 
              S0     = as.mcmc(S0), 
              M0     = as.mcmc(M0), 
              a      = as.mcmc(A),
              data   = z, 
              Alpha  = as.mcmc(Alphacon),
              DENS   = DENS,
              time   = t2-t1,
              HP     = hyperpar))
}

# 
# ###############################################################################################################
# 

SIMU2  <- function(z,J,nsim=10000,burn=10000,thin=1,HP,SS=1234){
  res1 <- Nollik_w4_BNP(z = z,   hyperpar = HP, NSIM = nsim,
                        BURN = burn,
                        J = J,
                        THIN = thin,
                        verbose.step = 5,
                        sigA = .5, fixedA = T,
                        setseed = SS,
                        optThresh = .4 )
  return(res1)
}