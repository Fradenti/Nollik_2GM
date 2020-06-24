# HP <- list(
#   KAPPA = 2,
#   aa    = 20,    bb   = 57,
#   a_rho = 1,    b_rho = 9,
#   ax    = 1,     bx   = 1, 
#   m01   = 0, V01   = 1, 
#   a01   = 3,     b01  = 1,
#   a_phi0 = 10, b_phi0 = 10, 
#   m_phi0 = 0,  V_phi0 = 1/100
# )
# 
# 
# 
# 
# Rcpp::sourceCpp(here::here("Riprova_corretta_W3_v4BNP.cpp"))
# source(here::here("Riprova_correttaW3_v4BNP3.R"))
  
  ### nota per futuro i dati definiscili tutti fuori
      library(plyr); 
      library(dbplyr)
      library(tidyverse)
      library(doParallel); library(foreach)
      ###########################################################################
       # stopCluster(cl)
        cl <- makeCluster(24)
        registerDoParallel(cl)
        
        clusterEvalQ(cl,expr = {
          Rcpp::sourceCpp(here::here("Riprova_corretta_W3_v4BNP.cpp"))
          source(here::here("Riprova_correttaW3_v4BNP3.R"))
          
          Ndata=1000
          set.seed(12345)
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
          ZZ <- replicate(50, c(rnorm(Ndata*.90,0,sqrt(1.5)),rnorm(Ndata*.05,5,1),rnorm(Ndata*.05,-5,1)))}
        )
        
        #########################################################################
      
      ############################################################## Scenario 1
        SCEN=1
        S1 <-parSapply(cl, 1:50, function(i,...) { 
          a = SIMU2(z = ZZ[,i],J = 30,nsim =5000,burn = 50000,thin = 30,HP = HP,SS = (1321*i))
          #      a = SIMU2(z = ZZ[1,,i],J = 30,nsim =10000,burn = 50000,thin = 50,HP = HP,SS = (1321*i))
        })
    save(S1,file="NS1BNPshort.RData")
    ############################################################## Scenario 2
        stopCluster(cl)
    
    
    
    
    
    
    
    
    cl2 <- makeCluster(24)
    registerDoParallel(cl2)
    library(parallel)
    clusterEvalQ(cl2,expr = {
      Rcpp::sourceCpp(here::here("Riprova_corretta_W3_v4BNP.cpp"))
      source(here::here("Riprova_correttaW3_v4BNP3.R"))      
        Ndata=1000
        set.seed(12345)
        HP <- list(
          KAPPA = 2,
          aa    = 20,    bb   = 57,
          a_rho = 1,    b_rho = 9,
          ax    = 1,     bx   = 1, 
          m01   = 0, V01   = 1, 
          a01   = 3,     b01  = 100,
          a_phi0 = 10, b_phi0 = 10, 
          m_phi0 = 0,  V_phi0 = 1/100,
          alpha_a=1,alpha_b=1
        )
        
        ZZ <- replicate(50, c(rnorm(Ndata*.90,0,1),rnorm(Ndata*.05,3,1),rnorm(Ndata*.05,-5,1.5) ) )}
      )
    ############################################################## Scenario 1
    SCEN=2
    S2 <-parSapply(cl2, 1:50, function(i,...) { 
      a = SIMU2(z = ZZ[,i],J = 30,nsim =5000,burn = 50000,thin = 30,HP = HP,SS = (1321*i))
      #      a = SIMU2(z = ZZ[1,,i],J = 30,nsim =10000,burn = 50000,thin = 50,HP = HP,SS = (1321*i))
    })
      save(S2,file="NS2BNPshort.RData")
    ############################################################## Scenario 3
    
    
    
    ############################################################## Scenario 4
    # Scenario Efron 6.1,7.1
      stopCluster(cl2)
      cl3 <- makeCluster(20)
      registerDoParallel(cl3)
      clusterEvalQ(cl3,expr = {
        Rcpp::sourceCpp(here::here("Riprova_corretta_W3_v4BNP.cpp"))
        source(here::here("Riprova_correttaW3_v4BNP3.R"))
        
        
        Ndata=1000
        set.seed(12345)
        HP <- list(
          KAPPA = 2,
          aa    = 20,    bb   = 57,
          a_rho = 1,    b_rho = 9,
          ax    = 1,     bx   = 1, 
          m01   = 0, V01   = 1, 
          a01   = 3,     b01  = 100,
          a_phi0 = 10, b_phi0 = 10, 
          m_phi0 = 0,  V_phi0 = 1/100,
          alpha_a=1,alpha_b=1
        )
        
        ZZ <- replicate(50, {
          sapply(1:Ndata,function(x) {       
            if(runif(1)<.9) {  c( rnorm(1,0,1),0) } else {
              c( rnorm(1,rnorm(1,-3,1),1),1)  }}  )})
      })
      ############################################################## Scenario 1
      SCEN=3
      S3 <-parSapply(cl3, 1:50, function(i,...) { 
        a = SIMU2(z = ZZ[1,,i],J = 30,nsim =5000,burn = 50000,thin = 30,HP = HP,SS = (1321*i))
        #      a = SIMU2(z = ZZ[1,,i],J = 30,nsim =10000,burn = 50000,thin = 50,HP = HP,SS = (1321*i))
        list(a, b =ZZ[2,,i] )
      })
    save(S3,file="NS3BNPshort.RData")
    
    ############################################################## Scenario 5
    # We generate z i ∼ N (δ i , 1), i = 1, . . . , N = 1000. 950 of the δ i were 0. The other 50 were drawn
    # (once and for all) from a Unif(2, 4) distribution.
    # Scenario 5
    # Scenario Murialidharan 6.1,7.1
    stopCluster(cl3)
    cl4 <- makeCluster(20)
    registerDoParallel(cl4)
    clusterEvalQ(cl4,expr = {
      Rcpp::sourceCpp(here::here("Riprova_corretta_W3_v4BNP.cpp"))
      source(here::here("Riprova_correttaW3_v4BNP3.R"))
      
        Ndata=1000
        set.seed(12345)
        HP <- list(
          KAPPA = 2,
          aa    = 20,    bb   = 57,
          a_rho = 1,    b_rho = 9,
          ax    = 1,     bx   = 1, 
          m01   = 0, V01   = 1, 
          a01   = 3,     b01  = 100,
          a_phi0 = 10, b_phi0 = 10, 
          m_phi0 = 0,  V_phi0 = 1/100,
          alpha_a=1,alpha_b=1
        )
        
        
        ZZ <- replicate(50, {
          sapply(1:Ndata,function(x) {       
            uu <- runif(1) 
            if(uu<.90) {
              c( rnorm(1,0,1),0) 
            } else if(uu<.95 ){
              c( rnorm(1,runif(1,2,4),1),1)  }
            else {
              c( rnorm(1,runif(1,-4,-2),1),1)}
          }  )
      })})
    ############################################################## Scenario 1
    SCEN=4
    S4 <-parSapply(cl4, 1:50, function(i,...) { 
      a = SIMU2(z = ZZ[1,,i],J = 30,nsim =5000,burn = 50000,thin = 30,HP = HP,SS = (1321*i))
      #      a = SIMU2(z = ZZ[1,,i],J = 30,nsim =10000,burn = 50000,thin = 50,HP = HP,SS = (1321*i))
      list(a, b =ZZ[2,,i] )
    })
    save(S4,file="NS4BNPshort.RData")
