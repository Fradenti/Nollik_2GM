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




  Rcpp::sourceCpp(here::here("W3/Riprova_corretta_W3_v3.cpp"))
  source(here::here("W3/Riprova_correttaW3_v3.R"))
  
  ### nota per futuro i dati definiscili tutti fuori
  library(plyr); 
  library(dbplyr)
  library(tidyverse)
  library(doParallel); library(foreach)
  ###########################################################################
    stopCluster(cl)
    cl <- makeCluster(24)
    registerDoParallel(cl)
    
    clusterEvalQ(cl,expr = {
      Rcpp::sourceCpp(here::here("W3/Riprova_corretta_W3_v3.cpp"))
      source(here::here("W3/Riprova_correttaW3_v3.R"))
      
      
      Ndata=1000
      set.seed(12345)
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
      ZZ <- replicate(50, c(rnorm(Ndata*.90,0,sqrt(1.5)),rnorm(Ndata*.05,5,1),rnorm(Ndata*.05,-5,1)))}
    )
    
    #########################################################################
  
  ############################################################## Scenario 1
    SCEN=1
    S1 <-parSapply(cl, 1:50, function(i,...) { 
      SIMU2(ZZ[,i],nsim =5000,burn = 50000,thin = 30,HP = HP,SS = (1321*i))})
    save(S1,file="NS1.RData")
    ############################################################## Scenario 2
        stopCluster(cl)
    
    
    
    
    
    
    
    
    cl2 <- makeCluster(20)
    registerDoParallel(cl2)
    library(parallel)
    clusterEvalQ(cl2,expr = {
      Rcpp::sourceCpp(here::here("W3/Riprova_corretta_W3_v3.cpp"))
      source(here::here("W3/Riprova_correttaW3_v3.R"))
      
        Ndata=1000
        set.seed(12345)
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
        ZZ <- replicate(50, c(rnorm(Ndata*.90,0,1),rnorm(Ndata*.05,3,1),rnorm(Ndata*.05,-5,1.5) ) )}
      )
    ############################################################## Scenario 1
    SCEN=2
    S2 <-parSapply(cl2, 1:50, function(i,...) { 
      SIMU2(ZZ[,i],nsim =5000,burn = 50000,thin=30,HP = HP,SS = 1321*i)})
    save(S2,file="NS2.RData")
    ############################################################## Scenario 3
    
    
    
    ############################################################## Scenario 4
    # Scenario Efron 6.1,7.1
      stopCluster(cl2)
      cl3 <- makeCluster(20)
      registerDoParallel(cl3)
      clusterEvalQ(cl3,expr = {
        Rcpp::sourceCpp("Riprova_corretta_W3_v3.cpp")
        source("Riprova_correttaW3_v3.R")
        
        
        Ndata=1000
        set.seed(12345)
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
        ZZ <- replicate(50, {
          sapply(1:Ndata,function(x) {       
            if(runif(1)<.9) {  c( rnorm(1,0,1),0) } else {
              c( rnorm(1,rnorm(1,-3,1),1),1)  }}  )})
      })
      ############################################################## Scenario 1
      SCEN=3
      S3 <-parSapply(cl3, 1:50, function(i,...) { 
        a   = SIMU2(ZZ[1,,i],nsim =5000,burn = 50000,COR =2,thin=30,HP = HP,SS = 1321*i)
        list(a, b =ZZ[2,,i] )
      })
    save(S3,file="NS3.RData")
    
    ############################################################## Scenario 5
    # We generate z i ∼ N (δ i , 1), i = 1, . . . , N = 1000. 950 of the δ i were 0. The other 50 were drawn
    # (once and for all) from a Unif(2, 4) distribution.
    # Scenario 5
    # Scenario Murialidharan 6.1,7.1
    stopCluster(cl3)
    cl4 <- makeCluster(20)
    registerDoParallel(cl4)
    clusterEvalQ(cl4,expr = {
      Rcpp::sourceCpp("Riprova_corretta_W3_v3.cpp")
      source("Riprova_correttaW3_v3.R")
      
        Ndata=1000
        set.seed(12345)
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
      a = SIMU2(ZZ[1,,i],nsim =5000,burn = 50000,COR =3,thin=30,HP = HP,SS = 1321*i)
      list(a, b =ZZ[2,,i] )
    })
    save(S4,file="Nw3S4.RData")
