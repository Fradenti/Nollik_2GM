###################
# Post-process hierarchical functions
###################

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


##############################################################################################################
f_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon[,j])) * dnorm(x,mean(res1$M0),sqrt(mean(res1$S0))) +  
    mean(res1$RHOcon[,j])  * wd_mixmom_theta(z = x,theta = mean(res1$THcon[,j]),
                                             m1 =    mean(res1$M1con),
                                             m2 =    mean(res1$M2con),
                                             s1 =    mean(res1$S1con),
                                             s2 =    mean(res1$S2con),
                                             k  =    round(mean(res1$Kcon)),
                                             a  =    mean(res1$a))
}

pif0_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon[,j]))*dnorm(x,mean(res1$M0),sqrt(mean(res1$S0)))
}

locdfr_postj <- function(x,res1,j) pif0_postj(x,res1,j)/f_postj(x,res1,j)
locdfr_post_kj <- function(x,res1,kkk,j) pif0_postj(x,res1,j)/f_postj(x,res1,j)-kkk


##############################################################################################################
##############################################################################################################
f_postji <- function(x,res1,j,i){
  (1-res1$RHOcon[i,j])*dnorm(x,res1$M0[i],sqrt(res1$S0[i])) +  
    res1$RHOcon[i,j] * wd_mixmom_theta(z = x,
                                       theta = res1$THcon[i,j],
                                       m1 =     res1$M1con[i],
                                       m2 =     res1$M2con[i],
                                       s1 =     res1$S1con[i],
                                       s2 =     res1$S2con[i],
                                       k  =     res1$Kcon[i],
                                       a  =     res1$a[i])
}

pif0_postji <- function(x,res1,j,i){
  (1-res1$RHOcon[i,j])*dnorm(x,(res1$M0[i]),sqrt(res1$S0[i]))
}

locdfr_postji <- function(x,res1,j,i) pif0_postji(x,res1,j,i)/f_postji(x,res1,j,i)
locdfr_post_kji <- function(x,res1,kkk,j,i) pif0_postj(x,res1,j,i)/f_postj(x,res1,j,i)-kkk




MCMC_locfdr_extractor <- function(res,J){
  MCMC_locfdr <- matrix(NA,nsim,100)
  for(i in 1:nsim){
    x <- seq(-10,10,length.out = 100)
    MCMC_locfdr[i,] <- locdfr_postji(x,res1 = res,j = J,i=i)
    cat(paste(i,"\n"))
  }
  return(MCMC_locfdr)
}


##################################################### BAYESIAN FDR -- > pick the threshold
vannucci <- function(k,p,alpha){
  fdr <- sum((1-p)*(p>k))/sum(p>k) 
  fdr  - alpha
}
choose_threshold <- function(prob,AL){
  uniroot(function(x) vannucci(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
}
# beware: my function picks as threshold the last NON-valid obs
# Michele's, Simon's ecc's one picks the first valid one

acc <- function(preds,actual) mean(preds==actual)