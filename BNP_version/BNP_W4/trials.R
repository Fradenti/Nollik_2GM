###  Cosa osservato di problematico:
#    Se J troppo piccolo salta, una volta aumentato J il problema scompare
#    Come salta? il parametro a schizza in alto e la f1 scompare

set.seed(19922598)
z <- c(rnorm(200,0,1), rnorm(50,-2), rnorm(50,2,4))
HP <- list(
  KAPPA = 2,
  aa    = 20,    bb   = 57,
  a_rho = 1,    b_rho = 9,
  ax    = 1,     bx   = 1, 
  m01   = 0, V01   = 100, 
  a01   = 1,     b01  = 1,
  a_phi0 = 10, b_phi0 = 10, 
  m_phi0 = 0,  V_phi0 = 1/100,
  alpha_a=1,alpha_b=1
)
hist(z)
d <- z#[1:1000]
res1 <- Nollik_w4_BNP(z = d,   
                      hyperpar = HP, 
                      NSIM = 5000,
                      BURN = 5000,
                      J = 20,
                      THIN = 1,
                      verbose.step = 50,
                      sigA = .1, 
                      fixedA = T,
                      setseed = 1, # 
                      optThresh = .4 )
plot(res1$data)
plot(res1$DENS[,2])
plot(res1$DENS[,3])
plot((c(res1$RHOcon)),type="l")

hist(res1$data,freq = F,breaks = 40)
lines(res1$DENS[,2]~res1$DENS[,1])
hist(res1$data,freq = F,breaks = 40)
lines(res1$DENS[,3]~res1$DENS[,1])
hist(res1$data,freq = F,breaks = 25)
lines(res1$DENS[,4]~res1$DENS[,1])
hist(res1$data,freq = F,breaks = 40)
lines(res1$DENS[,5]~res1$DENS[,1])

plot((c(res1$RHOcon)),type="l")
plot((c(res1$S0)),type="l")
plot((c(res1$M0)),type="l")
plot((c(res1$a)))
plot(((res1$a)))
plot(((res1$Alpha)))

res1$Xcon[,1999]
res1$Picon[2000,]

plot(res1$Picon[,1])
plot(res1$Picon[500,][[1]])

plot(t(res1$MuSigArray[,,900]))

