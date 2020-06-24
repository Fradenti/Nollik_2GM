functions{
  
//  real Norm_pdf(real x, real m, real s2){
//    return   exp(-(x-m)*(x-m)/(2*s2))/sqrt(2*pi()*s2);
//  }
  
  real w(real x, real csi){
    return    1 - exp( - exp( 2 * (log(x*x)-log(csi*csi) ))); //k=2
  }
  
  real per_integrale(real x, real csi){
    return   exp( - exp( 2 * (log(x*x)-log(csi*csi) ))); //k=2
  }
  
  
  real integrand(real x, real xc, real[] theta, real[] x_r, int[] x_i){
    real logVAL = log(per_integrale( x, theta[6])) + 
                  log_sum_exp( log( 1-theta[5] ) + normal_lpdf( x | theta[1], sqrt(theta[2])) ,  
                  log(theta[5])                  + normal_lpdf( x | theta[3], sqrt(theta[4])));
    return  exp(logVAL);
  }
  
  
  real f1_lpdf(real z, real[] theta, real Kcal){
    
    return log(w(z,theta[6])) - 
            log(1-Kcal) + 
            log_sum_exp( log(1-theta[5]) +  normal_lpdf( z | theta[1], sqrt(theta[2]) ) ,  
                         log(theta[5]) + normal_lpdf( z |  theta[3], sqrt(theta[4]))) ;
  }

}


data {
  real a_rho;   real b_rho;
  real a_alpha; real b_alpha;
  real m_f0;    real <upper=0> m_f1;  real <lower=0> m_f2;
  real V_f0;    real V_f1;            real V_f2;
  real a_IG0;   real b_IG0;
  real a_IG1;   real b_IG1;
  real a_IG2;   real b_IG2;
  real a_CSI;   real b_CSI;
  int<lower=1> N;          // number of data points
  real z[N];               // observations
  real a;
  real b;
 }

transformed data {
  real x_r[0];
  int  x_i[0];
}

parameters {
  real <lower=0, upper=1> alpha;   
  real <lower=0, upper=1> rho;
  real mu0;     real <lower = 0> sig0;
  real mu1;
  real mu2;
  real <lower = 0> sig1;
  real <lower = 0> sig2;
  real csi;
  }

transformed parameters{
  real theta[6];
  theta[1] = mu1;
  theta[3] = mu2;
  theta[2] = sig1;
  theta[4] = sig2;
  theta[5] = alpha;
  theta[6] = csi;
   }

model {
  real  Kcal; 
  
  alpha  ~ beta(a_alpha,b_alpha);
  rho    ~ beta(a_rho,b_rho);
  
  csi    ~ inv_gamma(a_CSI, b_CSI);
  sig0   ~ inv_gamma(a_IG0, b_IG0); 
  sig1   ~ inv_gamma(a_IG1, b_IG1); 
  sig2   ~ inv_gamma(a_IG2, b_IG2); 
  
  mu0    ~ normal(m_f0, sqrt(V_f0 * sig0) ); 
  mu1    ~ normal(m_f1, sqrt(V_f1 * theta[2]) ); 
  mu2    ~ normal(m_f2, sqrt(V_f2 * theta[4]) );

  Kcal =  integrate_1d( integrand, a, b, theta, x_r, x_i, 0.001);

  for (n in 1:N) {
    target +=  log_sum_exp(log(1-rho) + normal_lpdf( z[n] | mu0 , sqrt(sig0)) , 
                           log(rho)   + f1_lpdf( z[n] | theta, Kcal));
    }
  
}

