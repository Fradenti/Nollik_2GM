// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppNumerical.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
using namespace Numer;
using namespace Rcpp;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::colvec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::colvec mvrnormArma_simpler(int n, arma::colvec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return mu +  arma::chol(sigma).t() * Y.t();
}


/* 1st block: Function for numeric integration to compute the expected value \mathcal{K} at the denominator
 * The weight is w_3 = 1 - exp( -(x/a)^2 )
 */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class NormConst_part1: public Func
{
private:
  double m;
  double s;
  double a;
  int k;
public:
  NormConst_part1(double m_,  double s_, double a_, int k_) : m(m_), s(s_), a(a_), k(k_) {}
  double operator()(const double& x) const
  {
    return exp( - exp( k * (log(x*x)-log(a*a) ))) * R::dnorm(x,m,sqrt(s),0) ;
  }
};

// [[Rcpp::export]]
double integrate_normconst1( double m, double s, double a, int k, double low, double upp)
{
  NormConst_part1 f(m,s,a,k);
  double err_est;
  int err_code;
  const double res = integrate(f, low, upp, err_est, err_code);
  return res;
}

// [[Rcpp::export]]
double Cost_xW3 (int k, double m, double s, double a, double low=-50, double upp=50){
  return 1-integrate_normconst1(m,s,a,k,low,upp);
}

// [[Rcpp::export]]
double Cost_xW3_Exact (double m, double s, double a){  //la si usa solo se 
  return 1 - exp( -(m*m)/(2*s+a*a))/sqrt(s * (2/(a*a)+1/(s)) );
}


/*Easy inverse gamma density*/
// [[Rcpp::export]]
double dinvgamma_cpp( double x, double a, double b){
  return  R::dgamma(1/x,a,1/b,0)/(x*x);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Likelihood function for single observation. It needs the two alternative constants as input.
 */
// [[Rcpp::export]]
double wd_mixmom_cpp2(double z, int x, double m1, double m2,
                     int k, double s1, double s2, double a, bool logscale,
                     arma::colvec cost){
  double g = ( 1 - exp( - exp( k * ( log(z*z)-log(a*a) ))) )  * 
    ( (1-x) * R::dnorm(z,m1,sqrt(s1),0) + x * R::dnorm(z,m2,sqrt(s2),0 ) );
  double RR = log(g) - log(cost(x));
  if(logscale){
    return(RR);
  }else{
    return(exp(RR));
  }}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * MH on NIG1 and NIG2
 */

/* 1st NIG: on the negatives: logposterior*/
// [[Rcpp::export]]
double log_NIGls_post_cpp1(double m, double ls, double m0, double k0, double a0, double b0, int n10, int k,double a){
  return log(dinvgamma_cpp(exp(ls),a0,b0)) + 
         R::dnorm(m, m0, sqrt(exp(ls)/k0), 1) + 
         (ls) - n10 * log(Cost_xW3( k, m, exp(ls),a )) 
         + log(m<0);
}

/* 2nd NIG: on the positives: logposterior*/
// [[Rcpp::export]]
double log_NIGls_post_cpp2(double m, double ls, double m0, double k0, double a0, double b0, int n11, int k,double a){
  return log(dinvgamma_cpp(exp(ls),a0,b0)) + 
         R::dnorm(m, m0, sqrt(exp(ls)/k0), 1) + 
         (ls) - n11 * log(Cost_xW3( k, m, exp(ls),a )) 
         + log(m>0);
}


/* 1st NIG: on the positives: MH step*/
// [[Rcpp::export]]
arma::colvec MH_step_m_ls_cpp1(arma::colvec previous_ms, 
                               double m0, double k0, double a0, double b0, 
                               arma::mat SIGMA, int n10 , int k, double a){
  
  double p_ls = log(previous_ms[1]);
  arma::colvec previous_mls(2); 
  previous_mls[0]= previous_ms[0]; previous_mls[1]=p_ls;
  arma::colvec proposal = mvrnormArma(1, previous_mls, SIGMA).t();
  arma::colvec uno(1); uno.fill(1); 
  arma::colvec zer(1); zer.fill(0);
  
  arma::colvec proposal2 = arma::join_cols(proposal,uno);
  arma::colvec previous2 = arma::join_cols(previous_mls,zer);
  
  double logratio = log_NIGls_post_cpp1(proposal[0],     proposal[1],     m0,  k0, a0, b0, n10, k, a) - 
                    log_NIGls_post_cpp1(previous_mls[0], previous_mls[1], m0 , k0, a0, b0, n10, k, a);
  if( R::runif(0,1) < exp(logratio) ){
    proposal2[1] = exp(proposal2[1]);
    return(proposal2);
  }else{
    previous2[1] = exp(previous2[1]);
    return(previous2);
  }
}

/* 2nd NIG: on the positives: MH step*/
// [[Rcpp::export]]
arma::colvec MH_step_m_ls_cpp2(arma::colvec previous_ms, 
                               double m0, double k0, double a0, double b0, 
                               arma::mat SIGMA, int n11 , int k, double a){
  
  double p_ls = log(previous_ms[1]);
  arma::colvec previous_mls(2); 
  previous_mls[0]= previous_ms[0]; previous_mls[1]=p_ls;
  arma::colvec proposal = mvrnormArma(1, previous_mls, SIGMA).t();
  arma::colvec uno(1); uno.fill(1); 
  arma::colvec zer(1); zer.fill(0);
  
  arma::colvec proposal2 = arma::join_cols(proposal,uno);
  arma::colvec previous2 = arma::join_cols(previous_mls,zer);
  
  double logratio = log_NIGls_post_cpp2(proposal[0], proposal[1], m0,  k0, a0, b0, n11, k, a) - 
                    log_NIGls_post_cpp2(previous_mls[0], previous_mls[1], m0 , k0, a0, b0, n11, k, a);
  if( R::runif(0,1) < exp(logratio) ){
    proposal2[1] = exp(proposal2[1]);
    return(proposal2);
  }else{
    previous2[1] = exp(previous2[1]);
    return(previous2);
  }
}


/* 1st NIG: on the positives: UPDATE FUNCTION*/
// [[Rcpp::export]]
arma::colvec update_m1s1_cpp(arma::colvec z, arma::colvec x, double m0j, double V0j, double a0j, double b0j,
                             arma::colvec lambda, double prevMj, double prevSj, arma::mat SIG, int k, double a){
  
  arma::colvec MS(2); MS(0) = prevMj; MS(1) = prevSj;
  arma::uvec ind_lambda1_xj = find(lambda==1 && x==0);
  arma::colvec selected_z  = z.elem(ind_lambda1_xj);
  int                  n1j = selected_z.n_rows;
  double Vn = 1/(n1j+1/V0j),
    an = a0j + n1j/2, mn=0.0, bn=0.0;
  
  if(n1j == 0){
    mn = Vn * (m0j/V0j);
    bn = b0j;
  }else{
    mn = Vn * (m0j/V0j + n1j * mean(selected_z));
    bn = b0j + .5*( m0j*m0j/V0j + accu(selected_z%selected_z) - mn*mn/Vn   );
  }
  return  MH_step_m_ls_cpp1(MS,
                            mn, 1/Vn,
                            an, bn,
                            SIG, n1j, 
                            k,a);
}

/* 2nd NIG: on the positives: UPDATE FUNCTION*/
// [[Rcpp::export]]
arma::colvec update_m2s2_cpp(arma::colvec z, arma::colvec x, double m0j, double V0j, double a0j, double b0j,
                             arma::colvec lambda, double prevMj, double prevSj, arma::mat SIG, int k, double a){
  
  arma::colvec MS(2); MS(0) = prevMj; MS(1) = prevSj;
  arma::uvec ind_lambda1_xj = find(lambda==1 && x==1);
  arma::colvec selected_z  = z.elem(ind_lambda1_xj);
  int                  n1j = selected_z.n_rows;
  
  double Vn = 1/(n1j+1/V0j),
    an = a0j + n1j/2, mn=0.0, bn=0.0;
  
  if(n1j == 0){
    mn = Vn * (m0j/V0j);
    bn = b0j;
  }else{
    mn = Vn * (m0j/V0j + n1j * mean(selected_z));
    bn = b0j + .5*( m0j*m0j/V0j + accu(selected_z%selected_z) - mn*mn/Vn   );
  }
  return  MH_step_m_ls_cpp2(MS,
                            mn, 1/Vn,
                            an, bn,
                            SIG, n1j, 
                            k,a);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// MH on parameter a of w3
// [[Rcpp::export]]
double LOG_post_la2(double la,  
                    double m1, double s1, 
                    double m2, double s2,  
                    int k, arma::colvec sub_z, 
                    double aa, double bb,
                    int   n10, int n11){
  arma::colvec costs(2);
  costs(0) = log(Cost_xW3( k, m1, s1, exp(la) ));
  costs(1) = log(Cost_xW3( k, m2, s2, exp(la) ));
      return log( dinvgamma_cpp(exp(la),aa,bb) ) + 
             accu( log( 1 - exp( -pow( sub_z/exp(la), 2*k  )))) + 
             (la) - n10*costs(0) - n11*costs(1) ;
}


// [[Rcpp::export]]
arma::colvec MH_step_a_cpp2(double previous_a, 
                            double m1, double s1, double m2, double s2, 
                            double aa, double bb,
                            double sigma_a, arma::colvec sub_z, int k,int n11,int n12){
  
  arma::colvec prev_a(2), prop_a(2); 
  prev_a[0] = log(previous_a); prev_a[1]=0.0;
  prop_a[0] = R::rnorm(prev_a[0], sqrt(sigma_a)); prop_a[1] = 1;
  
  double logratio = LOG_post_la2(prop_a[0], m1,s1, m2,s2, k,sub_z, aa, bb, n11,n12 )-
                    LOG_post_la2(prev_a[0], m1,s1, m2,s2, k,sub_z, aa, bb ,n11,n12 );
  
  if( R::runif(0,1) < exp(logratio) ){
    prop_a[0] = exp(prop_a[0]);
    return(prop_a);
  }else{
    prev_a[0] = exp(prev_a[0]);
    return(prev_a);
  }
}







/*UPDATE RHO*/
//#3
// [[Rcpp::export]]
double update_rho_cpp(arma::colvec lambda, double a_rho, double b_rho, int n){
  return R::rbeta(a_rho+accu(lambda),n-accu(lambda)+b_rho);
}
//################################
//################################
//#4
/*UPDATE THETA*/
// [[Rcpp::export]]
double update_theta_cpp(arma::colvec x, double ax, double bx){
  int ngamma1 = accu(x != -1);
  return R::rbeta(ax + accu(x==1), ngamma1 - accu(x==1) + bx);
}
//################################
//################################
//#7
/*UPDATE JOINTLY UPDATE lambda, x*/
// [[Rcpp::export]]
arma::mat update_lambda_x_cpp_parallel(arma::colvec z, double rho, double theta,
                                       int k, int n, 
                                       double m1, double m2, double m0,
                                       double s1, double s2, double s0,
                                       arma::colvec possible_label123,  
                                       double a,  int cores){
  arma::mat newLambdaGamma(n,2);
  arma::colvec costs(2);
  /*I need to compute them again, since a has been updated*/
  costs(0) = Cost_xW3( k, m1, s1, a );
  costs(1) = Cost_xW3( k, m2, s2, a );
  omp_set_num_threads(cores);
  #pragma omp parallel for schedule(static)
  for(int i=0; i<n ; i++){
    arma::colvec probs(3);
    probs[0] = (1-rho) * R::dnorm(z[i],m0,sqrt(s0),0);
    probs[1] = rho * wd_mixmom_cpp2(z[i],0, m1, m2, k, s1, s2,a, 0, costs) * (1-theta);
    probs[2] = rho * wd_mixmom_cpp2(z[i],1, m1, m2, k, s1, s2,a, 0, costs) * (theta);
    int INT = RcppArmadillo::sample(possible_label123, 1, TRUE, probs)[0];
    
    if(INT==1){
      newLambdaGamma(i,0) = 0;
      newLambdaGamma(i,1) = -1;
    }else if(INT==2){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 0;
    }else if(INT==3){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 1;
    }
    
  }
  return(newLambdaGamma);
}



/* UPDATE M0S0 */

// [[Rcpp::export]]
double rt_cpp(double nu, double lambda){
  double TAA;
  TAA = R::rnorm(lambda,1.0) / exp(.5*log( R::rchisq(nu)/nu ));
  return(TAA);
}



// [[Rcpp::export]]
arma::colvec update_m0s0_cpp(arma::colvec sub_z0, double m_phi0, double V_phi0, double a_phi0, double b_phi0){
  arma::colvec res(2);
  int   n_phi0 = sub_z0.n_rows;
  
  double Vn = 1/(n_phi0+1/V_phi0),
    an = a_phi0 + n_phi0/2, mn=0.0, bn=0.0;
  
  if(n_phi0 == 0){
    mn = Vn * (m_phi0/V_phi0);
    bn = b_phi0;
  }else{
    mn = Vn * (m_phi0/V_phi0 + n_phi0 * mean(sub_z0)); 
    bn = b_phi0 + .5*( m_phi0*m_phi0/V_phi0 + accu(sub_z0%sub_z0) - mn*mn/Vn   );
  }
  
  res(1) = 1/ rgamma(1,an,1/(bn))[0];
  res(0) = rt_cpp(2*an, 0.0) * pow( ( bn/((1/Vn) * an)) , 0.5) + mn;
  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////777




//NESTED CASE



/*UPDATE JOINTLY UPDATE lambda, x*/
/*
//#7
// [[Rcpp::export]]
arma::mat update_lambda_x_HIER(arma::colvec z, arma::colvec rho, 
                               arma::colvec zg,
                               double theta,
                                       int k, int n, 
                                       double m1, double m2, double m0,
                                       double s1, double s2, double s0,
                                       arma::colvec possible_label123,  
                                       double a,  int cores){
  arma::mat newLambdaGamma(n,2);
  arma::colvec costs(2);
  //I need to compute them again, since a has been updated/
  costs(0) = Cost_xW3( k, m1, s1, a );
  costs(1) = Cost_xW3( k, m2, s2, a );
  omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
  for(int i=0; i<n ; i++){
    arma::colvec probs(3);
    probs[0] = (1-rho[zg[i]-1]) * R::dnorm(z[i],m0,sqrt(s0),0);
    probs[1] = rho[zg[i]-1] * wd_mixmom_cpp2(z[i],0, m1, m2, k, s1, s2,a, 0, costs) * (1-theta);
    probs[2] = rho[zg[i]-1] * wd_mixmom_cpp2(z[i],1, m1, m2, k, s1, s2,a, 0, costs) * (theta);
    int INT = RcppArmadillo::sample(possible_label123, 1, TRUE, probs)[0];
    
    if(INT==1){
      newLambdaGamma(i,0) = 0;
      newLambdaGamma(i,1) = -1;
    }else if(INT==2){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 0;
    }else if(INT==3){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 1;
    }
    
  }
  return(newLambdaGamma);
}
*/



// [[Rcpp::export]]
arma::colvec update_rho_cpp_nest(arma::colvec lambda, 
                                 arma::colvec z_g, 
                                 int nG, arma::colvec freq_each_g,
                                 double a_rho, double b_rho, int n){ //sto supponendo stessa prior per tutti i gruppi
  arma::colvec newRho(nG);
  for(int j=0; j<nG; j++){
    arma::uvec ind = find(z_g==(j+1));
    newRho(j) = R::rbeta( a_rho+accu(lambda.elem(ind)), freq_each_g[j]-accu(lambda.elem(ind))+b_rho);
  }
  return newRho;
}













////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////ù//
// Suppose now also theta varies : each group allocates a different proportion of relevant obeservation to the 1st or 2nd component

//#7
/*UPDATE JOINTLY UPDATE lambda, x*/
// [[Rcpp::export]]
arma::mat update_lambda_x_HIER2(arma::colvec z, 
                               arma::colvec rho, 
                               arma::colvec zg,
                               arma::colvec theta,
                               int k, int n, 
                               double m1, 
                               double m2, 
                               double m0,
                               double s1, 
                               double s2, 
                               double s0,
                               arma::colvec possible_label123,  
                               double a,  int cores){
  arma::mat newLambdaGamma(n,2);
  arma::colvec costs(2);
  /*I need to compute them again, since a has been updated*/
  costs(0) = Cost_xW3( k, m1, s1, a );
  costs(1) = Cost_xW3( k, m2, s2, a );
  omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
  for(int i=0; i<n ; i++){
    arma::colvec probs(3);
    probs[0] = (1-rho[zg[i]-1]) * R::dnorm(z[i],m0,sqrt(s0),0);
    probs[1] = rho[zg[i]-1] * wd_mixmom_cpp2(z[i],0, m1, m2, k, s1, s2,a, 0, costs) * (1-theta[zg[i]-1]);
    probs[2] = rho[zg[i]-1] * wd_mixmom_cpp2(z[i],1, m1, m2, k, s1, s2,a, 0, costs) * (theta[zg[i]-1]);
    int INT = RcppArmadillo::sample(possible_label123, 1, TRUE, probs)[0];
    
    if(INT==1){
      newLambdaGamma(i,0) = 0;
      newLambdaGamma(i,1) = -1;
    }else if(INT==2){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 0;
    }else if(INT==3){
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = 1;
    }
    
  }
  return(newLambdaGamma);
}

// [[Rcpp::export]]
arma::colvec update_theta_cpp_nest(arma::colvec x, 
                                 arma::colvec z_g, 
                                 int nG, arma::colvec freq_each_g,
                                 double ax, double bx){ //sto supponendo stessa prior per tutti i gruppi
  arma::colvec newX(nG);
  for(int j=0; j<nG; j++){
    int n0      = accu(z_g==(j+1) && x == 0);
    int n1      = accu(z_g==(j+1) && x == 1);
    newX(j) = R::rbeta( ax+n1, n0+bx);
  }
  return newX;
}
