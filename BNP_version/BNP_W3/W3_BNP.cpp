// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppNumerical.h>

using namespace Numer;
using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double log_weight_function(double x, double a, int k) {
  return - exp( k * (log(x*x)-log(a*a) )) ;  // qui non assegno 1-, è solo parte di weight function
}
// [[Rcpp::export]]
arma::colvec log_weight_function_vec(arma::colvec x, double a, int k) {
  return - exp( k * (log(x%x)-log(a*a) ));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::colvec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
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
    return exp(log_weight_function(x,a,k)) * R::dnorm(x,m,sqrt(s),0) ;
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
double wd_mixmom_cppBNP(double z, int x, // x goes from 1 to J 
                        arma::mat MuSig2,
                        int k, double a, bool logscale,
                        arma::colvec cost){
  double g = log( 1 - exp( log_weight_function(z,a,k) ))  + 
              R::dnorm(z, MuSig2(0,x-1), sqrt(MuSig2(1,x-1)),1);
  double RR = (g) - log(cost(x-1));
  if(logscale){
    return(RR);
  }else{
    return(exp(RR));
  }}

/*
// [[Rcpp::export]]
arma::colvec wd_mixmom_cppBNP_vec_OLDBUGGED(arma::colvec z,  // x goes from 1 to J 
                        arma::mat MuSig2,int J,arma::colvec pis,
                        int k, double a){
  arma::colvec g(z.n_elem); g.fill(0);
  for(int j=0;j<J;j++){
    
  g +=  arma::normpdf(z,MuSig2(0,j),sqrt(MuSig2(1,j))) * 
        pis(j) / Cost_xW3(k, MuSig2(0,j), MuSig2(1,j), a);
  }
    return(( 1 - exp( log_weight_function_vec(z,a,k) ))  % g);
}
*/
// [[Rcpp::export]]
arma::colvec wd_mixmom_cppBNP_vec_NEW(arma::colvec z,  // x goes from 1 to J 
                                  arma::mat MuSig2,int J,arma::colvec pis,
                                  int k, double a){
  arma::colvec g1(z.n_elem); g1.fill(0);
  double g2=0.0;
  for(int j=0;j<J;j++){
    
    g1 +=  arma::normpdf(z,MuSig2(0,j),sqrt(MuSig2(1,j))) * pis(j) ;
    g2 +=  Cost_xW3(k, MuSig2(0,j), MuSig2(1,j), a) * pis(j) ;
  
  }
  
  return(( 1 - exp( log_weight_function_vec(z,a,k) ))  % g1/g2);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * MH on NIG_x
 */

/// [[Rcpp::export]]
double log_NIGls_post_cpp1(double m, double ls, double m0, double k0, double a0, double b0, int n10, int k,double a){
  return log(dinvgamma_cpp(exp(ls),a0,b0)) + 
         R::dnorm(m, m0, sqrt(exp(ls)/k0), 1) + 
         (ls) - n10 * log(Cost_xW3( k, m, exp(ls),a ));
}


/* 1st NIG: on the positives: MH step*/
// [[Rcpp::export]]
arma::colvec MH_step_m_ls_cpp(arma::colvec previous_ms, 
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


/* 1st NIG: on the positives: UPDATE FUNCTION*/
// [[Rcpp::export]]
arma::mat update_msBNP_cpp(arma::colvec z, arma::colvec x, double m0j, double V0j, double a0j, double b0j,
                             arma::colvec lambda, arma::mat prevMuSig2, arma::cube SIG, int k, double a, int J){
  
  arma::mat updatedMuSig2(3,J);
  for(int jj=0; jj<J; jj++){
  arma::uvec ind_lambda1_xj = find(x==(jj+1));//find(lambda==1 && x==(j+1));
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
  
  updatedMuSig2.col(jj) = MH_step_m_ls_cpp( prevMuSig2.col(jj),
                            mn, 1/Vn,
                            an, bn,
                            SIG.slice(jj), n1j, 
                            k,a);
  }
  return updatedMuSig2;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// MH on parameter a of w3
// [[Rcpp::export]]
double LOG_post_la2(double la,  
                    arma::mat MuSig2,  //(2xJ)
                    int k, arma::colvec sub_z, 
                    double aa, double bb,
                    arma::colvec nj, int J){
  arma::colvec lcosts_nj(J);
  for(int j=0; j<J; j++){
  if(nj(j)==0){ 
    lcosts_nj(j)=0; 
    }else{
    lcosts_nj(j) = log(Cost_xW3( k, MuSig2(0,j), MuSig2(1,j), exp(la) ))*nj(j);
    }}
      return log( dinvgamma_cpp(exp(la),aa,bb) ) +
        // attenzione qui, non ho cambiato a generica log_weight_function per comodità
             accu( log( 1 - exp( -pow( sub_z/exp(la), 2*k  )))) + 
             (la) - accu(lcosts_nj);
}


// [[Rcpp::export]]
arma::colvec MH_step_a_cpp2(double previous_a, arma::mat MuSig2, int J, arma::colvec nj,
                            double aa, double bb,
                            double sigma_a, arma::colvec sub_z, int k){
  
  arma::colvec prev_a(2), prop_a(2); 
  prev_a[0] = log(previous_a); prev_a[1]=0.0;
  prop_a[0] = R::rnorm(prev_a[0], sqrt(sigma_a)); prop_a[1] = 1;
  
  double logratio = LOG_post_la2(prop_a[0], MuSig2, k,sub_z, aa, bb, nj,J )-
                    LOG_post_la2(prev_a[0], MuSig2, k,sub_z, aa, bb, nj,J );
  
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
/*UPDATE PIS*/
// [[Rcpp::export]]
arma::colvec UPD_Pi_v_z(arma::colvec zj, 
                         int J, double alpha){
   arma::colvec v_z(J), pi_z(J);
   for(int j=0; j<J; j++){
     v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
   }
   return v_z;  
   }

// [[Rcpp::export]]
arma::colvec UPD_Pi_v_z2(arma::colvec zj, 
                        int J, double alpha){
  arma::colvec v_z(J);
  for(int j=0; j<(J-1); j++){
    v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
  }
  v_z[J-1]=1.; 
  return v_z; 
}



  
// [[Rcpp::export]] 
arma::colvec SB_given_u2(arma::colvec V) {
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2   =arma::shift(cumprod(mV),+1); 
  mV2(0) = 1;
  return(V%mV2);
}
  
  
  
  
  
//################################
//################################
//#7
/*UPDATE JOINTLY UPDATE lambda, x*/
// [[Rcpp::export]]
arma::mat update_lambda_x_cpp_parallel(arma::colvec z, double rho, 
                                       arma::colvec pis,
                                       int k, int n, 
                                       arma::mat MuSig2,
                                       double m0,
                                       double s0, int J,
                                       arma::colvec possible_label0J,  
                                       double a){
  arma::mat newLambdaGamma(n,2);
  arma::colvec costs(J);
  
  
  /*I need to compute them again, since a has been updated*/
  for(int ll=0; ll<J ; ll++){
  costs(ll) = Cost_xW3( k, MuSig2(0,ll), MuSig2(1,ll), a );
  }

  for(int i=0; i<n ; i++){
    arma::colvec probs(J+1), probs2(J+1);
    int LAB =-1;
    probs[0] = log(1-rho) + R::dnorm(z[i],m0,sqrt(s0),1);

    for(int jj=1; jj<(J+1) ; jj++){
      probs[jj] = log(rho) + wd_mixmom_cppBNP(z[i],jj, MuSig2, k,a, 1, costs) + log( pis(jj-1) );
    }
  
    if(arma::is_finite(max(probs))){
      probs2 = exp( probs - max(probs));
      LAB    = RcppArmadillo::sample(possible_label0J, 1, TRUE, probs2)[0];
      }else{
      LAB    = RcppArmadillo::sample(possible_label0J, 1, TRUE)[0];  
    }
    //Rcout << 3;
    if(LAB==0){
      newLambdaGamma(i,0) = 0;
      newLambdaGamma(i,1) = 0;
    }else{ 
      newLambdaGamma(i,0) = 1;
      newLambdaGamma(i,1) = LAB;
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

/////////////////////////////////   //////////////////////////////////////////////////////////////////////////777

