//#########################################################
//#                                                       #
//#            Bliss method : rcpp code                   #
//#                                                       # 
//#########################################################
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <string>
#include <iostream>
#include <vector>
#include <cstring> 
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//#############################################################################
//############################ basic functions  ###############################
//#############################################################################

// square 
inline double sq(double x){
  return x * x;
}

// The R function : ginv (generalized matrix inversion using SVD decomposition)
// [[Rcpp::export]]
mat ginv_cpp (mat & X, double tol){
  int p;
  mat u;
  vec s;
  mat v;
  
  svd(u,s,v,X);
  p = s.size();
  
  tol   = tol * s(0);
  mat S = zeros<mat>(p,p);
  
  for( unsigned i=0 ; i<p; ++i){
    if( s(i) > tol ) S(i,i) = 1/s(i);
  }
  
  return( v * (S * trans(u)) );
}

// Termwise product 
vec termwise_product (vec & v, vec & u){
  vec res = v;  
  for(int i=0 ; i<v.size() ; ++i){
    res(i) = res(i) * u(i);
  }
  return(res);
}

// Weighted sample in 0:n-1.
int sample_weight (vec proba){
  if(sum(proba) == 0)   proba = ones<vec>(proba.n_rows) ;
  vec proba_cum = cumsum(proba)/sum(proba) ;
  
  int ret = 0;
  double u = R::runif(0,1);
  while (ret <= proba_cum.n_rows and u > proba_cum(ret)) {
    ret++;
  }
  return ret;
}

// Vector weighted sample in 0:n-1.
vec sample_weight (int n, vec proba){
  vec res = zeros<vec>(n);
  for (unsigned i = 0 ; i < n ; ++i) {
    res(i) = sample_weight(proba);
  }
  return res;
}

// Return the vector vec[-k].
vec vec_drop_k(vec & vecteur, int k){
  unsigned n = vecteur.n_rows;  
  if (k==0   && n > 1) return vecteur.subvec(1,n-1);
  if (k==n-1 && n > 1) return vecteur.subvec(0,n-2);
  vec res = zeros<vec>(n-1);
  
  res.subvec( 0  , k-1 ) = vecteur.subvec(0,k-1);
  res.subvec( k  , n-2 ) = vecteur.subvec(k+1,n-1);
  return res;
}

// Return the matrix mat[,-k].
mat mat_drop_col_k(mat & matrix, int k){
  unsigned n = matrix.n_rows;  
  unsigned p = matrix.n_cols;  
  if(k < 0  || k > (p-1)) return matrix;
  if(k==0   && p>1) return matrix.cols(1,p-1);
  if(k==p-1 && p>1) return matrix.cols(0,p-2);
  
  mat res = zeros<mat>(n,p-1);
  res.cols( 0 , k-1 ) = matrix.cols(0   , k-1);
  res.cols( k , p-2 ) = matrix.cols(k+1 , p-1);
  
  return res;
}

// seq.
vec sequence(int a,int b,double by){
  int range = floor((b-a)/by + 1) ;
  vec res = zeros<vec>(range);
  for(int i=0 ; i<range ; i++){
    res(i) = a + i*by;
  }
  return res;
}

// Compute the square root matrix using the SVD decomposition 
mat sqrt_mat (mat & X){
  int p;
  mat u;
  vec s;
  mat v;
  
  svd(u,s,v,X);
  p = s.size();
  
  mat S = zeros<mat>(p,p);
  
  for( unsigned i=0 ; i<p; ++i){
    S(i,i) = sqrt(s(i));
  }
  
  return( u * (S * trans(v)) );
}

// Simulate from a multidimensional gaussian.
vec mvrnormArma(vec mu, mat VarCovar, double sigma_sq) {
  int ncols = VarCovar.n_cols;
  vec Y = randn<vec>(ncols); 
  VarCovar = chol(VarCovar); 
  
  return  mu + sqrt(sigma_sq) * trans(trans(Y) * VarCovar);
}

// Compute a trapezoidal approximation of area under curve.
// [[Rcpp::export]]
double integrate_trapeze (vec & x, vec & y){
  vec diff_x = x.subvec(1,x.size()-1) - x.subvec(0,x.size()-2);
  vec cumu_y = y.subvec(1,y.size()-1) + y.subvec(0,y.size()-2);
  return sum( termwise_product (diff_x  , cumu_y ) )/2 ;
}

//  Compute the norm of a vector.
double norm_fct(vec & x,vec & y){
  vec tmp = zeros<vec>(x.size());
  for(int i=0 ; i<tmp.size() ; ++i){
    tmp(i) = sq( y(i) );
  }
  double res;
  res = sqrt(integrate_trapeze(x,tmp));
  
  return res;
}

// Find the first element of a vector which equals to n.
int which_first (NumericVector & v, int n){
  for(int i=0 ; i<v.size() ; ++i){
    if( v(i)==n ) return i;
  }
  return -1;
}

// Concatenate two vectors.
vec concatenate (vec & v, vec & u){
  vec res = zeros<vec>( v.size()+u.size() ) ;
  
  res.subvec( 0        ,   v.size()-1 ) = v;
  res.subvec( v.size() , res.size()-1 ) = u;
  return res;
}

// Compute a moving average on the vector v. 
vec moving_average_cpp (vec & v, int range){
  int n = v.size();
  vec res = zeros<vec>(n) ;
  int b_inf;
  int b_sup;
  
  for( int i=0; i<n; ++i){
    if(i-range < 0  ) b_inf = 0   ; else b_inf = i-range;
    if(i+range > n-1) b_sup = n-1 ; else b_sup = i+range;
    res(i) = mean( v.subvec(b_inf,b_sup) )  ;
  }
  
  return res;  
}

// Use to compute an uniform function, see function beta_build.
vec uniform_cpp (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  vec index = sequence(m-l,m+l,1);
  int tmp;
  for(int i=0 ; i<index.size() ; ++i){
    tmp = index(i);
    if( (tmp <= p) && (tmp >= 1) ){
      res(index(i)-1) = 1;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Beta_build in cpp.
// [[Rcpp::export]]
vec beta_build_cpp (vec & beta_star, vec & m, vec & l, vec & grid, 
                    int p, int K, std::string basis, mat & scale_ml ){
  vec res = zeros<vec>(p) ;
  
  if(basis == "uniform"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i)/scale_ml( m(i)-1 , l(i)-1 ) *
        uniform_cpp(m(i),l(i),grid);
    }    
  }
  return res;
}

// Compute the functions beta_i for each iteration i.
// [[Rcpp::export]]
mat compute_beta_functions_cpp (mat &  trace, int p, int K, vec & grid, 
                                std::string basis, mat & scale_ml){
  mat res = zeros<mat>(trace.n_rows,p);
  vec beta_star;
  vec m   ;
  vec l   ;
  vec tmp ;
  
  for(int i=0 ; i<res.n_rows ; ++i){
    tmp  = trans(trace.row(i))     ;
    beta_star = tmp.subvec(0,K-1)    ;
    m    = tmp.subvec(K,2*K-1)   ;
    l    = tmp.subvec(2*K,3*K-1) ;
    
    res.row(i) = trans(beta_build_cpp(beta_star,m,l,grid,p,K,basis,scale_ml)) ;
  }
  return res ;
}

// Extract an element from a cube.
double cube_extract(NumericVector & cube, int x , int y, int z, vec & dims){
  double res;
  res = cube[x + y*dims(0) + z*dims(1)*dims(0)];
  return res;
}

// Compute all the alternative for the value of the intergral for all m and l. 
arma::cube potential_intervals_List(List & X_list, List & grids,vec & l_max_vec, 
                                     CharacterVector & basis_vec, int q){ 
  mat X = as<mat>(X_list[q]);
  vec grid = as<vec>(grids[q]);
  int l_max = l_max_vec(q);
  std::string basis = as<std::string>(basis_vec(q));
  
  int n = X.n_rows ;
  int p = X.n_cols ;
  vec tub;
  
  arma::cube res(p,l_max,n+1);
  vec tmp;
  vec X_tmp;
  vec tmp2;
  if(basis == "uniform"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  
  // pour plus de clarete, je recommence une nouvelle double boucle plutot 
  // que de copier ce bout de code dans chacun des if. A changer ?
  for(int i=0 ; i<p ; ++i){
    for(int j=0 ; j<l_max ; ++j){
      // normalize by the scale \hat{s}_k
      tub = res.tube(i,j);
      tub = tub.subvec(0,n-1) ;
      res(i,j,n) = stddev( tub ); 
      for( int k=0 ; k<n ; ++k){
        res(i,j,k) = res(i,j,k) / res(i,j,n);
      }
    }
  }
  return res;
}

// Compute the matrix V (for a Ridge Zellner prior) 
// (for Q functional covaribles)                             
mat compute_Sigma_inv (int Q, vec & K, double g, mat & X_tilde, int sum_K,
                       mat & lambda_id0){                                       
  mat Sigma_inv = zeros<mat>(sum_K+1,sum_K+1);
  mat lambda_id = lambda_id0 ;
  
  mat X_tilde_temp = mat_drop_col_k(X_tilde,0);  
  mat u;
  vec s;
  mat v;
  svd(u,s,v,X_tilde_temp);
  
  Sigma_inv(0,0) = 1/lambda_id(0,0); 
  for( unsigned i=1 ; i<sum_K+1; ++i){
    lambda_id(i,i) = lambda_id(i,i) * max(s);
  }
  
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    Sigma_inv.submat(1+count,1+count,K(q)+count,K(q)+count) = 
      ( trans(X_tilde.cols(1+count,K(q)+count)) *
      X_tilde.cols(1+count,K(q)+count) +
      // ( trans( X_tilde.submat(0,1+count,n_y-1,K(q)+count)) *
      // X_tilde.submat(0,1+count,n_y-1,K(q)+count) +
      lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
    count = count + K(q);
  }
  
  return Sigma_inv;
}

// Compute SSE from the matrix W_sigma_inv
double compute_SSE(vec & y, mat & X_tilde, vec & beta_tilde, 
                   mat & mat_precision){
  double SSE = dot( y - X_tilde * beta_tilde, mat_precision * 
                    (y - X_tilde * beta_tilde));
  return(SSE);
}

// Compute the precision matrix W 
mat compute_W(List & weights, vec & n_e, unsigned E, int n){
  mat res = zeros<mat>(n,n);
  int count = 0;
  for( unsigned e=0 ; e<E ; ++e){
    vec weights_temp = weights[e];
    for( unsigned i=0 ; i<n_e(e) ; ++i){
      res( count + i , count + i ) = weights_temp(i) ;
    }
    count = count + n_e(e);
  }
  return(res);
}

// Extract a subvector from the cube all_intervals with a m_k and a l_k.
vec all_intervals_extract (NumericVector & all_intervals, int mk , int lk, 
                           vec & dims) {
  vec res = zeros<vec>(dims(2));
  for (int i = 0; i < dims(2); i++) {
    res(i) = cube_extract(all_intervals, mk - 1, lk - 1, i, dims);
  }
  return res;
}

// Extract a submatrix from the cube all_intervals with the vectors m and l.
mat extraire_X_tilde(vec & m, vec & l, NumericVector & all_intervals, 
                     vec & dims){
  int K = m.size();  
  mat res = ones<mat>(dims(2), K + 1);
  for(int i=0; i<K ; i++){
    res.col(i+1) = all_intervals_extract(all_intervals,m(i),l(i),dims);
  }  
  return res;
}  

// Compute the loss function for a proposal d
double loss_cpp (vec & d, vec & grid, vec & posterior_expe){
  vec tmp  = d-posterior_expe ;
  
  return sq(norm_fct(grid, tmp) ); 
}

// Compute the decrease of the Temperature
double cooling_cpp (int i, double Temp){  
  double res;
  res = Temp / log( ( i / 10)*10 + exp(1));
  return res;
}
//#############################################################################
//############################ Auxiliary functions ############################
//#############################################################################
// Update the parameter m_k
int mk_update_List (int count, int k, vec & y_bar, vec & beta_tilde, 
                    double sigma_sq,vec & m_q, vec & l_q, mat & X_tilde_bar, 
                    mat & X_tilde, NumericVector & all_intervals_q, 
                    vec & all_intervals_dims_q, vec & m_alternative_q, 
                    vec & phi_m_q, int p_q, int Q,vec & K, double g, int sum_K, 
                    mat & lambda_id0, mat & W) {                                    
  double aux;
  vec aux2 = zeros<vec>(p_q);
  vec aux3 = zeros<vec>(p_q);
  mat X_tilde_bar_temp = X_tilde_bar;
  mat Sigma_inv_temp;
  
  // Compute the probabilities
  vec probs = ones<vec>(p_q);
  vec X_tilde_bar_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    X_tilde_bar_qki = all_intervals_extract(all_intervals_q,m_alternative_q(i),
                                            l_q(k),all_intervals_dims_q);
    
    X_tilde_bar_temp.col(count + k + 1) = X_tilde_bar_qki;
    Sigma_inv_temp = compute_Sigma_inv (Q,K, g, X_tilde,sum_K,lambda_id0) ;  
      
    aux2(i) = compute_SSE(y_bar,X_tilde_bar_temp, beta_tilde, W) ;
    aux2(i) = aux2(i) + dot( beta_tilde , Sigma_inv_temp * beta_tilde );
    aux2(i) = aux2(i) / (2*sigma_sq) ;
    aux3(i) = sqrt( det(Sigma_inv_temp) );
  }  
  
  double min_aux = min(aux2);
  for(int  i=0 ; i<p_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_m_q(i);
  }
  // Simulate a mk
  int mk = sample_weight(probs) + 1 ;
  
  return mk;
}

// Update the parameter l_k
int lk_update_List (int count, int k, vec & y_bar, vec & beta_tilde,
                    double sigma_sq,vec & m_q, vec & l_q, mat & X_tilde_bar, 
                    mat & X_tilde, NumericVector & all_intervals_q, 
                    vec & all_intervals_dims_q,vec & l_alternative_q, 
                    vec & phi_l_q, int lmax_q, int Q,vec K, double g, int sum_K,
                    mat & lambda_id0, mat & W) {
  double aux;
  vec aux2 = zeros<vec>(lmax_q);
  vec aux3 = zeros<vec>(lmax_q);
  mat X_tilde_bar_temp = X_tilde_bar;
  mat Sigma_inv_temp;
  
  // Compute the probabilities
  vec probs = ones<vec>(lmax_q);
  vec X_tilde_bar_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<lmax_q ; ++i){
    X_tilde_bar_qki = all_intervals_extract(all_intervals_q,m_q(k),
                                        l_alternative_q(i),
                                        all_intervals_dims_q);
    
    X_tilde_bar_temp.col(count + k + 1) = X_tilde_bar_qki;
    Sigma_inv_temp = compute_Sigma_inv (Q,K, g, X_tilde,sum_K,lambda_id0) ;  
    
    aux2(i) = compute_SSE(y_bar,X_tilde_bar_temp, beta_tilde, W) ;
    aux2(i) = aux2(i) + dot( beta_tilde , Sigma_inv_temp * beta_tilde );
    aux2(i) = aux2(i) / (2*sigma_sq) ;
    aux3(i) = sqrt( det(Sigma_inv_temp) );
  }  
  
  double min_aux = min(aux2);
  for(int  i=0 ; i<lmax_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  int lk = sample_weight(probs) + 1 ;
  
  return lk;
}

// update the parameter sigma_sq
double sigma_sq_update_List (vec & y_bar, vec & beta_tilde, mat & W, 
                             mat & Sigma_beta_tilde_inv, mat & X_tilde_bar,
                             double a_star, double b) {                                    
  double b_star = b + 0.5 * compute_SSE(y_bar,X_tilde_bar,beta_tilde,W) +
    0.5 * dot(beta_tilde, Sigma_beta_tilde_inv * beta_tilde);
  double res = 1. / (R::rgamma(a_star, 1/b_star) );
  
  return res ;
}

// update the parameter beta_star
vec beta_tilde_update (mat & Sigma_beta_tilde_inv, vec & mu_beta_tilde, 
                       double sigma, double tol) {                  
  return mvrnormArma( ginv_cpp(Sigma_beta_tilde_inv,tol) * mu_beta_tilde , 
                      ginv_cpp(Sigma_beta_tilde_inv,tol), sigma); 
}

// function to update the matrix X_tilde with the new intervals
mat update_X_tilde (int Q, vec & K, List & all_intervals, 
                    List & all_intervals_dims, List & m, List & l, 
                    mat & X_tilde){
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    for(unsigned k=0 ; k<K(q) ; ++k) {    
      vec m_temp = m[q] ;
      vec l_temp = l[q] ;
      NumericVector all_intervals_temp = all_intervals[q];
      vec all_intervals_dims_temp = all_intervals_dims[q];
      
      X_tilde.col(k+1+count) = all_intervals_extract(all_intervals_temp,
                  m_temp(k),l_temp(k),all_intervals_dims_temp);
    }
    count = count + K(q);
  }
  
  return X_tilde;
}
//##############################################################################
//#################### Gibbs Sampler and Simulated Annealing ###################
//##############################################################################

// Perform the Gibbs Sampler algorithm for the Bliss model
// returned values; trace and param.
// trace : a matrix, with the different parameters in columns and the 
// iterations in rows. 
// The different parameters are : beta_star_1, m1, l1, beta_star_2, m2, 
//l2, ..., beta_star_Q, mQ, lQ, mu, sigma2.
// Hence the trace has (iter + 1) rows (+1 for the initailisation), and 
// 3*(K1+K2+...KQ)+2 columns. Indeed, 
// beta_star_1, m1 and l1 are of length K, ..., beta_star_Q, mQ and lQ are of 
// length KQ.
// [[Rcpp::export]]
List Bliss_multiple_WL_cpp (int Q, int iter, List & grids, bool posterior, 
                            
                            vec & y_obs, List & X_obs, vec & y_expert, 
                            List & X_expert, vec & n_e, 
                            List & weights, vec & average_weights,
                            
                            vec & K, vec & lmax, List & probs_m, 
                            List & probs_l, double lambda, mat & V_tilde, 
                            double tol, CharacterVector & basis) {                  
  std::cout << "Gibbs Sampler: " <<  std::endl;
  std::cout << "\t Initialization." <<  std::endl;
  
  // Compute the value of n and the p's                                         
  int n_y = as<mat>(X_obs[0]).n_rows ;
  int n = sum(n_e)  ;
  int E = n_e.size();
  
  vec p = zeros<vec>(Q) ;
  for(int i=0 ; i<Q ; ++i){
    p(i) = as<mat>(X_obs[i]).n_cols; 
  }
  
  // if prior == false : merge observed data and pseudo data
  vec y_bar;
  List X_bar(Q);
  if(posterior == false){
    y_bar = y_expert;
    X_bar = X_expert;
  }else{
    y_bar = concatenate(y_obs,y_expert);
    for(unsigned q=0 ; q<Q ; ++q){
      mat temp          = zeros<mat>(n,p(q));
      mat X_obs_temp    = X_obs[q];
      mat X_expert_temp = X_expert[q];
      
      temp.rows(0,n_y-1) = X_obs_temp;
      temp.rows(n_y,n-1) = X_expert_temp;
      X_bar[q] = temp;
    }
  }
  y_bar = y_bar - mean(y_bar) ;                                                 // centrage a faire ici ?
  
  // Compute projection of the x_i on all the intervals
  std::cout << "\t Compute the decompositions (expert)." <<  std::endl;
  List scale_ml(Q);                                                             // scale_ml is used to normalize the predictors
  List all_intervals(Q);                                                        // will be contain all the projections
  List all_intervals_dims(Q);                                                   // will be contain the dim of the all_intervals's
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (X_bar, grids, lmax, basis,q); 
    scale_ml[q]      = temp.slice(n);
    
    temp = temp.subcube(0,0,0,p(q)-1,lmax(q)-1,n-1);
    all_intervals[q] = temp;
    
    vec temp2 = zeros<vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = lmax(q) ;
    temp2(2) = n    ;
    all_intervals_dims[q] = temp2;
  }
  
  // Compute projection of the x_i on all the intervals (for the observed x)
  std::cout << "\t Compute the decompositions (data)." <<  std::endl;
  List scale_ml_obs(Q);                                                         // scale_ml is used to normalize the predictors
  List all_intervals_obs(Q);                                                    // will be contain all the projections
  List all_intervals_dims_obs(Q);                                               // will be contain the dim of the all_intervals's
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (X_obs, grids, lmax, basis,q); 
    scale_ml[q]      = temp.slice(n_y);
    
    temp = temp.subcube(0,0,0,p(q)-1,lmax(q)-1,n_y-1);
    all_intervals_obs[q] = temp;
    
    vec temp2 = zeros<vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = lmax(q) ;
    temp2(2) = n_y     ;
    all_intervals_dims_obs[q] = temp2;
  }
  
  // Compute the matrix of lambda for the Ridge penalty of the Rigde Zellner prior...
  int sum_K = sum(K);
  mat lambda_id0  = zeros<mat>(sum_K+1,sum_K+1) ;                               
  lambda_id0(0,0) = 100*var(y_obs);               
  for( unsigned i=1 ; i<sum_K+1; ++i){                            
    lambda_id0(i,i) = lambda ;
  }
  
  // Determine the start point
  std::cout << "\t Determine the start point." <<  std::endl;
  double a = 0.1;
  double b = 0.1;
  double g = y_obs.size();
  
  double sigma_sq           ;
  vec beta_tilde            ;
  mat Sigma_beta_tilde_inv  ;
  
  double a_star = a + 0.5*( dot(n_e,average_weights) + sum_K + 1) ;                            
  
  vec mu_beta_tilde ;
  mat Sigma_inv     ;
  
  bool success = false ;
  mat R                ;
  mat test             ;
  
  List m(Q) ;
  List l(Q) ;
  
  mat X_tilde_bar = ones<mat>(n  ,sum_K+1) ;  
  mat X_tilde     = ones<mat>(n_y,sum_K+1) ; 
  mat W = compute_W(weights,n_e,E,n)       ;
  
  // Try to determine a starting point which not leads to a non-invertible matrix problem
  while(success == false){ 
    // Initialization of sigma_sq
    sigma_sq  = var(y_bar) ;
    
    // Initialization of the middle and length of the intervals
    for( unsigned q=0 ; q<Q ; ++q){  
      vec probs_m_temp = probs_m[q];
      m[q]         = sample_weight(K(q),probs_m_temp) + 1 ; 
      vec probs_l_temp = probs_l[q];
      l[q]         = sample_weight(K(q),probs_l_temp) + 1 ;
    }
    
    // Initialize X_tilde_bar and X_tilde
    int count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      for(unsigned k=0 ; k<K(q) ; ++k) {    
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        NumericVector all_intervals_temp     = all_intervals[q];
        NumericVector all_intervals_temp_obs = all_intervals_obs[q];
        vec all_intervals_dims_temp     = all_intervals_dims[q];
        vec all_intervals_dims_temp_obs = all_intervals_dims_obs[q];
        
        X_tilde_bar.col(k+1+count) = all_intervals_extract(all_intervals_temp,
                        m_temp(k),l_temp(k),all_intervals_dims_temp);
        
        X_tilde.col(k+1+count) = all_intervals_extract(all_intervals_temp_obs,
                    m_temp(k),l_temp(k),all_intervals_dims_temp_obs);
      }
      count = count + K(q);
    }
    
    // Initialize the current Sigma_inv matrix (which depend on the intervals)
    Sigma_inv = compute_Sigma_inv (Q,K, g, X_tilde,sum_K,lambda_id0) ;  
    
    // Check if there is a non-invertible matrix problem
    Sigma_beta_tilde_inv = Sigma_inv + trans(X_tilde_bar) * W * X_tilde_bar;
    test            = ginv_cpp(Sigma_beta_tilde_inv,tol)    ; 
    success         = accu(abs(test)) != 0                  ;
  }
  // Initialization of beta_tilde
  beta_tilde = mvrnormArma( zeros<vec>(sum_K+1) , ginv_cpp(Sigma_inv,tol) ,
                            sigma_sq) ; 
  
  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*sum_K+2);  
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    vec m_temp = m[q] ;
    vec l_temp = l[q] ;
    trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
      trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
    trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
    trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;
    
    trace(0,3*sum_K  ) = beta_tilde(0) ;  
    trace(0,3*sum_K+1) = sigma_sq ;  
    count = count + K(q) ;
  }
  
  // Initialize some variable used in the Gibbs loop
  List l_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    l_alternative[q] = sequence(1,lmax(q),1);
  }
  List m_alternative(Q);
  for( int q=0 ; q<Q ; ++q){
    m_alternative[q] = sequence(1,p(q),1);
  }
  
  // The MWG loop
  std::cout << "\t Start the Gibbs loop." <<  std::endl;
  for(unsigned i=1  ; i < iter+1 ; ++i ) {      
    // Progress
    // std::cout << "\t " << i << std::endl;
    if( i % (iter / 10)  == 0) 
      std::cout << "\t " << i / (iter / 100) << "%" << std::endl;
    
    // update sigma_sq  
    sigma_sq = sigma_sq_update_List(y_bar,beta_tilde,W,Sigma_beta_tilde_inv,
                                    X_tilde_bar,a_star,b);
    
    // update m  
    count = 0 ; 
    // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( unsigned q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k 
      vec m_q = m[q];
      vec l_q = l[q];
      int p_q = p(q);
      NumericVector all_intervals_q = all_intervals[q];
      vec all_intervals_dims_q      = all_intervals_dims[q];
      vec m_alternative_q = sequence(1,p_q,1) ; 
      vec probs_m_q       = probs_m[q];
      
      for(int k=0 ; k<K(q) ; ++k){
        // update m_k
        m_q(k) = mk_update_List(count,k,y_bar,beta_tilde,sigma_sq,m_q,l_q,      
            X_tilde_bar,X_tilde,all_intervals_q,all_intervals_dims_q,
            m_alternative_q, probs_m_q,p_q, Q,K,g,sum_K,lambda_id0,  W);
      }
      
      // update the value "X_tilde_bar" and "X_tilde"
      X_tilde_bar = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,   
                                   X_tilde_bar);
      X_tilde     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                   m,l,X_tilde);
      
      // Update the m_q value
      m[q] = m_q;
      // Update count 
      count = count + K(q);
    }
    
    // update l
    count = 0 ; 
    // count is used to browse some vec/mat when p(q) is not constant wrt q.
    for( unsigned q=0 ; q<Q ; ++q ){
      // Compute some quantities which do not vary with k 
      vec m_q = m[q];
      vec l_q = l[q];
      int lmax_q = lmax(q);
      NumericVector all_intervals_q = all_intervals[q];
      vec all_intervals_dims_q      = all_intervals_dims[q];
      vec l_alternative_q = sequence(1,lmax_q,1) ; 
      vec probs_l_q         = probs_l[q];
      
      for(int k=0 ; k<K(q) ; ++k){
        // update l_k
        l_q(k) = lk_update_List(count,k,y_bar,beta_tilde,sigma_sq,m_q,l_q,      
            X_tilde_bar,X_tilde,all_intervals_q,all_intervals_dims_q,
            l_alternative_q, probs_l_q,lmax_q, Q,K,g,sum_K,lambda_id0,W); 
        
        // update the value "X_tilde_bar" and "X_tilde"
        X_tilde_bar = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,   
                                     X_tilde_bar);
        X_tilde     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                     m,l,X_tilde);
      }
      
      // Update the l_q value
      l[q] = l_q;
      // Update count 
      count = count + K(q);
    }
    
    // update the value "Sigma_inv" (only after the updating of the m's and l's)
    Sigma_inv = compute_Sigma_inv (Q,K, g, X_tilde,sum_K,lambda_id0) ;  
    
    // update the matrix Sigma_beta_tilde (only after the updating of 
    // the m's and l's)   
    Sigma_beta_tilde_inv = Sigma_inv + trans(X_tilde_bar) * W * X_tilde_bar;
    test            = ginv_cpp(Sigma_beta_tilde_inv,tol)    ; 
    success         = accu(abs(test)) != 0                  ;
    if(success){
      // update the beta_tilde   
      mu_beta_tilde = trans(X_tilde_bar) * W * y_bar;
      beta_tilde = beta_tilde_update(Sigma_beta_tilde_inv,mu_beta_tilde ,sigma_sq,
                                     tol) ; 
      
      // update the matrix trace
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        vec m_temp = m[q] ;
        vec l_temp = l[q] ;
        trace.row(i).subvec( 3*count        , 3*count+  K(q)-1) = 
          trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
        trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1) = trans(m_temp);
        trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1) = trans(l_temp);
        
        trace(i,3*sum_K  ) = beta_tilde(0) ;  
        trace(i,3*sum_K+1) = sigma_sq ;  
        count = count + K(q) ;
      }
    }else{ //... go back to the beginning of the updating process.
      i     = i - 1 ;
      count = 0;
      for( unsigned q=0 ; q<Q ; ++q){
        beta_tilde.subvec( 1+count , K(q)+count ) =
          trans(trace.row(i).subvec( 3*count , 3*count+  K(q)-1))  ;
        m[q] = trans(trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)) ;
        l[q] = trans(trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)) ;
        
        beta_tilde(0) = trace(i,3*sum_K  ) ;  
        sigma_sq     = trace(i,3*sum_K+1) ;  
        count = count + K(q) ;
      }
      
      // update the value "X_tilde_bar" and "X_tilde"
      X_tilde_bar = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,   
                                   X_tilde_bar);
      X_tilde     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                   m,l,X_tilde);
      
      // update the value "Sigma_inv" 
      Sigma_inv = compute_Sigma_inv (Q,K, g, X_tilde,sum_K,lambda_id0) ;  
    }
  }
  
  // return the trace and the parameters
  std::cout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]=trace,
                       _["param"]=List::create(_["all_intervals"]=all_intervals,
                                          _["scale_ml"]=scale_ml));
}

// Perform the Simulated Annealing algorithm to minimize the loss function
// [[Rcpp::export]]
List Bliss_Simulated_Annealing_cpp (int iter, mat beta_functions, vec & grid, 
                                    int burnin, double Temp,int k_max, 
                                    int l_max, int dm, int dl, 
                                    int p,std::string basis, mat & scale_ml){
  std::cout << "Simulated Annealing:" <<  std::endl;
  // Initialization
  std::cout << "\t Initialization." <<  std::endl;
  int N = beta_functions.n_rows;
  vec posterior_expe = zeros<vec>(p); 
  vec posterior_var  = zeros<vec>(p);
  for(int i=0 ; i<p ; ++i){
    posterior_expe(i) = mean(beta_functions.col(i));
    posterior_var(i)  =  var(beta_functions.col(i));
  }
  
  vec probs;
  int k;
  vec m;
  vec l;
  vec beta_star;
  vec d;
  double d_loss;
  int boundary_min;
  int boundary_max;
  vec difference;
  
  vec d_tmp;
  double d_loss_tmp;
  double proba_acceptance;
  double u;
  int j;
  double Temperature;
  vec beta_star_tmp;
  vec m_tmp;
  vec l_tmp;
  int k_tmp;
  int accepted;
  vec choice_prob_interval; 
  vec choice_prob;
  int choice;
  int interval_min;
  int interval_max;
  double var_beta_star;
  vec boundaries_min;
  vec boundaries_max;
  vec boundaries;
  int new_m;
  int new_l;
  double new_beta_star;
  
  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*k_max+3);
  
  // Determine the start point
  std::cout << "\t Determine the start point." <<  std::endl;
  probs = ones<vec>(k_max);
  k      = sample_weight( probs )+1;
  m      = zeros<vec>(k); 
  l      = zeros<vec>(k);
  beta_star   = zeros<vec>(k);
  
  probs = ones<vec>(p);
  m(0)   = sample_weight( probs )+1;
  probs = ones<vec>(l_max);
  l(0)   = sample_weight( probs )+1;  
  
  boundary_min = m(0)-l(0)-1;
  boundary_max = m(0)+l(0)-1;
  if(boundary_min < 0   ) boundary_min = 0   ;
  if(boundary_max > p-1 ) boundary_max = p-1 ;  
  beta_star(0) = mean(posterior_expe.subvec( boundary_min , boundary_max ));
  d = beta_build_cpp(beta_star,m,l,grid,p,1,basis,scale_ml);
  
  if(k > 1){
    for(int i=1 ; i<k ; ++i ){
      // Compute the difference ...
      difference = abs(posterior_expe - d);      
      // ... and its smoothed version. 
      difference = moving_average_cpp(difference,4);
      
      // Which intervals are possible ?
      for(int o=0 ; o<i ; ++o){
        if( m(o) - l(o) -1 > 0  ) boundary_min = m(o) - l(o) -1; else boundary_min = 1;
        if( m(o) + l(o) +1 < p+1) boundary_max = m(o) + l(o) +1; else boundary_max = p; 
        for(int z=boundary_min; z < boundary_max+1 ; ++z){
          difference(z-1) = 0;  
        }
      }
      
      // Simulate an interval
      if(sum(difference) > 0){
        vec boundaries;
        vec boundaries_min;
        vec boundaries_max;
        int boundary_min;
        int boundary_max;
        
        // Simulate a m
        m(i) = sample_weight(difference) +1;
        
        // Simulate a l
        boundaries_max = zeros<vec>(i);
        boundaries_min = zeros<vec>(i);
        
        boundaries_min = abs(m.subvec(0,i-1) - l.subvec(0,i-1) - m(i))-1 ;
        boundaries_max = abs(m.subvec(0,i-1) + l.subvec(0,i-1) - m(i))-1 ;
        
        boundaries = zeros<vec>(2*i+1);
        boundaries(0) = l_max;
        boundaries.subvec(1  ,i  ) = boundaries_min;
        boundaries.subvec(1+i,i+i) = boundaries_max;
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);
        
        if(boundary_max < 1) boundary_max = 1;      
        
        l(i) = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        // Simulate a beta_star (from the smoothed difference)
        if( m(i) - l(i) -1 > 0  ) boundary_min = m(i) - l(i) -1; else boundary_min = 1;
        if( m(i) + l(i) +1 < p+1) boundary_max = m(i) + l(i) +1; else boundary_max = p; 
        beta_star(i) = mean( difference.subvec(boundary_min-1 , boundary_max-1) );
        // Compute the function with these intervals
        d = beta_build_cpp(beta_star,m,l,grid,p,i+1,basis,scale_ml);
      }else{ 
        // sortir de la boucle avec une bonne valeur de k
        unsigned i_tmp = i; 
        i = k ;
        k = i_tmp; 
      }
    }      
  }    
  
  // Compute the first function with K intervals (and its loss)
  d      = beta_build_cpp(beta_star,m,l,grid,p,k,basis,scale_ml);
  d_loss = loss_cpp(d,grid,posterior_expe);    
  
  // Update the trace with the start point
  trace.row(0).subvec( 0       ,           k-1) = trans(beta_star.subvec(0,k-1)) ;
  trace.row(0).subvec( k_max   , k_max   + k-1) = trans(m.subvec(0,k-1))         ;
  trace.row(0).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec(0,k-1))         ;
  trace(0,3*k_max)   = 1      ;
  trace(0,3*k_max+1) = k      ;
  trace(0,3*k_max+2) = d_loss ;      
  
  // Start the loop
  std::cout << "\t Start the loop." <<  std::endl;      
  for(int i=0 ; i<iter ; ++i){
    Temperature = cooling_cpp(i,Temp);
    // Progress
    if( (i+1) % (iter / 10)  == 0) 
      Rcpp::Rcout << "\t " << (i+1) / (iter / 100) << "%" << std::endl;
    // Initialize the proposal
    beta_star_tmp = beta_star;
    m_tmp    = m   ;
    l_tmp    = l   ;
    k_tmp    = k   ;
    accepted  = 0   ;
    choice_prob_interval = ones<vec>(k);    
    
    // Choose a move 
    choice_prob = ones<vec>(5);
    if(k == k_max) choice_prob(3) = choice_prob(3) -1;
    if(k == 1    ) choice_prob(4) = choice_prob(4) -1;
    
    choice = sample_weight( choice_prob ) + 1;
    // change a beta_star
    if(choice == 1){        
      // choose an interval
      j = sample_weight(choice_prob_interval);
      
      // Simulate a new beta_star_k
      interval_min = m(j)-l(j) -1;
      interval_max = m(j)+l(j) -1;
      if(interval_min < 0   ) interval_min = 0   ;
      if(interval_max > p-1 ) interval_max = p-1 ;
      
      var_beta_star = mean(posterior_var.subvec( interval_min , interval_max));    
      beta_star_tmp(j) = R::rnorm( beta_star(j), sqrt(var_beta_star) );
    }
    // change a m_k
    if(choice == 2){      
      // choose an interval
      j = sample_weight(choice_prob_interval);
      
      // Simulate a new m_k
      if(k > 1){
        boundaries_max = zeros<vec>(k);
        boundaries_max.subvec(0,k-2) = vec_drop_k(m,j).subvec(0,k-2) - 
          vec_drop_k(l,j).subvec(0,k-2) - l(j)-1 ;
        boundaries_max(k-1)          = p ;
        boundaries_max               = sort(boundaries_max);
        boundary_max                 = boundaries_max(boundaries_max.size()-1);
        
        for(int o=0 ; o<boundaries_max.size() ; ++o){
          if(m(j) < boundaries_max(o) ){
            boundary_max = boundaries_max(o) ;
            break;
          }
        }
        
        boundaries_min = zeros<vec>(k);
        boundaries_min.subvec(0,k-2) = vec_drop_k(m,j).subvec(0,k-2) + 
          vec_drop_k(l,j).subvec(0,k-2) + l(j)+1 ;
        boundaries_min(k-1)          = 1 ;
        boundaries_min               = sort(boundaries_min);
        boundary_min                 = boundaries_min(0);
        
        for(int o=boundaries_max.size()-1 ; o>=0 ; --o){
          if(m(j) > boundaries_min(o) ){
            boundary_min = boundaries_min(o) ;
            break;
          }
        }        
      }else{
        boundary_max = m(j) + dm ;
        boundary_min = m(j) - dm ;
        if(boundary_max > p) boundary_max = p ;
        if(boundary_min < 1) boundary_min = 1 ;
      }
      
      if(boundary_max - boundary_min + 1 > 0){
        probs = ones<vec>( boundary_max - boundary_min + 1 ) ;
        m_tmp(j)  = sample_weight( probs ) + boundary_min ;     
      }
    }
    // change a l_k
    if(choice == 3){
      // choose an interval
      j = sample_weight(choice_prob_interval);
      
      // Simulate a new l_k
      if(k > 1){
        boundaries_max = zeros<vec>(k-1);
        boundaries_min = zeros<vec>(k-1);
        boundaries_max = abs(vec_drop_k(m,j) + vec_drop_k(l,j) - m(j))-1 ;
        boundaries_min = abs(vec_drop_k(m,j) - vec_drop_k(l,j) - m(j))-1 ;
        
        boundaries = zeros<vec>(2*k-1);
        boundaries(0) = dl;
        boundaries.subvec(1,k-1)   = boundaries_min.subvec(0,k-2);
        boundaries.subvec(k,2*(k-1)) = boundaries_max.subvec(0,k-2);
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);
        
        if(boundary_max > 1){
          l_tmp(j) = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        } 
      }
    }
    // birth
    if(choice == 4){    
      // compute the difference ...
      difference = posterior_expe - d;      
      // ... and its smoothed version 
      difference = moving_average_cpp(difference,4);
      
      // Which intervals are possible ?
      for(int o=0 ; o<k ; ++o){
        if( m(o) - l(o) -1 > 0  ) boundary_min = m(o) - l(o) -1; else
          boundary_min = 1;
        if( m(o) + l(o) +1 < p+1) boundary_max = m(o) + l(o) +1; else 
          boundary_max = p; 
        for(int z=boundary_min; z < boundary_max+1 ; ++z){
          difference(z-1) = 0;  
        }
      }
      if(sum(abs(difference)) > 0){
        // update k
        k_tmp = k+1;
        // Simulate a new m
        new_m = sample_weight(abs(difference)) +1;
        m_tmp = zeros<vec>(k_tmp);
        m_tmp.subvec(0,k-1) = m;
        m_tmp(k_tmp-1)      = new_m;
        // Simulate a new l 
        boundaries_max = zeros<vec>(k_tmp-1);
        boundaries_min = zeros<vec>(k_tmp-1);
        
        boundaries_min = abs(m - l - new_m)-1 ;
        boundaries_max = abs(m + l - new_m)+1 ;
        
        boundaries = zeros<vec>(2*k+1);
        boundaries(0) = l_max;
        boundaries.subvec(1  ,k  ) = boundaries_min;
        boundaries.subvec(1+k,k+k) = boundaries_max;
        boundaries = sort(boundaries);
        boundary_max = boundaries(0);
        
        new_l = sample_weight( ones<vec>(boundary_max) ) + 1 ;
        l_tmp = zeros<vec>(k_tmp);
        l_tmp.subvec(0,k-1) = l;
        l_tmp(k_tmp-1)      = new_l;
        
        // Simulate a new beta_star (from the smoothed difference)        
        if( new_m - new_l -1 > 0  ) boundary_min = new_m - new_l -1; else boundary_min = 1;
        if( new_m + new_l +1 < p+1) boundary_max = new_m + new_l +1; else boundary_max = p; 
        new_beta_star = mean( difference.subvec(boundary_min-1 , boundary_max-1) );
        beta_star_tmp = zeros<vec>(k_tmp);
        beta_star_tmp.subvec(0,k_tmp-2) = beta_star;
        beta_star_tmp(k_tmp-1)          = new_beta_star;
      } 
    }       
    // death
    if(choice == 5){
      // Choose an interval to drop
      j = sample_weight(choice_prob_interval);
      
      // Drop the interval
      k_tmp = k-1;
      beta_star_tmp = zeros<vec>(k_tmp);
      m_tmp    = zeros<vec>(k_tmp);
      l_tmp    = zeros<vec>(k_tmp);
      
      beta_star_tmp = vec_drop_k(beta_star,j);
      m_tmp    = vec_drop_k(m   ,j);
      l_tmp    = vec_drop_k(l   ,j); 
    }      
    
    // Compute the acceptance probability
    d_tmp      = beta_build_cpp(beta_star_tmp,m_tmp,l_tmp,grid,p,k_tmp,basis,
                                scale_ml);
    d_loss_tmp = loss_cpp(d_tmp,grid,posterior_expe);
    
    proba_acceptance = exp( -( d_loss_tmp-d_loss )/ Temperature );      

    // Accept/reject
    u = R::runif(0,1) ;      
    if(u < proba_acceptance){
      beta_star = zeros<vec>(k_tmp)     ;
      l = zeros<vec>(k_tmp)     ;
      m = zeros<vec>(k_tmp)     ;
      accepted  = 1             ;
      beta_star = beta_star_tmp.subvec(0,k_tmp-1) ;
      m         = m_tmp.subvec(0,k_tmp-1)         ;
      l         = l_tmp.subvec(0,k_tmp-1)         ;
      k         = k_tmp         ;
      d         = d_tmp         ;
      d_loss    = d_loss_tmp    ;
    }
    
    // Update the trace
    trace.row(i+1).subvec( 0       ,           k-1) = trans(beta_star.subvec( 0,k-1)) ;
    trace.row(i+1).subvec( k_max   , k_max   + k-1) = trans(m.subvec( 0,k-1))         ;
    trace.row(i+1).subvec( k_max*2 , k_max*2 + k-1) = trans(l.subvec( 0,k-1))         ;
    trace(i+1,3*k_max  ) = accepted ;
    trace(i+1,3*k_max+1) = k        ; 
    trace(i+1,3*k_max+2) = d_loss   ;      
  }
  
  // Return the result
  std::cout << "\t Return the result." <<  std::endl;  
  return  List::create(_["trace"]         =trace,
                       _["posterior_expe"]=posterior_expe,
                       _["posterior_var"] =posterior_var);
}
