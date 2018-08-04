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
#include <Rcpp.h>
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

// Return the matrix mat[,-k].
mat mat_drop_row_k(mat & matrix, int k){
  unsigned n = matrix.n_rows;  
  unsigned p = matrix.n_cols;  
  if(k < 0  || k > (n-1)) return matrix;
  if(k==0   && n>1) return matrix.rows(1,n-1);
  if(k==p-1 && n>1) return matrix.rows(0,n-2);
  
  mat res = zeros<mat>(n-1,p);
  res.rows( 0 , k-1 ) = matrix.rows(0   , k-1);
  res.rows( k , n-2 ) = matrix.rows(k+1 , n-1);
  
  return res;
}

// useless
// // Return a sub vector.
// vec sub_vector(vec & v, vec & index){
//   vec res = zeros<vec>(index.size());
//   
//   for( unsigned i=0 ; i < index.size() ; ++i){
//     res(i) = v( index(i) );
//   }
//   
//   return(res);
// }

// seq.
vec sequence(int a,int b,double by){
  int range = floor((b-a)/by + 1) ;
  vec res = zeros<vec>(range);
  for(int i=0 ; i<range ; i++){
    res(i) = a + i*by;
  }
  return res;
}

// seq.
vec sequence2(double a, double b, unsigned le){
  double by = (b-a)/ le ;
  vec res = zeros<vec>(le);
  for(int i=0 ; i<le ; i++){
    res(i) = a + i*by;
  }
  return res;
}

// Extract an element from a cube.
double cube_extract(NumericVector & cube, int x , int y, int z, vec & dims){
  double res;
  res = cube[x + y*dims(0) + z*dims(1)*dims(0)];
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
  VarCovar = chol(VarCovar); // XXXXXXXXXXXX
  
  return  mu + sqrt(sigma_sq) * trans(trans(Y) * VarCovar);
}

// Compute a trapezoidal approximation of area under curve.
// [[Rcpp::export]]
double integrate_trapeze (vec & x, vec & y){
  vec diff_x = vec_drop_k(x,0) - vec_drop_k(x,x.size()-1);
  vec cumu_y = vec_drop_k(y,0) + vec_drop_k(y,y.size()-1);
  return sum( termwise_product (diff_x  , cumu_y ) )/2 ;
}

//  Compute the norm of the function y on the evaluation grid x.
double norm_fct(vec & x,vec & y){
  vec tmp = zeros<vec>(x.size());
  for(int i=0 ; i<tmp.size() ; ++i){
    tmp(i) = sq( y(i) );
  }
  double res;
  res = sqrt(integrate_trapeze(x,tmp));
  
  return res;
}

//  Compute the norm of the function y on the evaluation grid x weighted by g
double weighted_norm_fct(vec & x,vec & y, vec & g){
  vec tmp = zeros<vec>(x.size());
  for(int i=0 ; i<tmp.size() ; ++i){
    tmp(i) = sq( y(i) ) * g(i);
  }
  double res;
  // res = sqrt(integrate_trapeze(x,tmp));
  res = integrate_trapeze(x,tmp);
  
  return res;
}

// Use to compute an uniform function, see function beta_build.
// [[Rcpp::export]]
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

// Use to compute an uniform function, see function beta_build.
vec uniform_cpp_unnormalized (int m, int l, vec & grid){
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
  return res;
}

// Use to compute a triangular function, see function beta_build.
vec triangular_cpp (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1 - i/l_double ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a triangular function, see function beta_build.
vec triangular_cpp_unnormalized (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1 - i/l_double ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  return res;
}

// Use to compute a gaussian function, see function beta_build.
vec gaussian_cpp (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  exp( - 9*sq(i/l_double)/2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function beta_build.
vec gaussian_cpp_unnormalized (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  exp( - 9*sq(i/l_double)/2) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  return res;
}

// Use to compute a gaussian function, see function beta_build.
vec Epanechnikov_cpp (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1-sq(i/l_double) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
  double res_norm = norm_fct(grid,res);
  res = res / res_norm ;
  return res;
}

// Use to compute a gaussian function, see function beta_build.
vec Epanechnikov_cpp_unnormalized (int m, int l, vec & grid){
  int p = grid.size();
  vec res = zeros<vec>(p);
  int tmp;
  double l_double = l;
  double tmp2;
  for(int i=0 ; i<l ; ++i){
    tmp  = m - i ;
    tmp2 =  1-sq(i/l_double) ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
    tmp = m + i ;
    if( (tmp <= p) && (tmp >= 1) ){
      res(tmp-1) = tmp2 ;
    }
  }
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
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i) * uniform_cpp_unnormalized(m(i),l(i),grid);
    }    
  }
  if(basis == "triangular"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i)/scale_ml( m(i)-1 , l(i)-1 )  * 
        triangular_cpp(m(i),l(i),grid);
    }    
  }
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i) * triangular_cpp_unnormalized(m(i),l(i),grid);
    }    
  }
  if(basis == "gaussian"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i)/scale_ml( m(i)-1 , l(i)-1 )  * 
        gaussian_cpp(m(i),l(i),grid);
    }    
  }
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i) * gaussian_cpp_unnormalized(m(i),l(i),grid);
    }    
  }
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i)/scale_ml( m(i)-1 , l(i)-1 )  *
        Epanechnikov_cpp(m(i),l(i),grid);
    }    
  }
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<K ; ++i){
      res = res + beta_star(i) * Epanechnikov_cpp_unnormalized(m(i),l(i),grid);
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

// Compute all the alternative for the value of the intergral for all m and l. 
arma::cube potential_intervals (mat & X, vec & grid, int l_max, 
                                std::string basis){
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
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "triangular"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "gaussian"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp_unnormalized(i+1,j+1,grid); 
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
  if(basis == "uniform_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = uniform_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "triangular"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "triangular_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = triangular_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "gaussian"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "gaussian_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = gaussian_cpp_unnormalized(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "Epanechnikov"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp(i+1,j+1,grid); 
          X_tmp = trans(X.row(k)) ;
          tmp2  = termwise_product( X_tmp , tmp) ;     
          
          res(i,j,k) = integrate_trapeze(grid, tmp2 );
        }
      }
    }    
  }  
  if(basis == "Epanechnikov_unnormalized"){
    for(int i=0 ; i<p ; ++i){
      for(int j=0 ; j<l_max ; ++j){
        for( int k=0 ; k<n ; ++k){
          tmp   = Epanechnikov_cpp_unnormalized(i+1,j+1,grid); 
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

// Compute the matrix V (for a Ridge Zellner prior) 
// (for Q functional covaribles)                    
mat compute_W_inv_RZ_List (int Q, vec & K, double g, mat & X_tilde, int sum_K,
                       mat & lambda_id0){                                       
  mat W_inv = zeros<mat>(sum_K+1,sum_K+1);
  mat lambda_id = lambda_id0 ;
  
  // mat X_tilde_temp = X_tilde.submat(1,1,sum_K,sum_K) ;
  mat X_tilde_temp = mat_drop_col_k(X_tilde,0);  
  mat u;
  vec s;
  mat v;
  svd(u,s,v,X_tilde_temp);
  
  W_inv(0,0) = 1/lambda_id(0,0); 
  for( unsigned i=1 ; i<sum_K+1; ++i){
    lambda_id(i,i) = lambda_id(i,i) * max(s);
  }
  
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    W_inv.submat(1+count,1+count,K(q)+count,K(q)+count) = 
      ( trans(X_tilde.cols(1+count,K(q)+count)) *
      X_tilde.cols(1+count,K(q)+count) +
      lambda_id.submat(1+count,1+count,K(q)+count,K(q)+count) )  /g;
    count = count + K(q);
  }
  
  return W_inv;
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


// Compute beta_s for a covariate q from the parameters : beta_q, m_q and l_q
vec compute_beta_s (int p, int K, vec & beta, vec & m, vec & l){
  vec res = zeros<vec>(p) ;
  
  // for each intervals
  for(int k=0 ; k<K ; ++k){
    vec indexes;
    // compute indexes
    if( m(k) - l(k) < 1){
      indexes = sequence( 0, m(k) + l(k)-1 , 1);
    }
    if( m(k) + l(k) > p){
      indexes = sequence( m(k) - l(k)-1, p-1 , 1); 
    }
    if( m(k) - l(k) < 1  && m(k) + l(k) > p ){
      indexes = sequence( 0, p-1 , 1); 
    }
    if( m(k) - l(k) >= 1 && m(k) + l(k) <= p){
      indexes = sequence( m(k) - l(k)-1, m(k) + l(k)-1 , 1);
    }
    // for each index
    for(int j=0 ; j<indexes.size() ; ++j ){
      if( beta(k) > 0 )
        res( indexes(j) ) = 1;
      if( beta(k) < 0 )
        res( indexes(j) ) = -1;
    }
  }
  
  return(res) ;
}

// Compute the distance between two beta_s_q 
double dist_beta_s (vec & grid, vec & beta_s1, vec & beta_s2, vec & conf){                  
  double res; 
  vec fct = (beta_s1 - beta_s2);
  res = weighted_norm_fct(grid, fct,conf);
  return(res);
}

//#############################################################################
//############################ Auxiliary functions ############################
//#############################################################################

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

// Update the parameter m_k
int mk_update_List (int count, int k, vec & Y, vec & beta_tilde, double sigma_sq, 
                vec & m_q, vec & l_q, mat & X_tilde, 
                NumericVector & all_intervals_q, vec & all_intervals_dims_q,
                vec & m_alternative_q, vec & eta_tilde,
                vec & phi_m_q, int p_q, std::string prior_beta, int Q,          
                vec & K, double g, int sum_K, mat & lambda_id0, int q,
                vec & beta_s_expert_q, vec & grid_q, double tau, vec & conf_q){
  double aux;
  vec aux2 = zeros<vec>(p_q);
  vec aux3 = zeros<vec>(p_q);
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);  // mis des + 1
  mat X_tilde_temp = X_tilde;
  mat W_inv_temp;
  double beta_s_distance;
  vec beta_s;
  vec m_q_temp = m_q;
  vec beta = beta_tilde.subvec(1+count, K(q)+count );
  
  // Compute the probabilities
  vec probs = ones<vec>(p_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<p_q ; ++i){
    m_q_temp(k) = m_alternative_q(i);
    beta_s = compute_beta_s( p_q, K(q), beta, m_q_temp, l_q );
    beta_s_distance = dist_beta_s(grid_q,beta_s,beta_s_expert_q,conf_q);
      
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_alternative_q(i),
                                        l_q(k),all_intervals_dims_q);
    
    X_tilde_temp.col(count + k + 1) = X_tilde_qki;
    aux2(i) = dot( Y - X_tilde_temp * beta_tilde ,
         Y - X_tilde_temp * beta_tilde)  /(2*sigma_sq) ;
    
    aux3(i) = 1 ; // if prior is not the Rigde Zellner prior
    if(prior_beta == "Ridge_Zellner"){
      W_inv_temp = compute_W_inv_RZ_List(Q,K,g,X_tilde_temp,sum_K,lambda_id0);
      
      aux2(i) = aux2(i) + 
        dot( beta_tilde , W_inv_temp * beta_tilde ) / (2*sigma_sq) +
        tau * beta_s_distance ;
      aux3(i) = sqrt( det(W_inv_temp) );
    }
  }  
  
  double min_aux = min(aux2);
  for(int  i=0 ; i<p_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    // aux3(i)  = aux3(i) / max_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_m_q(i) ;
  }
  // Simulate a mk
  int mk = sample_weight(probs) + 1 ;
  
  return mk;
  }

// Update the parameter l_k
int lk_update_List (int count, int k, vec & Y, vec & beta_tilde, double sigma_sq,
                vec & m_q, vec & l_q, mat & X_tilde, 
                NumericVector & all_intervals_q, vec & all_intervals_dims_q,
                vec & l_alternative_q, vec & eta_tilde, 
                vec & phi_l_q, int lmax_q, int p_q, std::string prior_beta, int Q,         
                vec & K, double g, int sum_K, mat & lambda_id0, int q,
                vec & beta_s_expert_q, vec & grid_q, double tau,  vec & conf_q){
  double aux;
  vec aux2 = zeros<vec>(lmax_q);
  vec aux3 = zeros<vec>(lmax_q);
  mat X_tilde_mqk = mat_drop_col_k(X_tilde,count + k + 1);  // mis des + 1
  mat X_tilde_temp = X_tilde;
  mat W_inv_temp;
  double beta_s_distance;
  vec beta_s;
  vec l_q_temp = l_q;
  vec beta = beta_tilde.subvec(1+count, K(q)+count );
  
  // Compute the probabilities
  vec probs = ones<vec>(lmax_q);
  vec X_tilde_qki = zeros<vec>(all_intervals_dims_q(2)) ;
  for(int  i=0 ; i<lmax_q ; ++i){
    l_q_temp(k) = l_alternative_q(i);
    beta_s = compute_beta_s( p_q, K(q), beta, m_q, l_q_temp );
    beta_s_distance = dist_beta_s(grid_q,beta_s,beta_s_expert_q,conf_q);
    
    X_tilde_qki = all_intervals_extract(all_intervals_q,m_q(k),
                                        l_alternative_q(i),
                                        all_intervals_dims_q);
    
    X_tilde_temp.col(count + k + 1) = X_tilde_qki;
    aux2(i) = dot( Y - X_tilde_temp * beta_tilde , 
         Y - X_tilde_temp * beta_tilde)  /(2*sigma_sq) ;
    
    aux3(i) = 1 ; // if prior is not the Rigde Zellner prior
    if(prior_beta == "Ridge_Zellner"){
      W_inv_temp = compute_W_inv_RZ_List(Q,K,g,X_tilde_temp,sum_K,lambda_id0);
      
      aux2(i) = aux2(i) +
        dot( beta_tilde , W_inv_temp * beta_tilde ) / (2*sigma_sq) +
        tau * beta_s_distance ;
      aux3(i) = sqrt( det(W_inv_temp) );
    }
  }  
  
  double min_aux = min(aux2);
  for(int  i=0 ; i<lmax_q ; ++i){
    aux2(i)  = aux2(i) - min_aux;
    // aux3(i)  = aux3(i) / max_aux;
    probs(i) = aux3(i) * exp( - aux2(i) ) * phi_l_q(i);
  }
  // Simulate a lk
  int lk = sample_weight(probs) + 1 ;
  
  return lk;
}


// update the parameter sigma_sq
double sigma_sq_update_List (vec & Y, vec & beta_tilde, vec & eta_tilde,        
                         mat & W_inv, mat & X_tilde, double a, double b, 
                         int n, int sum_K) {                                    
  double a_tmp     = sum_K+n+1 ;
  double a_star    = a + a_tmp/2 ;
  
  vec Y_tmp        = Y - X_tilde * beta_tilde ;
  double Y_tmp2    = dot(Y_tmp, Y_tmp) ;                                    
  vec beta_tilde_tmp     = beta_tilde - eta_tilde ;
  double beta_tilde_tmp2 = dot(beta_tilde_tmp, W_inv * beta_tilde_tmp) ;
  
  double b_star    = b + 0.5*( Y_tmp2 + beta_tilde_tmp2);
  double res = 1. / (R::rgamma(a_star, 1/b_star) );
  
  return res ;
}

// Compute beta_s from the parameters beta, m and l.
List compute_beta_s_List (vec & p_vec, vec & K, int Q, vec & beta_tilde,          
                         List & m_List, List & l_List){
  List res(Q);
  
  // for each covariables 
  int count = 0;
  for( unsigned q=0 ; q<Q ; ++q){
    vec m = m_List[q];
    vec l = l_List[q];
    int p = p_vec(q) ;
    vec beta = beta_tilde.subvec(1+count, K(q)+count );
    
    vec beta_s_q = compute_beta_s(p,K(q),beta,m,l);

    res[q] = beta_s_q;
    
    count = count + K(q);
  }
  
  // return 
  return(res);
}


// Compute the distance between two beta_s (List) 
vec dist_beta_s_List (List & grids, List & beta_s1_List, List & beta_s2_List,   
                      int Q, List & confidence){
  vec res = zeros<vec>(Q); 
  
  for( unsigned q=0 ; q<Q ; ++q){
    vec grid    = grids[q];
    vec beta_s1 = beta_s1_List[q];
    vec beta_s2 = beta_s2_List[q];
    vec conf_q  = confidence[q];
    res(q) = dist_beta_s(grid,beta_s1,beta_s2,conf_q);
  }
  
  return(res);
}


// Compute a proposal for beta 
vec compute_beta_proposal (vec & beta, vec & eta_tilde, mat & W_inv,
                           double sigma_sq, std::string proposal_distribution, 
                           double tol, double rho){
  vec beta_prop;
  
  if(proposal_distribution == "random_walk"){
    int K = beta.size();
    mat VarCovar = zeros<mat>(K,K);
    for(unsigned k=0 ; k<K ; ++k){
      VarCovar(k,k) = rho;
    }
    vec epsilon = randn<vec>(K); 
    epsilon(0) = 0;
    epsilon = VarCovar * epsilon;
    
    beta_prop = beta + epsilon;
  }
  
  return(beta_prop);
}

// Compute the probability of acceptance alpha
List compute_prob_accep_List (vec & beta_tilde, vec & beta_prop,              
                                mat & X_tilde, double sigma_sq, vec & y, 
                                double mu, List & grids, List & beta_s, 
                                List & beta_s_prop, List & beta_s_expert, int Q, 
                                std::string proposal_distribution, 
                                mat & W_inv, double tau, List & confidence){
  double res;
  double res1;
  double res2;
  vec vec_unity = ones<vec>(y.size());   
  
  double aux;
  // double SSE_mmu;
  
  // std::cout << trans(beta_prop) << std::endl;
  
  vec beta = vec_drop_k(beta_tilde,0);
  vec beta_prop_temp = vec_drop_k(beta_prop,0);
  mat X = mat_drop_col_k(X_tilde,0);
  
  mat W_inv_temp;
  W_inv_temp = W_inv.submat(1,1,W_inv.n_rows-1,W_inv.n_cols-1);
  
  vec aux_y = y - mu * vec_unity;
  mat aux_mat;
  
  aux_mat = trans(X) * X + W_inv_temp;
  vec d_beta_s      = dist_beta_s_List(grids,beta_s     , beta_s_expert,Q,confidence);
  vec d_beta_s_prop = dist_beta_s_List(grids,beta_s_prop, beta_s_expert,Q,confidence);
  
  aux = dot( beta_prop_temp - beta , aux_mat * (beta_prop_temp + beta)) -
    2 * dot( beta_prop_temp - beta , trans(X) * aux_y);
  
  res  = exp( - aux / (2*sigma_sq) - tau * sum(d_beta_s_prop - d_beta_s));
  res1 = aux / (2*sigma_sq);
  res2 = sum(d_beta_s_prop - d_beta_s);
    
  if(res>1) res=1; 
  
  return( List::create(_["alpha"]      = res,
                       _["SSE_diff"]  = res1, 
                       _["dist_diff"] = res2));
}

// update the parameter beta_star with a Metropolis step.
// Take into account of the prior "beta_s"
List beta_tilde_update_beta_s_List (vec & beta_tilde, vec & Y, double sigma_sq, 
                               vec & eta_tilde, mat & W_inv, mat & X_tilde, 
                               double tol,List & beta_s_expert, List & beta_s,
                               std::string proposal_distribution, vec & K, 
                               vec & p,List & m, List & l, double mu, 
                               List & grids, int Q, double rho, double tau,
                               List & confidence){
  double u; 
  unsigned accepted=0;
  List res_alpha;
  // Simulate a proposal 
  vec beta_prop = compute_beta_proposal(beta_tilde, eta_tilde, W_inv, sigma_sq, 
                                        proposal_distribution, tol, rho);

  // Compute beta_s_prop
  List beta_s_prop = compute_beta_s_List(p,K,Q,beta_prop,m,l);

  // Compute the probability of acceptance
  res_alpha = compute_prob_accep_List(beta_tilde, beta_prop, X_tilde, sigma_sq, Y,
                                  mu, grids, beta_s, beta_s_prop, beta_s_expert,
                                  Q,proposal_distribution,W_inv, tau,confidence) ;
  
  double alpha     = res_alpha[0];
  double SSE_diff  = res_alpha[1];
  double dist_diff = res_alpha[2];
  double RSS;
  
  // Simulate a choice
  u = R::runif(0,1) ;
  if(u < alpha){
    beta_tilde = beta_prop  ;
    beta_s     = beta_s_prop;
    accepted   = 1;
  }
  
  RSS = 0.5 / sigma_sq * 
    dot(Y - X_tilde * beta_tilde, Y - X_tilde * beta_tilde);
  
  // return 
  return( List::create(_["beta_tilde"] = beta_tilde,
                       _["beta_s"]     = beta_s, 
                       _["alpha"]      = alpha,
                       _["accepted"]   = accepted,
                       _["SSE_diff"]   = SSE_diff,
                       _["dist_diff"]  = dist_diff,
                       _["RSS"]        = RSS)
  ); 
}


// Update the matrix X_tilde after an update of a m_k or a l_k
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

// Update the intercept 
double mu_update( vec & y, mat & X_tilde, vec & beta_tilde, double sigma_sq, 
                  int n, double v_0){
  double res; 
  mat X = mat_drop_col_k(X_tilde, 0);
  vec beta = vec_drop_k(beta_tilde,0);
  
  double mu_expe = sum(y - X * beta) / (n + 1/v_0); 
  double mu_sd   = sqrt( sigma_sq / (n + 1/v_0) );
  
  res = R::rnorm(mu_expe,mu_sd);
  return(res);
}

// update tau 
double tau_update (double lambda_tau, List & beta_s_expert, List & beta_s, 
                   int Q, List & grids,List & confidence){
  vec d_beta_s = dist_beta_s_List(grids , beta_s , beta_s_expert , Q,confidence) ;
  double lambda = sum(d_beta_s) + lambda_tau;
  
  double res = 1 * R::rgamma(1, 1/ lambda );
  return(res);
}

// A function which compute the weights of the Importance Sampling approximation
mat compute_weights_IS (vec & y, mat & trace, int Q, vec & K, int sum_K, 
                        List & all_intervals, List & all_intervals_dims){
  unsigned iter = trace.n_rows;
  unsigned count ;
  vec trace_tmp ;
  
  unsigned n  = y.size();
  mat weights = zeros<mat>(iter,n);
  mat X_tilde = ones<mat>(n,sum_K+1) ;
  
  vec beta_tilde = zeros<vec>(sum_K+1)  ;
  List m(Q)       ;
  List l(Q)       ;
  double sigma_sq ;
  
  for( unsigned t=0 ; t<iter ; ++t){
    count = 0;
    // posterior value of theta 
    beta_tilde(0) = trace(t,3*sum_K);
    sigma_sq = trace(t,3*sum_K+1);
    for( unsigned q=0 ; q<Q ; ++q){
      trace_tmp = trans(trace.row(t).subvec(3*count , 3*count+  K(q)-1));
      beta_tilde.subvec(1+count,K(q)+count) = trace_tmp ;
      
      trace_tmp = trans(trace.row(t).subvec( 3*count+  K(q) , 3*count+2*K(q)-1));
      m[q] = trace_tmp;
      
      trace_tmp = trans(trace.row(t).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1));
      l[q] = trace_tmp;
      
      count = count + K(q);
    }
    // Compute X
    X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                             X_tilde);
    
    // Compute the weights
    weights.row(t) = trans(exp( ( y - X_tilde * beta_tilde ) / (2*sigma_sq) ));  
  }
  
  return(weights);
}

// a function which compute the utility from the IS weights
double compute_utility(mat & weights, int n, int T){
  double approx_IS;
  double approx_MC;
  double res;
  
  for(unsigned i=0 ; i<n ; ++i){
    for(unsigned t=0 ; t<T ; ++t){
      approx_IS = approx_IS + 1./ weights(t,i);
    }
    approx_IS = log( 1./T * approx_IS);
    approx_MC = approx_MC + approx_IS;
  }
  
  res = 1./n * approx_MC ;
  
  return(res);
}
//##############################################################################
//#################### Gibbs Sampler and Simulated Annealing ###################
//##############################################################################

// Perfom the loop of the MWG algorithm
mat MWG_loop (unsigned iter, vec & beta_tilde, double sigma_sq, List & m, 
               List & l, mat & X_tilde, mat & W_inv, List & beta_s_List, 
               double tau, unsigned print,
               int n, double v_0, vec & Y, vec & eta_tilde, 
               double a, double b, double sum_K, unsigned Q, 
               vec & p, List & all_intervals, List & all_intervals_dims,
               List & probs_m, std::string & prior_beta, vec & K, 
               mat & lambda_id0, vec & lmax, List & probs_l, 
               double tol, double g, List & beta_s_expert_List, List & confidence,
               std::string & proposal_distribution, List & grids, double rho,
               std::string & is_tau, double lambda_tau){
  double mu           ;
  List res_MH         ;
  double alpha        ;
  unsigned accepted   ; 
  double d_beta_s     ;
  double RSS          ;
  
  // Initialize the matrix trace
  mat trace = zeros<mat>(iter+1,3*sum_K+9);  
  unsigned count = 0 ;
  for( unsigned q=0 ; q<Q ; ++q){
    vec m_temp = m[q] ;
    vec l_temp = l[q] ;
    trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
      trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
    trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
    trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;
    count = count + K(q) ;
  }
  trace(0,3*sum_K  ) = beta_tilde(0) ;  
  trace(0,3*sum_K+1) = sigma_sq      ;  
  trace(0,3*sum_K+2) = tau           ;
  
  for(unsigned i=1  ; i < iter+1 ; ++i ) {      
    // Progress
    // std::cout << "\t " << i << std::endl;
    if( print == 1){
      if( (i % (iter / 10)  == 0)) 
            std::cout << "\t " << i / (iter / 100) << "%" << std::endl;
    }
    
    // update mu
    beta_tilde(0) = mu_update(Y, X_tilde, beta_tilde, sigma_sq, n, v_0);
    
    // update sigma_sq
    sigma_sq = sigma_sq_update_List(Y,beta_tilde,eta_tilde,W_inv,X_tilde,a,b,
                                    n,sum_K) ;
    
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
      
      vec beta_s_expert_q = beta_s_expert_List[q];
      vec grid_q = grids[q];
      vec conf_q = confidence[q];
      
      for(int k=0 ; k<K(q) ; ++k){
        // update m_k
        // if(prior_beta == "Ridge_Zellner"){
        m_q(k) = mk_update_List(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
            all_intervals_q,all_intervals_dims_q,m_alternative_q,
            eta_tilde,probs_m_q,p_q, prior_beta,Q,K,g,sum_K,lambda_id0,q,
            beta_s_expert_q,grid_q, tau,conf_q);
          
        // update the value "X_tilde"
        X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                                 X_tilde);
        
        // Do not need to update W_inv for now (because the "mk_update_List"
        // function does not require this matrix as an input. It have to be
        // internal computed for each possibility of m_k)
      }
      
      // update beta_s_List;
      vec beta_q = beta_tilde.subvec(1+count, K(q)+count );
      beta_s_List[q] = compute_beta_s(p_q,K(q),beta_q,m_q,l_q);
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
      int p_q = p(q);
      NumericVector all_intervals_q = all_intervals[q];
      vec all_intervals_dims_q      = all_intervals_dims[q];
      vec l_alternative_q = sequence(1,lmax_q,1) ;
      vec probs_l_q         = probs_l[q];
      
      vec beta_s_expert_q = beta_s_expert_List[q];
      vec grid_q = grids[q];
      vec conf_q = confidence[q];
      
      for(int k=0 ; k<K(q) ; ++k){
        // update l_k
        // if(prior_beta == "Ridge_Zellner"){
        l_q(k) = lk_update_List(count,k,Y,beta_tilde,sigma_sq,m_q,l_q,X_tilde,
            all_intervals_q,all_intervals_dims_q,l_alternative_q,
            eta_tilde,probs_l_q,lmax_q, p_q,prior_beta,Q,K,g,sum_K,lambda_id0,q,
            beta_s_expert_q,grid_q, tau,conf_q);
        
        
        vec beta_s_expert_q = beta_s_expert_List[q];
        vec grid_q = grids[q];
        // update the value "X_tilde"
        X_tilde = update_X_tilde(Q,K,all_intervals,all_intervals_dims,m,l,
                                 X_tilde);
        
        // Do not need to update W_inv for now (because the "mk_update_List"
        // function does not require this matrix as an input. It have to be
        // internal computed for each possibility of l_k)
      }
      
      // update beta_s_List;
      vec beta_q = beta_tilde.subvec(1+count, K(q)+count );
      beta_s_List[q] = compute_beta_s(p_q,K(q),beta_q,m_q,l_q);
      // Update the m_q value
      l[q] = l_q;
      // Update count
      count = count + K(q);
    }
    
    // update the value "W_inv" (only after the updating of the m's and l's)
    // if(prior_beta == "Ridge_Zellner")
    W_inv = compute_W_inv_RZ_List (Q,K, g, X_tilde,sum_K,lambda_id0) ;
    
    // update beta : Metropolis step
    mu = beta_tilde(0);
    res_MH =
      beta_tilde_update_beta_s_List(beta_tilde,Y,sigma_sq, eta_tilde, W_inv,
                                    X_tilde,tol, beta_s_expert_List,beta_s_List,
                                    proposal_distribution, K, p, m,l, mu, grids,
                                    Q,rho,tau,confidence);
    vec  res_MH_1 = res_MH[0];
    List res_MH_2 = res_MH[1];
    alpha = res_MH[2];
    accepted = res_MH[3];
    
    double SSE_diff  = res_MH[4];
    double dist_diff = res_MH[5];
    
    RSS         = res_MH[6];
    beta_tilde  = res_MH_1;
    beta_s_List = res_MH_2;
    
    if(is_tau == "random"){
      // update tau 
      d_beta_s = sum(dist_beta_s_List(grids , beta_s_List , beta_s_expert_List ,
                                      Q,confidence)) ;
      tau = tau_update(lambda_tau,beta_s_expert_List,beta_s_List,Q,grids,
                       confidence);
    }
    
    // update the matrix trace
    count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      vec m_temp = m[q] ;
      vec l_temp = l[q] ;
      trace.row(i).subvec( 3*count        , 3*count+  K(q)-1) =
        trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
      trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1) = trans(m_temp);
      trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1) = trans(l_temp);
      
      count = count + K(q) ;
    }
    trace(i,3*sum_K  ) = beta_tilde(0) ;
    trace(i,3*sum_K+1) = sigma_sq      ;
    trace(i,3*sum_K+2) = tau           ;
    trace(i,3*sum_K+3) = alpha         ;
    trace(i,3*sum_K+4) = accepted      ;
    trace(i,3*sum_K+5) = SSE_diff      ;
    trace(i,3*sum_K+6) = dist_diff     ;
    trace(i,3*sum_K+7) = d_beta_s      ;
    trace(i,3*sum_K+8) = RSS           ;
  }
  
  // Return the trace
  return(trace);
}


// Perform a MWG to sample from the posterior distribution.
// The metropolis step is for update the beta parameter. 
// [[Rcpp::export]]
List Bliss_MWG_multiple_cpp (int Q, vec & Y, List & X, int iter,  
                               List & grids, vec & K, vec & lmax, 
                               vec & eta_tilde, double a, double b, 
                               List & probs_m, List & probs_l, 
                               std::string prior_beta, 
                               double g, double lambda, mat & V_tilde, 
                               double tol,   CharacterVector & basis,
                               List & beta_s_expert_List, List & confidence,
                               std::string proposal_distribution, double rho, 
                               unsigned iter_tau, vec & tau_vec,
                               std::string choose_tau_vec,
                               std::string & is_tau, double lambda_tau,
                               unsigned tau_vec_size){
  std::cout << "Metropolis within Gibbs: " <<  std::endl;
  std::cout << "\t Initialization." <<  std::endl;
  
  //########## Initialization of some quantities
  // Compute the value of n and the p's                                         
  int n = as<mat>(X[0]).n_rows ;
  vec p = zeros<vec>(Q)        ;
  for(int i=0 ; i<Q ; ++i){
    p(i) = as<mat>(X[i]).n_cols; 
  }
  
  // Compute projection of the x_i on all the intervals
  List scale_ml(Q)          ; // scale_ml is used to normalize the predictors
  List all_intervals(Q)     ; // will be contain all the projections
  List all_intervals_dims(Q); // will be contain the dim of the all_intervals's
  for( int q=0 ; q<Q ; ++q){
    arma::cube temp = potential_intervals_List (X, grids, lmax, basis,q) ;
    scale_ml[q]      = temp.slice(n);
    
    temp = temp.subcube(0,0,0,p(q)-1,lmax(q)-1,n-1);
    all_intervals[q] = temp;
    
    vec temp2 = zeros<vec>(3);
    temp2(0) = p(q)    ;
    temp2(1) = lmax(q) ;
    temp2(2) = n    ;
    all_intervals_dims[q] = temp2;
  }
  
  
  // Compute the matrix of lambda for the Ridge penalty of Rigde Zellner prior
  int sum_K = sum(K);
  mat lambda_id0  = zeros<mat>(sum_K+1,sum_K+1) ;
  // if(prior_beta == "Ridge_Zellner"){                                            
  lambda_id0(0,0) = 100*var(Y);
  for( unsigned i=1 ; i<sum_K+1; ++i){                            
    lambda_id0(i,i) = lambda ;
    // }
  }
  
  // Some objects
  vec beta_tilde  ;
  double sigma_sq ;
  double d_beta_s ;
  mat W_inv       ;
  double tau      ;
  
  mat X_tilde = ones<mat>(n,sum_K+1) ; 
  double v_0  = V_tilde(0,0)         ;
  
  List m(Q)            ;
  List l(Q)            ;
  List res_MH(2)       ;
  List beta_s_List(Q)  ;
  List m_alternative(Q);
  List l_alternative(Q);
  double RSS;
  List res_CV;
  
  for( int q=0 ; q<Q ; ++q){
    l_alternative[q] = sequence(1,lmax(q),1);
  }
  for( int q=0 ; q<Q ; ++q){
    m_alternative[q] = sequence(1,p(q),1);
  }
  
  //########## Determine the start point
  std::cout << "\t Determine the start point." <<  std::endl;
  // Initialization of sigma_sq
  sigma_sq  = var(Y) ;
  
  // Initialization of the middle and length of the intervals
  for( unsigned q=0 ; q<Q ; ++q){  
    vec probs_m_temp = probs_m[q];
    m[q]         = sample_weight(K(q),probs_m_temp) + 1 ; 
    vec probs_l_temp = probs_l[q];
    l[q]         = sample_weight(K(q),probs_l_temp) + 1 ;
  }
  
  // Initialize the current X_tilde matrix (which depend on the intervals)
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
  
  // Initialize the current W_inv matrix (which depend on the intervals)
  // if(prior_beta == "Ridge_Zellner")
  W_inv = compute_W_inv_RZ_List(Q,K, g, X_tilde,sum_K,lambda_id0) ;
  
  // Initialization of beta_tilde (without taking account of beta_s)
  beta_tilde = mvrnormArma( eta_tilde , ginv_cpp(W_inv,tol) , sigma_sq) ; 
  
  // Compute beta_s_List
  beta_s_List = compute_beta_s_List(p,K,Q,beta_tilde,m,l);
  
  // update tau (if tau is chosen to be random)
  if(is_tau == "random"){
    tau = 0;
  }
  
  //########## Determine a value for rho
  if(rho == 0){
    std::cout << "\t Determine rho." <<  std::endl;
    int success = 0;
    int count2  = 0;
    rho = 1;
    unsigned iter_rho = 1e2;
    mat trace_rho;
    
    while( (success==0) && (count2 < 20) ){
      std::cout << "\t \t if rho = "<< rho << " ?";
      vec    accepteds ;
      vec    beta_tilde_temp = beta_tilde ;
      double tau_temp        = tau;
      double sigma_sq_temp   = sigma_sq ;
      List   m_temp          = m ;
      List   l_temp          = l ;
      mat    X_tilde_temp    = X_tilde ;
      mat    W_inv_temp      = W_inv ;        
      List   beta_s_List_temp = beta_s_List;
      
      trace_rho = MWG_loop(iter_rho,beta_tilde_temp,sigma_sq_temp,m_temp,l_temp,
                           X_tilde_temp,W_inv_temp,beta_s_List_temp,tau_temp,0,
                           n,v_0,Y,eta_tilde,a,b,sum_K,Q,p,all_intervals,
                           all_intervals_dims,probs_m,prior_beta,K,lambda_id0,
                           lmax,probs_l,tol,g,beta_s_expert_List,confidence,
                           proposal_distribution,grids,rho,is_tau,lambda_tau);
      RSS       = mean(trace_rho.col(3*sum_K+8));
      accepteds = trace_rho.col(3*sum_K+4);
      
      // Check the values of alphas 
      if( mean(accepteds) > 0.5 ){
        std::cout << " Then mean acceptance rate = " << mean(accepteds) << ". So, check another value." <<  std::endl;
        rho = rho * 2;
      }
      if( mean(accepteds) < 0.3 ){
        std::cout << " Then mean acceptance rate = " << mean(accepteds) << ". So, check another value." <<  std::endl;
        rho = rho / 2;
      }
      if( (mean(accepteds) < 0.5) && (mean(accepteds) > 0.3) ){
        std::cout << " Then mean acceptance rate = " << mean(accepteds) << ". So, we fixe rho to "<< rho << "." <<  std::endl;
        success = 1;
      }
      
      count2 = count2 +1;
    }
  }
  
  //########## Determine a value for tau
  if(is_tau != "random"){
    if(choose_tau_vec == "data_driven"){
      tau_vec = sequence2(0,RSS,tau_vec_size);
    }
    vec utility(tau_vec.size());
    List weights_IS(tau_vec.size());
    if(tau_vec.size() > 1){
      std::cout << "\t Determine tau." <<  std::endl;
      // std::cout << "\t \t tau_vec = "<< trans(tau_vec) <<  std::endl;
      std::cout << "\t \t tau_vec = (" ;
      for(unsigned t=0 ; t<tau_vec.size() ; ++t){
        std::cout << " " << tau_vec(t) << " ";
      }
      std::cout << ")" << std::endl;
      
      // load start point 
      vec beta_tilde_tmp   = beta_tilde ;
      double sigma_sq_tmp  = sigma_sq ;
      List m_tmp           = m ;
      List l_tmp           = l ;
      mat X_tilde_tmp      = X_tilde ;
      mat W_inv_tmp        = W_inv ;
      List beta_s_List_tmp = beta_s_List ;
      
      // initialization of some objects
      int burnin_tau = iter_tau/10;
      List trace_tau(tau_vec.size()); 
      List var_weights(tau_vec.size());
      
      mat trace_tmp2 ; 
      mat trace_tmp ; 
      mat weights_tmp;
      vec var_weights_vec(n);
      
      
      // For each tau_k :
      for(unsigned t=0 ; t<tau_vec.size() ; ++t){
        std::cout << "\t \t Compute the utility for tau = " << tau_vec(t) << "." << std::endl;
        // Simulate from the posterior distribution
        trace_tmp2 = MWG_loop(iter_tau,
                              beta_tilde_tmp,sigma_sq_tmp,m_tmp,l_tmp,X_tilde_tmp,
                              W_inv_tmp,beta_s_List_tmp,
                              tau_vec(t),0,
                              n,v_0,Y,eta_tilde, a, b, sum_K, Q, p, all_intervals,
                              all_intervals_dims, probs_m, prior_beta, K, 
                              lambda_id0, lmax, probs_l, tol, g,
                              beta_s_expert_List, confidence, 
                              proposal_distribution, grids, rho,
                              is_tau,lambda_tau);
        trace_tmp = trace_tmp2.rows(burnin_tau,trace_tmp2.n_rows-1);
        trace_tau[t] = trace_tmp; 
        // Compute the weights 
        weights_tmp = compute_weights_IS(Y, trace_tmp,Q, K, sum_K, all_intervals,
                                         all_intervals_dims);
        weights_IS[t] = weights_tmp;
        // Compute the variance of the weights
        for(unsigned i=0 ; i<n ; ++i){
          var_weights_vec(i) = var( weights_tmp.col(i)) ;
        }
        var_weights[t] = var_weights_vec;
        
        // Compute the expected utility of tau_k
        utility(t) = compute_utility(weights_tmp,n,iter_tau-burnin_tau);
      }
      // Choose tau_star
      unsigned decision = 0; 
      for(unsigned t=1 ; t<tau_vec.size() ; ++t){
        if( utility(t) > utility(decision) ) decision = t;
      }
      tau = tau_vec(decision);
      std::cout << "\t \t tau is fixed to " << tau <<  std::endl;
    }
    res_CV = List::create(_["utility"]=utility,
                          _["weights_IS"]=weights_IS,
                          _["tau_vec"]=tau_vec,
                          _["tau"]=tau);
  }
  
  //########## The MWG loop
  std::cout << "\t Start the MWG loop." <<  std::endl;
  mat trace = MWG_loop(iter,beta_tilde,sigma_sq,m,l,X_tilde,W_inv,beta_s_List, 
                       tau,1,
                       n,v_0,Y,eta_tilde, a, b, sum_K, Q, p, all_intervals,
                       all_intervals_dims, probs_m, prior_beta, K, lambda_id0, 
                       lmax, probs_l, tol, g, beta_s_expert_List, confidence,
                       proposal_distribution, grids, rho,is_tau,lambda_tau);
  //########## return the trace and the parameters
  std::cout << "\t Return the result." <<  std::endl;
  return  List::create(_["trace"]=trace, 
                       _["IS_LOO"]=res_CV,
                       _["param"]=List::create(_["a"]=a,
                                          _["b"]=b,
                                          _["phi_m"]=probs_m,
                                          _["phi_l"]=probs_l,
                                          _["K"]=K,
                                          _["eta_tilde"]=eta_tilde,
                                          _["l_max"]=lmax,
                                          _["all_intervals"]=all_intervals,
                                          _["grids"]=grids,
                                          _["scale_ml"]=scale_ml
                       ));
  
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
