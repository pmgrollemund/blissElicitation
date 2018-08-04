

colnames(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace)
hist(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"beta_1_q_1"],xlab="",ylab="",main="",nclass=50)
hist(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),1],xlab="",ylab="",main="",nclass=50)

hist(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"m_1_q_1"],xlab="",ylab="",main="",nclass=50)
table(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"m_1_q_1"])
hist(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),4],xlab="",ylab="",main="",nclass=50)
table(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),4])

hist(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"l_1_q_1"],xlab="",ylab="",main="",nclass=50)
table(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"l_1_q_1"])
table(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),7])


hist(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"accepted"],xlab="",ylab="",main="",nclass=50)


hist(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"sigma_sq"],xlab="",ylab="",main="",nclass=50,xlim=c(0,1200),freq=F)
hist(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),11],xlab="",ylab="",main="",nclass=50,xlim=c(0,1200),freq = F)

plot(density(res_Bliss_Pernes_bliss$res.Gibbs_Sampler$trace[-(1:2e3),11]),main="",xlab="",ylab="")
lines(density(res_Bliss_Pernes_MWG$res.Gibbs_Sampler$trace[-(1:2e3),"sigma_sq"]),col=2)











// Initialization of beta_tilde
beta_tilde = mvrnormArma( zeros<vec>(sum_K+1) , ginv_cpp(Sigma_inv,tol) , sigma_sq) ; 

// Initialize the matrix trace
mat trace = zeros<mat>(iter+1,3*sum_K+5);  
int count = 0;
for( unsigned q=0 ; q<Q ; ++q){
  vec m_temp = m[q] ;
  vec l_temp = l[q] ;
  trace.row(0).subvec( 3*count        , 3*count+  K(q)-1) =
    trans(beta_tilde.subvec( 1+count , K(q)+count )) ;
  trace.row(0).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)   = trans(m_temp) ;
  trace.row(0).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)   = trans(l_temp) ;
  
  trace(0,3*sum_K  ) = beta_tilde(0) ;  
  trace(0,3*sum_K+1) = sigma_sq      ;  
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

// // Determine a value for rho
rho = sigma_sq/3; // pour determiner en utilisant une boucle gibbs voir brouillon

// The MWG loop
std::cout << "\t Start the MWG loop." <<  std::endl;
for(unsigned i=1  ; i < iter+1 ; ++i ) {      
  // Progress
  // std::cout << "\t " << i << std::endl;
  if( i % (iter / 10)  == 0) 
    std::cout << "\t " << i / (iter / 100) << "%" << std::endl;
  
  // update sigma_sq  
  res_update = sigma_sq_update_RWM(sigma_sq,rho,weights,average_weights,n_e,
                                   E,y_bar,X_tilde_bar,beta_tilde,coef1,
                                   Sigma_inv) ; 
  res_update_1 = res_update[0];
  sigma_sq = res_update_1;
  res_update_1 = res_update[1];
  alpha    = res_update_1;
  res_update_1 = res_update[2];
  accepted = res_update_1;
  res_update_1 = res_update[3];
  sigma_sq_prop = res_update_1;
  
  
  // // Compute the matrix W
  W = compute_W(weights,average_weights,sigma_sq,n_e,E,n);
  
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
  Sigma_beta_tilde_inv = Sigma_inv/sigma_sq + 
    trans(X_tilde_bar) * W * X_tilde_bar;
  
  if(success){
    // update the beta_tilde   
    mu_beta_tilde = trans(X_tilde_bar) * W * y_bar;
    beta_tilde = beta_tilde_update(Sigma_beta_tilde_inv,mu_beta_tilde ,tol) ; 
    
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
      trace(i,3*sum_K+1) = sigma_sq      ;  
      count = count + K(q) ;
    }
    
    trace(i,3*sum_K+2) = alpha         ;
    trace(i,3*sum_K+3) = accepted      ;
    
    trace(i,3*sum_K+4) = sigma_sq_prop ;  
  }else{ //... go back to the beginning of the updating process.
    i     = i - 1 ;
    count = 0;
    for( unsigned q=0 ; q<Q ; ++q){
      beta_tilde.subvec( 1+count , K(q)+count ) =
        trans(trace.row(i).subvec( 3*count , 3*count+  K(q)-1))  ;
      m[q] = trans(trace.row(i).subvec( 3*count+  K(q) , 3*count+2*K(q)-1)) ;
      l[q] = trans(trace.row(i).subvec( 3*count+2*K(q) , 3*count+3*K(q)-1)) ;
      
      beta_tilde(0) = trace(i,3*sum_K  ) ;  
      sigma_sq      = trace(i,3*sum_K+1) ;  
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




###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""
###################################"""""""


if(rho == 0){
  std::cout << "\t Determine rho." <<  std::endl;
  int success2 = 0;
  int count2 = 0;
  rho = sigma_sq/3;
  unsigned iter_rho = 1e2;
  
  while( (success2==0) && (count2 < 20) ){
    std::cout << "\t \t if rho = "<< rho << " ?";
    vec    accepteds        = zeros<vec>(iter_rho);
    vec    beta_tilde_temp  = beta_tilde ;
    double sigma_sq_temp    = sigma_sq ;
    List   m_temp           = m ;
    List   l_temp           = l ;
    mat    X_tilde_temp     = X_tilde ;
    mat    X_tilde_bar_temp = X_tilde_bar ;
    mat    Sigma_inv_temp   = Sigma_inv;        
    mat    W_temp           = W;
    
    double sigma_sq_save ;
    vec beta_tilde_save  ;
    List m_save          ;
    List l_save          ;
    mat X_tilde_bar_save ;
    mat D_w_sigma_save   ;
    mat W_inv_save       ;
    mat X_tilde_save     ;
    mat Sigma_inv_save   ;
    mat W_save           ;
    
    for(unsigned i=0  ; i < iter_rho ; ++i ) {
      // Save parameters
      sigma_sq_save    = sigma_sq_temp;
      X_tilde_bar_save = X_tilde_bar_temp;
      X_tilde_save     = X_tilde_temp;
      m_save           = m_temp;
      l_save           = l_temp;
      beta_tilde_save  = beta_tilde_temp;
      Sigma_inv_save   = Sigma_inv_temp;
      W_save           = W_temp;
      
      // update sigma_sq (RWM)
      res_update = sigma_sq_update_RWM(sigma_sq_temp,rho,weights,
                                       average_weights,n_e,E,y_bar,X_tilde_bar_temp,
                                       beta_tilde_temp,coef1,Sigma_inv_temp) ; 
      res_update_1 = res_update[0];
      sigma_sq_temp = res_update_1;
      res_update_1 = res_update[2];
      accepteds(i) = res_update_1;
      
      // // Compute the matrix W
      W_temp = compute_W(weights,average_weights,sigma_sq_temp,n_e,E,n);
      
      // update m  
      count = 0 ; 
      // count is used to browse some vec/mat when p(q) is not constant wrt q.
      for( unsigned q=0 ; q<Q ; ++q ){
        // Compute some quantities which do not vary with k 
        vec m_q = m_temp[q];
        vec l_q = l_temp[q];
        int p_q = p(q);
        NumericVector all_intervals_q = all_intervals[q];
        vec all_intervals_dims_q      = all_intervals_dims[q];
        vec m_alternative_q = sequence(1,p_q,1) ; 
        vec probs_m_q       = probs_m[q];
        
        for(int k=0 ; k<K(q) ; ++k){
          // update m_k
          m_q(k) = mk_update_List(count,k,y_bar,beta_tilde_temp,sigma_sq_temp,
                                  m_q,l_q,X_tilde_bar_temp,X_tilde_temp,all_intervals_q,
                                  all_intervals_dims_q,m_alternative_q, probs_m_q,p_q, Q,K,g,
                                  sum_K,lambda_id0,  W_temp);
        }
        
        // update the value "X_tilde_bar" and "X_tilde"
        X_tilde_bar_temp = update_X_tilde(Q,K,all_intervals,all_intervals_dims,
                                          m_temp,l_temp,X_tilde_bar_temp);
        X_tilde_temp     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                          m_temp,l_temp,X_tilde_temp);
        
        // Update the m_q value
        m_temp[q] = m_q;
        // Update count 
        count = count + K(q);
      }
      
      // update l
      count = 0 ; 
      // count is used to browse some vec/mat when p(q) is not constant wrt q.
      for( unsigned q=0 ; q<Q ; ++q ){
        // Compute some quantities which do not vary with k 
        vec m_q = m_temp[q];
        vec l_q = l_temp[q];
        int lmax_q = lmax(q);
        NumericVector all_intervals_q = all_intervals[q];
        vec all_intervals_dims_q      = all_intervals_dims[q];
        vec l_alternative_q = sequence(1,lmax_q,1) ; 
        vec probs_l_q         = probs_l[q];
        
        for(int k=0 ; k<K(q) ; ++k){
          // update l_k
          l_q(k) = lk_update_List(count,k,y_bar,beta_tilde_temp,sigma_sq_temp,
                                  m_q,l_q,X_tilde_bar_temp,X_tilde_temp,all_intervals_q,
                                  all_intervals_dims_q,l_alternative_q, probs_l_q,lmax_q, Q,K,g,
                                  sum_K,lambda_id0,W_temp); 
          
          // update the value "X_tilde_bar" and "X_tilde"
          X_tilde_bar_temp = update_X_tilde(Q,K,all_intervals,all_intervals_dims,
                                            m_temp,l_temp,X_tilde_bar_temp);
          X_tilde_temp     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                            m_temp,l_temp,X_tilde_temp);
        }
        
        // Update the l_q value
        l[q] = l_q;
        // Update count 
        count = count + K(q);
      }
      
      // update the value "Sigma_inv_temp" (only after the updating of the m's and l's)
      Sigma_inv_temp = compute_Sigma_inv (Q,K, g, X_tilde_temp,sum_K,lambda_id0) ;  
      
      // update the matrix Sigma_beta_tilde (only after the updating of 
                                             // the m's and l's)   
      Sigma_beta_tilde_inv = Sigma_inv_temp/sigma_sq_temp + 
        trans(X_tilde_bar_temp) * W_temp * X_tilde_bar_temp;
      
      if(success){
        // update the beta_tilde   
        mu_beta_tilde = trans(X_tilde_bar_temp) * W_temp * y_bar;
        beta_tilde_temp = beta_tilde_update(Sigma_beta_tilde_inv,mu_beta_tilde ,tol) ; 
      }else{ //... go back to the beginning of the updating process.
        i     = i - 1 ;
        sigma_sq_temp    = sigma_sq_save;
        X_tilde_bar_temp = X_tilde_bar_save;
        X_tilde_temp     = X_tilde_save;
        m_temp           = m_save;
        l_temp           = l_save;
        beta_tilde_temp  = beta_tilde_save;
        mat Sigma_inv_temp = Sigma_inv_save   ;
        mat W_temp         = W_save           ;
      }
      
      // update the value "X_tilde_bar" and "X_tilde"
      X_tilde_bar_temp = update_X_tilde(Q,K,all_intervals,all_intervals_dims,
                                        m_temp,l_temp,X_tilde_bar_temp);
      X_tilde_temp     = update_X_tilde(Q,K,all_intervals_obs,all_intervals_dims_obs,
                                        m_temp,l_temp,X_tilde_temp); 
    }
    
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
      success2 = 1;
    }
    
    count2 = count2 +1;
  }
}