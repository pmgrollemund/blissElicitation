#########################################################
#                                                       #
#    Bliss method : MCMC and optimization algorithms    #
#                                                       # 
#########################################################
################################# ----
# Bliss_MWG_multiple
################################# ----
# description : Perform the Gibbs Sampler to sample from the posterior 
#               distribution of the Bliss model.
# value : a list containing :
#             trace : a matrix. Each row is an iteration of the Gibbs Sampler. 
#                The colums are: beta1star,m1,l1,...,betaQstar,mQ,lQ,mu,sigma2.
#             param : a list containing a, b, V_tilde, K, eta_tilde, l_max, grids and 
#                         all_intervals (see the function potential_intervals).
# argument 1 : data, a list containing 1) the number of covariates, 2) the functions x_qi(t) observed at
#              grids of time points and 3) the outcome values y_i.
# argument 2 : param, a list containing :
#              iter : an integer, the number of iterations of the
#                 Gibbs sampler algorithm.
#              grids : a list of numerical vectors, the qth vector is the grid of observation points for the qth covariate.
#              K : a vector of integers, hyperparameters of the Bliss model, the number of intervals for the Q covariates.
#              l_max : a vector of integers, hyperparameters of the Bliss model. Beware, give the index 
#                 corresponding to the value of the maximum width of the intervals (optional)
#              basis : a vector of characters among : "uniform" (default), 
#                 "epanechnikov", "gauss" and "triangular". This indicates the 
#                 shapes of the Q coefficient functions on the intervals. (optional)
#              eta_tilde : a numerical vector of length (1+sum(K)), hyperparameter of the Bliss model. 
#                 By default, eta is the vector (0,...,0). (optional)
#              V_tilde : a matrix of dimension (1+sum(K))*(1+sum(K)), hyperparameter 
#                 of the Bliss model. (optional) (nasty code)
#              a : a nonnegative value, hyperparameter of the Bliss model. By 
#                 default, a = 0.1. (optional)
#              b : a nonnegative value, hyperparameter of the Bliss model. By 
#                 default, b = 0.1. (optional)
#              g : hyperparameter of the Bliss model,  a nonnegative value, 
#                the coefficient of the Zellner prior.
#              phi_m : a list of numerical vectors. The priors of the mq, q=1,...,Q. If not specified, a uniform 
#                 distribution is used for phi_m[[q]].
#              phi_l : a list of numerical vectors. The priors of the lq, q=1,...,Q. If not specified, a uniform 
#                 distribution is used for phi_l[[q]].
#              phi_l_mean : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the mean of the Gamma prior of lq.
#              phi_l_sd : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the standard deviation of the Gamma prior of lq.
#              prior_beta : a character string, which indicates the prior on the
#                beta_star[q]. The possible values are : 
#                1) "diag" (default) for a diagonal matrix prior,
#                2) "Ridge_Zellner" for the Ridge Zellner prior,
#              G : a matrix, the precision matrix (inverse of covariance) of 
#                 y | \mu, \beta, \sigma, m, l.
#              elicited : a boolean value, indicates if there is some 
#                 elicited data in y nor x.
#              conf : a numerical value, the trust that you have in the 
#                 expert's opinions. 
#              n0 : a positive integer, the number of non-elicited observations.

# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# 
# res_Bliss_Gibbs_Sampler <- Bliss_Gibbs_Sampler(data,param)
# K       <- param$K 
# theta_1 <- res_Bliss_Gibbs_Sampler$trace[1,]
# theta_1

# beta_star <- theta_1[1+1:K]
# m         <- theta_1[1:K+K+2]
# l         <- theta_1[1:K+2*K+2]
# beta_1    <- beta_build_cpp(beta_star,m,l,param$grid,param$p,K,"uniform")
# plot(param$grid,beta_1,type="l")
Bliss_MWG_multiple <- function(data,param){
  # Initialize
  x_mult <- data$x_mult
  y <- data$y
  Q <- length(x_mult)
  # load objects
  iter   <- param[["iter"]]
  grids   <- param$grids
  grids_l <- list()
  for (q in 1:Q){
    grids_l[[q]] <- (grids[[q]] - grids[[q]][1])[-1]
  }
  K      <- param$K
  p      <- param[["p"]]
  
  # load optional objects
  l_max   <- param[["l_max"]]
  eta_tilde     <- param$eta_tilde
  V_tilde <- param[["V_tilde"]]
  lambda  <- param[["lambda"]]
  g       <- param[["g"]]
  a       <- param$a
  b       <- param[["b"]]
  phi_m   <- param[["phi_m"]]
  phi_l   <- param[["phi_l"]]
  basis   <- param[["basis"]]
  phi_l_mean <- param[["phi_l_mean"]]
  phi_l_sd   <- param[["phi_l_sd"]]
  prior_beta <- param$prior_beta
  beta_s_expert <- param[["beta_s_expert"]]
  conf_expert <- param[["conf_expert"]]
  beta_s_dist <- param[["beta_s_dist"]]
  MH_proposal <- param[["MH_proposal"]]
  rho         <- param[["rho"]]
  
  is_tau      <- param[["is_tau"]]
  lambda_tau  <- param[["lambda_tau"]]
  iter_tau    <- param[["iter_tau"]]
  tau_vec     <- param[["tau_vec"]]
  tau_vec_size     <- param[["tau_vec_size"]]
  
  
  
  # Initialize the necessary unspecified objects 
  if( is.null(beta_s_expert)){
    stop("Please specify the elicited signed function.")
  }else{
    if (!is.list(beta_s_expert))
      stop("beta_s_expert has to be list.")
  }
  if(is.null(MH_proposal)) MH_proposal <- "random_walk"
  p <- numeric()
  for(q in 1:Q){
    p[q] <- length(grids[[q]])
  }
  if(is.null(prior_beta))
    stop("Please specify a value for the vector prior_beta.")
  if(is.null(K)) stop("Please specify a value for the vector K.")
  if(!is.null(K)){
    for (q in 1:Q){
      if(is.na(K[q])) 
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(eta_tilde))     eta_tilde <- rep(0,1+sum(K))
  if(!is.null(eta_tilde)){
    for (q in 1:Q){
      if(is.na(eta_tilde[q])) 
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(a))       a       <- 1e-1
  if(is.null(b))       b       <- 1e-1
  if(is.null(basis)){
    basis <- character()
    for (q in 1:Q){
      basis[q] <- "uniform"
    }
  }  
  if(!is.null(basis)){
    for (q in 1:Q){
      if(is.na(basis[q])){basis[q] <- "uniform"}
    }
  }      
  if (is.null(phi_l_mean)) phi_l_mean <- rep(NA,Q)
  if (is.null(phi_l_sd)) phi_l_sd <- rep(NA,Q)
  for(q in 1:Q){
    if(!is.null(phi_l[[q]]) && is.character(phi_l[[q]]) && 
       phi_l[[q]] != "Gamma") 
      stop("The qth component of phi_l should be a numeric vector or 'Gamma'.")
    if(!is.null(phi_l[[q]]) && is.character(phi_l[[q]]) &&
       phi_l[[q]] == "Gamma"){
      
      if(is.na(phi_l_mean[q])) phi_l_mean[q] <- 
          diff(range(grids[[q]]))/5 + grids[[q]][1]
      if(is.na(phi_l_sd[q]))   phi_l_sd[q]   <- 
          diff(range(grids[[q]]))/5 
      phi_l[[q]] <- prior_l(phi_l_mean[q]/K[q],phi_l_sd[q]/K[q],grids_l[[q]])   
    }
  }
  if(is.null(phi_l)){
    phi_l <- list()
    for (q in 1:Q){
      if(is.null(l_max)) l_max <- floor(p/5)
      if(!is.null(l_max) & is.na(l_max[q])){l_max[q] <- floor(p[q]/5)}
      phi_l[[q]] <- rep(1/l_max[q],l_max[q])
    }  
  }
  
  for (q in 1:Q){
    l_max[q] <- length(phi_l[[q]])
  }
  if(!is.null(phi_m)){
    for (q in 1:Q){
      if(is.na(phi_m[[q]])){phi_m[[q]] <- rep(1/p[q],p[q])}
    }    
  }  
  if(is.null(phi_m)){
    phi_m <- list()
    for (q in 1:Q){
      phi_m[[q]] <- rep(1/p[q],p[q])
    }  
  }
  if(is.null(prior_beta)) prior_beta <- "diag"
  if(is.character(prior_beta) & prior_beta != "diag" &
     prior_beta != "Ridge_Zellner") 
    stop("The prior for beta should be 'diag' or 'Ridge_Zellner'.")
  if (prior_beta == "diag"){
    V_tilde <- diag(1+sum(K)) 
    V_tilde[1,1] <- 1e2 * var(y)
    indice_q <- 1
    for(q in 1:Q){
      var_x <- min(apply(x_mult[[q]],2,var))
      V_tilde[(indice_q+1):(indice_q+K[q]),(indice_q+1):(indice_q+K[q])] <- 
        diag(K[q]) * 1e2 * var(y)/var_x
      indice_q <- indice_q+K[q]
    }
    if(is.null(g))      g <- length(y)
    if(is.null(lambda)) lambda <- 5
  }
  if (prior_beta == "Ridge_Zellner"){
    if(is.null(g))       g <- length(y)
    if(is.null(lambda))  lambda <- 5
    if(is.null(V_tilde)) V_tilde <- diag(1+sum(K))
    
  }
  if(is.null(rho)) rho <- 0
  
  if(is.null(conf_expert)){
    conf_expert <- list()[1:Q]
    for(q in 1:Q){
      conf_expert[[q]] <- list()[1:length(beta_s_expert[[q]])]
      for(j in 1:length(beta_s_expert[[q]])){
        conf_expert[[q]][[j]] <- rep(1,length(grids[[q]])) #/ length(beta_s_expert[[q]])  
      }
    }
  } 
  
  if(is.null(beta_s_dist)) beta_s_dist <- "L2"
  L_max <- 0
  beta_s_expert_average <- list()[1:Q]
  
  # beta_s moyen quand la distance est la distance L2
  if(beta_s_dist == "L2"){
    conf_tilde <- list()
    # Define the L2 distance 
    dist <- function(grid,bs1,bs2,g){
      integrate_trapeze(grid, (bs1-bs2)^2 * g  )
    }
    for(q in 1:Q){
      # Compute sum of confidences
      conf_tilde[[q]] <- 0*conf_expert[[q]][[1]]
      for(j in 1:length(beta_s_expert[[q]])){
        conf_tilde[[q]] <- conf_tilde[[q]] + conf_expert[[q]][[j]]
      }
      
      # Compute the average beta_s of the experts
      beta_s_expert_average[[q]] <- rep(0,length(grids[[q]]))
      
      # if( sum(conf_tilde[[q]]) > 0 ){
      #   for(j in 1:length(beta_s_expert[[q]])){
      #     beta_s_expert_average[[q]] <- beta_s_expert_average[[q]] +
      #       beta_s_expert[[q]][[j]] * conf_expert[[q]][[j]] / conf_tilde[[q]]
      #   }
      # }
      
      for(j in 1:length(beta_s_expert[[q]])){
        for( t in 1:length(grids[[q]]) ){
          if( conf_tilde[[q]][t] != 0 ){
            beta_s_expert_average[[q]][t] <- beta_s_expert_average[[q]][t] +
              beta_s_expert[[q]][[j]][t] * conf_expert[[q]][[j]][t] / conf_tilde[[q]][t]
          }
        }
      }
      
      
      # # Compute the maximum of the distance 
      # beta_s_opppose <- beta_s_expert_average[[q]]
      # beta_s_opppose[ beta_s_opppose >  0 ] <- 2
      # beta_s_opppose[ beta_s_opppose <= 0 ] <- 1
      # beta_s_opppose[ beta_s_opppose == 2 ] <- -1
      # 
      # L_max <- L_max + dist(grids[[q]],beta_s_opppose,
      #                       beta_s_expert_average[[q]],conf_tilde)
    }
  }
  
  # beta_s moyen quand la distance est la distance 0-1
  if(beta_s_dist == "0-1"){
    # Define the 0-1 distance 
    dist <- function(grid,bs1,bs2,g){
      integrate_trapeze(grid, as.numeric(bs1 != bs2) * g  )
    }
    
    for(q in 1:Q){
      # Compute sum of confidences
      conf_tilde <- 0*conf_expert[[q]][[1]]
      for(j in 1:length(beta_s_expert[[q]])){
        conf_tilde <- conf_tilde + conf_expert[[q]][[j]]
      }
      
      # Compute the average beta_s of the experts
      beta_s_expert_average[[q]] <- rep(0,length(grids[[q]]))
      for(t in 1:length(grids[[q]])){
        S_1  <- 0
        S_0  <- 0
        S_m1 <- 0
        
        for(j in 1:length(beta_s_expert[[q]])){
          if(beta_s_expert[[q]][[j]][t] == 1){
            S_0  <- S_0  + 1*conf_expert[[q]][[j]][t]
            S_m1 <- S_m1 + 4*conf_expert[[q]][[j]][t]
          }
          if(beta_s_expert[[q]][[j]][t] == 0){
            S_1  <- S_1  + 1*conf_expert[[q]][[j]][t]
            S_m1 <- S_m1 + 1*conf_expert[[q]][[j]][t]
          }
          if(beta_s_expert[[q]][[j]][t] == -1){
            S_1  <- S_1  + 4*conf_expert[[q]][[j]][t]
            S_0  <- S_0  + 1*conf_expert[[q]][[j]][t]
          }
        }
        
        if(S_1 < S_0  && S_1 < S_m1)
          beta_s_expert_average[[q]][t] <- 1
        if(S_0 < S_1  && S_0 < S_m1)
          beta_s_expert_average[[q]][t] <- 0
        if(S_m1 < S_0 && S_m1 < S_1)
          beta_s_expert_average[[q]][t] <- -1
      }
      
      # Compute the maximum of the distance 
      beta_s_opppose <- beta_s_expert_average[[q]]
      beta_s_opppose[ beta_s_opppose == 1 ] <- 2
      beta_s_opppose[ beta_s_opppose == -1 ] <- 1
      beta_s_opppose[ beta_s_opppose == 2 ] <- -1
      beta_s_opppose[ beta_s_opppose == 0 ] <- 1
      
      L_max <- L_max + dist(grids[[q]],beta_s_opppose,
                            beta_s_expert_average[[q]],conf_tilde)
    }
  }
  
  if(is.null(iter_tau)) iter_tau <- 1e3
  if(is.null(tau_vec)){
    tau_vec <- rep(0,2)
    choose_tau_vec <- "data_driven"
  }else{
    choose_tau_vec <- "user"
  }
  if(length(tau_vec) == 1) 
    choose_tau_vec <- "user"
  
  if(is.null(lambda_tau)) lambda_tau <- 0
  if(is.null(is_tau))     is_tau     <- "random"
  if(is.null(tau_vec_size)) tau_vec_size <- 10
  
  # Perfome the Metropolis within Gibbs and return the result.
  res <- Bliss_MWG_multiple_cpp(Q,y,x_mult,iter,grids,K,l_max,eta_tilde,a,b,
                                phi_m,phi_l,prior_beta,g,lambda,V_tilde,
                                tol=sqrt(.Machine$double.eps),basis,
                                beta_s_expert_average,conf_tilde,MH_proposal,
                                rho,iter_tau,tau_vec,choose_tau_vec,
                                is_tau,lambda_tau,tau_vec_size)
  res$beta_s_expert_average <- beta_s_expert_average
  res$conf_tilde <- conf_tilde
  trace_names <- NULL
  for(q in 1:Q){
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("beta",k,"q",q,sep="_"))
    }
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("m",k,"q",q,sep="_"))
    }
    for(k in 1:K[q]){
      trace_names <- c(trace_names,paste("l",k,"q",q,sep="_"))
    }
  }
  colnames(res$trace) <- c(trace_names ,"mu","sigma_sq","tau","alpha",
                           "accepted","SSE_diff","dist_diff","d_beta_s","RSS" )
  return(res)
}

################################# ----
# Bliss_Gibbs_Sampler_multiple
################################# ----
# description : Perform the Gibbs Sampler to sample from the posterior 
#               distribution of the Bliss model.
# value : a list containing :
#             trace : a matrix. Each row is an iteration of the Gibbs Sampler. 
#                The colums are: beta1star,m1,l1,...,betaQstar,mQ,lQ,mu,sigma2.
#             param : a list containing a, b, V_tilde, K, eta_tilde, l_max, grids and 
#                         all_intervals (see the function potential_intervals).
# argument 1 : data, a list containing 1) the number of covariates, 2) the functions x_qi(t) observed at
#              grids of time points and 3) the outcome values y_i.
# argument 2 : param, a list containing :
#              iter : an integer, the number of iterations of the
#                 Gibbs sampler algorithm.
#              grids : a list of numerical vectors, the qth vector is the grid of observation points for the qth covariate.
#              K : a vector of integers, hyperparameters of the Bliss model, the number of intervals for the Q covariates.
#              l_max : a vector of integers, hyperparameters of the Bliss model. Beware, give the index 
#                 corresponding to the value of the maximum width of the intervals (optional)
#              basis : a vector of characters among : "uniform" (default), 
#                 "epanechnikov", "gauss" and "triangular". This indicates the 
#                 shapes of the Q coefficient functions on the intervals. (optional)
#              eta_tilde : a numerical vector of length (1+sum(K)), hyperparameter of the Bliss model. 
#                 By default, eta is the vector (0,...,0). (optional)
#              V_tilde : a matrix of dimension (1+sum(K))*(1+sum(K)), hyperparameter 
#                 of the Bliss model. (optional) (nasty code)
#              a : a nonnegative value, hyperparameter of the Bliss model. By 
#                 default, a = 0.1. (optional)
#              b : a nonnegative value, hyperparameter of the Bliss model. By 
#                 default, b = 0.1. (optional)
#              g : hyperparameter of the Bliss model,  a nonnegative value, 
#                the coefficient of the Zellner prior.
#              phi_m : a list of numerical vectors. The priors of the mq, q=1,...,Q. If not specified, a uniform 
#                 distribution is used for phi_m[[q]].
#              phi_l : a list of numerical vectors. The priors of the lq, q=1,...,Q. If not specified, a uniform 
#                 distribution is used for phi_l[[q]].
#              phi_l_mean : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the mean of the Gamma prior of lq.
#              phi_l_sd : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the standard deviation of the Gamma prior of lq.
#              prior_beta : a character string, which indicates the prior on the
#                beta_star[q]. The possible values are : 
#                1) "diag" (default) for a diagonal matrix prior,
#                2) "Ridge_Zellner" for the Ridge Zellner prior,
#              G : a matrix, the precision matrix (inverse of covariance) of 
#                 y | \mu, \beta, \sigma, m, l.
#              elicited : a boolean value, indicates if there is some 
#                 elicited data in y nor x.
#              conf : a numerical value, the trust that you have in the 
#                 expert's opinions. 
#              n0 : a positive integer, the number of non-elicited observations.

# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# 
# res_Bliss_Gibbs_Sampler <- Bliss_Gibbs_Sampler(data,param)
# K       <- param$K 
# theta_1 <- res_Bliss_Gibbs_Sampler$trace[1,]
# theta_1

# beta_star <- theta_1[1+1:K]
# m         <- theta_1[1:K+K+2]
# l         <- theta_1[1:K+2*K+2]
# beta_1    <- beta_build_cpp(beta_star,m,l,param$grid,param$p,K,"uniform")
# plot(param$grid,beta_1,type="l")
Bliss_Gibbs_Sampler_multiple <- function(data,param){
  # Initialize
  x_mult <- data$x_mult
  y <- data$y
  Q <- length(x_mult)
  # load objects
  iter   <- param[["iter"]]
  grids   <- param$grids
  grids_l <- list()
  for (q in 1:Q){
    grids_l[[q]] <- (grids[[q]] - grids[[q]][1])[-1]
  }
  K      <- param$K
  p      <- param[["p"]]
  
  # load optional objects
  l_max   <- param[["l_max"]]
  eta_tilde     <- param$eta_tilde
  V_tilde <- param[["V_tilde"]]
  lambda  <- param[["lambda"]]
  g       <- param[["g"]]
  a       <- param$a
  b       <- param[["b"]]
  phi_m   <- param[["phi_m"]]
  phi_l   <- param[["phi_l"]]
  basis   <- param[["basis"]]
  phi_l_mean <- param[["phi_l_mean"]]
  phi_l_sd   <- param[["phi_l_sd"]]
  prior_beta <- param$prior_beta
  G          <- param[["G"]]
  elicited   <- param[["elicited"]]
  conf       <- param[["conf"]]
  n0         <- param[["n0"]]
  
  # Initialize the necessary unspecified objects 
  p <- numeric()
  for(q in 1:Q){
    p[q] <- length(grids[[q]])
  }
  if(is.null(prior_beta))
    stop("Please specify a value for the vector prior_beta.")
  if(is.null(K)) stop("Please specify a value for the vector K.")
  if(!is.null(K)){
    for (q in 1:Q){
      if(is.na(K[q])) 
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(eta_tilde))     eta_tilde <- rep(0,1+sum(K))
  if(!is.null(eta_tilde)){
    for (q in 1:Q){
      if(is.na(eta_tilde[q])) 
        stop("Please specify a value for all components of the vector K.")
    }
  }
  if(is.null(a))       a       <- 1e-1
  if(is.null(b))       b       <- 1e-1
  if(is.null(basis)){
    basis <- character()
    for (q in 1:Q){
      basis[q] <- "uniform"
    }
  }  
  if(!is.null(basis)){
    for (q in 1:Q){
      if(is.na(basis[q])){basis[q] <- "uniform"}
    }
  }      
  if (is.null(phi_l_mean)) phi_l_mean <- rep(NA,Q)
  if (is.null(phi_l_sd)) phi_l_sd <- rep(NA,Q)
  for(q in 1:Q){
    if(!is.null(phi_l[[q]]) && is.character(phi_l[[q]]) && 
       phi_l[[q]] != "Gamma") 
      stop("The qth component of phi_l should be a numeric vector or 'Gamma'.")
    if(!is.null(phi_l[[q]]) && is.character(phi_l[[q]]) &&
       phi_l[[q]] == "Gamma"){
      
      if(is.na(phi_l_mean[q])) phi_l_mean[q] <- 
          diff(range(grids[[q]]))/5 + grids[[q]][1]
      if(is.na(phi_l_sd[q]))   phi_l_sd[q]   <- 
          diff(range(grids[[q]]))/5 
      phi_l[[q]] <- prior_l(phi_l_mean[q]/K[q],phi_l_sd[q]/K[q],grids_l[[q]])   
    }
  }
  if(is.null(phi_l)){
    phi_l <- list()
    for (q in 1:Q){
      if(is.null(l_max)) l_max <- floor(p/5)
      if(!is.null(l_max) & is.na(l_max[q])){l_max[q] <- floor(p[q]/5)}
      phi_l[[q]] <- rep(1/l_max[q],l_max[q])
    }  
  }
  
  for (q in 1:Q){
    l_max[q] <- length(phi_l[[q]])
  }
  if(!is.null(phi_m)){
    for (q in 1:Q){
      if(is.na(phi_m[[q]])){phi_m[[q]] <- rep(1/p[q],p[q])}
    }    
  }  
  if(is.null(phi_m)){
    phi_m <- list()
    for (q in 1:Q){
      phi_m[[q]] <- rep(1/p[q],p[q])
    }  
  }
  if(is.null(prior_beta)) prior_beta <- "diag"
  if(is.character(prior_beta) & prior_beta != "diag" &
     prior_beta != "Ridge_Zellner") 
    stop("The prior for beta should be 'diag' or 'Ridge_Zellner'.")
  if (prior_beta == "diag"){
    V_tilde <- diag(1+sum(K)) 
    V_tilde[1,1] <- 1e2 * var(y)
    indice_q <- 1
    for(q in 1:Q){
      var_x <- min(apply(x_mult[[q]],2,var))
      V_tilde[(indice_q+1):(indice_q+K[q]),(indice_q+1):(indice_q+K[q])] <- 
        diag(K[q]) * 1e2 * var(y)/var_x
      indice_q <- indice_q+K[q]
    }
    if(is.null(g))      g <- length(y)
    if(is.null(lambda)) lambda <- 5
  }
  if (prior_beta == "Ridge_Zellner"){
    if(is.null(g))       g <- length(y)
    if(is.null(lambda))  lambda <- 5
    if(is.null(V_tilde)) V_tilde <- diag(1+sum(K))
    #     V_tilde[1,1] <- 1e2 * var(y)
    #     V_tilde[-1,-1] <- diag(sum(K))  * 5
    
  }
  if(is.null(G)) G <- diag(length(y))
  if(is.null(elicited)) elicited <- FALSE
  if(elicited) conf <- compute_conf(n0,G) else conf <- 0
  
  # Perfome the Gibbs Sampler and return the result.
  res <- Bliss_Gibbs_Sampler_multiple_cpp(Q,y,x_mult,iter,grids,K,l_max,
                                          eta_tilde,a,b,phi_m,phi_l,
                                          prior_beta,g,lambda,V_tilde,
                                          tol=sqrt(.Machine$double.eps),
                                          basis,G,conf)
  
  return(res)
}

################################# ----
# Bliss_Simulated_Annealing
################################# ----
# description : Perform the Simulated Annealing algorithm to determine the
#               minimum of the posterior expectation loss, i.e. the 
#               Bliss estimate.
# value : a list containing :
#             Bliss_estimate : a numerical vector, corresponding the Bliss estimate.
#             posterior_expe : a numerical vector, which is posterior 
#                 expectation of beta(t), for t fixed and for each t in the 
#                 grid ot time points.
#             posterior_var : a numerical vector, which is posterior 
#                 expectation of beta(t), for t fixed and for each t in the 
#                 grid ot time points.
#             trace : a matrix. Each row is an iteration of the Simulated 
#                 Annealing algorithm.
#             argmin : an integer, which is the index of the iteration 
#                 minimizing the loss.
# argument 1 : res.Gibbs_Sampler, a list. This is the result of the 
#              Bliss_Gibbs_Sampler function.
# argument 2 : param, a list containing :
#                 burnin : an integer, the number of iterations to drop from 
#                     the Gibbs sample. (optional)
#                 iter : an integer, the number of iteration of 
#                     the Gibbs Sampler algorithm.
#                 iter_sann : an integer, the number of iteration of 
#                     the Simulated Annealing algorithm. (optional)
#                 Temp : a non negative value, the initial temperature for the 
#                     cooling function of the Simulated Annealing. (optional)
#                 k_max : an integer, the maximum number of intervals.
#                 l_max : an integer, the maximum value for the parameter l. (optional)
#                 basis : a character vectors, used to compute the coefficient 
#                     function beta, see the function beta_build. (optional)
# argument 3 : scale_ml, a matrix returned by the function "Bliss_Gibbs_Sampler".
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# 
# res.Simulated_Annealing <- Bliss_Simulated_Annealing(beta_functions,param)
# 
# ylim <- range(c(res.Simulated_Annealing$Bliss_estimate,
#                 res.Simulated_Annealing$posterior_expe))
# plot(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",ylim=ylim)
# lines(param$grid,res.Simulated_Annealing$posterior_expe,lty=2)
Bliss_Simulated_Annealing <- function(beta_functions,param,scale_ml){  
  # Initialize
  grid <- param$grid
  iter <- param[["iter"]]
  p    <- param$p
  
  # load optional objects
  Temp      <- param$Temp
  k_max     <- param$k_max
  iter_sann <- param[["iter_sann"]]
  burnin    <- param[["burnin"]]
  l_max     <- param[["l_max_sann"]]
  basis     <- param[["basis"]]
  
  # Initialize the necessary unspecified objects   
  if(is.null(Temp))      Temp      <- 1000
  if(is.null(k_max))     k_max     <- 5
  if(is.null(iter_sann)) iter_sann <- 1e5
  if(is.null(burnin))    burnin    <- floor(iter/5)
  if(is.null(l_max))     l_max     <- floor(p/5)  
  if(is.null(basis))     basis     <- "uniform"
  
  # Check if the burnin value is correct.
  if(iter <= burnin+1){
    burnin <- floor(iter/5)
    cat("Burnin is too large. New burnin : ",burnin,"\n")
  } 
  
  # dm is related to the random walk. A new m_k' is chosen in [m_k - dm , m_k + dm].
  dm <- floor(p/5)+1
  # dl is related to the random walk. A new l_k' is chosen in [l_k - dl , l_k + dl].
  dl <- floor(l_max/2)+1
  
  # Compute the Simulated Annealing algorithm (3 times to find a suitable value
  # for the initial temperature)
  res_sann_list <- list()
  
  # first time
  cat("First Simulated Annealing. \n")
  res_sann_list[[1]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                      grid,burnin,Temp,k_max,
                                                      l_max,dm,dl,p,basis,
                                                      scale_ml) 
  # Revise the initial temperature
  Temp <- min(abs(range(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])
                  - median(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)])))
  # second time
  cat("Second Simulated Annealing. \n")
  res_sann_list[[2]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                      grid,burnin,Temp,k_max,
                                                      l_max,dm,dl,p,basis,
                                                      scale_ml) 
  # Revise the initial temperature
  Temp <- min(abs(range(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])
                  - median(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)])))
  # third time
  cat("Third Simulated Annealing. \n")
  res_sann_list[[3]] <- Bliss_Simulated_Annealing_cpp(iter_sann,beta_functions,
                                                      grid,burnin,Temp,k_max,
                                                      l_max,dm,dl,p,basis,
                                                      scale_ml) 
  # Comparison and selection  
  mins      <- c(min(res_sann_list[[1]]$trace[,ncol(res_sann_list[[1]]$trace)]),
                 min(res_sann_list[[2]]$trace[,ncol(res_sann_list[[2]]$trace)]),
                 min(res_sann_list[[3]]$trace[,ncol(res_sann_list[[3]]$trace)]))
  index <- which(mins==min(mins))
  res_sann <- res_sann_list[[index]]
  
  # Determine the estimate  
  index <- which(res_sann$trace[,ncol(res_sann$trace)] %in% 
                   min(res_sann$trace[,ncol(res_sann$trace)]))[1]
  
  argmin <- res_sann$trace[index,]
  res_k  <- argmin[length(argmin)-1]
  
  # Compute the estimate
  beta_star <- argmin[1:res_k]
  m         <- argmin[1:res_k+  k_max]
  l         <- argmin[1:res_k+2*k_max]
  k         <- argmin[3*k_max+2]  
  estimate  <- beta_build_cpp(beta_star,m,l,grid,p,k,basis,scale_ml)
  
  return(list(Bliss_estimate = estimate,
              posterior_expe = res_sann$posterior_expe,
              posterior_var  = res_sann$posterior_var,
              trace          = res_sann$trace,
              argmin         = argmin))
}
