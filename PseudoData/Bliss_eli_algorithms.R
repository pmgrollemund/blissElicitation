#########################################################
#                                                       #
#    Bliss method : MCMC and optimization algorithms    #
#                                                       # 
#########################################################
################################# ----
# Bliss_multiple_WL
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
Bliss_multiple_WL <- function(data,param){
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
  n_y    <- length(y) 
  
  # load optional objects
  l_max   <- param[["l_max"]]
  V_tilde <- param[["V_tilde"]]
  lambda  <- param[["lambda"]]
  phi_m   <- param[["phi_m"]]
  phi_l   <- param[["phi_l"]]
  basis   <- param[["basis"]]
  phi_l_mean <- param[["phi_l_mean"]]
  phi_l_sd   <- param[["phi_l_sd"]]
  
  conf       <- param[["conf"]]
  rho_mat    <- param[["rho_mat"]]
  y_expert   <- param[["y_expert"]]
  x_expert   <- param[["x_expert"]]
  posterior  <- param[["posterior"]]
  
  
  # Initialize the necessary unspecified objects 
  p <- numeric()
  for(q in 1:Q){
    p[q] <- length(grids[[q]])
  }
  if(is.null(K)) stop("Please specify a value for the vector K.")
  if(!is.null(K)){
    for (q in 1:Q){
      if(is.na(K[q])) 
        stop("Please specify a value for all components of the vector K.")
    }
  }
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
  
  # For the matrix of the Rifge Zellner prior
  g <- length(y)
  V_tilde <- diag(1+sum(K))
  if(is.null(lambda))  lambda <- 5
  
  ##### part : expert's information
  if(is.null(posterior)) posterior <- TRUE
  if(is.null(y_expert)){
    y_expert <- y
    x_expert_final <- x_mult
    
    weights <- list(rep(1,n_y),rep(0,n_y))                                      ## changer pour mettre qu'une donnees ? ou rien mettre ?
    res_weight <- weights
    average_weights <- c(1,0)
    posterior <- TRUE
    n_e <- c(n_y,n_y)
  }else{
    # Compute the weights
    res_weight <- compute_weight(conf,rho_mat,n_y)  
    weights    <- res_weight$weights 
    average_weights <- NULL
    # Centrage des y_expert
    for(e in 1:length(weights) ){ 
      if(!is.list(y_expert)) 
        y_expert <- list(y_expert)
      if(!is.list(weights)) 
        weights <- list(weights)
        
      average_weights <- c(average_weights,mean(weights[[e]])) 
    }
    
    # if you want the expert's prior or posterior #### faire plus propre pour avoir le prior expert
    if(posterior){
      if(is.list(y_expert))
        n_e <- c(n_y, sapply(y_expert,length)) else
          n_e <- c(n_y, length(y_expert))
      weights <- c( list(rep(1,n_y)), weights)
      average_weights <- c(1,average_weights)
    }else{
      if(is.list(y_expert))
        n_e <- sapply(y_expert,length) else
          n_e <- length(y_expert)
    }
    
    # Centrage des x_experts
    # for(q in 1:Q){
    #   for(e in 1:length(y_expert))
    #     x_expert[[q]][[e]] <- apply(x_expert[[q]][[e]],2,
    #                                 function(v) v - mean(v))                    
    # }
    
    x_expert_final <- list()
    length(x_expert_final) <- Q
    for(q in 1:Q){
      for(e in 1:length(y_expert)){
        x_expert_final[[q]] <- rbind(x_expert_final[[q]],x_expert[[q]][[e]])
        x_expert_final[[q]] <- apply(x_expert_final[[q]],2,
                                     function(v) v - mean(v))
      }
    }
  }
  
  y_expert <- unlist(y_expert)
  # Perfome the Gibbs Sampler and return the result.
  res <- Bliss_multiple_WL_cpp(Q,iter, grids, posterior, 
                                
                                y, x_mult, unlist(y_expert), x_expert_final,
                                n_e,weights, average_weights, 
                                
                                K, l_max, phi_m, phi_l, lambda, V_tilde, 
                                tol=sqrt(.Machine$double.eps),basis)
  
  
  res$res_weight <- res_weight
  res$n_e <- n_e
  res$y_bar <- c(y,unlist(y_expert))
  res$x_bar <- x_expert_final
  res$average_weights <- average_weights
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
  colnames(res$trace) <- c(trace_names ,"mu","sigma_sq")
  
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
