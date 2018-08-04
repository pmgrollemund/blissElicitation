
#########################################################
#                                                       #
#              Bliss method : main script               #
#                                                       # 
#########################################################
#### Install the required packages ----
if(!suppressMessages(suppressWarnings(require(Rcpp)))){
  cat("Install the 'Rcpp' package. \n")
  install.packages("Rcpp", repos = "http://ftp.igh.cnrs.fr/pub/CRAN/")
}
if(!suppressMessages(suppressWarnings(require(RcppArmadillo)))){
  cat("Install the 'RcppArmadillo' package. \n")
  install.packages("RcppArmadillo",repos = "http://ftp.igh.cnrs.fr/pub/CRAN/")
}
if(!suppressMessages(suppressWarnings(require(MASS)))){
  cat("Install the 'MASS' package. \n")
  install.packages("MASS", repos = "http://ftp.igh.cnrs.fr/pub/CRAN/")
}
if(!suppressMessages(suppressWarnings(require(ggplot2)))){
  cat("Install the 'ggplot2' package. \n")
  install.packages("ggplot2", repos = "http://ftp.igh.cnrs.fr/pub/CRAN/")
}
if(!suppressMessages(suppressWarnings(require(RColorBrewer)))){
  cat("Install the 'RColorBrewer' package. \n")
  install.packages("RColorBrewer", repos = "http://ftp.igh.cnrs.fr/pub/CRAN/")
}
#### Load the packages and the R scripts ----
# Rccp functions
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(ggplot2)
library(RColorBrewer)
# R script
sourceCpp("./Bliss_eli_multiple_rcpp.cpp")
source("./Bliss_eli_algorithms.R")
source("./Bliss_eli_basic_Functions.R")
source("./Bliss_simulate.R")

################################# ----
# Bliss_multiple
################################# ----
# description : Perform the Bliss method to obtain an interpretable estimate 
#               of the Q coefficient functions of the functional linear 
#               regression model with Q functional covariates.
# value : return a list containing :
#             Bliss_estimate : a list of numerical vectors, the Bliss estimates : we obtain Q estimates, one for each covariable.
#             posterior_density_estimate : a list of lists containing the estimates of the
#                 posterior density. (obtained with the kde2d function). Firts level of the list: Q lists, for the Q covariates. Second
#                 level: for each covariate we have a list of three components, that is the  x and y coordinates of the grid points, and z
#                 a matrix of the estimated density: rows correspond to the value of x, columns to the value of y. 
#             beta_functions : a list of matrices. For the qth covariate, beta_functions[[q]] is a matrix where each row is a
#                 function beta_qi associated to the iteration i of the Gibbs sampler.
#             res.Simulated_Annealing : a list of lists: res.Simulated_Annealing[[q]] is the result of the 
#                 function Bliss_Simulated_Annealing applied for the qth covariate.
#             res.Gibbs_Sampler : a list of lists: res.Gibbs_Sampler[[q]] is the the result of the function 
#                 Bliss_Gibbs_Sampler for the qth covariate.
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
#              phi_l : a list of numerical vectors and/or characters. The priors of the lq, q=1,...,Q. If not specified, a uniform 
#                 distribution is used for phi_l[[q]].
#              phi_l_mean : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the mean of the Gamma prior of lq.
#              phi_l_sd : a Q vector of numerical values. if "phi_l[[q]]" is "Gamma", phi_l_mean[q]
#                corresponds to the standard deviation of the Gamma prior of lq.
#              prior_beta : a character string, which indicates the prior on the
#                beta_star[q]. The possible values are : 
#                1) "diag" (default) for a diagonal matrix prior,
#                2) "Ridge_Zellner" for the Ridge Zellner prior,
#              burnin : an integer, the number of iteration to drop of the 
#                 Gibbs sampler. (optional)
#              thin : an integer, used to thin the Gibbs sample to compute an 
#                 estimate of the posterior density of beta(t). (optional)
#              lims.kde : a Qx2 matrix. lims.kde[q,] corresponds to the parameter (yl,yu) for the representation 
#                 of the posterior density of the qth covariable. There are limits of the y-axis for the function kde2d. (optional)
#              n : an integer, the number of grid points in each direction, for the kde2d function. 
#                 See function kde2d. (optional)
#              iter_sann : an integer, the number of iteration of the 
#                 Simulated Annealing. (optional)
#              Temp : a vector of nonnegative values, the Q initial temperatures for the 
#                 cooling function of the Q Simulated Annealings. (optional)
#              k_max : a vector of integers, k_max[q] is the maximum number of intervals for the  
#                 function beta_q(t) at each iteration. (optional)
#              cols : a vector of colors for the function image. 
#                 Only if plot=TRUE. (optional)
#              new_grids : a list of Q numerical vectors. If new_grids is not NULL, the 
#                 coefficient functions beta_q(t) at each iteration are computed on these grids 
#                 (only) to plot a graphical representation of the posterior
#                  distribution.
#              h1 : a vector of numerical values which are the Q bandwidths of the kernel 
#                 density estimation for the t-axis (optional). h1[q] the bandwidth of the kernel density estimation for the qth covariate.
#              ylim : a Qx2 matrix, the qth line gives the limits for the y-axis for the plotting function 
#                 "image_Bliss" for the qth covariate.
#              main : a vector of characters, main[q] is for the plotting function 
#                 "image_Bliss" for the qth covariate.
#              n_chains : number of chains to do in the Gibbs sampler.
#              G : a matrix, the precision matrix (inverse of covariance) of 
#                 y | \mu, \beta, \sigma, m, l.
#              elicited : a boolean value, indicates if there is some 
#                 elicited data in y nor x.
#              conf : a numerical value, the trust that you have in the 
#                 expert's opinions. 
#              n0 : a positive integer, the number of non-elicited observations.
# argument 3 : plot, a logical value. If it is TRUE, the Bliss estimates and 
#                 the estimates of the posterior density are plotted. (optional)
# example : see the script Bliss_demo.
Bliss_multiple   <- function(data,param,plot=FALSE){
  # preprocessing
  Q <- length(data$x_mult)
  p  <- numeric()
  for (q in 1:Q){
    # data$x_mult[[q]]   <- apply(data$x_mult[[q]],2,
    #                             function(vect) vect - mean(vect))
    p[q] <- length(param$grids[[q]])
  }
  param$p <- p
    
  # How many chains i have to do ?
  n_chains <- param[["n_chains"]]
  if(is.null(n_chains)){
    n_chains <- 1
    param[["n_chains"]] <- n_chains
  } 
  
  # Initialize the list "chains"
  chains <- list()
  # Each chain :
  for(j in 1:n_chains){
    cat("Chain ",j,": \n",sep="")
    # Initialize the list "chains[[j]]"
    chains[[j]] <- list()
    # Execute the Gibbs Sampler algorithm to obtain a posterior sample
    chains[[j]]$res.Gibbs_Sampler <- Bliss_multiple_WL(data,param)
    
    # Compute the functions beta_i for each iteration i of the Gibbs sample.
    chains[[j]]$beta_functions <- 
      compute_beta_functions_mult(chains[[j]]$res.Gibbs_Sampler,param)
    
    # Estimate the density of the posterior sample of functions beta_i
    chains[[j]]$posterior_density_estimate <- list()
    for(q in 1:Q){
      param_density <- list(grid= param$grids[[q]],
                            iter= param$iter,
                            p   = param[["p"]][q],
                            n        = param[["n"]],
                            thin     = param$thin,
                            burnin   = param[["burnin"]],
                            lims.kde = param$lims.kde[[q]],
                            h1       = param$h1,
                            new_grid = param[["new_grid"]]
                            )
      chains[[j]]$posterior_density_estimate[[q]] <- 
        density_estimation(chains[[j]]$beta_functions[[q]],param_density)
    }
    
  }
  
  # Choose a chain 
  j <- sample(n_chains,1)
  
  res.Gibbs_Sampler          <- chains[[j]]$res.Gibbs_Sampler
  beta_functions             <- chains[[j]]$beta_functions
  posterior_density_estimate <- chains[[j]]$posterior_density_estimate
  
  # Execute the Simulated Annealing algorithm to obtain an estimate
  res.Simulated_Annealing <- list()
  for(q in 1:Q){
    beta_functions_tmp <- beta_functions[[q]]
    scale_tmp <- res.Gibbs_Sampler$param$scale_ml[[q]]
    param_sann <- list( grid = param$grids[[q]],
                        iter = param[["iter"]],
                        p    = param$p[q],
                        Temp = param$Temp[q],
                        k_max = param$k_max[q],
                        iter_sann = param[["iter_sann"]],
                        burnin    = param[["burnin"]],
                        l_max     = param[["l_max_sann"]][q],
                        basis     = param[["basis"]][q])
    
    res.Simulated_Annealing[[q]] <- 
      Bliss_Simulated_Annealing(beta_functions_tmp,param_sann,scale_tmp)
  }
  
  # Plot the Bliss estimate and the posterior density estimate
  if(plot){
    for(q in 1:Q){
      image_Bliss(posterior_density_estimate[[q]],param)
      lines(param$grid,res.Simulated_Annealing[[q]]$Bliss_estimate,
            type="s",lwd=2)
      #lines(param$grid,res.Simulated_Annealing$posterior_expe,type="l",lty=2)
    }
  }
   
  # Do not return the list "chains" if n_chains is 1.
  if(n_chains == 1) chains <- NULL
  
  # The object to return
  Bliss_estimate <- list()
  for(q in 1:Q){
    Bliss_estimate[[q]] <- res.Simulated_Annealing[[q]]$Bliss_estimate
  }
  res <- list(Bliss_estimate             = Bliss_estimate,
              posterior_density_estimate = posterior_density_estimate,
              beta_functions             = beta_functions,
              res.Simulated_Annealing    = res.Simulated_Annealing,
              res.Gibbs_Sampler          = res.Gibbs_Sampler,
              chains                     = chains)
    
  return(res)
}