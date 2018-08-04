#########################################################
#                                                       #
#           Bliss method : basic Functions              #
#                                                       # 
#########################################################
################################# ----
# plot_beta_s_posterior
################################# ----
# description : 
# value : 
# argument 1 : 
# argument 2 : 
# argument 3 : 
# details  : 
# example :
plot_beta_s_posterior <- function(res_bliss_mult,param,which=1,...){
  # Initialize some objects
  Q <- length(res_bliss_mult$beta_s)
  posterior_mean_beta_s <- list()
  posterior_0 <- list()
  posterior_p <- list()
  posterior_n <- list()
  col_assignation_0  <- list()
  col_assignation_p  <- list()
  col_assignation_n  <- list()
  
  # load graphical parameter
  have_to_plot <- param$plot
  cols <- param$cols
  
  if(is.null(cols)) cols <- rev(heat.colors(1e2))
  if(is.null(have_to_plot)) have_to_plot <- TRUE 
  
  breaks <- seq(0,1,le=length(cols)+1)
  
  # loop on number of covariates
  for(q in 1:Q){
    # load the grid of the q^th covariate
    grid <- param$grids[[q]]
    
    # compute the posterior frequencies
    posterior_0[[q]] <- apply(res_bliss_mult$beta_s[[q]],2,function(v)
      sum(v==0)/length(v))
    posterior_p[[q]] <- apply(res_bliss_mult$beta_s[[q]],2,function(v)
      sum(v==1)/length(v))
    posterior_n[[q]] <- apply(res_bliss_mult$beta_s[[q]],2,function(v)
      sum(v==-1)/length(v))
    
    # Compute the poseterior mean
    posterior_mean_beta_s[[q]] <- apply(res_bliss_mult$beta_s[[q]],2,mean)
    
    # assign colors to posterior frequencies
    col_assignation_0 <- rep(0,length(posterior_0[[q]]))
    col_assignation_p <- rep(0,length(posterior_p[[q]]))
    col_assignation_n <- rep(0,length(posterior_n[[q]]))
    for(i in 1:length(posterior_0[[q]])){
      col_assignation_0[i] <- in_interval(posterior_0[[q]][i],breaks)
      col_assignation_p[i] <- in_interval(posterior_p[[q]][i],breaks)
      col_assignation_n[i] <- in_interval(posterior_n[[q]][i],breaks)
    }
    
    if(have_to_plot==T){
      if(which==1){
        plot(grid,param$grids[[1]],type="n",ylim= c(-1,1),xlab="",ylab="",...)
        points(grid,rep(0,length(posterior_0[[q]])),col=cols[col_assignation_0 ],pch=15)
        points(grid,rep(1,length(posterior_p[[q]])),col=cols[col_assignation_p ],pch=15)
        points(grid,rep(-1,length(posterior_n[[q]])),col=cols[col_assignation_n ],pch=15)
        lines(grid,posterior_mean_beta_s[[q]])
      }
      if(which==2){
        plot(grid,posterior_0[[q]],type="l",xlab="",ylab="",ylim=c(0,1),...)
        lines(grid,posterior_p[[q]],col=2,...)
        lines(grid,posterior_n[[q]],col=3,...)
      }
      if(which==3){
        lines(grid,posterior_0[[q]],...)
        lines(grid,posterior_p[[q]],col=2,...)
        lines(grid,posterior_n[[q]],col=3,...)
      }
      if(which==4){
        plot(grid,posterior_0[[q]],type="n",xlab="",ylab="",ylim=c(0,1),...)
        lines(grid,posterior_p[[q]],col=2,type="n",...)
        lines(grid,posterior_n[[q]],col=3,type="n",...)
        points(grid,posterior_0[[q]],col=1,pch="0",...)
        points(grid,posterior_p[[q]],col=2,pch="p",...)
        points(grid,posterior_n[[q]],col=3,pch="n",...)
      }
    }
  }
  res <- list(posteriors=list(posterior_0,
                              posterior_p,
                              posterior_n),
              col_assignation = list(col_assignation_0,
                                     col_assignation_p,
                                     col_assignation_n),
              posterior_mean_beta_s = posterior_mean_beta_s,
              breaks = breaks,
              cols = cols)
  if(have_to_plot == FALSE)
    return(res)
}
in_interval <- function(value,breaks){
  if(value == 0) return(1) else return( tail(which(value > breaks),1) )
}
################################# ----
# compute_beta_s_mult
################################# ----
# description : 
# value : 
# argument 1 : 
# argument 2 : 
# argument 3 : 
# details  : 
# example :
compute_beta_s_mult <- function(beta_functions_list){
  cat("Compute the functions beta_s_i. \n")
  Q <- length(beta_functions_list)
  beta_s_list <- list()
  length(beta_s_list) <- Q
  for(q in 1:Q){
    beta_functions <- beta_functions_list[[q]]
    beta_s_list[[q]] <- compute_beta_s(beta_functions)
  }
  
  return(beta_s_list)
}

compute_beta_s <- function(beta_functions){
  n <- nrow(beta_functions)
  beta_s <- matrix(0,n,ncol(beta_functions))
  
  for(i in 1:n){
    beta_s[i,] <- as.numeric( beta_functions[i,] > 0) - 
      as.numeric( beta_functions[i,] < 0)
  }
  
  return(beta_s)
}
################################# ----
# compute_conf
################################# ----
# description : 
# value : 
# argument 1 : 
# argument 2 : 
# argument 3 : 
# details  : 
# example :
compute_conf <- function(n0,G){
  sub_G <- G[ -(1:n0) , -(1:n0) ]
  res   <- sum( sub_G ) / n0
  return(res)
}
################################# ----
# finer_grid
################################# ----
# description : Compute a curve on a new finer gird.
# value : a numerical vector, the curve evaluated on the new grid.
# argument 1 : curve, a numerical vector.
# argument 2 : grid, a numerical vector, the former grid.
# argument 3 : grid, a numerical vector.
# details  : This is nasty code.
# example :
# grid <- seq(0,1,l=1e1)
# new_grid <- seq(0,1,l=1e2)
# curve <- 3*grid^2 + sin(grid*2*pi)
# 
# plot(grid,curve,type="o")
# lines(new_grid,finer_grid(curve,grid,new_grid),type="o",col=makeTransparent(2))
finer_grid <- function(curve,grid,new_grid){
  res <- rep(0,length(new_grid))
  for(i in 1:(length(grid)-1)){
    index <- new_grid %between% grid[i:(i+1)]
    res[index] <- curve[i] + 0:(sum(index)-1) / 
      (sum(index)-1) * (curve[i+1]-curve[i])
  }
  return(res)
}
################################# ----
# prior_l
################################# ----
# description : Compute the probability function of the Gamma prior on l.
# value : a numerical vector, which is the prability function on "grid_l".
# argument 1 : m, a positive value, the mean of the Gamma prior.
# argument 2 : sd, a nonnegative value, the standard deviation of the Gamma prior.
# argument 3 : grid_l : a numerical value : the discrete support of the 
#                 parameter l.
# example :
# grid_l <- seq(0,5,l=100)
# f <- prior_l(3,1,grid_l)
# plot(grid_l,f,type="h",xlab="",ylab="")
#
# f <- prior_l(1,0.5,grid_l)
# plot(grid_l,f,type="h",xlab="",ylab="")
#
# f <- prior_l(1,5,grid_l)
# plot(grid_l,f,type="h",xlab="",ylab="")
prior_l <- function(m,s,grid_l){
  # Compute the scale and the rate 
  alpha <- m^2/s^2
  beta  <- m/s^2
  
  # Compute the probability function on the grid
  step <- diff(grid_l)[1] / 2
  probs <- pgamma(grid_l + step ,shape=alpha,rate=beta) -
    pgamma(grid_l - step ,shape=alpha,rate=beta)
  
  return(probs)
}

################################# ----
# makeTransparent
################################# ----
# description : Make transparent some color.
# value : a character string coding for a color.
# argument 1 : ... must be a color. Numerical or character string.
# argument 2 : alpha, a numerical value between 0 and 1, corresponding to the 
#               transparency of the returned color.
# details : Thanks to Ricardo Oliveros-Ramos for the function on the web :
#               http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
# example :
# cols <- makeTransparent(2:6,0.3) 
# res_hist <- hist(rnorm(1e4,0,1),nclass=2e2,border = 0,col = cols[1],
#                  xlim=c(-2,6),xlab="",main = "")
# hist(rnorm(1e4,1,1),border = 0,col = cols[2],add=T,nclass=2e2)
# hist(rnorm(1e4,2,1),border = 0,col = cols[3],add=T,nclass=2e2)
# hist(rnorm(1e4,3,1),border = 0,col = cols[4],add=T,nclass=2e2)
# hist(rnorm(1e4,4,1),border = 0,col = cols[5],add=T,nclass=2e2)
makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

################################# ----
# diagnostics
################################# ----
# description : Perform some diagnostics for the chains resulting of the 
#                 Gibbs Sampler algorithm.
# value : a list containing the diagnostics which can be plotted with the 
#                 function "plot_diagnostics".
# argument 1 : chains, a list containing :
#                 res.Gibbs_Sampler, a list resulting of the function 
#                   Bliss_Gibbs_Sampler.
#                 beta_functions, a matrix which is the result of the 
#                   function compute_beta_functions. 
#               posterior_density_estimate, a list which is the result of the
#                   function density_estimation.
# argument 2 : param, a list containing :
#                 iter : the number of iterations.
#                 burnin : the number of iterations to drop.
#                 K : the number of intervals of the beta_i's functions.
#                 p : the number of time points.
#                 ts : a vector. The sample of time points t_j used for compute
#                   some diagnostics of \beta(t_j). (optional)
#                 l_ts : an integer, the number of time points to sample, if
#                   ts is not specified. (optional)
#                 lag_max : an integer, the maximal lag when compute the 
#                   autocorrelation of the trace.
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# 
# XXXXXXXXXXXXXXXX BULU XXXXXXXXXXXXXXXXXXXXXXXXXXX
diagnostics <- function(chains,param){
  cat("Compute some diagnostics on the posterior sample. \n")
  # load objects
  n_chains  <- param[["n_chains"]]
  burnin    <- param[["burnin"]]
  K         <- param[["K"]]
  p         <- length(param[["grid"]])
  iter      <- param[["iter"]]
  # load optional objects
  ts      <- param[["ts"]]  
  l_ts    <- param[["l_ts"]]  
  lag_max <- param[["lag_max"]]
  # Initialize the necessary unspecified objects 
  if(is.null(l_ts))    l_ts    <- 10
  if(is.null(ts))      ts      <- floor(seq(1,p,l=l_ts))
  if(is.null(lag_max)) lag_max <- min(100,floor(iter/50))
  l_ts <- length(ts)
  # Initialize
  n_class     <- min(1e3 , floor( (iter-burnin+1)/100 ) )
  lags        <- 1:lag_max
  DF_mu       <-  data.frame()
  DF_sigma    <-  data.frame()
  DF_beta     <-  list()
  trace_mu    <- NULL
  trace_sigma <- NULL
  trace_beta  <- NULL
  autocorr_lag_mu    <- NULL
  autocorr_lag_sigma <- NULL
  autocorr_lag_beta  <- NULL
  ylim_mu        <- NULL 
  ylim_sigma     <- NULL  
  ylim_beta      <- rep(NA,2*l_ts)  
  dim(ylim_beta) <- c(2,l_ts)
  hist_mu       <- list()
  hist_sigma    <- list()
  hist_beta     <- list()
  density_mu    <- list()
  density_sigma <- list()
  density_beta  <- list()
  length(DF_beta)       <- l_ts
  length(hist_mu)       <- length(n_chains)
  length(hist_sigma)    <- length(n_chains)
  length(hist_beta)     <- l_ts
  length(density_mu)    <- length(n_chains)
  length(density_sigma) <- length(n_chains)
  length(density_beta)  <- l_ts
  
  #### diagnostic for mu and sigma_sq 
  # Gather the different chains
  for(j in 1:n_chains){
    trace_chain <- chains[[j]]$res.Gibbs_Sampler$trace[-(1:burnin),]
    trace_mu    <- cbind(trace_mu   ,trace_chain[,1])
    trace_sigma <- cbind(trace_sigma,trace_chain[,1+K+1])
    
    DF_mu <- rbind(DF_mu,data.frame(chain = paste("chain_",j,sep=""),
                                    obs   =  trace_chain[,1]))
    DF_sigma <- rbind(DF_sigma,data.frame(chain = paste("chain_",j,sep=""),
                                          obs   =  trace_chain[,1+K+1]))
  }
  
  # Autocorr mu
  for(j in 1:n_chains){
    n_iter <- nrow(trace_mu)
    autocorr_lag_chain <- NULL
    for(l in lags){
      indice     <- 1:(n_iter-l)
      indice_lag <- 1:(n_iter-l) + l 
      
      autocorr_lag_chain <- c(autocorr_lag_chain,
                              cor(trace_mu[indice,j],
                                  trace_mu[indice_lag,j]))
      
    }
    autocorr_lag_mu <- cbind(autocorr_lag_mu,autocorr_lag_chain)
  }
  
  # Autocorr sigma
  for(j in 1:n_chains){
    n_iter <- nrow(trace_sigma)
    autocorr_lag_chain <- NULL
    for(l in lags){
      indice     <- 1:(n_iter-l)
      indice_lag <- 1:(n_iter-l) + l 
      
      autocorr_lag_chain <- c(autocorr_lag_chain,
                              cor(trace_sigma[indice,j],
                                  trace_sigma[indice_lag,j]))
      
    }
    
    autocorr_lag_sigma <- cbind(autocorr_lag_sigma,autocorr_lag_chain)
  }
  
  # Compute the histograms of mu and sigma_sq
  breaks_mu    <- seq(min(DF_mu$obs),max(DF_mu$obs),l=n_class+1)
  breaks_sigma <- seq(min(DF_sigma$obs),max(DF_sigma$obs),l=n_class+1)
  step_mu      <- diff(breaks_mu)[1]
  step_sigma   <- diff(breaks_sigma)[1]
  for(j in 1:n_chains){
    hist_mu[[j]] <- hist(DF_mu$obs[DF_mu$chain==paste("chain_",j,sep="")],
                         plot = FALSE,breaks=breaks_mu)
    density_mu[[j]] <- density(DF_mu$obs[DF_mu$chain==paste("chain_",j,sep="")])
    ylim_mu <- range(ylim_mu, hist_mu[[j]]$density,density_mu[[j]]$y)
    
    hist_sigma[[j]] <- hist(DF_sigma$obs[DF_sigma$chain==paste("chain_",j,sep="")], 
                            plot = FALSE,breaks=breaks_sigma)
    density_sigma[[j]] <- density(DF_sigma$obs[DF_sigma$chain==paste("chain_",j,sep="")])
    ylim_sigma <- range(ylim_sigma, hist_sigma[[j]]$density,density_sigma[[j]]$y)
  }
  
  #### diagnostic for beta(t) (for some t)
  breaks_beta <- rep(0,(n_class+1)*l_ts)
  dim(breaks_beta) <- c(n_class+1, l_ts)
  trace_beta <- rep(0,l_ts*(iter-burnin+1)*n_chains)
  dim(trace_beta) <- c(l_ts,iter-burnin+1,n_chains)
  autocorr_lag_beta      <- rep(0,l_ts*lag_max*n_chains)
  dim(autocorr_lag_beta) <- c(l_ts,n_chains,lag_max)
  for(o in 1:l_ts){
    for(j in 1:n_chains){
      trace_chain <- chains[[j]]$beta_functions[-(1:burnin),ts[o]]
      trace_beta[o,,j] <- trace_chain
      
      DF_beta[[o]] <- rbind(DF_beta[[o]],data.frame(chain = paste("chain_",j,sep=""),
                                                    obs   =  trace_chain))
    }
    
    #  autocorrelation of beta
    for(j in 1:n_chains){
      n_iter <- nrow(as.matrix(trace_beta[o,,]))
      autocorr_lag_chain <- NULL
      for(l in lags){
        indice     <- 1:(n_iter-l)
        indice_lag <- 1:(n_iter-l) + l 
        
        autocorr_lag_chain <- c(autocorr_lag_chain,
                                cor(trace_beta[o,indice,j],
                                    trace_beta[o,indice_lag,j]))
        
      }
      autocorr_lag_beta[o,j,] <- autocorr_lag_chain
    }
    # Compute the histograms and the kernel density estimations 
    breaks_beta[,o] <- seq(min(DF_beta[[o]]$obs),max(DF_beta[[o]]$obs),l=n_class+1)
    step_beta   <- diff(breaks_beta[,o])[1]
    for(j in 1:n_chains){
      hist_beta[[o]][[j]] <- hist(DF_beta[[o]]$obs[DF_beta[[o]]$chain==paste("chain_",j,sep="")],
                                  plot = FALSE,breaks=breaks_beta[,o])
      density_beta[[o]][[j]] <- density(DF_beta[[o]]$obs[DF_beta[[o]]$chain==paste("chain_",j,sep="")])
      ylim_beta[,o] <- range(ylim_beta[,o], hist_beta[[o]][[j]]$density,
                             density_beta[[o]][[j]]$y,na.rm = T)
    }
  }
  
  return( list(DF_mu=DF_mu, trace_mu=trace_mu, autocorr_lag_mu=autocorr_lag_mu,
               hist_mu=hist_mu, density_mu=density_mu, ylim_mu=ylim_mu,
               
               DF_sigma=DF_sigma, trace_sigma=trace_sigma, 
               autocorr_lag_sigma=autocorr_lag_sigma,hist_sigma=hist_sigma,
               density_sigma=density_sigma, ylim_sigma=ylim_sigma,
               
               DF_beta=DF_beta, trace_beta=trace_beta, 
               autocorr_lag_beta=autocorr_lag_beta,hist_beta=hist_beta, 
               density_beta=density_beta, ylim_beta=ylim_beta,
               breaks_beta=breaks_beta,
               
               lags=lags, ts=ts))
}

################################# ----
# plot_diagnostics
################################# ----
# description : Plot the diagnostics computed by the function "diagnostics".
# value : nothing.
# argument 1 : res_diagnostics, a list containing the result of the 
#                 function "diagnostics".
# argument 2 : param, a list containing :
#                 n_chains : the number of chains.
# argument 3 : chain, an integer, corresponding to the index of the chain to 
#                 those diagnositcs have to be plotted. If chain is NULL 
#                 (defaut), the diagnostics of all the chains are plotted on 
#                 the same plots. (optional)
# argument 4 : which_plot, a string character indicating which parameter is of 
#                 interest. The possibe value : 'mu', 'sigma_sq' or 'beta'.  
# argument 5 : time, a numerical value belonging to the "ts" vector (option
#                 of the function "diagnostics"). (optional)
#         
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# 
# XXXXXXXXXXXXXXXX BULU XXXXXXXXXXXXXXXXXXXXXXXXXXX
plot_diagnostics <- function(res_diagnostics,param,chain=NULL,which_plot=NULL,
                             time=NULL){
  if(is.null(which_plot)) stop("Please specify which diagnostics you want : 'mu', 'sigma_sq' or 'beta' ?")
  if(which_plot == "beta"){
    if(is.null(time)){
      time <- sample(res_diagnostics$ts,1)
      cat("You didn't specify which time point have to be taken into account ",
          "for the diagnostics. As an example, here is the diagnostics for time = ",
          time,". \n",sep="")
    }
    if(!( time %in% res_diagnostics$ts)) stop("'time' must belong to 'ts'.")
  } 
  
  # load objects
  if(which_plot == "mu"){
    DF_mu           <- res_diagnostics$DF_mu
    trace_mu        <- res_diagnostics$trace_mu
    autocorr_lag_mu <- res_diagnostics$autocorr_lag_mu
    hist_mu         <- res_diagnostics$hist_mu
    density_mu      <- res_diagnostics$density_mu
    ylim_mu         <- res_diagnostics$ylim_mu
    n_class <- length(hist_mu[[1]]$breaks)-1
  }
  if(which_plot == "sigma_sq"){
    DF_sigma           <- res_diagnostics$DF_sigma
    trace_sigma        <- res_diagnostics$trace_sigma
    autocorr_lag_sigma <- res_diagnostics$autocorr_lag_sigma
    hist_sigma         <- res_diagnostics$hist_sigma
    density_sigma      <- res_diagnostics$density_sigma
    ylim_sigma         <- res_diagnostics$ylim_sigma
    n_class <- length(hist_sigma[[1]]$breaks)-1
  }
  if(which_plot == "beta"){
    DF_beta           <- res_diagnostics$DF_beta
    trace_beta        <- res_diagnostics$trace_beta
    autocorr_lag_beta <- res_diagnostics$autocorr_lag_beta
    hist_beta         <- res_diagnostics$hist_beta
    density_beta      <- res_diagnostics$density_beta
    ylim_beta         <- res_diagnostics$ylim_beta
    breaks_beta       <- res_diagnostics$breaks_beta
    n_class <- length(hist_beta[[1]][[1]]$breaks)-1
  }
  
  # Initialize
  n_chains     <- param[["n_chains"]]
  if(n_chains < 5){
    cols_chains <- 1:n_chains+1
  } else cols_chains  <- colorRampPalette(brewer.pal(n_chains,"Spectral"))(n_chains)
  cols_chains2 <- makeTransparent(cols_chains,alpha=1/(max(5,n_chains)))
  lags    <- res_diagnostics$lags
  ts      <- res_diagnostics$ts
  l_ts    <- length(ts)
  
  if(is.numeric(chain)){ # For a specific chain or ...
    if(which_plot == "mu"){
      # Trace of mu
      readline(paste("Press [enter] to plot the trace of mu for the chain ",chain,": ",sep=""))
      plot(trace_mu[,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",
           ylab="",main=paste("Trace of mu for the Chain ",chain,sep=""))
      # Autocorrelation of mu
      readline(paste("Press [enter] to plot the autocorrelation of mu for the chain ",chain,":",sep=""))
      plot(lags,autocorr_lag_mu[,chain],type="h",xlab="lag",ylab="correlation",
           col=cols_chains[chain],
           main=paste("Autocorrelation of mu for the Chain ",chain,sep=""))
      # Histogram of the empirical posterior sample of mu
      readline(paste("Press [enter] to plot the histogram of the posterior sample of mu for the chain",chain,":",sep=""))
      DF_mu_tmp <- DF_mu$obs[DF_mu$chain == paste("chain",1,sep="_")]
      hist(DF_mu_tmp,nclass=n_class,col=cols_chains[chain],border=0,xlab="",
           ylab="density",main=paste("Histogram of the posterior sample of mu \n for the chain ",
                                     chain,sep="" ))
    }
    if(which_plot == "sigma_sq"){
      # Trace of sigma_sq
      readline(paste("Press [enter] to plot the trace of sigma_sq for the chain" ,chain,":",sep=""))
      plot(trace_sigma[,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",
           ylab="",main=paste("Trace of sigma for the Chain ",chain,sep=""))
      # Autocorrelation of sigma_sq
      readline(paste("Press [enter] to plot the autocorrelation of sigma_sq for the chain",chain,":",sep=""))
      plot(lags,autocorr_lag_sigma[,chain],type="h",xlab="lag",ylab="correlation",
           col=cols_chains[chain],
           main=paste("Autocorrelation of sigma for the Chain ",chain,sep=""))
      # Histogram of the empirical posterior distribution of sigma_sq
      readline(paste("Press [enter] to plot the histogram of the posterior sample of sigma_sq for the chain",chain,":",sep=""))
      DF_sigma_tmp <- DF_sigma$obs[DF_sigma$chain == paste("chain",1,sep="_")]
      hist(DF_sigma_tmp,nclass=n_class,col=cols_chains[chain],border=0,xlab="",
           ylab="density",main=paste("Histogram of the posterior sample of sigma \n for the chain ",
                                     chain,sep="" ))
    }
    if(which_plot == "beta"){
      ts_index <- which(ts == time)
      # Trace of beta( time )
      readline(paste("Press [enter] to plot the trace of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
      plot(trace_beta[ts_index,,chain],type="l",lty=1,col=cols_chains[chain],xlab="Iterations",ylab="",
           main=paste("Trace of beta(",round(data$grid[time],2),") for the Chain ",chain,sep=""))
      # Autocorrelation of beta( time )
      readline(paste("Press [enter] to plot the autocorrelation of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
      plot(lags,autocorr_lag_beta[ts_index,chain,],type="h",xlab="lag",ylab="correlation",col=cols_chains[chain],
           main=paste("Autocorrelation of beta(",round(data$grid[time],2),") for the Chain ",chain,sep=""))
      # Histogram of the empirical posterior distribution of beta( time )
      readline(paste("Press [enter] to plot the histogram of the posterior sample of beta(",round(data$grid[time],2),") for the chain",chain,":",sep=""))
      DF_beta_tmp <- DF_beta[[ts_index]]$obs[DF_beta[[ts_index]]$chain == paste("chain",1,sep="_")]
      hist(DF_beta_tmp,nclass=n_class,col=cols_chains[chain],border=0,
           main=paste("Histogram of the posterior sample of beta(",
                      round(data$grid[time],2),") \n for the chain",chain,sep="")
           ,xlab="",ylab="density")
    }
  }else{ # ... or all the chain on the same plot.
    matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
    layout(matrix_layout)
    if(which_plot == "mu"){
      # Trace of mu
      readline(paste("Press [enter] to plot the trace of mu :",sep=""))
      matplot(trace_mu,type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
              main="Trace of mu")
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      # Autocorrelation of mu
      readline(paste("Press [enter] to plot the autocorrelation of mu :",sep=""))
      plot(lags,type="n",ylim=range(autocorr_lag_mu),xlab="lag",ylab="correlation",
           main="Autocorrelation of mu ")
      for(j in 1:n_chains){
        lines(lags+j/(n_chains+1),autocorr_lag_mu[,j],type="h",col=cols_chains[j])
      }
      
      breaks_mu <- hist_mu[[1]]$breaks
      step_mu   <- diff( breaks_mu )[1]
      # Histograms of the empirical posterior distribution of mu
      readline(paste("Press [enter] to plot the histogram of the posterior sample of mu :",sep=""))
      ggplot(DF_mu, aes(x=obs, fill=chain)) +
        geom_histogram(binwidth=5*step_mu, colour="black", position="dodge") +
        scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains) 
      
      readline(paste("Press [enter] to plot another histogram of the posterior sample of mu :",sep=""))
      matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
      layout(matrix_layout)
      plot(1,type="n",xlab="",main="mu for different chains",axes = FALSE,ylab="",
           xlim=range(breaks_mu),ylim=ylim_mu)
      axis(1) ; axis(2)
      for(j in 1:n_chains){
        hist_mu[[j]]$counts <- hist_mu[[j]]$density
        lines(hist_mu[[j]],border = 0,col = cols_chains2[j])
        lines(density_mu[[j]],col=cols_chains[j])
      }
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
    }
    if(which_plot == "sigma_sq"){
      # Trace of sigma_sq
      readline(paste("Press [enter] to plot the trace of sigma_sq :",sep=""))
      matplot(trace_sigma,type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
              main="Trace of sigma")
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
      # Autocorrelation of sigma_sq
      readline(paste("Press [enter] to plot the autocorrelation of sigma_sq :",sep=""))
      plot(lags,type="n",ylim=range(autocorr_lag_sigma),xlab="lag",ylab="correlation",
           main="Autocorrelation of mu ")
      for(j in 1:n_chains){
        lines(lags+j/(n_chains+1),autocorr_lag_sigma[,j],type="h",col=cols_chains[j])
      }
      
      breaks_sigma <- hist_sigma[[1]]$breaks
      step_sigma   <- diff( breaks_sigma )[1]
      # Histograms of the empirical posterior distribution of sigma_sq
      readline(paste("Press [enter] to plot the histogram of the posterior sample of sigma_sq :",sep=""))
      ggplot(DF_sigma, aes(x=obs, fill=chain)) +
        geom_histogram(binwidth=5*step_sigma, colour="black", position="dodge") +
        scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains)
      
      readline(paste("Press [enter] to plot another histogram of the posterior sample of sigma_sq :",sep=""))
      matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
      layout(matrix_layout) 
      plot(1,type="n",xlab="",main="sigma for different chains",axes = FALSE,ylab="",
           xlim=range(breaks_sigma),ylim=ylim_sigma)
      axis(1) ; axis(2)
      for(j in 1:n_chains){
        hist_sigma[[j]]$counts <- hist_sigma[[j]]$density
        lines(hist_sigma[[j]],border = 0,col = cols_chains2[j])
        lines(density_sigma[[j]],col=cols_chains[j])
      }
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
      
    }
    if(which_plot == "beta"){
      ts_index <- which(ts == time)
      # Trace of beta( time )
      readline(paste("Press [enter] to plot the trace of beta(",round(data$grid[time],2),"):",sep=""))
      matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
      layout(matrix_layout)
      matplot(trace_beta[ts_index,,],type="l",lty=1,col=cols_chains,xlab="Iterations",ylab="",
              main=paste("Trace of beta(",round(data$grid[time],2),")",sep=""))
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
      
      # Autocorrelation of beta( time )
      readline(paste("Press [enter] to plot the autocorrelation of beta(",round(data$grid[time],2),"):",sep=""))
      matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
      layout(matrix_layout)
      plot(lags,type="n",ylim=range(autocorr_lag_beta[ts_index,,]),xlab="lag",ylab="correlation",
           main=paste("Autocorrelation of beta(",round(data$grid[time],2),") ",sep=""))
      for(j in 1:n_chains){
        lines(lags+j/(n_chains+1),autocorr_lag_beta[ts_index,j,],type="h",col=cols_chains[j])
      }
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
      
      step_beta   <- diff( breaks_beta[,ts_index] )[1]
      # Histograms of the empirical posterior distribution of beta( time )
      readline(paste("Press [enter] to plot the histogram of the posterior sample of beta(",round(data$grid[time],2),"):",sep=""))
      ggplot(DF_beta[[ts_index]], aes(x=obs, fill=chain)) +
        geom_histogram(binwidth=5*step_beta, colour="black", position="dodge") +
        scale_fill_manual(breaks=paste("chain_",1:n_chains,sep=""), values=cols_chains) 
      
      readline(paste("Press [enter] to plot another trace of beta(",round(data$grid[time],2),"):",sep=""))
      matrix_layout <- matrix(c(rep(1,3),2), 1, 4, byrow = TRUE)
      layout(matrix_layout)
      plot(1,type="n",xlab="",main=paste("beta(",round(data$beta_function[time],2),
                                         ") for different chains",sep=""),axes = FALSE,ylab="",
           xlim=range(breaks_beta[,ts_index]),ylim=ylim_beta[,ts_index])
      axis(1) ; axis(2)
      for(j in 1:n_chains){
        hist_beta[[ts_index]][[j]]$counts <- hist_beta[[ts_index]][[j]]$density
        lines(hist_beta[[ts_index]][[j]],border = 0,col = cols_chains2[j])
        lines(density_beta[[ts_index]][[j]],col=cols_chains[j])
      }
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
      points(x=rep(0.1,n_chains),y=seq(0,1,l=n_chains+2)[-c(1,n_chains+2)],pch=16,col=cols_chains,cex=2)
      text(x=rep(0.3,n_chains), y = seq(0,1,l=n_chains+2)[-c(1,n_chains+2)], 
           labels = rev(paste("chain_",1:n_chains,sep="")),cex=1)
      layout(1)
      
    }
  }
}

################################# ----
# autocorr 
################################# ----
# description : Compute the autocorrelation of the sample ( x_i(t) )_{i=1,...,n}.
# value : a symmetric matrix.
# argument 1 : data, a list containing 1) x, the functions x_i(t) and 
#                 2) grid, the grid of time points.       
# argument 2 : plot, a logical value. If it is TRUE, an image (heat map) of the
#                 autocorreraltion matrix is plotted. (optional)
# example : 
#### Autocorrelation of the function x_i(t)
# param <- list(n=50,p=100,beta_type="smooth")
# data <- sim(param)
# res_autocorr <- autocorr(data)
# cols <- rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(50))
# image(res_autocorr,col=cols)
#
#### Autocorrelation of the function beta_j(t).
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
#
# beta_functions_autocorr <- autocorr(list(grid = data$grid, x = beta_functions))
# image(beta_functions_autocorr)
autocorr <- function(data,plot=F){
  # Initialize
  x     <- data$x
  x.cor <- matrix(0,ncol(x),ncol(x))
  
  # Compute the correlations
  for(i in 1:ncol(x)){
    for(j in i:ncol(x)){
      x.cor[i,j] <- cor(x[,i],x[,j])
      x.cor[j,i] <- cor(x[,i],x[,j])
    }
  }
  
  # Plot the autocorrelation ?
  if(plot) image(x.cor)
  
  # Return the result 
  return(x.cor)
}

################################# ----
# between
################################# ----
# description : Check if a value is in an interval.
# value : a logical value.
# argument 1 : value, a numerical value.
# argument 2 : interval, a numerical vector of lenght 2 : (low,up).
# example : 
# 1 %between% c(0,2)
# 2 %between% c(0,2)
# 3 %between% c(0,2)
"%between%" <- function(value,interval){
  (value >= interval[1]) & (value <= interval[2])
}

################################# ----
# compute_beta_functions
################################# ----
# description : Compute the function beta_i(t) for each iteration i of the Gibbs
#               sample.
# value : return a matrix. Each row is a function observed on the grid of
#         time points.
# argument 1 : res.Gibbs_Sampler, a list (provided by the function 
#              Bliss_Gibbs_Sampler).
# argument 2 : param (optional),a list containing 1) p, the number of 
#              time points, 2) K, the number of intervals, 3) grid, the grid
#              of time points and 4) basis, a character vector.
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
# beta_functions <- compute_beta_functions(res.Gibbs_Sampler,param)
# indexes <- sample(1e4,1e2,replace=F)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
# matplot(param$grid,t(beta_functions[indexes,]),type="l",lty=1,col=cols)
compute_beta_functions <- function(res.Gibbs_Sampler,param){
  cat("Compute the functions beta_i. \n")
  # Initialize parameters
  K    <- param$K
  grid <- param$grid
  p    <- param$p
  
  # Initialize parameters
  basis <- param[["basis"]]
  if(is.null(basis)) basis <- "uniform"
  
  # Compute the functions beta_i for each iteration i of the Gibbs Sampler
  beta_functions <- compute_beta_functions_cpp(res.Gibbs_Sampler$trace,p,K,grid,
                                               basis,res.Gibbs_Sampler$param$scale_ml)
  # Return the functions
  return(beta_functions)
}

compute_beta_functions_mult <- function(res.Gibbs_Sampler,param){
  cat("Compute the functions beta_i. \n")
  # Initialize parameters
  K     <- param[["K"]]
  grids <- param$grids
  p     <- param[["p"]]
  
  # Initialize parameters
  basis <- param[["basis"]]
  if(is.null(basis)) basis <- rep("uniform",length(K))
  
  # Compute the functions beta_i for each iteration i of the Gibbs Sampler and for each covariable
  beta_functions <- list()
  count <- 0
  for(q in 1:length(K)){
    trace_tmp <- res.Gibbs_Sampler$trace[,(1+count):(count+3*K[q])]
    
    beta_functions[[q]] <- compute_beta_functions_cpp(trace_tmp,p[q],K[q],grids[[q]],
                                                      basis[q],res.Gibbs_Sampler$param$scale_ml[[q]])
    count <- count + 3*K[q]
  }
  # Return the functions
  return(beta_functions)
}

################################# ----
# corr_matrix 
################################# ----
# description : Compute an autocorrelation matrix to simulate functions x_i(t).
# value : a symmetric matrix.
# argument 1 : diagonal, a numerical vector corresponding to the diagonal of 
#                 the final matrix.
# argument 2 : ksi, a "coefficient of correlation". See the article Bliss, 
#                 Section 3.1 for more details.
# example : 
#### Test 1 : weak autocorrelation
# ksi     <- 1
# diagVar <- abs(rnorm(100,50,5))
# Sigma   <- corr_matrix(diagVar,ksi^2)
# persp(Sigma)
#### Test 2 : strong autocorrelation
# ksi     <- 0.2
# diagVar <- abs(rnorm(100,50,5))
# Sigma   <- corr_matrix(diagVar,ksi^2)
# persp(Sigma)
corr_matrix <- function(diagonal,ksi){  
  # Initialize 
  p <- length(diagonal)
  res <- diag(diagonal)
  
  # Compute the correlation matrix 
  for(i in 1:p-1){
    for(j in (i+1):p){
      res[i,j] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])    
      res[j,i] <- exp(-ksi*(i-j)^2/p)*sqrt(res[i,i]*res[j,j])    
    }
  }
  
  # return the matrix
  return(res)
}

################################# ----
# image_Bliss 
################################# ----
# description : Plot the representation of the estimated posterior density.
# argument 1 : posterior_density_estimate, a list. The result of the function
#                 density_estimation.
# argument 2 : param, a list containing 1) "cols", a vector of colors 
#                 for the function image (optional), 2) col_scale, a 
#                 character vector, 3) ylim, a numerical two-vector (optional) 
#                 and 4) main, a character string. 
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
#
# param$cols <- colorRampPalette(brewer.pal(9,"Reds"))(1e2)
# image_Bliss(posterior_density_estimate,param)
# lines(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",lwd=2)
# lines(param$grid,data$beta_function,col=3,lwd=2)
#
# param$cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(1e2)
# image_Bliss(posterior_density_estimate,param)
# lines(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",lwd=2)
# lines(param$grid,data$beta_function,col=3,lwd=2)
#
# param$cols <- heat.colors(12)
# param$col_scale <- "quantile"
# image_Bliss(posterior_density_estimate,param)
# lines(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",lwd=2)
# lines(param$grid,data$beta_function,col=3,lwd=2)
#
# param$cols <- terrain.colors(12)
# image_Bliss(posterior_density_estimate,param)
# lines(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",lwd=2)
# lines(param$grid,data$beta_function,col=3,lwd=2)
#
# param$cols <- topo.colors(12)
# image_Bliss(posterior_density_estimate,param)
# lines(param$grid,res.Simulated_Annealing$Bliss_estimate,type="l",lwd=2)
# lines(param$grid,data$beta_function,col=3,lwd=2)
image_Bliss <- function(posterior_density_estimate,param){
  cols      <- param$cols
  col_scale <- param$col_scale
  ylim      <- param$ylim
  main      <- param$main
  if(is.null(cols)){
    #cols <- colorRampPalette(brewer.pal(9,"Spectral"))(1e2)
    cols <- rev(heat.colors(100))
  } 
  nbre_cols <-length(cols)
  if(is.null(col_scale)) col_scale <- "regular"
  if(col_scale == "regular"){
    breaks <- seq(min(as.vector(posterior_density_estimate$res.kde2d$z)),
                  max(as.vector(posterior_density_estimate$res.kde2d$z)),
                  length=nbre_cols+1)
  } 
  if(col_scale == "quantile"){
    breaks <- quantile(as.vector(posterior_density_estimate$res.kde2d$z),
                       seq(0,1,length=nbre_cols+1))
  }
  if(is.null(ylim)) ylim <- range(posterior_density_estimate$res.kde2d$y)
  if(is.null(main)) main <- "BLiSS estimates"
  image(posterior_density_estimate$res.kde2d,col=cols,breaks = breaks,
        main=main,ylim=ylim,useRaster = TRUE)
}

################################# ----
# density_estimation
################################# ----
# description : Compute a graphical representation of the posterior 
#               distribution of beta.
# details : The sample is thinned in order to reduce the number of points and 
#           so the time of the computation of the function kde2d. 
# value : return the result of the kde2d function, i.e. an estimate of the  
#         posterior density on a two-dimensional grid.
# argument 1 : beta_functions, a list (provided by the function 
#              compute_beta_functions).
# argument 2 : param (optional), a list containing 1) burnin, an interger which
#                 is the number of iteration to drop (optional), 2) thin, a
#                 numerical value (optional), 3) lims.kde, a numerical vector 
#                 (yl,yu) (optional), 4) n, an integer related to the 
#                 precision of the heat map (optional), 5) h1, a numerical
#                 value which is the bandwidth of the kernel density estimation
#                 for the t-axis (optional) and 6) new_grid, a numerical 
#                 vector.
# example :
# load("./Example.RData")
# sourceCpp("./Bliss_rcpp.cpp")
#
# density_estimate <- density_estimation(beta_functions,param)
# image_Bliss(density_estimate,param) 
#
# param$thin <- 1e2
# density_estimate <- density_estimation(beta_functions,param)
# image_Bliss(density_estimate,param)
#
# param$lims.kde <- c(-3,8)
# density_estimate <- density_estimation(beta_functions,param)
# image_Bliss(density_estimate,param)
density_estimation <- function(beta_functions,param){
  cat("Compute the estimation of the posterior density.\n")
  # Initialize
  grid <- param$grid
  iter <- param[["iter"]]
  p    <- param$p
  
  # load optional objects
  n        <- param[["n"]]
  thin     <- param$thin
  burnin   <- param[["burnin"]]
  lims.kde <- param$lims.kde
  h1       <- param$h1
  new_grid <- param[["new_grid"]]
  
  # Initialize the necessary unspecified objects 
  max_points <- 1e5       
  if(!is.null(new_grid)) p      <- length(new_grid)
  if(is.null(n))         n      <- 512
  if(is.null(burnin))    burnin <- floor(iter/5)  
  if(is.null(thin))      thin   <- floor((iter-burnin)*p/max_points)
  
  # Check if the burnin isn't too large.
  if(2*burnin > iter){
    burnin <- floor(iter/5)
    cat("\t Burnin is too large. New burnin : ",burnin,".\n")
  }
  
  # Thin the sample of beta functions   
  thin_min   <- max(1,floor((iter-burnin)*p/max_points))
  if(thin <  thin_min){
    cat("\t 'thin = ",thin,"' is too small. Now, thin = ",
        thin_min,".\n",sep="")
    thin <- thin_min
  } 
  cat("\t Thin the sample.\n")
  beta_functions <- beta_functions[seq(1+burnin,iter,by=thin),]
  
  
  # Compute the functions beta_i on the new grid (if claimed).
  if(!is.null(new_grid)){
    beta_functions_save <- beta_functions
    beta_functions <- matrix(0,nrow(beta_functions),p)
    cat("Compute the coefficient functions on the new grid.\n")
    for(i in 1:nrow(beta_functions)){
      beta_functions[i,] <- finer_grid(beta_functions_save[i,],
                                       grid,new_grid) 
    }
    param$old_grid <- grid
    param$grid     <- new_grid
    grid           <- new_grid
  }
  
  # Filter the functions with a range too large which make a
  # noninterpretable graphical representation
  cat("\t Drop the extreme values.\n")
  t_beta <- rep(grid,nrow(beta_functions))
  y_beta <- as.vector(t(beta_functions))
  
  t_beta <- t_beta[y_beta %between% quantile(y_beta,c(0.01,0.99))]
  y_beta <- y_beta[y_beta %between% quantile(y_beta,c(0.01,0.99))]
  
  # If a window is given to plot the posterior density.
  if(!is.null(lims.kde)){
    index   <- which(y_beta %between% lims.kde)
    t_beta  <- t_beta[index]
    y_beta  <- y_beta[index]
  }
  
  # If there are too much of y_beta=0 (more than 50%), the default bandwidth of the 
  # density estimation of the kde function is 0, using the defaut function 
  # bw.nrd0. Indeed, with this function the bandwidth is :
  #     bw = 0.9 * min( \hat{\sigma} , IQR /1.34 ) * n^{-1/5}
  # So, if there are to much of y_beta=0, we choose to reduce them to 25%.
  if(sum(y_beta==0)/length(y_beta) > 0.25){    
    cat("\t Bandwidth untractrable. Some points are dropped.\n")
    res.table <- table(t_beta[y_beta==0])    
    # remove all the points which y=0
    t_beta <- t_beta[y_beta!=0]
    y_beta <- y_beta[y_beta!=0]    
    Toadd <- length(y_beta)*0.25/(1-0.25)
    
    # Add 25% of them for some correct t_beta
    t_beta_0 <- sample(as.numeric(names(res.table)),Toadd,prob=res.table,replace=T)
    y_beta_0 <- rep(0,Toadd)    
    t_beta <- c(t_beta,t_beta_0)
    y_beta <- c(y_beta,y_beta_0)
  }
  
  cat("\t Perform the 'kde2d' function.\n")
  h2 <- bandwidth.nrd(y_beta)
  if(is.null(h1)) h1 <- bandwidth.nrd(t_beta) 
  
  points <- cbind(t_beta,y_beta)
  if(is.null(lims.kde)){
    res.kde2d <- kde2d(x=t_beta,y=y_beta,n=n,h=c(h1,h2))
  }else{
    res.kde2d <- kde2d(x=t_beta,y=y_beta,lims=c(range(t_beta),lims.kde),
                       n=n,h=c(h1,h2))
  }  
  return(list(beta_functions = beta_functions,
              points         = points,
              thin           = thin,
              res.kde2d      = res.kde2d))
  cat("\t Done.\n")
}

################################# ----
# Fourier_basis_build 
################################# ----
# description : Define a Fourier basis to simulate functions x_i(t).
# value : a matrix.
# argument 1 : grid, a numerical vector.
# argument 2 : dim, a numerical value. It corresponds to dim(basis)/2.
# argument 3 : per, a numerical value which corresponds to the period of the 
#                 sine and cosine functions.
# example : see the function sim_functions (Bliss_simulate).
Fourier_basis_build <- function(grid,dim,per=2*pi){
  sapply(grid,function(x) c(cos(2*pi*x*(1:dim)/per),sin(2*pi*x*(1:dim)/per) )  ) 
}

################################# ----
# random_walk
################################# ----
# description : Compute a random walk. (gaussian)
# value : a matrix where each row is a random walk.
# argument 1 : n, an integer, the number of random walk.
# argument 2 : p, an integer, the length of the random walks.
# argument 3 : mu, a numerical vector, the average random walk.
# argument 4 : sigma, a numerical value which is the standard deviation of the 
#                 gaussian distribution used to compute the random walk.
# argument 5 : start, a numerical vector which is the initial value of 
#                 the random walks. (optional)
# example : see the function sim_functions (Bliss_simulate).
random_walk <- function(n,p,mu,sigma,start=rep(0,n)){
  res <- matrix(0,n,p)  
  for(i in 1:n){
    add     <- rnorm(p,mu,sigma)
    res[i,] <- cumsum(add)
    res[i,] <- start[i] + res[i,] 
  }  
  res <- start + res
  return(res)  
}

################################# ----
# sigmoid 
################################# ----
# description : Compute a sigmoid function.
# details : Used to simulate a coefficient function or functions x_i(t).
# value : a numerical vector.
# argument 1 : x, a numerical vector, a grid of points. 
# argument 2 : asym, the value of the asymptote of the sigmoid function. (optional)
# argument 3 : v, a numerical value which is related to the slope 
#                 at the origin. (optional)
# example : 
### Test 1 :
# x <- seq(-7,7,0.1)
# y <- sigmoid(x)
# plot(x,y,type="l",main="Sigmoid function")
### Test 2 :
# x  <- seq(-7,7,0.1)
# y  <- sigmoid(x)
# y2 <- sigmoid(x,asym=0.5)
# y3 <- sigmoid(x,v   =  5)
# plot(x,y,type="l",main="Other sigmoid functions")
# lines(x,y2,col=2)
# lines(x,y3,col=3)
sigmoid <- function(x,asym=1,v=1){
  (asym^-1 + exp(-v*x))^-1
}

################################# ----
# sigmoid_sharp 
################################# ----
# description : Compute a sharp function from the sigmoid function 
# Details : Used to simulate a coefficient function or functions x_i(t).
# value : a numerical vector.
# argument 1 : x, a numerical vector, a grid of points. 
# argument 2 : loc, a numerical value, the instant of the sharp. (optional)
# argument 3 : ... Arguments to be passed to the function sigmoid. (optional)
# example : 
### Test 1 :
# x <- seq(-7,7,0.1)
# y <- sigmoid_sharp(x)
# plot(x,y,type="l",main="Sharp sigmoid")
### Test 2 :
# x  <- seq(-7,7,0.1)
# y  <- sigmoid_sharp(x,loc=3)
# y2 <- sigmoid_sharp(x,loc=3,asym=0.5)
# y3 <- sigmoid_sharp(x,loc=3,v   =  5)
# plot(x,y,type="l",main="Other sharp sigmoids")
# lines(x,y2,col=2)
# lines(x,y3,col=3)
sigmoid_sharp <- function(x,loc=0,...){
  4*(sigmoid(x-loc,...) * sigmoid(-x+loc,...))
}
# 4 should be replace by (a+1)^2 such that the maximum of the curve 
# provided by sigmoid_sharp is 1. 