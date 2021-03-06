#########################################################
#                                                       #
#          Bliss method : simulate dataset              #
#                                                       # 
#########################################################

################################# ----
# sim
################################# ----
# description : Simulate a dataset for the functional linear model. 
# value : a list containing :
#             y : a numerical vector, the outcomes y_i.
#             x : a numerical vector, the functional covariates on a grid of
#                 times points (grid). 
#             beta_function : a numerical vector, the coefficient function on 
#                 a grid of times points (grid).
#             grid : a numerical vector, which is the grid of observation points.
# argument 1 : param, a list containing :
#                 n : an integer, the number of functions.
#                 p : an integer, the number of observation times.
#                 b_inf : a numerical value, the lower boundary of the support. (optional)
#                 b_sup : a numerical value, the upper boundary of the support. (optional)
#                 mu : a numerical value, the intercept of the model. (optional)
#                 r : a nonnegative value, the signal to noise ratio. (optional)
#                 link : a function, the link function to simulate data from a GFLM. (optional)
#                 grid : a numerical vector, the grid of observation points.
#                 beta_type : a character vector. It indicates the shape of 
#                     the coefficient function. (optional)
#                 x_type : a character vector. It indicates the shape of the 
#                     functions x_i(t). (optional)
#                 dim : a numerical value, the dimension of the Fourier basis, 
#                     if "type" is "Fourier" or "Fourier2". (optional)
#                 ksi : a numerical value which is a "coefficient of correlation", 
#                     see the Bliss article Section 3.1 for more details.
#                 diagVar : a numerical vector, the diagonal of the 
#                     autocorrelation matrix of the functions x_i(t).
# example :
### Dataset 1
# param <- list(n=50,p=100,beta_type="smooth")
# data <- sim(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(50)
# par(mfrow=c(2,1))
# matplot(data$grid,t(data$x),type="l",lty=1,col=cols)
# plot(data$grid,data$beta_function,type="l")
# abline(h=0,lty=2,col="gray")
# par(mfrow=c(1,1))
### Dataset 1
# param <- list(n=10,p=200,beta_type="simple",b_inf=-2.12,b_sup=3.14)
# data <- sim(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)
# par(mfrow=c(2,1))
# matplot(data$grid,t(data$x),type="l",lty=1,col=cols)
# plot(data$grid,data$beta_function,type="l")
# abline(h=0,lty=2,col="gray")
# par(mfrow=c(1,1))
sim <- function(param){      
  cat("Simulation of the data.\n")
  # load objects 
  n     <- param$n
  p     <- param$p 
  
  # load optional objects 
  b_inf <- param$b_inf
  b_sup <- param$b_sup
  mu    <- param$mu
  r     <- param$r
  link  <- param$link
  
  # Initialize the necessary unspecified objects 
  if(is.null(b_inf)) b_inf <- 0
  if(is.null(b_sup)) b_sup <- 1
  if(is.null(mu))    mu    <- 1
  if(is.null(r))     r     <- 5
  if(is.null(link))  link  <- function(x,b) return(x*b)
  
  # Deduce objects
  param$grid <- seq(b_inf,b_sup,length=p)
  
  # Simulated the functions x_i(t) on the grid.
  cat("\t Simulate functions x_i(t).\n")
  x <- sim_functions(param)
  
  # Choose a coefficient function beta
  cat("\t Choose a coefficient function.\n")
  beta_function <- choose_beta(param)
  
  # Compute the expectation of the outcomes y  
  cat("\t Compute the outcomes y_i.\n")
  x_beta <- matrix(0,nrow(x),ncol(x))
  y_expe <- vector("numeric",n)
  for(i in 1:n){
    x_beta[i,] <- link( x[i,] , beta_function ) 
    y_expe[i]  <- mu + integrate_trapeze(param$grid,x_beta[i,])
  }   
  
  # Compute the error \varepsilon
  err <- rnorm(n,0,1)
  err <- sd(y_expe) * err / (sd(err) * sqrt(r))
  
  # Compute the outcomes y_i
  y <- y_expe + err
  
  # Return the data.
  return(list(y             = y,
              expe          = y_expe,
              x             = x,
              beta_function = beta_function,
              grid          = param$grid))
}

################################# ----
# sim_multiple
################################# ----
# description : Simulate a dataset for the functional linear model. 
# value : a list containing :
#             y : a numerical vector, the outcomes y_i.
#             x_mult : a list of matrices, the qth matrice representing the qth functional covariate on a grid of
#                 times points (grid[[q]]). 
#             beta_function_mult : a list of numerical vectors, the qth vector representing the coefficient function of 
#                 the qth covariate on a grid of times points (grid[[q]]).
#             grids : a list of numerical vectors, the qth vector is the grid of observation points for the qth covariate.
# argument 1 : param, a list containing :
#                 n : an integer, the number of observations.
#                 p : a vector of integers, the qth component is the number of observation times for the qth covariate.
#                 b_inf : a vector of numerical values, the qth component is the lower boundary of the support for the 
#                      observation times for the qth covariate. (optional)
#                 b_sup : a vector of numerical values, the qth component is the upper boundary of the support for the 
#                      observation times for the qth covariate. (optional)
#                 mu : a numerical value, the intercept of the model. (optional)
#                 r : a nonnegative value, the signal to noise ratio. (optional)
#                 links : a list of functions, the qth element of the list is a link function to simulate data from a GFLM. (optional)
#                 grids : a list of numerical vectors, the qth vector is the grid of observation points for the qth covariate
#                      (optional, can be deduced from b_inf,b_sup and p).
#                 beta_types : a character vector. It indicates the shapes of 
#                     the coefficient functions. If we have Q covariates, Q shapes can be indicated.
#                 x_types : a character vector. It indicates the shapes of the 
#                     functions x_qi(t). If we have Q covariates, Q shapes can be indicated. (optional)
#                 dim : vector of numerical values, the dimensions of the Fourier basis, if "type" is "Fourier" or "Fourier2". 
#                      One value for each covariate. (optional)
#                      If the qth shape is not "Fourier" ou "Fourier2", the associated component should be specified as NA (optional)
#                 ksi : vector of numerical values which are "coefficients of correlation" (one value for each covariate), 
#                     see the Bliss article Section 3.1 for more details. (optional)
#                 diagVar : a list of numerical vectors, the qth vector is the diagonal of the 
#                     autocorrelation matrix of the functions x_qi(t). (optional)
# example :
# param <- list(n=10,p=c(200,100),beta_types=c("simple","smooth"),b_inf=c(0,-2.12),b_sup=c(1,3.14))
# data <- sim_multiple(param)
# data$y
# data$expe
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)
# q=2
# par(mfrow=c(2,1))
# matplot(data$grids[[q]],t(data$x_mult[[q]]),type="l",lty=1,col=cols)
# plot(data$grids[[q]],data$beta_function_mult[[q]],type="l")
# abline(h=0,lty=2,col="gray")
# par(mfrow=c(1,1))

sim_multiple <- function(param){      
  cat("Simulation of the data.\n")
  # load objects 
  n     <- param[['n']]
  p     <- param[['p']] 
  
  # load optional objects 
  b_inf <- param[['b_inf']]
  b_sup <- param[['b_sup']]
  mu    <- param[['mu']]
  r     <- param[['r']]
  links  <- param[['links']]
  grids  <- param[['grids']]
  beta_types <- param[['beta_types']]
  x_types <- param[['x_types']]
  dim <- param[['dim']]
  ksi <- param[['ksi']]
  diagVar <- param[['diagVar']]
  
  # Initialize the necessary unspecified objects 
  Q <- length(p)
  if(is.null(b_inf)) b_inf <- rep(0,Q)
  if(is.null(b_sup)) b_sup <- rep(1,Q)
  if(is.null(mu))    mu    <- 1
  if(is.null(r))     r     <- 5
  if(is.null(links)) {
      links  <- list()
      for (q in 1:Q){links[[q]] <- function(x,b) return(x*b)}
  }
  
  # Deduce objects
  if(is.null(grids)) { 
    grids  <- list()
    for (q in 1:Q){grids[[q]] <- seq(b_inf[q],b_sup[q],length=p[q])}
    param[['grids']] <- grids 
  }
  if(!is.null(grids)) {
    coherence <- TRUE
    for (q in 1:Q){ coherence <- coherence * length(grids[[q]])==p[q] }
    if(coherence == FALSE) stop("The lengths of the grids (parameter grids) should correspond to the number of observation times (parameter p).")
  }
  
  # Simulated the functions x_qi(t) on the grids.
  cat("\t Simulate functions x_qi(t).\n")
  x_mult <- list()
  for (q in 1:Q){
    # dim[q] est soit NULL, soit NA, soit c'est un entier. Il faut gérer ces cas sachant que sim_functions marche pour un parametre dim NULL ou entier mais pas NA.
    if(is.null(dim[q])) {param_sim_functions <- list(n=n,p=p[q],grid=grids[[q]],x_type=x_types[q],dim=NULL,ksi=ksi[q],diagVar=diagVar[[q]])}
    if(!is.null(dim[q])) {
      if(is.na(dim[q])) {
        param_sim_functions <- list(n=n,p=p[q],grid=grids[[q]],x_type=x_types[q],dim=NULL,ksi=ksi[q],diagVar=diagVar[[q]])
      }else{
        param_sim_functions <- list(n=n,p=p[q],grid=grids[[q]],x_type=x_types[q],dim=dim[q],ksi=ksi[q],diagVar=diagVar[[q]])
      }
    }
    x_mult[[q]] <- sim_functions(param_sim_functions)
  }

  # Choose a coefficient function beta
  cat("\t Choose a coefficient function.\n")
  beta_function_mult <- list()
  for (q in 1:Q){
    param_choose_beta <- list(p=p[q],grid=grids[[q]],beta_type=beta_types[q])
    beta_function_mult[[q]] <- choose_beta(param_choose_beta)
  }
  
  # Compute the expectation of the outcomes y  
  cat("\t Compute the outcomes y_i.\n")
  y_expe <- rep(mu,n)
  for(i in 1:n){
    for(q in 1:Q){
      x_beta <- links[[q]]( x_mult[[q]][i,] , beta_function_mult[[q]] ) 
      y_expe[i]  <- y_expe[i] + integrate_trapeze(grids[[q]],x_beta)     
    }
  }   
  
  # Compute the error \varepsilon
  err <- rnorm(n,0,1)
  err <- sd(y_expe) * err / (sd(err) * sqrt(r))
  
  # Compute the outcomes y_i
  y <- y_expe + err
  
  # Return the data.
  return(list(y             = y,
              #expe          = y_expe,
              x_mult             = x_mult,
              beta_function_mult = beta_function_mult,
              grids          = grids))
}

################################# ----
# sim_functions
################################# ----
# description : Simulate the functions x_qi(t) for a given q.
# details : Several shape are available for the functions x_qi(t).
# value : a matrix. Each row is a function x_qi(t) (i=1,...,n and q is fixed).
# argument 1 : param, a list containing :
#                 n : an integer, the number of functions.
#                 p : an integer, the number of observation times.
#                 grid : a numerical vector, the grid of observation times.
#                 x_type : a character vector, the shape of the 
#                     functions x_i(t). (optional)
#                 dim : a numerical value, the dimension of the Fourier basis, 
#                     if "type" is "Fourier" or "Fourier2". (optional)
#                 ksi : a numerical value, a "coefficient of correlation", 
#                     see the Bliss article Section 3.1 for more details.
#                 diagVar : a numerical vector, the diagonal of the 
#                     autocorrelation matrix of the functions x_i(t).
# example :
#### Fourier 
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="Fourier")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### Fourier2
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="Fourier2")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### random_walk
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_walk")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### random_sharp
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="random_sharp")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### uniform
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="uniform")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### gaussian
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="gaussian")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
#### mvgauss 
# param <- list(n=15,p=100,grid=seq(0,1,length=100),x_type="mvgauss")
# x <- sim_functions(param)
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(15)
# matplot(param$grid,t(x),type="l",lty=1,col=cols)
sim_functions <- function(param){  
  # load objects 
  n <- param$n
  p <- param$p
  grid <- param$grid
  
  # load optional objects 
  type    <- param$x_type
  dim     <- param$dim
  ksi     <- param$ksi
  diagVar <- param$diagVar
  
  # Initialize the necessary unspecified objects 
  if(is.null(type))    type    <- "mvgauss"
  if(is.null(dim))     dim     <- 4
  if(is.null(ksi))     ksi     <- 1
  if(is.null(diagVar)) diagVar <- abs(rnorm(p,1,1/10))
  
  # Deduce objects
  by <- diff(grid)[1]
  
  # Simulate the functions x_i(t)
  if(type == "Fourier"){
    # Set a Fourier basis 
    Fourier_basis <- Fourier_basis_build(grid = grid,
                                         dim  = dim,
                                         per  = 1.5*(max(grid)-min(grid)))    
    # Choose the coefficients 
    a_n <- mvrnorm(n,(dim:1)/dim,
                   diag((dim:1)/(50*dim)))
    b_n <- mvrnorm(n,(dim:1)/dim,
                   diag((dim:1)/(50*dim)))
    
    # Compute the functions x_i(t)
    x <- a_n %*% Fourier_basis[1:dim,] + b_n %*% Fourier_basis[(dim+1):(2*dim),]    
  }
  if(type == "Fourier2"){
    # Set a Fourier basis 
    Fourier_basis <- Fourier_basis_build(grid = grid,
                                         dim  = dim,
                                         per  = 1.5*(max(grid)-min(grid)))   
    # Choose the coefficients 
    a_n <- runif(n*dim,-3,3)
    dim(a_n) <- c(n,dim)
    b_n <- runif(n*dim,-3,3)
    dim(b_n) <- c(n,dim)
    
    # Determiner les courbes 
    x <- a_n %*% Fourier_basis[1:dim,] + b_n %*% Fourier_basis[(dim+1):(2*dim),]
  }
  if(type == "random_walk"){
    start <- rnorm(n,0,2)
    x <- random_walk(n,p,0,1,start)
  }
  if(type == "random_sharp"){    
    locs <- runif(n*2,grid[1],tail(grid,1))
    dim(locs) <- c(n,2)
    
    asyms <- runif(n*2,1,5)
    dim(asyms) <- c(n,2)
    
    vs <- runif(n*2, 1/(4*by), 1/(3*by) )
    dim(vs) <- c(n,2)
    
    s <- sample(c(-1,1),2*n,replace=T)
    dim(s) <- c(n,2)
    
    x <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- s[i,1] * sigmoid_sharp(grid,locs[i,1],asyms[i,1],vs[i,1]) +
        s[i,2] * sigmoid_sharp(grid,locs[i,2],asyms[i,2],vs[i,2])
    }
  }
  if(type == "uniform"){
    x <- matrix(0,n,p)
    for(j in 1:p){
      x[,j] <- runif(n,-5,5)
    }
  }
  if(type == "gaussian"){
    x <- matrix(0,n,p)
    for(j in 1:p){
      x[,j] <- rnorm(n,0,4)
    }
  }
  if(type == "mvgauss"){   
    mu      <- (1:p-p/2)^2/(p^2/4)
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- mvrnorm(1,mu,Sigma)  
    }
  }
  if(type == "mvgauss_different_scale"){   
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 10 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- mvrnorm(1,mu,Sigma)  
    }
  }
  if(type == "mvgauss_different_scale2"){   
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 100 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- mvrnorm(1,mu,Sigma)  
    }
  }
  if(type == "mvgauss_different_scale3"){   
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[1:floor(p/3)] <- 1000 * diagVar[1:floor(p/3)]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- mvrnorm(1,mu,Sigma)  
    }
  }
  if(type == "mvgauss_different_scale4"){   
    mu      <- (1:p-p/2)^2/(p^2/4)
    diagVar[floor(2*p/3):p] <- 100 * diagVar[floor(2*p/3):p]
    Sigma   <- corr_matrix(diagVar,ksi^2)
    x       <- matrix(0,n,p)
    for(i in 1:n){
      x[i,] <- mvrnorm(1,mu,Sigma)  
    }
    
  }
  
  # Return the functions
  return(x)
}

################################# ----
# choose_beta
################################# ----
# description : Compute a coefficient function for the function linear model.
# details : Several shape are available for the coefficient function.
# value : a numerical vector which is the observation of the coefficient 
#             function at the given grid of times points: "grid".  
# argument 1 : param, a list containing 1) p, a numerical value, 2) grid, 
#                 a numerical vector and 3) beta_type, a character vector.
# example :
#### smooth
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="smooth")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="l")
#### random_smooth
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="random_smooth")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="l")
#### simple
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="simple")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="s")
#### simple2
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="simple2")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="s")
#### random_simple
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="random_simple")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="s")
#### sinusoid
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="sinusoid")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="l")
#### flat_sinusoid
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="flat_sinusoid")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="l")
#### sharp
# param <- list(p=100,grid=seq(0,1,length=100),beta_type="sharp")
# beta_function <- choose_beta(param)
# plot(param$grid,beta_function,type="l")
choose_beta <- function(param){
  # load objects 
  p    <- param$p
  grid <- param$grid
  type <- param$beta_type
  
  # Compute a "standard" grid on (0,1).
  grid2 <- (grid - min(grid))/ (max(grid) - min(grid))
  # Choose a function beta
  if(type == "smooth"){
    beta_function <- 5*exp(-((grid2-0.25)*20)^2) +
      2*exp(-((grid2-0.75)*20)^2) - 
      2*exp(-((grid2-0.5)*20)^2)
  } 
  if(type == "random_smooth"){
    beta_function <- runif(1,-5,5)*exp(-((grid2-runif(1,0,1))*20)^2) +
      runif(1,-5,5)*exp(-((grid2-runif(1,0,1))*20)^2) + 
      runif(1,-5,5)*exp(-((grid2-runif(1,0,1))*20)^2)
  } 
  if(type == "simple"){
    beta_function <- rep(0,p)
    beta_function[round(p/10):round(3*p/10)] <- 3
    beta_function[round(4.5*p/10):round(5.5*p/10)] <- 4
    beta_function[round(8*p/10):round(9.5*p/10)] <- -1
  } 
  if(type == "simple_K10"){
    beta_function <- rep(0,p)
    beta_function[round(0.5*p/10):round(2*p/10)]   <- 1 + beta_function[round(0.5*p/10):round(2*p/10)]
    beta_function[round(p/10):round(2*p/10)]       <- 2 + beta_function[round(p/10):round(2*p/10)]
    beta_function[round(0.8*p/10):round(1.7*p/10)] <- 1 + beta_function[round(0.8*p/10):round(1.7*p/10)]
    beta_function[round(4.5*p/10):round(7*p/10)]   <- 2 + beta_function[round(4.5*p/10):round(7*p/10)]
    beta_function[round(5*p/10):round(7*p/10)]     <- 1 + beta_function[round(5*p/10):round(7*p/10)]
    beta_function[round(5*p/10):round(6*p/10)]     <- 2 + beta_function[round(5*p/10):round(6*p/10)]
    beta_function[round(8*p/10):round(9.5*p/10)]   <- -0.5 + beta_function[round(8*p/10):round(9.5*p/10)]
    beta_function[round(8*p/10):round(10*p/10)]    <- -1 + beta_function[round(8*p/10):round(10*p/10)]
    beta_function[round(8*p/10):round(9.5*p/10)]   <- -1 + beta_function[round(8*p/10):round(9.5*p/10)]
    beta_function[round(8.8*p/10):round(9.5*p/10)] <- -0.5 + beta_function[round(8.8*p/10):round(9.5*p/10)]
  } 
  if(type == "simple2"){
    beta_function <- rep(0,p)
    beta_function[round(p/10):round(3*p/10)] <- 3
    beta_function[round(8*p/10):round(9.5*p/10)] <- -1
  } 
  if(type == "simple3"){
    beta_function <- rep(0,p)
    beta_function[round(p/10):round(7*p/10)] <- 3
    beta_function[round(8*p/10):round(9.5*p/10)] <- -1
  } 
  if(type == "random_simple"){
    beta_function <- rep(0,p)
    boundaries <- sort(sample(1:p,6))
    beta_function[boundaries[1]:boundaries[2]] <- runif(1,-5,5)
    beta_function[boundaries[3]:boundaries[4]] <- runif(1,-5,5)
    beta_function[boundaries[5]:boundaries[6]] <- runif(1,-5,5)
  } 
  if(type == "sinusoid"){
    beta_function <- sin(grid2 * 2* pi)
  }
  if(type == "flat_sinusoid"){
    beta_function <- rep(0,p)
    flat          <- round(p/3)
    beta_function[1:flat]  <- sin(10/(p-flat+10)               * 2* pi)
    beta_function[flat:p]  <- sin((10:(p-flat+10))/(p-flat+10) * 2* pi)
    beta_function          <- beta_function * sigmoid(1:p-flat)
  }
  if(type == "sharp"){
    beta_function <- rep(0,p)  
    shift         <- max(grid) - min(grid)
    
    beta_function <- beta_function +
      2 * sigmoid_sharp(grid,min(grid) + 0.2 * shift,v=100,asym=1) -
      3 * sigmoid_sharp(grid,min(grid) + 0.6 * shift,v=100,asym=1)  
  }
  
  # Return the chosen function
  return(beta_function)
}