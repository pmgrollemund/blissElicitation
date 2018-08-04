#########################################################
#                                                       #
#              Bliss eli method : example               #
#                                                       # 
######################################################### ----
source("./Bliss_eli_multiple_WL.R")
#### Simulate observed data ----
Q <- 2
param <- list(n=50,p=c(10,20),beta_types=c("simple","smooth"),b_inf=c(0,0),
              b_sup=c(1,1))
data <- sim_multiple(param)
#### Apply Bliss_eli without experts' information ----
param     <- list(iter=1e4,burnin=2e3,K=c(3,3),grids=data$grids,
                  phi_l=list("Gamma","Gamma"))
res_Bliss_mult_WL <- Bliss_multiple(data,param)
#### Simulate pseudo data ----
E <- 2
n_e <- c(20,20)
# simulate the expert responses as responses of another model
param <- list(n=n_e[1],p=c(10,20),beta_types=c("smooth","smooth"),b_inf=c(0,0),
              b_sup=c(1,1))
data_expert_1 <- sim_multiple(param)
param <- list(n=n_e[2],p=c(10,20),beta_types=c("sinusoid","simple"),b_inf=c(0,0),
              b_sup=c(1,1))
data_expert_2 <- sim_multiple(param)

y_expert <- list(matrix(data_expert_1$y,1,n_e[1],
                        byrow = T),
                 matrix(data_expert_2$y,1,n_e[2],
                        byrow = T))
x_expert <- list(list(data_expert_1$x_mult[[1]],
                      data_expert_2$x_mult[[1]]),
                 list(data_expert_1$x_mult[[2]],
                      data_expert_2$x_mult[[2]]))
conf <- list(matrix(rep(1,n_e[1]),1,n_e[1],
                    byrow = T),
             matrix(rep(1,n_e[2]),1,n_e[2],
                    byrow = T))
rho_mat <- diag(E)
#### Apply Bliss_eli with experts' information ----
param     <- list(iter=1e4,burnin=2e3,K=c(3,3),grids=data$grids,
                  phi_l=list("Gamma","Gamma"),y_expert=y_expert, 
                  x_expert = x_expert,conf=conf,
                  rho_mat = rho_mat)
param$posterior <- TRUE
res_Bliss_mult_WL_posterior <- Bliss_multiple(data,param)
param$posterior <- FALSE
res_Bliss_mult_WL_prior <- Bliss_multiple(data,param)
#### Plot the results ----
param$main <- ""
par(mfrow=c(2,1))
for(q in 1:Q){
  param$ylim <- range( c(res_Bliss_mult_WL$posterior_density_estimate[[q]]$res.kde2d$y,
                         data$beta_function_mult[[q]])) + c (-0.5,0.5)
  image_Bliss(res_Bliss_mult_WL$posterior_density_estimate[[1]],param)
  lines(param$grids[[q]],
        res_Bliss_mult_WL$res.Simulated_Annealing[[q]]$posterior_expe,lty=2)
  
  image_Bliss(res_Bliss_mult_WL_posterior$posterior_density_estimate[[q]],param)
  lines(param$grids[[q]],
        res_Bliss_mult_WL_prior$res.Simulated_Annealing[[q]]$posterior_expe,
        lty=2,col=4)
  lines(param$grids[[q]],
        res_Bliss_mult_WL_posterior$res.Simulated_Annealing[[q]]$posterior_expe)
  lines(param$grids[[q]],
        res_Bliss_mult_WL$res.Simulated_Annealing[[q]]$posterior_expe,lty=2,col=5)
}