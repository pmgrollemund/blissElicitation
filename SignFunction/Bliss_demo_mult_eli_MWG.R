# setwd("~/Documents/These/Code/Methodes/BLiSS/Bliss_multiple")
#########################################################
#                                                       #
#              Bliss method : examples                  #
#                                                       # 
#########################################################
setwd("~/Documents/These/Code/Methodes/BLiSS/Bliss_multiple_elicitation/MWG/")
source("./Bliss_eli_multiple_MWG.R")
setwd("~/Documents/")
################################# ----
# Dataset test 1 : simple and smooth coefficient functions 
################################# ----
#### simulate the dataset ----
Q <- 1
param <- list(n=100,p=c(100),beta_types=c("simple"),b_inf=c(0),
              b_sup=c(1))
data <- sim_multiple(param)
#### Apply the Bliss method ----
param     <- list(iter=1e3,burnin=1e1,K=c(3),grids=data$grids,
                  prior_beta="Ridge_Zellner",phi_l=list("Gamma"),
                  conf_expert  =list(list(c(rep(0.5,30),rep(1,20),rep(0.1,30),
                                            rep(1,20)), #expert 1
                                          c(rep(0.5,20),rep(0.1,50),rep(0.2,30)), #expert 2
                                          c(rep(0.3,40),rep(0.4,20),rep(0.5,40)))), #expert 3
                  beta_s_expert=list(list(c(rep(0,30),rep(1,20),rep(-1,30),
                                            rep(0,20)), #expert 1
                                       c(rep(0,20),rep(1,50),rep(-1,30)), #expert 2
                                       c(rep(0,40),rep(1,20),rep(-1,40)))), #expert 3
                  tau=0.1)
res_Bliss_mult <- Bliss_multiple(data,param)
##### Plot the results ----
for(q in 1:Q){
  image_Bliss(res_Bliss_mult$posterior_density_estimate[[q]],param)
  # the Bliss estimate
  lines(param$grids[[q]],res_Bliss_mult$Bliss_estimate[[q]],type="s",lwd=2) 
  # the posterior expection of beta(t)
  lines(param$grids[[q]],res_Bliss_mult$res.Simulated_Annealing[[q]]$posterior_expe,type="l",lty=2)
  # the "true" coefficient function for this simulated dataset 
  lines(param$grids[[q]],data$beta_function[[q]],col=3,lwd=2)
  # and the prior beta_s
  lines(param$grids[[q]],res_Bliss_mult$res.Gibbs_Sampler$beta_s_expert_average[[q]],col=4,lty=2,type="s",lwd=2)
}

head(res_Bliss_mult$res.Gibbs_Sampler$trace,50)

boxplot(res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)-1],xlab="alpha",ylab="",main="")

sum(res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)])

res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)-1]

plot(0:param$iter,res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)-1])

hist(res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)],xlab="accepted",ylab="",main="")

mean(res_Bliss_mult$res.Gibbs_Sampler$trace[
  ,ncol(res_Bliss_mult$res.Gibbs_Sampler$trace)])
#### save the results ----
save.image("./Example_1_MWG.RData")
rm(list=ls())









################################# ----
# Dataset test 2 : une seule variable
################################# ----
Q <- 2
#### simulate the dataset ----
param2 <- list(n=1000,p=c(100,100),beta_types=c("simple","smooth"),b_inf=c(0,0),
               b_sup=c(1,1))
data2 <- sim_multiple(param2)
#### Apply the Bliss method ----
param2     <- list(iter=5e4,burnin=2e3,K=c(3,3),grids=data2$grids,
                   prior_beta="Ridge_Zellner",phi_l=list("Gamma","Gamma"),
                   beta_s_expert=list(c(rep(0,30),rep(1,20),rep(-1,30),
                                        rep(0,20)),
                                      c(rep(0,30),rep(1,20),rep(-1,30),
                                        rep(0,20))
                                      )
                   )
res_Bliss_mult2 <- Bliss_multiple(data2,param2)
##### Plot the results ----
for(q in 1:Q){
  image_Bliss(res_Bliss_mult2$posterior_density_estimate[[q]],param2)
  # the Bliss estimate
  lines(param2$grids[[q]],res_Bliss_mult2$Bliss_estimate[[q]],type="s",lwd=2) 
  # the posterior expection of beta(t)
  lines(param2$grids[[q]],res_Bliss_mult2$res.Simulated_Annealing[[q]]$posterior_expe,type="l",lty=2)
  # the "true" coefficient function for this simulated dataset 
  lines(param2$grids[[q]],data2$beta_function[[q]],col=3,lwd=2)
}
#### save the results ----
save.image("./Example_2.RData")
rm(list=ls())

################################# ----
# Dataset test 3 : MH saut normaux 
################################# ----
#### simulate the dataset ----
Q <- 1
param <- list(n=1000,p=c(100),beta_types=c("simple"),b_inf=c(0),
              b_sup=c(1))
data <- sim_multiple(param)
#### Apply the Bliss method ----
param     <- list(iter=5e4,burnin=2e3,K=c(3),grids=data$grids,
                  prior_beta="Ridge_Zellner",phi_l=list("Gamma"),
                  beta_s_expert=list(c(rep(0,30),rep(1,20),rep(-1,30),
                                       rep(0,20))),
                  proposal_distribution = "random_walk")
res_Bliss_mult <- Bliss_multiple(data,param)
##### Plot the results 
for(q in 1:Q){
  image_Bliss(res_Bliss_mult$posterior_density_estimate[[q]],param)
  # the Bliss estimate
  lines(param$grids[[q]],res_Bliss_mult$Bliss_estimate[[q]],type="s",lwd=2) 
  # the posterior expection of beta(t)
  lines(param$grids[[q]],res_Bliss_mult$res.Simulated_Annealing[[q]]$posterior_expe,type="l",lty=2)
  # the "true" coefficient function for this simulated dataset 
  lines(param$grids[[q]],data$beta_function[[q]],col=3,lwd=2)
}
#### save the results ----
save.image("./Example_3.RData")
rm(list=ls())

################################# ----
# Dataset test 3 : MH saut normaux 
################################# ----
#### simulate the dataset ----
Q <- 1
param <- list(n=100,p=c(100),beta_types=c("simple"),b_inf=c(0),
              b_sup=c(1))
data <- sim_multiple(param)
#### Apply the Bliss method ----
param     <- list(iter=1e5,burnin=2e3,K=c(3),grids=data$grids,
                  prior_beta="Ridge_Zellner",phi_l=list("Gamma"),
                  beta_s_expert=list(c(rep(0,30),rep(1,20),rep(-1,30),
                                       rep(0,20))),
                  proposal_distribution = "random_walk")
res_Bliss_mult <- Bliss_multiple(data,param)
##### Plot the results 
for(q in 1:Q){
  image_Bliss(res_Bliss_mult$posterior_density_estimate[[q]],param)
  # the Bliss estimate
  lines(param$grids[[q]],res_Bliss_mult$Bliss_estimate[[q]],type="s",lwd=2) 
  # the posterior expection of beta(t)
  lines(param$grids[[q]],res_Bliss_mult$res.Simulated_Annealing[[q]]$posterior_expe,type="l",lty=2)
  # the "true" coefficient function for this simulated dataset 
  lines(param$grids[[q]],data$beta_function[[q]],col=3,lwd=2)
}
#### save the results ----
save.image("./Example_4.RData")
rm(list=ls())