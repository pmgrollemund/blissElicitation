#########################################################
#                                                       #
#       Bliss method : elicitation conf to weights      #
#                                                       # 
#########################################################
#### Test (2 experts) ----
library(mvtnorm)
n <- 5
rho <- 0.8
conf <- c(rep(0.8,n),rep(0.7,n))
# Covariance matrix
V <- diag(x = 1,2*n)
# diag(V[1:5,6:10]) <- rho
# diag(V[6:10,1:5]) <- rho
for(i in 1:n){
  V[i,n+i] = rho * sqrt(V[i,i]*V[n+i,n+i])
  V[n+i,i] = rho * sqrt(V[i,i]*V[n+i,n+i])
}
sum(diag(V)) / t(rep(1,2*n)) %*% V %*% rep(1,2*n) * (2* n)
sum(diag(V)) / t(rep(1,2*n)) %*% V %*% rep(1,2*n) * (2* sum(conf)/2)

View(V)

0.8 * 4.166667 / n
0.7 * 4.166667 / n

# Correlation matrix
C <- diag(x = rep(1,2*n),2*n)
for(i in 1:n){
  C[i,n+i] = rho 
  C[n+i,i] = rho 
}
sum(diag(C)) / t(rep(1,2*n)) %*% C %*% rep(1,2*n) * (2*n)
sum(diag(C)) / t(rep(1,2*n)) %*% C %*% rep(1,2*n) * (2* sum(conf)/2)
#### Test experts and data ----
n <- 10
E <- 4
rho <- diag(1,E)
for(j in 1:E){
  for(k in 1:E){
    # Expert 1
    if(j == 1 & k == 2){
      rho[j,k] = 0.8
      rho[k,j] = rho[j,k] 
    } 
    if(j == 1 & k == 3){
      rho[j,k] = 0.8
      rho[k,j] = rho[j,k] 
    } 
    if(j == 1 & k == 4){
      rho[j,k] = 0.8
      rho[k,j] = rho[j,k] 
    } 
    # Expert 2
    if(j == 2 & k == 3){
      rho[j,k] = 0
      rho[k,j] = rho[j,k] 
    } 
    if(j == 2 & k == 4){
      rho[j,k] = 0
      rho[k,j] = rho[j,k] 
    } 
    # Expert 3
    if(j == 3 & k == 4){
      rho[j,k] = 0
      rho[k,j] = rho[j,k] 
    } 
  }
}
  
C <- diag(x=1,n*E)
for(j in 1:E){
  for(k in 1:E){
    for(i in 1:n){
      C[ (j-1)*n + i , (k-1)*n + i] = rho[j,k] 
      C[ (k-1)*n + i , (j-1)*n + i] = rho[j,k]  
    }
  }
}

# View(C)
sum(diag(C)) / t(rep(1,E*n)) %*% C %*% rep(1,E*n) * (E*n)
#### Test pour Meili (nouvelle formule) ----
n <- 5
E <- 3
rho <- diag(1,E)
rho[1,2] <- rho[2,1] <- 1
rho[1,3] <- rho[3,1] <- 0
rho[2,3] <- rho[3,2] <- 0

compute_ESS <- function(n,E,conf,rho){
  ESS <- 0
  ESS_k <- rep(0,E)
  weights <- matrix(0,ncol = n, nrow=E) 
  
  rho_k_sq <- rep(0,E)
  for(k in 1:E){
    rho_k_sq[k] <- sum( rho[k,-k]^2)
    ESS_k[k]    <- sum(conf[k,]) / (1 + rho_k_sq[k]) 
    weights[k,] <- conf[k,] * ESS_k[k] / sum(conf[k,])
  }
  ESS <- sum(ESS_k)
  
  comparison <- data.frame( conf , apply(conf,1,sum), rep("->",E) , weights, ESS_k) 
  colnames(comparison) <- c( paste("c",1:n,sep=""),"sum_conf", "", paste("w",1:n,sep=""), "ESS_k")
  comparison
  
  return(list(ESS=ESS,
              ESS_k=ESS_k,
              weights = weights,
              comparison=comparison))
}
# test 1 :
conf <- rbind( c(1,1,1,1,1), 
               c(1,1,1,1,1),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison

# test 2 :
conf <- rbind( c(0.5,0.5,0.5,0.5,0.5), 
               c(1,1,1,1,1),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison

# test 3 :
conf <- rbind( c(1,1,1,1,1), 
               c(0.1,0.2,0.3,0.4,0.5),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison

# test 4 :
rho <- diag(1,E)
rho[1,2] <- rho[2,1] <- 0.5
rho[1,3] <- rho[3,1] <- 0
rho[2,3] <- rho[3,2] <- 0


conf <- rbind( c(1,1,1,1,1), 
               c(1,1,1,1,1),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison

# test 5 :
rho <- diag(1,E)
rho[1,2] <- rho[2,1] <- 0.75
rho[1,3] <- rho[3,1] <- 0
rho[2,3] <- rho[3,2] <- 0


conf <- rbind( c(1,1,1,1,1), 
               c(1,1,1,1,1),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison

# test 6 :
rho <- diag(1,E)
rho[1,2] <- rho[2,1] <- 0.75
rho[1,3] <- rho[3,1] <- 0.75
rho[2,3] <- rho[3,2] <- 0.75


conf <- rbind( c(1,1,1,1,1), 
               c(1,1,1,1,1),
               c(1,1,1,1,1))
res_ESS <- compute_ESS(n,E,conf,rho)

res_ESS$ESS
res_ESS$comparison
