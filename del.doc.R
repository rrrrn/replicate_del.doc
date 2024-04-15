source("site_update_param.R")
library(foreach)
library(doParallel)

del.doc <- function(max_iter=2, tol=1e-4, Omega_list_1, Omega_list_2, beta_init, t_init, N, X_list, Y_list){
  ## beta_init: K x d matrix
  ## t_init: K x d matrix
  num_cores <- detectCores()-1
  registerDoParallel(num_cores)
  K = dim(beta_init)[1]
  d = dim(beta_init)[2]
  # init a and b
  a <- b <- apply(beta_init, 2, mean)
  # init u and u2
  u <- u2 <- matrix(rep(0, K*d), nrow=K)
  
  beta_new <- matrix(0, K, d)
  t_new <- matrix(0, K, d)
  
  diff <- Inf
  iter <- 0
  while (iter < max_iter & diff > tol){
    iter <- iter + 1
    for(i in 1:K){
      res <- site_update_param(beta_init = beta_init[i,], t_init = t_init[i,], Omega1 = Omega_list_1[[i]], Omega2 = Omega_list_2[[i]], u = u[i,], u2 = u2[i,], 
                               X_list[[i]], Y_list[[i]], a = a, b= b, N = N)
      beta_new[i,] <- res$beta
      t_new[i,] <- res$t
    }
    
    a_sum <- rep(0, d); b_sum <- rep(0, d); a_denom <- matrix(0, d, d); b_denom <- matrix(0, d, d)
    for(i in 1:K){
      a_sum <- a_sum + Omega_list_1[[i]] %*% matrix(beta_new[i,], ncol=1)
      b_sum <- b_sum + Omega_list_2[[i]] %*% matrix(t_new[i,], ncol=1)
      a_denom <- a_denom + Omega_list_1[[i]]
      b_denom <- b_denom + Omega_list_2[[i]]
    }
    
    a <- (solve(a_denom) %*% a_sum)[,1]
    b <- (solve(b_denom) %*% b_sum)[,1]
    
    print(a)
    
    for(i in 1:K){
      u[i,] <- beta_new[i,] - a + u[i,]
      u2[i,] <- t_new[i,] - b + u2[i,]
    }
    
    # check convergence
    if(max(abs(beta_new - beta_init)) < tol & max(abs(t_new - t_init)) < tol){
      break
    }
    
    beta_init <- beta_new
    t_init <- t_new
    print(paste("iter: ", iter))
  }
  stopImplicitCluster()
  return(list(beta=beta_new, t=t_new, a=a, b=b, u=u, u2=u2))
}