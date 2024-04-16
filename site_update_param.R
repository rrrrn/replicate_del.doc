library(stats)
library(nloptr)
library(quadprog)
fn <- function(beta, t, Omega, u, x, y, a, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  z <- (y - x %*% beta)[,1] ## n x 1 vector
  h_fn <- sum(log(1-(x*z)%*%t))/N
  p_fn <- 0.5 * t(beta - a + u) %*% Omega %*% (beta - a + u)
  return(h_fn + p_fn)
}

grad_fn <- function(beta, t, Omega, u, x, y, a, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  ## t: d x 1 vector
  z <- (y - x %*% beta)[,1] ## n x 1 matrix
  denom <- (1-(x*z)%*%t) ## n x 1 vector
  num <- diag(x %*% t(x)) ## n x 1 vector
  
  grad_h_fn <- sum(num/denom) * t/N
  grad_p_fn <- Omega %*% (beta - a + u)
  
  return(grad_h_fn + grad_p_fn)
}

fn_t <- function(beta, t, Omega, u2, z, x, y, b, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  h_fn <- sum(log(1-(x*z)%*%t))/N
  p_fn <- 0.5 * t(t - b + u2) %*% Omega %*% (t - b + u2)
  return(-h_fn + p_fn)
}

grad_fn_t <- function(beta, t, Omega, u2, z, x, y, b, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  ## t: d x 1 vector
  n <- dim(x)[1]
  d <- dim(x)[2]
  num = (x*z) ## N x d matrix
  denom <- (1-num%*%t) ## N x 1 vector
  denom <- matrix(rep(denom, d), n, d, byrow = FALSE) ## N x d matrix
  
  grad_h_fn <- apply(num/denom, 2, sum)/N ## d x 1 vector
  grad_p_fn <- Omega %*% (t - b + u2)
  
  return(grad_h_fn + grad_p_fn)
}


site_update_param <- function(beta_init, t_init, Omega1, Omega2, u, u2, x, y, a, b, N){
  ## use optim to find minimizer of beta
  res <- optim(beta_init, fn = function(beta) fn(beta, t = t_init, Omega = Omega1, 
                                                 u = u, x = x, y = y, a = a, N = N), 
               gr = function(beta) grad_fn(beta, t = t_init, Omega = Omega1, 
                                           u = u, x = x, y = y, a = a, N = N), 
               method = "BFGS")
  beta_new <- (res$par)
  
  ## use optim to find minimizer of t
  z = (y - x %*% beta_new)[,1]
  num = x*z ## N x d matrix
  
  eval_g <- function(t){
    return(list("constraints" = num%*%t-1, "jacobian" = num))
  }
  
  ## project t onto feasible set
  Dmat <- diag(dim(x)[2])
  dvec <- t_init
  Amat <- -t(num)
  epsilon <- 1e-6
  bvec <- -rep(1, dim(x)[1]) + epsilon
  t_init <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
  
  res_t <- constrOptim(t_init, f = function(t) fn_t(beta = beta_new, t = t, Omega = Omega2, 
                                                    u2 = u2, z = z, x = x, y = y, b = b, N = N), 
                       ui = -num, ci = -rep(1, dim(x)[1]), grad = function(t) grad_fn_t(beta = beta_new, t = t, 
                                                                                       Omega = Omega2, u2 = u2, z = z, x = x, y = y, b = b, N = N))
  t_new <- (res_t$par)
  return(list(beta = beta_new, t = t_new))
}