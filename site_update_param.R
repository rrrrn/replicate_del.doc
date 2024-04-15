library(stats)
fn <- function(beta, t, Omega, u, x, y, a, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  z <- diag((y - x %*% beta)[,1]) ## n x n matrix
  h_fn <- sum(log(1-z%*%(x)%*%t))/N
  p_fn <- 0.5 * t(beta - a + u) %*% Omega %*% (beta - a + u)
  return(h_fn + p_fn)
}

grad_fn <- function(beta, t, Omega, u, x, y, a, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  ## t: d x 1 vector
  z <- diag((y - x %*% beta)[,1]) ## N x N matrix
  denom <- (1-z%*%(x)%*%t) ## N x 1 vector
  num <- diag(x %*% t(x))
  
  grad_h_fn <- sum(num/denom) * t/N
  grad_p_fn <- Omega %*% (beta - a + u)
  
  return(grad_h_fn + grad_p_fn)
}

fn_t <- function(beta, t, Omega, u2, z, x, y, b, N){
  ## x: n x d matrix
  ## y: n x 1 vector
  ## beta: d x 1 vector
  h_fn <- sum(log(1-z%*%(x)%*%t))/N
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
  denom <- (1-z%*%(x)%*%t) ## N x 1 vector
  denom <- matrix(rep(denom, d), n, d, byrow = FALSE) ## N x d matrix
  num <- z%*%(x) ## N x d matrix
  
  grad_h_fn <- apply(num/denom, 2, sum)/N
  grad_p_fn <- Omega %*% (t - b + u2)
  
  return(grad_h_fn + grad_p_fn)
}

box_constraint_t <- function(beta, t, x, y){
  z = diag((y - x %*% beta)[,1])
  num = z%*%(x) ## N x d matrix
  bounds = find_bounds(num, rep(1, dim(num)[1]))
  return(bounds)
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
  z = diag((y - x %*% beta_new)[,1])
  num = z%*%(x) ## N x d matrix
  res_t <- constrOptim(t_init, f = function(t) fn_t(beta = beta_new, t = t, Omega = Omega2, 
                                                     u2 = u2, z = z, x = x, y = y, b = b, N = N), 
                       grad = function(t) grad_fn_t(beta = beta_new, t = t, Omega = Omega2, 
                                                 u2 = u2, z = z, x = x, y = y, b = b, N = N), 
                       ui = -num, ci = -rep(1, dim(num)[1]), method = "BFGS")
  t_new <- (res_t$par)
  return(list(beta = beta_new, t = t_new))
}