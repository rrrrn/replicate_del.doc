site_update_p <- function(p, t, a, u, Omega, x, y, z){
  p0 = p
  optim( p0, # start point
        fn = function(z) { func(y, A, z, lambda) }, # objective function
        gr = function(z) { gradient(y, A, z, lambda) }, # gradient
        method = "BFGS", # specify method
        control=list(maxit=max.iter) )
}