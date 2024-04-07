library(rockchalk)
## hyperparameters
d = 5; ## number of covariates
N = 1000 ## number of observations
K = 20 ## number of sites

## covariance matrix
Sigma <- matrix(0, d, d)
for(i in 1:5){
  for(j in 1:5){
    Sigma[i,j] <- Sigma[j,i] <- 0.5^(abs(i-j))
  }
}

## data generation
beta <- c(2, 0.5, 4, sqrt(6), -3)
for(i in 1:K){
  sim_data <- data.frame()
  x <- mvrnorm(N, rep(0, d), Sigma)
  Ep <- rnorm(1, 0, 5)
  epsilon <- rnorm(N, Ep, 1)
  y <- x %*% beta + epsilon
  sim_data <- rbind(sim_data, data.frame(y, x))
  write.csv(sim_data, paste("dataset/sim_data_", i, ".csv", sep = ""), row.names = FALSE)
}
