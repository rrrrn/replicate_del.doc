if(!require(momentfit)){
  install.packages("momentfit")
  library(momentfit)
}

if(!dir.exists("dataset/beta")){
  dir.create("dataset/beta")
}
if(!dir.exists("dataset/t_value")){
  dir.create("dataset/t_value")
}

## calculate local EL estimates as initial values for integrative EL computation
filelist <- list.files("dataset/sim_data", full.names = TRUE)
beta_init <- matrix(0, length(filelist), 5)
t_init <- matrix(0, length(filelist), 5)
for(i in 1:length(filelist)){
  df <- read.csv(filelist[i])
  results <- gel4(y~-1+X1+X2+X3+X4+X5,~-1+X1+X2+X3+X4+X5, data = df, gelType = "EL", vcov = "iid")
  beta_init[i,] <- results@theta
  t_init[i,] <- results@lambda
}

beta_init <- data.frame(beta_init)
colnames(beta_init) <- c("X1", "X2", "X3", "X4", "X5")
t_init <- data.frame(t_init)
colnames(t_init) <- c("X1", "X2", "X3", "X4", "X5")
write.csv(beta_init, "dataset/beta/beta_init.csv", row.names = FALSE)
write.csv(t_init, "dataset/t_value/t_init.csv", row.names = FALSE)


## follow the same procedure, compute centralized EL estimates for all datasets combined
df_integrative <- data.frame()
for(i in 1:length(filelist)){
  df <- read.csv(filelist[i])
  df_integrative <- rbind(df_integrative, df)
}
df_integrative <- data.frame(df_integrative)

results <- gel4(y~-1+X1+X2+X3+X4+X5,~-1+X1+X2+X3+X4+X5, data = df_integrative, gelType = "EL", vcov = "iid")
beta_central <- results@theta
t_central <- results@lambda
beta_central <- t(data.frame(beta_central))
t_central <- t(data.frame(t_central))

write.csv(beta_central, "dataset/beta/beta_central.csv", row.names = FALSE)
write.csv(t_central, "dataset/t_value/t_central.csv", row.names = FALSE)



