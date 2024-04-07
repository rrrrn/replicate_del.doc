if(!require(momentfit)){
  install.packages("momentfit")
  library(momentfit)
}

## calculate local EL estimates as initial values for integrative EL computation
filelist <- list.files("dataset", full.names = TRUE)
beta_init <- matrix(0, length(filelist), 5)
for(i in 1:length(filelist)){
  df <- read.csv(filelist[i])
  results <- gel4(y~-1+X1+X2+X3+X4+X5,~X1+X2+X3+X4+X5, data = df, gelType = "EL", vcov = "iid")
  beta_init[i,] <- coef(results)
}

beta_init <- data.frame(beta_init)
colnames(beta_init) <- c("X1", "X2", "X3", "X4", "X5")
write.csv(beta_init, "dataset/beta/beta_init.csv", row.names = FALSE)


## follow the same procedure, compute centralized EL estimates for all datasets combined
df_integrative <- data.frame()
for(i in 1:length(filelist)){
  df <- read.csv(filelist[i])
  df_integrative <- rbind(df_integrative, df)
}
df_integrative <- data.frame(df_integrative)

results <- gel4(y~-1+X1+X2+X3+X4+X5,~X1+X2+X3+X4+X5, data = df_integrative, gelType = "EL", vcov = "iid")
beta_central <- coef(results)
beta_central <- t(data.frame(beta_central))

write.csv(beta_central, "dataset/beta/beta_central.csv", row.names = FALSE)


