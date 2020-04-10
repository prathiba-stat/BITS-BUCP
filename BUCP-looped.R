    #words that follow the # symbol are comments and will not be executed
    #THIS PROGRAM FIT A BUCP MODEL WITH UNKNOWN CHANGE-POINTS for multiple subjects
    #Please refer BUCP.R for detailed comments on the syntax
install.packages("runjags")
library(runjags)
setwd("~/Desktop/Research/Trauma-SCED/Plots")
source('plot_SSD.R')
filename <- 'Au_Figure2_SelfBlame'
data <- as.matrix(read.csv(paste0(filename, ".csv"), header = T))
full.results <- c()
P <- 2
for (i in 1:(ncol(data)/2)){
  Tb <- length(which(data[ ,(2*i - 1)] == 0))
  Tt <- length(which(data[ ,(2*i - 1)] == 1))
  T <- Tb + Tt
  y <- matrix(NA, 2, max(Tb, Tt))
  y[1,] <- data[1:Tb, 2*i]
  y[2,] <- data[(Tb+1):T, 2*i]
  plot_SSD(y, paste0(colnames(data)[2*i], "-SSDplot.jpg"))
  yall <- c(data[,2])
  beta1 <- mean(data[1:Tb, 2*i], na.rm = TRUE)
  beta2 <- mean(data[(Tb + 1):T, 2*i], na.rm = TRUE)
  
  BUCP.model1 <- "model {
      temp[1] <- 0 
      yhat[1] <- beta[1, 1]
       for (i in 2:T) {
         dummy[i] <- step(CP - i)
         temp[i] <- dummy[i] * beta[1, 1] + (1 - dummy[i]) * beta[2, 1] 
         yhat[i] <- ifelse(CP == i + 1, temp[i], 
                          temp[i] + rho * (y[i - 1] - temp[i - 1]))
        y[i] ~ dnorm(yhat[i], tau)
     }
      y[1] ~ dnorm(yhat[1], tau) #equation 1 for baseline
      for (i in 1:P){
        beta[i, 1] ~ dnorm(mu[i], 1)
        mu[i] ~ dnorm(5, .05)
      }
   es <- (beta[2, 1] - beta[1, 1])*sqrt(tau)
   sigma ~ dunif(0.1, 5)
   tau <- pow(sigma, -2)
   rho ~ dunif(-1, 1)
   CP ~ dcat(pi)
  
}"
############# end of model definition##########################


  results.BUCP <- autorun.jags(
    model = BUCP.model1, #the model is BITS.model1 defined above
    data = list(y = yall, P = P, T = T, 
                pi = c(rep(0, 3), rep(1/(T - 6), (T - 6)), rep(0, 3))), 
    monitor = c("CP", "beta", "sigma", "rho", "es"), 
    n.chains = 4,
    startsample = 30000,
    inits = function() { 
      list(
        beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
        sigma = runif(1, 0.1, 5),
        rho = runif(1, -1, 1),
        CP = sample(seq(1:21),1, replace = T)
      )
    }, 
    method = "rjparallel"
  )
  
  
  results.BUCP$draws <- combine.mcmc(results.BUCP$mcmc)
  
  results.BUCP
  full.results <- rbind(full.results, 
                        c(colnames(data)[2*i], rep(NA, 10)),
                        results.BUCP$summaries)
  jpeg(paste0(colnames(data)[2*i], "-Change-Point.jpg"), quality = 100)
  plot(results.BUCP, vars = "CP", plot.type = c("trace", "histogram"))
  dev.off()
  print(i)
}
write.csv(full.results, paste0(filename, "-BUCPresults.csv"))

