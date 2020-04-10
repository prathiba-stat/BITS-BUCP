#words that follow the # symbol are comments and will not be executed
#THIS PROGRAM FIT A BITS MODEL WITH KNOWN CHANGE-POINTS AND LOOPED FOR MULTIPLE SUBJECTS
#Please refer BITS.R for more detailed information about the syntax
install.packages("runjags")
library(runjags)
setwd("~/Desktop/Research/Trauma-SCED/Plots")
source('plot_SSD.R')
source('plots.R')
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
  beta1 <- rowMeans(y, na.rm = TRUE)[1]
  beta2 <- rowMeans(y, na.rm = TRUE)[2]
  
  BITS.model1 <- "model {
    yhat[1, 1] <- beta[1, 1]
    yhat[2, 1] <- beta[2, 1]
     for (i in 2:Tb) {
        yhat[1, i] <- beta[1, 1]
        y[1, i] ~ dnorm(yhat[1, i] + rho * (y[1, (i - 1)] - yhat[1, (i - 1)]), tau)
     }
     for (i in 2:Tt) {
        yhat[2, i] <- beta[2, 1]
        y[2, i] ~ dnorm(yhat[2, i] + rho * (y[2, (i - 1)] - yhat[2, (i - 1)]), tau)
     }
      y[1, 1] ~ dnorm(yhat[1, 1], tau) 
      y[2, 1] ~ dnorm(yhat[2, 1], tau) 
      
      for (i in 1:P){
        beta[i, 1] ~ dnorm(mu[i], 1)
        mu[i] ~ dnorm(5, .05)
      }
    es <- (beta[2, 1] - beta[1, 1])*sqrt(tau)
    sigma ~ dunif(0.1, 5)
    tau <- pow(sigma, -2)
    rho ~ dunif(-1, 1)
}"
  ############# end of model definition##########################
  
  
  results <- autorun.jags(
    model = BITS.model1, 
    data = list(y = y, P = P, Tb = Tb, Tt = Tt), 
    monitor = c("beta", "sigma", "rho", "es"), 
    n.chains = 4,
    startsample = 30000,
    inits = function() { 
      list(
        beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
        sigma = runif(1, 0.1, 5),
        rho = runif(1, -1, 1)
      )
    }, 
    method = "rjparallel"
  )
  results$draws <- combine.mcmc(results$mcmc)
  full.results <- rbind(full.results, 
                        c(colnames(data)[2*i], rep(NA, 10)),
                        results$summaries)
  
  #plot the posterior of the effect size 
  samples <- combine.mcmc(results$mcmc)
  es <- samples[, "es"]
  
  #compVal is the centre of the value around which the ROPE would be drawn
  #since this is a test of effect size we will have the upper limit as the default
  #the default is any value higher than (compVal - ropeRad)
  plots(es, compVal = 1, ropeRad = 0.5, maintitle = "effect size", 
        HDImass = .95, plotname = paste0(colnames(data)[2*i], "-effect size rope"))
  
  jpeg(paste0(colnames(data)[2*i], "-BITS.jpg"), quality = 100)
  plot(results, layout = c(5, 2), plot.type = c("trace", "histogram"))
  dev.off()
  print(i)
}
write.csv(full.results, paste0(filename, "-BITSresults.csv"))

