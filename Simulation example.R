setwd("~/Desktop/Research/Trauma-SCED/Plots")
source("plot_SSD.R")
library(gsarima)
baseline.mean<-60
difference <- -10
n <- 8
ar <- 0.2
T <- 2*n
base.X=matrix(c(rep(baseline.mean, n+length(ar))), ncol=1)
ybase.sim <- garsim(n=(n+length(ar)), phi=ar, beta=c(1), link= "identity",
                    family= "gaussian", X=base.X)
ybase<-ybase.sim[(1+length(ar)):(n+length(ar))]
#acf(tsy, 1, plot = FALSE)
treatment.mean<-baseline.mean + difference
treat.X=matrix(c(rep(treatment.mean, n+length(ar))), ncol=1)
ytreat.sim <- garsim(n=(n+length(ar)), phi=ar, beta=c(1), link= "identity",
                     family= "gaussian", X=treat.X)
ytreat<-ytreat.sim[(1+length(ar)):(n+length(ar))]
data <- c(ybase, ytreat)
plot_SSD(rbind(ybase, ytreat), "simulated.jpg")
data <- cbind((1:(2*n)), data)
beta1 <- mean(data[1:n, 2], na.rm = TRUE)
beta2 <- mean(data[(n + 1):T, 2], na.rm = TRUE)
yall <- c(data[,2])
P <- 2

#Model is defined below:
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
        mu[i] ~ dnorm(50, .1)
      }
es <- (beta[1, 1] - beta[2, 1])*sqrt(tau)
   sigma ~ dgamma(1, 1)
   tau <- pow(sigma, -2)
   rho ~ dunif(-1, 1)
  CP ~ dcat(pi)
  
}"
############# end of model definition##########################

#Begin running the model with the data

results.BUCP <- autorun.jags(
  model = BUCP.model1, 
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
jpeg(paste0(colnames(data)[2], "simulation-gamma-Change-Point.jpg"), quality = 100)
plot(results.BUCP, vars = "CP", plot.type = c("trace", "histogram"))
dev.off()
write.csv(results.BUCP$summaries, "simulationresults-gamma.csv")

