    #words that follow the # symbol are comments and will not be executed
    #THIS PROGRAM FIT A BUCP MODEL WITH UNKNOWN CHANGE-POINTS
install.packages("runjags") #this command installs the runjags package
library(runjags) #loads runjags into the environment
install.packages('rjags')
library(rjags)
setwd("~/Desktop/Research/Trauma-SCED/Plots") #sets the working directory
source("plot_SSD.R") #This is a prewritten program to plot SSD data

    #data is a matrix that contains all data
filename <- 'Au_Figure2_SelfBlame'
data <- as.matrix(read.csv(paste0(filename, ".csv"), header = T))
    #y is the dependent variable with observations in each phase. Each phase is denoted by separate rows
    #P is assigned the value 2, i.e. the number of phases
P <- 2

    #Baseline has Tb observations. This is obtained using the syntax below
Tb <- length(which(data[,1] == 0))

    #Treatment phase has Tt observations which is obtained using the syntax below
Tt <- length(which(data[,1] == 1))

    #T is the total number of observations/time-points
T <- Tb + Tt

    #yall is the DV in the form of a vector
yall <- c(data[,2])
y <- matrix(NA, 2, max(Tb, Tt))
y[1,] <- data[1:Tb, 2]
y[2,] <- data[(Tb+1):T, 2]
plot_SSD(y, paste0(colnames(data)[2], "-SSDplot.jpg"))
    #the following two lines compute the means of the two phases
    #beta1 is the mean of the first Tb observations and beta2 - the rest
    #these betas will be used as starting values in autorunjags
beta1 <- mean(data[1:Tb, 2], na.rm = TRUE)
beta2 <- mean(data[(Tb + 1):T, 2], na.rm = TRUE)

    #Model is defined below:
BUCP.model1 <- "model {
    #this is just a temporary variable to fill in the temp matrix
      temp[1] <- 0 
      
    #yhat for time-point 1 is set to the intercept or the estimated mean of the first phase
    #this is because yhat is the expected value of y and the expected value of y is the mean
      yhat[1] <- beta[1, 1]
      
    #for loop begins at 2 and runs till Tb (for baseline phase)
       for (i in 2:T) {
    # create a dummy variable to indicate the phase
         dummy[i] <- step(CP - i)
    #the expected value of temp in the baseline is beta[1, 1] if it belongs to phase 1
    #it is beta[2, 1] if it belongs to phase 2
         temp[i] <- dummy[i] * beta[1, 1] + (1 - dummy[i]) * beta[2, 1] 
    #yhat is assigned the temp value including the autocorrelation, rho except for the first value in phase 2
    #thus yhat is the expected value at a given time-point
         yhat[i] <- ifelse(CP == i + 1, temp[i], 
                          temp[i] + rho * (y[i - 1] - temp[i - 1]))
    #y is drawn from a distribution with expected value yhat 
    #and autocorrelation rho which is multiplied by the error of the 
    #previous time-point (i - 1)
    #tau is the precision = 1/(standard deviation) = 1/sigma
        y[i] ~ dnorm(yhat[i], tau)
     }
      y[1] ~ dnorm(yhat[1], tau) #equation 1 for baseline
#Prior specifications
#For both phases
      for (i in 1:P){
#the intercepts for baseline and treatment phases are drawn 
#from distributions with means mu[1] and mu[2], respectively 
# and precision 1
        beta[i, 1] ~ dnorm(mu[i], 1)
        
#because mu is the number of vocal responses in five-minute intervals, 
#the expected value is set at 40 with a standard deviation of 20 (precision = 1/20 = .05)
#change the values inside () to reflect your belief or literature
#You can set the lower limit of mu to 0 by adding I(0, ) after dnorm(40, .05)
        mu[i] ~ dnorm(5, .05)
      }
#standard deviation can uniformly vary between 0.1 to 5. 
es <- (beta[2, 1] - beta[1, 1])*sqrt(tau)
#Change values inside the 
#parentheses to reflect your belief or literature
   sigma ~ dunif(0.1, 5)
#tau is sigma^(-2)
   tau <- pow(sigma, -2)
#autocorrelation can vary uniformly between -1 and 1
   rho ~ dunif(-1, 1)
#change-point CP can vary from 3 to T - 3
  CP ~ dcat(pi)
  
}"
############# end of model definition##########################

#Begin running the model with the data

results.BUCP <- autorun.jags(
  #autorun.jags automatically runs the model until convergence is indicated
  model = BUCP.model1, #the model is BITS.model1 defined above
  data = list(y = yall, P = P, T = T, 
              pi = c(rep(0, 3), rep(1/(T - 6), (T - 6)), rep(0, 3))), #input data are the y vector of observations,
  #the total number of time-points T, the number of baseline observations Tb, and phases P
  monitor = c("CP", "beta", "sigma", "rho", "es"), 
  n.chains = 4,#these are the parameters to be monitored i.e., checked for convergence and obtained
  startsample = 30000,
  #burn-in/discard the first 30000 parameter estimates to ensure that 
  #there is no effect of the starting values on the retained estimates
  inits = function() { #initialize/assign starting values
    #change the specs to see how starting values affect the estimates
    #before and after burning-in. Once burned-in they should not.
    list(
      beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
      #intercept starting values around the phase means
      sigma = runif(1, 0.1, 5),
      #standard deviation can be any value between -1 and 5
      rho = runif(1, -1, 1),
      #autocorrelations between -1 and 1
      CP = sample(seq(1:21),1, replace = T)
    )
  }, 
  method = "rjparallel"
  #run the chains in parallel
)


# combine all chains into a single chain for convenience
results.BUCP$draws <- combine.mcmc(results.BUCP$mcmc)

results.BUCP
jpeg(paste0(colnames(data)[2], "-Change-Point.jpg"), quality = 100)
plot(results.BUCP, vars = "CP", plot.type = c("trace", "histogram"))
dev.off()
write.csv(results.BUCP$summaries, paste0(filename, "-BUCPresults.csv"))

