#please download JAGS and R before you run these analyses

#https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
#https://www.r-project.org/
#Rstudio makes it a bit more convenient to run the scripts
#https://rstudio.com/
#words that follow the # symbol are comments and will not be executed
#THIS PROGRAM FIT A BITS MODEL WITH KNOWN CHANGE-POINTS
install.packages("runjags") #this command installs the runjags package

library(runjags) #loads runjags into the environment

install.packages('rjags')
library(rjags)
setwd("~/Desktop/Research/Trauma-SCED/Plots") #sets the working directory
source('plot_SSD.R') #loads the file that plots single case data
source('plots.R') #loads program that plots posterior and ROPE
source('plot_SSD_line.R')

filename <- 'Au_Figure2_SelfBlame' #change the filename for your own analysis on this line

data <- as.matrix(read.csv(paste0(filename, ".csv"), header = T))#y is the dependent variable with observations in each phase. Each phase is denoted by separate rows

#P is assigned the value 2, i.e. the number of phases
P <- 2

#Baseline has Tb observations. This is obtained using the syntax below
Tb <- length(which(data[ ,1] == 0))

#Treatment phase has Tt observations which is obtained using the syntax below
Tt <- length(which(data[ ,1] == 1))

#T is the total number of observations/time-points
T <- Tb + Tt

#create an empty matrix y to fill in baseline data on row 1 and treatment data on row 2
y <- matrix(NA, 2, max(Tb, Tt))
y[1,] <- data[1:Tb, 2]
y[2,] <- data[(Tb+1):T, 2]

#plot and save the single case data against time as a jpeg file
plot_SSD(y, paste0(colnames(data)[2], "-SSDplot.jpg"))

#plot and save the single case data with line of best fit against time as a jpeg file
plot_SSD_line(y, paste0(colnames(data)[2], "-SSDplotline.jpg"))

#the following two lines compute the means of the two phases
#beta1 is the mean of the first Tb observations and beta2 - the rest
#these betas will be used as starting values in autorunjags
beta1 <- rowMeans(y, na.rm = TRUE)[1]
beta2 <- rowMeans(y, na.rm = TRUE)[2]

#JAGS Model is defined below:
BITS.model1 <- "model {
     
      y[1, 1] ~ dnorm(yhat[1, 1], tau) #equation 1 for baseline
      y[2, 1] ~ dnorm(yhat[2, 1], tau) #equation 1 for treatment
#yhat for time-point 1 is set to the intercept or the estimated mean of the first phase
#this is because yhat is the expected value of y and the expected value of y is the mean
    yhat[1, 1] <- beta[1, 1] #equation 4 for baseline
    
#similarly the the first time-point of phase 2
#is assigned to the estimated mean of the second phase
    yhat[2, 1] <- beta[2, 1] #equation 4 for treatment
    
#for loop begins at 2 and runs till Tb (for baseline phase)
     for (i in 2:Tb) {

#y is drawn from a distribution with expected value yhat 
#and autocorrelation rho which is multiplied by the error of the 
#previous time-point (i - 1)
#tau is the precision = 1/(standard deviation) = 1/sigma
#equation 2 for baseline
        y[1, i] ~ dnorm(yhat[1, i] + rho * (y[1, (i - 1)] - yhat[1, (i - 1)]), tau)

#the expected value of yhat in the baseline is beta[1, 1]
        yhat[1, i] <- beta[1, 1] #equation 4 for baseline
     }

#the second for loop begins at 2 and runs till Tt (for treatment phase)
#same syntax applied as baseline except with treatment data
     for (i in 2:Tt) {
#equation 2 for treatment
        y[2, i] ~ dnorm(yhat[2, i] + rho * (y[2, (i - 1)] - yhat[2, (i - 1)]), tau)
     
        yhat[2, i] <- beta[2, 1] #equation 4 for treatment
     }


#Prior specifications for both phases
      for (i in 1:P){
    
#the intercepts for baseline and treatment phases are drawn 
#from distributions with means mu[1] and mu[2], respectively 
# and precision 1
        beta[i, 1] ~ dnorm(mu[i], 1)
        
#because data for this are ranging from 2 to 7
#the expected value is set at 5 with a standard deviation of 20 (precision = 1/20 = .05)
#change the values inside () to reflect your belief or literature
#You can set the lower limit of mu to 0 by adding I(0, ) after dnorm(5, .05) if need be
        mu[i] ~ dnorm(5, .05)
      }
      
#effect size is computed as standardized absolute mean difference
es <- (beta[1, 1] - beta[2, 1])*sqrt(tau)

#standard deviation can uniformly vary between 0.1 to 5. 
#Change values inside the 
#parentheses to reflect your belief or literature
   sigma ~ dunif(0.1, 5)

#tau is sigma^(-2)
   tau <- pow(sigma, -2)

#autocorrelation can vary uniformly between -1 and 1
   rho ~ dunif(-1, 1)
}"
############# end of model definition##########################

#Begin running the model with the data

results <- autorun.jags(
  #autorun.jags automatically runs the model until convergence is indicated
  model = BITS.model1, #the model is BITS.model1 defined above
  data = list(y = y, P = P, Tb = Tb, Tt = Tt), #input data are the y vector of observations,
  #the number of baseline observations Tb, and phases P
  
  monitor = c("beta", "sigma", "rho", "es"), 
  n.chains = 4,#these are the parameters to be monitored i.e., checked for convergence and obtained
  #burn-in/discard the first 30000 parameter estimates to ensure that 
  #there is no effect of the starting values on the retained estimates
  startsample = 30000,
  inits = function() { #initialize/assign starting values
    #change the specs to see how starting values affect the estimates
    #before and after burning-in. Once burned-in they should not.
    list(
      #intercept starting values around the phase means
      beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
      #standard deviation can be any value between -1 and 5
      sigma = runif(1, 0.1, 5),
      #autocorrelations between -1 and 1
      rho = runif(1, -1, 1)
    )
  }, 
  method = "rjparallel"   #run the chains in parallel
)


# combine all chains into a single chain for convenience
results$draws <- combine.mcmc(results$mcmc)

#examine the results
results

#plot the posterior of the effect size 
samples <- combine.mcmc(results$mcmc)
es <- samples[, "es"]

#compVal is the centre of the value around which the ROPE would be drawn
#since this is a test of effect size we will have the upper limit as the default
#the default is any value higher than (compVal - ropeRad)
plots(es, compVal = 1, ropeRad = 0.5, maintitle = "effect size", 
      HDImass = .95, plotname = paste0(colnames(data)[2], "-effect size rope"))

#plot and save the trace and histograms of the posteriors in a jpeg file
jpeg(paste0(colnames(data)[2], "-BITS.jpg"), quality = 100)
plot(results, layout = c(5, 2), plot.type = c("trace", "histogram"))
dev.off()
#write the summary of the posterior distributions to a csv file
write.csv(results$summaries, paste0(filename, "-BITSresults.csv"))
