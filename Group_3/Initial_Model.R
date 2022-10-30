#Libraries
library(coda)
library(rjags)
library(ecoforecastR)
library(lubridate)

#Data
chlorophyll <- #Data for chlorophyll a from new site
biomass <- #Data for biomass from new site
#Temp <- #Weather data for temperature from new site
O2 <- #Dissolved oxygen data from new site

data <- list(chlorophyll = y, biomass = x1, O2 = x2, n = length(y))
  
#Model
model1 <- "
model{

  b ~ dmnorm(b0, Vb)
  S ~ dgamma(s1, s2)

  for(i in 1:n){
    mu[i] <- b[1] + b[2]*x1[i] + b[3]*x2[i] ##Process Model
    y[i] ~ dnorm(mu[i], S)    ##Data Model
  }
}
"

#Priors
data$b0 <- as.vector(c(0,0,0)) #Assume mean of each parameter is 0
data$Vb <- solve(diag(10000,3)) #Assume variance of data is 10000
data$s1 <- 0.1
data$s2 <- 0.1

#Initial Conditions
nchain = 3
inits <- list()
for(i in 1:nchain){
  inits[[i]] <- list(b = rnorm(2,0,5), S = runif(1,1/200,1/20))
}
#These three chains have random initial conditions (from the JAGS assignment)
#Idk if there are better choices

#Run JAGS
j.model <- jags.model(file = textConnection(model1),
                      data = data,
                      inits = inits,
                      n.chains = 3)

jags.out <- coda.samples(model = j.model,
                         variable.names = c("b","S"),
                         n.iter = 5000)


#Evaluate Outputs
plot(jags.out)

#Assess Convergence
gelman.diag(jags.out)
#Should be under 1.05 (The closer to 1 the better)

BGR <- gelman.plot(jags.out)
#Remove points from when plot is >1.05

burnin = 500                                   ## determine convergence
#Change burnin based on what the gelman plot shows
jags.burn <- window(jags.out, start = burnin)  ## remove burn-in
plot(jags.burn)                                ## check diagnostics post burn-in

#Go back and make sure you have enough samples for your analysis after burn-in
#Should be >5000

#Check autocorrelation
acfplot(jags.burn)
#Fast decay to 0 indicates stronger independence

effectiveSize(jags.burn)


#MCMC Statistics
summary(jags.burn)
#SD is uncertainty about the parameter
#SE is indicator of precision of results
#(Declines asymptotically, as MCMC length increases)
#Time-Series SE is the one that matters more because it takes into account
#Autocorrelation

out <- as.matrix(jags.burn)
#Pairwise scatter plots & Correlation
pairs(out)	## pairs plot to evaluate parameter correlation
cor(out)    ## correlation matrix among model parameters

