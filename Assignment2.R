#############################
# Author: Floris Holstege, Markus Mueller
# Goal: Computations for assignment 2 of econometrics I 
############################



# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plm, 
               lmtest, 
               matlib,
               reshape2,
               dplyr)

# Ensure that it is reproducable 
set.seed(1234567)


########
# Question 1
########


# Generate two independent variables, normally distributed
mX1 <- rnorm(100,mean = 0, sd = 1)
mX2 <- rnorm(100,mean = 1, sd = 2)

# Parameters 
a = 0.9
b1 = 0.5
b2 = -0.2
c1 = 0.01
c2 = 0.04

# Dependent variable generated 
y = a + (b1 * mX1) + (b2 * mX2)

# dataset with simulated data
dfSimulated_data <- data.frame(y = numeric(), X1 = numeric(), X2 = numeric())

for(i in 1:100){
  
  # Generate two independent variables, normally distributed
  mX1 <- rnorm(100,mean = 0, sd = 1)
  mX2 <- rnorm(100,mean = 1, sd = 2)
  
  # Dependent variable generated 
  y = a + (b1 * mX1) + (b2 * mX2)
  
  # save for this generated sample
  dfSample = data.frame(y = y, X1 = mX1, X2 = mX2)
  
  # add to overal dataset
  dfSimulated_data = rbind(dfSimulated_data, dfSample)
  
}

lmSimulated <- lm(y ~ X1 + X2, data = dfSimulated_data)

# compute SE with constant var
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="const"))

# compute white SE
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="HC3"))

# White standard errors are slightly smaller, but probably too small of a difference to infer anything. 
# The residuals are homoskedastic (see histogram below, and also logically) so smaller SE not caused by ehtereoskedasticity
hist(lmSimulated$residuals, xlim = c(-1e-13, 1e-13))


# Parameters for WLS and FLS
a0 = 2
b0 = 2
c1 = 0.01
c2 = 0.04
s2 = 0.25

# weights for WLS (given param)
vZ = rgamma(nrow(dfSimulated_data), scale = a0, shape = a0)
vWeights = s2 * exp((c1 * vZ^2) +(c2 * vZ))

# weights for FWLS (est. param)
residualsOLS <- residuals(lmSimulated)
VarEst <- lm(lmSimulated$residuals^2 ~ X1 + X2, data = dfSimulated_data)$fittes.values

# WLS
lmSimulated_WLS <- lm(y ~ X1 + X2, data = dfSimulated_data, weights=vWeights)

# FWLS 
lmSimualted_FWLS <- lm(y ~ X1 + X2, data = dfSimulated_data, weights=VarEst)

# similar SE 
summary(lmSimualted_FWLS)
summary(lmSimulated_WLS)


nSample = 100
nSimulations = 100

# weights for WLS (given param)
vZ = rgamma(nSample, scale = a0, shape = a0)
vWeights = s2 * exp((c1 * vZ^2) +(c2 * vZ))

summary(lmSimulated)$coefficients[,2]


# function to simulate results and compare OLS, WLS, FLS 
simulationResults<- function(nSimulations, nSample, vWeights){
  
  # dataframe with results 
  results <- data.frame(Estimator = character(), 
                        a = numeric(), 
                        b1 = numeric(), 
                        b2 = numeric(),
                        a_se = numeric(),
                        b1 = numeric(),
                        b2 = numeric())
  iIndexResults = 1
  
  
  # create n simulations
  for(i in 1:nSimulations){
    
    # Generate two independent variables, normally distributed
    mX1 <- rnorm(nSample,mean = 0, sd = 1)
    mX2 <- rnorm(nSample,mean = 1, sd = 2)
    
    # Dependent variable generated 
    vY = a + (b1 * mX1) + (b2 * mX2)
    
    # OLS 
    lmOLS <- lm(vY ~ mX1 + mX2)
    SE_OLS <- summary(lmOLS)$coefficients[,2]
    
    # WLS
    lmWLS <- lm(vY ~ mX1 + mX2, weights = vWeights)
    SE_WLS <- summary(lmWLS)$coefficients[,2]
    
    # weights for FWLS (est. param)
    residualsOLS <- residuals(lmOLS)
    VarEst <- lm(lmOLS$residuals^2 ~ mX1 + mX2)$fittes.values
    
    # FWLS 
    lmFWLS <- lm(vY ~ mX1 + mX2, data = dfSimulated_data, weights=VarEst)
    SE_lmFWLS <- summary(lmFWLS)$coefficients[,2]
    
    # add to the df
    results[iIndexResults,-1] <- c(lmOLS$coefficients, SE_OLS)
    results[iIndexResults + 1, -1] <- c(lmWLS$coefficients, SE_WLS)
    results[iIndexResults + 2, -1] <- c(lmFWLS$coefficients, SE_FWLS)
    
    results$Estimator[iIndexResults] <- "OLS"
    results$Estimator[iIndexResults + 1] <- "WLS"
    results$Estimator[iIndexResults + 2] <- "FWLS"
    
    iIndexResults = iIndexResults + 3
    
    
    
    
  }
  
  return(results)
  
}


dfResultSim <- simulationResults(nSimulations, nSample, vWeights)

# results: FWLS and OLS consistent, WLS not because of the gamma distribution is random each time (but only slight divergencese)
resultsMelted <- melt(dfResultSim, id.vars = "Estimator")  





#### 
# Question 4
####
  
