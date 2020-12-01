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
               dplyr,
               estimatr,
               ggplot2,
               scales,
               Hmisc, 
               GGally,
               stargazer,
               AER)

# Ensure that it is reproducable 
set.seed(1234567)


########
# Question 1A
########

# Generate two independent variables, normally distributed
mX1 <- rnorm(100,mean = 1, sd = 2)
mX2 <- rnorm(100,mean = 0, sd = 1)

# Parameters 
a = 0.9
b1 = 0.5
b2 = -0.2

# Parameters for sigma2_i
a0 = 2
b0 = 2
c1 = 0.01
c2 = 0.04
sigma2 = 0.25


########
# Question 1B
########


# dataset where simulated data will be stored
dfSimulated_data <- data.frame(y = numeric(), X1 = numeric(), X2 = numeric())

vErrors <- c()
vSigma2_i <- c()
vZ_saved <- c()

# generate 100 samples of 100, according to the model
for(i in 1:100){
  
  # Generate two independent variables, normally distributed
  mX1 <- rnorm(100,mean = 0, sd = 1)
  mX2 <- rnorm(100,mean = 1, sd = 2)
  
  # Dependent variable generated 
  y = a + (b1 * mX1) + (b2 * mX2)
  
  # create errors
  # talked to walter - I think we need to use sigma2_i here, otherwise no heteroskadisticity 
  vZ = rgamma(100, scale = a0, shape = a0)
  sigma2_i = sigma2 * exp((c1* vZ^2)+ (c2 * vZ))
  errors = rnorm(100,0,sd = sigma2_i)
  
  # Dependent variable generated 
  y = a + (b1 * mX1) + (b2 * mX2) + errors
  
  vErrors <- c(vErrors, errors)
  vSigma2_i <- c(vSigma2_i, sigma2_i)
  vZ_saved <- c(vZ_saved, vZ)
  
  # save for this generated sample
  dfSample = data.frame(y = y, X1 = mX1, X2 = mX2)
  
  # add to overall dataset
  dfSimulated_data = rbind(dfSimulated_data, dfSample)
  
}

#
plot(dfSimulated_data$X1,vErrors, main = "X1 and the errors", ylim = c(-10,10), xlab = "X1")
plot(dfSimulated_data$X2, vErrors, main = "X2 and the errors", ylim = c(-10,10), xlab = "X2")

# test the standard OLS on the simualted data
lmSimulated <- lm(y ~ X1 + X2, data = dfSimulated_data)

# compute SE with constant var
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="const"))

# compute white SE
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="HC1"))

# White standard errors are slightly bigger, but probably too small of a difference to infer anything. 
# The residuals are normal (see histogram below, and also logically) so smaller SE not caused by heteroskedasticity
hist(lmSimulated$residuals, main = "Residuals with OLS", xlab="Residuals", breaks = 1000, xlim = c(-10,10))


########
# Question 1C
########

# weight for WLS
WLS_weight <- exp((vZ_saved^2 *c1 )+ (vZ_saved * c2))

# step 3: apply these weights in WLS
lmSimulated_WLS <- lm(y ~ X1+ X2, data = dfSimulated_data, weights=1/vSigma2_i)

## For the FWLS, we use the residuals of the OLS to estimate s2, c1, and c2. 

# get residuals from OLS estimate
residualsOLS <- residuals(lmSimulated)

# model the variables in a dataframe( intercept = sigma2)
model_sigma2_i <- data.frame( vZ_saved^2, vZ_saved)
colnames(model_sigma2_i) <- c("vZ.2", "vZ")

# Take the log of the residuals squared, and see how well these 
VarEst <- lm(log(residualsOLS^2) ~ vZ.2 + vZ ,data = model_sigma2_i )
c1_est <- VarEst$coefficients[2]
c2_est <- VarEst$coefficients[3]

FWLS_weight <- exp((vZ_saved^2*c1_est) + (vZ_saved * c2_est))

# FWLS - uses est. sigma2_ i in weights
lmSimulated_FWLS <- lm(y ~ X1 + X2, data = dfSimulated_data, weights =1/exp(VarEst$fitted.values))

# similar SE - slightly bigger for the FWLS
summary(lmSimulated_WLS)
summary(lmSimulated_FWLS)

########
# Question 1D
########


# function to simulate results and compare OLS, WLS, FLS 
simulationResults<- function(nSimulations, nSample, a, b1,b2,c1,c2, a0,b0,sigma2){
  
  # dataframe with results 
  results <- data.frame(Estimator = character(), 
                        a = numeric(), 
                        b1 = numeric(), 
                        b2 = numeric(),
                        a_se = numeric(),
                        b1_se = numeric(),
                        b2_se = numeric())
  iIndexResults = 1
  
  
  # create n simulations
  for(i in 1:nSimulations){
    
    # Generate two independent variables, normally distributed
    mX1 <- rnorm(nSample,mean = 0, sd = 1)
    mX2 <- rnorm(nSample,mean = 1, sd = 2)
    
    # create errors
    vZ = rgamma(nSample, scale = a0, shape = b0)
    sigma2_i = sigma2 * exp((c1* vZ^2)+ (c2 * vZ))
    errors = rnorm(nSample,0,sd = sigma2_i)
    
    # add to y
    # Dependent variable generated 
    vY = a + (b1 * mX1) + (b2 * mX2) + errors
    
    # OLS 
    lmOLS <- lm(vY ~ mX1 + mX2)
    SE_OLS <- summary(lmOLS)$coefficients[,2]
    
    # weight for WLS
    WLS_weight <- exp((vZ^2 *c1 )+ (vZ * c2))
    
    # for WLS: repeat previous steps, take sigma2_i by which the errors are generated
    lmWLS <- lm(vY ~ mX1 + mX2, weights=1/WLS_weight)
    SE_WLS <- summary(lmWLS)$coefficients[,2]
    
    # for FWLS: again, repeat previous steps
    # model the variables in a dataframe( intercept = sigma2)
    model_sigma2_i <- data.frame( vZ^2, vZ)
    
    # get residuals from OLS model
    residualsOLS <- residuals(lmOLS)
    
    # Take the log of the residuals squared, and see how well these 
    VarEst <- lm(log(residualsOLS^2) ~ vZ.2 + vZ ,data = model_sigma2_i )
    c1_est <- VarEst$coefficients[2]
    c2_est <- VarEst$coefficients[3]
    
    FWLS_weight <- exp((vZ^2*c1_est) + (vZ * c2_est))
    
    # FWLS - uses est. sigma2_ i in weights
    lmFWLS <- lm(vY ~ mX1 + mX2, weights = 1/FWLS_weight)
    SE_FWLS <- summary(lmFWLS)$coefficients[,2]

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

###
# Question 1D
#
###

# 100 simulations of a sample of 100
nSample = 10000
nSimulations = 100

# simulate results
dfResultSim <- simulationResults(nSimulations, nSample, a, b1,b2,c1,c2,a0,b0, sigma2)

# results: FWLS and OLS consistent, WLS not because of the gamma distribution is random each time (but only slight divergencese)
resultsMelted <- melt(dfResultSim, id.vars = "Estimator")  

# gather results for SE
resultsSE <- resultsMelted %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# summarize for mean in plots
resultsSE_summarized <- resultsSE %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value))
resultsSE_summarized

# create lots to compare the standard errors
ggplot(resultsSE, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsSE_summarized, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")

# gather results for coefficients
resultsParam <- resultsMelted %>% filter(variable %in% c("a", "b1", "b2"))

# check in dataframe
resultsParam_summarized <- resultsParam %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value), sd = sd(value))
resultsParam_summarized

# create lots to compare the coefficients
ggplot(resultsParam, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsParam_summarized, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")


###
# Question 1E
#
###

### under alternative parameters
c1_0 = 0
c2_0 = 0

dfResultSim_0 <- simulationResults(nSimulations, nSample, a, b1,b2,c1_0,c2_0,a0,b0, sigma2)
resultsMelted_0 <-  melt(dfResultSim_0, id.vars = "Estimator")

# gather results for SE
resultsSE_0 <- resultsMelted_0 %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# summarize for mean in plots
resultsSE_summarized_0 <- resultsSE_0 %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value))
resultsSE_summarized_0

# create lots to compare the standard errors
ggplot(resultsSE_0, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsSE_summarized_0, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")

# gather results for coefficients
resultsParam_0 <- resultsMelted_0 %>% filter(variable %in% c("a", "b1", "b2"))

# check in dataframe
resultsParam_summarized_0 <- resultsParam_0 %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value), sd = sd(value))
resultsParam_summarized_0

# create lots to compare the coefficients
ggplot(resultsParam_0, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsParam_summarized_0, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")



###################
# Question 4A
###################


dfFamily <- read.csv("Data/DataAS21.csv")

# Independent = family size, dependent = family income
lm_Famincome_size <- lm(familyincome ~ familysize,data = dfFamily)
summary(lm_Famincome_size) 

# Independent = family size, dependent = weeks worked
lm_worked_size <- lm(weeksworked ~ familysize, data = dfFamily)
summary(lm_worked_size)

# Independent = family size, dependent = hours per week
lm_hours_size <- lm(hoursperweek ~ familysize,data = dfFamily)
summary(lm_hours_size)

# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -2 hours worked less per week

# Independent = family size, dependent = labor income of the mother
lm_Labincome_size <- lm( laborincome ~ familysize,data = dfFamily)
summary(lm_Labincome_size)

# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -2889 in the labor income of the mother

# output table for latex document
stargazer(lm_Famincome_size, lm_worked_size, lm_hours_size, lm_Labincome_size)


###################
# Question 4B, C
###################
dfFamily <- read.csv("DataAS21.csv")

# Correlations to check if instruments are related to the other independent variables
# Twin variable seems only mildly (0.1) correlated with family size, age child is slightly better (0.3). Same sex and education seem both unrelated to family size(0)
mCor <- round(cor(dfFamily[,-1]),2)
ggcorr(dfFamily[,-1], nbreaks=8, palette='RdGy', label=TRUE, label_size=5, label_color='white')

# get residuals for the regressions
residuals_Famincome_size <- residuals(lm_Famincome_size)
residuals_worked_size <- residuals(lm_worked_size)
residuals_hours_size <- residuals(lm_hours_size)
residuals_Labincome_size <- residuals(lm_Labincome_size)


# check correlation of residuals to any of the potential instruments
dfResiduals <- data.frame(residuals_Famincome_size, residuals_worked_size, residuals_hours_size, residuals_Labincome_size)
dfVar_Residuals <- cbind(dfFamily, dfResiduals)
mCor_var_residuals <-round(cor(dfVar_Residuals), 2)
ggcorr(dfVar_Residuals[,-1], nbreaks=8, palette='RdGy', label=TRUE, label_size=5, label_color='white')
# It does nto seem like either twin or age child has a problematic correlation with the residuals of any of the regressions


# two stage least squares with instruments twin and age child
lm_twin_familySize <- lm(familysize ~ twin, data = dfFamily)
lm_ageChild_familySize <- lm(familysize ~agechild, data = dfFamily)

explained_familySize_twin <- lm_twin_familySize$fitted.values
explained_familySize_ageChild <- lm_ageChild_familySize$fitted.values


#re-estimate models with the twin instrument variable (explained part)
lm_Famincome_twinIV <- ivreg(familyincome ~ familysize | twin, data = dfFamily)
summary(lm_Famincome_twinIV, diagnostics=TRUE)

lm_worked_twinIV <- ivreg(weeksworked ~ familysize | twin, data = dfFamily)
summary(lm_worked_twinIV, diagnostics=TRUE)

lm_hours_twinIV <- ivreg(hoursperweek ~ familysize | twin, data = dfFamily)
summary(lm_worked_twinIV, diagnostics=TRUE)

lm_Labincome_twinIV <- ivreg(laborincome ~ familysize | twin, data = dfFamily)
summary(lm_Labincome_twinIV, diagnostics=TRUE)

colnames(dfFamily)

# output table for latex document
stargazer(lm_Famincome_twinIV, lm_worked_twinIV, lm_hours_twinIV, lm_Labincome_twinIV)

# re-restimate models with the age child instrument variable (explained part)
lm_Famincome_ageChildIV <- ivreg(familyincome ~ familysize | samesex, data = dfFamily)
summary(lm_Famincome_ageChildIV, diagnostics=TRUE)

lm_worked_ageChildIV <- ivreg(weeksworked ~ familysize | samesex, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)

lm_hours_ageChildIV <- ivreg(hoursperweek ~ familysize | samesex, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)

lm_Labincome_ageChildIV <- ivreg(laborincome ~ familysize | samesex, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)

# output table for latex document
stargazer(lm_Famincome_ageChildIV, lm_worked_ageChildIV, lm_hours_ageChildIV, lm_Labincome_ageChildIV)


##### alternative code - shows what ivreg does 
#lm_Famincome_twin <- lm(familyincome ~ explained_familySize_twin, data = dfFamily)
#summary(lm_Famincome_twin, diagnostics=TRUE)


#lm_worked_twin <- lm(weeksworked ~ explained_familySize_twin, data = dfFamily)
#summary(lm_worked_twin, diagnostics=TRUE)
# 
# lm_hours_twin <- lm(hoursperweek ~ explained_familySize_twin, data = dfFamily)
# summary(lm_hours_twin, diagnostics=TRUE)
# 
# lm_Labincome_twin <- lm(laborincome ~ explained_familySize_twin, data = dfFamily)
# summary(lm_Labincome_twin, diagnostics=TRUE)



