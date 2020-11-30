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
# Question 1A, 1B
########

# Generate two independent variables, normally distributed
mX1 <- rnorm(100,mean = 0, sd = 1)
mX2 <- rnorm(100,mean = 1, sd = 2)

# Parameters 
a = 0.9
b1 = 0.5
b2 = -0.2


# Dependent variable generated 
y = a + (b1 * mX1) + (b2 * mX2)

# dataset with simulated data
dfSimulated_data <- data.frame(y = numeric(), X1 = numeric(), X2 = numeric())

# generate 100 samples of 100, according to the model
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

# test the standard OLS on the simualted data
lmSimulated <- lm(y ~ X1 + X2, data = dfSimulated_data)

# compute SE with constant var
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="const"))

# compute white SE
coeftest(lmSimulated, vcov = vcovHC(lmSimulated, type="HC3"))


# White standard errors are slightly smaller, but probably too small of a difference to infer anything. 
# The residuals are homoskedastic (see histogram below, and also logically) so smaller SE not caused by heteroskedasticity
hist(lmSimulated$residuals, xlim = c(-1e-13, 1e-13), main = "Residuals with OLS", xlab="Residuals")



########
# Question 1C
########

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

# get residuals from OLS estimate
residualsOLS <- residuals(lmSimulated)

# regress these with Z = (X1, X2)
VarEst <- lm(residualsOLS^2 ~ X1 + X2 ,  data = dfSimulated_data)
summary(VarEst)$coefficients
mean(VarEst$fitted.values)

# WLS
lmSimulated_WLS <- lm(y ~ X1+ X2, data = dfSimulated_data, weights=vWeights)


# FWLS 
lmSimualted_FWLS <- lm(y ~ X1 + X2, data = dfSimulated_data, weights=VarEst$fitted.values)

# similar SE 
summary(lmSimulated_WLS)$coefficient
summary(lmSimualted_FWLS)$coefficient


hist(lmSimulated_WLS$residuals, xlim = c(-1e-13, 1e-13), main = "Residuals with WLS", xlab = "Residuals")
hist(lmSimualted_FWLS$residuals, xlim = c(-1e-13, 1e-13), main = "Residuals with FWLS", xlab = "Residuals")

nSample = 100
nSimulations = 100

# weights for WLS (given param)
vZ = rgamma(nSample, scale = a0, shape = a0)
vWeights = s2 * exp((c1 * vZ^2) +(c2 * vZ))



# function to simulate results and compare OLS, WLS, FLS 
simulationResults<- function(nSimulations, nSample, vWeights){
  
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


dfResultSim <- simulationResults(nSimulations, nSample, vWeights)

# results: FWLS and OLS consistent, WLS not because of the gamma distribution is random each time (but only slight divergencese)
resultsMelted <- melt(dfResultSim, id.vars = "Estimator")  

# gather results for SE
resultsSE <- resultsMelted %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# summarize for mean in plots
resultsSE_summarized <- resultsSE %>%
  group_by(Estimator, variable) %>%
  summarize(avg = mean(value))

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
  summarize(avg = mean(value), sd = sd(value))


### under alternative parameters
c1_0 = 0
c2_0 = 0

vWeights_0 = s2 * exp((c1_0 * vZ^2) +(c2_0 * vZ))


dfResultSim_0 <- simulationResults(nSimulations, nSample, vWeights_0)
resultsMelted_0 <-  melt(dfResultSim_0, id.vars = "Estimator")

# gather results for SE
resultsSE_0 <- resultsMelted_0 %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# summarize for mean in plots
resultsSE_summarized_0 <- resultsSE_0 %>%
  group_by(Estimator, variable) %>%
  summarize(avg = mean(value))

# create lots to compare the standard errors
ggplot(resultsSE_0, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 8, position="identity")+
  geom_vline( data = resultsSE_summarized, mapping = aes(xintercept = avg))+
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
stargazer(lm_Famincome_size)
# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -3179 in family income

# Independent = family size, dependent = weeks worked
lm_worked_size <- lm(weeksworked ~ familysize, data = dfFamily)
summary(lm_worked_size)
stargazer(lm_worked_size)
# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -4 weeks worked less in the year for the mother

# Independent = family size, dependent = hours per week
lm_hours_size <- lm(hoursperweek ~ familysize,data = dfFamily)
summary(lm_hours_size)
stargazer(lm_hours_size)
# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -2 hours worked less per week

# Independent = family size, dependent = labor income of the mother
lm_Labincome_size <- lm( laborincome ~ familysize,data = dfFamily)
summary(lm_Labincome_size)
stargazer(lm_Labincome_size)
# Statistically significant at 5% - one unit of change in family size (e.g. one more or less family member) leads on average to -2889 in the labor income of the mother




###################
# Question 4B, C
###################

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
summary(lm_worked_twinIV, diagnostics=TRUE)

# re-restimate models with the age child instrument variable (explained part)
lm_Famincome_ageChildIV <- ivreg(familyincome ~ familysize | ageChild, data = dfFamily)
summary(lm_Famincome_ageChildIV, diagnostics=TRUE)

lm_worked_ageChildIV <- ivreg(weeksworked ~ familysize | ageChild, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)

lm_hours_ageChildIV <- ivreg(hoursperweek ~ familysize | ageChild, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)

lm_Labincome_ageChildIV <- ivreg(laborincome ~ familysize | ageChild, data = dfFamily)
summary(lm_worked_ageChildIV, diagnostics=TRUE)



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



