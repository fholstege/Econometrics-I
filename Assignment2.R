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
# Parameters 
########

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

# Alternative C1, C2
c1_0 = 0
c2_0 = 0

####
# Key function: simulates results given parameters
#####

# function to simulate results and compare OLS, WLS, FLS 
simulationResults<- function(nSimulations, nSample, a, b1,b2,c1,c2, a0,b0,sigma2){
  
  # dataframe with results 
  results <- data.frame(Estimator = character(), 
                        a = numeric(), 
                        b1 = numeric(), 
                        b2 = numeric(),
                        a_se = numeric(),
                        b1_se = numeric(),
                        b2_se = numeric(), 
                        a_white_se = numeric(),
                        b1_white_se = numeric(),
                        b2_white_se = numeric())
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
    White_SE_OLS <- coeftest(lmOLS, vcov = vcovHC(lmOLS, type="HC1"))[,2]
    
    # weight for WLS
    WLS_weight <- exp((vZ^2 *c1 )+ (vZ * c2))
    
    # for WLS: repeat previous steps, take sigma2_i by which the errors are generated
    lmWLS <- lm(vY ~ mX1 + mX2, weights=1/WLS_weight)
    SE_WLS <- summary(lmWLS)$coefficients[,2]
    White_SE_WLS <- coeftest(lmWLS, vcov = vcovHC(lmWLS, type="HC1"))[,2]
    
    
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
    White_SE_FWLS <- coeftest(lmWLS, vcov = vcovHC(lmWLS, type="HC1"))[,2]

    # add to the df
    results[iIndexResults,-1] <- c(lmOLS$coefficients, SE_OLS, White_SE_OLS)
    results[iIndexResults + 1, -1] <- c(lmWLS$coefficients, SE_WLS, White_SE_WLS)
    results[iIndexResults + 2, -1] <- c(lmFWLS$coefficients, SE_FWLS, White_SE_FWLS)
    
    results$Estimator[iIndexResults] <- "OLS"
    results$Estimator[iIndexResults + 1] <- "WLS"
    results$Estimator[iIndexResults + 2] <- "FWLS"
    
    iIndexResults = iIndexResults + 3
    
    
  }
  
  return(results)
  
}

###
# Results per question
###

# 100 simulations of a sample of 100
nSample = 100
nSimulations = 100
dfResultSim <- simulationResults(nSimulations, nSample, a, b1,b2,c1,c2,a0,b0, sigma2)
resultsMelted <- melt(dfResultSim, id.vars = "Estimator")  

# gather results for SE
resultsSE <- resultsMelted %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# gather results for coefficients
resultsParam <- resultsMelted %>% filter(variable %in% c("a", "b1", "b2"))

# summarize results for SE 
resultsSE_summarized <- resultsSE %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value))
xtable(resultsSE_summarized %>% arrange(variable))

# Summarize results for Param
resultsParam_summarized <- resultsParam %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value), sd = sd(value))
xtable(resultsParam_summarized %>% arrange(variable))

# 1B: Compare white standard errors and constant variance standard errors in OLS
resultsCompareSE <- resultsMelted %>% filter(variable %in% c("a_white_se", "b1_white_se", "b2_white_se", "a_se", "b1_se", "b2_se"), Estimator %in% "OLS")
resultsCompareSE

# summarise in table
resultsCompareSE_summarized <- resultsCompareSE %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value))
resultsCompareSE_summarized
xtable(resultsCompareSE_summarized %>% arrange(variable))


# 1C:  Comparison of WLS and FWLS
Param_WLS_FWLS <- resultsParam_summarized %>% filter(Estimator %in% c("WLS", "FWLS"))
SE_WLS_FWLS <- resultsSE_summarized %>% filter(Estimator %in% c("WLS", "FWLS"))
colnames(SE_WLS_FWLS)[3] <- "avg_SE"
WLS_FWLS_CompareTable <- cbind(Param_WLS_FWLS, SE_WLS_FWLS[,ncol(SE_WLS_FWLS)])
xtable(WLS_FWLS_CompareTable)

# 1D: Create plots for SE and parameter comparison
# create lots to compare the standard errors
ggplot(resultsSE, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsSE_summarized, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")

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
# Repeat above analysis under different parameters
###

dfResultSim_0 <- simulationResults(nSimulations, nSample, a, b1,b2,c1_0,c2_0,a0,b0, sigma2)
resultsMelted_0 <-  melt(dfResultSim_0, id.vars = "Estimator")

# gather results for SE
resultsSE_0 <- resultsMelted_0 %>% filter(variable %in% c("a_se", "b1_se", "b2_se"))

# summarize SE results 
resultsSE_summarized_0 <- resultsSE_0 %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value))
xtable(resultsSE_summarized_0 %>% arrange(variable))

# create lots to compare the standard errors
ggplot(resultsSE_0, aes(x = value, fill = Estimator))+
  geom_histogram(alpha = 0.5, bins = 10, position="identity")+
  geom_vline( data = resultsSE_summarized_0, mapping = aes(xintercept = avg))+
  facet_wrap(~Estimator + variable, scale = "free_x")+
  scale_x_continuous(breaks =  pretty_breaks())+
  theme_bw()+
  labs(y = "Count", x = "Value")

# gather results for parameters
resultsParam_0 <- resultsMelted_0 %>% filter(variable %in% c("a", "b1", "b2"))

# summarize Param results 
resultsParam_summarized_0 <- resultsParam_0 %>%
  group_by(Estimator, variable) %>%
  summarise(avg = mean(value), sd = sd(value))
xtable(resultsParam_summarized_0 %>% arrange(variable))

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

# load data
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

# Independent = family size, dependent = labor income of the mother
lm_Labincome_size <- lm( laborincome ~ familysize,data = dfFamily)
summary(lm_Labincome_size)

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


#re-estimate models with the twin instrument variable (explained part)
lm_Famincome_twinIV <- ivreg(familyincome ~ familysize | twin, data = dfFamily)
summary(lm_Famincome_twinIV, diagnostics=TRUE)

lm_worked_twinIV <- ivreg(weeksworked ~ familysize | twin, data = dfFamily)
summary(lm_worked_twinIV, diagnostics=TRUE)

lm_hours_twinIV <- ivreg(hoursperweek ~ familysize | twin, data = dfFamily)
summary(lm_hours_twinIV, diagnostics=TRUE)

lm_Labincome_twinIV <- ivreg(laborincome ~ familysize | twin, data = dfFamily)
summary(lm_Labincome_twinIV, diagnostics=TRUE)

# output table for latex document
stargazer(lm_Famincome_twinIV, lm_worked_twinIV, lm_hours_twinIV, lm_Labincome_twinIV)

# re-restimate models with the age child instrument variable (explained part)
lm_Famincome_samesexIV <- ivreg(familyincome ~ familysize | samesex, data = dfFamily)
summary(lm_Famincome_samesexIV, diagnostics=TRUE)

lm_worked_samesexIV <- ivreg(weeksworked ~ familysize | samesex, data = dfFamily)
summary(lm_worked_samesexIV , diagnostics=TRUE)

lm_hours_samesexIV <- ivreg(hoursperweek ~ familysize | samesex, data = dfFamily)
summary(lm_hours_samesexIV, diagnostics=TRUE)

lm_Labincome_samesexIV <- ivreg(laborincome ~ familysize | samesex, data = dfFamily)
summary(lm_Labincome_samesexIV, diagnostics=TRUE)

# output table for latex document
stargazer(lm_Famincome_samesexIV,lm_worked_samesexIV,lm_hours_samesexIV,lm_Labincome_samesexIV)


###################
# Question 4E - Hausman test for each model
###################
######################################
# Z: twin, Y: familyincome
######################################
# step (1) regress y on x
u1 <- resid(lm_Famincome_size)

# step (2) regress endogenous x on z
Res2 <- lm(familysize ~ twin, dfFamily)
u2 <- resid(Res2)

# step (3) regress residual of step 1 on x and residual of step 2
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)

n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared

# test statistic
H <- n3*Rsq

# pvalue
1-pchisq(H,1) 

######################################
# Z: twin, Y: weeks worked
######################################
u1 <- resid(lm_worked_size)
Res2 <- lm(familysize ~ twin, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: twin, Y: hours worked
######################################
u1 <- resid(lm_hours_size)
Res2 <- lm(familysize ~ twin, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: twin, Y: mothers income
######################################
u1 <- resid(lm_Labincome_size)
Res2 <- lm(familysize ~ twin, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: samesex, Y: familyincome
######################################
u1 <- resid(lm_Famincome_size)
Res2 <- lm(familysize ~ samesex, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: samesex, Y: weeks worked
######################################
u1 <- resid(lm_worked_size)
Res2 <- lm(familysize ~ samesex, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: samesex, Y: hours worked
######################################
u1 <- resid(lm_hours_size)
Res2 <- lm(familysize ~ samesex, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 

######################################
# Z: samesex, Y: mothers income
######################################
u1 <- resid(lm_Labincome_size)
Res2 <- lm(familysize ~ samesex, dfFamily)
u2 <- resid(Res2)
Res3 <- lm(u1 ~ familysize + u2, dfFamily)
summary(Res3)
n3 <- nobs(Res3)
Rsq <- summary(Res3)$r.squared
H <- n3*Rsq
1-pchisq(H,1) 
