

# Put here the packages required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(stargazer,
               lmtest,
               matlib)

# load data, get summary statistics
df <- read.csv("Data/DataAS1.csv")
summary(df)

### Question 4.a

# define first model, according to question 4.a.1
model_1 <- lm(birthweight ~ age + unmarried + educ ,data = df)
summary(model_1)

# define the second model, according to question 4.a.2
model_2 <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)
summary(model_2)

### Question 4.b

# looking at the summary of model_2, the alcohol and drinks variables are statistically insignificant (e.g. cannot reject H:0)
# but the smoker variable is statistically significant.

# get table of the coefficients and their significance, showing t-values
stargazer(model_2, report=('vc*t'))

#Answer: can reject the null hypothesis that they are all equal and zero

### Question 4.C

# test for non-linearity in model_2 - ^2 and ^3
yest <- model_2$fitted.values
yest2 <- as.matrix(yest^2)
yest3 <- as.matrix(yest^3)

# add these fitted values squared, to the power three, to the original model
test_linear <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks + yest2 + yest3, data = df)
summary(test_linear)

# get residuals for the F-test
RSS_REST <- t(resid(test_linear))%*%resid(test_linear)
RSS_Orig <- t(resid(model_2)) %*% resid(model_2)

# define F-test
Ftest <- function(k, g, n, RSS_test, RSS_compare){
  
  Ftest = ((RSS_compare - RSS_test)/g)/(RSS_test/(n-k))

  Pval <- 1-pf(Ftest, g, n-k)
  return(Pval)
}

# calc p-val for the F-test - likely to have non-linear dynamics
pval_linear <- Ftest(k=6, g = 2, n = nrow(df), RSS_REST, RSS_Orig)
pval_linear

#Answer: p-value is 9%, e.g. a 9% chance that there are non-linear dynamics that are not captured. Cannot reject H0 at 5% significance

### Question 4.D + E

# test for non-linear dynamics now the dependent variable is logged
model_2_Log <- lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)
summary(model_2_Log)

# create model with fitted values, squared and to the power three, logged
test_linear_logged <-  lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks + yest2 + yest3, data = df)
summary(test_linear_logged)

# get residuals and perform F-test
RSS_REST_Log <-  t(resid(test_linear_logged))%*%resid(test_linear_logged)
RSS_Orig_Log <- t(resid(model_2_Log)) %*% resid(model_2_Log)

pval_linear_logged <- Ftest(k=6, g=2, n=nrow(df), RSS_REST_Log, RSS_Orig_Log)
pval_linear_logged

#Answer: 
# - The model gets an adjusted R-squared of 4.8%. This is worse than the model with the dependent variable unlogged (5.5%)
# - The same variables remain statistically significant
# - the F-test has a higher p-value, showing that logging the dependent variable makes it more likely to be linear
# - Question: how to interpret AIC? weird that its lower for logged model, despite having lower adjusted R-Squared

### compare linear y and log y model
summary(model_2)
summary(model_2_Log)
AIC(model_2)
AIC(model_2_Log)


### Question 5: compute ML 

# get dependent and independent variables
mY <- as.matrix(df$birthweight)
mX_model1 <- as.matrix(data.frame(1, df$age, df$unmarried, df$educ))

# define log-likelihood function
LL <- function(mB, mY, mX, sigma2){
  n      <- nrow(mX)
  e <- mY - (mX %*% mB)
  
  logL <- .5*n*log(2*pi)-.5*n*log(sigma2)-((t(e)%*%e)/(2*sigma2))
  return(-logL)

}

# set initial beta's
Beta <- c(1,1,1,1)
mBeta <- as.matrix(Beta)

# optimize log likelihood function
ML_1 <- optim(mBeta,LL,method="BFGS",hessian=T,mY=mY,mX=mX_model1, sigma2 = 1)

# has the same beta's as the lm() 
ML_1$par

# Question 6: GMM estimate
est.Beta.GMM <- inv(t(mX_model1) %*%mX_model1) %*% t(mX_model1) %*% mY
est.Beta.GMM

