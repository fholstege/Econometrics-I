# Authors: Floris Holstege, Markus Mueller
# Goal: Code for econometrics I, assignment I
# Date: 17/11/2020

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(stargazer,
               lmtest,
               matlib,
               gmm)

# load data, get summary statistics
df <- read.csv("Data/DataAS1.csv")
summary(df)

### Question 4.a

# define first model, according to question 4.a.1
model_1 <- lm(birthweight ~ age + unmarried + educ ,data = df)

# get summary of coefficients and standard errors, create latex table
summary(model_1)
stargazer(model_1)

# define the second model, according to question 4.a.2
model_2 <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)

# get summary of coefficients and standard errors, create latex table
summary(model_2)
stargazer(model_2)


### Question 4.b

# define F-test. Returns p-value of said test
Ftest <- function(k, g, n, RSS_test, RSS_compare){
  Ftest = ((RSS_compare - RSS_test)/g)/(RSS_test/(n-k))
  Pval <- 1-pf(Ftest, g, n-k)
  return(Pval)
}

# residual sum of squares for both model 1 and 2
RSS_1 <- t(resid(model_1)) %*% resid(model_1)
RSS_2 <-  t(resid(model_2)) %*% resid(model_2)

# get result of F-test for comparing model 1 and 2
Ftest_model1_2 <- Ftest(3,3,nrow(df), RSS_2, RSS_1)
Ftest_model1_2

### Question 4.C

# Create non-linear transformations of the Yhat in model 2
yest <- model_2$fitted.values
yest2 <- as.matrix(yest^2)
yest3 <- as.matrix(yest^3)

# add these non-linear tranformations of Yhat to the original model
test_linear <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks + yest2 + yest3, data = df)
summary(test_linear)

# get residuals for the F-test
RSS_REST <- t(resid(test_linear))%*%resid(test_linear)
RSS_Orig <- t(resid(model_2)) %*% resid(model_2)

# calc p-val for the F-test - likely to have non-linear dynamics
Ftest_linear <- Ftest(k=6, g = 2, n = nrow(df), RSS_REST, RSS_Orig)
Ftest_linear

### Question 4.D

# Test for non-linear dynamics now the dependent variable is logged
model_2_Log <- lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)
summary(model_2_Log)

# get Yhat from logged model, and create non-linear transformations of this
yest_log <- model_2_Log$fitted.values
yest2_log <- as.matrix(yest_log^2)
yest3_log <- as.matrix(yest_log^3)

# create model with with non-linear transformations of Yhat
test_linear_logged <-  lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks + yest2_log + yest3_log, data = df)
summary(test_linear_logged)

# get residuals
RSS_REST_Log <-  t(resid(test_linear_logged))%*%resid(test_linear_logged)
RSS_Orig_Log <- t(resid(model_2_Log)) %*% resid(model_2_Log)

# perform F-test for logged model
Ftest_linear_logged <- Ftest(k=6, g=2, n=nrow(df), RSS_REST_Log, RSS_Orig_Log)
Ftest_linear_logged_2 <- ((RSS_Orig_Log-RSS_REST_Log)/2)/(RSS_REST_Log/(3000-6))
Ftest_linear_logged

### Question 4.E

### compare linear y and log y model in terms of R2, Adjusted R2
summary(model_2)
summary(model_2_Log)

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

# has approximately the same beta's as the lm() 
ML_1$par

# manually calculate sigma2
sigma2_ML1<- t(mY-mX_model1%*%ML_1$par)%*%(mY-mX_model1%*%ML_1$par)/nrow(mX_model1) 
sigma2_ML1


# Question 6: GMM estimate

# add X and Y together in matrix
mXY <- cbind(mX_model1,mY)

# define sample moments 
g <- function(a, mXY){
  # each x is orthogonal to residual
  m1 <- mXY[,1]*(mXY[,5]-a[1]*mXY[,1]-a[2]*mXY[,2]-a[3]*mXY[,3]-a[4]*mXY[,4])
  m2 <- mXY[,2]*(mXY[,5]-a[1]*mXY[,1]-a[2]*mXY[,2]-a[3]*mXY[,3]-a[4]*mXY[,4])
  m3 <- mXY[,3]*(mXY[,5]-a[1]*mXY[,1]-a[2]*mXY[,2]-a[3]*mXY[,3]-a[4]*mXY[,4])
  m4 <- mXY[,4]*(mXY[,5]-a[1]*mXY[,1]-a[2]*mXY[,2]-a[3]*mXY[,3]-a[4]*mXY[,4])
  
  f  <- cbind(m1,m2,m3,m4)
  return(f)
}

# estimate via gmm package
GMM <- gmm(g, mXY, c(0,0,0,0),type="iterative", wmatrix = "optimal")

# get gmm coeficients
GMM$coefficients

# manually calculate sigma2
sigma2_GMM<- t(mY-mX_model1%*%GMM$coefficients)%*%(mY-mX_model1%*%GMM$coefficients)/nrow(mX_model1)
sigma2_GMM

