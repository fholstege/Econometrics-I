

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

# get summary of coefficients and standard errors, create latex table
summary(model_1)
stargazer(model_1)

# define the second model, according to question 4.a.2
model_2 <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)

# get summary of coefficients and standard errors, create latex table
summary(model_2)

### Question 4.b

# define F-test
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

# add these non-linear tranformations to the original model
test_linear <- lm(birthweight ~ age + unmarried + educ + alcohol + smoker + drinks + yest2 + yest3, data = df)
summary(test_linear)

# get residuals for the F-test
RSS_REST <- t(resid(test_linear))%*%resid(test_linear)
RSS_Orig <- t(resid(model_2)) %*% resid(model_2)


# calc p-val for the F-test - likely to have non-linear dynamics
Ftest_linear <- Ftest(k=6, g = 2, n = nrow(df), RSS_REST, RSS_Orig)
Ftest_linear

anova(model_1, model_2)


### Question 4.D

# Test for non-linear dynamics now the dependent variable is logged
model_2_Log <- lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks ,data = df)
summary(model_2_Log)

yest_log <- model_2_Log$fitted.values
yest2_log <- as.matrix(yest_log^2)
yest3_log <- as.matrix(yest_log^3)

# create model with fitted values, squared and to the power three, logged
test_linear_logged <-  lm(log(birthweight) ~ age + unmarried + educ + alcohol + smoker + drinks + yest2_log + yest3_log, data = df)
summary(test_linear_logged)

# get residuals and perform F-test
RSS_REST_Log <-  t(resid(test_linear_logged))%*%resid(test_linear_logged)
RSS_Orig_Log <- t(resid(model_2_Log)) %*% resid(model_2_Log)

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
  

    return(-logL)
}

# log-likelihood funtion to additionally estimate sigma2
LL2 <- function(par, mY, mX){
  n <- nrow(mX)
  k <- ncol(mX)
  mB <- par[1:k]
  sigma2 <- par[k+1]
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

# calculate associated standard errors from fisher information
observed_FI <- solve(ML_1$hessian)
se_1 <- sqrt(diag(observed_FI))

# calculate t-test statistics
t_1 <- ML_1$par/se_1

# manually calculate sigma2
sigma2_ML1<- t(mY-mX_model1%*%ML_1$par)%*%(mY-mX_model1%*%ML_1$par)/nrow(mX_model1)
sigma2_ML1


# additinally estimate sigma2
ML_2 <- optim(c(1,1,1,1,1),LL2,method="BFGS",hessian=T,mY=mY,mX=mX_model1)
ML_2$par


# calculate associated standard errors from fisher information
observed_FI <- solve(ML_2$hessian)
se_2 <- sqrt(diag(observed_FI))

# calculate t-test statistics
t_1 <- ML_2$par/se_2




# construct fake lm object to input into stargazer
dep_var <- "birthweight"
indep_vars <- c("Intercept", "age", "unmarried", "educ", "sigma2")
f <- as.formula(paste(dep_var, " ~ 0 + ", paste(indep_vars, collapse = " + ")))

# construct stargazer table based on calculated statistics
fake_data <- as.data.frame(cbind(1, 1, df))
names(fake_data)[1] <- "Intercept"
names(fake_data)[2] <- "sigma2"
stargazer(lm(f, data = fake_data),
          type = "latex",
          coef = list(ML_2$par), # supply self-estimated coefficients
          se = list(se_2), # supply self-estimated standard errors
          t.auto = TRUE, # calculate t values based on supplied coef and se
          p.auto = TRUE, # calculate p value based on supplied coef and se
          omit.stat = "all",
          align = TRUE,
          digits = 3,
          title = "Results of the Maximum Likelihood Estimation")



# Question 6: GMM estimate
est.Beta.GMM <- inv(t(mX_model1) %*%mX_model1) %*% t(mX_model1) %*% mY
est.Beta.GMM


# specify moment conditions

mXY <- cbind(mX_model1,mY)

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
library(gmm)
GMM <- gmm(g, mXY, c(1,1,1,1),type="iterative", wmatrix = "optimal")

# gmm coefs
GMM$coefficients

# manually calculate sigma2
sigma2_GMM<- t(mY-mX_model1%*%GMM$coefficients)%*%(mY-mX_model1%*%GMM$coefficients)/nrow(mX_model1)
sigma2_GMM

