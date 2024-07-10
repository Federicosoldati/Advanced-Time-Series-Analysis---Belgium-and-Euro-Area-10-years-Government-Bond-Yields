# R code for Advance Time Series Analysis assignment 
# Student: Federico Soldati - r0924528


library(CADFtest)
library(forecast)
library(vars)

##### Univariate TSA: Explorative analysis #####


BE <- read.table("Bonds/BE.csv",header = TRUE,sep=",")
BE_ts <- ts(BE[,2],frequency = 4,start=c(1970,1),end = c(1995,4))

ts.plot(BE_ts,col="red")

acf(BE_ts)

max.lag <- round(sqrt(length(BE_ts)))
Box.test(BE_ts, lag = max.lag, type = "Ljung-Box")

CADFtest(BE_ts,type="drift",criterion="BIC",max.lag.y=max.lag)

# The series is not stationary

dBE_ts <- diff(BE_ts) 

ts.plot(dBE_ts,col="red")

CADFtest(dBE_ts,type="drift",criterion="BIC",max.lag.y=max.lag)

hist(dBE_ts, breaks = 30)
boxplot(dBE_ts)

Box.test(dBE_ts, lag = max.lag, type = "Ljung-Box")

# The series is now stationary


# We'll analyse dBE_ts during the univariate analysis

monthplot(BE_ts)

# It appears there is no seasonal effect

TREND <- 1:104
Q1 <- rep(c(1,0,0,0),26)
Q2 <- rep(c(0,1,0,0),26)
Q3 <- rep(c(0,0,1,0),26)
Q4 <- rep(c(0,0,0,1),26)

fit <- lm(BE_ts ~ TREND + Q1 + Q2 + Q3 + Q4)
summary(fit)

# Where Q4 is omitted to avoid perfect multi-collinearity. 

ts.plot(fit$residuals)

CADFtest(fit$residuals,max.lag.y = max.lag,type="drift",criterion = "BIC")

# The model is crearly not valid. It is not surpring since the series has no trend or seasonality, 
# so dummys quarters have no effect. The residuals are also not white noise.

##### Univariate TSA: Fit Arima #####

par(mfrow=c(2,1))

acf(diff(BE_ts))
pacf(diff(BE_ts))

# It's possible to see one significant correlation in the acf and one in pacf. 
# It has been choosed to fit and confront two different models, a MA(1) and a AR(1).

fit_ar <- arima(BE_ts, order = c(1,1,0), seasonal = list(order = c(0,0,0)))
summary(fit_ar)
par(mfrow=c(1,1))
plot(fit_ar$residuals)
par(mfrow=c(2,1))
acf(fit_ar$residuals)
pacf(fit_ar$residuals)

abs(fit_ar$coef/sqrt(diag(fit_ar$var.coef)))
Box.test(fit_ar$residuals, lag = max.lag, type = "Ljung-Box")

# The model is valid

fit_ma <- arima(BE_ts, order = c(0,1,1), seasonal = list(order = c(0,0,0)))
summary(fit_ma)
par(mfrow=c(1,1))
plot(fit_ma$residuals)
par(mfrow=c(2,1))
acf(fit_ma$residuals)
pacf(fit_ma$residuals)

abs(fit_ma$coef/sqrt(diag(fit_ma$var.coef)))
Box.test(fit_ma$residuals, lag = max.lag, type = "Ljung-Box")

# Both the models are valid, as they both show only one correlation significant 
# (which is ok since the significance level of the CI)

ts.plot(fit_ar$residuals)

par(mfrow=c(1,1))
qqnorm(fit_ar$residuals)
qqline(fit_ar$residuals)

# Residuals are homoscedastic, is not necessary to fit a garch model

# Comparing the two models,using AIC and SIC

AIC(fit_ar);AIC(fit_ma)
AIC(fit_ar,k=log(length(BE_ts)));AIC(fit_ma,k=log(length(BE_ts)))

# Both the measures indicates the AR(1) model as better 

##### Univariate TSA: Forecasting ##### 

myforecastAR<-predict(fit_ar,n.ahead=6)
expected<-myforecastAR$pred

# The confidence bounds of the 95% prediction interval:
lower<-myforecastAR$pred-qnorm(0.975)*myforecastAR$se
upper<-myforecastAR$pred+qnorm(0.975)*myforecastAR$se

par(mfrow=c(1,1))
plot.ts(BE_ts,xlim=c(1994,1997.2),ylim=c(2.5,12))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")


myforecastMA<-predict(fit_ma,n.ahead=6)
expected<-myforecastMA$pred


# The  of the 95% prediction interval:
lower<-myforecastMA$pred-qnorm(0.975)*myforecastMA$se
upper<-myforecastMA$pred+qnorm(0.975)*myforecastMA$se


plot.ts(BE_ts,xlim=c(1994,1997.2),ylim=c(2.5,12))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")


##### Univariate TSA: Confronting the two ARIMA models  ##### 

# Now it is used an expanding-window approach to forecast values of BE_ts for 1-period.
# The loops are necessary to compare the forecasts for the ARIMA(1,1,0) model and the ARIMA(0,1,1) model.

y<-BE_ts
S=round(0.75*length(y))
h=1
error1.h<-c()
for (i in S:(length(y)-h)){
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h)){
  mymodel.sub<-arima(y[1:i], order = c(0,1,1),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}

# To evaluate the forecasting performance of the two models the Mean Absolute Error (MAE) is calculated:

MAE1 <- mean(abs(error1.h))
MAE2 <- mean(abs(error2.h))

MAE1;MAE2

# The AR(1) model presents a better (lower) Mean Absolute Error

# Now a Diebold-Mariano test is performed, to see if the forecast performance of the two models is significantly different

dm.test(error1.h,error2.h,h=h,power=1)

# P-value = 0.953 > 5%, thus H0 is not reject and it can be concluded that the forecast
# performance of the two models, using the absolute value loss, is not significantly different.

# The same is now done using the squared value loss.

error1.h<-c()
for (i in S:(length(y)-h)){
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h<-c()
for (i in S:(length(y)-h)){
  mymodel.sub<-arima(y[1:i], order = c(0,1,1),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
MSE1 <- mean(abs(error1.h^2))
MSE2 <- mean(abs(error2.h^2))

MSE1;MSE2

# The MA(1) obtains the lowest MSE.


dm.test(error1.h,error2.h,h=h,power=2)


# Diebold-Mariano test:
# P-value = 0.793 > 5%, thus H0 is not reject and it can be concluded  that the forecast performance
# of the two models, using the squared value loss, is not significantly different.

##### Multivariate TSA: Cointegration test #####

EU <- read.table("Bonds/Eu.csv",header = TRUE,sep=",")
EU_ts <- ts(EU[,2],frequency = 4,start=c(1970,1),end = c(1995,4))

ts.plot(BE_ts,EU_ts,col=(c("red","blue")))

ts.plot(BE_ts-EU_ts)

max.lag <- round(sqrt(length(EU_ts)))
CADFtest(EU_ts,type="drift",criterion="BIC",max.lag.y=max.lag)

# EU_ts is not stationary

dEU_ts <- diff(EU_ts)

ts.plot(dEU_ts,col="blue")

CADFtest(dEU_ts,type="drift",criterion="BIC",max.lag.y=max.lag)

# The series in differences is stationary
# BE_ts and UE_ts are I(1)

# Cointegration test

fit_ci <- lm(BE_ts ~ EU_ts)

# Saving the residuals:

res_fit_ci <- fit_ci$residuals

# Unit root test on the residuals to check their stationarity 

CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max.lag)

# Test statistics: -3.6296, which is smaller that the Engle-Granger ADF test statistics for one explanatory variable −3.41. 
# H0 of no cointegration is reject and it is possible to conclude that BE_ts and EU_ts are cointegrated.


##### Multivariate TSA: Error-Correction Model #####

# Construction and estimation of an error-correction model (ECM):
ECT <- res_fit_ci[-length(res_fit_ci)]
fit_ecm <- lm(dBE_ts ~ dEU_ts + ECT)

# Checking the validity of the model:
summary(fit_ecm)
max.lag <- round(sqrt(length(fit_ecm$residuals)))

Box.test(fit_ecm$residuals, lag = max.lag, type="Ljung-Box")

# The Q-test on the residual has p-value = 0.07 > 5% H0 is not rejected: the model is valid. 
# However, the decision to reject H0 in the test could be discussed given the low p-value.

# R2 = 53.05%, thus 53% of the variance of dBE_ts is explained by the regressors in the ECM model
# The overall F-statistics has p-value = 0.00 < 5%, thus H0 is rejected and the regressors are jointly significant.-
# Also, all the variables are significant (apart the intercept)

# The absolute value of the estimated coefficient of ECT is the speed of adjustment towards
# the long run equilibrium   

##### Multivariate TSA: ADLM(2) #####


# Estimation of an Autoregressive Dynamic model of order 2 for dBE_ts (ADLM(2))
n <- length(dBE_ts)
lag <- 2
dBE.0 <- dBE_ts[(lag+1):n]
dEU.0 <- dEU_ts[(lag+1):n]
dBE.1 <- dBE_ts[lag:(n-1)]
dEU.1 <- dEU_ts[lag:(n-1)]
dBE.2 <- dBE_ts[(lag-1):(n-2)]
dEU.2 <- dEU_ts[(lag-1):(n-2)]

fit_adlm <- lm(dBE.0 ~ dBE.1+dBE.2+dEU.1+dEU.2)

ts.plot(fit_adlm$residuals)
acf(fit_adlm$residuals)
Box.test(fit_adlm$residuals, lag = max.lag, type="Ljung-Box")

# We note 1 significant autocorrelation lag, but the model is valid

summary(fit_adlm)

# R-squared:  0.1594, 15.94% of the variance of dBE.0 is explained by the model 

# The overall F-statistics has p-value = 0.002 < 5%, it is possible to conclude that the regressors are jointly significant  
# The only significant variable is dEU.1, that has a positive value.
# We can infer that the value of dBE.0 depends more on the previous value of dEU than on the previous value of dBE

# Testing for Granger causality comparing the ADLM(2) with the model without lagged explanatory variables:

fit_adlm_nox <- lm(dBE.0 ~ dBE.1+dBE.2)
anova(fit_adlm,fit_adlm_nox)

# P-value = 0.04577 < 5%, thus we do reject H0 of no Granger Causality. We conclude
# that dEU has a incremental explanatory power in predicting dBE

##### Multivariate TSA: VAR #####

data<-data.frame(dBE_ts,dEU_ts)
names(data)<-c("dBE","dEU")
VARselect(data,lag.max=10,type="const")

# The maximum number of lags is specified to be 10, a constant is included in each equation.
# The order of the VAR model selected by Schwarz information criterion is 1.

# Estimating the VAR model

fit_varautom<-VAR(data,type="const",p=1)
summary(fit_varautom)

# The plots the correlograms and cross-correlogram of the residual are creted to check the model. 

varautom_residuals<-resid(fit_varautom)

acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,2])
Box.test(varautom_residuals[,1], lag = max.lag, type="Ljung-Box")
Box.test(varautom_residuals[,2], lag = max.lag, type="Ljung-Box")

# The two model is validated by the Box-Ljung test and by the correlograms and the cross-correlogram


# Impulse response functions:

irf_var<-irf(fit_varautom,ortho=FALSE,boot=TRUE)
plot(irf_var)

# The IRFs provide an easy way to interpret the estimated coefficients of the VAR model. 
# Given aunitary impulse in dBE at time t, it is not observed any significant response. 
# Given a unitary impulse in dEU at time t, we observe a significant positive 
# response of dBE at time t + 1.

# forecasting using the var model

myforecastVAR<-predict(fit_varautom,n.ahead=6)
expected_var<- myforecastVAR$fcst$dBE
expected[1:6] <- expected_var[,1]
# The confidence bounds of the 95% prediction interval:
lower_var<-expected_var[,2]
upper_var<-expected_var[,3]
lower[1:6] <- lower_var
upper[1:6] <- upper_var

par(mfrow=c(1,1))
plot.ts(BE_ts,xlim=c(1994,1997.2),ylim=c(2.5,12))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")

