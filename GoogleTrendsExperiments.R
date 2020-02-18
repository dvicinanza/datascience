#install.packages("gtrendsR")
#install.packages("reshape2")
#install.packages("RApiDatetime")

library(gtrendsR)
library(reshape2)
library(zoo)

#Sys.setenv(TZ="Europe/London")

keyword=c("eduroam")
#geo=c("NL","GB","US")
geo=c("US")

#restricted to geo variable indicated above
res <- gtrends(keyword = keyword, geo = geo)

#worldwide search
#res <- gtrends(keyword = keyword)

plot(res)

rel_queries <- res$related_queries$value
#rel_queries <- res$related_topics$value
head(rel_queries)

x=c()
for (i in 1:length(rel_queries)){
  x=cbind(x,gtrends(rel_queries[i], geo = geo)$interest_over_time$hits)
}
y=zoo(res$interest_over_time$hits,as.Date(res$interest_over_time$date,"%d/%m/%Y"))
x=as.data.frame(x)
names(x)=c(res$related_queries$value[1:length(rel_queries)])

#View(x)

##To plot all the correlators
#require(ggplot2)
#require(reshape2)
#df <- data.frame(time =as.Date(res$interest_over_time$date,"%d/%m/%Y"),
#                 x=x)
#df <- melt(df ,  id.vars = 'time', variable.name = 'series')

## plot on same grid, each series colored differently -- 
## good if the series have same scale
#ggplot(df, aes(time,value)) + geom_line(aes(colour = series))

## or plot on different plots
#ggplot(df, aes(time,value)) + geom_line() + facet_grid(series ~ .)


#removing columns with duplicate name
temp=x; temp <- temp[, !duplicated(colnames(temp))]
x=temp


#Bayesian analysis using BSTS

library(bsts)

# set a few parameters
numiter <- 2000
npred <- 5
# describe state space model consisting of
# trend and seasonal components
ss <- AddLocalLinearTrend(list(),y)

#Change nseason to 12 if analysing monthly data, leave 52 for weekly
ss <- AddSeasonal(ss,y,nseasons=52)

# estimate the model using just autoregression
model <- bsts(y, ss, niter = 2000)
plot(model)
#plot(model,"comp")

# estimate the model using related series as posterior
fit2 <- bsts(y~ ., state.specification = ss, data = x, niter = 2000, expected.model.size=npred)
plot(fit2)
# Probability of inclusion of top predictors (p > .15)
plot(fit2,"coef",inc=0.01)
a=plot(fit2,"coef",inc=0.01)
top_correlators=as.vector(names(a$inclusion.prob[a$inclusion.prob>=0.01]))
top_correlators=gsub("`", '', top_correlators)
#plot(fit2,"comp")


#source("~/Documents/GEANT/DataScience_eduroam/oosf.R")

# http://people.ischool.berkeley.edu/~hal/Papers/2015/oosf.R
#######################################
# computes 1-step ahead forecast (H. R. Varian, google)
#######################################
library(dyn)

# y = dependent variable (zoo format)
# x = predictor variables (zoo format)
# k = start predicting at observation k

# uses lag 1
OutOfSampleForecast01 <- function(y,x,k) {
  if (is.zoo(x)==F | is.zoo(y)==F)  stop("x and y must be zoo objects")
  n <- length(y)
  d <- cbind(y,lag(y,-1),x)
  colnames(d) <- c("y","lagy",colnames(x))
  y.actual <- y.pred0 <- y.pred1 <- rep(NA,n)
  for (t in k:(n-1)) {
    reg0 <- lm(y~lagy,data=d[1:t,])
    reg1 <- lm(y~., data=d[1:t,])
    t1 <- t+1
    y.actual[t1] <-  d[t1:t1,]$y
    y.pred0[t1] <- predict(reg0,newdata=d[t1:t1,])
    y.pred1[t1] <- predict(reg1,newdata=d[t1:t1,])
  }
  z <- cbind(y.actual,y.pred0,y.pred1)[(k+1):n,]
  z <- zoo(z,index(y)[(k+1):n])
  return(z)}

# uses lag 1 and lag 12
OutOfSampleForecast12 <- function(y,x,k) {
  if (is.zoo(x)==F | is.zoo(y)==F) stop("x and y must be zoo objects")
  n <- length(y)
  d <- cbind(y,lag(y,-1),lag(y,-12),x)
  colnames(d) <- c("y","lagy.1","lagy.12",colnames(x))
  y.actual <- y.pred0 <- y.pred1 <- rep(NA,n)
  for (t in k:(n-1)) {
    reg0 <- lm(y~lagy.1+lagy.12,data=d[1:t,])
    reg1 <- lm(y~., data=d[1:t,])
    t1 <- t+1
    y.actual[t1] <-  d[t1:t1,]$y
    y.pred0[t1] <- predict(reg0,newdata=d[t1:t1,])
    y.pred1[t1] <- predict(reg1,newdata=d[t1:t1,])
  }
  z <- cbind(y.actual,y.pred0,y.pred1)[(k+1):n,]
  z <- zoo(z,index(y)[(k+1):n])
  return(z)}

# maeDates: computes mean absolute forecast error
# x= output of outOfSampleForecast01 or outOfSampleForecast12
# t0= index of start date
# t1= index of end date
# logged = TRUE if data is in logs, FALSE otherwise

MaeReport <- function(x,t0,t1,logged=TRUE) {
  if (is.zoo(x)==F)  stop("input must be a zoo object")
  if (missing(t0)) t0 <- start(x)
  if (missing(t1)) t1 <- end(x)
  dts <- window(x,start=t0,end=t1)
  plot(dts,plot.type="single",col=1:3)
  if (logged) {
    mae0 <- mean(abs(dts[,1]-dts[,2]))
    mae1 <- mean(abs(dts[,1]-dts[,3]))
  } else {
    mae0 <- mean(abs(dts[,1]-dts[,2])/dts[,1])
    mae1 <- mean(abs(dts[,1]-dts[,3])/dts[,1])
  }
  rpt <-   c(mae0,mae1,1-mae1/mae0)
  names(rpt) <- c("mae.base","mae.trends","mae.delta")
  return(rpt)
}





# choose top 10 predictors
top_correlators=rev(top_correlators)

top_correlators[1:10]

#x1 <- zoo(x[,cbind("connect to eduroam","eduroam wifi","unlv eduroam","eduroam ncsu")],as.Date(res$interest_over_time$date,"%d/%m/%Y"))
#x1 <- zoo(x[,cbind("connecting to eduroam","eduroam wifi","connect to eduroam","eduroam login")],as.Date(res$interest_over_time$date,"%d/%m/%Y"))
x1 <- zoo(x[,top_correlators[1:10]],as.Date(res$interest_over_time$date,"%d/%m/%Y"))

reg1 <- OutOfSampleForecast12(y,x1,k=52*2)
# mae.delta is the ratio of the trends MAE to the base MAE
MaeReport(reg1)




#Seasonal decomposition using Kalman Filtering
legend("topleft", 
       legend = c("actual", "base regression", "regression+correlators"), 
       col = c("black","red","green"))

decompose_eduroam = stats::decompose(ts(y,frequency = 52), "multiplicative")
plot(as.ts(decompose_eduroam$seasonal))
plot(as.ts(decompose_eduroam$trend))
plot(as.ts(decompose_eduroam$random))
plot(decompose_eduroam)
