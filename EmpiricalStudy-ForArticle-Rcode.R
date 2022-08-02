#Code for the empirical study
#Last update: 02.08.2022
library(R.matlab)
library(lubridate)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(reshape2)
library(trawl)
library(ambit)

#Read in the data
data <- readMat('A_20200504_data.mat')
str(data)

#Multiply data by 100 to convert to dollar cents
spread <-data$spr*100
plot(spread,type="l")
summary(spread)

#Since the minimum spread level in the data is one tick (one dollar cent), 
#we work on this time series minus one, i.e on $x_t = s_t-1$.
allx <- spread-1

summary(allx)

x <- data$y #spread data in cent (5s sampling)
n<-length(x)


plot(x, type="l")
summary(x)
acf(x)

#Plots using ggplot
theme_update(text = element_text(size=30)) #Increase font size in labels
#Plot as time series with correct times
time_seq <- seq(from = lubridate::mdy_hms("04-05-2020 10:30:00"),to = lubridate::mdy_hms("04-05-2020 16:00:00"), by = "5 sec")
spread.ts <- data.frame(time=time_seq,
                        value=x)
p <- ggplot(spread.ts, aes(x=time,y=value))+
  geom_line()+
  xlab("Time")+
  ylab("Cents")
p
ggsave("TS.eps", width = 20, height = 20, units = "cm")

#Plot the acf
my_acf <- acf(x, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab("Autocorrelation")
q
ggsave("ACF.eps", width = 20, height = 20, units = "cm")


# Fitting the trawl function
my_lag <-n-1
my_Delta <-1/12#sampling in 5s=1/12 min
esttrawlfct <- nonpar_trawlest(x, Delta=my_Delta, lag=my_lag)$a_hat
plot(esttrawlfct)

plot(esttrawlfct[1:100])
abline(h=0,col="red")

#print(esttrawlfct)

#Using ggplot
l_seq <- seq(from = 0,to = (my_lag-1), by = 1)
esttrawlfct.data <- data.frame(l=l_seq[1:31],
                               value=esttrawlfct[1:31])
p2 <- ggplot(esttrawlfct.data, aes(x=l,y=value))+
  geom_point(size=3)+
  xlab("l")+
  ylab(TeX("$\\hat{a}(\\cdot)$"))
p2
ggsave("hata.eps", width = 20, height = 20, units = "cm")


# Confidence intervals for the trawl function
#We plot the asymptotic variance on our $\Delta$-grid:
c4 <- c4est(x, my_Delta)  
av_vector <-numeric(my_lag)
for(i in 1:my_lag){
  av_vector[i] <- asymptotic_variance_est((i-1)*my_Delta, c4, varlevyseed=1, my_Delta, esttrawlfct)$v
}
#Note that av_vector[1] contains the value for i=0 etc
plot(av_vector)
plot(sqrt(av_vector/(my_Delta*n)))
plot(av_vector[1:31])

#Using ggplot
l_seq <- seq(from = 0,to = (my_lag-1), by = 1)
av.data <- data.frame(l=l_seq[1:31],
                      value=av_vector[1:31])
p4 <- ggplot(av.data, aes(x=l,y=value))+
  geom_point( size=3)+
  xlab("l")+
  ylab(TeX("$\\hat{\\sigma}_a^2(\\cdot)$"))
p4
ggsave("hatav.eps", width = 20, height = 20, units = "cm")

#We can now compute the upper and lower confidence bounds.
#$$
#  \hat a(t)\pm z_{\beta/2}\sqrt{\frac{\sigma^2_a(t)}{n\Delta_n}}
#$$
#
upper <-numeric(my_lag+1)
lower <-numeric(my_lag+1)
z<- qnorm(0.05/2,lower.tail = FALSE)
for(i in 0:(my_lag-1)){
  upper[i+1] <- esttrawlfct[i+1]+z*sqrt(av_vector[i+1]/(my_Delta*n))
  lower[i+1] <- esttrawlfct[i+1]-z*sqrt(av_vector[i+1]/(my_Delta*n))
}

tt<-(1:(my_lag+1))
plot(esttrawlfct)
lines(tt, upper, lty=3, col="red")
lines(tt, lower, lty=3, col="blue")


plot(esttrawlfct[1:20])
lines(tt[1:20], upper[1:20], lty=3, col="red")
lines(tt[1:20], lower[1:20], lty=3, col="blue") 

###With ggplot
l_seq <- seq(from = 0,to = (my_lag-1), by = 1)

xx<-l_seq[1:31]
y1<-upper[1:31]
y2<-lower[1:31]
y3<-esttrawlfct[1:31]


df <- data.frame(xx,y1,y2,y3)

mdf <- melt(data=df,id.vars="xx")

g5<-ggplot(mdf, aes( x=xx, y=value, colour=variable, group=variable, shape=variable )) + 
  geom_point(size=3, show.legend=FALSE)+#geom_line()+
  #geom_line(data = data.frame(x,y1,y2))+
  scale_color_manual(values=c("y1"="blue","y2"="blue","y3"="black")) +
  scale_shape_manual(values = c(1, 1, 19))+
  #scale_linetype_manual(values=c("y1"="dashed","y2"="dashed","y3"="solid"))+
  xlab("l")+
  ylab(TeX("$\\hat{a}(\\cdot)$"))+
  theme(legend.position = "none")

g5
ggsave("CIs.eps", width = 20, height = 20, units = "cm")

which(y3<0)
#We note that the trawl function estimate is negative for index 18, 22, 23, 25, 28, 31 etc.
#(i=17, 21, 22, 24 etc.).

##########################
#Forecasting
##########################
is_length <- 3000 #in sample length

my_lag <-is_length-1
my_Delta <-1/12 #sampling in 5s=1/12 min





hrange <- 20
noos <- n-is_length-hrange #=941 for is_length=3000
#Create forecast matrices:
CondMean <-matrix(0, nrow=noos, ncol=hrange)
CondMean_IV <-matrix(0, nrow=noos, ncol=hrange)
CondMean_linear <-matrix(0, nrow=noos, ncol=hrange)
NaiveFC <-matrix(0, nrow=noos, ncol=hrange)
NaiveMeanFC <-matrix(0, nrow=noos, ncol=hrange)

ActualValues <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection_linear <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference_linear <- matrix(0, nrow=noos, ncol=hrange)


for(tau in 1:noos){
  #tau<-1
  
  
  data <-x[tau:(is_length+tau-1)] # take new subsample of length is_length
  actualvalue <- x[(is_length+tau-1)]
  datamean <- mean(data) #sum(data)/(is_length) is also correct
  
  #Estimate trawl function
  esttrawlfct_fc <- nonpar_trawlest(data, Delta=my_Delta, lag=my_lag)$a_hat #estimated trawl   function for forecasting 
  
  #Differences between points for linear approximation (set sign to positive)
  diff_esttrawlfct_fc <- -diff(esttrawlfct_fc)
  
  #Estimate the Lebesgue measure of the trawl set (components)
  #lebA_fc <-sum(esttrawlfct_fc)*my_Delta
  lebA_fc <- LebA_est(data, my_Delta)
  
  #Estimate the Lebesgue measure of the trawl set (components)
  #using linear approximation between points
  
  lebA_linear_fc <-sum(esttrawlfct_fc[2:my_lag])*my_Delta   +sum(diff_esttrawlfct_fc[2:(my_lag-1)])*my_Delta*0.5
  
  
  for(h in 1:hrange){ #forecasting horizon loop in 1:hrange
    #h<-1
    #Estimate Leb(A intersection A_h)
    
    #lebAintersection_fc <-sum(esttrawlfct_fc[(h+1):my_lag])*my_Delta
    lebAintersection_fc <-LebA_slice_est(data, my_Delta, (h+1)*my_Delta)$LebAintersection 
    
    lebAintersection_linear_fc <-sum(esttrawlfct_fc[(h+2):my_lag])*my_Delta   +sum(diff_esttrawlfct_fc[(h+1):(my_lag-1)])*my_Delta*0.5
    #Estimate Leb(A \ A_h)
    
    lebAsetdiff_fc <-lebA_fc-lebAintersection_fc#sum(esttrawlfct_fc[1:h])*my_Delta
    
    
    lebAsetdiff_linear_fc <-lebA_linear_fc-lebAintersection_linear_fc
    
    #ActualValue
    ActualValues[tau,h]<- x[(is_length+tau-1+h)]
    #Compute forecast(s)
    #naive
    NaiveFC[tau,h] <- actualvalue
    NaiveMeanFC[tau,h] <- datamean
    
    #trawl formula
    
    
    Weight_Intersection[tau, h]<-lebAintersection_fc/lebA_fc
    Weight_SetDifference[tau,h]<-lebAsetdiff_fc/lebA_fc
    
    Weight_Intersection_linear[tau, h]<-lebAintersection_linear_fc/lebA_linear_fc
    Weight_SetDifference_linear[tau,h]<-lebAsetdiff_linear_fc/lebA_linear_fc
    
    CondMean[tau,h] <- actualvalue*Weight_Intersection[tau, h]+datamean*Weight_SetDifference[tau,h]
    
    CondMean_IV[tau,h]<-round(CondMean[tau,h])
    
    CondMean_linear[tau,h] <- actualvalue*Weight_Intersection_linear[tau, h]+datamean*Weight_SetDifference_linear[tau,h]
    
  } #end for h
  
}#end for tau


my_mse <-function(x, y){
  length<-length(x)
  return(sum((x-y)^2)/length)
}

my_mae <-function(x, y){
  length<-length(x)
  return(sum(abs(x-y))/length)
}

#Computing the error measures

CondMean_MAE <- numeric(hrange)
CondMean_MSE <- numeric(hrange)

CondMean_IV_MAE <- numeric(hrange)
CondMean_IV_MSE <- numeric(hrange)

CondMean_linear_MAE <- numeric(hrange)
CondMean_linear_MSE <- numeric(hrange)

Naive_MAE <- numeric(hrange)
Naive_MSE <- numeric(hrange)

NaiveMean_MAE <- numeric(hrange)
NaiveMean_MSE <- numeric(hrange)

for(h in 1:hrange){
  CondMean_MAE[h] <-my_mae(CondMean[,h], ActualValues[,h])
  CondMean_MSE[h] <-my_mse(CondMean[,h], ActualValues[,h])
  
  CondMean_IV_MAE[h] <-my_mae(CondMean_IV[,h], ActualValues[,h])
  CondMean_IV_MSE[h] <-my_mse(CondMean_IV[,h], ActualValues[,h])
  
  CondMean_linear_MAE[h] <-my_mae(CondMean_linear[,h], ActualValues[,h])
  CondMean_linear_MSE[h] <-my_mse(CondMean_linear[,h], ActualValues[,h])
  
  Naive_MAE[h] <-my_mae(NaiveFC[,h], ActualValues[,h])
  Naive_MSE[h] <-my_mse(NaiveFC[,h], ActualValues[,h])
  
  NaiveMean_MAE[h] <-my_mae(NaiveMeanFC[,h], ActualValues[,h])
  NaiveMean_MSE[h] <-my_mse(NaiveMeanFC[,h], ActualValues[,h])
}

print(CondMean_MAE/Naive_MAE)

print(CondMean_MSE/Naive_MSE)



print(CondMean_IV_MAE/Naive_MAE)

print(CondMean_IV_MSE/Naive_MSE)

print(CondMean_linear_MAE/Naive_MAE)

print(CondMean_linear_MSE/Naive_MSE)

print(CondMean_linear_MAE/CondMean_MAE)

print(CondMean_linear_MSE/CondMean_MSE)

#Plot the errors
plot(CondMean_MSE/Naive_MSE, type="h")

plot(CondMean_linear_MSE/Naive_MSE, type="h")

plot(CondMean_linear_MSE/CondMean_MSE, type="h")



plot(x[(is_length+1):n],type="l")
#lines(CondMean[,1])
lines(CondMean[,1],col="red")
lines(NaiveFC[,1],col="blue")
lines(NaiveMeanFC[,1],col="green")
lines(CondMean_linear[,1],col="brown")


plot(x[(is_length+1):((is_length+1)+100)],type="l")
#lines(CondMean[1:100,1])
lines(CondMean[1:100,1],col="red")
lines(NaiveFC[1:100,1],col="blue")
lines(NaiveMeanFC[1:100,1],col="green")
lines(CondMean_linear[1:100,1],col="brown")

plot(x[(is_length+1):((is_length+1)+100)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[1:100,10],col="red")
lines(NaiveFC[1:100,10],col="blue")
lines(NaiveMeanFC[1:100,10],col="green")
lines(CondMean_linear[1:100,10],col="brown")

plot(x[((is_length+1)+20):((is_length+1)+20+100)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[(1+20):(100+20),20],col="red")
lines(NaiveFC[(1+20):(100+20),20],col="blue")
lines(NaiveMeanFC[(1+20):(100+20),20],col="green")
lines(CondMean_linear[(1+20):(100+20),20],col="brown")

plot(x[((is_length+1)+20):((is_length+1)+20+30)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[(1+20):(30+20),20],col="red")
lines(NaiveFC[(1+20):(30+20),20],col="blue")
lines(NaiveMeanFC[(1+20):(30+20),20],col="green")
lines(CondMean_linear[(1+20):(30+20),20],col="brown")

## Comparing the nonparametric forecasts with a parametric forecast
##We fit a NegBin trawl model with supOU trawl function.
#We first fit such a model to the full sample.


data <-x

#Estimate trawl function parameters
esttrawlfct_fc <- fit_LMtrawl(data, Delta=my_Delta, GMMlag=10,plot=TRUE) #estimated trawl   function for forecasting 

my_alpha <- esttrawlfct_fc$alpha
my_alpha
my_H <- esttrawlfct_fc$H
my_H


#Estimate the Lebesgue measure of the trawl set (components)
lebA_fc_c1 <-my_alpha/(my_H-1)
lebA_fc_c1

#Estimate negaitve binomial parameters
negbinpar <- fit_marginalNB(data, lebA_fc_c1, plot=TRUE)

my_theta <- negbinpar$theta
my_m <- negbinpar$m
my_mtilde <- (1-my_theta)^2/my_theta

my_c <- my_theta*my_m/(1-my_theta)^2

lebA_fc <-lebA_fc_c1*my_c

my_theta
my_m
my_mtilde
my_c
lebA_fc


#Next we fit a NegBin trawl model with supOU trawl function
#to the data using the same rolling window as in the nonparametric set-up.


is_length <- 3000 #in sample length

my_lag <-is_length-1
my_Delta <-1/12#sampling in 5s=1/12 min




hrange <- 20
noos <- n-is_length-hrange #=941 for is_length=3000

#noos <- 4
#Create forecast matrices:
CondMean_p <-matrix(0, nrow=noos, ncol=hrange)
CondMean_IV_p <-matrix(0, nrow=noos, ncol=hrange)

ActualValues <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection_p <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference_p <- matrix(0, nrow=noos, ncol=hrange)

parvec_alpha <-numeric(noos)
parvec_H <- numeric(noos)
parvec_theta <- numeric(noos)
parvev_m <- numeric(noos)
parvec_mtilde <-numeric(noos)
parvec_c <- numeric(noos)
parvec_lebA <-numeric(noos)

set.seed(1)

for(tau in 1:noos){
  #tau<-1
  
  
  data <-x[tau:(is_length+tau-1)] # take new subsample of length is_length
  actualvalue <- x[(is_length+tau-1)]
  #datamean <- mean(data) #sum(data)/(is_length) is also correct
  
  
  #Estimate trawl function parameters
  esttrawlfct_fc <- fit_LMtrawl(data, Delta=my_Delta, GMMlag=10) #estimated trawl   function for forecasting 
  
  my_alpha <- esttrawlfct_fc$alpha
  my_H <- esttrawlfct_fc$H
  
  #Estimate the Lebesgue measure of the trawl set (components)
  lebA_fc_c1 <-my_alpha/(my_H-1)
  
  #Estimate negaitve binomial parameters
  negbinpar <- fit_marginalNB(data, lebA_fc_c1)
  
  my_theta <- negbinpar$theta
  my_m <- negbinpar$m
  my_mtilde <- (1-my_theta)^2/my_theta
  
  my_c <- my_theta*my_m/(1-my_theta)^2
  
  lebA_fc <-lebA_fc_c1*my_c
  
  parvec_alpha[tau] <-my_alpha
  parvec_H[tau] <- my_H
  parvec_theta[tau] <- my_theta
  parvev_m[tau] <- my_m
  parvec_mtilde[tau] <-my_mtilde
  parvec_c[tau] <- my_c
  parvec_lebA[tau] <-lebA_fc
  
  for(h in 1:hrange){ #forecasting horizon loop in 1:hrange
    #h<-1
    #Estimate Leb(A intersection A_h)
    
    lebAintersection_fc <-my_c*my_alpha/(my_H-1)*(1+h*my_Delta/my_alpha)^(1-my_H)
    
    
    #Estimate Leb(A \ A_h)
    
    lebAsetdiff_fc <-lebA_fc-lebAintersection_fc
    
    
    
    #ActualValue
    ActualValues[tau,h]<- x[(is_length+tau-1+h)]
    #Compute forecast(s)
    
    #trawl formula
    
    
    Weight_Intersection_p[tau, h]<-lebAintersection_fc/lebA_fc
    Weight_SetDifference_p[tau,h]<-lebAsetdiff_fc/lebA_fc
    
    
    CondMean_p[tau,h] <- actualvalue*Weight_Intersection_p[tau, h]+(1-my_theta)*lebA_fc*Weight_SetDifference_p[tau,h]
    
    CondMean_IV_p[tau,h]<-round(CondMean_p[tau,h])
    
    
    
  } #end for h
  
}#end for tau


#We briefly consider the parameter estimates for the forecasting exercise:
  

boxplot(parvec_alpha)
boxplot(parvec_H)

boxplot(parvec_theta)
boxplot(parvev_m)

boxplot(parvec_mtilde)
boxplot(parvec_c)

boxplot(parvec_lebA)

#Next we compute the error measures:

CondMean_p_MAE <- numeric(hrange)
CondMean_p_MSE <- numeric(hrange)

CondMean_IV_p_MAE <- numeric(hrange)
CondMean_IV_p_MSE <- numeric(hrange)



for(h in 1:hrange){
  CondMean_p_MAE[h] <-my_mae(CondMean_p[,h], ActualValues[,h])
  CondMean_p_MSE[h] <-my_mse(CondMean_p[,h], ActualValues[,h])
  
  CondMean_IV_p_MAE[h] <-my_mae(CondMean_IV_p[,h], ActualValues[,h])
  CondMean_IV_p_MSE[h] <-my_mse(CondMean_IV_p[,h], ActualValues[,h])
  
}

print(CondMean_p_MAE/Naive_MAE)

print(CondMean_p_MSE/Naive_MSE)



print(CondMean_IV_p_MAE/Naive_MAE)

print(CondMean_IV_p_MSE/Naive_MSE)

####
print(CondMean_p_MAE/CondMean_MAE)

print(CondMean_p_MSE/CondMean_MSE)

plot(x[(is_length+1):n],type="l")
#lines(CondMean[,1])
lines(CondMean[,1],col="red")
lines(NaiveFC[,1],col="blue")

lines(CondMean_p[,1],col="brown")


plot(x[(is_length+1):((is_length+1)+100)],type="l")
#lines(CondMean[1:100,1])
lines(CondMean[1:100,1],col="red")
lines(NaiveFC[1:100,1],col="blue")

lines(CondMean_p[1:100,1],col="brown")

plot(x[(is_length+1):((is_length+1)+100)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[1:100,10],col="red")
lines(NaiveFC[1:100,10],col="blue")

lines(CondMean_p[1:100,10],col="brown")

plot(x[((is_length+1)+20):((is_length+1)+20+100)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[(1+20):(100+20),20],col="red")
lines(NaiveFC[(1+20):(100+20),20],col="blue")

lines(CondMean_p[(1+20):(100+20),20],col="brown")

plot(x[((is_length+1)+20):((is_length+1)+20+30)],type="l")
#lines(CondMean[1:100,10])
lines(CondMean[(1+20):(30+20),20],col="red")
lines(NaiveFC[(1+20):(30+20),20],col="blue")

lines(CondMean_p[(1+20):(30+20),20],col="brown")

#Plot the errors with ggplot
ratio1 <- CondMean_MSE/CondMean_p_MSE
ratio1
df1<-data.frame(lag=seq(1,20,1), ratio=ratio1)

df1p <-ggplot(data=df1, aes(x=lag, y=ratio))+
  geom_bar(stat="identity")+
  ylim(0,1.16)+
  geom_hline(yintercept=1, lty=2, col='red')+
  xlab("Lag")+
  ylab("Ratio of MSEs")
df1p
ggsave("CondMean-nonpar-par-MSE-ratio.eps", width = 20, height = 20, units = "cm")

ratio2 <- CondMean_MAE/CondMean_p_MAE
ratio2
df2<-data.frame(lag=seq(1,20,1), ratio=ratio2)

df2p <-ggplot(data=df2, aes(x=lag, y=ratio))+
  geom_bar(stat="identity")+
  ylim(0,1.16)+
  geom_hline(yintercept=1, lty=2, col='red')+
  xlab("Lag")+
  ylab("Ratio of MAEs")
df2p
ggsave("CondMean-nonpar-par-MAE-ratio.eps", width = 20, height = 20, units = "cm")

