#This file contains the R code for the empirical study presented in the article
#"Nonparametric estimation for trawl processes: Theory and Applications" 
#by Orimar Sauri (Aalborg University) and Almut Veraart (Imperial College London).
#Last update: 09.08.2022
library(R.matlab)
library(lubridate)
library(ggplot2)
library(dplyr)
library(forecast)
library(latex2exp)
library(reshape2)
library(gtools)
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
n <-length(x)


plot(x, type="l")
summary(x)
acf(x)

#Plot using ggplot
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
my_lag <- n-1
my_Delta <- 1/12 #sampling in 5s=1/12 min
esttrawlfct <- nonpar_trawlest(x, Delta=my_Delta, lag=my_lag)$a_hat
plot(esttrawlfct)

plot(esttrawlfct[1:100])
abline(h=0,col="red")


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
# We plot the asymptotic variance on our $\Delta$-grid:
c4 <- c4est(x, my_Delta)  
av_vector <-numeric(my_lag)
for(i in 1:my_lag){
  av_vector[i] <- asymptotic_variance_est((i-1)*my_Delta, c4, varlevyseed=1, my_Delta, esttrawlfct)$v
}

#Note that av_vector[1] contains the value for i=0 etc.
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
  geom_point(size=3, show.legend=FALSE)+
  scale_color_manual(values=c("y1"="blue","y2"="blue","y3"="black")) +
  scale_shape_manual(values = c(1, 1, 19))+
  xlab("l")+
  ylab(TeX("$\\hat{a}(\\cdot)$"))+
  theme(legend.position = "none")

g5
ggsave("CIs.eps", width = 20, height = 20, units = "cm")

which(y3<0)
#We note that the trawl function estimate is negative for 
#index 18, 22, 23, 25, 28, 31 etc.
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
NaiveFC <-matrix(0, nrow=noos, ncol=hrange)

ActualValues <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference <- matrix(0, nrow=noos, ncol=hrange)


for(tau in 1:noos){
  
  
  data <-x[tau:(is_length+tau-1)] # take new subsample of length is_length
  actualvalue <- x[(is_length+tau-1)]
  datamean <- mean(data) 
  
  #Estimate trawl function
  esttrawlfct_fc <- nonpar_trawlest(data, Delta=my_Delta, lag=my_lag)$a_hat 
  
  #Estimate the Lebesgue measure of the trawl set (components)
  lebA_fc <- LebA_est(data, my_Delta)
  
  
  for(h in 1:hrange){ #forecasting horizon loop in 1:hrange
    
    
    slices <- LebA_slice_est(data, my_Delta, h*my_Delta)
    
    #Estimate Leb(A intersection A_h)
    lebAintersection_fc <-slices$LebAintersection 
    
    #Estimate Leb(A \ A_h)
    lebAsetdiff_fc <-slices$LebAsetdifference 
    
    #ActualValue
    ActualValues[tau,h]<- x[(is_length+tau-1+h)]
    
    #Compute forecast(s)
    #naive
    NaiveFC[tau,h] <- actualvalue
    
    #Trawl prediction formula:
    Weight_Intersection[tau, h]<-lebAintersection_fc/lebA_fc
    Weight_SetDifference[tau,h]<-lebAsetdiff_fc/lebA_fc
    
    CondMean[tau,h] <- actualvalue*Weight_Intersection[tau, h]+datamean*Weight_SetDifference[tau,h]
    
    
  } #end for h
  
}#end for tau


#Computing the error measures

CondMean_MAE <- numeric(hrange)
CondMean_MSE <- numeric(hrange)

Naive_MAE <- numeric(hrange)
Naive_MSE <- numeric(hrange)


for(h in 1:hrange){
  CondMean_MAE[h] <-my_mae(CondMean[,h], ActualValues[,h])
  CondMean_MSE[h] <-my_mse(CondMean[,h], ActualValues[,h])
  
  Naive_MAE[h] <-my_mae(NaiveFC[,h], ActualValues[,h])
  Naive_MSE[h] <-my_mse(NaiveFC[,h], ActualValues[,h])
  
}

print(CondMean_MAE/Naive_MAE)

print(CondMean_MSE/Naive_MSE)


#Plot the errors
plot(CondMean_MSE/Naive_MSE, type="h")




## Comparing the nonparametric forecasts with a parametric forecast
##We fit a NegBin trawl model with supOU trawl function.
#We first fit such a model to the full sample.
data <- x

#Estimate trawl function parameters
esttrawlfct_fc <- fit_LMtrawl(data, Delta=my_Delta, GMMlag=10,plot=TRUE) 

my_alpha <- esttrawlfct_fc$alpha
my_alpha
my_H <- esttrawlfct_fc$H
my_H


#Estimate the Lebesgue measure of the trawl set (components)
lebA_fc_c1 <-my_alpha/(my_H-1)
lebA_fc_c1

#Estimate negative binomial parameters
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
my_Delta <-1/12 #sampling in 5s=1/12 min




hrange <- 20
noos <- n-is_length-hrange #=941 for is_length=3000


#Create forecast matrices:
CondMean_p <-matrix(0, nrow=noos, ncol=hrange)
CondMean_IV_p <-matrix(0, nrow=noos, ncol=hrange)

ActualValues <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection_p <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference_p <- matrix(0, nrow=noos, ncol=hrange)

set.seed(1)

for(tau in 1:noos){
  
  data <-x[tau:(is_length+tau-1)] # take new subsample of length is_length
  actualvalue <- x[(is_length+tau-1)]
   
  #Estimate trawl function parameters
  esttrawlfct_fc <- fit_LMtrawl(data, Delta=my_Delta, GMMlag=10)  
  
  my_alpha <- esttrawlfct_fc$alpha
  my_H <- esttrawlfct_fc$H
  
  #Estimate the Lebesgue measure of the trawl set (components)
  lebA_fc_c1 <-my_alpha/(my_H-1)
  
  #Estimate negative binomial parameters
  negbinpar <- fit_marginalNB(data, lebA_fc_c1)
  
  my_theta <- negbinpar$theta
  my_m <- negbinpar$m
  my_mtilde <- (1-my_theta)^2/my_theta
  
  my_c <- my_theta*my_m/(1-my_theta)^2
  
  lebA_fc <-lebA_fc_c1*my_c
  
  
  for(h in 1:hrange){ #forecasting horizon loop in 1:hrange
    
    #Estimate Leb(A intersection A_h)
    lebAintersection_fc <-my_c*my_alpha/(my_H-1)*(1+h*my_Delta/my_alpha)^(1-my_H)
    
    #Estimate Leb(A \ A_h)
    lebAsetdiff_fc <-lebA_fc-lebAintersection_fc
    
    #ActualValue
    ActualValues[tau,h]<- x[(is_length+tau-1+h)]
    
    #Compute forecast(s)
    #Trawl prediction formula
    Weight_Intersection_p[tau, h]<-lebAintersection_fc/lebA_fc
    Weight_SetDifference_p[tau,h]<-lebAsetdiff_fc/lebA_fc
    
    
    CondMean_p[tau,h] <- actualvalue*Weight_Intersection_p[tau, h]+(1-my_theta)*lebA_fc*Weight_SetDifference_p[tau,h]
    
    
  } #end for h
  
}#end for tau


#Next we compute the error measures:

CondMean_p_MAE <- numeric(hrange)
CondMean_p_MSE <- numeric(hrange)

for(h in 1:hrange){
  CondMean_p_MAE[h] <-my_mae(CondMean_p[,h], ActualValues[,h])
  CondMean_p_MSE[h] <-my_mse(CondMean_p[,h], ActualValues[,h])
}

print(CondMean_p_MAE/Naive_MAE)

print(CondMean_p_MSE/Naive_MSE)

####
print(CondMean_p_MAE/CondMean_MAE)

print(CondMean_p_MSE/CondMean_MSE)





#########################
#Forecasting using the ACF as weights
######################
is_length <- 3000 #in sample length

my_lag <- is_length-1
my_Delta <- 1/12 #sampling in 5s=1/12 min





hrange <- 20
noos <- n-is_length-hrange #=941 for is_length=3000
#Create forecast matrices:
CondMean_simple <-matrix(0, nrow=noos, ncol=hrange)

ActualValues <- matrix(0, nrow=noos, ncol=hrange)

Weight_Intersection_simple <- matrix(0, nrow=noos, ncol=hrange)
Weight_SetDifference_simple <- matrix(0, nrow=noos, ncol=hrange)


for(tau in 1:noos){
  
  
  data <-x[tau:(is_length+tau-1)] # take new subsample of length is_length
  actualvalue <- x[(is_length+tau-1)]
  datamean <- mean(data) 
  
  
  for(h in 1:hrange){ #forecasting horizon loop in 1:hrange
    
    
  
    #ActualValue
    ActualValues[tau,h]<- x[(is_length+tau-1+h)]
    
    #Compute forecast(s)
    #Trawl prediction formula
    my_acf <- acf(data, plot=FALSE)$acf[h+1]
    Weight_Intersection_simple[tau, h]<-my_acf
    Weight_SetDifference_simple[tau,h]<-1-my_acf
    
    CondMean_simple[tau
                    ,h] <- actualvalue*Weight_Intersection_simple[tau, h]+datamean*Weight_SetDifference_simple[tau,h]
    
   
    
  } #end for h
  
}#end for tau



#Computing the error measures

CondMean_simple_MAE <- numeric(hrange)
CondMean_simple_MSE <- numeric(hrange)



for(h in 1:hrange){
  CondMean_simple_MAE[h] <-my_mae(CondMean_simple[,h], ActualValues[,h])
  CondMean_simple_MSE[h] <-my_mse(CondMean_simple[,h], ActualValues[,h])
  
}

print(CondMean_simple_MAE/Naive_MAE)

print(CondMean_simple_MSE/Naive_MSE)

print(CondMean_simple_MAE/CondMean_MAE)

print(CondMean_simple_MSE/CondMean_MSE)


#Run Diebold-Mariano tests for forecast comparisons
#Create error matrices:
CondMean_Error <- CondMean-ActualValues 
Naive_Error <- NaiveFC - ActualValues
CondMean_p_Error <- CondMean_p - ActualValues
CondMean_simple_Error <- CondMean_simple - ActualValues

#Nonparametric trawl versus naive
CM_Naive1 <-1+numeric(hrange)
for(h in 2:hrange){
  print(h)
  dm<-dm.test(Naive_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =1)$p.value
  print(dm)
  CM_Naive1[h]<-dm
}
CM_Naive1_stars <- stars.pval(CM_Naive1)

#Nonparametric trawl versus naive
CM_Naive2 <-1+numeric(hrange)
for(h in 2:hrange){
  print(h)
  dm<-dm.test(Naive_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =2)$p.value
  print(dm)
  CM_Naive2[h]<-dm
}
CM_Naive2_stars <- stars.pval(CM_Naive2)


#Nonparametric trawl versus parametric trawl
CM_CMp1 <-numeric(hrange)
for(h in 1:hrange){
  print(h)
  dm<-dm.test(CondMean_p_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =1)$p.value
  print(dm)
  CM_CMp1[h]<-dm
}
CM_CMp1_stars <- stars.pval(CM_CMp1)

#Nonparametric trawl versus naive
CM_CMp2 <-numeric(hrange)
for(h in 1:hrange){
  print(h)
  dm<-dm.test(CondMean_p_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =2)$p.value
  print(dm)
  CM_CMp2[h]<-dm
}
CM_CMp2_stars <- stars.pval(CM_CMp2)





#Nonparametric trawl versus acf-based nonparametric trawl
#Power 1
#Nonparametric trawl versus parametric trawl
CM_CMsimple1 <-numeric(hrange)
for(h in 1:hrange){
  print(h)
  dm<-dm.test(CondMean_simple_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =1)$p.value
  print(dm)
  CM_CMsimple1[h]<-dm
}
CM_CMsimple1_stars <- stars.pval(CM_CMsimple1)

#Nonparametric trawl versus naive
CM_CMsimple2 <-numeric(hrange)
for(h in 1:hrange){
  print(h)
  dm<-dm.test(CondMean_simple_Error[,h], CondMean_Error[,h] , alternative="greater", h=h, power =2)$p.value
  print(dm)
  CM_CMsimple2[h]<-dm
}
CM_CMsimple2_stars <- stars.pval(CM_CMsimple2)




#####All error plots
#Plot the errors with ggplot
ratio1 <- CondMean_MSE/Naive_MSE
ratio1


df1<-data.frame(lag=seq(1,20,1), ratio=ratio1)

df1p <-ggplot(data=df1, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  #ylim(-1,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MSEs - 1")
df1p

label1.df <- data.frame(lag=seq(1,20,1), ratio=1+rep(0.03,20))

df1p+geom_text(data=label1.df, label=CM_Naive2_stars,size=10, angle=90)
ggsave("CondMean-Naive-MSE-ratio.eps", width = 20, height = 20, units = "cm")

ratio2 <- CondMean_MAE/Naive_MAE
ratio2
df2<-data.frame(lag=seq(1,20,1), ratio=ratio2)

df2p <-ggplot(data=df2, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  geom_bar(stat="identity")+
  #ylim(0,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MAEs - 1")
df2p



df2p+geom_text(data=label1.df, label=CM_Naive1_stars,size=10, angle=90)
ggsave("CondMean-Naive-MAE-ratio.eps", width = 20, height = 20, units = "cm")


################################

#Plot the errors with ggplot
ratio1 <- CondMean_MSE/CondMean_p_MSE
ratio1
df1<-data.frame(lag=seq(1,20,1), ratio=ratio1)

df1p <-ggplot(data=df1, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  #ylim(-1,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MSEs - 1")
df1p

label1.df <- data.frame(lag=seq(1,20,1), ratio=1+rep(0.03,20))

df1p+geom_text(data=label1.df, label=CM_CMp2_stars,size=10, angle=90)
ggsave("CondMean-nonpar-par-MSE-ratio.eps", width = 20, height = 20, units = "cm")

ratio2 <- CondMean_MAE/CondMean_p_MAE
ratio2
df2<-data.frame(lag=seq(1,20,1), ratio=ratio2)

df2p <-ggplot(data=df2, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  geom_bar(stat="identity")+
  #ylim(0,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MAEs - 1")
df2p



df2p+geom_text(data=label1.df, label=CM_CMp1_stars,size=10, angle=90)
ggsave("CondMean-nonpar-par-MAE-ratio.eps", width = 20, height = 20, units = "cm")

######################################################

#Plot the errors with ggplot
ratio1 <- CondMean_MSE/CondMean_simple_MSE
ratio1
df1<-data.frame(lag=seq(1,20,1), ratio=ratio1)

df1p <-ggplot(data=df1, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  #ylim(-1,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MSEs - 1")
df1p

label1.df <- data.frame(lag=seq(1,20,1), ratio=1+rep(0.03,20))

df1p+geom_text(data=label1.df, label=CM_CMsimple2_stars,size=10, angle=90)

ggsave("CondMean-simple-MSE-ratio.eps", width = 20, height = 20, units = "cm")

ratio2 <- CondMean_MAE/CondMean_simple_MAE
ratio2
df2<-data.frame(lag=seq(1,20,1), ratio=ratio2)

df2p <-ggplot(data=df2, aes(x=lag, y=ratio-1, fill=ifelse(ratio-1>0,"+","-")), show.legend=FALSE)+
  scale_fill_manual(values=c("blue","red"), name=" ") +
  geom_bar(stat="identity")+
  #ylim(0,1.2)+
  geom_hline(yintercept=0, lty=2, col='black')+
  xlab("Lag")+
  ylab("Ratio of MAEs - 1")
df2p



df2p+geom_text(data=label1.df, label=CM_CMsimple2_stars,size=10, angle=90)
ggsave("CondMean-simple-MAE-ratio.eps", width = 20, height = 20, units = "cm")



