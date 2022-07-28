#Code for the empirical study
#Last update: 28.07.2022
library(R.matlab)
library(lubridate)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(reshape2)
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
esttrawlfct <- nonpar_trawlest(x, Delta=my_Delta, lag=my_lag) 
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