#Code for the empirical study
#Last update: 28.07.2022
library(R.matlab)

#Read in the data
data <- readMat('A_20200504_data.mat')
str(data)

#Multiply data by 100 to convert to dollar cents
spread <-data$spr*100
plot(spread,type="l")

summary(spread)