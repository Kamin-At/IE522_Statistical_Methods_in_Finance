#install.packages("dplyr")
# install.packages("moments")

library(dplyr)
library(moments)

data = read.csv("C:/Users/kamin/OneDrive/Documents/min/mfe/uiuc/Fall2023/Stats_for_Fin/Week1/assignment/TSLA.csv")
print(data)
data["lag_adj_close"] = dplyr::lag(data$Adj.Close, n= 1)
data["log_ret"] = with(data, log(data$Adj.Close/data$lag_adj_close))
hist(data$log_ret, breaks = 30)
print(paste("skewness:",skewness(data$log_ret, na.rm = TRUE)))
print(paste("kurtosis:",kurtosis(data$log_ret, na.rm = TRUE)))