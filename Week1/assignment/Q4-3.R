#install.packages("dplyr")

library(dplyr)

data = read.csv("C:/Users/kamin/OneDrive/Documents/min/mfe/uiuc/Fall2023/Stats_for_Fin/Week1/assignment/TSLA.csv")
data["lag_adj_close"] = dplyr::lag(data$Adj.Close, n= 1)
log_ret = with(data, diff(log(data$Adj.Close)))
data["log_ret"] = with(data, log(data$Adj.Close/data$lag_adj_close))
data["log_ret_lag1"] = dplyr::lag(data$log_ret, n= 1)
print(data)
print("Test for corr(log_ret, log_ret_lag1):")
print(cor.test(data$log_ret, data$log_ret_lag1, use = "complete.obs"))
plot(data$log_ret, data$log_ret_lag1)