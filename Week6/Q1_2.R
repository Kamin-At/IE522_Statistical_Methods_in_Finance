library(moments)

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week1/assignment/TSLA.csv")

log_ret = diff(log(data$Adj.Close))
kurtosis(log_ret)