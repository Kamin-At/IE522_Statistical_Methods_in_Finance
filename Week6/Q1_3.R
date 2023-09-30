library(moments)

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week1/assignment/TSLA.csv")

log_ret = diff(log(data$Adj.Close))

num_trial = 5000
num_sample = length(log_ret)#1257

sample_kurtoses = rep(0, num_trial)

for (i in 1: num_trial){
  xstar = sample(log_ret, num_trial, replace = TRUE)
  sample_kurtoses[i] = kurtosis(xstar)
}

hist(sample_kurtoses, 30, freq = FALSE)