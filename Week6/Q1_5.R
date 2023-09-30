library(moments)
library(comprehenr)

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week1/assignment/TSLA.csv")

log_ret = diff(log(data$Adj.Close))

num_trial = 5000
num_sample = length(log_ret)#1257

sample_kurtoses = rep(0, num_trial)

for (i in 1: num_trial){
  xstar = sample(log_ret, num_trial, replace = TRUE)
  sample_kurtoses[i] = kurtosis(xstar)
}

cat("P(T>=6) = ",length(to_list(for (i in sample_kurtoses) if (i >= 6) i ))/length(sample_kurtoses))