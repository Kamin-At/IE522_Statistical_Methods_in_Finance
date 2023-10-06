library(moments)

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week6/TSLANVDA.csv")
log_ret_tsla = diff(log(data$TSLAAdjClose), lag=1)

shapiro.test(log_ret_tsla)
#For  Shapiro-Wilk test, the P-value is much less than 0.05, meaning that we accept the alternative hypothesis that the TSLA stock returns did not come from a normal distribution.
jarque.test(log_ret_tsla)
#For  Jarque-Bera test, the P-value is much less than 0.05, meaning that we accept the alternative hypothesis that the TSLA stock returns did not come from a normal distribution.

qqnorm(log_ret_tsla, xlab = "Theoretical Normal quantile", ylab = "Sample quantile")
qqline(log_ret_tsla)
# We can also clearly see that the TSLA return data did not follow a normal distribution from the normal Q-Q plot.