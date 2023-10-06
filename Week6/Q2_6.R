library(moments)

mu0 = 0
alpha = 0.05

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week6/TSLANVDA.csv")
log_ret_tsla = diff(log(data$TSLAAdjClose), lag=1)
log_ret_nvda = diff(log(data$NVDAAdjClose), lag=1)
D = log_ret_tsla - log_ret_nvda
se = sd(D)/sqrt(length(D))

lb = mean(D) - se*qnorm(p = 1 - alpha)
cat("The lower bound of D_bar:",lb)

# So, the confidence interval is $$\bar{D] \in [-0.00290, \inf)