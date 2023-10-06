library(moments)

mu0 = 0

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week6/TSLANVDA.csv")
log_ret_tsla = diff(log(data$TSLAAdjClose), lag=1)
log_ret_nvda = diff(log(data$NVDAAdjClose), lag=1)
D = log_ret_tsla - log_ret_nvda
se = sd(D)/sqrt(length(D))
cat("mu_D:", mean(D))
Z0 = (mean(D) - mu0)/se
cat("Z0:", Z0)
pval = pnorm(Z0)
cat("p-value:", pval)

hist(D, 40, freq = FALSE)
# The p-value of the one sided test is 0.2521, and the test has its z value of $$Z = -0.6678$$.