library(moments)

mu0 = 0
alpha = 0.05

data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week6/TSLANVDA.csv")
log_ret_tsla = diff(log(data$TSLAAdjClose), lag=1)
log_ret_nvda = diff(log(data$NVDAAdjClose), lag=1)
D = log_ret_tsla - log_ret_nvda
se = sd(D)/sqrt(length(D))

Z0 = (mean(D) - mu0)/se

pval = pnorm(Z0)
pval
cat("p-value:", pval)
critical_value = qnorm(p=alpha)
cat("mu_D:", mean(D))
cat("Z0:", Z0)
cat("critical_value Z of D_bar:", critical_value)

# With 95% confidence level, we cannot reject the null hypothesis $$H_0$$ as $$\mu_D > critical_value of \mu_D$$ 