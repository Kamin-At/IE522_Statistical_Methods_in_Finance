data = read.csv("C:/min/coding_project/IE522_Statistical_Methods_in_Finance/Week1/assignment/TSLA.csv")

log_ret = diff(log(data$Adj.Close))
tmp_ecdf = ecdf(log_ret)
plot(tmp_ecdf, xlab="log returns of TSLA stock (x)", ylab="Fn(x)")