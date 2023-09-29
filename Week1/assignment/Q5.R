data = c(1,2,3,4,5)
given_quantile = 0.7
print("data")
print(data)

C = given_quantile * (length(data) - 1) + 1
C_lower = floor(C)
C_upper = ceiling(C)
gamma_remainder = C - C_lower
solution = data[C_lower] * (1 - gamma_remainder) + data[C_upper] * (gamma_remainder)
print(paste(paste0("Q(", given_quantile, ") ="), solution))
print(paste(paste0("F(Q(", given_quantile, ")) ="), length(data[data<=solution])/length(data)))