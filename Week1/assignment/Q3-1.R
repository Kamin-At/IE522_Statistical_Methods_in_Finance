data = c(1,2,4,6,9,12)
given_quantile = 0.65
print("data")
print(data)

C = given_quantile * (length(data) - 1) + 1
C_lower = floor(C)
C_upper = ceiling(C)
gamma_remainder = C - C_lower
solution = data[C_lower] * (1 - gamma_remainder) + data[C_upper] * (gamma_remainder)
print(paste("solution:", solution))