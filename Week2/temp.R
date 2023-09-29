library(ISLR)
library(VGAM)
library(moments)


x = seq(0, 1, by = 0.01)

tmp = rlaplace(1000)
tmp2 = rlaplace(1000,location = 0,scale = 5)
tmp3 = rlaplace(1000,location = 0,scale = 0.1)
qqplot(tmp,tmp2,xlim=c(-6,6),ylim=c(-6,6))
qqline(tmp2,xlim=c(-6,6),ylim=c(-6,6))

qqplot(tmp,tmp3,xlim=c(-6,6),ylim=c(-6,6))
qqline(tmp3,xlim=c(-6,6),ylim=c(-6,6))
